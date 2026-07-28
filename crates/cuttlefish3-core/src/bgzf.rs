//! Parallel decompression of block-structured (BGZF) gzip input.
//!
//! A plain gzip member is one LZ77 stream whose deflate block boundaries are
//! neither byte-aligned nor indexed, so it cannot be split without speculative
//! decoding. BGZF instead concatenates independent gzip members and records each
//! member's compressed length in an extra-field subfield, so members can be
//! located by reading headers alone and inflated concurrently.
//!
//! [`ParallelBgzfReader`] exposes that as a plain [`Read`], producing exactly the
//! bytes a serial decoder would. Parsing therefore stays unchanged, and only the
//! decompression step — the actual bottleneck on large inputs — is parallel.

use flate2::read::GzDecoder;
use std::io::{self, Read};
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};

/// Fixed portion of a gzip header, up to and including `XLEN`.
const GZIP_FIXED_HEADER: usize = 12;
/// `FLG.FEXTRA`, which BGZF always sets.
const FEXTRA: u8 = 0x04;
/// Decompressed bytes in a BGZF block never exceed 64 KiB.
const MAX_BLOCK_PAYLOAD: usize = 64 * 1024;
/// Blocks in flight per worker; bounds memory to workers * this * 64 KiB.
const BLOCKS_IN_FLIGHT: usize = 4;

/// Returns the total size of the BGZF block whose header starts `header`.
///
/// Returns `None` when the bytes are not a BGZF block header, which is how a
/// plain gzip stream is distinguished from a block-structured one.
fn bgzf_block_size(header: &[u8]) -> Option<usize> {
    if header.len() < GZIP_FIXED_HEADER {
        return None;
    }
    // Magic, DEFLATE method, and an extra field are all required.
    if header[0] != 0x1f || header[1] != 0x8b || header[2] != 0x08 || header[3] & FEXTRA == 0 {
        return None;
    }
    let extra_len = u16::from_le_bytes([header[10], header[11]]) as usize;
    if header.len() < GZIP_FIXED_HEADER + extra_len {
        return None;
    }
    let extra = &header[GZIP_FIXED_HEADER..GZIP_FIXED_HEADER + extra_len];
    // Scan subfields for `BC`, which carries the block's total size minus one.
    let mut offset = 0usize;
    while offset + 4 <= extra.len() {
        let slen = u16::from_le_bytes([extra[offset + 2], extra[offset + 3]]) as usize;
        if extra[offset] == b'B' && extra[offset + 1] == b'C' && slen == 2 {
            if offset + 6 > extra.len() {
                return None;
            }
            let bsize = u16::from_le_bytes([extra[offset + 4], extra[offset + 5]]) as usize;
            return Some(bsize + 1);
        }
        offset += 4 + slen;
    }
    None
}

/// Returns whether `head` begins a BGZF stream.
pub fn is_bgzf(head: &[u8]) -> bool {
    bgzf_block_size(head).is_some()
}

/// Number of leading bytes [`is_bgzf`] needs.
pub const PROBE_BYTES: usize = 64;

/// A [`Read`] that inflates BGZF blocks on a worker pool and returns the
/// decompressed stream in order.
pub struct ParallelBgzfReader {
    /// One ordered channel per worker; blocks are dealt round-robin, so reading
    /// the channels in the same rotation restores the original block order.
    outputs: Vec<Receiver<io::Result<Vec<u8>>>>,
    next: usize,
    current: Vec<u8>,
    position: usize,
    finished: bool,
    workers: Vec<std::thread::JoinHandle<()>>,
    dispatcher: Option<std::thread::JoinHandle<()>>,
}

impl ParallelBgzfReader {
    /// Starts a reader over `source`, inflating with `workers` threads.
    ///
    /// `source` must be positioned at the start of the stream.
    pub fn new<R>(mut source: R, workers: usize) -> Self
    where
        R: Read + Send + 'static,
    {
        let workers = workers.max(1);
        let mut block_txs = Vec::with_capacity(workers);
        let mut outputs = Vec::with_capacity(workers);
        let mut handles = Vec::with_capacity(workers);

        for _ in 0..workers {
            let (block_tx, block_rx) = sync_channel::<Vec<u8>>(BLOCKS_IN_FLIGHT);
            let (out_tx, out_rx) = sync_channel::<io::Result<Vec<u8>>>(BLOCKS_IN_FLIGHT);
            block_txs.push(block_tx);
            outputs.push(out_rx);
            handles.push(std::thread::spawn(move || {
                while let Ok(block) = block_rx.recv() {
                    let mut payload = Vec::with_capacity(MAX_BLOCK_PAYLOAD);
                    let result = GzDecoder::new(block.as_slice())
                        .read_to_end(&mut payload)
                        .map(|_| payload);
                    if out_tx.send(result).is_err() {
                        break;
                    }
                }
            }));
        }

        let dispatcher = std::thread::spawn(move || {
            dispatch_blocks(&mut source, &block_txs);
        });

        Self {
            outputs,
            next: 0,
            current: Vec::new(),
            position: 0,
            finished: false,
            workers: handles,
            dispatcher: Some(dispatcher),
        }
    }

    /// Pulls the next decompressed block, preserving stream order.
    fn advance(&mut self) -> io::Result<bool> {
        if self.finished {
            return Ok(false);
        }
        let slot = self.next % self.outputs.len();
        self.next += 1;
        match self.outputs[slot].recv() {
            Ok(Ok(payload)) => {
                self.current = payload;
                self.position = 0;
                Ok(true)
            }
            Ok(Err(error)) => {
                self.finished = true;
                Err(error)
            }
            // A closed channel means the dispatcher stopped handing out blocks.
            Err(_) => {
                self.finished = true;
                Ok(false)
            }
        }
    }
}

/// Splits `source` into whole BGZF blocks and deals them round-robin.
fn dispatch_blocks<R: Read>(source: &mut R, block_txs: &[SyncSender<Vec<u8>>]) {
    let mut next = 0usize;
    loop {
        let mut header = vec![0u8; GZIP_FIXED_HEADER];
        match read_exact_or_eof(source, &mut header) {
            Ok(0) => return,
            Ok(n) if n < GZIP_FIXED_HEADER => return,
            Ok(_) => {}
            Err(_) => return,
        }
        let extra_len = u16::from_le_bytes([header[10], header[11]]) as usize;
        header.resize(GZIP_FIXED_HEADER + extra_len, 0);
        if source.read_exact(&mut header[GZIP_FIXED_HEADER..]).is_err() {
            return;
        }
        let Some(total) = bgzf_block_size(&header) else {
            return;
        };
        if total <= header.len() {
            return;
        }
        let mut block = header;
        block.resize(total, 0);
        let filled = GZIP_FIXED_HEADER + extra_len;
        if source.read_exact(&mut block[filled..]).is_err() {
            return;
        }
        if block_txs[next % block_txs.len()].send(block).is_err() {
            return;
        }
        next += 1;
    }
}

fn read_exact_or_eof<R: Read>(source: &mut R, buffer: &mut [u8]) -> io::Result<usize> {
    let mut filled = 0;
    while filled < buffer.len() {
        match source.read(&mut buffer[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => continue,
            Err(error) => return Err(error),
        }
    }
    Ok(filled)
}

impl Read for ParallelBgzfReader {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        while self.position == self.current.len() {
            if !self.advance()? {
                return Ok(0);
            }
        }
        let take = (self.current.len() - self.position).min(out.len());
        out[..take].copy_from_slice(&self.current[self.position..self.position + take]);
        self.position += take;
        Ok(take)
    }
}

impl Drop for ParallelBgzfReader {
    fn drop(&mut self) {
        // Dropping the receivers unblocks the workers, which in turn lets the
        // dispatcher's sends fail and its thread exit.
        self.outputs.clear();
        if let Some(dispatcher) = self.dispatcher.take() {
            let _ = dispatcher.join();
        }
        for worker in self.workers.drain(..) {
            let _ = worker.join();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    /// Builds a BGZF stream by emitting each chunk as its own gzip member with
    /// the `BC` extra-field subfield BGZF requires.
    fn encode_bgzf(payload: &[u8], block_size: usize) -> Vec<u8> {
        let mut out = Vec::new();
        for chunk in payload.chunks(block_size.max(1)) {
            let mut deflated = Vec::new();
            let mut encoder = GzEncoder::new(&mut deflated, Compression::default());
            encoder.write_all(chunk).unwrap();
            encoder.finish().unwrap();
            // Re-wrap the member with an extra field carrying its total size.
            let body = &deflated[10..];
            let total = 12 + 6 + body.len();
            out.extend_from_slice(&[0x1f, 0x8b, 0x08, FEXTRA, 0, 0, 0, 0, 0, 0xff]);
            out.extend_from_slice(&6u16.to_le_bytes());
            out.extend_from_slice(b"BC");
            out.extend_from_slice(&2u16.to_le_bytes());
            out.extend_from_slice(&((total - 1) as u16).to_le_bytes());
            out.extend_from_slice(body);
        }
        out
    }

    #[test]
    fn detects_bgzf_and_rejects_plain_gzip() {
        let bgzf = encode_bgzf(b"ACGTACGTACGT", 4);
        assert!(is_bgzf(&bgzf));

        let mut plain = Vec::new();
        let mut encoder = GzEncoder::new(&mut plain, Compression::default());
        encoder.write_all(b"ACGTACGTACGT").unwrap();
        encoder.finish().unwrap();
        assert!(!is_bgzf(&plain));
    }

    #[test]
    fn parallel_reader_reproduces_the_serial_stream() {
        let payload: Vec<u8> = (0..200_000u32)
            .map(|index| b"ACGT"[(index % 4) as usize])
            .collect();
        let bgzf = encode_bgzf(&payload, 8 * 1024);
        for workers in [1usize, 2, 8] {
            let mut decoded = Vec::new();
            ParallelBgzfReader::new(std::io::Cursor::new(bgzf.clone()), workers)
                .read_to_end(&mut decoded)
                .unwrap();
            assert_eq!(decoded, payload, "mismatch with {workers} worker(s)");
        }
    }

    #[test]
    fn parallel_reader_handles_an_empty_stream() {
        let mut decoded = Vec::new();
        ParallelBgzfReader::new(std::io::Cursor::new(Vec::new()), 4)
            .read_to_end(&mut decoded)
            .unwrap();
        assert!(decoded.is_empty());
    }
}
