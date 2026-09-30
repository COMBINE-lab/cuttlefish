//! lz4-compressed block streams for intermediates written once and read back
//! sequentially.
//!
//! A stream is a sequence of blocks, each an 8-byte header -- the raw length
//! and the stored length, both `u32` little-endian -- followed by the stored
//! bytes. The stored length equals the raw length when a block did not shrink
//! and was kept uncompressed; otherwise the bytes are an lz4 block. A writer
//! picks its block size, at most [`MAX_BLOCK_BYTES`]; the reader buffers one
//! block. The compression scratch is per thread, so a writer holds only its
//! pending block -- no more than the `BufWriter` it replaces.
//!
//! The format is private to a build: files are written and read by the same
//! binary within one run, so it carries no magic or version.

use std::io::{self, Read, Write};

/// Default raw bytes per block. Large enough for lz4 to find the repetition
/// in fixed-width records, small enough that a writer per bucket stays cheap.
pub(crate) const BLOCK_BYTES: usize = 256 * 1024;

/// Largest block a reader accepts.
pub(crate) const MAX_BLOCK_BYTES: usize = 1024 * 1024;

/// Whether new writers compress their blocks. Each build sets it from its
/// own parameters ([`crate::intermediates::apply_params`]) before opening any
/// stream; on is the default for streams opened outside a build. With it
/// off, blocks are stored raw and the format -- and the reader -- are
/// unchanged.
static COMPRESS_BLOCKS: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(true);

pub(crate) fn set_compress_blocks(compress: bool) {
    COMPRESS_BLOCKS.store(compress, std::sync::atomic::Ordering::Relaxed);
}

thread_local! {
    static COMPRESSED: std::cell::RefCell<Vec<u8>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// Writes one block of `raw` -- compressed, or as is if that is no smaller
/// or compression is off.
fn write_block(inner: &mut impl Write, raw: &[u8], compress: bool) -> io::Result<()> {
    debug_assert!(!raw.is_empty() && raw.len() <= MAX_BLOCK_BYTES);
    if !compress {
        let mut header = [0u8; HEADER_BYTES];
        header[..4].copy_from_slice(&(raw.len() as u32).to_le_bytes());
        header[4..].copy_from_slice(&(raw.len() as u32).to_le_bytes());
        inner.write_all(&header)?;
        return inner.write_all(raw);
    }
    COMPRESSED.with_borrow_mut(|stored| {
        stored.resize(lz4_flex::block::get_maximum_output_size(raw.len()), 0);
        let compressed = lz4_flex::block::compress_into(raw, stored)
            .map_err(|error| io::Error::other(error.to_string()))?;
        let payload = if compressed < raw.len() {
            &stored[..compressed]
        } else {
            raw
        };
        let mut header = [0u8; HEADER_BYTES];
        header[..4].copy_from_slice(&(raw.len() as u32).to_le_bytes());
        header[4..].copy_from_slice(&(payload.len() as u32).to_le_bytes());
        inner.write_all(&header)?;
        inner.write_all(payload)
    })
}

const HEADER_BYTES: usize = 8;

/// Buffers writes into blocks and writes each one compressed.
///
/// Call [`Write::flush`] before dropping: it writes the final, partial block.
/// Dropping without it loses that block, as dropping a `BufWriter` whose
/// flush failed would.
pub(crate) struct Lz4BlockWriter<W: Write> {
    inner: W,
    raw: Vec<u8>,
    block_bytes: usize,
    compress: bool,
}

impl<W: Write> Lz4BlockWriter<W> {
    pub(crate) fn new(inner: W) -> Self {
        Self::with_block_bytes(inner, BLOCK_BYTES)
    }

    pub(crate) fn with_block_bytes(inner: W, block_bytes: usize) -> Self {
        assert!((1..=MAX_BLOCK_BYTES).contains(&block_bytes));
        Self {
            inner,
            raw: Vec::new(),
            block_bytes,
            compress: COMPRESS_BLOCKS.load(std::sync::atomic::Ordering::Relaxed),
        }
    }

    /// Overrides the process-wide setting for this writer.
    #[cfg(test)]
    pub(crate) fn with_compression(mut self, compress: bool) -> Self {
        self.compress = compress;
        self
    }

    /// Writes `bytes` straight out as blocks, without copying them into the
    /// pending block: for callers that already batch their own writes.
    pub(crate) fn write_blocks(&mut self, bytes: &[u8]) -> io::Result<()> {
        self.write_pending()?;
        for block in bytes.chunks(self.block_bytes) {
            write_block(&mut self.inner, block, self.compress)?;
        }
        Ok(())
    }

    /// The underlying writer, after writing any pending block.
    #[cfg(test)]
    pub(crate) fn into_inner(mut self) -> W {
        self.write_pending().expect("pending block written");
        self.inner
    }

    fn write_pending(&mut self) -> io::Result<()> {
        if self.raw.is_empty() {
            return Ok(());
        }
        write_block(&mut self.inner, &self.raw, self.compress)?;
        self.raw.clear();
        Ok(())
    }
}

impl<W: Write> Write for Lz4BlockWriter<W> {
    fn write(&mut self, mut bytes: &[u8]) -> io::Result<usize> {
        let written = bytes.len();
        while !bytes.is_empty() {
            let take = (self.block_bytes - self.raw.len()).min(bytes.len());
            self.raw.extend_from_slice(&bytes[..take]);
            bytes = &bytes[take..];
            if self.raw.len() == self.block_bytes {
                self.write_pending()?;
            }
        }
        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.write_pending()?;
        self.inner.flush()
    }
}

/// Reads a stream written by [`Lz4BlockWriter`].
pub(crate) struct Lz4BlockReader<R: Read> {
    inner: R,
    raw: Vec<u8>,
    position: usize,
    stored: Vec<u8>,
}

impl<R: Read> Lz4BlockReader<R> {
    pub(crate) fn new(inner: R) -> Self {
        Self {
            inner,
            raw: Vec::new(),
            position: 0,
            stored: Vec::new(),
        }
    }

    /// Loads the next block; `false` at a clean end of stream.
    fn next_block(&mut self) -> io::Result<bool> {
        let mut header = [0u8; HEADER_BYTES];
        let mut filled = 0;
        while filled < HEADER_BYTES {
            match self.inner.read(&mut header[filled..])? {
                0 if filled == 0 => return Ok(false),
                0 => return Err(io::ErrorKind::UnexpectedEof.into()),
                n => filled += n,
            }
        }
        let raw_len = u32::from_le_bytes(header[..4].try_into().expect("u32")) as usize;
        let stored_len = u32::from_le_bytes(header[4..].try_into().expect("u32")) as usize;
        // The writer never emits an empty block, and stores a block
        // compressed only when that is smaller. Accepting either shape would
        // let zero-filled or torn tails pass for a clean end of stream.
        if raw_len == 0 || raw_len > MAX_BLOCK_BYTES || stored_len > raw_len {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "malformed lz4 block header",
            ));
        }
        self.raw.resize(raw_len, 0);
        self.position = 0;
        if stored_len == raw_len {
            self.inner.read_exact(&mut self.raw)?;
            return Ok(true);
        }
        self.stored.resize(stored_len, 0);
        self.inner.read_exact(&mut self.stored)?;
        let decoded = lz4_flex::block::decompress_into(&self.stored, &mut self.raw)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))?;
        if decoded != raw_len {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "lz4 block decoded to the wrong length",
            ));
        }
        Ok(true)
    }
}

impl<R: Read> Read for Lz4BlockReader<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if out.is_empty() {
            return Ok(0);
        }
        while self.position == self.raw.len() {
            if !self.next_block()? {
                return Ok(0);
            }
        }
        let take = (self.raw.len() - self.position).min(out.len());
        out[..take].copy_from_slice(&self.raw[self.position..self.position + take]);
        self.position += take;
        Ok(take)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn round_trip(data: &[u8], chunk: usize) -> Vec<u8> {
        let mut out = Vec::new();
        // Compressed and raw blocks read back through the same reader.
        for compress in [true, false] {
            let mut writer = Lz4BlockWriter::new(Vec::new()).with_compression(compress);
            for piece in data.chunks(chunk.max(1)) {
                writer.write_all(piece).unwrap();
            }
            writer.flush().unwrap();
            if !compress {
                assert_eq!(
                    writer.inner.len(),
                    data.len() + data.len().div_ceil(BLOCK_BYTES) * HEADER_BYTES
                );
            }
            out.clear();
            Lz4BlockReader::new(&writer.inner[..])
                .read_to_end(&mut out)
                .unwrap();
            assert_eq!(out, data);
        }
        out
    }

    #[test]
    fn streams_round_trip_across_block_boundaries() {
        let mut state = 7u64;
        // Compressible records, then incompressible noise, then empty.
        let mut data: Vec<u8> = (0..3 * BLOCK_BYTES + 123)
            .map(|i| if i % 8 < 5 { (i / 8 % 251) as u8 } else { 0 })
            .collect();
        data.extend((0..BLOCK_BYTES + 7).map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state as u8
        }));
        for chunk in [1, 7, 4096, BLOCK_BYTES, 3 * BLOCK_BYTES] {
            assert_eq!(round_trip(&data, chunk), data, "chunk {chunk}");
        }
        assert!(round_trip(&[], 1).is_empty());
    }

    #[test]
    fn direct_blocks_and_small_blocks_interleave() {
        let data: Vec<u8> = (0..300_000u32).map(|i| (i % 7 * 3) as u8).collect();
        let mut writer = Lz4BlockWriter::with_block_bytes(Vec::new(), 4096);
        writer.write_all(&data[..10]).unwrap();
        writer.write_blocks(&data[10..200_000]).unwrap();
        writer.write_all(&data[200_000..200_001]).unwrap();
        writer.write_blocks(&data[200_001..]).unwrap();
        writer.flush().unwrap();
        let mut out = Vec::new();
        Lz4BlockReader::new(&writer.inner[..])
            .read_to_end(&mut out)
            .unwrap();
        assert_eq!(out, data);
    }

    #[test]
    fn headers_the_writer_never_emits_are_errors() {
        let header = |raw: u32, stored: u32| {
            let mut bytes = raw.to_le_bytes().to_vec();
            bytes.extend_from_slice(&stored.to_le_bytes());
            bytes.resize(bytes.len() + stored as usize, 0);
            bytes
        };
        for bad in [
            header(0, 0),
            header(16, 17),
            header(MAX_BLOCK_BYTES as u32 + 1, 8),
            vec![0u8; 64],
        ] {
            let error = Lz4BlockReader::new(&bad[..])
                .read_to_end(&mut Vec::new())
                .unwrap_err();
            assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        }
    }

    #[test]
    fn truncated_streams_are_errors() {
        let mut writer = Lz4BlockWriter::new(Vec::new());
        writer.write_all(&vec![3u8; 10_000]).unwrap();
        writer.flush().unwrap();
        let bytes = writer.inner;
        for cut in [1, HEADER_BYTES, bytes.len() - 1] {
            let mut out = Vec::new();
            assert!(
                Lz4BlockReader::new(&bytes[..cut])
                    .read_to_end(&mut out)
                    .is_err(),
                "cut at {cut}"
            );
        }
    }
}
