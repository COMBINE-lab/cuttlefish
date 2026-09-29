//! lz4-compressed block streams for intermediates written once and read back
//! sequentially.
//!
//! A stream is a sequence of blocks, each an 8-byte header -- the raw length
//! and the stored length, both `u32` little-endian -- followed by the stored
//! bytes. The stored length equals the raw length when a block did not shrink
//! and was kept uncompressed; otherwise the bytes are an lz4 block. Blocks hold
//! at most [`BLOCK_BYTES`] raw bytes, so both ends buffer one block.
//!
//! The format is private to a build: files are written and read by the same
//! binary within one run, so it carries no magic or version.

use std::io::{self, Read, Write};

/// Raw bytes per block. Large enough for lz4 to find the repetition in
/// fixed-width records, small enough that a writer per bucket stays cheap.
pub(crate) const BLOCK_BYTES: usize = 256 * 1024;

const HEADER_BYTES: usize = 8;

/// Buffers writes into blocks and writes each one compressed.
///
/// Call [`Write::flush`] before dropping: it writes the final, partial block.
/// Dropping without it loses that block, as dropping a `BufWriter` whose
/// flush failed would.
pub(crate) struct Lz4BlockWriter<W: Write> {
    inner: W,
    raw: Vec<u8>,
    stored: Vec<u8>,
}

impl<W: Write> Lz4BlockWriter<W> {
    pub(crate) fn new(inner: W) -> Self {
        Self {
            inner,
            raw: Vec::with_capacity(BLOCK_BYTES),
            stored: Vec::new(),
        }
    }

    fn write_block(&mut self) -> io::Result<()> {
        if self.raw.is_empty() {
            return Ok(());
        }
        self.stored
            .resize(lz4_flex::block::get_maximum_output_size(self.raw.len()), 0);
        let compressed = lz4_flex::block::compress_into(&self.raw, &mut self.stored)
            .map_err(|error| io::Error::other(error.to_string()))?;
        let payload = if compressed < self.raw.len() {
            &self.stored[..compressed]
        } else {
            &self.raw[..]
        };
        let mut header = [0u8; HEADER_BYTES];
        header[..4].copy_from_slice(&(self.raw.len() as u32).to_le_bytes());
        header[4..].copy_from_slice(&(payload.len() as u32).to_le_bytes());
        self.inner.write_all(&header)?;
        self.inner.write_all(payload)?;
        self.raw.clear();
        Ok(())
    }
}

impl<W: Write> Write for Lz4BlockWriter<W> {
    fn write(&mut self, mut bytes: &[u8]) -> io::Result<usize> {
        let written = bytes.len();
        while !bytes.is_empty() {
            let take = (BLOCK_BYTES - self.raw.len()).min(bytes.len());
            self.raw.extend_from_slice(&bytes[..take]);
            bytes = &bytes[take..];
            if self.raw.len() == BLOCK_BYTES {
                self.write_block()?;
            }
        }
        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.write_block()?;
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
        if raw_len > BLOCK_BYTES || stored_len > lz4_flex::block::get_maximum_output_size(raw_len) {
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
        let mut writer = Lz4BlockWriter::new(Vec::new());
        for piece in data.chunks(chunk.max(1)) {
            writer.write_all(piece).unwrap();
        }
        writer.flush().unwrap();
        let mut out = Vec::new();
        Lz4BlockReader::new(&writer.inner[..])
            .read_to_end(&mut out)
            .unwrap();
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
