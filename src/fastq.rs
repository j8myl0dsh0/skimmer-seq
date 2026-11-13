//! FASTQ record parsing and writing.
//!
//! One record is four lines: `@name`, sequence, `+`, quality. Parsing validates
//! the shape of every record before any of it reaches the filter, so a truncated
//! file fails loudly instead of silently emitting half a record.

use std::fmt;
use std::io::{self, BufRead, Write};

/// A single FASTQ record. Sequences stay as bytes; nothing is assumed about
/// UTF-8 or about the alphabet beyond `ACGTN`-ish content.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    /// Header text after the leading `@`.
    pub name: String,
    /// Raw sequence bytes.
    pub seq: Vec<u8>,
    /// Raw Phred+33 quality bytes, same length as `seq`.
    pub qual: Vec<u8>,
}

/// Why a record could not be read.
#[derive(Debug)]
pub enum ReadError {
    /// The underlying reader failed.
    Io(io::Error),
    /// One of the four lines was missing.
    Truncated { line: usize },
    /// The header line did not start with `@`.
    BadHeader { line: usize },
    /// Sequence and quality lengths disagree, which makes the rest of the file
    /// ambiguous.
    LengthMismatch { name: String, seq: usize, qual: usize },
}

impl fmt::Display for ReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReadError::Io(e) => write!(f, "read failed: {e}"),
            ReadError::Truncated { line } => write!(f, "record truncated at line {line}"),
            ReadError::BadHeader { line } => write!(f, "expected '@' header at line {line}"),
            ReadError::LengthMismatch { name, seq, qual } => write!(
                f,
                "record {name:?}: sequence is {seq} bp but quality is {qual} bp"
            ),
        }
    }
}

impl std::error::Error for ReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ReadError::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for ReadError {
    fn from(e: io::Error) -> Self {
        ReadError::Io(e)
    }
}

/// Pull-based record reader over any buffered byte source.
pub struct Reader<R: BufRead> {
    inner: R,
    line: usize,
}

impl<R: BufRead> Reader<R> {
    pub fn new(inner: R) -> Self {
        Reader { inner, line: 0 }
    }

    /// Reads one line without its terminator. Returns the byte count so the
    /// caller can detect a clean EOF on an empty first line.
    fn line(&mut self, buf: &mut Vec<u8>) -> Result<usize, ReadError> {
        buf.clear();
        let n = self.inner.read_until(b'\n', buf)?;
        self.line += 1;
        if buf.last() == Some(&b'\n') {
            buf.pop();
        }
        if buf.last() == Some(&b'\r') {
            buf.pop();
        }
        Ok(n)
    }

    /// Reads the next record. `Ok(None)` at clean EOF; partial records are an
    /// error.
    pub fn next_record(&mut self) -> Result<Option<Record>, ReadError> {
        let mut header = Vec::new();
        if self.line(&mut header)? == 0 {
            return Ok(None);
        }
        if !header.starts_with(b"@") {
            return Err(ReadError::BadHeader { line: self.line });
        }

        let mut seq = Vec::new();
        if self.line(&mut seq)? == 0 {
            return Err(ReadError::Truncated { line: self.line });
        }

        let mut plus = Vec::new();
        if self.line(&mut plus)? == 0 || !plus.starts_with(b"+") {
            return Err(ReadError::Truncated { line: self.line });
        }

        let mut qual = Vec::new();
        if self.line(&mut qual)? == 0 {
            return Err(ReadError::Truncated { line: self.line });
        }

        if seq.len() != qual.len() {
            return Err(ReadError::LengthMismatch {
                name: String::from_utf8_lossy(&header[1..]).into_owned(),
                seq: seq.len(),
                qual: qual.len(),
            });
        }

        Ok(Some(Record {
            name: String::from_utf8_lossy(&header[1..]).trim().to_string(),
            seq,
            qual,
        }))
    }
}

/// Writes a record back out in FASTQ form with a bare `+` separator.
pub fn write_record<W: Write>(out: &mut W, rec: &Record) -> io::Result<()> {
    out.write_all(b"@")?;
    out.write_all(rec.name.as_bytes())?;
    out.write_all(b"\n")?;
    out.write_all(&rec.seq)?;
    out.write_all(b"\n+\n")?;
    out.write_all(&rec.qual)?;
    out.write_all(b"\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn round_trips_a_record() {
        let input = b"@r1\nACGT\n+\nIIII\n";
        let mut reader = Reader::new(Cursor::new(&input[..]));
        let rec = reader.next_record().unwrap().unwrap();
        assert_eq!(rec.name, "r1");
        assert_eq!(rec.seq, b"ACGT".to_vec());

        let mut out = Vec::new();
        write_record(&mut out, &rec).unwrap();
        assert_eq!(out, input.to_vec());
    }

    #[test]
    fn stops_cleanly_at_eof() {
        let mut reader = Reader::new(Cursor::new(&b""[..]));
        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn rejects_length_mismatch() {
        let mut reader = Reader::new(Cursor::new(&b"@r1\nACGT\n+\nII\n"[..]));
        match reader.next_record() {
            Err(ReadError::LengthMismatch { seq, qual, .. }) => {
                assert_eq!((seq, qual), (4, 2));
            }
            other => panic!("expected length mismatch, got {other:?}"),
        }
    }

    #[test]
    fn rejects_missing_quality_line() {
        let mut reader = Reader::new(Cursor::new(&b"@r1\nACGT\n+\n"[..]));
        assert!(matches!(
            reader.next_record(),
            Err(ReadError::Truncated { .. })
        ));
    }
}
