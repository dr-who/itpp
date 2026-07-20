//! Little-endian binary primitives + LEB128 varints for the container.

/// Cursor over a byte slice with checked reads.
pub struct Reader<'a> {
    data: &'a [u8],
    pub pos: usize,
}

#[derive(Debug)]
pub enum FormatError {
    Truncated,
    BadMagic,
    BadVersion(u16),
    Utf8,
    MissingSection(&'static str),
}

impl std::fmt::Display for FormatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FormatError::Truncated => write!(f, "unexpected end of container"),
            FormatError::BadMagic => write!(f, "not an ITPP container (bad magic)"),
            FormatError::BadVersion(v) => write!(f, "unsupported container version {v}"),
            FormatError::Utf8 => write!(f, "invalid UTF-8 in a string field"),
            FormatError::MissingSection(s) => write!(f, "required section {s} missing"),
        }
    }
}

impl std::error::Error for FormatError {}

pub type Result<T> = std::result::Result<T, FormatError>;

impl<'a> Reader<'a> {
    #[must_use]
    pub fn new(data: &'a [u8]) -> Self {
        Reader { data, pos: 0 }
    }

    pub fn seek(&mut self, pos: usize) {
        self.pos = pos;
    }

    pub fn remaining(&self) -> usize {
        self.data.len().saturating_sub(self.pos)
    }

    pub fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        if self.pos + n > self.data.len() {
            return Err(FormatError::Truncated);
        }
        let s = &self.data[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }

    pub fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }

    pub fn u16(&mut self) -> Result<u16> {
        Ok(u16::from_le_bytes(self.take(2)?.try_into().unwrap()))
    }

    pub fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }

    pub fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }

    pub fn tag(&mut self) -> Result<[u8; 4]> {
        Ok(self.take(4)?.try_into().unwrap())
    }

    /// Unsigned LEB128.
    pub fn varint(&mut self) -> Result<u64> {
        let mut result: u64 = 0;
        let mut shift = 0;
        loop {
            let byte = self.u8()?;
            result |= u64::from(byte & 0x7F) << shift;
            if byte & 0x80 == 0 {
                break;
            }
            shift += 7;
            if shift >= 64 {
                return Err(FormatError::Truncated);
            }
        }
        Ok(result)
    }

    /// varint length prefix + bytes.
    pub fn blob(&mut self) -> Result<&'a [u8]> {
        let n = self.varint()? as usize;
        self.take(n)
    }

    pub fn string(&mut self) -> Result<String> {
        let b = self.blob()?;
        std::str::from_utf8(b).map(str::to_string).map_err(|_| FormatError::Utf8)
    }
}

/// Append helpers on a byte vec.
pub trait Writer {
    fn put_u8(&mut self, v: u8);
    fn put_u16(&mut self, v: u16);
    fn put_u32(&mut self, v: u32);
    fn put_u64(&mut self, v: u64);
    fn put_tag(&mut self, tag: &[u8; 4]);
    fn put_varint(&mut self, v: u64);
    fn put_blob(&mut self, bytes: &[u8]);
    fn put_string(&mut self, s: &str);
}

impl Writer for Vec<u8> {
    fn put_u8(&mut self, v: u8) {
        self.push(v);
    }
    fn put_u16(&mut self, v: u16) {
        self.extend_from_slice(&v.to_le_bytes());
    }
    fn put_u32(&mut self, v: u32) {
        self.extend_from_slice(&v.to_le_bytes());
    }
    fn put_u64(&mut self, v: u64) {
        self.extend_from_slice(&v.to_le_bytes());
    }
    fn put_tag(&mut self, tag: &[u8; 4]) {
        self.extend_from_slice(tag);
    }
    fn put_varint(&mut self, mut v: u64) {
        loop {
            let mut byte = (v & 0x7F) as u8;
            v >>= 7;
            if v != 0 {
                byte |= 0x80;
            }
            self.push(byte);
            if v == 0 {
                break;
            }
        }
    }
    fn put_blob(&mut self, bytes: &[u8]) {
        self.put_varint(bytes.len() as u64);
        self.extend_from_slice(bytes);
    }
    fn put_string(&mut self, s: &str) {
        self.put_blob(s.as_bytes());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn varint_roundtrip() {
        for v in [0u64, 1, 127, 128, 300, 16384, u32::MAX as u64, u64::MAX] {
            let mut buf = Vec::new();
            buf.put_varint(v);
            let mut r = Reader::new(&buf);
            assert_eq!(r.varint().unwrap(), v);
        }
    }

    #[test]
    fn blob_and_string() {
        let mut buf = Vec::new();
        buf.put_string("hello");
        buf.put_blob(&[1, 2, 3]);
        let mut r = Reader::new(&buf);
        assert_eq!(r.string().unwrap(), "hello");
        assert_eq!(r.blob().unwrap(), &[1, 2, 3]);
    }
}
