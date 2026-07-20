//! Nucleotide sequence storage.
//!
//! We keep sequences as raw ASCII bytes (uppercased) rather than a packed 2-bit form at the
//! type level: losslessness demands we round-trip `N`, IUPAC ambiguity codes, and anything
//! else the input contains. The 2-bit packing is a *coder* (a baseline), not the storage
//! model — see `itpp-codec`.

/// A nucleotide sequence: ASCII, uppercased on construction.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct Sequence(Vec<u8>);

impl Sequence {
    #[must_use]
    pub fn from_bytes(mut bytes: Vec<u8>) -> Self {
        for b in &mut bytes {
            b.make_ascii_uppercase();
        }
        Sequence(bytes)
    }

    #[must_use]
    #[allow(clippy::should_implement_trait)] // infallible convenience ctor, not the fallible FromStr
    pub fn from_str(s: &str) -> Self {
        Sequence::from_bytes(s.as_bytes().to_vec())
    }

    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        &self.0
    }

    #[must_use]
    pub fn into_bytes(self) -> Vec<u8> {
        self.0
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

/// Complement a single base, preserving case-folded IUPAC semantics for the common cases.
/// Unknown bytes map to `N` on the reverse strand.
#[must_use]
pub fn complement_base(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'U' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' => b'N',
        // IUPAC ambiguity codes
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        _ => b'N',
    }
}

/// Reverse-complement a byte slice.
#[must_use]
pub fn reverse_complement(bytes: &[u8]) -> Vec<u8> {
    bytes.iter().rev().map(|&b| complement_base(b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_basic() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AACG"), b"CGTT");
        assert_eq!(reverse_complement(b"NNGC"), b"GCNN");
    }

    #[test]
    fn uppercases() {
        assert_eq!(Sequence::from_str("acgt").as_bytes(), b"ACGT");
    }
}
