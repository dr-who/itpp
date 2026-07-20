//! A minimal FASTA reader (enough to pull a backbone or haplotype sequences).

use itpp_core::Sequence;

/// Parse FASTA text into `(name, sequence)` records. The name is the id up to the first
/// whitespace on the header line.
#[must_use]
pub fn parse(text: &str) -> Vec<(String, Sequence)> {
    let mut records = Vec::new();
    let mut name: Option<String> = None;
    let mut seq: Vec<u8> = Vec::new();
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(n) = name.take() {
                records.push((n, Sequence::from_bytes(std::mem::take(&mut seq))));
            }
            name = Some(header.split_whitespace().next().unwrap_or("").to_string());
        } else {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if let Some(n) = name.take() {
        records.push((n, Sequence::from_bytes(seq)));
    }
    records
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_two_records() {
        let text = ">chr6 MHC\nACGT\nACGT\n>alt\nTTTT\n";
        let r = parse(text);
        assert_eq!(r.len(), 2);
        assert_eq!(r[0].0, "chr6");
        assert_eq!(r[0].1.as_bytes(), b"ACGTACGT");
        assert_eq!(r[1].1.as_bytes(), b"TTTT");
    }
}
