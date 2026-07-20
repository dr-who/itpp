//! Reference points the graph model must beat. These are *not* the reported result — they
//! calibrate it. See `docs/design.md`.

/// Order-0 Shannon entropy in bits per symbol (redundancy-blind). Cheap sanity baseline only.
#[must_use]
pub fn shannon_order0_bits(data: &[u8]) -> f64 {
    if data.is_empty() {
        return 0.0;
    }
    let mut hist = [0u64; 256];
    for &b in data {
        hist[b as usize] += 1;
    }
    let n = data.len() as f64;
    let mut bits = 0.0;
    for &c in &hist {
        if c == 0 {
            continue;
        }
        let p = c as f64 / n;
        bits -= p * p.log2();
    }
    bits
}

/// Total bits to store `data` as naive 2-bit-packed ACGT, charging 8 bits for every base
/// outside `{A,C,G,T}` (an escape). This is the "2.0 bits/base" reference for pure ACGT.
#[must_use]
pub fn twobit_total_bits(data: &[u8]) -> u64 {
    let mut bits: u64 = 0;
    for &b in data {
        match b {
            b'A' | b'C' | b'G' | b'T' => bits += 2,
            _ => bits += 8,
        }
    }
    bits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shannon_uniform_dna_is_two_bits() {
        let data = b"ACGTACGTACGTACGT";
        let h = shannon_order0_bits(data);
        assert!((h - 2.0).abs() < 1e-9, "got {h}");
    }

    #[test]
    fn twobit_counts_escapes() {
        assert_eq!(twobit_total_bits(b"ACGT"), 8);
        assert_eq!(twobit_total_bits(b"ACGTN"), 8 + 8);
    }
}
