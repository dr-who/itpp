//! DNA → protein translation for the "3-mers become proteins" panel.
//! Standard genetic code (NCBI transl_table 1). `*` = stop, `X` = unknown/ambiguous.

/// Translate a single codon (3 uppercase ASCII bases) to a one-letter amino acid.
#[must_use]
pub fn translate_codon(c: &[u8]) -> char {
    if c.len() != 3 {
        return 'X';
    }
    match [c[0].to_ascii_uppercase(), c[1].to_ascii_uppercase(), c[2].to_ascii_uppercase()] {
        [b'T', b'T', b'T'] | [b'T', b'T', b'C'] => 'F',
        [b'T', b'T', b'A'] | [b'T', b'T', b'G'] => 'L',
        [b'C', b'T', _] => 'L',
        [b'A', b'T', b'T'] | [b'A', b'T', b'C'] | [b'A', b'T', b'A'] => 'I',
        [b'A', b'T', b'G'] => 'M',
        [b'G', b'T', _] => 'V',
        [b'T', b'C', _] => 'S',
        [b'C', b'C', _] => 'P',
        [b'A', b'C', _] => 'T',
        [b'G', b'C', _] => 'A',
        [b'T', b'A', b'T'] | [b'T', b'A', b'C'] => 'Y',
        [b'T', b'A', b'A'] | [b'T', b'A', b'G'] => '*',
        [b'C', b'A', b'T'] | [b'C', b'A', b'C'] => 'H',
        [b'C', b'A', b'A'] | [b'C', b'A', b'G'] => 'Q',
        [b'A', b'A', b'T'] | [b'A', b'A', b'C'] => 'N',
        [b'A', b'A', b'A'] | [b'A', b'A', b'G'] => 'K',
        [b'G', b'A', b'T'] | [b'G', b'A', b'C'] => 'D',
        [b'G', b'A', b'A'] | [b'G', b'A', b'G'] => 'E',
        [b'T', b'G', b'T'] | [b'T', b'G', b'C'] => 'C',
        [b'T', b'G', b'A'] => '*',
        [b'T', b'G', b'G'] => 'W',
        [b'C', b'G', _] => 'R',
        [b'A', b'G', b'T'] | [b'A', b'G', b'C'] => 'S',
        [b'A', b'G', b'A'] | [b'A', b'G', b'G'] => 'R',
        [b'G', b'G', _] => 'G',
        _ => 'X',
    }
}

/// Split `dna` into consecutive codons from `frame` (0,1,2) and translate to a protein string.
#[must_use]
pub fn translate(dna: &[u8], frame: usize) -> String {
    let mut out = String::new();
    let mut i = frame;
    while i + 3 <= dna.len() {
        out.push(translate_codon(&dna[i..i + 3]));
        i += 3;
    }
    out
}

/// The list of 3-mers (as strings) for `dna` from `frame`, for the codon panel.
#[must_use]
pub fn codons(dna: &[u8], frame: usize) -> Vec<String> {
    let mut out = Vec::new();
    let mut i = frame;
    while i + 3 <= dna.len() {
        out.push(String::from_utf8_lossy(&dna[i..i + 3]).into_owned());
        i += 3;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_codons() {
        assert_eq!(translate_codon(b"ATG"), 'M'); // start
        assert_eq!(translate_codon(b"TAA"), '*'); // stop
        assert_eq!(translate_codon(b"TGG"), 'W');
        assert_eq!(translate_codon(b"GCA"), 'A');
        assert_eq!(translate_codon(b"aaa"), 'K'); // lowercase ok
    }

    #[test]
    fn translate_orf() {
        // ATG GCA TAA -> M A *
        assert_eq!(translate(b"ATGGCATAA", 0), "MA*");
        // frame shift
        assert_eq!(translate(b"AATGGCATAA", 1), "MA*");
    }

    #[test]
    fn codon_split() {
        assert_eq!(codons(b"ATGGCATAA", 0), vec!["ATG", "GCA", "TAA"]);
        assert_eq!(codons(b"ATGG", 0), vec!["ATG"]); // trailing partial dropped
    }
}
