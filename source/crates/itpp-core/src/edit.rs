//! Edit scripts: how a repeat *instance* or a private *delta* differs from a reference
//! (a dictionary consensus, or a graph node). Storing `(reference_id, edit-script)` instead
//! of raw bases is the factoring that removes duplicated components from the total.

/// A single edit operation transforming a reference sequence into a target.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EditOp {
    /// Copy `len` bases straight from the reference.
    Copy { len: u32 },
    /// Substitute the next `bases.len()` reference bases with `bases`.
    Sub { bases: Vec<u8> },
    /// Insert `bases` (not present in the reference).
    Ins { bases: Vec<u8> },
    /// Delete `len` reference bases.
    Del { len: u32 },
}

/// An ordered list of [`EditOp`]s.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct EditScript(pub Vec<EditOp>);

impl EditScript {
    #[must_use]
    pub fn apply(&self, reference: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        let mut r = 0usize;
        for op in &self.0 {
            match op {
                EditOp::Copy { len } => {
                    let end = (r + *len as usize).min(reference.len());
                    out.extend_from_slice(&reference[r..end]);
                    r = end;
                }
                EditOp::Sub { bases } => {
                    out.extend_from_slice(bases);
                    r = (r + bases.len()).min(reference.len());
                }
                EditOp::Ins { bases } => out.extend_from_slice(bases),
                EditOp::Del { len } => r = (r + *len as usize).min(reference.len()),
            }
        }
        out
    }

    /// Number of literal bases carried in the script (the part that is genuinely novel and
    /// must be entropy-coded; `Copy`/`Del` are near-free positional metadata).
    #[must_use]
    pub fn literal_bases(&self) -> usize {
        self.0
            .iter()
            .map(|op| match op {
                EditOp::Sub { bases } | EditOp::Ins { bases } => bases.len(),
                _ => 0,
            })
            .sum()
    }
}

/// A trivial anchored diff: emit `Copy` over the common prefix, one `Sub`/`Ins`/`Del` block
/// for the middle, and `Copy` over the common suffix. Good enough to model SNP/indel-level
/// divergence of a repeat instance from its consensus for v0; a real aligner replaces this.
#[must_use]
pub fn anchored_diff(reference: &[u8], target: &[u8]) -> EditScript {
    let mut ops = Vec::new();
    let max_pre = reference.len().min(target.len());
    let mut pre = 0;
    while pre < max_pre && reference[pre] == target[pre] {
        pre += 1;
    }
    let mut suf = 0;
    while suf < (max_pre - pre)
        && reference[reference.len() - 1 - suf] == target[target.len() - 1 - suf]
    {
        suf += 1;
    }
    if pre > 0 {
        ops.push(EditOp::Copy { len: pre as u32 });
    }
    let ref_mid = &reference[pre..reference.len() - suf];
    let tgt_mid = &target[pre..target.len() - suf];
    if !ref_mid.is_empty() && !tgt_mid.is_empty() {
        // Substitute the overlapping middle, then insert/delete the remainder.
        let common = ref_mid.len().min(tgt_mid.len());
        ops.push(EditOp::Sub { bases: tgt_mid[..common].to_vec() });
        if tgt_mid.len() > common {
            ops.push(EditOp::Ins { bases: tgt_mid[common..].to_vec() });
        } else if ref_mid.len() > common {
            ops.push(EditOp::Del { len: (ref_mid.len() - common) as u32 });
        }
    } else if !tgt_mid.is_empty() {
        ops.push(EditOp::Ins { bases: tgt_mid.to_vec() });
    } else if !ref_mid.is_empty() {
        ops.push(EditOp::Del { len: ref_mid.len() as u32 });
    }
    if suf > 0 {
        ops.push(EditOp::Copy { len: suf as u32 });
    }
    EditScript(ops)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn apply_roundtrip_snp() {
        let reference = b"ACGTACGTAC";
        let target = b"ACGTTCGTAC"; // single sub at pos 4
        let script = anchored_diff(reference, target);
        assert_eq!(script.apply(reference), target);
    }

    #[test]
    fn apply_roundtrip_ins_del() {
        let reference = b"ACGTACGT";
        for target in [b"ACGTGGACGT".to_vec(), b"ACGT".to_vec(), b"ACGTACGTAAA".to_vec()] {
            let script = anchored_diff(reference, &target);
            assert_eq!(script.apply(reference), target, "target={target:?}");
        }
    }
}
