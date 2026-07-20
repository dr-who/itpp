//! Entropy coders and the description-length ledger for ITPP.
//!
//! - [`range`] — an LZMA-style binary range coder (the engine that produces real bits).
//! - [`model::SeqCoder`] — an adaptive order-k context model over it (PPM-style), the main
//!   sequence coder; "achieved bits" = real output length.
//! - [`baseline`] — 2-bit and order-0 Shannon reference points to beat.
//! - [`ledger`] — per-section accounting → total bits → **bits/char**.
//!
//! Kept dependency-free and plain-data so Python/C bindings stay cheap. External coders
//! (`xz`/`zstd`/`bzip2`) are invoked from the CLI layer and folded into the same [`ledger`].

pub mod baseline;
pub mod ledger;
pub mod model;
pub mod range;

pub use ledger::{CoderResult, Ledger, SectionLedger};
pub use model::SeqCoder;

/// Convenience: run the internal coders on `data` and return their achieved bits, labelled.
/// The CLI adds external cross-checks (xz/zstd/bzip2) to the same section.
#[must_use]
pub fn measure_internal(data: &[u8], orders: &[usize]) -> Vec<CoderResult> {
    let mut out = Vec::new();
    out.push(CoderResult { coder: "2bit".into(), bits: baseline::twobit_total_bits(data) });
    for &k in orders {
        let bits = SeqCoder::new(k).cost_bits(data);
        out.push(CoderResult { coder: format!("order-{k}"), bits });
    }
    out
}
