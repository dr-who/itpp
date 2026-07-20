//! Adaptive order-k context-model coder (PPM-style) for nucleotide byte strings.
//!
//! The alphabet is learned from the data (distinct bytes) so the coder is fully general and
//! lossless — it round-trips `N`, IUPAC codes, anything. For DNA the alphabet is ~ACGT(N),
//! and the order-k model captures the local dependence that pushes real DNA below 2 bits/base.
//!
//! Contexts (the previous `k` symbols) are held in a `HashMap`, allocated lazily, so memory
//! scales with the number of *observed* contexts rather than `alphabet^k`.

use crate::range::{
    bits_for, bittree_decode, bittree_encode, RangeDecoder, RangeEncoder, PROB_INIT,
};
use std::collections::HashMap;

/// A self-describing compressed blob: header (order, alphabet, length) + range-coded payload.
///
/// Layout: `k:u8 | alphabet_len:u8 | alphabet bytes | orig_len:u64(LE) | payload`.
pub struct SeqCoder {
    order: usize,
}

impl SeqCoder {
    #[must_use]
    pub fn new(order: usize) -> Self {
        // Keep the packed context key within 56 bits: order * bits_per_symbol <= 56.
        SeqCoder { order }
    }

    /// Compress `data`. The output is standalone (carries its own header).
    #[must_use]
    pub fn encode(&self, data: &[u8]) -> Vec<u8> {
        // Build the alphabet.
        let mut seen = [false; 256];
        for &b in data {
            seen[b as usize] = true;
        }
        let alphabet: Vec<u8> = (0u16..256).filter(|&b| seen[b as usize]).map(|b| b as u8).collect();
        let mut index = [0u32; 256];
        for (i, &b) in alphabet.iter().enumerate() {
            index[b as usize] = i as u32;
        }
        let nbits = bits_for(alphabet.len().max(1));
        let bps = nbits; // bits per symbol used to pack the context key
        let order = effective_order(self.order, bps);
        let ctx_mask: u64 = if order == 0 { 0 } else { (1u64 << (bps * order as u32)) - 1 };

        let mut header = Vec::new();
        header.push(order as u8);
        header.extend_from_slice(&(alphabet.len() as u16).to_le_bytes());
        header.extend_from_slice(&alphabet);
        header.extend_from_slice(&(data.len() as u64).to_le_bytes());

        if alphabet.len() <= 1 {
            // All bytes identical (or empty): payload is empty, header suffices.
            return header;
        }

        let mut enc = RangeEncoder::new();
        let mut contexts: HashMap<u64, Vec<u16>> = HashMap::new();
        let mut ctx: u64 = 0;
        for &b in data {
            let sym = index[b as usize];
            let probs = contexts.entry(ctx).or_insert_with(|| vec![PROB_INIT; 1 << nbits]);
            bittree_encode(&mut enc, probs, nbits, sym);
            ctx = ((ctx << bps) | u64::from(sym)) & ctx_mask;
        }
        header.extend_from_slice(&enc.finish());
        header
    }

    /// Decompress a blob produced by [`SeqCoder::encode`]. `order` in `self` is ignored — the
    /// blob is self-describing.
    #[must_use]
    pub fn decode(blob: &[u8]) -> Vec<u8> {
        let order = blob[0] as usize;
        let alpha_len = u16::from_le_bytes([blob[1], blob[2]]) as usize;
        let mut pos = 3;
        let alphabet = &blob[pos..pos + alpha_len];
        pos += alpha_len;
        let orig_len =
            u64::from_le_bytes(blob[pos..pos + 8].try_into().expect("len field")) as usize;
        pos += 8;

        if alpha_len <= 1 {
            let fill = alphabet.first().copied().unwrap_or(b'N');
            return vec![fill; orig_len];
        }

        let nbits = bits_for(alpha_len);
        let bps = nbits;
        let ctx_mask: u64 = if order == 0 { 0 } else { (1u64 << (bps * order as u32)) - 1 };

        let mut dec = RangeDecoder::new(&blob[pos..]);
        let mut contexts: HashMap<u64, Vec<u16>> = HashMap::new();
        let mut ctx: u64 = 0;
        let mut out = Vec::with_capacity(orig_len);
        for _ in 0..orig_len {
            let probs = contexts.entry(ctx).or_insert_with(|| vec![PROB_INIT; 1 << nbits]);
            let sym = bittree_decode(&mut dec, probs, nbits);
            out.push(alphabet[sym as usize]);
            ctx = ((ctx << bps) | u64::from(sym)) & ctx_mask;
        }
        out
    }

    /// Achieved cost in bits for `data` under this coder (i.e. `encode(data).len() * 8`).
    #[must_use]
    pub fn cost_bits(&self, data: &[u8]) -> u64 {
        self.encode(data).len() as u64 * 8
    }
}

/// Clamp the requested order so the packed context key fits in 56 bits.
fn effective_order(order: usize, bits_per_symbol: u32) -> usize {
    if bits_per_symbol == 0 {
        return 0;
    }
    let max_order = (56 / bits_per_symbol) as usize;
    order.min(max_order)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip(order: usize, data: &[u8]) {
        let coder = SeqCoder::new(order);
        let blob = coder.encode(data);
        assert_eq!(SeqCoder::decode(&blob), data, "order={order} len={}", data.len());
    }

    #[test]
    fn roundtrip_dna() {
        let mut data = Vec::new();
        let bases = b"ACGT";
        for i in 0..5000u32 {
            data.push(bases[((i.wrapping_mul(2654435761)) >> 29) as usize % 4]);
        }
        for k in [0, 2, 4, 8, 12] {
            roundtrip(k, &data);
        }
    }

    #[test]
    fn roundtrip_edge_cases() {
        roundtrip(4, b"");
        roundtrip(4, b"A");
        roundtrip(4, b"AAAAAAAA");
        roundtrip(4, b"ACGTNNNNRYSWacgt"); // mixed alphabet incl. N/IUPAC/lowercase
    }

    #[test]
    fn repetitive_compresses_below_two_bits() {
        // A tandem repeat should code well under 2 bits/base with context.
        let unit = b"ACGTACGTTT";
        let mut data = Vec::new();
        for _ in 0..1000 {
            data.extend_from_slice(unit);
        }
        let coder = SeqCoder::new(8);
        let bits = coder.cost_bits(&data);
        let bpc = bits as f64 / data.len() as f64;
        assert!(bpc < 1.0, "expected < 1 bit/base on tandem repeat, got {bpc}");
    }
}
