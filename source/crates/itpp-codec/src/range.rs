//! A binary range coder (LZMA-style) with adaptive bit probabilities, plus a bit-tree for
//! multi-bit symbols. This is the engine behind every "achieved bits" number the ledger
//! reports: we encode real data and measure the real output length.
//!
//! The design is the well-known LZMA rc: 11-bit probabilities, 5-bit adaptation, `1<<24`
//! renormalization threshold, with the canonical carry-propagating `shift_low`.

const TOP: u32 = 1 << 24;
pub const PROB_BITS: u32 = 11;
pub const PROB_TOTAL: u32 = 1 << PROB_BITS; // 2048
pub const PROB_INIT: u16 = (PROB_TOTAL / 2) as u16; // 1024
const MOVE_BITS: u32 = 5;

/// Encoder: feed bits, get a byte stream.
pub struct RangeEncoder {
    low: u64,
    range: u32,
    cache: u8,
    cache_size: u64,
    out: Vec<u8>,
}

impl Default for RangeEncoder {
    fn default() -> Self {
        Self::new()
    }
}

impl RangeEncoder {
    #[must_use]
    pub fn new() -> Self {
        RangeEncoder { low: 0, range: 0xFFFF_FFFF, cache: 0, cache_size: 1, out: Vec::new() }
    }

    fn shift_low(&mut self) {
        if (self.low >> 32) != 0 || self.low < 0xFF00_0000 {
            let mut temp = self.cache;
            loop {
                self.out.push((temp as u64).wrapping_add(self.low >> 32) as u8);
                temp = 0xFF;
                self.cache_size -= 1;
                if self.cache_size == 0 {
                    break;
                }
            }
            self.cache = (self.low >> 24) as u8;
        }
        self.cache_size += 1;
        self.low = (self.low << 8) & 0xFFFF_FFFF;
    }

    /// Encode one bit under an adaptive probability `prob` (probability that bit == 0).
    pub fn encode_bit(&mut self, prob: &mut u16, bit: u32) {
        let bound = (self.range >> PROB_BITS) * u32::from(*prob);
        if bit == 0 {
            self.range = bound;
            *prob += ((PROB_TOTAL - u32::from(*prob)) >> MOVE_BITS) as u16;
        } else {
            self.low += u64::from(bound);
            self.range -= bound;
            *prob -= *prob >> MOVE_BITS;
        }
        while self.range < TOP {
            self.range <<= 8;
            self.shift_low();
        }
    }

    /// Finish the stream. Must be called exactly once.
    #[must_use]
    pub fn finish(mut self) -> Vec<u8> {
        for _ in 0..5 {
            self.shift_low();
        }
        self.out
    }
}

/// Decoder: mirror of [`RangeEncoder`].
pub struct RangeDecoder<'a> {
    code: u32,
    range: u32,
    data: &'a [u8],
    pos: usize,
}

impl<'a> RangeDecoder<'a> {
    #[must_use]
    pub fn new(data: &'a [u8]) -> Self {
        let mut d = RangeDecoder { code: 0, range: 0xFFFF_FFFF, data, pos: 0 };
        // First emitted byte is always 0 (initial cache); skip it, then load 4 code bytes.
        d.next_byte();
        for _ in 0..4 {
            d.code = (d.code << 8) | u32::from(d.next_byte());
        }
        d
    }

    fn next_byte(&mut self) -> u8 {
        let b = self.data.get(self.pos).copied().unwrap_or(0);
        self.pos += 1;
        b
    }

    pub fn decode_bit(&mut self, prob: &mut u16) -> u32 {
        let bound = (self.range >> PROB_BITS) * u32::from(*prob);
        let bit;
        if self.code < bound {
            self.range = bound;
            *prob += ((PROB_TOTAL - u32::from(*prob)) >> MOVE_BITS) as u16;
            bit = 0;
        } else {
            self.code -= bound;
            self.range -= bound;
            *prob -= *prob >> MOVE_BITS;
            bit = 1;
        }
        while self.range < TOP {
            self.range <<= 8;
            self.code = (self.code << 8) | u32::from(self.next_byte());
        }
        bit
    }
}

/// Number of bits needed to index `n` distinct symbols (min 1).
#[must_use]
pub fn bits_for(n: usize) -> u32 {
    if n <= 1 {
        1
    } else {
        usize::BITS - (n - 1).leading_zeros()
    }
}

/// Encode `symbol` (`< 1<<nbits`) MSB-first through a bit-tree of `1<<nbits` probabilities.
pub fn bittree_encode(enc: &mut RangeEncoder, probs: &mut [u16], nbits: u32, symbol: u32) {
    let mut m: u32 = 1;
    for i in (0..nbits).rev() {
        let bit = (symbol >> i) & 1;
        enc.encode_bit(&mut probs[m as usize], bit);
        m = (m << 1) | bit;
    }
}

/// Inverse of [`bittree_encode`].
pub fn bittree_decode(dec: &mut RangeDecoder, probs: &mut [u16], nbits: u32) -> u32 {
    let mut m: u32 = 1;
    for _ in 0..nbits {
        let bit = dec.decode_bit(&mut probs[m as usize]);
        m = (m << 1) | bit;
    }
    m - (1 << nbits)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_bit_roundtrip() {
        let bits = [0u32, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0];
        let mut enc = RangeEncoder::new();
        let mut p = PROB_INIT;
        for &b in &bits {
            enc.encode_bit(&mut p, b);
        }
        let data = enc.finish();
        let mut dec = RangeDecoder::new(&data);
        let mut p2 = PROB_INIT;
        for &b in &bits {
            assert_eq!(dec.decode_bit(&mut p2), b);
        }
    }

    #[test]
    fn bittree_roundtrip() {
        let syms: Vec<u32> = (0..500).map(|i| (i * 7 + 3) % 4).collect();
        let nbits = bits_for(4);
        let mut enc = RangeEncoder::new();
        let mut probs = vec![PROB_INIT; 1 << nbits];
        for &s in &syms {
            bittree_encode(&mut enc, &mut probs, nbits, s);
        }
        let data = enc.finish();
        let mut dec = RangeDecoder::new(&data);
        let mut probs2 = vec![PROB_INIT; 1 << nbits];
        for &s in &syms {
            assert_eq!(bittree_decode(&mut dec, &mut probs2, nbits), s);
        }
    }

    #[test]
    fn bits_for_values() {
        assert_eq!(bits_for(1), 1);
        assert_eq!(bits_for(2), 1);
        assert_eq!(bits_for(4), 2);
        assert_eq!(bits_for(5), 3);
        assert_eq!(bits_for(16), 4);
    }
}
