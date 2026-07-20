//! A deterministic synthetic pangenome, so the whole meter runs end-to-end with no genomics
//! toolchain or downloads. It is MHC-like in structure — a shared backbone, SNP/indel bubbles
//! that individuals choose between, and polymorphic repeat insertions drawn from a small
//! consensus library (Alu/LINE/SINE/satellite/ERV) — which is exactly the redundancy the
//! dictionary layer should factor out. Everything is seeded, so runs are reproducible.

use itpp_core::{
    DictClass, DictEntry, Dictionary, Graph, NodeId, RepeatInstance, Sequence, Step, Strand, Walk,
};

/// Knobs for the generator.
#[derive(Debug, Clone)]
pub struct SynthParams {
    pub haplotypes: usize,
    pub backbone_blocks: usize,
    pub block_len: usize,
    /// Probability a backbone block has an alternate allele (SNP bubble).
    pub snp_rate: f64,
    pub insertion_sites: usize,
    pub num_consensus: usize,
    pub consensus_len: usize,
    pub seed: u64,
}

impl Default for SynthParams {
    fn default() -> Self {
        SynthParams {
            haplotypes: 8,
            backbone_blocks: 60,
            block_len: 400,
            snp_rate: 0.35,
            insertion_sites: 12,
            num_consensus: 4,
            consensus_len: 300,
            seed: 42,
        }
    }
}

/// SplitMix64 — tiny, dependency-free, deterministic.
struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n
    }
    fn chance(&mut self, p: f64) -> bool {
        (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64) < p
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[(self.next_u64() & 3) as usize]
    }
    fn random_seq(&mut self, len: usize) -> Vec<u8> {
        (0..len).map(|_| self.base()).collect()
    }
    /// Copy with a few point substitutions (models divergence from a consensus/reference).
    fn mutate(&mut self, seq: &[u8], subs: usize) -> Vec<u8> {
        let mut out = seq.to_vec();
        for _ in 0..subs {
            if out.is_empty() {
                break;
            }
            let pos = self.below(out.len() as u64) as usize;
            out[pos] = self.base();
        }
        out
    }
}

const CLASSES: [DictClass; 5] =
    [DictClass::Alu, DictClass::Line, DictClass::Sine, DictClass::Satellite, DictClass::Erv];

/// Build the synthetic [`Graph`] (segments, edges, walks, dictionary, repeat instances).
#[must_use]
pub fn generate(p: &SynthParams) -> Graph {
    let mut rng = Rng(p.seed);
    let mut g = Graph::new();
    let mut next_id: NodeId = 0;

    // 1. Repeat/ERV consensus library.
    let mut consensi: Vec<Vec<u8>> = Vec::new();
    for j in 0..p.num_consensus {
        let seq = rng.random_seq(p.consensus_len);
        g.dictionary.entries.push(DictEntry {
            id: j as u32,
            class: CLASSES[j % CLASSES.len()],
            name: format!("fam{j}"),
            consensus: Sequence::from_bytes(seq.clone()),
        });
        consensi.push(seq);
    }

    // 2. Backbone blocks + optional alternate alleles.
    let mut backbone_ids: Vec<NodeId> = Vec::with_capacity(p.backbone_blocks);
    let mut alt_of: Vec<Option<NodeId>> = Vec::with_capacity(p.backbone_blocks);
    for _ in 0..p.backbone_blocks {
        let seq = rng.random_seq(p.block_len);
        let id = next_id;
        next_id += 1;
        g.add_segment(id, Sequence::from_bytes(seq.clone()));
        backbone_ids.push(id);
        if rng.chance(p.snp_rate) {
            let alt_seq = rng.mutate(&seq, 3.max(p.block_len / 100));
            let alt_id = next_id;
            next_id += 1;
            g.add_segment(alt_id, Sequence::from_bytes(alt_seq));
            alt_of.push(Some(alt_id));
        } else {
            alt_of.push(None);
        }
    }

    // 3. Polymorphic repeat insertions between blocks.
    // site -> (insertion node id, index of block it follows)
    let mut insertions: Vec<(NodeId, usize)> = Vec::new();
    for _ in 0..p.insertion_sites {
        let after = rng.below(p.backbone_blocks as u64) as usize;
        let fam = rng.below(p.num_consensus.max(1) as u64) as usize;
        let inst_seq = rng.mutate(&consensi[fam], 6.max(p.consensus_len / 40));
        let id = next_id;
        next_id += 1;
        g.add_segment(id, Sequence::from_bytes(inst_seq));
        g.instances.insert(id, RepeatInstance { dict_id: fam as u32, strand: Strand::Forward });
        insertions.push((id, after));
    }

    // 4. Backbone walk (the reference spine).
    let backbone_steps: Vec<Step> =
        backbone_ids.iter().map(|&n| Step { node: n, strand: Strand::Forward }).collect();
    let backbone_name = "CHM13_backbone".to_string();
    add_chain_edges(&mut g, &backbone_steps);
    g.add_walk(Walk { name: backbone_name.clone(), steps: backbone_steps });
    g.backbone = Some(backbone_name);

    // 5. Haplotype walks: choose ref/alt per block, carry a subset of insertions.
    for h in 0..p.haplotypes {
        let mut steps = Vec::new();
        for i in 0..p.backbone_blocks {
            let node = match alt_of[i] {
                Some(alt) if rng.chance(0.5) => alt,
                _ => backbone_ids[i],
            };
            steps.push(Step { node, strand: Strand::Forward });
            for &(ins_id, after) in &insertions {
                if after == i && rng.chance(0.4) {
                    steps.push(Step { node: ins_id, strand: Strand::Forward });
                }
            }
        }
        add_chain_edges(&mut g, &steps);
        g.add_walk(Walk { name: format!("sample{h:02}#1"), steps });
    }

    g
}

fn add_chain_edges(g: &mut Graph, steps: &[Step]) {
    for pair in steps.windows(2) {
        let e = itpp_core::Edge {
            from: pair[0].node,
            from_strand: pair[0].strand,
            to: pair[1].node,
            to_strand: pair[1].strand,
        };
        if !g.edges.contains(&e) {
            g.add_edge(e);
        }
    }
}

/// Convenience alias so callers can hold the full library type.
pub type Consensi = Dictionary;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generates_consistent_graph() {
        let p = SynthParams { haplotypes: 6, backbone_blocks: 20, ..Default::default() };
        let g = generate(&p);
        assert_eq!(g.walks.len(), 7); // backbone + 6 samples
        assert!(g.backbone.is_some());
        assert!(!g.dictionary.entries.is_empty());
        assert!(!g.instances.is_empty());
        // every walk spells without dangling node references
        for w in &g.walks {
            assert!(g.spell(w).is_some(), "walk {} has a dangling node", w.name);
        }
    }

    #[test]
    fn deterministic() {
        let p = SynthParams::default();
        assert_eq!(generate(&p), generate(&p));
    }
}
