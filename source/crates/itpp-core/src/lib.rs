//! Core types for the ITPP two-part-code pangenome.
//!
//! A pangenome is a [`Graph`] of sequence [`Segment`]s joined by [`Edge`]s. A haplotype
//! (an individual's genome, or one of its two copies) is a [`Walk`]: an ordered list of
//! oriented segments. The whole point of the representation is that segments shared by many
//! haplotypes are stored once, and repeat families collapse to a [`Dictionary`] consensus
//! plus per-copy edit scripts — so the *total* description length shrinks as duplication is
//! factored out (see `docs/design.md`).
//!
//! These types are deliberately plain data (no async, no external deps) so the later
//! Python (PyO3) and C-ABI bindings are cheap.

use std::collections::BTreeMap;

pub mod edit;
pub mod seq;

pub use edit::{EditOp, EditScript};
pub use seq::{reverse_complement, Sequence};

/// Segment identifier. Stable within a container.
pub type NodeId = u64;

/// Orientation of a segment as traversed by a walk or referenced by an edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    #[must_use]
    pub fn flip(self) -> Self {
        match self {
            Strand::Forward => Strand::Reverse,
            Strand::Reverse => Strand::Forward,
        }
    }

    /// GFA orientation sign.
    #[must_use]
    pub fn as_sign(self) -> char {
        match self {
            Strand::Forward => '+',
            Strand::Reverse => '-',
        }
    }

    #[must_use]
    pub fn from_sign(c: char) -> Option<Self> {
        match c {
            '+' => Some(Strand::Forward),
            '-' => Some(Strand::Reverse),
            _ => None,
        }
    }
}

/// A single sequence node in the graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Segment {
    pub id: NodeId,
    pub seq: Sequence,
}

/// A directed adjacency between two oriented segment ends (GFA `L` line semantics).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Edge {
    pub from: NodeId,
    pub from_strand: Strand,
    pub to: NodeId,
    pub to_strand: Strand,
}

/// One oriented step of a walk.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Step {
    pub node: NodeId,
    pub strand: Strand,
}

/// A haplotype: an ordered traversal of the graph. This is the per-individual "path through
/// the graph" — encoding it (given the shared graph) is the residual side of the two-part code.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Walk {
    pub name: String,
    pub steps: Vec<Step>,
}

/// Biological class of a dictionary entry (the repeat/ERV library).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DictClass {
    Alu,
    Line,
    Sine,
    Satellite,
    Erv,
    Other,
}

impl DictClass {
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            DictClass::Alu => "Alu",
            DictClass::Line => "LINE",
            DictClass::Sine => "SINE",
            DictClass::Satellite => "Satellite",
            DictClass::Erv => "ERV",
            DictClass::Other => "Other",
        }
    }
}

/// A consensus sequence for a repeat/ERV family. Instances elsewhere reference it as
/// `(id, orientation, edit-script)` instead of storing the bases again.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DictEntry {
    pub id: u32,
    pub class: DictClass,
    pub name: String,
    pub consensus: Sequence,
}

/// The repeat/ERV dictionary.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Dictionary {
    pub entries: Vec<DictEntry>,
}

impl Dictionary {
    #[must_use]
    pub fn get(&self, id: u32) -> Option<&DictEntry> {
        self.entries.iter().find(|e| e.id == id)
    }
}

/// A segment that is an *instance* of a dictionary family: instead of storing its bases, we
/// store `(dict_id, strand)` and recover the sequence by diffing against the consensus. The
/// per-copy edit script is derived at measure time (see `itpp-codec`), so this stays plain data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RepeatInstance {
    pub dict_id: u32,
    pub strand: Strand,
}

/// The pangenome graph plus the walks through it.
///
/// `backbone` names the walk designated as the reference spine (its sequence anchors the
/// coordinate space). Segments are keyed by id for stable iteration and delta coding.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Graph {
    pub segments: BTreeMap<NodeId, Segment>,
    pub edges: Vec<Edge>,
    pub walks: Vec<Walk>,
    pub dictionary: Dictionary,
    /// Segments known to be instances of a dictionary family (the repeat/ERV factoring).
    pub instances: BTreeMap<NodeId, RepeatInstance>,
    /// Non-human spans flagged by a contamination screen, kept out of the human bits/char.
    pub contamination: Vec<Contam>,
    pub backbone: Option<String>,
}

/// A contamination call: a span of a sample's walk attributed to a non-human taxon.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contam {
    pub sample: String,
    pub node: NodeId,
    pub taxon: String,
}

impl Graph {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_segment(&mut self, id: NodeId, seq: Sequence) {
        self.segments.insert(id, Segment { id, seq });
    }

    pub fn add_edge(&mut self, edge: Edge) {
        self.edges.push(edge);
    }

    pub fn add_walk(&mut self, walk: Walk) {
        self.walks.push(walk);
    }

    /// Total bases across all *distinct* segments (the shared sequence content).
    #[must_use]
    pub fn segment_bases(&self) -> usize {
        self.segments.values().map(|s| s.seq.len()).sum()
    }

    /// Total bases across all walks (the denominator of cohort bits/char — every haplotype's
    /// full length, counting shared segments once per traversal).
    #[must_use]
    pub fn walk_bases(&self) -> usize {
        self.walks.iter().map(|w| self.walk_len(w)).sum()
    }

    fn walk_len(&self, w: &Walk) -> usize {
        w.steps
            .iter()
            .map(|s| self.segments.get(&s.node).map_or(0, |seg| seg.seq.len()))
            .sum()
    }

    /// Reconstruct the full nucleotide sequence spelled by a walk. This is the inverse used
    /// by `itpp verify` to prove losslessness.
    ///
    /// Returns `None` if the walk references a segment not in the graph.
    #[must_use]
    pub fn spell(&self, walk: &Walk) -> Option<Sequence> {
        let mut out = Vec::new();
        for step in &walk.steps {
            let seg = self.segments.get(&step.node)?;
            match step.strand {
                Strand::Forward => out.extend_from_slice(seg.seq.as_bytes()),
                Strand::Reverse => out.extend(reverse_complement(seg.seq.as_bytes())),
            }
        }
        Some(Sequence::from_bytes(out))
    }

    #[must_use]
    pub fn walk_by_name(&self, name: &str) -> Option<&Walk> {
        self.walks.iter().find(|w| w.name == name)
    }

    /// The walk designated as the backbone spine, if any.
    #[must_use]
    pub fn backbone_walk(&self) -> Option<&Walk> {
        self.backbone.as_deref().and_then(|n| self.walk_by_name(n))
    }

    /// Set of segment ids traversed by the backbone walk (the coordinate anchor).
    #[must_use]
    pub fn backbone_nodes(&self) -> std::collections::BTreeSet<NodeId> {
        self.backbone_walk()
            .map(|w| w.steps.iter().map(|s| s.node).collect())
            .unwrap_or_default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spell_forward_and_reverse() {
        let mut g = Graph::new();
        g.add_segment(1, Sequence::from_str("ACGT"));
        g.add_segment(2, Sequence::from_str("TTAA"));
        let w = Walk {
            name: "h1".into(),
            steps: vec![
                Step { node: 1, strand: Strand::Forward },
                Step { node: 2, strand: Strand::Reverse },
            ],
        };
        // ACGT + revcomp(TTAA)=TTAA
        assert_eq!(g.spell(&w).unwrap().as_bytes(), b"ACGTTTAA");
    }

    #[test]
    fn bases_accounting() {
        let mut g = Graph::new();
        g.add_segment(1, Sequence::from_str("ACGTACGT")); // 8
        g.add_segment(2, Sequence::from_str("GG")); // 2
        g.add_walk(Walk {
            name: "a".into(),
            steps: vec![
                Step { node: 1, strand: Strand::Forward },
                Step { node: 2, strand: Strand::Forward },
            ],
        });
        g.add_walk(Walk {
            name: "b".into(),
            steps: vec![Step { node: 1, strand: Strand::Forward }],
        });
        assert_eq!(g.segment_bases(), 10); // distinct content
        assert_eq!(g.walk_bases(), 18); // 10 + 8
    }
}
