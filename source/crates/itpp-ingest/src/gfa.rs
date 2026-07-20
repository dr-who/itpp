//! Minimal GFA v1 reader/writer: enough of `S` / `L` / `P` / `W` to ingest an HPRC
//! minigraph-cactus graph and to round-trip our own synthetic graphs.
//!
//! Segment names are interned to [`NodeId`]s in first-appearance order, so the resulting
//! [`Graph`] is self-consistent regardless of the original naming scheme.

use itpp_core::{Edge, Graph, NodeId, Sequence, Step, Strand, Walk};
use std::collections::HashMap;

#[derive(Default)]
struct Interner {
    map: HashMap<String, NodeId>,
    next: NodeId,
}

impl Interner {
    fn intern(&mut self, name: &str) -> NodeId {
        if let Some(&id) = self.map.get(name) {
            return id;
        }
        let id = self.next;
        self.next += 1;
        self.map.insert(name.to_string(), id);
        id
    }
}

/// Parse GFA text into a [`Graph`]. The backbone walk is chosen by name heuristic (see
/// [`pick_backbone`]) or defaults to the first walk.
#[must_use]
pub fn parse(text: &str) -> Graph {
    let mut g = Graph::new();
    let mut interner = Interner::default();
    // rGFA (minigraph) reference tags: node -> (rank, offset, contig). rank 0 = the reference
    // spine. When a graph has no P/W lines (as minigraph output), we synthesize the backbone
    // walk from the rank-0 segments in offset order.
    let mut rgfa: Vec<(NodeId, i64, String)> = Vec::new();
    let mut has_walk = false;

    for line in text.lines() {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "S" if f.len() >= 3 => {
                let id = interner.intern(f[1]);
                g.add_segment(id, Sequence::from_str(f[2]));
                if let Some((rank, off, sn)) = parse_rgfa_tags(&f[3..]) {
                    if rank == 0 {
                        rgfa.push((id, off, sn));
                    }
                }
            }
            "L" if f.len() >= 5 => {
                if let (Some(fs), Some(ts)) =
                    (Strand::from_sign(sign(f[2])), Strand::from_sign(sign(f[4])))
                {
                    let from = interner.intern(f[1]);
                    let to = interner.intern(f[3]);
                    g.add_edge(Edge { from, from_strand: fs, to, to_strand: ts });
                }
            }
            "P" if f.len() >= 3 => {
                let steps = parse_path_steps(f[2], &mut interner);
                g.add_walk(Walk { name: f[1].to_string(), steps });
                has_walk = true;
            }
            "W" if f.len() >= 7 => {
                // W  sample  hap  seqid  start  end  walk
                let name = format!("{}#{}#{}", f[1], f[2], f[3]);
                let steps = parse_walk_steps(f[6], &mut interner);
                g.add_walk(Walk { name, steps });
                has_walk = true;
            }
            _ => {}
        }
    }

    // minigraph rGFA has no P/W lines: build the backbone spine from rank-0 segments.
    if !has_walk && !rgfa.is_empty() {
        rgfa.sort_by_key(|(_, off, _)| *off);
        let name = rgfa.first().map(|(_, _, sn)| sn.clone()).unwrap_or_else(|| "reference".into());
        let steps: Vec<Step> =
            rgfa.iter().map(|(id, _, _)| Step { node: *id, strand: Strand::Forward }).collect();
        g.add_walk(Walk { name: format!("{name}#0#ref"), steps });
    }

    g.backbone = pick_backbone(&g);
    g
}

/// Extract `(SR rank, SO offset, SN contig)` from rGFA optional tag fields, if present.
fn parse_rgfa_tags(tags: &[&str]) -> Option<(i32, i64, String)> {
    let mut sr = None;
    let mut so = None;
    let mut sn = String::new();
    for t in tags {
        if let Some(v) = t.strip_prefix("SR:i:") {
            sr = v.parse().ok();
        } else if let Some(v) = t.strip_prefix("SO:i:") {
            so = v.parse().ok();
        } else if let Some(v) = t.strip_prefix("SN:Z:") {
            sn = v.to_string();
        }
    }
    Some((sr?, so?, sn))
}

fn sign(s: &str) -> char {
    s.chars().next().unwrap_or('+')
}

/// `P` segment list: `12+,13-,14+`.
fn parse_path_steps(list: &str, interner: &mut Interner) -> Vec<Step> {
    list.split(',')
        .filter(|t| !t.is_empty())
        .filter_map(|t| {
            let (name, sign) = t.split_at(t.len().saturating_sub(1));
            let strand = Strand::from_sign(sign.chars().next()?)?;
            Some(Step { node: interner.intern(name), strand })
        })
        .collect()
}

/// `W` walk string: `>12>13<14`.
fn parse_walk_steps(walk: &str, interner: &mut Interner) -> Vec<Step> {
    let mut steps = Vec::new();
    let mut chars = walk.char_indices().peekable();
    while let Some((_, c)) = chars.next() {
        let strand = match c {
            '>' => Strand::Forward,
            '<' => Strand::Reverse,
            _ => continue,
        };
        let mut name = String::new();
        while let Some(&(_, nc)) = chars.peek() {
            if nc == '>' || nc == '<' {
                break;
            }
            name.push(nc);
            chars.next();
        }
        if !name.is_empty() {
            steps.push(Step { node: interner.intern(&name), strand });
        }
    }
    steps
}

/// Heuristic: pick a reference-spine walk, preferring GRCh38 (so genomic coordinates match
/// annotation databases like dbSNP/ClinVar), then CHM13/T2T, then any reference-ish name,
/// else the first walk. Hints are tried in priority order.
#[must_use]
pub fn pick_backbone(g: &Graph) -> Option<String> {
    const HINTS: [&str; 6] = ["grch38", "chm13", "t2t", "reference", "backbone", "_ref"];
    for hint in HINTS {
        if let Some(w) = g.walks.iter().find(|w| w.name.to_ascii_lowercase().contains(hint)) {
            return Some(w.name.clone());
        }
    }
    g.walks.first().map(|w| w.name.clone())
}

/// Serialize a [`Graph`] back to GFA v1 (segments, links, and paths as `P` lines). Segment
/// names are the interned numeric ids, so a re-import is identical.
#[must_use]
pub fn write(g: &Graph) -> String {
    let mut out = String::from("H\tVN:Z:1.0\n");
    for seg in g.segments.values() {
        out.push_str(&format!(
            "S\t{}\t{}\n",
            seg.id,
            std::str::from_utf8(seg.seq.as_bytes()).unwrap_or("")
        ));
    }
    for e in &g.edges {
        out.push_str(&format!(
            "L\t{}\t{}\t{}\t{}\t0M\n",
            e.from,
            e.from_strand.as_sign(),
            e.to,
            e.to_strand.as_sign()
        ));
    }
    for w in &g.walks {
        let list: Vec<String> =
            w.steps.iter().map(|s| format!("{}{}", s.node, s.strand.as_sign())).collect();
        out.push_str(&format!("P\t{}\t{}\t*\n", w.name, list.join(",")));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_gfa() {
        let text = "H\tVN:Z:1.0\n\
                    S\t1\tACGTACGT\n\
                    S\t2\tTTTT\n\
                    S\t3\tGGGGCCCC\n\
                    L\t1\t+\t2\t+\t0M\n\
                    L\t2\t+\t3\t+\t0M\n\
                    P\tCHM13_backbone\t1+,2+,3+\t*\n\
                    P\tsampleA\t1+,3+\t*\n";
        let g = parse(text);
        assert_eq!(g.segments.len(), 3);
        assert_eq!(g.walks.len(), 2);
        assert_eq!(g.backbone.as_deref(), Some("CHM13_backbone"));
        // spelling the backbone concatenates the three segments
        let bb = g.backbone_walk().unwrap();
        assert_eq!(g.spell(bb).unwrap().as_bytes(), b"ACGTACGTTTTTGGGGCCCC");
        // re-parse of our own writer is structurally identical
        let g2 = parse(&write(&g));
        assert_eq!(g.segments.len(), g2.segments.len());
        assert_eq!(g.walks.len(), g2.walks.len());
        assert_eq!(g.spell(bb).unwrap(), g2.spell(g2.backbone_walk().unwrap()).unwrap());
    }

    #[test]
    fn parse_rgfa_synthesizes_backbone() {
        // minigraph-style rGFA: no P/W lines, reference marked by SR:i:0 + SO offsets.
        let text = "S\ts1\tACGT\tLN:i:4\tSN:Z:chr6_alt\tSO:i:0\tSR:i:0\n\
                    S\ts2\tGGGG\tLN:i:4\tSN:Z:chr6_alt\tSO:i:4\tSR:i:0\n\
                    S\ts3\tTT\tLN:i:2\tSR:i:1\n\
                    L\ts1\t+\ts2\t+\t0M\n\
                    L\ts1\t+\ts3\t+\t0M\n";
        let g = parse(text);
        assert_eq!(g.segments.len(), 3);
        assert_eq!(g.walks.len(), 1, "backbone synthesized from rank-0 segments");
        let bb = g.backbone_walk().unwrap();
        assert_eq!(bb.steps.len(), 2, "only rank-0 segments on the spine");
        // spelled in SO order: s1 then s2
        assert_eq!(g.spell(bb).unwrap().as_bytes(), b"ACGTGGGG");
        assert!(g.backbone.as_deref().unwrap().contains("chr6_alt"));
    }

    #[test]
    fn parse_w_lines() {
        let text = "S\t1\tAAAA\nS\t2\tCCCC\nW\tHG002\t1\tchr6\t0\t8\t>1>2\n";
        let g = parse(text);
        assert_eq!(g.walks.len(), 1);
        assert_eq!(g.walks[0].name, "HG002#1#chr6");
        assert_eq!(g.spell(&g.walks[0]).unwrap().as_bytes(), b"AAAACCCC");
    }
}
