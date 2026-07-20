//! The genome-browser engine over an ITPP container.
//!
//! Given a nucleotide query (e.g. `ATCG`) it finds every occurrence across the pangenome's
//! segments (both strands), and for each one extracts the **local subgraph** — the nodes
//! before and after the hit and the bubbles between them — laid out left→right so the UI can
//! draw a tube-map. Backbone nodes carry a real genomic coordinate for mouse-over. It also
//! translates the matched context into 3-mers and protein for the codon panel.
//!
//! Everything here is plain Rust and unit-tested natively; `lib.rs` wraps it for WASM.

use crate::codon;
use itpp_core::{seq::reverse_complement, Graph, NodeId, Strand};
use std::collections::{BTreeMap, BTreeSet};

/// One node placed in the tube-map.
pub struct LNode {
    pub id: NodeId,
    pub x: i32,
    pub lane: i32,
    pub len: usize,
    pub seq: String,
    pub is_hit: bool,
    pub backbone: bool,
    pub pos: String,
}

/// One match of the query and its local graph context.
pub struct Hit {
    pub node: NodeId,
    pub offset: usize,
    pub strand: Strand,
    pub pos: String,
    pub haplotypes: Vec<String>,
    pub nodes: Vec<LNode>,
    pub edges: Vec<(NodeId, NodeId)>,
    pub codons: Vec<String>,
    pub protein: String,
}

pub struct Browser {
    graph: Graph,
    backbone_nodes: BTreeSet<NodeId>,
    backbone_off: BTreeMap<NodeId, usize>, // backbone node -> cumulative bases before it
    chrom: String,
    start_coord: i64,
    succ: BTreeMap<NodeId, Vec<NodeId>>,
    pred: BTreeMap<NodeId, Vec<NodeId>>,
    node_walks: BTreeMap<NodeId, Vec<usize>>,
}

impl Browser {
    #[must_use]
    pub fn new(graph: Graph) -> Self {
        let backbone_nodes = graph.backbone_nodes();

        // cumulative offset of each backbone node along the spine + parse genomic coord
        let mut backbone_off = BTreeMap::new();
        let mut cum = 0usize;
        if let Some(bw) = graph.backbone_walk() {
            for step in &bw.steps {
                backbone_off.entry(step.node).or_insert(cum);
                cum += graph.segments.get(&step.node).map_or(0, |s| s.seq.len());
            }
        }
        let (chrom, start_coord) = parse_region(graph.backbone.as_deref().unwrap_or(""));

        // adjacency + node→walks from the observed walks (real haplotype context)
        let mut succ: BTreeMap<NodeId, Vec<NodeId>> = BTreeMap::new();
        let mut pred: BTreeMap<NodeId, Vec<NodeId>> = BTreeMap::new();
        let mut node_walks: BTreeMap<NodeId, Vec<usize>> = BTreeMap::new();
        for (wi, w) in graph.walks.iter().enumerate() {
            for pair in w.steps.windows(2) {
                push_unique(&mut succ, pair[0].node, pair[1].node);
                push_unique(&mut pred, pair[1].node, pair[0].node);
            }
            for s in &w.steps {
                let e = node_walks.entry(s.node).or_default();
                if e.last() != Some(&wi) {
                    e.push(wi);
                }
            }
        }

        Browser { graph, backbone_nodes, backbone_off, chrom, start_coord, succ, pred, node_walks }
    }

    /// Human-readable region descriptor for the pulldown, e.g. `chr6:31,825,251 (MHC-C4)`.
    #[must_use]
    pub fn region(&self) -> String {
        self.graph.backbone.clone().unwrap_or_else(|| "(no backbone)".into())
    }

    #[must_use]
    pub fn n_segments(&self) -> usize {
        self.graph.segments.len()
    }

    #[must_use]
    pub fn n_haplotypes(&self) -> usize {
        self.graph.walks.len()
    }

    /// Find the query across all segments (both strands); return up to `max_hits` hits, each
    /// with a local subgraph laid out within `radius` graph-steps up/downstream.
    #[must_use]
    pub fn query(&self, raw: &str, max_hits: usize, radius: usize) -> Vec<Hit> {
        let q: Vec<u8> = raw.trim().bytes().map(|b| b.to_ascii_uppercase()).collect();
        if q.is_empty() {
            return Vec::new();
        }
        let rc = reverse_complement(&q);
        let mut hits = Vec::new();

        for seg in self.graph.segments.values() {
            let hay = seg.seq.as_bytes();
            for (strand, needle) in [(Strand::Forward, &q), (Strand::Reverse, &rc)] {
                // reverse strand only distinct if not a palindrome
                if strand == Strand::Reverse && rc == q {
                    continue;
                }
                let mut from = 0;
                while let Some(rel) = find_sub(&hay[from..], needle) {
                    let offset = from + rel;
                    hits.push(self.build_hit(seg.id, offset, strand, radius, radius));
                    if hits.len() >= max_hits {
                        return hits;
                    }
                    from = offset + 1;
                }
            }
        }
        hits
    }

    /// Run [`Browser::query`] and serialize the results to JSON for the UI. Each node carries
    /// `x` (column) and `backbone`, so the UI can collapse the variant lanes of a column into a
    /// single glyph and expand them on click.
    #[must_use]
    pub fn query_json(&self, raw: &str, max_hits: usize, radius: usize) -> String {
        let hits = self.query(raw, max_hits, radius);
        let total = self.count_matches(raw);
        let q_up: String = raw.trim().to_ascii_uppercase();
        let mut s = String::from("{");
        s.push_str(&format!("\"region\":{},", jstr(&self.region())));
        s.push_str(&format!("\"chrom\":{},", jstr(&self.chrom)));
        s.push_str(&format!("\"query\":{},", jstr(&q_up)));
        s.push_str(&format!("\"n_total\":{},", total));
        s.push_str(&format!("\"n_hits\":{},", hits.len()));
        s.push_str("\"hits\":[");
        for (i, h) in hits.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str(&hit_to_json(h));
        }
        s.push_str("]}");
        s
    }

    /// Build the JSON for the local subgraph around one already-found match, with independent
    /// upstream/downstream radii. Powers the ◀/▶ "extend context" buttons. Returns `{}` if the
    /// node isn't in the graph.
    #[must_use]
    pub fn context_json(
        &self,
        node: NodeId,
        offset: usize,
        strand_sign: &str,
        left: usize,
        right: usize,
    ) -> String {
        if !self.graph.segments.contains_key(&node) {
            return "{}".into();
        }
        let strand = if strand_sign == "-" { Strand::Reverse } else { Strand::Forward };
        let off = offset.min(self.graph.segments[&node].seq.len().saturating_sub(1));
        let hit = self.build_hit(node, off, strand, left, right);
        format!("{{\"hit\":{}}}", hit_to_json(&hit))
    }

    fn build_hit(
        &self,
        node: NodeId,
        offset: usize,
        strand: Strand,
        left: usize,
        right: usize,
    ) -> Hit {
        // BFS upstream (pred) and downstream (succ) collecting depth per node.
        let mut depth: BTreeMap<NodeId, i32> = BTreeMap::new();
        depth.insert(node, 0);
        self.bfs(node, left, true, &mut depth);
        self.bfs(node, right, false, &mut depth);

        // group by x (depth) then assign lanes: backbone at lane 0, others alternate.
        let mut by_x: BTreeMap<i32, Vec<NodeId>> = BTreeMap::new();
        for (&n, &d) in &depth {
            by_x.entry(d).or_default().push(n);
        }
        let mut nodes = Vec::new();
        for (&x, ids) in &by_x {
            let mut ordered = ids.clone();
            ordered.sort_by_key(|n| (!self.backbone_nodes.contains(n), *n));
            let mut alt = 1;
            for n in ordered {
                let lane = if self.backbone_nodes.contains(&n) {
                    0
                } else {
                    let l = if alt % 2 == 1 { (alt + 1) / 2 } else { -(alt / 2) };
                    alt += 1;
                    l
                };
                let seg = &self.graph.segments[&n];
                nodes.push(LNode {
                    id: n,
                    x,
                    lane,
                    len: seg.seq.len(),
                    seq: preview(seg.seq.as_bytes()),
                    is_hit: n == node,
                    backbone: self.backbone_nodes.contains(&n),
                    pos: self.node_pos(n),
                });
            }
        }

        // edges internal to the collected node set
        let present: BTreeSet<NodeId> = depth.keys().copied().collect();
        let mut edges = Vec::new();
        for (&a, outs) in &self.succ {
            if present.contains(&a) {
                for &b in outs {
                    if present.contains(&b) {
                        edges.push((a, b));
                    }
                }
            }
        }

        // haplotypes through the hit node
        let haplotypes: Vec<String> = self
            .node_walks
            .get(&node)
            .map(|idx| idx.iter().map(|&i| sample_name(&self.graph.walks[i].name)).collect())
            .unwrap_or_default();

        // codon/protein panel: translate the matched context on the hit node
        let seg = &self.graph.segments[&node];
        let window_end = (offset + 60).min(seg.seq.len());
        let context = &seg.seq.as_bytes()[offset..window_end];
        let codons = codon::codons(context, 0);
        let protein = codon::translate(context, 0);

        Hit {
            node,
            offset,
            strand,
            pos: self.hit_pos(node, offset),
            haplotypes,
            nodes,
            edges,
            codons,
            protein,
        }
    }

    /// Lay out the WHOLE graph in genomic coordinates for the zoomable overview ("see the
    /// world"): the backbone is a straight left→right axis (each node at its cumulative bp
    /// offset, lane 0); variants sit at their nearest backbone anchor, packed into lanes
    /// above/below. `x` is a bp offset from the region start; the UI scales/pans it.
    #[must_use]
    pub fn overview_json(&self) -> String {
        let mut pos: BTreeMap<NodeId, (i64, i32)> = BTreeMap::new();
        for (&n, &off) in &self.backbone_off {
            pos.insert(n, (off as i64, 0));
        }
        // variants sorted by anchor x, packed into alternating lanes by interval
        let mut variants: Vec<(NodeId, i64, usize)> = self
            .graph
            .segments
            .values()
            .filter(|s| !self.backbone_nodes.contains(&s.id))
            .map(|s| (s.id, self.anchor_x(s.id), s.seq.len()))
            .collect();
        variants.sort_by_key(|v| v.1);
        let mut lane_end: Vec<i64> = Vec::new();
        for (id, x, len) in &variants {
            let li = (0..lane_end.len()).find(|&li| lane_end[li] <= *x).unwrap_or_else(|| {
                lane_end.push(0);
                lane_end.len() - 1
            });
            lane_end[li] = x + (*len as i64).max(1) + 24;
            let lane = if li % 2 == 0 { (li / 2 + 1) as i32 } else { -((li / 2 + 1) as i32) };
            pos.insert(*id, (*x, lane));
        }

        let span = self.backbone_off.values().copied().max().unwrap_or(0)
            + self
                .graph
                .backbone_walk()
                .and_then(|w| w.steps.last())
                .and_then(|s| self.graph.segments.get(&s.node))
                .map_or(0, |s| s.seq.len());

        let mut s = String::from("{");
        s.push_str(&format!("\"chrom\":{},", jstr(&self.chrom)));
        s.push_str(&format!("\"start\":{},", self.start_coord));
        s.push_str(&format!("\"span\":{},", span));
        s.push_str(&format!("\"n_nodes\":{},", self.graph.segments.len()));
        s.push_str("\"nodes\":[");
        let mut first = true;
        for seg in self.graph.segments.values() {
            let (x, lane) = pos.get(&seg.id).copied().unwrap_or((0, 0));
            if !first {
                s.push(',');
            }
            first = false;
            // full sequence so the UI can render actual nucleotides at base-level zoom; capped
            // per node to bound payload (long nodes show their first BASES_CAP bases).
            const BASES_CAP: usize = 4000;
            let full = seg.seq.as_bytes();
            let seq_out = std::str::from_utf8(&full[..full.len().min(BASES_CAP)]).unwrap_or("");
            s.push_str(&format!(
                "{{\"id\":{},\"x\":{},\"lane\":{},\"len\":{},\"bb\":{},\"kind\":{},\"seq\":{}}}",
                seg.id,
                x,
                lane,
                seg.seq.len(),
                self.backbone_nodes.contains(&seg.id),
                jstr(self.node_kind(seg.id, full)),
                jstr(seq_out)
            ));
        }
        s.push_str("],\"edges\":[");
        for (i, e) in self.graph.edges.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str(&format!("[{},{}]", e.from, e.to));
        }
        s.push_str("]}");
        s
    }

    /// Coarse node type for colouring, computed from our own code (no RepeatMasker yet):
    /// backbone spine, tandem / low-complexity repeat, large interspersed-repeat or SV
    /// candidate (Alu-sized and up), indel, or SNP. Biological family calls (Alu/LINE/…)
    /// come later from the annotation layer.
    fn node_kind(&self, id: NodeId, seq: &[u8]) -> &'static str {
        if self.backbone_nodes.contains(&id) {
            return "backbone";
        }
        if is_tandem(seq) {
            return "tandem";
        }
        match seq.len() {
            n if n >= 250 => "large",
            n if n <= 4 => "snp",
            _ => "indel",
        }
    }

    /// x (bp offset) of the nearest backbone anchor to a variant node.
    fn anchor_x(&self, v: NodeId) -> i64 {
        // nearest backbone upstream → its end; else nearest downstream → its start; else 0
        for (adj, downstream) in [(&self.pred, false), (&self.succ, true)] {
            let mut seen: BTreeSet<NodeId> = BTreeSet::new();
            let mut frontier = vec![v];
            for _ in 0..8 {
                let mut next = Vec::new();
                for n in frontier.drain(..) {
                    if let Some(neigh) = adj.get(&n) {
                        for &m in neigh {
                            if let Some(&off) = self.backbone_off.get(&m) {
                                let len = self.graph.segments.get(&m).map_or(0, |s| s.seq.len());
                                return if downstream { off as i64 } else { off as i64 + len as i64 };
                            }
                            if seen.insert(m) {
                                next.push(m);
                            }
                        }
                    }
                }
                frontier = next;
            }
        }
        0
    }

    /// Node ids (and coords) of up to `max` matches of the query, for highlighting on the
    /// overview and building "jump to" targets.
    #[must_use]
    pub fn matches_json(&self, raw: &str, max: usize) -> String {
        let q: Vec<u8> = raw.trim().bytes().map(|b| b.to_ascii_uppercase()).collect();
        let mut s = String::from("{\"matches\":[");
        if q.is_empty() {
            s.push_str("]}");
            return s;
        }
        let rc = reverse_complement(&q);
        let both = rc != q;
        let mut count = 0;
        let mut first = true;
        'outer: for seg in self.graph.segments.values() {
            let hay = seg.seq.as_bytes();
            for needle in std::iter::once(&q).chain(if both { Some(&rc) } else { None }) {
                if find_sub(hay, needle).is_some() {
                    if !first {
                        s.push(',');
                    }
                    first = false;
                    s.push_str(&format!(
                        "{{\"node\":{},\"pos\":{}}}",
                        seg.id,
                        jstr(&self.hit_pos(seg.id, 0))
                    ));
                    count += 1;
                    if count >= max {
                        break 'outer;
                    }
                    break; // one entry per segment is enough to highlight it
                }
            }
        }
        s.push_str(&format!("],\"n_total\":{}}}", self.count_matches(raw)));
        s
    }

    /// Total occurrences of the query (both strands) across all segments — so the UI can say
    /// "showing 40 of N". Short queries (e.g. a 3-mer) match thousands of times.
    #[must_use]
    pub fn count_matches(&self, raw: &str) -> usize {
        let q: Vec<u8> = raw.trim().bytes().map(|b| b.to_ascii_uppercase()).collect();
        if q.is_empty() {
            return 0;
        }
        let rc = reverse_complement(&q);
        let both = rc != q;
        let mut total = 0;
        for seg in self.graph.segments.values() {
            let hay = seg.seq.as_bytes();
            total += count_sub(hay, &q);
            if both {
                total += count_sub(hay, &rc);
            }
        }
        total
    }

    fn bfs(&self, start: NodeId, radius: usize, upstream: bool, depth: &mut BTreeMap<NodeId, i32>) {
        let adj = if upstream { &self.pred } else { &self.succ };
        let mut frontier = vec![start];
        for step in 1..=radius as i32 {
            let mut next = Vec::new();
            for n in frontier.drain(..) {
                if let Some(neigh) = adj.get(&n) {
                    for &m in neigh {
                        let d = if upstream { -step } else { step };
                        if let std::collections::btree_map::Entry::Vacant(e) = depth.entry(m) {
                            e.insert(d);
                            next.push(m);
                        }
                    }
                }
            }
            frontier = next;
            if frontier.is_empty() {
                break;
            }
        }
    }

    fn node_pos(&self, n: NodeId) -> String {
        if let Some(&off) = self.backbone_off.get(&n) {
            format!("{}:{}", self.chrom, commas(self.start_coord + off as i64))
        } else {
            "(variant allele)".into()
        }
    }

    fn hit_pos(&self, n: NodeId, offset: usize) -> String {
        if let Some(&off) = self.backbone_off.get(&n) {
            let c = self.start_coord + off as i64 + offset as i64;
            format!("{}:{}", self.chrom, commas(c))
        } else {
            format!("variant node {n} +{offset}")
        }
    }
}

// ---- helpers ----------------------------------------------------------------

fn jstr(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            _ => out.push(c),
        }
    }
    out.push('"');
    out
}

fn hit_to_json(h: &Hit) -> String {
    let mut s = String::from("{");
    s.push_str(&format!("\"node\":{},", h.node));
    s.push_str(&format!("\"offset\":{},", h.offset));
    s.push_str(&format!("\"strand\":\"{}\",", h.strand.as_sign()));
    s.push_str(&format!("\"pos\":{},", jstr(&h.pos)));
    s.push_str(&format!("\"n_haplotypes\":{},", h.haplotypes.len()));
    s.push_str("\"haplotypes\":[");
    for (i, hp) in h.haplotypes.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str(&jstr(hp));
    }
    s.push_str("],");
    s.push_str(&format!("\"protein\":{},", jstr(&h.protein)));
    s.push_str("\"codons\":[");
    for (i, c) in h.codons.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str(&jstr(c));
    }
    s.push_str("],\"nodes\":[");
    for (i, n) in h.nodes.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str(&format!(
            "{{\"id\":{},\"x\":{},\"lane\":{},\"len\":{},\"seq\":{},\"hit\":{},\"backbone\":{},\"pos\":{}}}",
            n.id, n.x, n.lane, n.len, jstr(&n.seq), n.is_hit, n.backbone, jstr(&n.pos)
        ));
    }
    s.push_str("],\"edges\":[");
    for (i, (a, b)) in h.edges.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str(&format!("[{a},{b}]"));
    }
    s.push_str("]}");
    s
}

fn push_unique(m: &mut BTreeMap<NodeId, Vec<NodeId>>, k: NodeId, v: NodeId) {
    let e = m.entry(k).or_default();
    if !e.contains(&v) {
        e.push(v);
    }
}

fn find_sub(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || needle.len() > hay.len() {
        return None;
    }
    hay.windows(needle.len()).position(|w| w == needle)
}

fn count_sub(hay: &[u8], needle: &[u8]) -> usize {
    if needle.is_empty() || needle.len() > hay.len() {
        return 0;
    }
    hay.windows(needle.len()).filter(|w| *w == needle).count()
}

/// True if the sequence is mostly periodic with a short period (homopolymer, di-/tri-…nucleotide
/// tandem repeat, low-complexity) — i.e. satellite-like repetition, detected without a library.
fn is_tandem(seq: &[u8]) -> bool {
    let n = seq.len();
    if n < 6 {
        return false;
    }
    for period in 1..=6usize {
        if n < period * 4 {
            break;
        }
        let matches = (period..n).filter(|&i| seq[i] == seq[i - period]).count();
        if matches as f64 / (n - period) as f64 > 0.82 {
            return true;
        }
    }
    false
}

fn preview(seq: &[u8]) -> String {
    if seq.len() <= 24 {
        String::from_utf8_lossy(seq).into_owned()
    } else {
        format!(
            "{}…{}",
            String::from_utf8_lossy(&seq[..12]),
            String::from_utf8_lossy(&seq[seq.len() - 8..])
        )
    }
}

/// `chm13#chr6:31825251-31908851` → (`chr6`, 31825251).
fn parse_region(name: &str) -> (String, i64) {
    // take the substring after the last '#', then split on ':' and '-'
    let tail = name.rsplit('#').next().unwrap_or(name);
    if let Some((chrom, rest)) = tail.split_once(':') {
        let start = rest.split('-').next().and_then(|s| s.parse::<i64>().ok()).unwrap_or(0);
        (chrom.to_string(), start)
    } else {
        (tail.to_string(), 0)
    }
}

/// `HG00438#1#JAHBCB…` → `HG00438#1` (sample + haplotype, dropping the contig coordinates).
fn sample_name(walk: &str) -> String {
    let parts: Vec<&str> = walk.split('#').collect();
    match parts.len() {
        0 => walk.to_string(),
        1 => parts[0].to_string(),
        _ => format!("{}#{}", parts[0], parts[1]),
    }
}

fn commas(n: i64) -> String {
    let neg = n < 0;
    let s = n.unsigned_abs().to_string();
    let b = s.as_bytes();
    let mut out = String::new();
    for (i, c) in b.iter().enumerate() {
        if i > 0 && (b.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(*c as char);
    }
    if neg {
        format!("-{out}")
    } else {
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use itpp_ingest::synth::{generate, SynthParams};

    fn browser() -> Browser {
        let g = generate(&SynthParams { haplotypes: 6, backbone_blocks: 20, seed: 3, ..Default::default() });
        Browser::new(g)
    }

    #[test]
    fn region_and_counts() {
        let b = browser();
        assert!(b.region().contains("backbone"));
        assert!(b.n_segments() > 20);
        assert_eq!(b.n_haplotypes(), 7);
    }

    #[test]
    fn query_finds_a_known_substring() {
        let b = browser();
        // pull a real 8-mer out of the first backbone segment and search for it
        let g = generate(&SynthParams { haplotypes: 6, backbone_blocks: 20, seed: 3, ..Default::default() });
        let first = g.backbone_walk().unwrap().steps[0].node;
        let sub = String::from_utf8(g.segments[&first].seq.as_bytes()[10..18].to_vec()).unwrap();
        let hits = b.query(&sub, 50, 3);
        assert!(!hits.is_empty(), "should find the 8-mer somewhere");
        let h = &hits[0];
        // the hit node is laid out at x=0, and appears in the node list
        assert!(h.nodes.iter().any(|n| n.is_hit && n.x == 0));
        // codon panel is populated
        assert!(!h.codons.is_empty());
    }

    #[test]
    fn empty_query_no_hits() {
        assert!(browser().query("", 10, 2).is_empty());
    }

    #[test]
    fn overview_lays_out_whole_graph() {
        let b = browser();
        let json = b.overview_json();
        assert!(json.contains("\"nodes\":["));
        assert!(json.contains("\"span\":"));
        assert!(json.contains("\"bb\":true")); // backbone nodes present
        // every segment is represented
        assert_eq!(json.matches("\"id\":").count(), b.n_segments());
        // matches highlighting returns node ids
        let g = generate(&SynthParams { haplotypes: 6, backbone_blocks: 20, seed: 3, ..Default::default() });
        let first = g.backbone_walk().unwrap().steps[0].node;
        let sub = String::from_utf8(g.segments[&first].seq.as_bytes()[10..18].to_vec()).unwrap();
        let m = b.matches_json(&sub, 100);
        assert!(m.contains("\"matches\":[") && m.contains("\"node\":"));
    }

    #[test]
    fn context_widens_the_window() {
        let b = browser();
        let g = generate(&SynthParams { haplotypes: 6, backbone_blocks: 20, seed: 3, ..Default::default() });
        let mid = g.backbone_walk().unwrap().steps[10].node;
        let sub = String::from_utf8(g.segments[&mid].seq.as_bytes()[5..13].to_vec()).unwrap();
        let hit = &b.query(&sub, 1, 2)[0];
        let narrow = hit.nodes.len();
        // extend far to both sides; the laid-out subgraph should grow (or stay, never shrink)
        let wide_json = b.context_json(hit.node, hit.offset, "+", 12, 12);
        let wide_nodes = wide_json.matches("\"id\":").count();
        assert!(wide_nodes >= narrow, "wide {wide_nodes} should be >= narrow {narrow}");
        assert!(wide_json.starts_with("{\"hit\":"));
        // unknown node → empty object
        assert_eq!(b.context_json(u64::MAX, 0, "+", 2, 2), "{}");
    }

    #[test]
    fn query_json_is_wellformed() {
        let b = browser();
        let g = generate(&SynthParams { haplotypes: 6, backbone_blocks: 20, seed: 3, ..Default::default() });
        let first = g.backbone_walk().unwrap().steps[0].node;
        let sub = String::from_utf8(g.segments[&first].seq.as_bytes()[10..18].to_vec()).unwrap();
        let json = b.query_json(&sub, 10, 3);
        assert!(json.starts_with('{') && json.ends_with('}'));
        assert!(json.contains("\"hits\":["));
        assert!(json.contains("\"backbone\":"));
        assert!(json.contains("\"x\":"));
        // balanced braces — cheap structural sanity
        let opens = json.matches('{').count();
        let closes = json.matches('}').count();
        assert_eq!(opens, closes, "unbalanced braces in json");
    }

    #[test]
    fn parse_region_works() {
        assert_eq!(parse_region("chm13#chr6:31825251-31908851"), ("chr6".into(), 31825251));
        assert_eq!(parse_region("CHM13_backbone"), ("CHM13_backbone".into(), 0));
    }

    #[test]
    fn sample_name_trims_contig() {
        assert_eq!(sample_name("HG00438#1#JAHBCB010:1-2"), "HG00438#1");
        assert_eq!(sample_name("grch38"), "grch38");
    }

    #[test]
    fn commas_format() {
        assert_eq!(commas(31825251), "31,825,251");
        assert_eq!(commas(0), "0");
    }
}
