//! The measurement harness: turn a [`Graph`] into a [`Ledger`] of achieved bits/char.
//!
//! The two-part code shows up as sections: `backbone` + `graph` + `dictionary` are the shared
//! model; `repeats` and `walks` are where factoring pays off. Each sequence section is coded
//! with the internal order-k model (real bits), optionally cross-checked against `xz`/`zstd`/
//! `bzip2`. The `repeats` section is measured **both** literally and dictionary-factored so the
//! ledger's best-of shows the win. Baselines answer "what would no graph cost?".

use itpp_codec::{measure_internal, Ledger, SectionLedger, SeqCoder};
use itpp_core::{edit::anchored_diff, seq::reverse_complement, Graph, NodeId, Strand};
use std::collections::BTreeSet;
use std::io::Write;
use std::process::{Command, Stdio};

pub struct MeasureOpts {
    pub orders: Vec<usize>,
    pub external: bool,
}

impl Default for MeasureOpts {
    fn default() -> Self {
        MeasureOpts { orders: vec![2, 8, 16], external: true }
    }
}

/// Compute the full ledger for a graph.
#[must_use]
pub fn measure(g: &Graph, opts: &MeasureOpts) -> Ledger {
    let total_bases = g.walk_bases() as u64;
    let mut led = Ledger::new(total_bases);

    let bb_nodes = g.backbone_nodes();
    let bb_ids: Vec<NodeId> =
        g.backbone_walk().map(|w| w.steps.iter().map(|s| s.node).collect()).unwrap_or_default();

    // --- shared model ---
    let bb_data = concat_ids(g, &dedup(&bb_ids));
    led.push_section(seq_section("backbone", &bb_data, opts));

    let other_ids: Vec<NodeId> = g
        .segments
        .keys()
        .copied()
        .filter(|id| !bb_nodes.contains(id) && !g.instances.contains_key(id))
        .collect();
    let other_data = concat_ids(g, &other_ids);
    led.push_section(seq_section("graph", &other_data, opts));

    let dict_data: Vec<u8> =
        g.dictionary.entries.iter().flat_map(|e| e.consensus.as_bytes().to_vec()).collect();
    led.push_section(seq_section("dictionary", &dict_data, opts));

    // --- repeats: literal vs dictionary-factored ---
    let mut inst_ids: Vec<NodeId> = g.instances.keys().copied().collect();
    inst_ids.sort_unstable();
    let inst_bases: u64 = inst_ids.iter().map(|id| g.segments[id].seq.len() as u64).sum();
    let mut rsec = SectionLedger::new("repeats", inst_bases);
    let inst_literal = concat_ids(g, &inst_ids);
    for cr in measure_internal(&inst_literal, &opts.orders) {
        rsec.add(format!("literal-{}", cr.coder), cr.bits);
    }
    if opts.external {
        for (n, b) in external_all(&inst_literal) {
            rsec.add(format!("literal-{n}"), b);
        }
    }
    rsec.add("dict-factored", factored_repeat_bits(g, &inst_ids, opts));
    led.push_section(rsec);

    // --- per-sample residual: the walks (allele choices through the graph) ---
    let walk_stream = encode_walks(g);
    let steps: u64 = g.walks.iter().map(|w| w.steps.len() as u64).sum();
    let mut wsec = SectionLedger::new("walks", steps);
    wsec.add("raw", walk_stream.len() as u64 * 8);
    for &k in &opts.orders {
        wsec.add(format!("order-{k}"), SeqCoder::new(k).cost_bits(&walk_stream));
    }
    if opts.external {
        for (n, b) in external_all(&walk_stream) {
            wsec.add(n, b);
        }
    }
    led.push_section(wsec);

    // --- baselines: what would "no graph" cost? ---
    add_baselines(g, &mut led, opts);
    led
}

fn dedup(ids: &[NodeId]) -> Vec<NodeId> {
    let set: BTreeSet<NodeId> = ids.iter().copied().collect();
    set.into_iter().collect()
}

fn concat_ids(g: &Graph, ids: &[NodeId]) -> Vec<u8> {
    let mut out = Vec::new();
    for id in ids {
        if let Some(seg) = g.segments.get(id) {
            out.extend_from_slice(seg.seq.as_bytes());
        }
    }
    out
}

/// Build a section by running the internal coders (+ optional externals) on `data`.
fn seq_section(name: &str, data: &[u8], opts: &MeasureOpts) -> SectionLedger {
    let mut s = SectionLedger::new(name, data.len() as u64);
    for cr in measure_internal(data, &opts.orders) {
        s.add(cr.coder, cr.bits);
    }
    if opts.external {
        for (n, b) in external_all(data) {
            s.add(n, b);
        }
    }
    s
}

/// Cost of storing repeat instances as `(dict_id, strand, edit-script)` against the consensus.
/// The consensus itself is charged once in the `dictionary` section, so here we only pay for
/// the novel edits plus a small per-instance header.
fn factored_repeat_bits(g: &Graph, inst_ids: &[NodeId], opts: &MeasureOpts) -> u64 {
    const PER_INSTANCE_HEADER_BITS: u64 = 48; // dict id + strand + op structure, generously
    let mut literal_stream = Vec::new();
    let mut count = 0u64;
    for &id in inst_ids {
        let Some(seg) = g.segments.get(&id) else { continue };
        let Some(inst) = g.instances.get(&id) else { continue };
        let Some(entry) = g.dictionary.get(inst.dict_id) else { continue };
        let consensus = match inst.strand {
            Strand::Forward => entry.consensus.as_bytes().to_vec(),
            Strand::Reverse => reverse_complement(entry.consensus.as_bytes()),
        };
        let script = anchored_diff(&consensus, seg.seq.as_bytes());
        // gather the literal (novel) bases the script carries
        for op in &script.0 {
            if let itpp_core::EditOp::Sub { bases } | itpp_core::EditOp::Ins { bases } = op {
                literal_stream.extend_from_slice(bases);
            }
        }
        count += 1;
    }
    let best_order = opts.orders.iter().copied().max().unwrap_or(8);
    let literal_bits = SeqCoder::new(best_order).cost_bits(&literal_stream);
    literal_bits + count * PER_INSTANCE_HEADER_BITS
}

/// Serialize the non-backbone walks to a byte stream (delta-coded node ids + strand) so a
/// coder can measure the residual cost of "which path each individual takes".
fn encode_walks(g: &Graph) -> Vec<u8> {
    let backbone = g.backbone.as_deref();
    let mut out = Vec::new();
    for w in &g.walks {
        if Some(w.name.as_str()) == backbone {
            continue;
        }
        let mut prev: i64 = 0;
        for s in &w.steps {
            let delta = s.node as i64 - prev;
            prev = s.node as i64;
            put_varint(&mut out, zigzag(delta));
            out.push(match s.strand {
                Strand::Forward => 0,
                Strand::Reverse => 1,
            });
        }
    }
    out
}

fn add_baselines(g: &Graph, led: &mut Ledger, opts: &MeasureOpts) {
    // Spell the whole cohort (every haplotype) once.
    let mut cohort = Vec::new();
    for w in &g.walks {
        if let Some(seq) = g.spell(w) {
            cohort.extend_from_slice(seq.as_bytes());
        }
    }
    // naive independent 2-bit packing — the "no model at all" reference
    led.add_baseline("2bit_percopy", itpp_codec::baseline::twobit_total_bits(&cohort));
    // order-0 Shannon of the cohort
    let h0 = itpp_codec::baseline::shannon_order0_bits(&cohort);
    led.add_baseline("shannon0", (h0 * cohort.len() as f64).round() as u64);
    // the strongest "no graph" reference: PPM-compress the concatenated haplotypes directly
    let best_order = opts.orders.iter().copied().max().unwrap_or(12);
    led.add_baseline(
        format!("order{best_order}_concat_nograph"),
        SeqCoder::new(best_order).cost_bits(&cohort),
    );
    // per-genome xz, independently (no shared model)
    if opts.external {
        let mut total = 0u64;
        let mut ok = true;
        for w in &g.walks {
            if let Some(seq) = g.spell(w) {
                match external_one("xz", &["-9", "-c", "-T", "1"], seq.as_bytes()) {
                    Some(bits) => total += bits,
                    None => {
                        ok = false;
                        break;
                    }
                }
            }
        }
        if ok {
            led.add_baseline("percopy_xz", total);
        }
    }
}

// ---- varint / zigzag helpers (local to the walk encoder) --------------------

fn zigzag(v: i64) -> u64 {
    ((v << 1) ^ (v >> 63)) as u64
}

fn put_varint(out: &mut Vec<u8>, mut v: u64) {
    loop {
        let mut byte = (v & 0x7F) as u8;
        v >>= 7;
        if v != 0 {
            byte |= 0x80;
        }
        out.push(byte);
        if v == 0 {
            break;
        }
    }
}

// ---- external compressor cross-checks ---------------------------------------

/// Run `data` through every available external compressor; skip ones not installed.
fn external_all(data: &[u8]) -> Vec<(String, u64)> {
    let mut out = Vec::new();
    for (name, cmd, args) in [
        ("xz", "xz", &["-9", "-c", "-T", "1"][..]),
        ("zstd", "zstd", &["-19", "-c", "-q"][..]),
        ("bzip2", "bzip2", &["-9", "-c"][..]),
    ] {
        if let Some(bits) = external_one(cmd, args, data) {
            out.push((name.to_string(), bits));
        }
    }
    out
}

/// Pipe `data` through `cmd args`, returning compressed size in bits. A writer thread avoids
/// the classic stdin/stdout pipe deadlock. Returns `None` if the tool is missing or fails.
fn external_one(cmd: &str, args: &[&str], data: &[u8]) -> Option<u64> {
    let mut child = Command::new(cmd)
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .ok()?;
    let mut stdin = child.stdin.take()?;
    let buf = data.to_vec();
    let writer = std::thread::spawn(move || {
        let _ = stdin.write_all(&buf);
        // stdin dropped here → EOF
    });
    let output = child.wait_with_output().ok()?;
    let _ = writer.join();
    if !output.status.success() {
        return None;
    }
    Some(output.stdout.len() as u64 * 8)
}

#[cfg(test)]
mod tests {
    use super::*;
    use itpp_ingest::synth::{generate, SynthParams};

    #[test]
    fn measure_produces_positive_bpc_and_beats_naive() {
        let g = generate(&SynthParams { haplotypes: 8, backbone_blocks: 40, ..Default::default() });
        // internal coders only (external tools may be absent in CI)
        let opts = MeasureOpts { orders: vec![2, 8, 16], external: false };
        let led = measure(&g, &opts);
        let bpc = led.bits_per_char();
        assert!(bpc > 0.0, "bpc should be positive");
        // The graph two-part code must beat naive 2-bit-per-base independent storage.
        let naive = led
            .baselines
            .iter()
            .find(|b| b.coder == "2bit_percopy")
            .map(|b| b.bits)
            .unwrap();
        assert!(
            led.total_bits() < naive,
            "graph total {} should beat 2bit percopy {}",
            led.total_bits(),
            naive
        );
    }

    #[test]
    fn dictionary_factoring_beats_literal_repeats() {
        let g = generate(&SynthParams::default());
        let opts = MeasureOpts { orders: vec![8, 16], external: false };
        let led = measure(&g, &opts);
        let repeats = led.sections.iter().find(|s| s.name == "repeats").unwrap();
        let best = repeats.best().unwrap();
        assert_eq!(best.coder, "dict-factored", "factoring should win on repeats");
    }
}
