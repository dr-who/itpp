//! The ITPP container: a chunked, indexed, versioned, **lossless** on-disk form of a
//! [`Graph`]. See `spec/itpp-format-v0.md`. Sequences are stored order-k range-coded (so the
//! file size is itself meaningful), the dictionary/instance annotations ride along for the
//! measure layer, and an index + footer make it seekable.
//!
//! v0 stores every segment's bases literally (compressed) — the dictionary *factoring* is a
//! measurement (see `itpp measure`), not yet the container's storage form. That keeps
//! losslessness trivial to guarantee; a later version can store instances as edit scripts.

pub mod io;

use io::{FormatError, Reader, Result, Writer};
use itpp_codec::SeqCoder;
use itpp_core::{
    Contam, DictClass, DictEntry, Edge, Graph, NodeId, RepeatInstance, Sequence, Step, Strand, Walk,
};

const MAGIC: &[u8; 4] = b"ITPP";
const FOOTER_MAGIC: &[u8; 4] = b"ITP1";
const VERSION: u16 = 0;
/// Context order used to range-code stored sequences.
const CONTAINER_ORDER: usize = 12;

const TAG_HEAD: [u8; 4] = *b"HEAD";
const TAG_DICT: [u8; 4] = *b"DICT";
const TAG_GRAF: [u8; 4] = *b"GRAF";
const TAG_SAMP: [u8; 4] = *b"SAMP";
const TAG_CTAM: [u8; 4] = *b"CTAM";
const TAG_INDX: [u8; 4] = *b"INDX";

fn class_to_u8(c: DictClass) -> u8 {
    match c {
        DictClass::Alu => 0,
        DictClass::Line => 1,
        DictClass::Sine => 2,
        DictClass::Satellite => 3,
        DictClass::Erv => 4,
        DictClass::Other => 5,
    }
}

fn class_from_u8(v: u8) -> DictClass {
    match v {
        0 => DictClass::Alu,
        1 => DictClass::Line,
        2 => DictClass::Sine,
        3 => DictClass::Satellite,
        4 => DictClass::Erv,
        _ => DictClass::Other,
    }
}

fn strand_to_u8(s: Strand) -> u8 {
    match s {
        Strand::Forward => 0,
        Strand::Reverse => 1,
    }
}

fn strand_from_u8(v: u8) -> Strand {
    if v == 0 {
        Strand::Forward
    } else {
        Strand::Reverse
    }
}

fn seq_coder() -> SeqCoder {
    SeqCoder::new(CONTAINER_ORDER)
}

// ---- section builders -------------------------------------------------------

fn build_head(g: &Graph) -> Vec<u8> {
    let mut p = Vec::new();
    p.put_u8(0); // enc = raw
    p.put_string(g.backbone.as_deref().unwrap_or(""));
    p.put_u8(CONTAINER_ORDER as u8);
    p
}

fn build_dict(g: &Graph) -> Vec<u8> {
    let coder = seq_coder();
    let mut p = Vec::new();
    p.put_u8(2); // enc = order-k
    p.put_varint(g.dictionary.entries.len() as u64);
    for e in &g.dictionary.entries {
        p.put_u32(e.id);
        p.put_u8(class_to_u8(e.class));
        p.put_string(&e.name);
        p.put_blob(&coder.encode(e.consensus.as_bytes()));
    }
    p
}

fn build_graf(g: &Graph) -> Vec<u8> {
    let coder = seq_coder();
    let mut p = Vec::new();
    // segments
    p.put_varint(g.segments.len() as u64);
    for seg in g.segments.values() {
        p.put_varint(seg.id);
        p.put_blob(&coder.encode(seg.seq.as_bytes()));
    }
    // edges
    p.put_varint(g.edges.len() as u64);
    for e in &g.edges {
        p.put_varint(e.from);
        p.put_u8(strand_to_u8(e.from_strand));
        p.put_varint(e.to);
        p.put_u8(strand_to_u8(e.to_strand));
    }
    // instances
    p.put_varint(g.instances.len() as u64);
    for (&node, inst) in &g.instances {
        p.put_varint(node);
        p.put_u32(inst.dict_id);
        p.put_u8(strand_to_u8(inst.strand));
    }
    p
}

fn build_samp(g: &Graph) -> Vec<u8> {
    let mut p = Vec::new();
    p.put_varint(g.walks.len() as u64);
    for w in &g.walks {
        p.put_string(&w.name);
        p.put_varint(w.steps.len() as u64);
        for s in &w.steps {
            p.put_varint(s.node);
            p.put_u8(strand_to_u8(s.strand));
        }
    }
    p
}

fn build_ctam(g: &Graph) -> Vec<u8> {
    let mut p = Vec::new();
    p.put_varint(g.contamination.len() as u64);
    for c in &g.contamination {
        p.put_string(&c.sample);
        p.put_varint(c.node);
        p.put_string(&c.taxon);
    }
    p
}

/// Serialize a [`Graph`] into an ITPP container.
#[must_use]
pub fn write_container(g: &Graph) -> Vec<u8> {
    let mut sections: Vec<([u8; 4], Vec<u8>)> = vec![
        (TAG_HEAD, build_head(g)),
        (TAG_DICT, build_dict(g)),
        (TAG_GRAF, build_graf(g)),
        (TAG_SAMP, build_samp(g)),
    ];
    let has_ctam = !g.contamination.is_empty();
    if has_ctam {
        sections.push((TAG_CTAM, build_ctam(g)));
    }

    let mut out = Vec::new();
    out.put_tag(MAGIC);
    out.put_u16(VERSION);
    let flags: u16 = if has_ctam { 0b01 } else { 0 };
    out.put_u16(flags);
    out.put_u64(0); // reserved

    let mut index: Vec<([u8; 4], u64, u64)> = Vec::new();
    for (tag, payload) in &sections {
        let offset = out.len() as u64;
        out.put_tag(tag);
        out.put_u64(payload.len() as u64);
        out.extend_from_slice(payload);
        index.push((*tag, offset, payload.len() as u64));
    }

    // INDX chunk
    let index_offset = out.len() as u64;
    let mut idx_payload = Vec::new();
    idx_payload.put_varint(index.len() as u64);
    for (tag, off, len) in &index {
        idx_payload.put_tag(tag);
        idx_payload.put_u64(*off);
        idx_payload.put_u64(*len);
    }
    out.put_tag(&TAG_INDX);
    out.put_u64(idx_payload.len() as u64);
    out.extend_from_slice(&idx_payload);

    // Footer
    out.put_u64(index_offset);
    out.put_tag(FOOTER_MAGIC);
    out.put_u32(0);
    out
}

// ---- reading ----------------------------------------------------------------

/// Parse an ITPP container back into a [`Graph`] (the exact graph that was written).
pub fn read_container(data: &[u8]) -> Result<Graph> {
    if data.len() < 16 + 16 {
        return Err(FormatError::Truncated);
    }
    let mut r = Reader::new(data);
    if r.tag()? != *MAGIC {
        return Err(FormatError::BadMagic);
    }
    let version = r.u16()?;
    if version != VERSION {
        return Err(FormatError::BadVersion(version));
    }
    let _flags = r.u16()?;
    let _reserved = r.u64()?;

    // Footer → index.
    let footer_off = data.len() - 16;
    r.seek(footer_off);
    let index_offset = r.u64()? as usize;
    if r.tag()? != *FOOTER_MAGIC {
        return Err(FormatError::BadMagic);
    }

    // Index.
    r.seek(index_offset);
    if r.tag()? != TAG_INDX {
        return Err(FormatError::MissingSection("INDX"));
    }
    let _idx_len = r.u64()?;
    let n = r.varint()?;
    let mut sections: Vec<([u8; 4], usize, usize)> = Vec::new();
    for _ in 0..n {
        let tag = r.tag()?;
        let off = r.u64()? as usize;
        let len = r.u64()? as usize;
        sections.push((tag, off, len));
    }

    let find = |want: [u8; 4]| sections.iter().find(|(t, _, _)| *t == want).copied();
    let payload = |off: usize| -> Result<&[u8]> {
        let mut cr = Reader::new(data);
        cr.seek(off);
        let _tag = cr.tag()?;
        let len = cr.u64()? as usize;
        let start = cr.pos;
        data.get(start..start + len).ok_or(FormatError::Truncated)
    };

    let mut g = Graph::new();

    // HEAD
    let (_, ho, _) = find(TAG_HEAD).ok_or(FormatError::MissingSection("HEAD"))?;
    {
        let mut hr = Reader::new(payload(ho)?);
        let _enc = hr.u8()?;
        let backbone = hr.string()?;
        let _order = hr.u8()?;
        g.backbone = if backbone.is_empty() { None } else { Some(backbone) };
    }

    // DICT
    if let Some((_, o, _)) = find(TAG_DICT) {
        let mut dr = Reader::new(payload(o)?);
        let _enc = dr.u8()?;
        let count = dr.varint()?;
        for _ in 0..count {
            let id = dr.u32()?;
            let class = class_from_u8(dr.u8()?);
            let name = dr.string()?;
            let consensus = Sequence::from_bytes(SeqCoder::decode(dr.blob()?));
            g.dictionary.entries.push(DictEntry { id, class, name, consensus });
        }
    }

    // GRAF
    let (_, go, _) = find(TAG_GRAF).ok_or(FormatError::MissingSection("GRAF"))?;
    {
        let mut gr = Reader::new(payload(go)?);
        let seg_count = gr.varint()?;
        for _ in 0..seg_count {
            let id: NodeId = gr.varint()?;
            let seq = Sequence::from_bytes(SeqCoder::decode(gr.blob()?));
            g.add_segment(id, seq);
        }
        let edge_count = gr.varint()?;
        for _ in 0..edge_count {
            let from = gr.varint()?;
            let from_strand = strand_from_u8(gr.u8()?);
            let to = gr.varint()?;
            let to_strand = strand_from_u8(gr.u8()?);
            g.add_edge(Edge { from, from_strand, to, to_strand });
        }
        let inst_count = gr.varint()?;
        for _ in 0..inst_count {
            let node = gr.varint()?;
            let dict_id = gr.u32()?;
            let strand = strand_from_u8(gr.u8()?);
            g.instances.insert(node, RepeatInstance { dict_id, strand });
        }
    }

    // SAMP
    let (_, so, _) = find(TAG_SAMP).ok_or(FormatError::MissingSection("SAMP"))?;
    {
        let mut sr = Reader::new(payload(so)?);
        let wcount = sr.varint()?;
        for _ in 0..wcount {
            let name = sr.string()?;
            let scount = sr.varint()?;
            let mut steps = Vec::with_capacity(scount as usize);
            for _ in 0..scount {
                let node = sr.varint()?;
                let strand = strand_from_u8(sr.u8()?);
                steps.push(Step { node, strand });
            }
            g.add_walk(Walk { name, steps });
        }
    }

    // CTAM (optional)
    if let Some((_, o, _)) = find(TAG_CTAM) {
        let mut cr = Reader::new(payload(o)?);
        let count = cr.varint()?;
        for _ in 0..count {
            let sample = cr.string()?;
            let node = cr.varint()?;
            let taxon = cr.string()?;
            g.contamination.push(Contam { sample, node, taxon });
        }
    }

    Ok(g)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_graph() -> Graph {
        let mut g = Graph::new();
        g.add_segment(0, Sequence::from_str("ACGTACGTACGTACGT"));
        g.add_segment(1, Sequence::from_str("TTGGCCAA"));
        g.add_segment(2, Sequence::from_str("ACGTACGTACGTACGA")); // alt of 0
        g.add_edge(Edge { from: 0, from_strand: Strand::Forward, to: 1, to_strand: Strand::Forward });
        g.dictionary.entries.push(DictEntry {
            id: 0,
            class: DictClass::Alu,
            name: "famA".into(),
            consensus: Sequence::from_str("ACGTACGTACGTACGT"),
        });
        g.instances.insert(2, RepeatInstance { dict_id: 0, strand: Strand::Forward });
        g.add_walk(Walk {
            name: "CHM13_backbone".into(),
            steps: vec![
                Step { node: 0, strand: Strand::Forward },
                Step { node: 1, strand: Strand::Forward },
            ],
        });
        g.add_walk(Walk {
            name: "s1".into(),
            steps: vec![
                Step { node: 2, strand: Strand::Forward },
                Step { node: 1, strand: Strand::Reverse },
            ],
        });
        g.contamination.push(Contam { sample: "s1".into(), node: 1, taxon: "E.coli".into() });
        g.backbone = Some("CHM13_backbone".into());
        g
    }

    #[test]
    fn container_roundtrip() {
        let g = sample_graph();
        let bytes = write_container(&g);
        let g2 = read_container(&bytes).expect("read");
        assert_eq!(g, g2);
        // spelled sequences survive
        for w in &g.walks {
            assert_eq!(g.spell(w), g2.spell(g2.walk_by_name(&w.name).unwrap()));
        }
    }

    #[test]
    fn rejects_garbage() {
        assert!(read_container(b"not an itpp file at all........................").is_err());
    }
}
