//! Importers and fixtures for ITPP.
//!
//! - [`gfa`] — minimal GFA v1 reader/writer (ingest an HPRC graph, or round-trip our own).
//! - [`fasta`] — a tiny FASTA reader.
//! - [`synth`] — a deterministic, MHC-like synthetic pangenome so the meter runs with no
//!   external toolchain.

pub mod fasta;
pub mod gfa;
pub mod synth;
