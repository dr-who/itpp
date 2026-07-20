//! ITPP genome browser — WASM interface.
//!
//! The heavy lifting (container parsing, query, subgraph layout, translation) lives in
//! [`browser`] and [`codon`] as plain, natively-tested Rust. This module is a thin
//! `wasm-bindgen` skin: load a container's bytes, then ask for JSON the page renders as a
//! tube-map. The same `itpp-format` code that writes containers reads them here in the browser.

pub mod browser;
pub mod codon;

use browser::Browser;
use wasm_bindgen::prelude::*;

/// A loaded pangenome ready to query, exposed to JavaScript.
#[wasm_bindgen]
pub struct GenomeBrowser {
    inner: Browser,
}

#[wasm_bindgen]
impl GenomeBrowser {
    /// Parse an ITPP container (its raw bytes) into a browser.
    #[wasm_bindgen(constructor)]
    pub fn new(data: &[u8]) -> Result<GenomeBrowser, JsError> {
        let graph = itpp_format::read_container(data).map_err(|e| JsError::new(&e.to_string()))?;
        Ok(GenomeBrowser { inner: Browser::new(graph) })
    }

    /// Region descriptor for the pulldown (the backbone path name).
    #[wasm_bindgen(getter)]
    pub fn region(&self) -> String {
        self.inner.region()
    }

    #[wasm_bindgen(getter)]
    pub fn segments(&self) -> usize {
        self.inner.n_segments()
    }

    #[wasm_bindgen(getter)]
    pub fn haplotypes(&self) -> usize {
        self.inner.n_haplotypes()
    }

    /// Search the pangenome for `query` and return the matches + local subgraphs as JSON.
    pub fn query(&self, query: &str, max_hits: usize, radius: usize) -> String {
        self.inner.query_json(query, max_hits, radius)
    }
}

/// wasm-bindgen smoke tests (run under the wasm test runner; skipped on native `cargo test`,
/// where the `browser`/`codon` modules are covered by ordinary unit tests).
#[cfg(all(test, target_arch = "wasm32"))]
mod wasm_tests {
    use super::*;
    use wasm_bindgen_test::*;

    fn container() -> Vec<u8> {
        let g = itpp_ingest::synth::generate(&itpp_ingest::synth::SynthParams {
            haplotypes: 4,
            backbone_blocks: 12,
            seed: 1,
            ..Default::default()
        });
        itpp_format::write_container(&g)
    }

    #[wasm_bindgen_test]
    fn load_and_query() {
        let b = GenomeBrowser::new(&container()).unwrap();
        assert!(b.segments() > 0);
        assert!(b.haplotypes() == 5);
        let json = b.query("ACGT", 5, 2);
        assert!(json.contains("\"hits\""));
    }
}
