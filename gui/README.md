# gui — ITPP genome browser (Rust → WASM)

Target UX: open `gui/genome-browser.html` in a browser, pick a chromosome/region from a
pulldown, type a nucleotide query (e.g. `ATCG`), and see:

- **one row per match** across the pangenome's haplotypes;
- a **left→right graph view** of the local context — divergences, insertions, deletions,
  parallel alleles (bubbles) — with mouse-over positional information;
- an **above/below panel**: the 3-mers of the match and their translation to protein, plus
  other views.

Data is loaded from **our own project only** (ITPP containers committed to / released from
`github.com/dr-who/itpp`). The `itpp-core` / `itpp-codec` / `itpp-format` crates are pure Rust
and compile to `wasm32-unknown-unknown`, so the same code that writes containers reads them in
the browser.

Status: scaffolding. Rendering stack + build/delivery are being finalized; this folder will hold
the wasm crate, the query/layout logic (with tests), and the self-contained `genome-browser.html`.
