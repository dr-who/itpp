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

## Tiled browser (streams a full chromosome)

`gui/tiled-browser.html` fetches a **multi-resolution tile pyramid** (built by `itpp-tile`) and
streams only the tiles in view — so it scales to a whole chromosome/genome, never loading more
than a screenful. Pure fetch+JS (no WASM). Serve the tiles and open it:

```sh
# generate tiles for a container
cargo run --release -p itpp-gui --bin itpp-tile -- --in database/mhc-c4.itpp --out gui/tiles/mhc-c4
# serve, then open http://localhost:8000/gui/tiled-browser.html?tiles=gui/tiles/mhc-c4/
python3 -m http.server 8000
# for chromosome-scale tiles on S3, serve them with:  rclone serve http exaba:sockjam-eedf/itpp/tiles --addr :8000
```

Zoomed out = variant-density + CNV markers; zoom in = graph nodes → proteins → nucleotides.
