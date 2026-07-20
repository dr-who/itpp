#!/usr/bin/env bash
# Build the self-contained gui/genome-browser.html.
#
# Compiles the Rust engine to WASM, generates wasm-bindgen glue (no-modules, so the page runs
# from file:// on double-click), and inlines the glue + WASM + default container as base64 into
# a single HTML file. Also builds a nodejs binding used by gui/tests/smoke.mjs.
set -euo pipefail
cd "$(dirname "$0")/.."   # repo root

CONTAINER="${1:-database/mhc-c4.itpp}"
ANNO="${2:-database/annotations/c4-clinvar.json}"

echo ">> cargo build --release --target wasm32-unknown-unknown"
cargo build -p itpp-gui --target wasm32-unknown-unknown --release

echo ">> wasm-bindgen (no-modules for the page, nodejs for tests)"
wasm-bindgen --target no-modules --no-typescript --out-dir gui/pkg \
  target/wasm32-unknown-unknown/release/itpp_gui.wasm
wasm-bindgen --target nodejs --no-typescript --out-dir gui/pkg-node \
  target/wasm32-unknown-unknown/release/itpp_gui.wasm

echo ">> inlining into gui/genome-browser.html (container: $CONTAINER, annotations: $ANNO)"
python3 - "$CONTAINER" "$ANNO" <<'PY'
import base64, sys, pathlib
root = pathlib.Path(".")
container, anno_path = sys.argv[1], sys.argv[2]
glue = (root/"gui/pkg/itpp_gui.js").read_text()
wasm_b64 = base64.b64encode((root/"gui/pkg/itpp_gui_bg.wasm").read_bytes()).decode()
cont_b64 = base64.b64encode(pathlib.Path(container).read_bytes()).decode()
anno_json = pathlib.Path(anno_path).read_text().strip() if pathlib.Path(anno_path).exists() else "[]"
html = (root/"gui/template.html").read_text()
html = html.replace("/*__GLUE__*/", glue)
html = html.replace("__WASM_B64__", wasm_b64)
html = html.replace("__CONTAINER_B64__", cont_b64)
html = html.replace("__ANNO_JSON__", anno_json)
out = root/"gui/genome-browser.html"
out.write_text(html)
print(f"   wrote {out} ({out.stat().st_size//1024} KB; wasm {len(wasm_b64)//1024}KB b64, container {len(cont_b64)//1024}KB b64, {anno_json.count('{')} annotations)")
PY

echo ">> done. open gui/genome-browser.html in a browser (double-click)."
