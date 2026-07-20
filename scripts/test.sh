#!/usr/bin/env bash
# Run the whole test suite: native Rust (all crates) + the WASM engine end-to-end via node.
set -euo pipefail
cd "$(dirname "$0")/.."

echo "== cargo test (native, all crates) =="
cargo test --workspace --quiet

echo "== clippy =="
cargo clippy --workspace --all-targets --quiet

echo "== wasm engine end-to-end (node) =="
if [ ! -d gui/pkg-node ]; then
  echo ".. building wasm bindings first"; bash gui/build.sh >/dev/null
fi
node gui/tests/smoke.mjs

echo "== ALL TESTS PASSED =="
