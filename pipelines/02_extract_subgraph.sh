#!/usr/bin/env bash
# Extract the MHC subgraph from the whole-genome minigraph-cactus GFA and emit data/graphs/mhc.gfa
# ready for `itpp import`. Needs vg + odgi on PATH (see scripts/setup-env.sh).
set -euo pipefail
cd "$(dirname "$0")/.."

for t in vg odgi; do
  command -v "$t" >/dev/null || { echo "!! '$t' not on PATH — run scripts/setup-env.sh" >&2; exit 1; }
done

GFA_GZ="$(ls data/raw/*mc*.gfa.gz 2>/dev/null | head -1 || true)"
[[ -n "$GFA_GZ" ]] || { echo "!! no minigraph-cactus GFA in data/raw — run 01_fetch_mhc.sh" >&2; exit 1; }
REGION="$(grep -P '^mhc_region\t' data/manifests/mhc.tsv | cut -f2)"   # e.g. chr6:28510120-33480577
mkdir -p data/graphs

echo ">> building vg + snarl indexes"
vg convert -g "$GFA_GZ" -p > data/graphs/whole.pg
vg index -x data/graphs/whole.xg data/graphs/whole.pg

echo ">> chunking region $REGION"
# Path name for GRCh38 chr6 in the HPRC graph is typically GRCh38#0#chr6.
PATHNAME="${MHC_PATH:-GRCh38#0#chr6}"
vg chunk -x data/graphs/whole.xg -p "${PATHNAME}:${REGION#chr6:}" -c 20 > data/graphs/mhc.vg

echo ">> exporting GFA with walks (W-lines)"
vg convert -f -W data/graphs/mhc.vg > data/graphs/mhc.gfa

echo ">> snarl decomposition (the sites each haplotype's walk chooses between)"
vg snarls data/graphs/mhc.vg > data/graphs/mhc.snarls || true
vg deconstruct -a -e -P "$PATHNAME" data/graphs/mhc.vg > data/graphs/mhc.vcf || true

echo "done -> data/graphs/mhc.gfa"
echo "next: cargo run --release -p itpp-cli -- import --gfa data/graphs/mhc.gfa --out data/mhc.itpp"
