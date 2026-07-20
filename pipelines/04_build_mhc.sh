#!/usr/bin/env bash
# Build a REAL full ~5 Mb MHC pangenome graph with minigraph (robust, single small tool) from
# the GRCh38 alternate MHC haplotypes, then import with itpp. minigraph emits rGFA (no P/W
# lines); our GFA importer reconstructs the backbone spine from the rank-0 (SR:i:0) segments.
#
# Approved one-time external fetch. Output: database/graphs/mhc.rgfa -> database/mhc.itpp.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p tools database/raw/mhc database/graphs

# --- minigraph (compile from source; tiny, needs only zlib + a C compiler) ---
if [[ ! -x tools/minigraph ]]; then
  echo ">> building minigraph"
  curl -L --fail -o /tmp/minigraph.tar.gz "https://github.com/lh3/minigraph/archive/refs/heads/master.tar.gz"
  tar -xzf /tmp/minigraph.tar.gz -C /tmp
  make -C /tmp/minigraph-master -j"$(nproc)"
  cp /tmp/minigraph-master/minigraph tools/minigraph
fi
tools/minigraph --version

# --- the 7 GRCh38 alternate MHC haplotypes (full ~4.7-5 Mb each), from UCSC hg38 ---
UCSC="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes"
ALTS=(chr6_GL000251v2_alt chr6_GL000252v2_alt chr6_GL000253v2_alt chr6_GL000254v2_alt \
      chr6_GL000255v2_alt chr6_GL000256v2_alt chr6_GL000250v2_alt)
for a in "${ALTS[@]}"; do
  f="database/raw/mhc/$a.fa.gz"
  [[ -f "$f" ]] || { echo ">> $a"; curl -L --fail -o "$f" "$UCSC/$a.fa.gz"; }
  gunzip -kf "$f"
done

# --- build the graph (first haplotype = reference spine) ---
echo ">> minigraph -cxggs (building the MHC graph)"
tools/minigraph -cxggs -t"$(nproc)" \
  database/raw/mhc/chr6_GL000251v2_alt.fa \
  database/raw/mhc/chr6_GL000252v2_alt.fa database/raw/mhc/chr6_GL000253v2_alt.fa \
  database/raw/mhc/chr6_GL000254v2_alt.fa database/raw/mhc/chr6_GL000255v2_alt.fa \
  database/raw/mhc/chr6_GL000256v2_alt.fa database/raw/mhc/chr6_GL000250v2_alt.fa \
  > database/graphs/mhc.rgfa
echo "   rGFA: $(grep -c '^S' database/graphs/mhc.rgfa) segments, $(grep -c '^L' database/graphs/mhc.rgfa) links"

# --- import + measure with itpp ---
cargo build --release -q
./target/release/itpp import --gfa database/graphs/mhc.rgfa --out database/mhc.itpp
./target/release/itpp stats  --in database/mhc.itpp
./target/release/itpp verify --in database/mhc.itpp
./target/release/itpp report --in database/mhc.itpp --dataset grch38-mhc-minigraph
echo ">> rebuild the browser against the full MHC:  bash gui/build.sh database/mhc.itpp"
