#!/usr/bin/env bash
# Build a REAL full ~5 Mb MHC pangenome graph (with per-haplotype paths) from the classic
# GRCh38 MHC reference haplotypes, using a single static `vg` binary — no conda, fits in a
# few hundred MB of disk. Output: database/graphs/mhc.gfa -> `itpp import`.
#
# This is the approved one-time external fetch. After it runs, the resulting container is
# committed to our project and everything else loads from our project only.
#
# Disk: needs ~1 GB free (vg binary + 9 haplotype FASTAs + working graph). Check `df -h .`.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p database/raw/mhc database/graphs tools

VG=tools/vg
if [[ ! -x "$VG" ]]; then
  echo ">> fetching static vg binary"
  curl -L --fail -o "$VG" "https://github.com/vgteam/vg/releases/latest/download/vg"
  chmod +x "$VG"
fi
"$VG" version

# The 8 GRCh38 alternate MHC haplotypes (APD/COX/DBB/MANN/MCF/QBL/SSTO) as alt scaffolds,
# plus the primary MHC window (PGF) sliced from chr6. All from UCSC hg38.
UCSC="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes"
ALTS=(chr6_GL000250v2_alt chr6_GL000251v2_alt chr6_GL000252v2_alt chr6_GL000253v2_alt \
      chr6_GL000254v2_alt chr6_GL000255v2_alt chr6_GL000256v2_alt)
for a in "${ALTS[@]}"; do
  f="database/raw/mhc/$a.fa.gz"
  [[ -f "$f" ]] || { echo ">> $a"; curl -L --fail -o "$f" "$UCSC/$a.fa.gz"; }
done
# primary MHC window (GRCh38 chr6:28,510,120-33,480,577)
if [[ ! -f database/raw/mhc/chr6_MHC_PGF.fa ]]; then
  echo ">> chr6 primary MHC window"
  curl -L --fail -o database/raw/mhc/chr6.fa.gz "$UCSC/chr6.fa.gz"
  "$VG" 2>/dev/null || true
  # slice with samtools if present, else awk fallback is left to the operator
  if command -v samtools >/dev/null; then
    gunzip -kf database/raw/mhc/chr6.fa.gz
    samtools faidx database/raw/mhc/chr6.fa chr6:28510120-33480577 \
      | sed '1s/.*/>chr6_MHC_PGF/' > database/raw/mhc/chr6_MHC_PGF.fa
  else
    echo "!! need samtools to slice the primary MHC window; skipping PGF (alts still build a graph)"
  fi
fi

echo ">> assembling all MHC haplotypes"
: > database/raw/mhc/all.fa
for f in database/raw/mhc/*.fa.gz; do gunzip -c "$f" >> database/raw/mhc/all.fa; done
[[ -f database/raw/mhc/chr6_MHC_PGF.fa ]] && cat database/raw/mhc/chr6_MHC_PGF.fa >> database/raw/mhc/all.fa

echo ">> vg msga: aligning haplotypes into a graph (this is the slow step)"
"$VG" msga -f database/raw/mhc/all.fa -t "$(nproc)" -b chr6_MHC_PGF > database/graphs/mhc.vg
"$VG" mod -U 10 database/graphs/mhc.vg > database/graphs/mhc.mod.vg   # normalize
"$VG" view database/graphs/mhc.mod.vg > database/graphs/mhc.gfa

echo ">> import + verify + measure with itpp"
cargo build --release -q
./target/release/itpp import --gfa database/graphs/mhc.gfa --out database/mhc.itpp
./target/release/itpp verify --in database/mhc.itpp --gfa database/graphs/mhc.gfa
./target/release/itpp report --in database/mhc.itpp --dataset grch38-mhc-full

echo ">> rebuild the browser against the full MHC"
bash gui/build.sh database/mhc.itpp
echo "done. commit database/mhc.itpp + the refreshed gui/genome-browser.html."
