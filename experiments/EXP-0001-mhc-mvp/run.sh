#!/usr/bin/env bash
# EXP-0001 — MHC milestone, reproducible end-to-end.
#
# Part A (runs anywhere, no genomics toolchain): a synthetic MHC-like pangenome exercises the
# whole meter and appends a bits/char record to results/metrics/ledger.jsonl.
# Part B (commented) is the identical flow on real HPRC MHC data once the toolchain is set up.
set -euo pipefail
cd "$(dirname "$0")/../.."   # repo root

cargo build --release
BIN=./target/release/itpp
WORK=experiments/EXP-0001-mhc-mvp

echo "===== Part A: synthetic MHC-like pangenome ====="
"$BIN" synth --out "$WORK/mhc.itpp" --gfa "$WORK/mhc.gfa" --haplotypes 10 --blocks 80 --seed 7
"$BIN" stats  --in "$WORK/mhc.itpp"
"$BIN" verify --in "$WORK/mhc.itpp" --gfa "$WORK/mhc.gfa"
"$BIN" report --in "$WORK/mhc.itpp" --dataset synth-mhc-seed7

echo
echo "ledger now has $(wc -l < results/metrics/ledger.jsonl) record(s)."
echo "watch bits/char fall across commits with:  jq .bits_per_char results/metrics/ledger.jsonl"

# ===== Part B: real HPRC MHC (needs scripts/setup-env.sh) =====
# conda activate itpp
# bash pipelines/01_fetch_mhc.sh
# bash pipelines/02_extract_subgraph.sh
# bash pipelines/03_annotate.sh
# "$BIN" import --gfa data/graphs/mhc.gfa --out data/mhc.itpp
# "$BIN" verify --in data/mhc.itpp --gfa data/graphs/mhc.gfa
# "$BIN" report --in data/mhc.itpp --dataset hprc-mhc-v1.1
