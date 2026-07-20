#!/usr/bin/env bash
# Run the producer annotations over the MHC region: repeats/ERV (RepeatMasker+Dfam),
# contamination (Kraken2), methylation (modkit). Outputs land in data/anno/ as BED/TSV that a
# future `itpp annotate` step folds into the container's DICT / CTAM / EPIG sections.
#
# M1 status: this emits the annotation files; the `itpp annotate` importer that attaches them
# to a container is the next increment. Each block no-ops with a message if its tool is absent.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p data/anno
REF="data/graphs/mhc_backbone.fa"   # produced by slicing the backbone; see 02_extract_subgraph.sh notes

# --- repeats + ERV ----------------------------------------------------------
if command -v RepeatMasker >/dev/null 2>&1 && [[ -f "$REF" ]]; then
  echo ">> RepeatMasker (Alu/LINE/SINE/Satellite + HERV/ERV via Dfam)"
  RepeatMasker -species human -pa 4 -dir data/anno "$REF"
  # data/anno/<ref>.out is the annotation table -> dictionary instances
else
  echo ".. skip RepeatMasker (tool or $REF missing) — provides DICT (repeat/ERV) layer"
fi

# --- contamination ----------------------------------------------------------
KDB="$(ls -d data/raw/k2_* 2>/dev/null | head -1 || true)"
if command -v kraken2 >/dev/null 2>&1 && [[ -n "$KDB" && -f "$REF" ]]; then
  echo ">> Kraken2 contamination screen"
  kraken2 --db "$KDB" --output data/anno/kraken2.out --report data/anno/kraken2.report "$REF"
else
  echo ".. skip Kraken2 (tool/DB/$REF missing) — provides CTAM (contamination) sidecar"
fi

# --- methylation ------------------------------------------------------------
ONT_BAM="${ONT_BAM:-data/raw/mhc.ont.bam}"   # an ONT alignment carrying MM/ML modified-base tags
if command -v modkit >/dev/null 2>&1 && [[ -f "$ONT_BAM" ]]; then
  echo ">> modkit pileup (methylation channel)"
  modkit pileup "$ONT_BAM" data/anno/methyl.bed --ref "$REF"
else
  echo ".. skip modkit (tool or ONT BAM missing) — provides EPIG (methylation) channel"
  echo "   sourcing matched ONT methylation for MHC is the known M1 data risk (see plan)."
fi

echo "annotations (where produced) in data/anno/"
