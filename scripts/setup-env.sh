#!/usr/bin/env bash
# Create the conda/mamba environment holding the upstream *producers* ITPP ingests.
# These tools are NOT reimplemented — they generate the GFA / annotations we feed to `itpp`.
#
# Usage:  bash scripts/setup-env.sh        # creates env `itpp`
#         conda activate itpp
set -euo pipefail

ENV_NAME="${ITPP_ENV:-itpp}"

if ! command -v mamba >/dev/null 2>&1 && ! command -v conda >/dev/null 2>&1; then
  echo "Need conda or mamba (install Miniforge: https://github.com/conda-forge/miniforge)" >&2
  exit 1
fi
SOLVER="$(command -v mamba || command -v conda)"

echo ">> creating '$ENV_NAME' with the pangenome producer toolchain"
"$SOLVER" create -y -n "$ENV_NAME" -c bioconda -c conda-forge \
  vg \
  odgi \
  gfatools \
  samtools \
  bcftools \
  minimap2 \
  repeatmasker \
  kraken2 \
  ont-modkit

cat <<EOF

Done. Next:
  conda activate $ENV_NAME
  bash pipelines/01_fetch_mhc.sh
  bash pipelines/02_extract_subgraph.sh
  bash pipelines/03_annotate.sh
  cargo run --release -p itpp-cli -- import --gfa data/graphs/mhc.gfa --out data/mhc.itpp
  cargo run --release -p itpp-cli -- measure --in data/mhc.itpp

Notes:
  * RepeatMasker needs a Dfam library (the bioconda package ships a minimal one; for the full
    human repeat set install Dfam and run its configuration).
  * Contamination for M1 uses Kraken2 (small DB). FCS-GX is more thorough but its DB is ~470 GB.
  * Methylation (modkit) needs an ONT sample with MM/ML tags; see pipelines/03_annotate.sh.
EOF
