#!/usr/bin/env bash
# Fetch the pinned MHC inputs listed in database/manifests/mhc.tsv into database/raw/.
# Idempotent: skips files already present. Requires curl.
set -euo pipefail
cd "$(dirname "$0")/.."

MANIFEST="database/manifests/mhc.tsv"
mkdir -p database/raw

fetch() { # key -> database/raw/<basename>
  local key="$1"
  local url
  url="$(grep -P "^${key}\t" "$MANIFEST" | cut -f2)"
  if [[ -z "$url" ]]; then echo "!! key '$key' not in $MANIFEST" >&2; return 1; fi
  case "$url" in
    chr6:*) echo ".. $key is a coordinate, not a download ($url)"; return 0;;
  esac
  local out="database/raw/$(basename "$url")"
  if [[ -f "$out" ]]; then echo "== have $out"; return 0; fi
  echo ">> $key -> $out"
  curl -L --fail --retry 3 -o "$out" "$url"
}

for key in hprc_mc_gfa grch38_ref kraken2_db dfam_min; do
  fetch "$key" || echo "!! fetch $key failed (continue)"
done

echo "done. MHC region = $(grep -P '^mhc_region\t' "$MANIFEST" | cut -f2)"
