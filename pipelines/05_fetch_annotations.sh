#!/usr/bin/env bash
# Fetch a dbSNP/ClinVar overlay track for a genomic window (GRCh38) from the Ensembl REST API,
# keep the clinically-significant (ClinVar-classified) variants, and write the small JSON the
# browser overlays: [{pos, id, sig, cons}] in backbone coordinates.
#
#   bash pipelines/05_fetch_annotations.sh 6 31972046 32055647 database/annotations/c4-clinvar.json
set -euo pipefail
cd "$(dirname "$0")/.."
CHR="${1:-6}"; START="${2:-31972046}"; END="${3:-32055647}"; OUT="${4:-database/annotations/c4-clinvar.json}"
mkdir -p "$(dirname "$OUT")"
tmp=$(mktemp -d)

# Ensembl caps region size, so fetch in ~28 kb chunks.
i=0
for (( s=START; s<END; s+=28000 )); do
  e=$(( s+28000 < END ? s+28000 : END ))
  curl -s --max-time 90 "https://rest.ensembl.org/overlap/region/human/${CHR}:${s}-${e}?feature=variation;content-type=application/json" -o "$tmp/c$i.json"
  i=$((i+1))
done

node -e '
const fs=require("fs"),dir=process.argv[1],out=process.argv[2];
let all=[];
for(const f of fs.readdirSync(dir)){ try{ all=all.concat(JSON.parse(fs.readFileSync(dir+"/"+f))); }catch(e){} }
const keep=all.filter(v=>Array.isArray(v.clinical_significance)&&v.clinical_significance.length
  && !(v.clinical_significance.length===1 && /^(benign|likely benign)$/i.test(v.clinical_significance[0])));
const anno=keep.map(v=>({pos:v.start,id:v.id,sig:v.clinical_significance.join("/"),cons:(v.consequence_type||"").replace(/_/g," ")}))
  .sort((a,b)=>a.pos-b.pos);
fs.writeFileSync(out, JSON.stringify(anno));
console.log(`${all.length} variants -> ${anno.length} clinically-significant -> ${out}`);
' "$tmp" "$OUT"
rm -rf "$tmp"
echo "rebuild browser:  bash gui/build.sh database/mhc-c4.itpp $OUT"
