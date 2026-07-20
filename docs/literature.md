# Literature & tools

## HPRC (humanpangenome.org) — the anchor
Builds a graph reference from ~90 telomere-to-telomere haplotype assemblies. Notably built
the same input set three ways, spanning the complexity/losslessness tradeoff:

- **Minigraph** — sparse SV-level backbone (~392K nodes / ~566K edges, 3.2 Gbp). Our "backbone".
- **Minigraph-Cactus** — base-level, near-lossless, practical default (~80M nodes / ~111M edges).
  Paper: Nat. Biotechnol. 2023, <https://www.nature.com/articles/s41587-023-01793-w>
- **PGGB** (wfmash → seqwish → smoothxg) — fully lossless, reference-free (~111M nodes).

## Graph ecosystem
- **vg** — variation graph toolkit. `vg deconstruct` decomposes into **snarls/bubbles** = the
  sites where a haplotype's walk diverges. `vcfbub` filters nested sites.
- **odgi** — graph analysis; a delta-encoded (id-delta) storage that gets **cheaper per
  haplotype past ~32** — direct empirical support for the total-bits objective.
  <https://academic.oup.com/bioinformatics/article/38/13/3319/6585331>
- **gfatools**, **samtools/bcftools**, **minimap2**.
- Snarl/superbubble decomposition in linear time (SPQR-tree framework), arXiv:2511.21919.
- Grammar-based compression of GFA paths, bioRxiv 2025.05.22.655470 — path (W/P) lines dominate
  GFA size; grammar/dictionary coding of walks is exactly our per-sample residual problem.

## Repeats / ERV (the dictionary layer)
- **RepeatMasker** + **Dfam** library + **RepeatModeler** + **ULTRA** — the T2T-CHM13 repeat
  pipeline (Alu, L1/LINE, SINE, satellite HORs; HERV/ERV families are Dfam-covered).
  T2T repeat/epigenetic state: Science 2022, doi:10.1126/science.abk3112.
- Higher-order repeat (HOR) detection in T2T-CHM13, PMC12385485.

## Contamination
- NCBI **FCS-GX** (thorough, DB ~470 GB) and **Kraken2** (light, small DB — used for M1).

## Epigenetics
- ONT modified-base calling → **modkit** → bedMethyl. A separate information channel, not
  sequence entropy; gets its own meter.

## Compression / information theory (the meter)
- **PPM** / order-k context models — our internal range coder.
- **NAF** (reference-free DNA archival), Bioinformatics 2019 — a strong reference-free baseline.
- **GeCo3** — context-mixing DNA compressor (neural mixing of context + substitution-tolerant
  models); a target to beat on the residual.
- General-purpose `xz` (LZMA), `zstd`, `bzip2` (BWT) — cross-check baselines.
