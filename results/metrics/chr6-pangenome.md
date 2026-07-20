# chr6 full-chromosome compression (HPRC v1.1-mc-grch38)

Extracted from the HPRC whole-genome minigraph-cactus graph (streamed from exaba S3).

- ~90 haplotypes, GRCh38 chr6 backbone = 170.8 Mb
- graph stores 175.8 Mb distinct sequence in 4,517,930 segments
- variation of 90 haplotypes adds only **5.0 Mb (+2.9%)** over one reference chromosome
- naive independent storage: 90 × 170.8 Mb = **15.37 Gb** (3.84 GB @2-bit)
- pangenome graph: **175.8 Mb** (43.9 MB @2-bit)
- **structural compression ≈ 87×**, bits/char ≈ **0.023** (2-bit; entropy coding lowers further)

Method: streamed the extracted chr6.gfa.gz, summed S-line sequence lengths (total + SR:i:0
backbone). The full entropy-coded `itpp measure` (real bits/char, backbone < 2 bits/base) is the
next step once the codec is scaled from Mb to 100s of Mb.
