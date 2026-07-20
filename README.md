# ITPP — Information-Theoretic Pangenomic Project

Model the human genome as a **graph plus per-individual walks**, and drive down its
**total description length**. The single number we care about:

```
bits/char = (real compressed size of the whole representation, in bits) / (bases described)
```

measured with actual coders (a PPM/order-k context model — our internal range coder — plus
LZ/BWT methods and `xz`/`zstd`/`bzip2` cross-checks), logged per commit so we can watch it
fall as the model improves. Minimizing total bits is the same thing as **removing duplicated
components** (shared sequence across individuals, repeat families, ERVs): the only way the
total drops is if redundancy was factored out.

See [`docs/design.md`](docs/design.md) for the two-part-code cost model and
[`spec/itpp-format-v0.md`](spec/itpp-format-v0.md) for the container format.

## Layout

```
crates/
  itpp-core     types: sequence, graph, walks, dictionaries, samples, container model
  itpp-codec    range coder + adaptive order-k model + multi-coder ledger
  itpp-format   the ITPP container: chunked, indexed, versioned, lossless round-trip
  itpp-ingest   GFA / FASTA importers + synthetic MHC-like fixture generator
  itpp-cli      the `itpp` binary
spec/           the format specification (versioned)
docs/           design + literature
pipelines/      thin wrappers around upstream producers (vg, RepeatMasker, Kraken2, modkit)
data/manifests/ pinned data URLs (downloads gitignored)
experiments/    reproducible runs; EXP-0001 = MHC MVP
results/metrics/ append-only bits/char ledger over commits
```

## Quick start (self-contained, no genomics toolchain needed)

```sh
cargo build --release
# generate a synthetic MHC-like pangenome (backbone + haplotype walks + repeat insertions)
cargo run --release -p itpp-cli -- synth --out /tmp/mhc.gfa --haplotypes 8 --seed 42
# import GFA -> ITPP container
cargo run --release -p itpp-cli -- import --gfa /tmp/mhc.gfa --out /tmp/mhc.itpp
# measure bits/char with per-layer breakdown vs baselines
cargo run --release -p itpp-cli -- measure --in /tmp/mhc.itpp
# prove losslessness (reconstruct every walk, byte-exact)
cargo run --release -p itpp-cli -- verify --in /tmp/mhc.itpp --gfa /tmp/mhc.gfa
```

Real HPRC MHC data uses the same commands on a GFA produced by
[`pipelines/`](pipelines/) — see [`experiments/EXP-0001-mhc-mvp`](experiments/EXP-0001-mhc-mvp).

## Status

Milestone 1 (MHC): format v0 + Rust reference implementation with working
`synth / import / measure / verify`. Upstream genomics tools are *data producers* we ingest,
not reimplemented. Python/C bindings come later on top of the `itpp-core`/`-codec`/`-format` API.
