# ITPP container format — v0 (draft)

Status: **draft / unstable**. Little-endian. All offsets are byte offsets from file start.

The container holds a two-part code: a **shared model** (backbone, graph, dictionaries) and
**per-sample residuals** (walks + private deltas), plus separate **contamination** and
**epigenetics** channels. It is **lossless**: every input haplotype is reconstructable
byte-exact. See [`../docs/design.md`](../docs/design.md) for the cost model.

## File skeleton

```
+-----------------+  offset 0
| FileHeader      |  16 bytes
+-----------------+
| Section chunks  |  variable, in any order
|  ...            |
+-----------------+
| Index chunk     |  (tag = "INDX") lists every section: tag, offset, len
+-----------------+
| Footer          |  16 bytes; points at the Index chunk
+-----------------+  EOF
```

### FileHeader (16 bytes)
| field        | type    | value |
|--------------|---------|-------|
| magic        | `[u8;4]`| `ITPP` |
| version      | `u16`   | `0` |
| flags        | `u16`   | bit0 = has-contamination, bit1 = has-epigenetics |
| reserved     | `u64`   | `0` |

### Chunk framing
Every section is a chunk:
```
tag: [u8;4]      # section kind, see below
len: u64         # payload length in bytes
payload: [u8;len]
```

### Footer (16 bytes, last in file)
| field        | type    | value |
|--------------|---------|-------|
| index_offset | `u64`   | byte offset of the `INDX` chunk |
| magic        | `[u8;4]`| `ITP1` |
| reserved     | `u32`   | `0` |

Readers seek to `EOF-16`, validate `ITP1`, read `index_offset`, parse `INDX`, then load
sections by tag. This makes the container seekable without scanning.

## Sections

Every section payload begins with a `u16` **encoding id** identifying how the bytes were
coded (so a reader/`verify` can invert it), followed by section-specific fields:

| enc id | meaning                                   |
|--------|-------------------------------------------|
| 0      | raw (uncompressed)                        |
| 1      | 2-bit packed ACGT (+ exceptions list)     |
| 2      | order-k adaptive range coder (PPM-style)  |
| 3      | external: xz / zstd / bzip2 (id in header)|

### `HEAD` — manifest
`enc=0`. Fields: format params, `k` (context order), sample count, sample names (length-
prefixed UTF-8), coordinate-system tag, and free-form provenance (producer tool versions,
source URLs, git commit) as length-prefixed key/value pairs.

### `BONE` — backbone spine
The reference path's sequence, order-k range-coded (`enc=2`). Header carries the original
base count and the model order `k`. Defines the coordinate anchor for node ids.

### `GRAF` — graph topology
Nodes and edges beyond the backbone. Per node: id (varint delta), sequence *source*
(backbone-span | dictionary-instance | literal), and literal sequence (range-coded) when
present. Edges: `(from, from_orient, to, to_orient)` as delta-coded varints.

### `DICT` — dictionaries
Repeat + ERV consensi. Per entry: family id, class (Alu/L1/SINE/Satellite/ERV/other),
consensus sequence (range-coded). Instances live in `GRAF`/`SAMP` as
`(dict_id, orientation, edit-script)`.

### `SAMP` — samples (per-individual residual)
Per sample: an ordered **walk** = list of oriented node ids (delta + range-coded), plus
**private deltas** = novel sequence not present as any node (range-coded). Reconstructing a
sample = concatenating node sequences along the walk, applying deltas.

### `CTAM` — contamination sidecar (optional)
Non-human spans: `(sample, start, len, taxon-label)`. Excluded from the human-model bits/char
denominator but stored for losslessness.

### `EPIG` — epigenetics channel (optional)
Per sample, per position: methylation state/probability, coded independently. Reported with
its **own** meter; never mixed into the sequence bits/char.

### `LEDG` — ledger footer
`enc=0`, JSON-ish key/values: per-section {bytes described, bits by coder, winning coder},
container total bits, total bases, `bits_per_char`, and baselines (2-bit, order-0 Shannon,
per-genome xz). This is the machine-readable measurement `itpp report` appends to
`results/metrics/ledger.jsonl`.

## Versioning
`version` in the FileHeader gates breaking changes. New section tags are ignored by older
readers if not required (a required-sections bitmap may be added pre-1.0). v0 makes **no**
stability promises.
