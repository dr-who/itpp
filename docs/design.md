# ITPP design: the two-part code and the bits/char meter

## The objective

We model a cohort of human genomes as a **graph plus per-individual walks**, and we
drive down the **total description length** of the whole representation. The metric is:

```
bits/char = (achieved compressed size of the whole container, in bits) / (total bases described)
```

We report it as **two absolute numbers plus their ratio**: what a chromosome/cohort costs to
encode today (the 2-bit reference, and the best no-graph compressor) versus the **total
information content** — the summed entropy of every ITPP component. e.g. *"chr6 MHC: 733,400
bits today → 96,208 bits (Σ components) = 7.6× smaller"*. The total is the measure; bits/char is
that total divided by bases.

"Achieved compressed size" is **real**: it is `compressed_bytes * 8` from actual coders
(a PPM / order-k context model — our internal range coder — plus LZ- and BWT-based
methods, and general-purpose `xz`/`zstd`/`bzip2` as cross-checks). Order-0 Shannon
entropy is kept **only** as a cheap lower-effort sanity baseline; it is not the number we
report as the result.

### Why *total* bits is the right objective

Minimizing the total compressed size is equivalent to **removing replicated / duplicated
components**. If two individuals share a megabase, a good representation stores it once; if
a million Alu copies derive from one ~300 bp consensus, a good representation stores the
consensus once plus a small per-copy edit script. Every reduction in total bits corresponds
to a redundancy we successfully factored out.

A crucial and counter-intuitive consequence: as we factor better, the **local per-character
entropy of the residual goes up** (what's left is closer to incompressible novelty), while
the **total** goes down. So we track total bits and bits/char, not the entropy of any one
stream — a rising residual entropy next to a falling total is the signature of success.

## The two-part code

```
D_total = L(model) + L(data | model)

L(model)      = backbone spine
              + graph topology (nodes/edges beyond the backbone)
              + dictionaries (repeat consensi: Alu, L1/LINE, SINE, satellite HORs; ERV/HERV)

L(data|model) = Σ_individuals [ walk through the graph (allele choices at snarls)
                              + private novel sequence (deltas) ]

separate channels (own meters, not mixed into the human-sequence bits/char):
              + epigenetics  (per-position methylation)
              + contamination sidecar (non-human spans, excluded from the human model)
```

`bits/char = D_total / Σ bases`. Logged per commit to `results/metrics/ledger.jsonl`;
`itpp report` renders the curve.

## How each part is coded (v0)

| Part                | Coder (v0)                                   | Notes |
|---------------------|----------------------------------------------|-------|
| Backbone spine      | internal order-k range coder (PPM-style)     | anchors node coordinate space |
| Node/segment seqs   | internal order-k range coder                 | only sequence *not* attributable to backbone or a dictionary |
| Repeat/ERV instances| `(dict_id, pos, edit-script)`, edits range-coded | consensus stored once in the dictionary |
| Walks (per sample)  | allele-choice stream, range-coded            | a haplotype = a walk = list of oriented node ids |
| Private deltas      | internal order-k range coder                 | novel sequence not in the graph |
| Cross-checks        | `xz`, `zstd`, `bzip2` (external, optional)   | sanity vs the internal coder; best-of is reported |

The **ledger** (`itpp-codec::Ledger`) records, per section: bytes described, achieved bits
by each available coder, the winning coder, and the resulting bits/char. The container total
is the sum; the winning per-section coder id is stored so `verify` can reconstruct.

## Correctness contract

The representation is **lossless**: `itpp verify` reconstructs every input haplotype from the
container and asserts byte-exact equality. Any mismatch is a hard failure (non-zero exit).
Coordinate mapping between GFA node space and the backbone is the main hazard and is what
`verify` guards.

## Baselines we compare against

- `2.0` bits/base — the naive 2-bit packing of ACGT.
- order-0 Shannon of the concatenated cohort — redundancy-blind reference point.
- per-genome `xz`/`zstd` of each individual independently — "no shared model" reference.
- the graph model total — should beat all of the above; the gap is the value of the pangenome.

## Roadmap of the metric

Each milestone adds a factoring layer and we expect bits/char to drop:
1. backbone + walks only (shared nodes already remove cross-individual duplication)
2. + repeat dictionary (Alu/LINE/SINE/satellite)
3. + ERV/HERV library
4. + contamination removed from the denominator/model
5. + better coders (higher-order context, mixing) squeezing the residual
