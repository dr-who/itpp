// End-to-end WASM test: load the real MHC-C4 container in Node via the wasm-bindgen nodejs
// binding, run a query, and assert the JSON is well-formed with real hits + coordinates.
//
//   node gui/tests/smoke.mjs
//
// Requires `gui/pkg-node` (built by gui/build.sh) and `database/mhc-c4.itpp`.
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";
import { createRequire } from "node:module";

const here = dirname(fileURLToPath(import.meta.url));
const root = join(here, "..", "..");
const require = createRequire(import.meta.url);

let failures = 0;
const ok = (cond, msg) => {
  if (cond) {
    console.log(`  ok   ${msg}`);
  } else {
    console.error(`  FAIL ${msg}`);
    failures++;
  }
};

const { GenomeBrowser } = require(join(root, "gui", "pkg-node", "itpp_gui.js"));
const bytes = readFileSync(join(root, "database", "mhc-c4.itpp"));

const b = new GenomeBrowser(bytes);
ok(b.segments === 1748, `segments = ${b.segments} (expect 1748)`);
ok(b.haplotypes === 90, `haplotypes = ${b.haplotypes} (expect 90)`);
ok(/chr6/.test(b.region), `region mentions chr6: ${b.region}`);

// A query that must exist in real human sequence.
const json = JSON.parse(b.query("ATG", 20, 3));
ok(json.n_hits > 0, `ATG has ${json.n_hits} hits`);
ok(json.hits.length > 0, "hits array populated");

const h = json.hits[0];
ok(h.nodes.some((n) => n.hit && n.x === 0), "hit node laid out at column x=0");
ok(h.nodes.length >= 1, `local subgraph has ${h.nodes.length} nodes`);
ok(Array.isArray(h.edges), "edges present");
ok(typeof h.protein === "string" && h.protein.length > 0, `protein: ${h.protein}`);
ok(h.codons.length > 0, `codons: ${h.codons.slice(0, 3).join(",")}…`);
ok(/chr6:|variant/.test(h.pos), `positional info: ${h.pos}`);
ok(h.n_haplotypes >= 1, `${h.n_haplotypes} haplotypes through the hit`);

// Empty query → no hits, still valid JSON.
const empty = JSON.parse(b.query("", 5, 2));
ok(empty.n_hits === 0, "empty query yields 0 hits");

// A longer real query still round-trips through the layout.
const longer = JSON.parse(b.query("GAATTC", 10, 4)); // EcoRI site, common
ok(Array.isArray(longer.hits), "longer query returns hits array");

// Overview: the whole graph laid out in genomic coordinates.
const ov = JSON.parse(b.overview());
ok(ov.n_nodes === 1748, `overview has ${ov.n_nodes} nodes`);
ok(ov.chrom === "chr6" && ov.start > 31_000_000, `overview coords ${ov.chrom}:${ov.start}`);
ok(ov.span > 80000, `span ${ov.span} bp`);
ok(ov.nodes.some((n) => n.bb) && ov.nodes.some((n) => !n.bb), "overview has backbone + variant nodes");
ok(ov.nodes.every((n) => typeof n.x === "number" && typeof n.lane === "number"), "every node has x + lane");
ok(ov.nodes.every((n) => typeof n.kind === "string"), "every node has a type/kind for colouring");
ok(new Set(ov.nodes.map((n) => n.kind)).size >= 2, "at least two node types present");
ok(Array.isArray(ov.edges) && ov.edges.length > 0, `overview has ${ov.edges.length} edges`);

// Matches: node ids to highlight/zoom-to.
const mm = JSON.parse(b.matches("CAGCAGCAGCAG", 500));
ok(mm.matches.length >= 1 && mm.matches[0].node !== undefined, "matches return node ids");
ok(mm.n_total >= mm.matches.length, "match total ≥ highlighted");

// Context (extend) around a real node returns a hit object.
const ctx = JSON.parse(b.context(ov.nodes[0].id, 0, "+", 6, 6));
ok(ctx.hit && Array.isArray(ctx.hit.nodes), "context returns a hit with nodes");

console.log(failures === 0 ? "\nALL WASM SMOKE TESTS PASSED" : `\n${failures} FAILURES`);
process.exit(failures === 0 ? 0 : 1);
