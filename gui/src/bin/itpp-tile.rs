//! `itpp-tile` — generate a multi-resolution tile pyramid from an ITPP container.
//!
//! Partitions the region by coordinate into levels (level 0 = whole span in one tile; each
//! level down doubles the resolution). Coarse levels are density aggregates (variant counts +
//! CNV-sized features); fine levels carry the graph nodes (with sequence at base resolution).
//! Each tile is a small JSON file, so a viewer streams only the tiles in view — the same data,
//! never loading more than the viewport. Point it at the S3 mount to store PB-scale pyramids.
//!
//!   itpp-tile --in FILE.itpp --out /mnt/exaba-itpp/tiles/mhc-c4
//!
//! Layout: <out>/manifest.json + <out>/L<level>/<tile>.json

use itpp_gui::browser::Browser;
use std::{env, fs};

const COARSE_ABOVE_BP: i64 = 50_000; // tiles wider than this are density aggregates
const SEQ_BELOW_BP: i64 = 8_000; // tiles narrower than this carry per-base sequence
const BASE_TILE_BP: i64 = 2_000; // finest tile width target
const MAX_NODES_PER_TILE: usize = 8_000;
const DENSITY_BINS: usize = 128;

fn flag<'a>(args: &'a [String], name: &str) -> Option<&'a str> {
    args.iter().position(|a| a == name).and_then(|i| args.get(i + 1)).map(String::as_str)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let input = flag(&args, "--in").expect("usage: itpp-tile --in FILE.itpp --out DIR");
    let outdir = flag(&args, "--out").expect("usage: itpp-tile --in FILE.itpp --out DIR");

    let bytes = fs::read(input).unwrap_or_else(|e| panic!("read {input}: {e}"));
    let graph = itpp_format::read_container(&bytes).unwrap_or_else(|e| panic!("parse {input}: {e}"));
    let b = Browser::new(graph);
    let span = b.span().max(1);

    // deepest level whose tiles are ~BASE_TILE_BP wide, capped so a chromosome stays ~<=32k
    // finest tiles (2^15) rather than exploding into hundreds of thousands of files
    const MAX_LEVEL_CAP: u32 = 15;
    let mut max_level = 0u32;
    while (span >> (max_level + 1)) > BASE_TILE_BP && max_level < MAX_LEVEL_CAP {
        max_level += 1;
    }

    // tile width per level (ceil), computed up front
    let tile_bp_per_level: Vec<i64> =
        (0..=max_level).map(|l| { let n: i64 = 1 << l; (span + n - 1) / n }).collect();

    // Write the manifest FIRST so a viewer can load and render the levels that are already
    // done while the rest stream in (progressive tiling, important for whole-genome runs).
    fs::create_dir_all(outdir).unwrap_or_else(|e| panic!("mkdir {outdir}: {e}"));
    let tile_bp_list =
        tile_bp_per_level.iter().map(|n| n.to_string()).collect::<Vec<_>>().join(",");
    let manifest = format!(
        "{{\"chrom\":\"{}\",\"start\":{},\"span\":{},\"max_level\":{},\"coarse_above_bp\":{},\"seq_below_bp\":{},\"tile_bp\":[{}]}}",
        b.chrom(), b.start(), span, max_level, COARSE_ABOVE_BP, SEQ_BELOW_BP, tile_bp_list
    );
    fs::write(format!("{outdir}/manifest.json"), manifest)
        .unwrap_or_else(|e| panic!("write manifest: {e}"));

    let mut total_tiles = 0u64;
    for level in 0..=max_level {
        let ntiles: i64 = 1 << level;
        let tile_bp = tile_bp_per_level[level as usize];
        let coarse = tile_bp > COARSE_ABOVE_BP;
        let leveldir = format!("{outdir}/L{level}");
        fs::create_dir_all(&leveldir).unwrap_or_else(|e| panic!("mkdir {leveldir}: {e}"));
        for t in 0..ntiles {
            let s0 = t * tile_bp;
            let s1 = (s0 + tile_bp).min(span);
            let json = if coarse {
                b.density_json(s0, s1, DENSITY_BINS)
            } else {
                b.window_json(s0, s1, tile_bp <= SEQ_BELOW_BP, MAX_NODES_PER_TILE)
            };
            fs::write(format!("{leveldir}/{t}.json"), json)
                .unwrap_or_else(|e| panic!("write tile L{level}/{t}: {e}"));
            total_tiles += 1;
        }
    }

    println!(
        "tiled {} → {outdir}: {} levels (0..{}), {} tiles, span {} bp on {}",
        input,
        max_level + 1,
        max_level,
        total_tiles,
        span,
        b.chrom()
    );
}
