//! `itpp` — the command-line interface.
//!
//!   itpp synth    --out FILE.itpp [--gfa FILE.gfa] [--haplotypes N] [--blocks N] [--seed S]
//!   itpp import   --gfa FILE.gfa --out FILE.itpp
//!   itpp measure  --in FILE.itpp [--json FILE] [--no-external]
//!   itpp verify   --in FILE.itpp [--gfa FILE.gfa]
//!   itpp stats    --in FILE.itpp
//!   itpp report   --in FILE.itpp [--ledger PATH] [--dataset NAME] [--no-external]

mod measure;

use itpp_core::Graph;
use measure::{measure as run_measure, MeasureOpts};
use std::process::ExitCode;

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("{USAGE}");
        return ExitCode::FAILURE;
    }
    let rest = &args[2..];
    let result = match args[1].as_str() {
        "synth" => cmd_synth(rest),
        "import" => cmd_import(rest),
        "measure" => cmd_measure(rest),
        "verify" => cmd_verify(rest),
        "stats" => cmd_stats(rest),
        "report" => cmd_report(rest),
        "-h" | "--help" | "help" => {
            println!("{USAGE}");
            Ok(())
        }
        other => Err(format!("unknown command '{other}'\n\n{USAGE}")),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

const USAGE: &str = "\
itpp — Information-Theoretic Pangenomic Project

USAGE:
  itpp synth   --out FILE.itpp [--gfa FILE.gfa] [--haplotypes N] [--blocks N]
               [--block-len N] [--seed S]
  itpp import  --gfa FILE.gfa --out FILE.itpp
  itpp measure --in FILE.itpp [--json FILE] [--no-external]
  itpp verify  --in FILE.itpp [--gfa FILE.gfa]
  itpp stats   --in FILE.itpp
  itpp report  --in FILE.itpp [--ledger PATH] [--dataset NAME] [--no-external]";

// ---- arg helpers ------------------------------------------------------------

fn flag<'a>(args: &'a [String], name: &str) -> Option<&'a str> {
    args.iter().position(|a| a == name).and_then(|i| args.get(i + 1)).map(String::as_str)
}

fn has_flag(args: &[String], name: &str) -> bool {
    args.iter().any(|a| a == name)
}

fn required<'a>(args: &'a [String], name: &str) -> Result<&'a str, String> {
    flag(args, name).ok_or_else(|| format!("missing required {name}"))
}

fn parse_or<T: std::str::FromStr>(args: &[String], name: &str, default: T) -> T {
    flag(args, name).and_then(|s| s.parse().ok()).unwrap_or(default)
}

fn load_container(path: &str) -> Result<Graph, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("reading {path}: {e}"))?;
    itpp_format::read_container(&bytes).map_err(|e| format!("parsing {path}: {e}"))
}

// ---- commands ---------------------------------------------------------------

fn cmd_synth(args: &[String]) -> Result<(), String> {
    let out = required(args, "--out")?;
    let params = itpp_ingest::synth::SynthParams {
        haplotypes: parse_or(args, "--haplotypes", 8),
        backbone_blocks: parse_or(args, "--blocks", 60),
        block_len: parse_or(args, "--block-len", 400),
        seed: parse_or(args, "--seed", 42),
        ..Default::default()
    };
    let g = itpp_ingest::synth::generate(&params);
    std::fs::write(out, itpp_format::write_container(&g))
        .map_err(|e| format!("writing {out}: {e}"))?;
    println!(
        "synthesized pangenome: {} segments, {} walks ({} haplotypes + backbone), {} dict families, {} repeat instances",
        g.segments.len(),
        g.walks.len(),
        params.haplotypes,
        g.dictionary.entries.len(),
        g.instances.len()
    );
    println!("wrote container: {out}");
    if let Some(gfa_path) = flag(args, "--gfa") {
        std::fs::write(gfa_path, itpp_ingest::gfa::write(&g))
            .map_err(|e| format!("writing {gfa_path}: {e}"))?;
        println!("wrote GFA: {gfa_path}");
    }
    Ok(())
}

fn cmd_import(args: &[String]) -> Result<(), String> {
    let gfa = required(args, "--gfa")?;
    let out = required(args, "--out")?;
    let text = std::fs::read_to_string(gfa).map_err(|e| format!("reading {gfa}: {e}"))?;
    let g = itpp_ingest::gfa::parse(&text);
    if g.segments.is_empty() {
        return Err(format!("{gfa} produced no segments — is it GFA?"));
    }
    std::fs::write(out, itpp_format::write_container(&g))
        .map_err(|e| format!("writing {out}: {e}"))?;
    println!(
        "imported {gfa}: {} segments, {} walks, backbone = {}",
        g.segments.len(),
        g.walks.len(),
        g.backbone.as_deref().unwrap_or("(none)")
    );
    println!("wrote container: {out}");
    Ok(())
}

fn measure_opts(args: &[String]) -> MeasureOpts {
    MeasureOpts { external: !has_flag(args, "--no-external"), ..Default::default() }
}

fn cmd_measure(args: &[String]) -> Result<(), String> {
    let input = required(args, "--in")?;
    let g = load_container(input)?;
    let led = run_measure(&g, &measure_opts(args));
    print!("{}", led.to_report());
    if let Some(json_path) = flag(args, "--json") {
        std::fs::write(json_path, led.to_json())
            .map_err(|e| format!("writing {json_path}: {e}"))?;
        println!("\nwrote metrics JSON: {json_path}");
    }
    Ok(())
}

fn cmd_verify(args: &[String]) -> Result<(), String> {
    let input = required(args, "--in")?;
    let bytes = std::fs::read(input).map_err(|e| format!("reading {input}: {e}"))?;
    let g = itpp_format::read_container(&bytes).map_err(|e| format!("parsing {input}: {e}"))?;

    // 1. container is self-consistent (re-serialize → re-parse is identical)
    let g2 = itpp_format::read_container(&itpp_format::write_container(&g))
        .map_err(|e| format!("re-parse: {e}"))?;
    if g != g2 {
        return Err("container failed self-consistency round-trip".into());
    }

    // 2. every walk spells without a dangling node reference
    let mut spelled_bases = 0u64;
    for w in &g.walks {
        match g.spell(w) {
            Some(s) => spelled_bases += s.len() as u64,
            None => return Err(format!("walk '{}' references a missing segment", w.name)),
        }
    }

    // 3. if a source GFA is given, spelled sequences must match byte-for-byte
    let mut checked = 0usize;
    if let Some(gfa) = flag(args, "--gfa") {
        let text = std::fs::read_to_string(gfa).map_err(|e| format!("reading {gfa}: {e}"))?;
        let src = itpp_ingest::gfa::parse(&text);
        for w in &g.walks {
            if let Some(sw) = src.walk_by_name(&w.name) {
                let a = g.spell(w).ok_or("spell failed")?;
                let b = src.spell(sw).ok_or("source spell failed")?;
                if a != b {
                    return Err(format!("walk '{}' does not match source GFA", w.name));
                }
                checked += 1;
            }
        }
    }

    println!(
        "VERIFY OK — {} walks spelled ({} bases), self-consistent{}",
        g.walks.len(),
        spelled_bases,
        if checked > 0 {
            format!(", {checked} matched against source GFA")
        } else {
            String::new()
        }
    );
    Ok(())
}

fn cmd_stats(args: &[String]) -> Result<(), String> {
    let input = required(args, "--in")?;
    let g = load_container(input)?;
    println!("── ITPP container: {input} ──");
    println!("backbone         : {}", g.backbone.as_deref().unwrap_or("(none)"));
    println!("segments         : {}", g.segments.len());
    println!("  segment bases  : {}", g.segment_bases());
    println!("edges            : {}", g.edges.len());
    println!("walks            : {}", g.walks.len());
    println!("  cohort bases   : {}", g.walk_bases());
    println!("dictionary       : {} families", g.dictionary.entries.len());
    println!("repeat instances : {}", g.instances.len());
    println!("contamination    : {} spans", g.contamination.len());
    Ok(())
}

fn cmd_report(args: &[String]) -> Result<(), String> {
    let input = required(args, "--in")?;
    let ledger_path = flag(args, "--ledger").unwrap_or("results/metrics/ledger.jsonl");
    let dataset = flag(args, "--dataset").unwrap_or("unnamed");
    let g = load_container(input)?;
    let mut led = run_measure(&g, &measure_opts(args));
    led.add_provenance("dataset", dataset);
    led.add_provenance("input", input);
    if let Ok(t) = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH) {
        led.add_provenance("unixtime", t.as_secs().to_string());
    }

    print!("{}", led.to_report());

    if let Some(parent) = std::path::Path::new(ledger_path).parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|e| format!("mkdir {parent:?}: {e}"))?;
        }
    }
    let mut line = led.to_json();
    line.push('\n');
    use std::io::Write;
    let mut f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(ledger_path)
        .map_err(|e| format!("opening {ledger_path}: {e}"))?;
    f.write_all(line.as_bytes()).map_err(|e| format!("appending {ledger_path}: {e}"))?;
    println!("\nappended ledger record → {ledger_path}");
    Ok(())
}
