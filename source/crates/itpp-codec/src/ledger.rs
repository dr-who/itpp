//! The description-length ledger: per-section achieved bits by each coder, the winner, and
//! the container total → **bits/char**. This is the machine-readable measurement appended to
//! `results/metrics/ledger.jsonl`.
//!
//! The denominator (`total_bases`) is the cohort's described content — every haplotype's full
//! length. Because shared segments and repeat families are stored once, the numerator (total
//! bits) is far smaller than storing each haplotype independently; the ratio is the payoff of
//! the pangenome, and it should fall as we factor out more duplication.

/// One coder's achieved cost on a section (or as a whole-cohort baseline).
#[derive(Debug, Clone)]
pub struct CoderResult {
    pub coder: String,
    pub bits: u64,
}

/// Achieved bits for one section of the container, across every coder tried.
#[derive(Debug, Clone)]
pub struct SectionLedger {
    pub name: String,
    /// Bases attributable to this section (for per-layer bits/char attribution).
    pub bases: u64,
    pub results: Vec<CoderResult>,
}

impl SectionLedger {
    #[must_use]
    pub fn new(name: impl Into<String>, bases: u64) -> Self {
        SectionLedger { name: name.into(), bases, results: Vec::new() }
    }

    pub fn add(&mut self, coder: impl Into<String>, bits: u64) {
        self.results.push(CoderResult { coder: coder.into(), bits });
    }

    /// The winning (minimum-bits) coder for this section.
    #[must_use]
    pub fn best(&self) -> Option<&CoderResult> {
        self.results.iter().min_by_key(|r| r.bits)
    }
}

/// The whole-container ledger.
#[derive(Debug, Clone, Default)]
pub struct Ledger {
    pub sections: Vec<SectionLedger>,
    /// Cohort bases described — the bits/char denominator.
    pub total_bases: u64,
    /// Whole-cohort reference points (2-bit, order-0 Shannon, per-genome xz, …).
    pub baselines: Vec<CoderResult>,
    /// Free-form provenance (commit, dataset, tool versions).
    pub provenance: Vec<(String, String)>,
}

impl Ledger {
    #[must_use]
    pub fn new(total_bases: u64) -> Self {
        Ledger { total_bases, ..Default::default() }
    }

    pub fn push_section(&mut self, s: SectionLedger) {
        self.sections.push(s);
    }

    pub fn add_baseline(&mut self, coder: impl Into<String>, bits: u64) {
        self.baselines.push(CoderResult { coder: coder.into(), bits });
    }

    pub fn add_provenance(&mut self, key: impl Into<String>, val: impl Into<String>) {
        self.provenance.push((key.into(), val.into()));
    }

    /// Sum of the winning coder over every section — the container total.
    #[must_use]
    pub fn total_bits(&self) -> u64 {
        self.sections.iter().filter_map(|s| s.best().map(|r| r.bits)).sum()
    }

    /// The north-star metric.
    #[must_use]
    pub fn bits_per_char(&self) -> f64 {
        if self.total_bases == 0 {
            return 0.0;
        }
        self.total_bits() as f64 / self.total_bases as f64
    }

    /// Machine-readable one-line JSON record for `results/metrics/ledger.jsonl`.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::from("{");
        s.push_str(&format!("\"total_bases\":{},", self.total_bases));
        s.push_str(&format!("\"total_bits\":{},", self.total_bits()));
        s.push_str(&format!("\"bits_per_char\":{:.6},", self.bits_per_char()));
        if let Some(r) = self.baseline_bits("2bit_percopy") {
            if self.total_bits() > 0 {
                s.push_str(&format!(
                    "\"reduction_vs_2bit\":{:.4},",
                    r as f64 / self.total_bits() as f64
                ));
            }
        }
        if let Some(n) = self.baseline_bits("concat_nograph") {
            if self.total_bits() > 0 {
                s.push_str(&format!(
                    "\"reduction_vs_nograph\":{:.4},",
                    n as f64 / self.total_bits() as f64
                ));
            }
        }
        s.push_str("\"sections\":[");
        for (i, sec) in self.sections.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let best = sec.best();
            s.push_str(&format!(
                "{{\"name\":{},\"bases\":{},\"best_coder\":{},\"best_bits\":{},\"coders\":{{",
                json_str(&sec.name),
                sec.bases,
                json_str(best.map_or("", |b| b.coder.as_str())),
                best.map_or(0, |b| b.bits),
            ));
            for (j, r) in sec.results.iter().enumerate() {
                if j > 0 {
                    s.push(',');
                }
                s.push_str(&format!("{}:{}", json_str(&r.coder), r.bits));
            }
            s.push_str("}}");
        }
        s.push_str("],\"baselines\":{");
        for (i, b) in self.baselines.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str(&format!("{}:{}", json_str(&b.coder), b.bits));
        }
        s.push_str("},\"provenance\":{");
        for (i, (k, v)) in self.provenance.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str(&format!("{}:{}", json_str(k), json_str(v)));
        }
        s.push_str("}}");
        s
    }

    /// A whole-cohort baseline whose name contains `needle`.
    #[must_use]
    pub fn baseline_bits(&self, needle: &str) -> Option<u64> {
        self.baselines.iter().find(|b| b.coder.contains(needle)).map(|b| b.bits)
    }

    /// Human-readable table for the terminal.
    #[must_use]
    pub fn to_report(&self) -> String {
        let mut out = String::new();

        // Headline: absolute TOTAL INFORMATION CONTENT vs how the same data is stored today.
        let itpp = self.total_bits();
        out.push_str("═══ TOTAL INFORMATION CONTENT ═════════════════════════\n");
        out.push_str(&format!("  bases described          : {:>16}\n", commas(self.total_bases)));
        if let Some(r) = self.baseline_bits("2bit_percopy") {
            out.push_str(&format!(
                "  reference  (2 bit/base)  : {:>16} bits  ({} bytes)\n",
                commas(r),
                commas(r / 8)
            ));
        }
        if let Some(n) = self.baseline_bits("concat_nograph") {
            out.push_str(&format!(
                "  no-graph PPM (best)      : {:>16} bits  ({} bytes)\n",
                commas(n),
                commas(n / 8)
            ));
        }
        out.push_str(&format!(
            "  ITPP total (Σ components): {:>16} bits  ({} bytes)   ← the measure\n",
            commas(itpp),
            commas(itpp / 8)
        ));
        if let Some(r) = self.baseline_bits("2bit_percopy") {
            if itpp > 0 {
                out.push_str(&format!(
                    "  → {:.2}× smaller than the 2-bit reference",
                    r as f64 / itpp as f64
                ));
            }
        }
        if let Some(n) = self.baseline_bits("concat_nograph") {
            if itpp > 0 {
                out.push_str(&format!(", {:.2}× smaller than no-graph PPM", n as f64 / itpp as f64));
            }
        }
        out.push_str(&format!("\n  = {:.4} bits/char\n", self.bits_per_char()));
        out.push_str("═══════════════════════════════════════════════════════\n\n");

        out.push_str("── per-component entropy ──────────────────────────────\n");
        out.push_str(&format!(
            "{:<22} {:>14} {:>14} {:>10}\n",
            "section", "bases", "bits", "bits/char"
        ));
        for sec in &self.sections {
            let bits = sec.best().map_or(0, |b| b.bits);
            let bpc = if sec.bases > 0 { bits as f64 / sec.bases as f64 } else { 0.0 };
            out.push_str(&format!(
                "{:<22} {:>14} {:>14} {:>10.4}\n",
                sec.name, sec.bases, bits, bpc
            ));
        }
        out.push_str(&"-".repeat(62));
        out.push('\n');
        out.push_str(&format!(
            "{:<22} {:>14} {:>14} {:>10.4}\n",
            "TOTAL",
            self.total_bases,
            self.total_bits(),
            self.bits_per_char()
        ));
        if !self.baselines.is_empty() {
            out.push_str("\nbaselines (whole cohort):\n");
            for b in &self.baselines {
                let bpc = if self.total_bases > 0 {
                    b.bits as f64 / self.total_bases as f64
                } else {
                    0.0
                };
                out.push_str(&format!("  {:<20} {:>10.4} bits/char\n", b.coder, bpc));
            }
        }
        out
    }
}

/// Group an integer with thousands separators for the human report.
fn commas(n: u64) -> String {
    let s = n.to_string();
    let bytes = s.as_bytes();
    let mut out = String::new();
    for (i, b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(*b as char);
    }
    out
}

fn json_str(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            _ => out.push(c),
        }
    }
    out.push('"');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn totals_and_bpc() {
        let mut l = Ledger::new(1000);
        let mut s = SectionLedger::new("backbone", 500);
        s.add("order-8", 800);
        s.add("2bit", 1000);
        l.push_section(s);
        let mut s2 = SectionLedger::new("walks", 500);
        s2.add("order-4", 200);
        l.push_section(s2);
        assert_eq!(l.total_bits(), 1000); // 800 (best) + 200
        assert!((l.bits_per_char() - 1.0).abs() < 1e-9);
        // JSON is well-formed enough to contain the headline number.
        assert!(l.to_json().contains("\"bits_per_char\":1.0"));
    }
}
