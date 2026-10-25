//! MultiQC custom content. MultiQC reads one section per custom-content
//! file, so the general-stats row and the two line graphs go in three
//! files: `<prefix>.general_mqc.json`, `<prefix>.insert_mqc.json` and
//! `<prefix>.coverage_mqc.json`. Percentages are shown as percentages
//! here (MultiQC convention), unlike the 0-1 fractions in the JSON report.

use serde_json::{Value, json};

use crate::report::json::Report;

/// Cumulative coverage is drawn up to this depth (or the maximum depth).
pub const CUMULATIVE_MAX_DEPTH: usize = 250;

fn pct(x: f64) -> f64 {
    100.0 * x
}

pub fn general(r: &Report) -> Value {
    let f = &r.coverage.filtered;
    let mut row = serde_json::Map::new();
    row.insert("mean_coverage".into(), json!(f.mean));
    if let Some(v) = f.pct_ge.get("30") {
        row.insert("pct_ge_30x".into(), json!(pct(*v)));
    }
    row.insert("duplicate_rate".into(), json!(pct(r.duplicates.rate)));
    row.insert("median_insert".into(), json!(r.insert_size.fr.median));
    if let Some(p) = r.uniformity.pmad_gc {
        row.insert("pmad_gc".into(), json!(p));
    }
    row.insert("qc_status".into(), json!(r.qc_status.overall));
    json!({
        "id": "warpqc_general",
        "section_name": "WarpQC",
        "plot_type": "generalstats",
        "pconfig": [
            {"mean_coverage": {"title": "Mean cov", "description": "WarpQC filtered mean depth (MAPQ, BQ, duplicates, overlap)", "suffix": "X", "format": "{:,.1f}", "min": 0}},
            {"pct_ge_30x": {"title": "≥ 30X", "description": "WarpQC: % of territory at filtered depth >= 30", "suffix": "%", "min": 0, "max": 100}},
            {"duplicate_rate": {"title": "Dups", "description": "WarpQC duplicate rate (Picard PERCENT_DUPLICATION)", "suffix": "%", "min": 0, "max": 100}},
            {"median_insert": {"title": "Insert", "description": "WarpQC median insert size, FR pairs", "suffix": " bp", "format": "{:,.0f}"}},
            {"pmad_gc": {"title": "PMAD", "description": "WarpQC GC-corrected PMAD (DRAGEN definition)", "format": "{:,.3f}"}},
            {"qc_status": {"title": "QC", "description": "WarpQC qc_status.overall"}}
        ],
        "data": { r.sample.clone(): row },
    })
}

pub fn insert_size(r: &Report) -> Value {
    let points: Vec<[u64; 2]> = r
        .insert_size
        .fr
        .histogram
        .iter()
        .enumerate()
        .filter(|(_, c)| **c > 0)
        .map(|(i, &c)| [i as u64, c])
        .collect();
    json!({
        "id": "warpqc_insert_size",
        "section_name": "WarpQC insert size",
        "description": "FR pairs, non-duplicate, from R2 |TLEN| (up to 2,000 bp).",
        "plot_type": "linegraph",
        "pconfig": {"id": "warpqc_insert_size_plot", "title": "WarpQC: insert size", "xlab": "Insert size (bp)", "ylab": "Pairs"},
        "data": { r.sample.clone(): points },
    })
}

/// % of territory at filtered depth >= x, for x = 0..=min(max depth, 250).
pub fn cumulative(hist: &[u64]) -> Vec<(u64, f64)> {
    let total: u64 = hist.iter().sum();
    if total == 0 {
        return vec![];
    }
    let mut at_least = total;
    let mut out = Vec::new();
    for (d, &c) in hist.iter().enumerate().take(CUMULATIVE_MAX_DEPTH + 1) {
        out.push((d as u64, 100.0 * at_least as f64 / total as f64));
        at_least -= c;
    }
    out
}

pub fn coverage(r: &Report) -> Value {
    let points: Vec<(u64, f64)> = cumulative(&r.coverage.filtered.histogram);
    json!({
        "id": "warpqc_cumulative_coverage",
        "section_name": "WarpQC cumulative coverage",
        "description": "Filtered depth over non-N analysis positions.",
        "plot_type": "linegraph",
        "pconfig": {"id": "warpqc_cumulative_coverage_plot", "title": "WarpQC: cumulative coverage", "xlab": "Depth (X)", "ylab": "% bases ≥ depth", "ymin": 0, "ymax": 100},
        "data": { r.sample.clone(): points },
    })
}

/// (suffix, content) for each file.
pub fn files(r: &Report) -> [(&'static str, Value); 3] {
    [
        (".general_mqc.json", general(r)),
        (".insert_mqc.json", insert_size(r)),
        (".coverage_mqc.json", coverage(r)),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cumulative_is_percent_at_or_above() {
        // 2 positions at 0, 1 at 1, 1 at 3.
        let c = cumulative(&[2, 1, 0, 1]);
        assert_eq!(c, vec![(0, 100.0), (1, 50.0), (2, 25.0), (3, 25.0)]);
        assert!(cumulative(&[]).is_empty());
    }

    #[test]
    fn general_stats_row() {
        let mut r = Report {
            sample: "S1".into(),
            ..Default::default()
        };
        r.coverage.filtered.mean = 31.5;
        r.coverage.filtered.pct_ge.insert("30".into(), 0.9);
        r.duplicates.rate = 0.08;
        r.uniformity.pmad_gc = Some(0.21);
        let g = general(&r);
        assert_eq!(g["plot_type"], "generalstats");
        let row = &g["data"]["S1"];
        assert_eq!(row["mean_coverage"], 31.5);
        assert_eq!(row["pct_ge_30x"], 90.0);
        assert_eq!(row["duplicate_rate"], 8.0);
        assert_eq!(row["pmad_gc"], 0.21);
    }
}
