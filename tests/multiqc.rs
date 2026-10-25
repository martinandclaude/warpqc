//! Runs the real MultiQC on WarpQC output when it is available (set
//! WARPQC_MULTIQC to its path, or have `multiqc` on PATH; CI installs it).

mod common;
use common::*;
use std::process::Command;

/// With WARPQC_MULTIQC set (as in CI) MultiQC must work; otherwise the test
/// skips when `multiqc` is not on PATH.
fn multiqc() -> Option<String> {
    let set = std::env::var("WARPQC_MULTIQC").ok();
    let exe = set.clone().unwrap_or_else(|| "multiqc".into());
    let ok = Command::new(&exe)
        .arg("--version")
        .output()
        .ok()
        .is_some_and(|o| o.status.success());
    if set.is_some() && !ok {
        panic!("WARPQC_MULTIQC={exe} does not run");
    }
    ok.then_some(exe)
}

#[test]
fn multiqc_reads_all_three_files() {
    let Some(mq) = multiqc() else {
        eprintln!("multiqc not found; skipping");
        return;
    };
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in");
    for (sample, fx) in [("S01", "insert.sam"), ("S02", "coverage.sam")] {
        let mut c = config(fx);
        c.sample = sample.into();
        let out = warpqc::pipeline::run(&c).unwrap();
        warpqc::report::write_all(
            &out.report,
            &input.join(sample),
            warpqc::report::Outputs::ALL,
        )
        .unwrap();
    }
    let outdir = dir.path().join("out");
    let st = Command::new(&mq)
        .args(["-q", "-f", "-o"])
        .arg(&outdir)
        .arg(&input)
        .status()
        .unwrap();
    assert!(st.success());
    let data = outdir.join("multiqc_data");
    let gs = std::fs::read_to_string(data.join("multiqc_general_stats.txt")).unwrap();
    let header = gs.lines().next().unwrap();
    for col in [
        "mean_coverage",
        "pct_ge_30x",
        "duplicate_rate",
        "median_insert",
        "qc_status",
    ] {
        assert!(header.contains(col), "{col} missing from {header}");
    }
    assert!(
        gs.lines().any(|l| l.starts_with("S01\t")) && gs.lines().any(|l| l.starts_with("S02\t")),
        "{gs}"
    );
    for plot in [
        "multiqc_warpqc_insert_size_plot.txt",
        "multiqc_warpqc_cumulative_coverage_plot.txt",
    ] {
        let t = std::fs::read_to_string(data.join(plot)).unwrap();
        assert_eq!(
            t.lines().filter(|l| l.starts_with("S0")).count(),
            2,
            "{plot}: {t}"
        );
    }
}
