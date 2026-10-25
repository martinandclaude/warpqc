mod common;
use common::*;
use std::process::Command;
use warpqc::report::json::Report;

fn thresholds(name: &str) -> std::path::PathBuf {
    root().join("tests/fixtures/thresholds").join(name)
}

#[test]
fn qc_status_from_thresholds() {
    let mut c = config("coverage.sam");
    c.thresholds = Some(thresholds("qc.yaml"));
    let (r, _) = run_cfg(&c);
    assert_eq!(get(&r, "/qc_status/overall").as_str(), Some("fail"));
    let checks = get(&r, "/qc_status/checks").as_array().unwrap();
    let st: Vec<&str> = checks
        .iter()
        .map(|c| c["status"].as_str().unwrap())
        .collect();
    assert_eq!(st, ["fail", "warn", "pass", "missing"]);
    close(&r, "/qc_status/checks/0/value", 103.0 / 780.0);
}

#[test]
fn typo_in_thresholds_fails_before_reading_input() {
    let dir = tempfile::tempdir().unwrap();
    // stdin is never written: the run must fail on the thresholds alone.
    let o = Command::new(env!("CARGO_BIN_EXE_warpqc"))
        .args([
            "run",
            "--sample",
            "S",
            "--contigs",
            "chr1,chr2,chrM",
            "--reference",
        ])
        .arg(reference())
        .arg("--thresholds")
        .arg(thresholds("typo.yaml"))
        .arg("--out")
        .arg(dir.path().join("S"))
        .arg(fixture("coverage.sam"))
        .output()
        .unwrap();
    assert_eq!(o.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&o.stderr).contains("meen"));
    assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
}

#[test]
fn report_round_trips_through_json() {
    let mut c = config("coverage.sam");
    c.thresholds = Some(thresholds("qc.yaml"));
    let out = warpqc::pipeline::run(&c).unwrap();
    let back = Report::from_json(&out.report.to_json()).unwrap();
    assert_eq!(back, out.report);
}

#[test]
fn reads_schema_1_0_reports_and_rejects_2_x() {
    let mut v: serde_json::Value = serde_json::to_value(
        warpqc::pipeline::run(&config("coverage.sam"))
            .unwrap()
            .report,
    )
    .unwrap();
    v.as_object_mut().unwrap().remove("uniformity");
    v["schema_version"] = "1.0".into();
    assert!(Report::from_json(&serde_json::to_vec(&v).unwrap()).is_ok());
    v["schema_version"] = "2.0".into();
    assert_eq!(
        Report::from_json(&serde_json::to_vec(&v).unwrap())
            .unwrap_err()
            .exit_code(),
        2
    );
}
