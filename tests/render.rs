mod common;
use common::*;
use std::process::Command;

fn warpqc(args: &[&std::ffi::OsStr]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_warpqc"))
        .args(args)
        .output()
        .unwrap()
}

#[test]
fn run_then_render_gives_identical_files() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a/S");
    let o = warpqc(&[
        "run".as_ref(),
        "--sample".as_ref(),
        "S".as_ref(),
        "--contigs".as_ref(),
        "chr1,chr2,chrM".as_ref(),
        "--bin-size".as_ref(),
        "20".as_ref(),
        "--reference".as_ref(),
        reference().as_os_str(),
        "--out".as_ref(),
        a.as_os_str(),
        fixture("coverage.sam").as_os_str(),
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let b = dir.path().join("b/S");
    std::fs::create_dir_all(b.parent().unwrap()).unwrap();
    let json = dir.path().join("a/S.qc.json");
    let o = warpqc(&[
        "render".as_ref(),
        "--json".as_ref(),
        json.as_os_str(),
        "--out".as_ref(),
        b.as_os_str(),
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    for suffix in [
        ".qc.html",
        ".general_mqc.json",
        ".insert_mqc.json",
        ".coverage_mqc.json",
    ] {
        let x = std::fs::read(dir.path().join(format!("a/S{suffix}"))).unwrap();
        let y = std::fs::read(dir.path().join(format!("b/S{suffix}"))).unwrap();
        assert_eq!(x, y, "{suffix} differs after render");
    }
    // Without thresholds the JSON is not rewritten.
    assert!(!dir.path().join("b/S.qc.json").exists());
}

#[test]
fn render_with_thresholds_rewrites_qc_status() {
    let dir = tempfile::tempdir().unwrap();
    let out = warpqc::pipeline::run(&config("coverage.sam")).unwrap();
    let prefix = dir.path().join("S");
    warpqc::report::write_all(&out.report, &prefix).unwrap();
    assert_eq!(out.report.qc_status.overall, "not_evaluated");
    let json = dir.path().join("S.qc.json");
    let qc = root().join("tests/fixtures/thresholds/qc.yaml");
    // Spec form: --out S.qc.html is accepted as the prefix S.
    let html = dir.path().join("S.qc.html");
    let o = warpqc(&[
        "render".as_ref(),
        "--json".as_ref(),
        json.as_os_str(),
        "--thresholds".as_ref(),
        qc.as_os_str(),
        "--out".as_ref(),
        html.as_os_str(),
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let r: serde_json::Value = serde_json::from_slice(&std::fs::read(&json).unwrap()).unwrap();
    assert_eq!(get(&r, "/qc_status/overall").as_str(), Some("fail"));
    assert!(
        std::fs::read_to_string(&html)
            .unwrap()
            .contains(r#"<span class="badge fail">fail</span>"#)
    );
    // A typo in the thresholds leaves everything as it was.
    let typo = root().join("tests/fixtures/thresholds/typo.yaml");
    let before = std::fs::read(&json).unwrap();
    let o = warpqc(&[
        "render".as_ref(),
        "--json".as_ref(),
        json.as_os_str(),
        "--thresholds".as_ref(),
        typo.as_os_str(),
    ]);
    assert_eq!(o.status.code(), Some(2));
    assert_eq!(std::fs::read(&json).unwrap(), before);
}

#[test]
fn render_prefix_rules() {
    use std::path::Path;
    use warpqc::report::render_prefix;
    assert_eq!(
        render_prefix(Path::new("x/S.qc.json"), None),
        Path::new("x/S")
    );
    assert_eq!(
        render_prefix(Path::new("x/S.qc.json"), Some(Path::new("y/T.qc.html"))),
        Path::new("y/T")
    );
    assert_eq!(
        render_prefix(Path::new("x/S.qc.json"), Some(Path::new("y/T"))),
        Path::new("y/T")
    );
}
