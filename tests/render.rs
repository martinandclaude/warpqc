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
        "--html".as_ref(),
        "--multiqc".as_ref(),
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
    warpqc::report::write_all(&out.report, &prefix, warpqc::report::Outputs::ALL).unwrap();
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

fn files_in(dir: &std::path::Path) -> Vec<String> {
    let mut v: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    v.sort();
    v
}

/// `run` writes only the JSON unless --html / --multiqc ask for more.
#[test]
fn run_writes_only_json_by_default() {
    for (flags, want) in [
        (vec![], vec!["S.qc.json"]),
        (vec!["--html"], vec!["S.qc.html", "S.qc.json"]),
        (
            vec!["--multiqc"],
            vec![
                "S.coverage_mqc.json",
                "S.general_mqc.json",
                "S.insert_mqc.json",
                "S.qc.json",
            ],
        ),
        (
            vec!["--html", "--multiqc"],
            vec![
                "S.coverage_mqc.json",
                "S.general_mqc.json",
                "S.insert_mqc.json",
                "S.qc.html",
                "S.qc.json",
            ],
        ),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("S");
        let reference = reference();
        let mut args: Vec<&std::ffi::OsStr> = vec![
            "run".as_ref(),
            "--sample".as_ref(),
            "S".as_ref(),
            "--contigs".as_ref(),
            "chr1,chr2,chrM".as_ref(),
            "--reference".as_ref(),
            reference.as_os_str(),
            "--out".as_ref(),
            out.as_os_str(),
        ];
        args.extend(flags.iter().map(std::ffi::OsStr::new));
        let input = fixture("insert.sam");
        args.push(input.as_os_str());
        let o = warpqc(&args);
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        assert_eq!(files_in(dir.path()), want, "flags {flags:?}");
    }
}

/// `render` writes both by default; a flag narrows it to that output.
#[test]
fn render_flags_select_outputs() {
    let out = warpqc::pipeline::run(&config("insert.sam")).unwrap();
    for (flag, want) in [
        (
            None,
            vec![
                "S.coverage_mqc.json",
                "S.general_mqc.json",
                "S.insert_mqc.json",
                "S.qc.html",
                "S.qc.json",
            ],
        ),
        (Some("--html"), vec!["S.qc.html", "S.qc.json"]),
        (
            Some("--multiqc"),
            vec![
                "S.coverage_mqc.json",
                "S.general_mqc.json",
                "S.insert_mqc.json",
                "S.qc.json",
            ],
        ),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("S");
        warpqc::report::write_all(&out.report, &prefix, warpqc::report::Outputs::NONE).unwrap();
        let json = dir.path().join("S.qc.json");
        let mut args: Vec<&std::ffi::OsStr> =
            vec!["render".as_ref(), "--json".as_ref(), json.as_os_str()];
        if let Some(f) = flag {
            args.push(f.as_ref());
        }
        let o = warpqc(&args);
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        assert_eq!(files_in(dir.path()), want, "flag {flag:?}");
    }
}
