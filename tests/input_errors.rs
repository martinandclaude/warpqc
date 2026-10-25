//! End-to-end through the binary: exit codes, and no report on failure.

mod common;
use common::*;
use std::process::{Command, Output, Stdio};

fn warpqc(args: &[&str], input: &str, stdin: Option<&[u8]>) -> (Output, tempfile::TempDir) {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("S");
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_warpqc"));
    cmd.args([
        "run",
        "--sample",
        "S",
        "--threads",
        "2",
        "--contigs",
        "chr1,chr2,chrM",
        "--reference",
    ])
    .arg(reference())
    .arg("--out")
    .arg(&out)
    .args(args)
    .arg(input)
    .stdin(if stdin.is_some() {
        Stdio::piped()
    } else {
        Stdio::null()
    })
    .stdout(Stdio::piped())
    .stderr(Stdio::piped());
    let mut child = cmd.spawn().unwrap();
    if let Some(bytes) = stdin {
        use std::io::Write;
        child.stdin.take().unwrap().write_all(bytes).unwrap();
    }
    (child.wait_with_output().unwrap(), dir)
}

fn assert_input_error(name: &str, needle: &str) {
    let (o, dir) = warpqc(&[], fixture(name).to_str().unwrap(), None);
    let err = String::from_utf8_lossy(&o.stderr);
    assert_eq!(o.status.code(), Some(2), "{name}: {err}");
    assert!(err.contains(needle), "{name}: {err}");
    assert_eq!(
        std::fs::read_dir(dir.path()).unwrap().count(),
        0,
        "{name}: output written"
    );
}

#[test]
fn coordinate_sorted_header() {
    assert_input_error("so_coordinate.sam", "coordinate-sorted");
}

#[test]
fn reappearing_read_name() {
    assert_input_error("reappearing_name.sam", "reappears");
}

#[test]
fn sq_mismatch() {
    assert_input_error("sq_mismatch.sam", "length 301");
}

#[test]
fn incomplete_templates() {
    assert_input_error("incomplete_templates.sam", "missing a primary mate");
}

#[test]
fn truncated_stdin_and_expect_reads() {
    let full = std::fs::read(fixture("insert.sam")).unwrap();
    // Cut mid-record.
    let (o, _) = warpqc(&[], "-", Some(&full[..full.len() - 7]));
    assert_eq!(o.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&o.stderr).contains("truncated"));
    // Whole stream, wrong expected count.
    let (o, dir) = warpqc(&["--expect-reads", "12"], "-", Some(&full));
    assert_eq!(o.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&o.stderr).contains("--expect-reads"));
    assert_eq!(std::fs::read_dir(dir.path()).unwrap().count(), 0);
    // Right count: report written from stdin.
    let (o, dir) = warpqc(&["--expect-reads", "11"], "-", Some(&full));
    assert_eq!(
        o.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
    let r: serde_json::Value =
        serde_json::from_slice(&std::fs::read(dir.path().join("S.qc.json")).unwrap()).unwrap();
    assert_eq!(u(&r, "/alignment/total_templates"), 11);
    assert_eq!(get(&r, "/inputs/input").as_str(), Some("-"));
    assert_eq!(get(&r, "/tool/name").as_str(), Some("warpqc"));
    assert_eq!(
        get(&r, "/qc_status/overall").as_str(),
        Some("not_evaluated")
    );
}

#[test]
fn unparseable_names_warn_and_continue() {
    let (r, _) = run_fixture("bad_names.sam");
    let w = get(&r, "/warnings").as_array().unwrap();
    assert!(
        w.iter().any(|w| w
            .as_str()
            .unwrap()
            .contains("optical duplicate detection is disabled")),
        "{w:?}"
    );
    assert!(get(&r, "/duplicates/optical_duplicate_pairs").is_null());
    assert_eq!(u(&r, "/duplicates/duplicate_pairs"), 1);
    assert_eq!(get(&r, "/inputs/lanes"), &serde_json::json!(["unknown"]));
    assert_eq!(u(&r, "/per_lane/unknown/templates"), 2);
}

#[test]
fn reference_mismatch_via_cache() {
    let dir = tempfile::tempdir().unwrap();
    let other = dir.path().join("other.fa");
    std::fs::write(&other, ">chr1\nACGT\n").unwrap();
    std::fs::write(dir.path().join("other.fa.fai"), "chr1\t4\t6\t4\t5\n").unwrap();
    let cache = dir.path().join("other.warpqc");
    let st = Command::new(env!("CARGO_BIN_EXE_warpqc"))
        .args(["prepare-ref", "--reference"])
        .arg(&other)
        .arg("--out")
        .arg(&cache)
        .status()
        .unwrap();
    assert!(st.success());
    let (o, _) = warpqc(
        &["--ref-cache", cache.to_str().unwrap()],
        fixture("insert.sam").to_str().unwrap(),
        None,
    );
    assert_eq!(o.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&o.stderr).contains("different reference"));
}

#[test]
fn prepare_ref_cache_gives_identical_report() {
    let dir = tempfile::tempdir().unwrap();
    let cache = dir.path().join("mini.warpqc");
    let st = Command::new(env!("CARGO_BIN_EXE_warpqc"))
        .args(["prepare-ref", "--reference"])
        .arg(reference())
        .arg("--out")
        .arg(&cache)
        .status()
        .unwrap();
    assert!(st.success());
    let (a, _) = run_cfg(&config("coverage.sam"));
    let mut c = config("coverage.sam");
    c.ref_cache = Some(cache);
    let (b, _) = run_cfg(&c);
    assert_eq!(a, b);
}
