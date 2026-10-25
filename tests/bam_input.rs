//! Every passing SAM fixture converted to BAM (with samtools) must give the
//! same report.

mod common;
use common::*;

#[test]
fn bam_matches_sam_for_every_fixture() {
    if std::process::Command::new("samtools")
        .arg("--version")
        .output()
        .is_err()
    {
        eprintln!("samtools not found; skipping");
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    for name in [
        "flags_mixed.sam",
        "dup_softclip_reverse.sam",
        "fragments.sam",
        "suppl_secondary.sam",
        "optical_chain.sam",
        "insert.sam",
        "coverage.sam",
        "bad_names.sam",
    ] {
        let bam = dir.path().join(name.replace(".sam", ".bam"));
        let st = std::process::Command::new("samtools")
            .args(["view", "-b", "-o"])
            .arg(&bam)
            .arg(fixture(name))
            .status()
            .unwrap();
        assert!(st.success(), "{name}");
        let (mut sam, sd) = run_cfg(&config(name));
        let mut c = config(name);
        c.input = Some(bam);
        let (mut b, bd) = run_cfg(&c);
        for r in [&mut sam, &mut b] {
            r["inputs"]["input"] = serde_json::Value::Null;
        }
        assert_eq!(sam, b, "{name}");
        assert_eq!(sd, bd, "{name}");
    }
}
