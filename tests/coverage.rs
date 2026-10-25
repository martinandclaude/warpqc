mod common;
use common::*;

#[test]
fn coverage_positions() {
    let (_, d) = run_fixture("coverage.sam");
    // O: overlap counted once; BQ-failing base in one mate covered by the other.
    for p in 11..=25 {
        let want = if p == 24 { 0 } else { 1 };
        assert_eq!(depth_at(&d, "chr2", p), (want, Some(1)), "chr2:{p}");
    }
    // Y: both mates fail BQ at 105.
    for p in 101..=110 {
        assert_eq!(depth_at(&d, "chr2", p), (if p == 105 { 0 } else { 1 }, Some(1)), "chr2:{p}");
    }
    // Z: MAPQ 10 -> raw only.
    assert_eq!(depth_at(&d, "chr2", 205), (0, Some(1)));
    // D1/D2: duplicates in raw only.
    assert_eq!(depth_at(&d, "chr2", 155), (1, Some(2)));
    assert_eq!(depth_at(&d, "chr2", 185), (1, Some(2)));
    // Deletion is not coverage.
    assert_eq!(depth_at(&d, "chr1", 305), (1, Some(1)));
    assert_eq!(depth_at(&d, "chr1", 306), (0, Some(0)));
    assert_eq!(depth_at(&d, "chr1", 307), (0, Some(0)));
    assert_eq!(depth_at(&d, "chr1", 308), (1, Some(1)));
    // QC-fail pair contributes nowhere.
    assert_eq!(depth_at(&d, "chr1", 105), (0, Some(0)));
    // chrM unpaired and mate-unmapped reads.
    assert_eq!(depth_at(&d, "chrM", 1), (1, Some(1)));
    assert_eq!(depth_at(&d, "chrM", 55), (1, Some(1)));
    // The decoy is not an analysis contig.
    assert!(d.iter().all(|(n, _, _)| n != "chrUn_decoy"));
}

#[test]
fn coverage_summary() {
    let (r, _) = run_fixture("coverage.sam");
    assert_eq!(u(&r, "/coverage/filtered/territory"), 780);
    close(&r, "/coverage/filtered/mean", 103.0 / 780.0);
    close(&r, "/coverage/raw/mean", 145.0 / 780.0);
    close(&r, "/coverage/filtered/pct_ge/1", 103.0 / 780.0);
    close(&r, "/coverage/raw/pct_ge/1", 125.0 / 780.0);
    close(&r, "/coverage/filtered/median", 0.0);
    assert!(get(&r, "/coverage/filtered/fold80").is_null());
    let h = get(&r, "/coverage/filtered/histogram").as_array().unwrap();
    assert_eq!((h[0].as_u64(), h[1].as_u64(), h.len()), (Some(677), Some(103), 2));
    close(&r, "/coverage/filtered/per_contig/chr1", 40.0 / 380.0);
    close(&r, "/coverage/filtered/per_contig/chr2", 43.0 / 300.0);
    close(&r, "/coverage/filtered/per_contig/chrM", 20.0 / 100.0);
    close(&r, "/coverage/filtered/autosomal_mean", 83.0 / 680.0);
    let e = |k: &str| u(&r, &format!("/coverage/exclusions/bases/{k}"));
    assert_eq!((e("counted"), e("mapq"), e("dupe"), e("baseq"), e("overlap"), e("capped")), (103, 20, 20, 5, 12, 0));
    close(&r, "/coverage/exclusions/mapq", 20.0 / 160.0);
    close(&r, "/coverage/exclusions/overlap", 12.0 / 160.0);
    close(&r, "/coverage/exclusions/total", 57.0 / 160.0);
}

#[test]
fn no_raw_drops_the_raw_array() {
    let mut c = config("coverage.sam");
    c.params.raw = false;
    let (r, d) = run_cfg(&c);
    assert!(get(&r, "/coverage/raw").is_null());
    assert!(d.iter().all(|(_, _, raw)| raw.is_none()));
}

/// 70,000 identical unpaired reads: raw depth saturates at 65,535; filtered
/// keeps one (the rest are fragment duplicates).
#[test]
fn raw_depth_saturates() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("deep.sam");
    let mut s = std::fs::read_to_string(fixture("header.txt")).unwrap();
    for i in 0..70_000 {
        s.push_str(&format!("r{i}:1:FCA:1:1101:{i}:1\t0\tchrM\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"));
    }
    std::fs::write(&path, s).unwrap();
    let mut c = config("coverage.sam");
    c.input = Some(path);
    let (r, d) = run_cfg(&c);
    assert_eq!(depth_at(&d, "chrM", 1), (1, Some(65_535)));
    assert_eq!(u(&r, "/duplicates/duplicate_fragments"), 69_999);
    assert_eq!(get(&r, "/coverage/raw/histogram").as_array().unwrap().len(), 65_536);
}
