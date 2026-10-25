mod common;
use common::*;

#[test]
fn reverse_trailing_softclip_groups_by_unclipped_5p() {
    // Members differ in clips and TLEN, so everything fed by the kept member
    // (coverage, insert size, bins, per-lane insert median) depends on thread
    // timing (locked decision 5). The duplicate counts below do not.
    let (r, _) = run_fixture_except(
        "dup_softclip_reverse.sam",
        &[
            "/coverage",
            "/insert_size",
            "/uniformity",
            "/per_lane/FCA:1/insert_median",
        ],
    );
    assert_eq!(u(&r, "/duplicates/pairs_examined"), 5);
    assert_eq!(u(&r, "/duplicates/duplicate_pairs"), 3);
    assert_eq!(u(&r, "/duplicates/duplicate_groups"), 1);
    assert_eq!(u(&r, "/duplicates/distinct_pair_signatures"), 2);
    close(&r, "/duplicates/rate", 0.6);
    assert_eq!(u(&r, "/duplicates/optical_duplicate_pairs"), 0);
}

#[test]
fn fragments_dedup_against_fragments_only() {
    let (r, _) = run_fixture("fragments.sam");
    assert_eq!(u(&r, "/duplicates/fragments_examined"), 4);
    assert_eq!(u(&r, "/duplicates/duplicate_fragments"), 2);
    assert_eq!(u(&r, "/duplicates/pairs_examined"), 1);
    assert_eq!(u(&r, "/duplicates/duplicate_pairs"), 0);
    close(&r, "/duplicates/rate", 2.0 / 6.0);
    assert!(
        get(&r, "/definitions/duplicates.fragments")
            .as_str()
            .unwrap()
            .contains("fragments only")
    );
}

#[test]
fn supplementary_and_secondary_are_ignored() {
    let (r, d) = run_fixture("suppl_secondary.sam");
    assert_eq!(u(&r, "/duplicates/pairs_examined"), 2);
    assert_eq!(u(&r, "/duplicates/duplicate_pairs"), 0);
    for p in 11..=20 {
        assert_eq!(depth_at(&d, "chr1", p), (1, Some(1)), "chr1:{p}");
    }
    for p in 101..=110 {
        assert_eq!(depth_at(&d, "chr1", p), (1, Some(1)), "chr1:{p}");
    }
}

#[test]
fn optical_chain_and_cross_lane() {
    let (r, _) = run_fixture("optical_chain.sam");
    assert_eq!(u(&r, "/duplicates/pairs_examined"), 7);
    assert_eq!(u(&r, "/duplicates/duplicate_pairs"), 5);
    assert_eq!(u(&r, "/duplicates/optical_duplicate_pairs"), 2);
    close(&r, "/duplicates/rate", 10.0 / 14.0);
    close(&r, "/duplicates/optical_rate", 4.0 / 14.0);
    close(&r, "/duplicates/non_optical_rate", 6.0 / 14.0);
    close(&r, "/duplicates/cross_lane_fraction", 0.8);
    // Picard estimate with 7 - 2 read pairs and 7 - 5 unique pairs.
    assert_eq!(u(&r, "/duplicates/est_library_size"), 2);
    assert_eq!(u(&r, "/per_lane/FCA:1/pairs_examined"), 6);
    // Fractional attribution: group 1 (n = 5) has 4 members in FCA:1 and 1 in
    // FCA:2, each counting 4/5; group 2 (n = 2) is all FCA:1, 1/2 each.
    close(&r, "/per_lane/FCA:1/duplicate_pairs", 4.0 * 0.8 + 2.0 * 0.5);
    close(&r, "/per_lane/FCA:1/duplicate_rate", 4.2 / 6.0);
    assert_eq!(u(&r, "/per_lane/FCA:1/optical_duplicate_pairs"), 2);
    assert_eq!(u(&r, "/per_lane/FCA:2/pairs_examined"), 1);
    close(&r, "/per_lane/FCA:2/duplicate_pairs", 0.8);
    assert_eq!(
        get(&r, "/inputs/lanes"),
        &serde_json::json!(["FCA:1", "FCA:2"])
    );
}

#[test]
fn optical_distance_is_a_parameter() {
    // At 1999 px the A-B and B-C links (2000 px per axis) break.
    let mut c = config("optical_chain.sam");
    c.params.optical_distance = 1999;
    let (r, _) = run_cfg(&c);
    assert_eq!(u(&r, "/duplicates/optical_duplicate_pairs"), 0);
}
