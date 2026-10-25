mod common;
use common::*;

#[test]
fn flags_mixed_matches_hand_counts_and_flagstat() {
    let (r, _) = run_fixture("flags_mixed.sam");
    let a = |p: &str| format!("/alignment/{p}");
    assert_eq!(u(&r, &a("total_records")), 17);
    assert_eq!(u(&r, &a("total_templates")), 8);
    assert_eq!(u(&r, &a("total_reads")), 15);
    assert_eq!(u(&r, &a("secondary_records")), 1);
    assert_eq!(u(&r, &a("supplementary_records")), 1);
    assert_eq!(u(&r, &a("qc_fail_reads")), 2);
    assert_eq!(u(&r, &a("primary_mapped")), 12);
    assert_eq!(u(&r, &a("paired_reads")), 14);
    close(&r, &a("primary_mapped_fraction"), 12.0 / 15.0);
    close(&r, &a("properly_paired_fraction"), 6.0 / 14.0);
    close(&r, &a("both_mates_mapped_fraction"), 10.0 / 14.0);
    close(&r, &a("singleton_fraction"), 1.0 / 14.0);
    assert_eq!(u(&r, &a("mate_diff_contig_mapq5_reads")), 3);
    close(&r, &a("chimeric_read_fraction"), 3.0 / 10.0);
    assert_eq!(u(&r, &a("chimeric_pairs_mapq5")), 1);
    close(&r, &a("chimeric_pair_fraction"), 1.0 / 5.0);
    let mq = get(&r, &a("mapq_histogram")).as_array().unwrap();
    assert_eq!(mq.len(), 61);
    assert_eq!((mq[60].as_u64(), mq[30].as_u64(), mq[20].as_u64(), mq[3].as_u64(), mq[0].as_u64()),
               (Some(8), Some(1), Some(1), Some(1), Some(1)));
    close(&r, &a("mapq0_fraction"), 1.0 / 12.0);
    close(&r, &a("mapq_ge20_fraction"), 10.0 / 12.0);
    close(&r, &a("mapq_ge60_fraction"), 8.0 / 12.0);
    close(&r, &a("off_target_fraction"), 1.0 / 12.0);
    close(&r, &a("soft_clipped_fraction"), 5.0 / 125.0);
    assert_eq!(u(&r, &a("aligned_bases")), 119);
    close(&r, &a("mismatch_rate"), 2.0 / 109.0);
    close(&r, &a("indel_rate_per_100"), 100.0 / 119.0);
    assert_eq!(u(&r, &a("reads_without_nm")), 1);
    assert_eq!(u(&r, &a("incomplete_templates")), 0);
}

/// The same counts as `samtools flagstat` (QC-pass + QC-fail), if samtools
/// is installed.
#[test]
fn flags_mixed_agrees_with_samtools_flagstat() {
    let Ok(out) = std::process::Command::new("samtools")
        .args(["flagstat", "-O", "tsv"])
        .arg(fixture("flags_mixed.sam"))
        .output()
    else {
        eprintln!("samtools not found; skipping");
        return;
    };
    let text = String::from_utf8(out.stdout).unwrap();
    let field = |label: &str| -> u64 {
        let l = text.lines().find(|l| l.split('\t').nth(2) == Some(label)).unwrap();
        let f: Vec<&str> = l.split('\t').collect();
        f[0].parse::<u64>().unwrap() + f[1].parse::<u64>().unwrap()
    };
    let (r, _) = run_fixture("flags_mixed.sam");
    assert_eq!(u(&r, "/alignment/total_records"), field("total (QC-passed reads + QC-failed reads)"));
    assert_eq!(u(&r, "/alignment/total_reads"), field("primary"));
    assert_eq!(u(&r, "/alignment/primary_mapped"), field("primary mapped"));
    assert_eq!(u(&r, "/alignment/paired_reads"), field("paired in sequencing"));
    assert_eq!(u(&r, "/alignment/mate_diff_contig_mapq5_reads"), field("with mate mapped to a different chr (mapQ>=5)"));
}
