mod common;
use common::*;
use warpqc::pipeline::RunConfig;

fn bins_config() -> RunConfig {
    let mut c = config("bins.sam");
    c.bin_size = 20;
    c.exclude_bed = Some(fixture("bins.bed"));
    c
}

#[test]
fn bin_counts_and_exclusions() {
    let one = bins_config();
    let (r, _) = run_cfg(&one);
    let mut many = bins_config();
    many.threads = 4;
    many.block_size = 64;
    let (r4, _) = run_cfg(&many);
    assert_eq!(r, r4, "1 and 4 threads disagree");

    let out = warpqc::pipeline::run(&one).unwrap();
    let chr1 = &out.bins[0];
    let chr2 = &out.bins[1];
    assert_eq!(out.bins.len(), 2, "only autosomes are binned");
    let ones = |v: &[u32]| {
        v.iter()
            .enumerate()
            .filter(|(_, c)| **c == 1)
            .map(|(i, _)| i)
            .collect::<Vec<_>>()
    };
    assert_eq!(chr1.0, "chr1");
    assert_eq!(ones(&chr1.1), vec![0, 1, 2, 3]);
    assert!(chr1.1.iter().all(|&c| c <= 1));
    assert_eq!(ones(&chr2.1), vec![1, 6, 8, 14]);
    assert!(chr2.1.iter().all(|&c| c <= 1));
    let unusable = |v: &[bool]| {
        v.iter()
            .enumerate()
            .filter(|(_, u)| !**u)
            .map(|(i, _)| i)
            .collect::<Vec<_>>()
    };
    assert_eq!(unusable(&chr1.2), vec![3, 10]);
    assert_eq!(unusable(&chr2.2), vec![0]);

    let un = |p: &str| format!("/uniformity/{p}");
    assert_eq!(u(&r, &un("bin_size")), 20);
    assert_eq!(u(&r, &un("bins_total")), 35);
    assert_eq!(u(&r, &un("bins_excluded_n")), 1);
    assert_eq!(u(&r, &un("bins_excluded_bed")), 2);
    assert_eq!(u(&r, &un("bins_used")), 32);
    // 7 used bins hold a read; with 25 empty bins the median is 0, so PMAD
    // is undefined.
    assert_eq!(u(&r, &un("zero_bins_dropped")), 25);
    assert!(get(&r, &un("pmad")).is_null());
    let w = get(&r, "/warnings").as_array().unwrap();
    assert!(
        w.iter().any(|w| w.as_str().unwrap().contains("chrZ")),
        "{w:?}"
    );
}

/// A dense synthetic sample: 2 autosomes x 100 bins, reads in every bin, so
/// PMAD, CV and the GC curve are all defined and reproducible.
#[test]
fn dense_sample_gives_defined_uniformity() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("dense.sam");
    let mut s = std::fs::read_to_string(fixture("header.txt")).unwrap();
    let mut n = 0;
    for (contig, len) in [("chr1", 400), ("chr2", 300)] {
        for start in (1..=len - 10).step_by(2) {
            // Reads at every second base: 10 per 20 bp bin.
            s.push_str(&format!(
                "r{n}:1:FCA:1:1101:{n}:1\t0\t{contig}\t{start}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
            ));
            n += 1;
        }
    }
    std::fs::write(&path, s).unwrap();
    let mut c = config("bins.sam");
    c.input = Some(path);
    c.bin_size = 20;
    let (r, _) = run_cfg(&c);
    let pm = get(&r, "/uniformity/pmad").as_f64().unwrap();
    assert_eq!(
        pm, 0.0,
        "every full bin has 10 reads except where reads stop"
    );
    assert!(get(&r, "/uniformity/cv").as_f64().is_some());
    assert!(
        !get(&r, "/uniformity/gc_curve")
            .as_array()
            .unwrap()
            .is_empty()
    );
    assert_eq!(
        get(&r, "/uniformity/profile/contigs")
            .as_array()
            .unwrap()
            .len(),
        2
    );
}
