//! Uniformity from binned read starts: PMAD (raw and GC-corrected), CV of
//! bin counts, the GC-bias curve with AT/GC dropout, a 1 Mb coverage profile;
//! and the Oexle/Mokry evenness score from the per-base depth histogram.
//!
//! PMAD follows the DRAGEN definition (locked decision 8):
//! v[i] = log2(count[i] / median count) over retained bins, zero-count bins
//! dropped before the log; d[i] = v[i] - v[i+1] within a contig only;
//! PMAD = median(|d - median(d)|), unscaled. Medians are htsjdk-style (mean
//! of the two middle values for an even count).

use serde::{Deserialize, Serialize};

use crate::metrics::bins::BinLayout;
use crate::reference::NO_GC;

/// GC strata with fewer bins than this are merged with a neighbour.
pub const MIN_STRATUM_BINS: usize = 100;
pub const PROFILE_WINDOW_BP: u64 = 1_000_000;

fn median(v: &mut [f64]) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    v.sort_by(f64::total_cmp);
    let n = v.len();
    Some(if n.is_multiple_of(2) {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    } else {
        v[n / 2]
    })
}

#[derive(Debug, Clone, PartialEq)]
pub struct Pmad {
    pub pmad: Option<f64>,
    pub median: f64,
    pub zero_dropped: u64,
}

/// PMAD from per-contig vectors of (normalised) counts of retained bins.
pub fn pmad(contigs: &[Vec<f64>]) -> Pmad {
    let mut all: Vec<f64> = contigs.iter().flatten().copied().collect();
    let m = median(&mut all).unwrap_or(0.0);
    let mut zero_dropped = 0;
    let mut d = Vec::new();
    for c in contigs {
        let mut prev: Option<f64> = None;
        for &x in c {
            if x <= 0.0 {
                zero_dropped += 1;
                continue;
            }
            let v = (x / m).log2();
            if let Some(p) = prev {
                d.push(p - v);
            }
            prev = Some(v);
        }
    }
    let pmad = if m > 0.0 {
        median(&mut d.clone()).and_then(|md| {
            let mut dev: Vec<f64> = d.iter().map(|x| (x - md).abs()).collect();
            median(&mut dev)
        })
    } else {
        None
    };
    Pmad {
        pmad,
        median: m,
        zero_dropped,
    }
}

/// Groups integer GC% values: walking up from 0, strata are added to the
/// current group until it holds `min` bins; a short tail joins the group
/// before it. Returns the group id per GC% (0..=100) and the group count.
pub fn gc_groups(bins_per_gc: &[usize; 101], min: usize) -> ([usize; 101], usize) {
    let mut group = [0usize; 101];
    let (mut g, mut in_group) = (0usize, 0usize);
    let mut last_closed: Option<usize> = None;
    for pct in 0..=100 {
        group[pct] = g;
        in_group += bins_per_gc[pct];
        if in_group >= min {
            last_closed = Some(g);
            g += 1;
            in_group = 0;
        }
    }
    if in_group > 0 || last_closed.is_none() {
        // Unclosed tail: merge into the last closed group, if any.
        if let Some(lc) = last_closed {
            for x in group.iter_mut() {
                if *x == g {
                    *x = lc;
                }
            }
            (group, lc + 1)
        } else {
            (group, 1)
        }
    } else {
        (group, g)
    }
}

/// Evenness score E (Mokry et al. 2010; Oexle 2016) from a depth histogram:
/// C = round(mean), D2 = positions with depth <= C,
/// E = 1 - (|D2| - sum(D2) / C) / N. As the Sequana implementation, the mean
/// is rounded half to even. None when C is 0.
pub fn evenness(hist: &[u64]) -> Option<f64> {
    let n: u64 = hist.iter().sum();
    if n == 0 {
        return None;
    }
    let mean = hist
        .iter()
        .enumerate()
        .map(|(d, &c)| d as f64 * c as f64)
        .sum::<f64>()
        / n as f64;
    let c = mean.round_ties_even();
    if c < 1.0 {
        return None;
    }
    let (mut count2, mut sum2) = (0u64, 0f64);
    for (d, &k) in hist.iter().enumerate().take(c as usize + 1) {
        count2 += k;
        sum2 += d as f64 * k as f64;
    }
    Some(1.0 - (count2 as f64 - sum2 / c) / n as f64)
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct GcPoint {
    pub gc: u8,
    pub bins: u64,
    pub reads: u64,
    /// (reads per bin at this GC) / (reads per bin overall).
    pub normalized_coverage: f64,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ProfileContig {
    pub name: String,
    /// Mean of count / median count over the window's retained bins.
    pub values: Vec<Option<f64>>,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Profile {
    pub window_bp: u64,
    pub contigs: Vec<ProfileContig>,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct UniformityReport {
    pub bin_size: u32,
    pub bins_total: u64,
    pub bins_excluded_n: u64,
    pub bins_excluded_bed: u64,
    pub bins_used: u64,
    pub median_bin_count: f64,
    pub zero_bins_dropped: u64,
    pub gc_strata: u64,
    /// Bins dropped from the GC-corrected PMAD because their stratum median is 0.
    pub gc_dropped_bins: u64,
    pub pmad: Option<f64>,
    pub pmad_gc: Option<f64>,
    pub cv: Option<f64>,
    pub evenness: Option<f64>,
    pub at_dropout: f64,
    pub gc_dropout: f64,
    pub gc_curve: Vec<GcPoint>,
    pub profile: Profile,
}

pub fn report(
    bins: &BinLayout,
    counts: &[u32],
    depth_hist: &[u64],
    min_stratum: usize,
) -> UniformityReport {
    let zero = vec![0u32; bins.len()];
    let counts = if counts.is_empty() { &zero[..] } else { counts };
    let mut excluded_n = 0;
    let mut excluded_bed = 0;
    let mut per_contig: Vec<Vec<(f64, u8)>> = Vec::new();
    for c in &bins.contigs {
        let mut v = Vec::with_capacity(c.n);
        let r = c.first..c.first + c.n;
        for ((&count, &gc), &bed) in counts[r.clone()]
            .iter()
            .zip(&bins.gc[r.clone()])
            .zip(&bins.bed_excluded[r])
        {
            if gc == NO_GC {
                excluded_n += 1;
            } else if bed {
                excluded_bed += 1;
            } else {
                v.push((count as f64, gc));
            }
        }
        per_contig.push(v);
    }

    let raw: Vec<Vec<f64>> = per_contig
        .iter()
        .map(|v| v.iter().map(|x| x.0).collect())
        .collect();
    let raw_p = pmad(&raw);

    // GC strata and their medians.
    let mut bins_per_gc = [0usize; 101];
    let mut reads_per_gc = [0u64; 101];
    for &(c, g) in per_contig.iter().flatten() {
        bins_per_gc[g as usize] += 1;
        reads_per_gc[g as usize] += c as u64;
    }
    let (group, n_groups) = gc_groups(&bins_per_gc, min_stratum);
    let mut members: Vec<Vec<f64>> = vec![Vec::new(); n_groups];
    for &(c, g) in per_contig.iter().flatten() {
        members[group[g as usize]].push(c);
    }
    let group_median: Vec<f64> = members
        .iter_mut()
        .map(|m| median(m).unwrap_or(0.0))
        .collect();
    let mut gc_dropped = 0;
    let corrected: Vec<Vec<f64>> = per_contig
        .iter()
        .map(|v| {
            v.iter()
                .filter_map(|&(c, g)| {
                    let m = group_median[group[g as usize]];
                    if m > 0.0 {
                        Some(c / m)
                    } else {
                        gc_dropped += 1;
                        None
                    }
                })
                .collect()
        })
        .collect();
    let gc_p = pmad(&corrected);

    let used: Vec<f64> = raw.iter().flatten().copied().collect();
    let n = used.len() as f64;
    let cv = (used.len() >= 2).then(|| {
        let mean = used.iter().sum::<f64>() / n;
        let sd = (used.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
        sd / mean
    });
    let cv = cv.filter(|x| x.is_finite());

    let total_bins: u64 = bins_per_gc.iter().map(|&b| b as u64).sum();
    let total_reads: u64 = reads_per_gc.iter().sum();
    let mut gc_curve = Vec::new();
    let (mut at, mut gcd) = (0.0, 0.0);
    for g in 0..=100 {
        if bins_per_gc[g] == 0 {
            continue;
        }
        let fb = bins_per_gc[g] as f64 / total_bins as f64;
        let fr = if total_reads > 0 {
            reads_per_gc[g] as f64 / total_reads as f64
        } else {
            0.0
        };
        if g <= 50 {
            at += (fb - fr).max(0.0);
        }
        if g >= 50 {
            gcd += (fb - fr).max(0.0);
        }
        gc_curve.push(GcPoint {
            gc: g as u8,
            bins: bins_per_gc[g] as u64,
            reads: reads_per_gc[g],
            normalized_coverage: if total_reads > 0 { fr / fb } else { 0.0 },
        });
    }

    let window = ((PROFILE_WINDOW_BP / bins.size as u64).max(1)) as usize;
    let m = raw_p.median;
    let profile = Profile {
        window_bp: window as u64 * bins.size as u64,
        contigs: bins
            .contigs
            .iter()
            .map(|c| ProfileContig {
                name: c.name.clone(),
                values: (c.first..c.first + c.n)
                    .collect::<Vec<_>>()
                    .chunks(window)
                    .map(|w| {
                        let vals: Vec<f64> = w
                            .iter()
                            .filter(|&&i| bins.usable(i))
                            .map(|&i| counts[i] as f64)
                            .collect();
                        (m > 0.0 && !vals.is_empty())
                            .then(|| vals.iter().sum::<f64>() / vals.len() as f64 / m)
                    })
                    .collect(),
            })
            .collect(),
    };

    UniformityReport {
        bin_size: bins.size,
        bins_total: bins.len() as u64,
        bins_excluded_n: excluded_n,
        bins_excluded_bed: excluded_bed,
        bins_used: total_bins,
        median_bin_count: raw_p.median,
        zero_bins_dropped: raw_p.zero_dropped,
        gc_strata: n_groups as u64,
        gc_dropped_bins: gc_dropped,
        pmad: raw_p.pmad,
        pmad_gc: gc_p.pmad,
        cv,
        evenness: evenness(depth_hist),
        at_dropout: at,
        gc_dropout: gcd,
        gc_curve,
        profile,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hand_computed_vector_with_zero_bin() {
        // median of [0,100,100,100,200,400] = 100; v = [0,1,0,(zero dropped),2,0]
        // d = [-1, 1, -2, 2]; median(d) = 0; |d| = [1,1,2,2] -> 1.5
        let p = pmad(&[vec![100.0, 200.0, 100.0, 0.0, 400.0, 100.0]]);
        assert_eq!(
            p,
            Pmad {
                pmad: Some(1.5),
                median: 100.0,
                zero_dropped: 1
            }
        );
    }

    #[test]
    fn no_difference_across_contig_boundary() {
        // Within contigs: d = [-1] and [2] -> median 0.5 -> |d - 0.5| = 1.5.
        // Joined, the boundary difference log2(200/400) = -1 would give 0.
        assert_eq!(
            pmad(&[vec![100.0, 200.0], vec![400.0, 100.0]]).pmad,
            Some(1.5)
        );
        assert_eq!(pmad(&[vec![100.0, 200.0, 400.0, 100.0]]).pmad, Some(0.0));
    }

    #[test]
    fn matches_independent_python_computation() {
        // Expected value from an independent Python computation:
        //   m = median(a + b); v = log2(x / m) for x > 0, per contig;
        //   d = [v[i] - v[i+1]]; pmad = median(|d - median(d)|)
        // with statistics.median throughout.
        let a = vec![
            212.0, 187.0, 240.0, 0.0, 199.0, 260.0, 175.0, 230.0, 205.0, 190.0,
        ];
        let b = vec![150.0, 310.0, 220.0, 205.0, 0.0, 0.0, 180.0, 195.0];
        let p = pmad(&[a, b]);
        assert!((p.pmad.unwrap() - PY_PMAD).abs() < 1e-12, "{p:?}");
        assert_eq!(p.zero_dropped, 3);
    }
    const PY_PMAD: f64 = 0.22510170859443404;

    #[test]
    fn gc_correction_removes_a_pure_gc_effect() {
        // GC 40 bins have 100 reads, GC 60 bins 200, in an irregular order
        // (40,40,60,40,60,60,40,60): raw d = [0,-1,1,-1,0,1,-1] -> PMAD 1;
        // after dividing by stratum medians every bin is 1 -> PMAD 0.
        let mut bins_per_gc = [0usize; 101];
        bins_per_gc[40] = 4;
        bins_per_gc[60] = 4;
        let (group, n) = gc_groups(&bins_per_gc, 2);
        assert_eq!(n, 2);
        assert_ne!(group[40], group[60]);
        let raw = vec![100.0, 100.0, 200.0, 100.0, 200.0, 200.0, 100.0, 200.0];
        assert_eq!(pmad(&[raw]).pmad, Some(1.0));
        let corrected = vec![1.0; 8];
        assert_eq!(pmad(&[corrected]).pmad, Some(0.0));
    }

    #[test]
    fn small_strata_merge_upwards_and_tail_merges_down() {
        let mut b = [0usize; 101];
        b[10] = 50;
        b[11] = 60;
        b[12] = 120;
        b[13] = 30;
        let (g, n) = gc_groups(&b, 100);
        assert_eq!(n, 2);
        assert_eq!((g[10], g[11], g[12], g[13], g[100]), (0, 0, 1, 1, 1));
        // Nothing reaches the minimum: one group.
        let mut few = [0usize; 101];
        few[40] = 5;
        assert_eq!(gc_groups(&few, 100).1, 1);
    }

    #[test]
    fn evenness_hand_computed() {
        // Depths 0,1,2,3,4 once each: mean 2, C 2, D2 = {0,1,2}:
        // E = 1 - (3 - 3/2) / 5 = 0.7.
        assert!((evenness(&[1, 1, 1, 1, 1]).unwrap() - 0.7).abs() < 1e-12);
        // Perfectly even coverage gives 1.
        assert_eq!(evenness(&[0, 0, 0, 0, 0, 10]), Some(1.0));
        assert_eq!(evenness(&[5]), None);
    }
}
