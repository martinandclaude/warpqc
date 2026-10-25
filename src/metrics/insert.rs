//! Insert-size histograms per pair orientation, following Picard
//! CollectInsertSizeMetrics: one count per pair from R2's |TLEN|, pairs with
//! both mates mapped on the same contig, non-duplicate, TLEN != 0.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::record::Aln;
use crate::summary::hist;

/// Inserts are binned densely up to this size and sparsely above it.
/// The JSON histogram stops at `REPORT_MAX`.
pub const MAX_INSERT: usize = 100_000;
pub const REPORT_MAX: usize = 2_000;
/// Picard's DEVIATIONS: mean and SD use inserts <= median + 10 * MAD.
pub const DEVIATIONS: f64 = 10.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    Fr = 0,
    Rf = 1,
    Tandem = 2,
}

/// htsjdk `SamPairUtil.getPairOrientation`, using both mates directly.
pub fn orientation(a: &Aln, b: &Aln) -> Orientation {
    if a.is_reverse() == b.is_reverse() {
        return Orientation::Tandem;
    }
    let (fwd, rev) = if a.is_reverse() { (b, a) } else { (a, b) };
    if fwd.pos < rev.end() {
        Orientation::Fr
    } else {
        Orientation::Rf
    }
}

#[derive(Debug, Clone, Default)]
pub struct InsertStats {
    /// Dense counts for inserts up to `MAX_INSERT` ...
    hist: [Vec<u64>; 3],
    /// ... and exact sparse counts above it, so order statistics stay exact.
    big: [BTreeMap<u64, u64>; 3],
}

impl InsertStats {
    pub fn add(&mut self, r1: &Aln, r2: &Aln) {
        if r2.tlen == 0 {
            return;
        }
        let o = orientation(r1, r2) as usize;
        let v = r2.tlen.unsigned_abs();
        if v as usize > MAX_INSERT {
            *self.big[o].entry(v).or_default() += 1;
            return;
        }
        let h = &mut self.hist[o];
        if h.is_empty() {
            h.resize(MAX_INSERT + 1, 0);
        }
        h[v as usize] += 1;
    }

    pub fn merge(&mut self, o: &InsertStats) {
        for k in 0..3 {
            if !o.hist[k].is_empty() {
                if self.hist[k].is_empty() {
                    self.hist[k].resize(MAX_INSERT + 1, 0);
                }
                for (a, b) in self.hist[k].iter_mut().zip(&o.hist[k]) {
                    *a += b;
                }
            }
            for (&v, &c) in &o.big[k] {
                *self.big[k].entry(v).or_default() += c;
            }
        }
    }

    fn bins(&self, k: usize) -> Vec<(u64, u64)> {
        let mut b = hist::sparse(&self.hist[k]);
        b.extend(self.big[k].iter().map(|(&v, &c)| (v, c)));
        b
    }

    pub fn report(&self) -> InsertReport {
        let bins: Vec<_> = (0..3).map(|k| self.bins(k)).collect();
        let total: u64 = bins.iter().map(|b| hist::count(b)).sum();
        let orient = |k: usize| OrientationReport::new(&bins[k], total);
        InsertReport {
            fr: orient(0),
            rf: orient(1),
            tandem: orient(2),
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct InsertReport {
    pub fr: OrientationReport,
    pub rf: OrientationReport,
    pub tandem: OrientationReport,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct OrientationReport {
    pub pairs: u64,
    pub fraction: f64,
    pub median: f64,
    pub mad: f64,
    pub mode: u64,
    pub min: u64,
    pub max: u64,
    pub p10: u64,
    pub p25: u64,
    pub p75: u64,
    pub p90: u64,
    /// Over inserts <= median + 10 MAD (Picard's upper-only trim).
    pub mean: f64,
    pub sd: f64,
    pub trim_width: u64,
    /// Counts for 0..=2000 bp (trailing zeros dropped) ...
    pub histogram: Vec<u64>,
    /// ... and above 2000 bp.
    pub histogram_overflow: u64,
}

impl OrientationReport {
    /// `h`: sorted (insert, count) bins.
    fn new(h: &hist::Bins, total: u64) -> OrientationReport {
        let pairs = hist::count(h);
        let fraction = if total > 0 {
            pairs as f64 / total as f64
        } else {
            0.0
        };
        if pairs == 0 {
            return OrientationReport {
                fraction,
                ..Default::default()
            };
        }
        let median = hist::median(h);
        let mad = hist::mad(h);
        let width = (median + DEVIATIONS * mad) as u64;
        let trimmed = hist::trim_by_width(h, width);
        let shown = hist::trim_by_width(h, REPORT_MAX as u64);
        let mut histogram = vec![0u64; shown.last().map_or(0, |b| b.0 as usize + 1)];
        for &(v, c) in shown {
            histogram[v as usize] = c;
        }
        OrientationReport {
            pairs,
            fraction,
            median,
            mad,
            mode: hist::mode(h),
            min: h[0].0,
            max: h[h.len() - 1].0,
            p10: hist::percentile(h, 0.10),
            p25: hist::percentile(h, 0.25),
            p75: hist::percentile(h, 0.75),
            p90: hist::percentile(h, 0.90),
            mean: hist::mean(trimmed),
            sd: hist::sd(trimmed),
            trim_width: width,
            histogram,
            histogram_overflow: pairs - hist::count(shown),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn upper_only_trim() {
        // 100 inserts at 300, 1 at 299 and 1 at 301: median 300, MAD 0 ->
        // width 300, so 301 and the 5000 outlier leave the mean but not the
        // order statistics; 299 stays because the trim is upper-only.
        let h = [(299, 1), (300, 100), (301, 1), (5000, 1)];
        let r = OrientationReport::new(&h, 103);
        assert_eq!((r.median, r.mad, r.trim_width), (300.0, 0.0, 300));
        assert!((r.mean - (300.0 * 100.0 + 299.0) / 101.0).abs() < 1e-9);
        assert_eq!(r.max, 5000);
        assert_eq!(r.histogram.len(), 302);
        assert_eq!(r.histogram_overflow, 1);
    }

    fn aln(pos1: i64, len: u32, reverse: bool) -> Aln {
        Aln {
            pos: pos1 - 1,
            cigar: vec![(crate::record::op::M, len)],
            flag: if reverse {
                crate::record::flag::REVERSE
            } else {
                0
            },
            ..Default::default()
        }
    }

    /// htsjdk >= 5.0.0 (Picard >= 3.5.0), symmetric in the two mates; the
    /// reviewer's two cases, which Picard 3.4.0 calls RF.
    #[test]
    fn orientation_matches_htsjdk_5() {
        // Forward 110-119, reverse 101-110: forward start 110 <= reverse end 110.
        let (f, r) = (aln(110, 10, false), aln(101, 10, true));
        assert_eq!(orientation(&f, &r), Orientation::Fr);
        assert_eq!(orientation(&r, &f), Orientation::Fr);
        // Dovetail: reverse 101-110, forward 105-114.
        let (r, f) = (aln(101, 10, true), aln(105, 10, false));
        assert_eq!(orientation(&r, &f), Orientation::Fr);
        // Forward 111-120 starts past the reverse end 110: RF.
        assert_eq!(
            orientation(&aln(111, 10, false), &aln(101, 10, true)),
            Orientation::Rf
        );
        assert_eq!(
            orientation(&aln(1, 10, true), &aln(50, 10, true)),
            Orientation::Tandem
        );
    }

    #[test]
    fn inserts_above_max_stay_in_order_statistics() {
        // Review finding: 300, 310 and 150,000 -> median 310, MAD 10 (Picard).
        let mut s = InsertStats::default();
        let r1 = Aln::default();
        let mut r2 = Aln {
            flag: crate::record::flag::REVERSE,
            cigar: vec![(crate::record::op::M, 10)],
            pos: 1_000,
            ..Default::default()
        };
        for t in [300i64, 310, 150_000] {
            r2.tlen = -t;
            s.add(&r1, &r2);
        }
        let fr = s.report().fr;
        assert_eq!(
            (fr.pairs, fr.median, fr.mad, fr.max),
            (3, 310.0, 10.0, 150_000)
        );
        assert_eq!(fr.histogram_overflow, 1);
    }
}
