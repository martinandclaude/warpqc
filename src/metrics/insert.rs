//! Insert-size histograms per pair orientation, following Picard
//! CollectInsertSizeMetrics: one count per pair from R2's |TLEN|, pairs with
//! both mates mapped on the same contig, non-duplicate, TLEN != 0.

use serde::Serialize;

use crate::record::Aln;
use crate::summary::hist;

/// Inserts are binned exactly up to this size; larger ones go to an
/// overflow count. The JSON histogram stops at `REPORT_MAX`.
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
    if fwd.pos <= rev.end() - 1 { Orientation::Fr } else { Orientation::Rf }
}

#[derive(Debug, Clone, Default)]
pub struct InsertStats {
    hist: [Vec<u64>; 3],
    overflow: [u64; 3],
}

impl InsertStats {
    pub fn add(&mut self, r1: &Aln, r2: &Aln) {
        if r2.tlen == 0 {
            return;
        }
        let o = orientation(r1, r2) as usize;
        let v = r2.tlen.unsigned_abs() as usize;
        if v > MAX_INSERT {
            self.overflow[o] += 1;
            return;
        }
        let h = &mut self.hist[o];
        if h.is_empty() {
            h.resize(MAX_INSERT + 1, 0);
        }
        h[v] += 1;
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
            self.overflow[k] += o.overflow[k];
        }
    }

    pub fn report(&self) -> InsertReport {
        let total: u64 = (0..3).map(|k| hist::count(&self.hist[k]) + self.overflow[k]).sum();
        let orient = |k: usize| OrientationReport::new(&self.hist[k], self.overflow[k], total);
        InsertReport {
            fr: orient(0),
            rf: orient(1),
            tandem: orient(2),
        }
    }
}

#[derive(Debug, Serialize)]
pub struct InsertReport {
    pub fr: OrientationReport,
    pub rf: OrientationReport,
    pub tandem: OrientationReport,
}

#[derive(Debug, Serialize, Default)]
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
    /// Inserts above `MAX_INSERT`; they enter the counts but not the
    /// order statistics.
    pub above_max: u64,
    /// Counts for 0..=2000 bp (trailing zeros dropped) ...
    pub histogram: Vec<u64>,
    /// ... and above 2000 bp.
    pub histogram_overflow: u64,
}

impl OrientationReport {
    fn new(h: &[u64], overflow: u64, total: u64) -> OrientationReport {
        let n = hist::count(h);
        let pairs = n + overflow;
        if n == 0 {
            return OrientationReport {
                pairs,
                fraction: if total > 0 { pairs as f64 / total as f64 } else { 0.0 },
                above_max: overflow,
                histogram_overflow: overflow,
                ..Default::default()
            };
        }
        let median = hist::median(h);
        let mad = hist::mad(h);
        let width = (median + DEVIATIONS * mad) as usize;
        let trimmed = hist::trim_by_width(h, width);
        let shown = &h[..hist::used_len(&h[..=REPORT_MAX.min(h.len() - 1)])];
        OrientationReport {
            pairs,
            fraction: pairs as f64 / total as f64,
            median,
            mad,
            mode: hist::mode(h),
            min: h.iter().position(|&c| c > 0).unwrap() as u64,
            max: (hist::used_len(h) - 1) as u64,
            p10: hist::percentile(h, 0.10),
            p25: hist::percentile(h, 0.25),
            p75: hist::percentile(h, 0.75),
            p90: hist::percentile(h, 0.90),
            mean: hist::mean(trimmed),
            sd: hist::sd(trimmed),
            trim_width: width as u64,
            above_max: overflow,
            histogram: shown.to_vec(),
            histogram_overflow: h.iter().skip(REPORT_MAX + 1).sum::<u64>() + overflow,
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
        let mut h = vec![0u64; MAX_INSERT + 1];
        h[300] = 100;
        h[299] = 1;
        h[301] = 1;
        h[5000] = 1;
        let r = OrientationReport::new(&h, 0, 103);
        assert_eq!((r.median, r.mad, r.trim_width), (300.0, 0.0, 300));
        assert!((r.mean - (300.0 * 100.0 + 299.0) / 101.0).abs() < 1e-9);
        assert_eq!(r.max, 5000);
        assert_eq!(r.histogram.len(), 302);
        assert_eq!(r.histogram_overflow, 1);
    }
}
