//! End-of-stream sweep over a depth array: depth histogram over non-N
//! analysis positions and per-contig sums, then the coverage summary.

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicU16, Ordering::Relaxed};

use rayon::prelude::*;
use serde::Serialize;

use crate::metrics::Layout;
use crate::reference::is_autosome;
use crate::summary::hist;

const CHUNK: u64 = 1 << 22;

pub struct DepthSweep {
    pub hist: Vec<u64>,
    /// Per analysis contig (in `layout.analysis` order): depth sum, territory.
    pub per_contig: Vec<(u64, u64)>,
}

pub fn sweep(depth: &[AtomicU16], layout: &Layout) -> DepthSweep {
    let mut tasks = Vec::new();
    for (k, &tid) in layout.analysis.iter().enumerate() {
        let off = layout.offset[tid as usize];
        let len = layout.lens[tid as usize];
        let mut s = 0;
        while s < len {
            tasks.push((k, off + s, off + (s + CHUNK).min(len)));
            s += CHUNK;
        }
    }
    let n = layout.analysis.len();
    let (hist, per_contig) = tasks
        .par_iter()
        .fold(
            || (vec![0u64; 1 << 16], vec![(0u64, 0u64); n]),
            |(mut h, mut pc), &(k, a, b)| {
                let (mut sum, mut terr) = (0u64, 0u64);
                for g in a..b {
                    if layout.is_n(g) {
                        continue;
                    }
                    let d = depth[g as usize].load(Relaxed);
                    h[d as usize] += 1;
                    sum += d as u64;
                    terr += 1;
                }
                pc[k].0 += sum;
                pc[k].1 += terr;
                (h, pc)
            },
        )
        .reduce(
            || (vec![0u64; 1 << 16], vec![(0u64, 0u64); n]),
            |(mut h1, mut p1), (h2, p2)| {
                for (a, b) in h1.iter_mut().zip(&h2) {
                    *a += b;
                }
                for (a, b) in p1.iter_mut().zip(&p2) {
                    a.0 += b.0;
                    a.1 += b.1;
                }
                (h1, p1)
            },
        );
    DepthSweep { hist, per_contig }
}

#[derive(Debug, Serialize)]
pub struct DepthReport {
    pub territory: u64,
    pub mean: f64,
    pub median: f64,
    pub sd: f64,
    pub mad: f64,
    /// Fraction of territory at depth >= threshold, keyed by threshold.
    pub pct_ge: BTreeMap<String, f64>,
    /// Mean / 20th-percentile depth, zero-depth positions included (Picard).
    /// Null when the 20th percentile is 0.
    pub fold80: Option<f64>,
    pub autosomal_mean: Option<f64>,
    pub per_contig: BTreeMap<String, f64>,
    /// Positions at each depth 0..=max observed.
    pub histogram: Vec<u64>,
}

pub fn report(s: &DepthSweep, layout: &Layout, thresholds: &[u32]) -> DepthReport {
    let h = &s.hist;
    let mean = hist::mean(h);
    let p20 = hist::percentile(h, 0.2);
    let mut per_contig = BTreeMap::new();
    let (mut asum, mut aterr) = (0u64, 0u64);
    for (k, &tid) in layout.analysis.iter().enumerate() {
        let name = &layout.names[tid as usize];
        let (sum, terr) = s.per_contig[k];
        per_contig.insert(name.clone(), if terr > 0 { sum as f64 / terr as f64 } else { 0.0 });
        if is_autosome(name) {
            asum += sum;
            aterr += terr;
        }
    }
    DepthReport {
        territory: hist::count(h),
        mean,
        median: hist::median(h),
        sd: hist::sd(h),
        mad: hist::mad(h),
        pct_ge: thresholds
            .iter()
            .map(|&t| (t.to_string(), hist::fraction_ge(h, t as usize)))
            .collect(),
        fold80: (p20 > 0).then(|| mean / p20 as f64),
        autosomal_mean: (aterr > 0).then(|| asum as f64 / aterr as f64),
        per_contig,
        histogram: h[..hist::used_len(h).max(1)].to_vec(),
    }
}
