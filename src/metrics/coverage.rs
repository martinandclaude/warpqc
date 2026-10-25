//! Per-base depth arrays and the exclusion counters.
//!
//! Filtered depth: primary, non-duplicate, not QC-fail, MAPQ >= min_mapq;
//! M/=/X bases with BQ >= min_bq (N read bases have BQ 0). A position covered
//! by both mates counts once if either mate's base passes.
//! Raw depth: primary, not QC-fail, any MAPQ, duplicates included, no BQ
//! filter; mate overlap counted once (mosdepth default with -F 2820).

use std::sync::atomic::{AtomicU16, Ordering::Relaxed};

use crate::metrics::{Layout, Params};
use crate::record::Aln;

pub struct Depth(Box<[AtomicU16]>);

impl Depth {
    pub fn zeroed(n: usize) -> Depth {
        // vec![0; n] comes from calloc, so untouched pages stay unmapped.
        let b: Box<[u16]> = vec![0u16; n].into_boxed_slice();
        let len = b.len();
        let p = Box::into_raw(b) as *mut AtomicU16;
        // SAFETY: AtomicU16 has the same size, alignment and bit validity as
        // u16, and the allocation came from a Box<[u16]> of the same length.
        Depth(unsafe { Box::from_raw(std::ptr::slice_from_raw_parts_mut(p, len)) })
    }

    /// Saturating increment; returns false if the position was already at
    /// 65,535. A compare-exchange loop, because a check followed by a
    /// fetch_add can wrap when threads race.
    #[inline]
    pub fn inc(&self, i: u64) -> bool {
        let a = &self.0[i as usize];
        let mut v = a.load(Relaxed);
        loop {
            if v == u16::MAX {
                return false;
            }
            match a.compare_exchange_weak(v, v + 1, Relaxed, Relaxed) {
                Ok(_) => return true,
                Err(now) => v = now,
            }
        }
    }

    pub fn as_slice(&self) -> &[AtomicU16] {
        &self.0
    }

    pub fn snapshot(&self) -> Vec<u16> {
        self.0.iter().map(|a| a.load(Relaxed)).collect()
    }
}

pub struct Coverage {
    pub filtered: Depth,
    pub raw: Option<Depth>,
}

impl Coverage {
    pub fn new(layout: &Layout, raw: bool) -> Coverage {
        let n = layout.total as usize;
        Coverage {
            filtered: Depth::zeroed(n),
            raw: raw.then(|| Depth::zeroed(n)),
        }
    }
}

/// Bases lost from filtered depth, counted at non-N analysis positions, in
/// Picard's order: MAPQ, then duplicate, then base quality, then overlap,
/// then the depth cap.
#[derive(Debug, Clone, Default)]
pub struct Exclusions {
    pub mapq: u64,
    pub dupe: u64,
    pub baseq: u64,
    pub overlap: u64,
    pub capped: u64,
}

impl Exclusions {
    pub fn merge(&mut self, o: &Exclusions) {
        self.mapq += o.mapq;
        self.dupe += o.dupe;
        self.baseq += o.baseq;
        self.overlap += o.overlap;
        self.capped += o.capped;
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Class {
    Pass,
    LowMapq,
    Duplicate,
}

pub fn classify(a: &Aln, duplicate: bool, params: &Params) -> Class {
    if a.mapq < params.min_mapq {
        Class::LowMapq
    } else if duplicate {
        Class::Duplicate
    } else {
        Class::Pass
    }
}

const RAW: u8 = 1;
const PASS: [u8; 2] = [2, 4];

/// Reusable per-worker scratch space for mate-overlap resolution.
#[derive(Default)]
pub struct Scratch {
    buf: Vec<u8>,
}

/// Adds one template's primary, non-QC-fail, mapped reads (at most two).
pub fn add_template(
    reads: &[(&Aln, Class)],
    cov: &Coverage,
    layout: &Layout,
    params: &Params,
    ex: &mut Exclusions,
    scratch: &mut Scratch,
) {
    debug_assert!(reads.len() <= 2);
    match reads {
        [(a, ca), (b, cb)]
            if layout.is_analysis(a.tid)
                && a.tid == b.tid
                && a.pos < b.end()
                && b.pos < a.end() =>
        {
            add_overlapping([(a, *ca), (b, *cb)], cov, layout, params, ex, scratch)
        }
        rs => {
            for (a, c) in rs {
                if layout.is_analysis(a.tid) {
                    add_single(a, *c, cov, layout, params, ex);
                }
            }
        }
    }
}

fn add_single(
    a: &Aln,
    class: Class,
    cov: &Coverage,
    layout: &Layout,
    params: &Params,
    ex: &mut Exclusions,
) {
    let off = layout.offset[a.tid as usize];
    let len = layout.lens[a.tid as usize] as i64;
    a.for_each_match_block(|r, q, l| {
        for i in 0..l as i64 {
            let pos = r + i;
            if pos < 0 || pos >= len {
                continue;
            }
            let g = off + pos as u64;
            if let Some(raw) = &cov.raw {
                raw.inc(g);
            }
            if layout.is_n(g) {
                continue;
            }
            match class {
                Class::LowMapq => ex.mapq += 1,
                Class::Duplicate => ex.dupe += 1,
                Class::Pass => {
                    if a.bq(q + i as usize) < params.min_bq {
                        ex.baseq += 1;
                    } else if !cov.filtered.inc(g) {
                        ex.capped += 1;
                    }
                }
            }
        }
    });
}

fn add_overlapping(
    reads: [(&Aln, Class); 2],
    cov: &Coverage,
    layout: &Layout,
    params: &Params,
    ex: &mut Exclusions,
    scratch: &mut Scratch,
) {
    let tid = reads[0].0.tid as usize;
    let off = layout.offset[tid];
    let len = layout.lens[tid] as i64;
    let start = reads[0].0.pos.min(reads[1].0.pos).max(0);
    let end = reads[0].0.end().max(reads[1].0.end()).min(len);
    if end <= start {
        return;
    }
    let buf = &mut scratch.buf;
    buf.clear();
    buf.resize((end - start) as usize, 0);
    for (k, (a, class)) in reads.iter().enumerate() {
        a.for_each_match_block(|r, q, l| {
            for i in 0..l as i64 {
                let pos = r + i;
                if pos < start || pos >= end {
                    continue;
                }
                let idx = (pos - start) as usize;
                buf[idx] |= RAW;
                if layout.is_n(off + pos as u64) {
                    continue;
                }
                match class {
                    Class::LowMapq => ex.mapq += 1,
                    Class::Duplicate => ex.dupe += 1,
                    Class::Pass => {
                        if a.bq(q + i as usize) < params.min_bq {
                            ex.baseq += 1;
                        } else {
                            buf[idx] |= PASS[k];
                        }
                    }
                }
            }
        });
    }
    for (idx, &s) in buf.iter().enumerate() {
        if s == 0 {
            continue;
        }
        let g = off + start as u64 + idx as u64;
        if let Some(raw) = &cov.raw {
            raw.inc(g);
        }
        let passes = (s & PASS[0] != 0) as u8 + (s & PASS[1] != 0) as u8;
        if passes > 0 && !cov.filtered.inc(g) {
            ex.capped += 1;
        }
        if passes == 2 {
            ex.overlap += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn saturates_at_u16_max() {
        let d = Depth::zeroed(2);
        for _ in 0..65_535 {
            assert!(d.inc(0));
        }
        assert!(!d.inc(0));
        assert!(!d.inc(0));
        assert_eq!(d.snapshot(), vec![65_535, 0]);
    }

    #[test]
    fn saturates_under_contention() {
        let d = Depth::zeroed(1);
        let ok = std::sync::atomic::AtomicU64::new(0);
        std::thread::scope(|s| {
            for _ in 0..8 {
                s.spawn(|| {
                    for _ in 0..10_000 {
                        if d.inc(0) {
                            ok.fetch_add(1, Relaxed);
                        }
                    }
                });
            }
        });
        assert_eq!(d.snapshot(), vec![65_535]);
        assert_eq!(ok.load(Relaxed), 65_535);
    }
}
