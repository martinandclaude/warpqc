//! Per-template metric updates. Shared state lives in `Shared`; every
//! histogram and counter is thread-local in `Local` and merged at the end.

pub mod align;
pub mod bins;
pub mod coverage;
pub mod dups;
pub mod insert;
pub mod optical;

use std::collections::HashMap;

use parking_lot::Mutex;

use crate::error::{Result, input};
use crate::input::group::Template;
use crate::readname::parse_illumina;
use crate::record::{Aln, flag};
use crate::reference::{AnalysisSet, Reference};
use align::AlignStats;
use coverage::{Class, Coverage, Exclusions, Scratch};
use dups::{DupTable, M_LOCATED, M_R1_REVERSE, Member, PairKey, frag_key, pack_end};
use insert::InsertStats;

#[derive(Debug, Clone)]
pub struct Params {
    pub min_mapq: u8,
    pub min_bq: u8,
    pub optical_distance: u32,
    pub raw: bool,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            min_mapq: 20,
            min_bq: 20,
            optical_distance: 2500,
            raw: true,
        }
    }
}

/// Header contigs and where the analysis contigs sit in the depth arrays.
pub struct Layout {
    /// Header order (tid order).
    pub names: Vec<String>,
    pub lens: Vec<u64>,
    /// Offset into the depth arrays per tid; `u64::MAX` for other contigs.
    pub offset: Vec<u64>,
    /// Analysis tids, in reference order.
    pub analysis: Vec<u32>,
    pub total: u64,
    /// One bit per depth-array position, set where the reference is N.
    nmask: Vec<u64>,
}

impl Layout {
    pub fn new(reference: &Reference, set: &AnalysisSet, tid_to_ref: &[usize]) -> Result<Layout> {
        if tid_to_ref.len() > u16::MAX as usize + 1 {
            return Err(input("more than 65,536 @SQ contigs are not supported"));
        }
        let names: Vec<String> = tid_to_ref
            .iter()
            .map(|&r| reference.contigs[r].name.clone())
            .collect();
        let lens: Vec<u64> = tid_to_ref
            .iter()
            .map(|&r| reference.contigs[r].len)
            .collect();
        let mut offset = vec![u64::MAX; names.len()];
        let mut analysis = Vec::new();
        let mut total = 0u64;
        for &ri in &set.contigs {
            let tid = tid_to_ref
                .iter()
                .position(|&r| r == ri)
                .expect("header covers the reference");
            offset[tid] = total;
            analysis.push(tid as u32);
            total += lens[tid];
        }
        let mut nmask = vec![0u64; (total as usize).div_ceil(64)];
        for &tid in &analysis {
            let off = offset[tid as usize];
            for &(a, b) in &reference.contigs[tid_to_ref[tid as usize]].n_runs {
                for g in off + a..off + b {
                    nmask[(g / 64) as usize] |= 1 << (g % 64);
                }
            }
        }
        Ok(Layout {
            names,
            lens,
            offset,
            analysis,
            total,
            nmask,
        })
    }

    #[inline]
    pub fn is_analysis(&self, tid: i32) -> bool {
        tid >= 0 && self.offset[tid as usize] != u64::MAX
    }

    #[inline]
    pub fn is_n(&self, g: u64) -> bool {
        self.nmask[(g / 64) as usize] >> (g % 64) & 1 == 1
    }
}

pub const UNKNOWN_LANE: &[u8] = b"unknown";

/// Lane keys (`flowcell:lane`) in the order first seen; the id is the index.
#[derive(Default)]
pub struct Lanes(Mutex<Vec<Vec<u8>>>);

impl Lanes {
    fn id(&self, key: &[u8]) -> Result<u8> {
        let mut v = self.0.lock();
        if let Some(i) = v.iter().position(|k| k == key) {
            return Ok(i as u8);
        }
        if v.len() == 255 {
            return Err(input("more than 255 flowcell:lane keys in the read names"));
        }
        v.push(key.to_vec());
        Ok((v.len() - 1) as u8)
    }

    pub fn names(&self) -> Vec<String> {
        self.0
            .lock()
            .iter()
            .map(|k| String::from_utf8_lossy(k).into_owned())
            .collect()
    }
}

pub struct Shared {
    pub layout: Layout,
    pub params: Params,
    pub library: u8,
    pub lanes: Lanes,
    pub dups: DupTable,
    pub coverage: Coverage,
    pub bins: bins::BinLayout,
}

#[derive(Debug, Clone, Default)]
pub struct DupCounts {
    pub pairs_examined: u64,
    pub dup_pairs: u64,
    pub frags_examined: u64,
    pub dup_frags: u64,
}

#[derive(Debug, Clone, Default)]
pub struct LaneStats {
    pub templates: u64,
    pub reads: u64,
    pub mapped: u64,
    pub pairs_examined: u64,
}

#[derive(Default)]
pub struct Local {
    pub align: AlignStats,
    pub insert: InsertStats,
    pub dup: DupCounts,
    pub excl: Exclusions,
    pub lanes: Vec<LaneStats>,
    pub unparsed_names: u64,
    /// Read starts per bin (see `bins`); empty until the first count.
    pub bin_counts: Vec<u32>,
    lane_cache: HashMap<Vec<u8>, u8>,
    scratch: Scratch,
}

impl Local {
    fn lane_id(&mut self, key: &[u8], shared: &Shared) -> Result<u8> {
        let id = match self.lane_cache.get(key) {
            Some(&id) => id,
            None => {
                let id = shared.lanes.id(key)?;
                self.lane_cache.insert(key.to_vec(), id);
                id
            }
        };
        if self.lanes.len() <= id as usize {
            self.lanes.resize(id as usize + 1, LaneStats::default());
        }
        Ok(id)
    }

    pub fn merge(&mut self, o: &Local) {
        self.align.merge(&o.align);
        self.insert.merge(&o.insert);
        self.dup.pairs_examined += o.dup.pairs_examined;
        self.dup.dup_pairs += o.dup.dup_pairs;
        self.dup.frags_examined += o.dup.frags_examined;
        self.dup.dup_frags += o.dup.dup_frags;
        self.excl.merge(&o.excl);
        if self.lanes.len() < o.lanes.len() {
            self.lanes.resize(o.lanes.len(), LaneStats::default());
        }
        for (a, b) in self.lanes.iter_mut().zip(&o.lanes) {
            a.templates += b.templates;
            a.reads += b.reads;
            a.mapped += b.mapped;
            a.pairs_examined += b.pairs_examined;
        }
        self.unparsed_names += o.unparsed_names;
        if !o.bin_counts.is_empty() {
            if self.bin_counts.is_empty() {
                self.bin_counts = o.bin_counts.clone();
            } else {
                for (a, b) in self.bin_counts.iter_mut().zip(&o.bin_counts) {
                    *a += b;
                }
            }
        }
    }
}

fn end_of(a: &Aln) -> u64 {
    pack_end(a.tid as u16, a.unclipped_5p(), a.is_reverse())
}

/// Updates every metric for one template.
pub fn process(t: &Template, sh: &Shared, lo: &mut Local) -> Result<()> {
    lo.align.add_template(t, &sh.layout);

    let name = parse_illumina(t.qname);
    let lane = match &name {
        Some(n) => lo.lane_id(n.lane_key, sh)?,
        None => {
            lo.unparsed_names += 1;
            lo.lane_id(UNKNOWN_LANE, sh)?
        }
    };
    {
        let ls = &mut lo.lanes[lane as usize];
        ls.templates += 1;
        for a in t.primaries() {
            ls.reads += 1;
            if a.is_mapped() {
                ls.mapped += 1;
            }
        }
    }

    let r1 = t.get(t.r1);
    let r2 = t.get(t.r2);
    // Duplicate verdicts, indexed like `primaries`: r1, r2, single.
    let mut dup = [false; 3];
    match (r1, r2) {
        (Some(a), Some(b)) if a.is_mapped() && b.is_mapped() => {
            let key = PairKey::new(sh.library, end_of(a), end_of(b));
            let mut member = Member {
                lane,
                flags: if a.is_reverse() { M_R1_REVERSE } else { 0 },
                ..Default::default()
            };
            if let Some(n) = name.filter(|n| n.tile <= u16::MAX as u32) {
                member.flags |= M_LOCATED;
                member.tile = n.tile as u16;
                member.x = n.x;
                member.y = n.y;
            }
            let first = sh.dups.insert_pair(key, member);
            lo.dup.pairs_examined += 1;
            lo.lanes[lane as usize].pairs_examined += 1;
            if first {
                if a.tid == b.tid {
                    lo.insert.add(a, b);
                }
            } else {
                lo.dup.dup_pairs += 1;
                dup = [true, true, false];
            }
        }
        _ => {
            for (k, a) in [r1, r2, t.get(t.single)].into_iter().enumerate() {
                if let Some(a) = a.filter(|a| a.is_mapped()) {
                    lo.dup.frags_examined += 1;
                    if !sh.dups.insert_frag(frag_key(sh.library, end_of(a))) {
                        lo.dup.dup_frags += 1;
                        dup[k] = true;
                    }
                }
            }
        }
    }

    let mut reads: [(&Aln, Class); 2] = [(&t.alns[0], Class::Pass); 2];
    let mut n = 0;
    for (k, a) in [r1, r2, t.get(t.single)].into_iter().enumerate() {
        if let Some(a) = a.filter(|a| a.is_mapped() && !a.has(flag::QC_FAIL))
            && n < 2
        {
            let class = coverage::classify(a, dup[k], &sh.params);
            if class == Class::Pass
                && let Some(i) = sh.bins.index(a.tid, a.pos)
            {
                if lo.bin_counts.is_empty() {
                    lo.bin_counts = vec![0; sh.bins.len()];
                }
                lo.bin_counts[i] += 1;
            }
            reads[n] = (a, class);
            n += 1;
        }
    }
    coverage::add_template(
        &reads[..n],
        &sh.coverage,
        &sh.layout,
        &sh.params,
        &mut lo.excl,
        &mut lo.scratch,
    );
    Ok(())
}
