//! Alignment stats: flagstat-style counts, MAPQ, clipping, mismatches,
//! indels, chimeras and off-target contigs.

use serde::{Deserialize, Serialize};

use crate::input::group::Template;
use crate::metrics::Layout;
use crate::record::{Aln, flag, op};

pub const MAPQ_BINS: usize = 61;

#[derive(Debug, Clone, Default)]
pub struct AlignStats {
    pub templates: u64,
    pub records: u64,
    pub primary_reads: u64,
    pub secondary: u64,
    pub supplementary: u64,
    pub qc_fail: u64,
    pub primary_mapped: u64,
    pub paired_reads: u64,
    pub both_mapped_reads: u64,
    pub properly_paired: u64,
    pub singletons: u64,
    /// Primary reads, mate mapped to a different contig (samtools `n_diffchr`).
    pub mate_diff_chr: u64,
    /// ... and the read's own MAPQ >= 5 (samtools `n_diffhigh`).
    pub mate_diff_chr_q5: u64,
    /// Templates with both primary mates mapped.
    pub both_mapped_pairs: u64,
    /// ... on different contigs, both mates MAPQ >= 5.
    pub chimeric_pairs_q5: u64,
    /// Primary mapped reads, MAPQ 0..=60; index 60 also holds > 60.
    pub mapq: Vec<u64>,
    pub off_target_reads: u64,
    pub query_bases: u64,
    pub soft_clipped_bases: u64,
    pub aligned_bases: u64,
    pub insertion_events: u64,
    pub deletion_events: u64,
    /// Aligned bases of reads that carry NM, and their mismatches.
    pub nm_aligned_bases: u64,
    pub mismatches: u64,
    pub reads_without_nm: u64,
}

impl AlignStats {
    pub fn merge(&mut self, o: &AlignStats) {
        macro_rules! add { ($($f:ident),*) => { $( self.$f += o.$f; )* } }
        add!(
            templates,
            records,
            primary_reads,
            secondary,
            supplementary,
            qc_fail,
            primary_mapped,
            paired_reads,
            both_mapped_reads,
            properly_paired,
            singletons,
            mate_diff_chr,
            mate_diff_chr_q5,
            both_mapped_pairs,
            chimeric_pairs_q5,
            off_target_reads,
            query_bases,
            soft_clipped_bases,
            aligned_bases,
            insertion_events,
            deletion_events,
            nm_aligned_bases,
            mismatches,
            reads_without_nm
        );
        self.mapq.resize(MAPQ_BINS, 0);
        for (a, b) in self.mapq.iter_mut().zip(o.mapq.iter()) {
            *a += b;
        }
    }

    fn add_primary(&mut self, a: &Aln, layout: &Layout) {
        self.primary_reads += 1;
        if a.has(flag::QC_FAIL) {
            self.qc_fail += 1;
        }
        let mapped = a.is_mapped();
        if a.has(flag::PAIRED) {
            self.paired_reads += 1;
            if mapped && a.has(flag::PROPER_PAIR) {
                self.properly_paired += 1;
            }
            if mapped && a.has(flag::MATE_UNMAPPED) {
                self.singletons += 1;
            }
            if mapped && !a.has(flag::MATE_UNMAPPED) {
                self.both_mapped_reads += 1;
                if a.mate_tid != a.tid {
                    self.mate_diff_chr += 1;
                    if a.mapq >= 5 {
                        self.mate_diff_chr_q5 += 1;
                    }
                }
            }
        }
        if !mapped {
            return;
        }
        self.primary_mapped += 1;
        self.mapq.resize(MAPQ_BINS, 0);
        self.mapq[(a.mapq as usize).min(MAPQ_BINS - 1)] += 1;
        if !layout.is_analysis(a.tid) {
            self.off_target_reads += 1;
        }
        let (mut query, mut soft, mut aligned, mut ins_b, mut del_b) =
            (0u64, 0u64, 0u64, 0u64, 0u64);
        for &(o, l) in &a.cigar {
            let l = l as u64;
            if op::consumes_query(o) {
                query += l;
            }
            match o {
                op::S => soft += l,
                op::I => {
                    ins_b += l;
                    self.insertion_events += 1;
                }
                op::D => {
                    del_b += l;
                    self.deletion_events += 1;
                }
                o if op::is_match(o) => aligned += l,
                _ => {}
            }
        }
        self.query_bases += query;
        self.soft_clipped_bases += soft;
        self.aligned_bases += aligned;
        match a.nm {
            Some(nm) => {
                self.nm_aligned_bases += aligned;
                self.mismatches += (nm as u64).saturating_sub(ins_b + del_b);
            }
            None => self.reads_without_nm += 1,
        }
    }

    pub fn add_template(&mut self, t: &Template, layout: &Layout) {
        self.templates += 1;
        for a in t.alns {
            self.records += 1;
            if a.has(flag::SECONDARY) {
                self.secondary += 1;
            } else if a.has(flag::SUPPLEMENTARY) {
                self.supplementary += 1;
            } else {
                self.add_primary(a, layout);
            }
        }
        if let (Some(r1), Some(r2)) = (t.get(t.r1), t.get(t.r2))
            && r1.is_mapped()
            && r2.is_mapped()
        {
            self.both_mapped_pairs += 1;
            if r1.tid != r2.tid && r1.mapq >= 5 && r2.mapq >= 5 {
                self.chimeric_pairs_q5 += 1;
            }
        }
    }
}

fn frac(n: u64, d: u64) -> f64 {
    if d == 0 { 0.0 } else { n as f64 / d as f64 }
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct AlignReport {
    pub total_templates: u64,
    pub total_records: u64,
    pub total_reads: u64,
    pub primary_mapped: u64,
    pub primary_mapped_fraction: f64,
    pub paired_reads: u64,
    pub both_mates_mapped_fraction: f64,
    pub properly_paired_fraction: f64,
    pub singleton_fraction: f64,
    pub secondary_records: u64,
    pub supplementary_records: u64,
    pub supplementary_fraction: f64,
    pub qc_fail_reads: u64,
    pub mapq_histogram: Vec<u64>,
    pub mapq0_fraction: f64,
    pub mapq_ge20_fraction: f64,
    pub mapq_ge60_fraction: f64,
    pub mate_diff_contig_mapq5_reads: u64,
    pub chimeric_read_fraction: f64,
    pub chimeric_pairs_mapq5: u64,
    pub chimeric_pair_fraction: f64,
    pub off_target_fraction: f64,
    pub soft_clipped_fraction: f64,
    pub mismatch_rate: Option<f64>,
    pub indel_rate_per_100: f64,
    pub aligned_bases: u64,
    pub reads_without_nm: u64,
}

impl AlignStats {
    pub fn report(&self) -> AlignReport {
        let mapped = self.primary_mapped;
        let mut mq = self.mapq.clone();
        mq.resize(MAPQ_BINS, 0);
        let ge = |t: usize| mq[t..].iter().sum::<u64>();
        AlignReport {
            total_templates: self.templates,
            total_records: self.records,
            total_reads: self.primary_reads,
            primary_mapped: mapped,
            primary_mapped_fraction: frac(mapped, self.primary_reads),
            paired_reads: self.paired_reads,
            both_mates_mapped_fraction: frac(self.both_mapped_reads, self.paired_reads),
            properly_paired_fraction: frac(self.properly_paired, self.paired_reads),
            singleton_fraction: frac(self.singletons, self.paired_reads),
            secondary_records: self.secondary,
            supplementary_records: self.supplementary,
            supplementary_fraction: frac(self.supplementary, self.primary_reads),
            qc_fail_reads: self.qc_fail,
            mapq_histogram: mq.clone(),
            mapq0_fraction: frac(mq[0], mapped),
            mapq_ge20_fraction: frac(ge(20), mapped),
            mapq_ge60_fraction: frac(ge(60), mapped),
            mate_diff_contig_mapq5_reads: self.mate_diff_chr_q5,
            chimeric_read_fraction: frac(self.mate_diff_chr_q5, self.both_mapped_reads),
            chimeric_pairs_mapq5: self.chimeric_pairs_q5,
            chimeric_pair_fraction: frac(self.chimeric_pairs_q5, self.both_mapped_pairs),
            off_target_fraction: frac(self.off_target_reads, mapped),
            soft_clipped_fraction: frac(self.soft_clipped_bases, self.query_bases),
            mismatch_rate: (self.nm_aligned_bases > 0)
                .then(|| frac(self.mismatches, self.nm_aligned_bases)),
            indel_rate_per_100: 100.0
                * frac(
                    self.insertion_events + self.deletion_events,
                    self.aligned_bases,
                ),
            aligned_bases: self.aligned_bases,
            reads_without_nm: self.reads_without_nm,
        }
    }
}
