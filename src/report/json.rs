//! The `<prefix>.qc.json` schema (v1.0). Fractions are 0-1 floats.

use std::collections::BTreeMap;

use serde::Serialize;

use crate::metrics::align::AlignReport;
use crate::metrics::coverage::Exclusions;
use crate::metrics::insert::{InsertReport, OrientationReport};
use crate::summary::depth::DepthReport;

pub const SCHEMA_VERSION: &str = "1.0";

#[derive(Debug, Serialize)]
pub struct Report {
    pub schema_version: &'static str,
    pub tool: Tool,
    pub sample: String,
    pub inputs: Inputs,
    pub parameters: Parameters,
    pub definitions: BTreeMap<&'static str, &'static str>,
    pub warnings: Vec<String>,
    pub alignment: Alignment,
    pub insert_size: InsertSize,
    pub duplicates: Duplicates,
    pub coverage: CoverageSection,
    pub per_lane: BTreeMap<String, Lane>,
    pub qc_status: QcStatus,
}

#[derive(Debug, Serialize)]
pub struct Tool {
    pub name: &'static str,
    pub version: &'static str,
    pub command: String,
}

#[derive(Debug, Serialize)]
pub struct Inputs {
    pub input: String,
    pub reference: String,
    /// MD5 over the ordered "name<TAB>length<TAB>M5" lines of the reference.
    pub reference_md5: String,
    pub lanes: Vec<String>,
}

#[derive(Debug, Serialize)]
pub struct Parameters {
    pub min_mapq: u8,
    pub min_bq: u8,
    pub optical_distance: u32,
    pub depth_thresholds: Vec<u32>,
    pub raw_depth: bool,
    pub library: String,
    pub analysis_contigs: Vec<String>,
}

#[derive(Debug, Serialize)]
pub struct Alignment {
    #[serde(flatten)]
    pub stats: AlignReport,
    pub incomplete_templates: u64,
}

#[derive(Debug, Serialize)]
pub struct InsertSize {
    /// FR is the headline orientation.
    #[serde(flatten)]
    pub fr: OrientationReport,
    pub orientation_fractions: BTreeMap<&'static str, f64>,
    pub rf: OrientationReport,
    pub tandem: OrientationReport,
}

impl From<InsertReport> for InsertSize {
    fn from(r: InsertReport) -> InsertSize {
        let orientation_fractions = [
            ("fr", r.fr.fraction),
            ("rf", r.rf.fraction),
            ("tandem", r.tandem.fraction),
        ]
        .into_iter()
        .collect();
        let mut rf = r.rf;
        let mut tandem = r.tandem;
        rf.histogram.clear();
        tandem.histogram.clear();
        InsertSize {
            fr: r.fr,
            orientation_fractions,
            rf,
            tandem,
        }
    }
}

#[derive(Debug, Serialize)]
pub struct Duplicates {
    /// Picard PERCENT_DUPLICATION.
    pub rate: f64,
    pub optical_rate: Option<f64>,
    pub non_optical_rate: Option<f64>,
    pub cross_lane_fraction: f64,
    /// Picard ESTIMATED_LIBRARY_SIZE (optical duplicates excluded).
    pub est_library_size: Option<u64>,
    pub pairs_examined: u64,
    pub duplicate_pairs: u64,
    pub optical_duplicate_pairs: Option<u64>,
    pub fragments_examined: u64,
    pub duplicate_fragments: u64,
    pub duplicate_groups: u64,
    pub distinct_pair_signatures: u64,
}

#[derive(Debug, Serialize)]
pub struct CoverageSection {
    pub filtered: DepthReport,
    pub raw: Option<DepthReport>,
    pub exclusions: ExclusionReport,
}

#[derive(Debug, Serialize)]
pub struct ExclusionReport {
    pub mapq: f64,
    pub dupe: f64,
    pub baseq: f64,
    pub overlap: f64,
    pub capped: f64,
    pub total: f64,
    pub bases: BTreeMap<&'static str, u64>,
}

impl ExclusionReport {
    /// Denominator as Picard: counted bases plus every exclusion.
    pub fn new(e: &Exclusions, counted: u64) -> ExclusionReport {
        let excluded = e.mapq + e.dupe + e.baseq + e.overlap + e.capped;
        let denom = (counted + excluded) as f64;
        let f = |x: u64| if denom > 0.0 { x as f64 / denom } else { 0.0 };
        ExclusionReport {
            mapq: f(e.mapq),
            dupe: f(e.dupe),
            baseq: f(e.baseq),
            overlap: f(e.overlap),
            capped: f(e.capped),
            total: f(excluded),
            bases: [
                ("counted", counted),
                ("mapq", e.mapq),
                ("dupe", e.dupe),
                ("baseq", e.baseq),
                ("overlap", e.overlap),
                ("capped", e.capped),
            ]
            .into_iter()
            .collect(),
        }
    }
}

#[derive(Debug, Serialize)]
pub struct Lane {
    pub templates: u64,
    pub reads: u64,
    pub mapped_fraction: f64,
    pub pairs_examined: u64,
    pub duplicate_pairs: u64,
    pub duplicate_rate: f64,
    pub optical_duplicate_pairs: Option<u64>,
    pub optical_rate: Option<f64>,
}

#[derive(Debug, Serialize)]
pub struct QcStatus {
    pub overall: &'static str,
    pub checks: Vec<serde_json::Value>,
}

pub fn definitions() -> BTreeMap<&'static str, &'static str> {
    [
        ("alignment.chimeric_read_fraction",
         "primary reads with the mate on another contig and own MAPQ >= 5 (samtools flagstat), over primary reads with both mates mapped"),
        ("alignment.chimeric_pair_fraction",
         "templates with mates on different contigs and both MAPQ >= 5, over templates with both mates mapped"),
        ("alignment.mismatch_rate",
         "(NM - inserted - deleted bases) / M,=,X bases, primary mapped reads carrying NM"),
        ("alignment.indel_rate_per_100", "insertion and deletion events per 100 M,=,X bases, primary mapped reads"),
        ("alignment.off_target_fraction", "primary mapped reads on non-analysis contigs / primary mapped reads"),
        ("insert_size",
         "one count per pair from R2 |TLEN|: both mates mapped, same contig, non-duplicate, TLEN != 0 (Picard); \
          mean and sd over inserts <= median + 10 MAD (Picard's upper-only trim); percentiles as htsjdk getPercentile"),
        ("duplicates.signature",
         "library + per mate (contig, unclipped 5' position, strand), ends in canonical order (Picard)"),
        ("duplicates.kept_member",
         "the first-seen member of a group is kept for coverage and insert size; Picard keeps the highest base-quality sum"),
        ("duplicates.fragments",
         "fragments (mate unmapped or unpaired) are compared with fragments only; Picard also marks fragments at pair-end positions"),
        ("duplicates.optical",
         "Picard: same lane (read group) and tile, |dx| <= d and |dy| <= d, split by R1 orientation, single linkage, k - 1 per cluster"),
        ("duplicates.optical_rate", "2 x optical duplicate pairs over the PERCENT_DUPLICATION denominator"),
        ("duplicates.cross_lane_fraction",
         "duplicate pairs in groups whose members span more than one flowcell:lane / all duplicate pairs"),
        ("per_lane.duplicates",
         "pairs only; in each group the member with the lowest (lane name, tile, x, y) is the keeper, the others count in their own lane"),
        ("coverage.filtered",
         "primary, non-duplicate, not QC-fail, MAPQ >= min_mapq; M,=,X bases with BQ >= min_bq and not N; mate overlap counted once. \
          Supplementary records excluded (Picard CollectWgsMetrics includes them); reads with an unmapped mate included \
          (Picard's default COUNT_UNPAIRED=false excludes them)"),
        ("coverage.raw",
         "primary, not QC-fail, any MAPQ, duplicates included, no BQ filter, M,=,X bases, mate overlap counted once (mosdepth default with -F 2820)"),
        ("coverage.territory", "non-N reference positions of the analysis contigs; depth saturates at 65,535 (Picard caps at 250)"),
        ("coverage.fold80", "mean / 20th-percentile depth with zero-depth positions included (Picard); null if that percentile is 0"),
        ("coverage.exclusions",
         "per base at non-N analysis positions, first matching reason in Picard order (MAPQ, duplicate, base quality, overlap, cap), \
          over counted + excluded bases. Picard charges whole reads for MAPQ and duplicates"),
    ]
    .into_iter()
    .collect()
}
