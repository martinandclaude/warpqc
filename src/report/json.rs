//! The `<prefix>.qc.json` schema (v1.0). Fractions are 0-1 floats.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::metrics::align::AlignReport;
use crate::metrics::coverage::Exclusions;
use crate::metrics::insert::{InsertReport, OrientationReport};
use crate::summary::depth::DepthReport;

pub const SCHEMA_VERSION: &str = "1.1";

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Report {
    pub schema_version: String,
    pub tool: Tool,
    pub sample: String,
    pub inputs: Inputs,
    pub parameters: Parameters,
    pub definitions: BTreeMap<String, String>,
    pub warnings: Vec<String>,
    pub alignment: Alignment,
    pub insert_size: InsertSize,
    pub duplicates: Duplicates,
    pub coverage: CoverageSection,
    /// Absent in schema 1.0 reports.
    #[serde(default)]
    pub uniformity: crate::summary::uniformity::UniformityReport,
    pub per_lane: BTreeMap<String, Lane>,
    pub qc_status: QcStatus,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Tool {
    pub name: String,
    pub version: String,
    pub command: String,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Inputs {
    pub input: String,
    pub reference: String,
    /// MD5 over the ordered "name<TAB>length<TAB>M5" lines of the reference.
    pub reference_md5: String,
    pub lanes: Vec<String>,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Parameters {
    pub min_mapq: u8,
    pub min_bq: u8,
    pub optical_distance: u32,
    pub depth_thresholds: Vec<u32>,
    pub raw_depth: bool,
    pub library: String,
    pub analysis_contigs: Vec<String>,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Alignment {
    #[serde(flatten)]
    pub stats: AlignReport,
    pub incomplete_templates: u64,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct InsertSize {
    /// FR is the headline orientation.
    #[serde(flatten)]
    pub fr: OrientationReport,
    pub orientation_fractions: BTreeMap<String, f64>,
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
        .map(|(k, v)| (k.to_string(), v))
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

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
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

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct CoverageSection {
    pub filtered: DepthReport,
    pub raw: Option<DepthReport>,
    pub exclusions: ExclusionReport,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct ExclusionReport {
    pub mapq: f64,
    pub dupe: f64,
    pub baseq: f64,
    pub overlap: f64,
    pub capped: f64,
    pub total: f64,
    pub bases: BTreeMap<String, u64>,
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
            .map(|(k, v)| (k.to_string(), v))
            .collect(),
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct Lane {
    pub templates: u64,
    pub reads: u64,
    pub mapped_fraction: f64,
    pub pairs_examined: u64,
    /// Fractional: each member of a group of n counts (n - 1)/n here.
    pub duplicate_pairs: f64,
    pub duplicate_rate: f64,
    pub optical_duplicate_pairs: Option<u64>,
    pub optical_rate: Option<f64>,
    /// FR pairs, non-duplicate (absent in schema 1.0).
    #[serde(default)]
    pub insert_median: Option<f64>,
    /// As alignment.mismatch_rate, for this lane (absent in schema 1.0).
    #[serde(default)]
    pub mismatch_rate: Option<f64>,
}

pub use crate::report::thresholds::QcStatus;

impl Report {
    /// Reads a report written by any 1.x version.
    pub fn from_json(bytes: &[u8]) -> crate::error::Result<Report> {
        let r: Report = serde_json::from_slice(bytes)
            .map_err(|e| crate::error::input(format!("not a WarpQC report: {e}")))?;
        if !r.schema_version.starts_with("1.") {
            return Err(crate::error::input(format!(
                "report schema {} is not supported (expected 1.x)",
                r.schema_version
            )));
        }
        Ok(r)
    }

    pub fn to_json(&self) -> Vec<u8> {
        let mut b = serde_json::to_vec_pretty(self).expect("report serialises");
        b.push(b'\n');
        b
    }

    pub fn apply_thresholds(&mut self, t: &crate::report::thresholds::Thresholds) {
        let v = serde_json::to_value(&*self).expect("report serialises");
        self.qc_status = t.evaluate(&v);
    }

    /// A report with the shape of a real one for this configuration (all
    /// values zero or null), used to check threshold paths before a run.
    pub fn skeleton(thresholds: &[u32], contigs: &[String], raw: bool) -> serde_json::Value {
        let depth = DepthReport {
            pct_ge: thresholds.iter().map(|t| (t.to_string(), 0.0)).collect(),
            per_contig: contigs.iter().map(|c| (c.clone(), 0.0)).collect(),
            ..Default::default()
        };
        let r = Report {
            coverage: CoverageSection {
                filtered: depth.clone(),
                raw: raw.then_some(depth),
                ..Default::default()
            },
            ..Default::default()
        };
        serde_json::to_value(&r).expect("report serialises")
    }
}

pub fn definitions() -> BTreeMap<String, String> {
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
          orientation as htsjdk >= 5.0.0 / Picard >= 3.5.0 (forward start <= reverse end; Picard 3.4.0 differs on dovetailed pairs); \
          mean and sd over inserts <= median + 10 MAD (Picard's upper-only trim); percentiles as htsjdk getPercentile"),
        ("duplicates.signature",
         "library + per mate (contig, unclipped 5' position, strand), ends in canonical order (Picard)"),
        ("duplicates.kept_member",
         "the first-seen member of a group is kept for coverage and insert size; Picard keeps the highest base-quality sum"),
        ("duplicates.fragments",
         "fragments (mate unmapped or unpaired) are compared with fragments only; Picard also marks fragments at pair-end positions"),
        ("duplicates.optical",
         "same lane (Picard: read group) and tile, |dx| <= d and |dy| <= d, split by R1 orientation as Picard; single linkage, k - 1 per cluster. Picard's keeper-based pass can count slightly fewer in chains (0.17 % fewer on a synthetic test set)"),
        ("duplicates.optical_rate", "2 x optical duplicate pairs over the PERCENT_DUPLICATION denominator"),
        ("duplicates.cross_lane_fraction",
         "duplicate pairs in groups whose members span more than one flowcell:lane / all duplicate pairs"),
        ("per_lane.duplicates",
         "pairs only; each member of a duplicate group of n counts (n - 1)/n duplicates in its own lane (expected value under a random keeper), so no lane is favoured and lanes sum to the total"),
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
    .map(|(k, v)| (k.to_string(), v.to_string()))
    .collect()
}
