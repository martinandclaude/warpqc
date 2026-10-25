//! `warpqc run`: reader thread -> N workers -> merge -> final sweep.

use std::collections::BTreeMap;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};

use crate::error::{IoContext, Result, input};
use crate::input::group::Grouper;
use crate::input::header::Header;
use crate::input::sam::{Blocker, LineSource, OrderStats, SamSource};
use crate::metrics::bins::BinLayout;
use crate::metrics::coverage::Coverage;
use crate::metrics::dups::DupTable;
use crate::metrics::optical::{DupAnalysis, analyze_group};
use crate::metrics::{Lanes, Layout, Local, Params, Shared, UNKNOWN_LANE, process};
use crate::record::TidMap;
use crate::reference::{AnalysisSet, Reference};
use crate::report::json::*;
use crate::summary::{depth, library, uniformity};

pub const DEFAULT_THRESHOLDS: [u32; 8] = [1, 5, 10, 15, 20, 30, 50, 100];
/// With more unparseable read names than this fraction, optical duplicate
/// detection is reported as disabled.
const MAX_UNPARSED_FOR_OPTICAL: f64 = 0.5;

#[derive(Debug, Clone)]
pub struct RunConfig {
    pub sample: String,
    /// `None` reads stdin.
    pub input: Option<PathBuf>,
    pub reference: PathBuf,
    pub ref_cache: Option<PathBuf>,
    pub contigs: Option<String>,
    pub threads: usize,
    pub params: Params,
    pub depth_thresholds: Vec<u32>,
    pub library: String,
    pub expect_reads: Option<u64>,
    pub block_size: usize,
    pub bin_size: u32,
    pub exclude_bed: Option<PathBuf>,
    pub thresholds: Option<PathBuf>,
    pub command: String,
    /// Keep copies of the depth arrays in the outcome (tests).
    pub keep_depth: bool,
}

impl RunConfig {
    pub fn new(sample: &str, reference: PathBuf, input: Option<PathBuf>) -> RunConfig {
        RunConfig {
            sample: sample.into(),
            input,
            reference,
            ref_cache: None,
            contigs: None,
            threads: 4,
            params: Params::default(),
            depth_thresholds: DEFAULT_THRESHOLDS.to_vec(),
            library: String::new(),
            expect_reads: None,
            block_size: 4 << 20,
            bin_size: crate::reference::DEFAULT_BIN_SIZE,
            exclude_bed: None,
            thresholds: None,
            command: String::new(),
            keep_depth: false,
        }
    }
}

pub struct Outcome {
    pub report: Report,
    /// Header tid order, analysis contigs only: (name, filtered, raw).
    pub depth: Vec<(String, Vec<u16>, Option<Vec<u16>>)>,
    /// Binned contigs: (name, read-start counts, usable).
    pub bins: Vec<(String, Vec<u32>, Vec<bool>)>,
}

impl<S: LineSource + ?Sized> LineSource for Box<S> {
    fn header(&self) -> &str {
        (**self).header()
    }
    fn next_line(&mut self, buf: &mut Vec<u8>) -> Result<bool> {
        (**self).next_line(buf)
    }
}

fn open_input(path: Option<&PathBuf>) -> Result<Box<dyn LineSource>> {
    let (mut r, name): (Box<dyn BufRead + Send>, String) = match path {
        None => (
            Box::new(BufReader::with_capacity(1 << 22, std::io::stdin())),
            "stdin".into(),
        ),
        Some(p) => (
            Box::new(BufReader::with_capacity(
                1 << 22,
                std::fs::File::open(p).ctx(format!("opening {}", p.display()))?,
            )),
            p.display().to_string(),
        ),
    };
    let head = r.fill_buf().ctx(format!("reading {name}"))?;
    if head.starts_with(&[0x1f, 0x8b]) {
        return crate::input::bam::open(r);
    }
    Ok(Box::new(SamSource::new(r)?))
}

pub fn run(cfg: &RunConfig) -> Result<Outcome> {
    let (reference, ref_warning) =
        Reference::load(&cfg.reference, cfg.ref_cache.as_deref(), cfg.bin_size)?;
    let set = AnalysisSet::resolve(&reference, cfg.contigs.as_deref())?;
    let src = open_input(cfg.input.as_ref())?;
    let header = Header::parse(src.header())?;
    let tid_to_ref = header.check(&reference)?;
    let layout = Layout::new(&reference, &set, &tid_to_ref)?;
    let tids = TidMap::new(&header.names());

    let mut warnings: Vec<String> = ref_warning.into_iter().collect();
    let mut bins = BinLayout::new(&layout, &reference, &tid_to_ref);
    if let Some(bed) = &cfg.exclude_bed
        && let Some(w) = bins.exclude_bed(bed, &layout)?
    {
        warnings.push(w);
    }
    if !set.missing_defaults.is_empty() {
        warnings.push(format!(
            "default analysis contigs not in the reference: {}",
            set.missing_defaults.join(", ")
        ));
    }

    let shared = Shared {
        coverage: Coverage::new(&layout, cfg.params.raw),
        layout,
        params: cfg.params.clone(),
        library: 0,
        lanes: Lanes::default(),
        dups: DupTable::new(cfg.expect_reads),
        bins,
    };

    let (order, local) = stream(src, &shared, &tids, cfg)?;

    if let Some(expect) = cfg.expect_reads
        && order.templates != expect
    {
        return Err(input(format!(
            "read {} templates but --expect-reads is {expect}; the stream was truncated or the count is wrong",
            order.templates
        )));
    }

    finish(shared, local, order, &reference, cfg, warnings)
}

fn stream(
    src: Box<dyn LineSource>,
    sh: &Shared,
    tids: &TidMap,
    cfg: &RunConfig,
) -> Result<(OrderStats, Local)> {
    let mut blocker = Blocker::new(src, cfg.block_size);
    let abort = AtomicBool::new(false);
    let threads = cfg.threads.max(1);
    let (tx, rx) = crossbeam_channel::bounded(threads * 2);
    std::thread::scope(|s| {
        let abort = &abort;
        let reader = s.spawn(move || -> Result<OrderStats> {
            let res = (|| {
                while !abort.load(Ordering::Relaxed) {
                    match blocker.next_block()? {
                        Some(b) => {
                            if tx.send(b).is_err() {
                                break;
                            }
                        }
                        None => break,
                    }
                }
                Ok(blocker.stats)
            })();
            if res.is_err() {
                abort.store(true, Ordering::Relaxed);
            }
            res
        });
        let workers: Vec<_> = (0..threads)
            .map(|_| {
                let rx = rx.clone();
                let mut tids = tids.clone();
                s.spawn(move || -> Result<Local> {
                    let mut lo = Local::default();
                    let mut g = Grouper::default();
                    for block in rx {
                        if abort.load(Ordering::Relaxed) {
                            continue;
                        }
                        if let Err(e) =
                            g.for_each_template(&block.data, &mut tids, |t| process(t, sh, &mut lo))
                        {
                            abort.store(true, Ordering::Relaxed);
                            return Err(e);
                        }
                    }
                    Ok(lo)
                })
            })
            .collect();
        drop(rx);
        let order = reader.join().expect("reader thread panicked");
        let mut merged = Local::default();
        let mut worker_err = None;
        for w in workers {
            match w.join().expect("worker thread panicked") {
                Ok(lo) => merged.merge(&lo),
                Err(e) => {
                    worker_err.get_or_insert(e);
                }
            }
        }
        let order = order?;
        if let Some(e) = worker_err {
            return Err(e);
        }
        Ok((order, merged))
    })
}

fn rate(n: u64, d: u64) -> f64 {
    if d == 0 { 0.0 } else { n as f64 / d as f64 }
}

fn finish(
    mut sh: Shared,
    lo: Local,
    order: OrderStats,
    reference: &Reference,
    cfg: &RunConfig,
    mut warnings: Vec<String>,
) -> Result<Outcome> {
    let lane_names = sh.lanes.names();
    let nl = lane_names.len();
    let d = cfg.params.optical_distance;
    let parts = sh
        .dups
        .for_each_group(|| DupAnalysis::new(nl), |acc, m| analyze_group(m, d, acc));
    let mut da = DupAnalysis::new(nl);
    for p in &parts {
        da.merge(p);
    }
    debug_assert_eq!(da.dup_pairs, lo.dup.dup_pairs);

    if order.incomplete > 0 {
        warnings.push(format!(
            "{} of {} templates are missing a primary mate (tolerated below {} %)",
            order.incomplete,
            order.templates,
            crate::input::sam::MAX_INCOMPLETE_FRACTION * 100.0
        ));
    }
    let optical_ok =
        (lo.unparsed_names as f64) <= MAX_UNPARSED_FOR_OPTICAL * lo.align.templates as f64;
    if lo.unparsed_names > 0 {
        warnings.push(format!(
            "{} of {} read names are not Illumina instrument:run:flowcell:lane:tile:x:y; they are pooled in lane \"{}\" \
             and excluded from optical duplicate detection{}",
            lo.unparsed_names,
            lo.align.templates,
            String::from_utf8_lossy(UNKNOWN_LANE),
            if optical_ok { "" } else { "; optical duplicate detection is disabled" }
        ));
    }

    let dc = &lo.dup;
    let denom = 2 * dc.pairs_examined + dc.frags_examined;
    let dup_rate = rate(2 * dc.dup_pairs + dc.dup_frags, denom);
    let optical_rate = optical_ok.then(|| rate(2 * da.optical, denom));
    let duplicates = Duplicates {
        rate: dup_rate,
        optical_rate,
        non_optical_rate: optical_rate.map(|o| dup_rate - o),
        cross_lane_fraction: rate(da.cross_lane_dup_pairs, da.dup_pairs),
        est_library_size: library::estimate_library_size(
            dc.pairs_examined - if optical_ok { da.optical } else { 0 },
            dc.pairs_examined - dc.dup_pairs,
        ),
        pairs_examined: dc.pairs_examined,
        duplicate_pairs: dc.dup_pairs,
        optical_duplicate_pairs: optical_ok.then_some(da.optical),
        fragments_examined: dc.frags_examined,
        duplicate_fragments: dc.dup_frags,
        duplicate_groups: da.groups,
        distinct_pair_signatures: sh.dups.distinct_pairs(),
    };

    let mut per_lane = BTreeMap::new();
    for (i, name) in lane_names.iter().enumerate() {
        let ls = lo.lanes.get(i).cloned().unwrap_or_default();
        per_lane.insert(
            name.clone(),
            Lane {
                templates: ls.templates,
                reads: ls.reads,
                mapped_fraction: rate(ls.mapped, ls.reads),
                pairs_examined: ls.pairs_examined,
                duplicate_pairs: da.lane_dups(i),
                duplicate_rate: if ls.pairs_examined > 0 {
                    da.lane_dups(i) / ls.pairs_examined as f64
                } else {
                    0.0
                },
                optical_duplicate_pairs: optical_ok.then_some(da.lane_optical[i]),
                optical_rate: optical_ok.then(|| rate(da.lane_optical[i], ls.pairs_examined)),
            },
        );
    }

    let layout = &sh.layout;
    let fsweep = depth::sweep(sh.coverage.filtered.as_slice(), layout);
    let counted: u64 = fsweep.per_contig.iter().map(|p| p.0).sum();
    let filtered = depth::report(&fsweep, layout, &cfg.depth_thresholds);
    let raw = sh.coverage.raw.as_ref().map(|r| {
        depth::report(
            &depth::sweep(r.as_slice(), layout),
            layout,
            &cfg.depth_thresholds,
        )
    });

    let uniformity = uniformity::report(
        &sh.bins,
        &lo.bin_counts,
        &fsweep.hist,
        uniformity::MIN_STRATUM_BINS,
    );

    let mut depth_copy = Vec::new();
    let mut bins_copy = Vec::new();
    if cfg.keep_depth {
        for c in &sh.bins.contigs {
            let r = c.first..c.first + c.n;
            bins_copy.push((
                c.name.clone(),
                if lo.bin_counts.is_empty() {
                    vec![0; c.n]
                } else {
                    lo.bin_counts[r.clone()].to_vec()
                },
                r.map(|i| sh.bins.usable(i)).collect(),
            ));
        }
        let f = sh.coverage.filtered.snapshot();
        let r = sh.coverage.raw.as_ref().map(|r| r.snapshot());
        for &tid in &layout.analysis {
            let (o, l) = (
                layout.offset[tid as usize] as usize,
                layout.lens[tid as usize] as usize,
            );
            depth_copy.push((
                layout.names[tid as usize].clone(),
                f[o..o + l].to_vec(),
                r.as_ref().map(|r| r[o..o + l].to_vec()),
            ));
        }
    }

    let mut lanes_sorted = lane_names.clone();
    lanes_sorted.sort();
    let report = Report {
        schema_version: SCHEMA_VERSION,
        tool: Tool {
            name: "warpqc",
            version: env!("CARGO_PKG_VERSION"),
            command: cfg.command.clone(),
        },
        sample: cfg.sample.clone(),
        inputs: Inputs {
            input: cfg
                .input
                .as_ref()
                .map_or("-".into(), |p| p.display().to_string()),
            reference: cfg.reference.display().to_string(),
            reference_md5: reference.digest(),
            lanes: lanes_sorted,
        },
        parameters: Parameters {
            min_mapq: cfg.params.min_mapq,
            min_bq: cfg.params.min_bq,
            optical_distance: cfg.params.optical_distance,
            depth_thresholds: cfg.depth_thresholds.clone(),
            raw_depth: cfg.params.raw,
            library: if cfg.library.is_empty() {
                cfg.sample.clone()
            } else {
                cfg.library.clone()
            },
            analysis_contigs: layout
                .analysis
                .iter()
                .map(|&t| layout.names[t as usize].clone())
                .collect(),
        },
        definitions: definitions(),
        warnings,
        alignment: Alignment {
            stats: lo.align.report(),
            incomplete_templates: order.incomplete,
        },
        insert_size: lo.insert.report().into(),
        duplicates,
        coverage: CoverageSection {
            filtered,
            raw,
            exclusions: ExclusionReport::new(&lo.excl, counted),
        },
        uniformity,
        per_lane,
        qc_status: QcStatus {
            overall: "not_evaluated",
            checks: vec![],
        },
    };
    Ok(Outcome {
        report,
        depth: depth_copy,
        bins: bins_copy,
    })
}
