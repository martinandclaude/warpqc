use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Args, Parser, Subcommand};

use warpqc::error::{Result, input};
use warpqc::metrics::Params;
use warpqc::pipeline::{DEFAULT_THRESHOLDS, RunConfig, run};
use warpqc::reference::{DEFAULT_BIN_SIZE, Reference};
use warpqc::report::{Outputs, write_all};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// WarpQC: whole-genome QC in one streaming pass over aligner output.
#[derive(Parser)]
#[command(name = "warpqc", version)]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
#[allow(clippy::large_enum_variant)] // parsed once at startup
enum Cmd {
    /// Compute QC metrics from name-grouped SAM/BAM and write the report files.
    Run(RunArgs),
    /// Scan a reference FASTA once and write the .warpqc cache.
    PrepareRef {
        /// Reference FASTA (with .fai).
        #[arg(long)]
        reference: PathBuf,
        /// Cache file to write.
        #[arg(long)]
        out: PathBuf,
        /// Bin size for the per-bin GC table (must match `run --bin-size`).
        #[arg(long, default_value_t = DEFAULT_BIN_SIZE)]
        bin_size: u32,
    },
    /// Rebuild the HTML and MultiQC files from a report (both, unless --html
    /// or --multiqc picks one), optionally re-evaluating qc_status with new
    /// thresholds (rewrites the JSON then).
    Render {
        /// Write only <prefix>.qc.html (with --multiqc: both).
        #[arg(long)]
        html: bool,
        /// Write only the MultiQC files (with --html: both).
        #[arg(long)]
        multiqc: bool,
        /// A <prefix>.qc.json written by `warpqc run`.
        #[arg(long)]
        json: PathBuf,
        /// Thresholds YAML to apply.
        #[arg(long)]
        thresholds: Option<PathBuf>,
        /// Output prefix [default: the JSON path without .qc.json].
        #[arg(long)]
        out: Option<PathBuf>,
    },
}

#[derive(Args)]
struct RunArgs {
    /// Sample name for the report.
    #[arg(long)]
    sample: String,
    /// Reference FASTA (with .fai); must match the @SQ lines.
    #[arg(long)]
    reference: PathBuf,
    /// Cache from `warpqc prepare-ref`, checked against the reference.
    #[arg(long)]
    ref_cache: Option<PathBuf>,
    /// Worker threads.
    #[arg(long, default_value_t = 8)]
    threads: usize,
    /// Analysis contigs: comma-separated list or a file with one per line
    /// [default: chr1-22, chrX, chrY, chrM].
    #[arg(long)]
    contigs: Option<String>,
    #[arg(long, default_value_t = 20)]
    min_mapq: u8,
    #[arg(long, default_value_t = 20)]
    min_bq: u8,
    /// Optical duplicate distance in pixels, per axis.
    #[arg(long, default_value_t = 2500)]
    optical_distance: u32,
    /// Do not keep the raw depth array (saves 2 bytes per analysis base).
    #[arg(long)]
    no_raw: bool,
    /// Depth thresholds for the "% of bases >= X" fields.
    #[arg(long, value_delimiter = ',')]
    depth_thresholds: Option<Vec<u32>>,
    /// Library name for the report [default: the sample name].
    #[arg(long)]
    library: Option<String>,
    /// Expected template (read pair) count; a mismatch fails the run.
    /// Also pre-sizes the duplicate table.
    #[arg(long)]
    expect_reads: Option<u64>,
    /// Warn/fail bounds per JSON path (YAML); sets qc_status.
    #[arg(long)]
    thresholds: Option<PathBuf>,
    /// Bin size for PMAD and the GC curve.
    #[arg(long, default_value_t = DEFAULT_BIN_SIZE)]
    bin_size: u32,
    /// BED of regions whose bins are left out of PMAD (blacklist, low mappability).
    #[arg(long)]
    exclude_bed: Option<PathBuf>,
    /// Also write <out>.qc.html.
    #[arg(long)]
    html: bool,
    /// Also write the MultiQC custom content (<out>.general_mqc.json,
    /// <out>.insert_mqc.json, <out>.coverage_mqc.json).
    #[arg(long)]
    multiqc: bool,
    /// Output prefix: writes <out>.qc.json, plus the --html / --multiqc files.
    #[arg(long)]
    out: PathBuf,
    /// Block size in bytes handed to each worker.
    #[arg(long, default_value_t = 4 << 20, hide = true)]
    block_size: usize,
    /// SAM or BAM path, or - for stdin.
    input: String,
}

fn run_cmd(a: RunArgs) -> Result<()> {
    if a.threads == 0 {
        return Err(input("--threads must be at least 1"));
    }
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(a.threads)
        .build_global();
    let mut cfg = RunConfig::new(
        &a.sample,
        a.reference,
        (a.input != "-").then(|| PathBuf::from(&a.input)),
    );
    cfg.ref_cache = a.ref_cache;
    cfg.contigs = a.contigs;
    cfg.threads = a.threads;
    cfg.params = Params {
        min_mapq: a.min_mapq,
        min_bq: a.min_bq,
        optical_distance: a.optical_distance,
        raw: !a.no_raw,
    };
    cfg.depth_thresholds = a
        .depth_thresholds
        .unwrap_or_else(|| DEFAULT_THRESHOLDS.to_vec());
    cfg.library = a.library.unwrap_or_default();
    cfg.expect_reads = a.expect_reads;
    cfg.block_size = a.block_size;
    cfg.thresholds = a.thresholds;
    cfg.bin_size = a.bin_size;
    cfg.exclude_bed = a.exclude_bed;
    cfg.command = std::env::args().collect::<Vec<_>>().join(" ");
    let out = run(&cfg)?;
    let outputs = Outputs {
        html: a.html,
        multiqc: a.multiqc,
    };
    write_all(&out.report, &a.out, outputs).map(|_| ())
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let res = match cli.cmd {
        Cmd::Run(a) => run_cmd(a),
        Cmd::PrepareRef {
            reference,
            out,
            bin_size,
        } => Reference::from_fasta(&reference, bin_size).and_then(|r| r.write_cache(&out)),
        Cmd::Render {
            json,
            thresholds,
            out,
            html,
            multiqc,
        } => {
            let outputs = if html || multiqc {
                Outputs { html, multiqc }
            } else {
                Outputs::ALL
            };
            warpqc::report::render(&json, thresholds.as_deref(), out.as_deref(), outputs)
                .map(|_| ())
        }
    };
    match res {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("warpqc: {e}");
            ExitCode::from(e.exit_code() as u8)
        }
    }
}
