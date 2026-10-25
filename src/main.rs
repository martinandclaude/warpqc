use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Args, Parser, Subcommand};

use warpqc::error::{IoContext, Result, input};
use warpqc::metrics::Params;
use warpqc::pipeline::{DEFAULT_THRESHOLDS, RunConfig, run};
use warpqc::reference::Reference;
use warpqc::report::write_atomic;

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
enum Cmd {
    /// Compute QC metrics from name-grouped SAM/BAM and write <out>.qc.json.
    Run(RunArgs),
    /// Scan a reference FASTA once and write the .warpqc cache.
    PrepareRef {
        /// Reference FASTA (with .fai).
        #[arg(long)]
        reference: PathBuf,
        /// Cache file to write.
        #[arg(long)]
        out: PathBuf,
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
    /// Output prefix: writes <out>.qc.json.
    #[arg(long)]
    out: PathBuf,
    /// Block size in bytes handed to each worker.
    #[arg(long, default_value_t = 4 << 20, hide = true)]
    block_size: usize,
    /// SAM or BAM path, or - for stdin.
    input: String,
}

fn with_suffix(prefix: &std::path::Path, suffix: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(suffix);
    PathBuf::from(s)
}

fn run_cmd(a: RunArgs) -> Result<()> {
    if a.threads == 0 {
        return Err(input("--threads must be at least 1"));
    }
    let _ = rayon::ThreadPoolBuilder::new().num_threads(a.threads).build_global();
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
    cfg.depth_thresholds = a.depth_thresholds.unwrap_or_else(|| DEFAULT_THRESHOLDS.to_vec());
    cfg.library = a.library.unwrap_or_default();
    cfg.expect_reads = a.expect_reads;
    cfg.block_size = a.block_size;
    cfg.command = std::env::args().collect::<Vec<_>>().join(" ");
    let out = run(&cfg)?;
    let json = with_suffix(&a.out, ".qc.json");
    if let Some(dir) = json.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir).ctx(format!("creating {}", dir.display()))?;
    }
    let mut bytes = serde_json::to_vec_pretty(&out.report).expect("report serialises");
    bytes.push(b'\n');
    write_atomic(&json, &bytes)
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let res = match cli.cmd {
        Cmd::Run(a) => run_cmd(a),
        Cmd::PrepareRef { reference, out } => Reference::from_fasta(&reference).and_then(|r| r.write_cache(&out)),
    };
    match res {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("warpqc: {e}");
            ExitCode::from(e.exit_code() as u8)
        }
    }
}
