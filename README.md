# WarpQC

WarpQC computes whole-genome sequencing QC in one streaming pass over aligner
output. It reads name-grouped SAM straight from the aligner on stdin, never sorts,
indexes or writes alignments, and writes one report per sample: JSON, a
self-contained HTML page and MultiQC custom content. Its job is to
replace the Picard, GATK and samtools QC steps in a FASTQ-delivery WGS pipeline.
The name is a nod to Picard.

This is **v0.2**. It covers:
- alignment stats, insert size, and duplicates (including optical and cross-lane);
- filtered and raw per-base depth;
- uniformity: PMAD raw and GC-corrected (DRAGEN definition), CV, evenness, the
  GC-bias curve and a 1 Mb coverage profile;
- QC thresholds, the HTML report, MultiQC output and `warpqc render`.

Sex check, mtDNA copy number, SNP fingerprint and contamination come in v0.3.
[SPEC.md](SPEC.md) is the design document and states every metric definition.

## Install

Static Linux x86_64 binaries are attached to each
[GitHub release](https://github.com/martinandclaude/warpqc/releases):

```bash
curl -LO https://github.com/martinandclaude/warpqc/releases/download/v0.1.0/warpqc-v0.1.0-x86_64-unknown-linux-musl.tar.gz
```

Every push to `main` also builds the static binary as a CI artifact.

To build from source:

```bash
cargo build --release
```

## Usage

Prepare the reference cache once per FASTA (it needs a `.fai` next to it):

```bash
warpqc prepare-ref --reference GRCh38_no_alt.fa --out GRCh38_no_alt.warpqc  # --bin-size 1000
```

Pipe the aligner straight in:

```bash
set -o pipefail
minibwa map -t 56 GRCh38_no_alt.fa <(cat S01_L*_R1.fastq.gz) <(cat S01_L*_R2.fastq.gz) \
  | warpqc run --sample S01 --reference GRCh38_no_alt.fa \
      --ref-cache GRCh38_no_alt.warpqc --threads 8 \
      --expect-reads 412345678 --out results/S01/S01 -
```

By default this writes only `results/S01/S01.qc.json`. Add `--html` for
`S01.qc.html` and `--multiqc` for the MultiQC files; either can also be made later
from the JSON with `warpqc render`. A SAM or BAM path can be given instead of
`-`, which is useful for tests.

| Option | Default | Meaning |
| --- | --- | --- |
| `--contigs` | chr1–22, X, Y, M | Analysis contigs: a comma list or a file with one per line |
| `--min-mapq` / `--min-bq` | 20 / 20 | Filters for filtered depth |
| `--optical-distance` | 2500 | Optical duplicate distance in pixels, per axis (Picard) |
| `--depth-thresholds` | 1,5,10,15,20,30,50,100 | Thresholds for `pct_ge` |
| `--no-raw` | off | Skip the raw depth array (saves 2 bytes per analysis base) |
| `--expect-reads` | none | Expected read pairs. A mismatch fails the run. Also pre-sizes the duplicate table |
| `--library` | sample name | Library name in the report |
| `--html` | off | Also write `<prefix>.qc.html` |
| `--multiqc` | off | Also write the three MultiQC custom-content files |
| `--thresholds` | none | YAML of warn/fail bounds per JSON path; sets `qc_status` (see below) |
| `--bin-size` | 1000 | Bin size for PMAD and the GC curve; match `prepare-ref --bin-size` |
| `--exclude-bed` | none | Regions whose bins are left out of PMAD (blacklist, low mappability) |

### Thresholds

```yaml
coverage.filtered.mean: {fail_below: 25, warn_below: 30}
coverage.filtered.pct_ge.20: {fail_below: 0.90, warn_below: 0.95}
duplicates.rate: {warn_above: 0.10, fail_above: 0.20}
uniformity.pmad_gc: {warn_above: 0.25, fail_above: 0.40}
```

How thresholds are applied:
- Keys are dot-separated paths into the JSON, and values are in JSON units (fractions, not percentages).
- Paths are checked before any input is read, so a typo fails at once with exit 2.
- A metric that comes out null counts as a warning.
- Failing QC never changes the exit code; read `qc_status.overall` (pass, warn or fail).

### Re-rendering

```bash
warpqc render --json results/S01/S01.qc.json --thresholds qc.yaml
```

This rebuilds the HTML and MultiQC files from the JSON (both by default;
`--html` or `--multiqc` picks one). With `--thresholds` it
also re-evaluates `qc_status` and rewrites the JSON, so thresholds set after the
first production samples can be applied to them without re-aligning.

## Input requirements

WarpQC needs **aligner order**: all records of one read name adjacent. It fails
with exit code 2 if:
- the header says `@HD SO:coordinate`;
- a read name reappears after its records were closed;
- more than 1 % of templates are missing a primary mate;
- the stream ends mid-record or mid-template;
- `@SQ` doesn't match the reference.

Read names are parsed as Illumina `instrument:run:flowcell:lane:tile:x:y` for
lanes and optical duplicates. Names that don't parse produce a warning, and the
run continues.

## Exit codes

- **0:** report written.
- **2:** input error (bad or coordinate-sorted input, truncation, reference mismatch).
- **3:** I/O or resource error.

The report is written to a temporary file and renamed at the end, so a failed run
leaves no report behind.

## Output

Each run writes `<prefix>.qc.json`; the HTML (`--html`) and MultiQC
(`--multiqc`) files are optional. The JSON is written first and the others are
rendered from it, so `warpqc render` can produce them at any time later.

- **`<prefix>.qc.json`** (schema 1.1, which reads 1.0 reports too). Fractions are 0–1 floats. Sections:
  - `alignment`, `insert_size`, `duplicates`
  - `coverage.filtered`, `coverage.raw` and `coverage.exclusions`
  - `uniformity`, `per_lane`, `qc_status`, `warnings`
  - a `definitions` block stating each metric's exact definition and where it differs from Picard
- **`<prefix>.qc.html`**: one self-contained page with no external requests and no JavaScript. It has:
  - a summary table with pass/warn/fail badges;
  - six plots: insert size, cumulative coverage, depth distribution, per-contig depth, GC bias and the genome profile;
  - a per-lane table;
  - the full JSON embedded in `<script type="application/json" id="warpqc-report">`.
- **`<prefix>.general_mqc.json`, `.insert_mqc.json`, `.coverage_mqc.json`**:
  MultiQC custom content. These are three files because MultiQC reads one section
  per custom-content file. Running `multiqc` over a run directory gives the
  cross-sample view, with WarpQC columns in General Statistics.

## Memory

At 30× human WGS, expect about:
- 6.2 GB for each depth array;
- 10–20 GB for the duplicate table;
- about 0.5 GB for the reference mask.

These are estimates still to be replaced by measured values (see SPEC.md).

## Tests

```bash
cargo test
```

The unit fixtures are hand-written SAM files in `tests/fixtures/sam/`. Each one
lists its hand-derived expected values in `@CO` header lines. Every fixture is run
with 1 worker and with 4 workers on tiny blocks, and the two reports must be
identical. Some tests also cross-check against `samtools` and run the real
MultiQC on the output when they are installed (set `WARPQC_MULTIQC` to point at
MultiQC); CI installs both.

## License

Apache-2.0, see [LICENSE](LICENSE).
