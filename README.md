# WarpQC

WarpQC computes whole-genome sequencing QC in one streaming pass over aligner
output. It reads name-grouped SAM straight from the aligner on stdin, never sorts,
indexes or writes alignments, and writes one JSON report per sample. Its job is to
replace the Picard, GATK and samtools QC steps in a FASTQ-delivery WGS pipeline.
The name is a nod to Picard.

This is **v0.1**: the core engine. It covers alignment stats, insert size,
duplicates (including optical and cross-lane), filtered and raw per-base depth,
and the JSON report. HTML, MultiQC output, PMAD/GC uniformity and QC thresholds
come in v0.2. [SPEC.md](SPEC.md) is the design document and states every metric
definition.

## Build

```bash
cargo build --release
```

The binary is `target/release/warpqc`. The production target is a static
`x86_64-unknown-linux-musl` build, made on Linux.

## Usage

Prepare the reference cache once per FASTA (it needs a `.fai` next to it):

```bash
warpqc prepare-ref --reference GRCh38_no_alt.fa --out GRCh38_no_alt.warpqc
```

Pipe the aligner straight in:

```bash
set -o pipefail
minibwa map -t 56 GRCh38_no_alt.fa <(cat S01_L*_R1.fastq.gz) <(cat S01_L*_R2.fastq.gz) \
  | warpqc run --sample S01 --reference GRCh38_no_alt.fa \
      --ref-cache GRCh38_no_alt.warpqc --threads 8 \
      --expect-reads 412345678 --out results/S01/S01 -
```

This writes `results/S01/S01.qc.json`. A SAM or BAM path can be given instead of
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

`<prefix>.qc.json` (schema 1.0) contains these sections:
- `alignment`
- `insert_size`
- `duplicates`
- `coverage.filtered`, `coverage.raw` and `coverage.exclusions`
- `per_lane`
- `warnings`
- a `definitions` block that states each metric's exact definition and where it differs from Picard.

Fractions are 0–1 floats.

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
identical. Some tests also cross-check against `samtools` when it is installed.
