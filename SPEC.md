# WarpQC — Rust WGS QC tool: spec & handoff

Oct 25, 2026 · @Martin

## Purpose and scope

WarpQC is one Rust binary that reads a sample's alignments straight from the aligner and writes one QC report per sample, as HTML and JSON. It replaces the GATK, Picard and samtools QC steps in the WGS service pipeline.

- **Context:** we sequence for groups that do their own analysis and deliver FASTQ. Alignment exists only to produce QC, so no BAM is kept.
- **Input:** minibwa SAM output for all lanes of one sample, streamed through a pipe in aligner order.
- **Output:** `<sample>.qc.json`, `<sample>.qc.html` and a MultiQC-compatible file.
- **Target data:** human WGS, NovaSeq X Plus, 2x151, PCR-free and PCR libraries, GRCh38.
- **Out of scope for v1:** variant calling, BAM/CRAM writing, exome or targeted metrics, and long reads.

## Locked decisions

These were agreed on 25 October 2026 and should not be reopened during the build without a reason.

| # | Decision | Consequence for the build |
| --- | --- | --- |
| 1 | Groups receive FASTQ only | No BAM is written or kept; the tool never needs an index |
| 2 | Aligner is minibwa | Input is minibwa SAM; any SAM/BAM in aligner order must also work |
| 3 | One job per sample, aligner piped into WarpQC | All lanes go through one minibwa run; lane comes from the read name |
| 4 | Streaming over name-grouped pairs, no sorting | No `samtools sort`, `fixmate` or `MC`/`ms` tags needed |
| 5 | Run-to-run coverage variation in the 4th decimal is acceptable | Single-pass duplicate exclusion; no deterministic two-pass mode in v0.1 |
| 6 | No BAM output option | No `--tee-bam`; problem samples are simply re-run |
| 7 | Coverage denominator = non-N genome | Default contigs chr1–22, X, Y, M; other contigs reported as a fraction of reads |
| 8 | PMAD follows the DRAGEN definition | log-scale normalized bin counts, median absolute deviation of adjacent-bin differences |

## Architecture

WarpQC is a single streaming pass over name-grouped pairs. It never sorts, never indexes and never writes alignments; all state lives in memory, about 24 GB for a 30× sample.

```
Lane FASTQs ──▶ minibwa map ──SAM on stdin──▶ Block reader ──▶ Workers × N
(R1, R2/lane)   (all lanes, one run)          (template-safe)    (one template at a time)
                                                                       │
                ┌──────────────────────────────┬───────────────────────┤
                ▼                              ▼                       ▼
        Duplicate table               Depth arrays            Thread-local stats
   (first-seen verdict per pair)   (filtered + raw, per base)  (insert, MAPQ, bins, lanes)
                │                              │                       │
                └──────────────┬───────────────┴───────────────────────┘
                               ▼
                         Final sweep ──▶ S01.qc.json ──▶ HTML + MultiQC
                 (depth histogram, PMAD, GC)  (source of truth)  (rendered from JSON)
```

The duplicate table is the hinge. Each worker asks it whether a pair's signature is new, and only first-seen pairs add to filtered depth and the insert-size histogram.

- **Shared state:** the duplicate table (sharded, mutex per shard) and two per-base depth arrays (atomic counters). Both are sized from the reference at startup.
- **Thread-local state:** every histogram and counter. These are merged once at end of stream.
- **End of stream:** one parallel sweep over the arrays and the duplicate table produces the summary numbers. The JSON is written first, and the HTML and MultiQC files are rendered from it.
- **Why it's possible:** both mates of a pair arrive together. That gives each pair its full signature, its mate overlap and its insert size without seeing any other record. Only cross-pair facts (duplicates, depth) need shared state.

## Input handling

WarpQC reads SAM (or BAM) from stdin or a path, in aligner order, and fails loudly on anything else.

**Accepted input**

- SAM text from stdin (`-`), the primary mode, or a SAM/BAM file path for tests.
- Records in aligner order: all records for one read name are adjacent. That covers the primary R1 and R2 plus any supplementary and secondary records.
- `@SQ` header lines must match the reference FASTA given with `--reference`, by name and length.
- A header with `@HD SO:coordinate` is rejected with a clear error. So is a read name that reappears after its group was closed: the tool detects this with a recent-names check (the last 65,536 closed names) and exits non-zero. `SO:queryname` and collated input keep templates together and are accepted.
- **Completeness check:** each paired template must contain both primary mates. Coordinate-sorted input without an `SO` tag shows up as nearly every template incomplete, which a recent-names window can miss. If more than 1% of templates are incomplete (checked every 10,000 templates and at end of input), the run fails with exit code 2; below that, the count is reported with a warning.
- **Truncation:** a final line without a newline, a final template missing a mate (checked even if earlier templates were incomplete), or a header with no records is an input error (exit 2). For BAM input, a stream that ends mid-record is also exit 2. `--expect-reads N` fails the run unless exactly N templates were read; it also pre-sizes the duplicate table.

**Pair grouping**

- Group consecutive records by QNAME into a template: one primary record per mate, plus supplementary and secondary records.
- Secondary records are counted in the alignment stats only. Supplementary records count towards stats and chimera rates, but not towards coverage, duplicates or insert size.
- Unpaired reads and pairs with one mate unmapped are counted and handled as single-end fragments, following Picard's rules.

**Lane and tile from the read name**

- Parse the Illumina name `instrument:run:flowcell:lane:tile:x:y`, optionally with an eighth UMI field and a `/1` `/2` suffix. The lane key is `flowcell:lane`.
- Names that don't parse are pooled into one "unknown" lane and excluded from optical detection, and the JSON carries a warning with the count. If more than half the templates have unparseable names, optical duplicate detection is reported as disabled (optical fields `null`). All other metrics still run.
- The library is a single library per sample by default. `--library` can override the name used in the report.

**Reference and contig set**

- `--reference GRCh38.fa` (with `.fai`) is required. It's used once at startup to build the non-N mask and the per-bin GC content (v0.2). In v0.1 it is read through the `.fai` directly; noodles is used only for BAM input.
- The default analysis contigs are chr1–22, chrX, chrY and chrM. `--contigs` takes a list or file to override this.
- Reads on other contigs (decoys, ALT, HLA, EBV, unplaced) are counted and reported as a fraction of mapped reads. They're excluded from coverage.
- The mask, GC table and contig list can be cached to a small sidecar file (`--ref-cache`) so startup takes seconds. The cache stores contig lengths, N runs and per-contig M5, and is checked against the MD5 of the FASTA's `.fai`; GC is added in v0.2.
- `reference_md5` in the JSON is the MD5 over the ordered `name<TAB>length<TAB>M5` lines, so it doesn't depend on FASTA line width.

## Metric definitions

Every metric is computed in the single streaming pass or in one sweep over the accumulators afterwards. Where a definition departs from Picard, the JSON says so in a `definitions` block.

### Alignment stats

- **Counts:** total templates, total reads, primary mapped %, both mates mapped %, properly paired %, singletons %, secondary and supplementary records, and QC-fail reads.
- **MAPQ:** a 0–60 histogram, plus fractions at MAPQ 0, ≥20 and ≥60.
- **Chimeras:** two numbers. `chimeric_read_fraction` is the samtools flagstat convention: primary reads whose mate is on another contig and whose *own* MAPQ is ≥5, over reads with both mates mapped. `chimeric_pair_fraction` counts templates with mates on different contigs and *both* MAPQ ≥5, over templates with both mates mapped. (The original text called the both-mates rule "the samtools convention"; samtools checks only the read's own MAPQ.)
- **Clipping and errors:** soft-clipped base fraction, mismatch rate ((NM − inserted and deleted bases) / aligned bases, over primary mapped reads that carry NM; minibwa always writes NM), and indel rate as insertion plus deletion *events* per 100 aligned bases (Picard PF\_INDEL\_RATE).
- **Off-target contigs:** the fraction of mapped reads on non-analysis contigs.

### Insert size

- **Pairs used:** both mates mapped, primary, same contig, not a duplicate, TLEN ≠ 0. Each pair is counted once, using R2's |TLEN| (Picard). No MAPQ filter. minibwa's TLEN is 5′-to-5′.
- **Orientation:** split into FR, RF and TANDEM as htsjdk `SamPairUtil.getPairOrientation` in htsjdk ≥ 5.0.0 / Picard ≥ 3.5.0: FR when the forward read's start ≤ the reverse read's end, using both mates' alignments (the MC tag in htsjdk). Picard 3.4.0 (htsjdk 4.2.0) used a strict `<` and R2's start + TLEN, and calls some dovetailed short-insert pairs RF. The concordance run must therefore use Picard ≥ 3.5.0. FR is the headline; the others are reported as fractions.
- **Histogram:** 1 bp resolution up to 2,000 bp, plus an overflow bin.
- **Reported:** median, MAD, mode, 10th/25th/75th/90th percentiles, and mean ± SD, all with htsjdk `Histogram` semantics (median of an even count = mean of the two middle values; percentile = first value whose cumulative fraction reaches p; SD with n − 1). Mean and SD use only inserts ≤ median + 10 × MAD: Picard trims from above only, so a long tail can't inflate the SD. Inserts are binned densely up to 100 kb and sparsely above it, so all statistics are exact and don't depend on the 2,000 bp report cut-off.

### Duplicates

- **Pair signature:** library, plus (contig, unclipped 5′ position, strand) for each mate. The two ends are ordered canonically by the full (contig, position, strand) tuple, so swapping R1 and R2 never changes the key, including when both ends share a position. This matches Picard's set criterion, where a same-position RF pair is forced to FR.
- **Unclipped 5′ position:** for a forward read, POS minus the leading soft and hard clips. For a reverse read, the alignment end plus the trailing clips (htsjdk `getUnclippedStart/End`). It can be negative or beyond the contig end, so it is stored with a 2³¹ bias.
- **Fragments:** reads whose mate is unmapped, and unpaired reads, get their own one-end signature table. In v0.1 they are deduplicated only against other fragments, not against pair ends; this is a documented deviation from Picard (under 1% of reads).
- **Counting:** a group of n members contributes n − 1 duplicates. The duplicate rate = (2 × duplicate pairs + duplicate fragments) / (2 × pairs examined + fragments examined), which is Picard's PERCENT\_DUPLICATION.
- **Optical:** as Picard MarkDuplicates, except that WarpQC uses true single linkage. Picard's keeper-based pass can count slightly fewer in chains (0.17 % fewer on a synthetic set). A group is first split by which end R1 is on (`orientationForOpticalDuplicates`; pairs from opposite strands of a molecule can't be optical copies). Within each part, members on the same flowcell, lane and tile are clustered by single linkage. Two members link when |dx| ≤ d **and** |dy| ≤ d, a per-axis box rather than a Euclidean distance, with d = `--optical-distance` (default 2500, the patterned-flowcell value). A cluster of size k contributes k − 1 optical duplicates. Groups over 300,000 members get no optical calls (Picard's `MAX_OPTICAL_DUPLICATE_SET_SIZE`).
- **Reported:** total, optical and non-optical duplicate rates; the same rates per lane; the cross-lane fraction; and estimated library size (Picard's Lander–Waterman estimate, optical duplicates excluded, reproduced step for step).
- **Optical rate:** 2 × optical duplicate pairs over the PERCENT\_DUPLICATION denominator.
- **Cross-lane fraction:** duplicate pairs (n − 1 per group) in groups whose members span more than one `flowcell:lane`, over all duplicate pairs.
- **Per-lane rates:** pairs only, attributed at end of stream so they don't depend on thread timing. A group of n has n − 1 duplicates but no natural keeper, so each member counts (n − 1)/n duplicates in its own lane: the expected value under a randomly chosen keeper. The lanes sum to the total, and no lane is favoured. A fixed keeper rule (for example the lowest lane name) would make the first lane look clean and the last lane dirty whenever most duplicates are cross-lane. The value is kept as integer member counts per group size, so it is exact and order-independent.
- **Which member counts:** the first-seen member of a group is kept for coverage and insert size; later members are excluded. Counts are deterministic, but which member is kept depends on thread timing. That's accepted (locked decision 5).

### Coverage

- **Filtered depth (headline):** primary, non-duplicate, not QC-fail, MAPQ ≥20. Only M, = and X bases count, and only those with BQ ≥20.
- **Mate overlap:** a position covered by both mates is counted once, if at least one mate's base passes BQ.
- **Raw depth:** primary mapped M, = and X bases, not QC-fail, duplicates included, with no MAPQ or BQ filter. Mate overlap is counted once, so raw depth matches mosdepth's default mode run with `-F 2820` (which excludes supplementary and QC-fail reads and keeps duplicates); mosdepth also counts overlap once and doesn't count deletions. It's on by default; `--no-raw` turns it off and saves 6.2 GB.
- **Supplementary and unpaired reads:** supplementary records never add depth. Reads whose mate is unmapped do count in filtered depth. Both differ from Picard CollectWgsMetrics, which keeps supplementary reads (it filters only secondary) and by default drops reads with an unmapped mate (`COUNT_UNPAIRED=false`).
- **Denominator:** non-N positions of the analysis contigs. Depth saturates at 65,535.
- **Reported:** mean, median, SD and MAD, and % of bases at ≥1, 5, 10, 15, 20, 30, 50 and 100× (the list is configurable with `--depth-thresholds`). Also the fold-80 penalty (mean / 20th-percentile depth, zero-depth positions included as in Picard; `null` when that percentile is 0), a per-contig mean, an autosomal mean, and the full depth histogram up to the highest depth seen.
- **Exclusion breakdown:** Picard-style fractions of bases lost to MAPQ, duplicates, BQ, overlap and the depth cap. These show at a glance *why* filtered coverage is lower than raw. Each base is counted at non-N analysis positions under the first reason that applies, in Picard's order (MAPQ, duplicate, BQ, overlap, cap), and the denominator is counted bases plus all excluded bases. One difference: Picard charges MAPQ and duplicate exclusions as whole-read aligned bases.

### Uniformity and PMAD

- **Bins:** 1,000 bp (`--bin-size`) on autosomes only, so sex-chromosome copy number doesn't inflate the result. A bin is skipped if it contains an N, or if it overlaps the optional `--exclude-bed` (low mappability or blacklist).
- **Bin count:** filtered (non-duplicate, MAPQ ≥20) reads whose alignment start falls in the bin.
- **Normalization:** v\[i\] = log2(count\[i\] / median count). Bins with a zero count are dropped from the log step, and the number dropped is reported.
- **GC-corrected variant:** before the log, divide each count by the median count of bins in its integer GC% stratum. Strata with fewer than 100 bins are merged with their neighbour.
- **Differences:** d\[i\] = v\[i\] − v\[i+1\] over consecutive retained bins *within* a contig. There is never a difference across a contig boundary.
- **PMAD:** Median(|d\[i\] − Median(d)|), with no 1.4826 scaling. Reported raw and GC-corrected. The GC-corrected value is the one to compare with DRAGEN.
- **Also reported:** CV of bin counts, the Oexle evenness score, and a GC-bias curve (normalized coverage per GC stratum) with Picard's AT and GC dropout.

### Later metrics (v0.3)

- **Sex check:** chrX/autosome and chrY/autosome depth ratios, giving an XX / XY / other call.
- **mtDNA copy number:** 2 × chrM mean / autosomal mean.
- **SNP fingerprint:** allele counts at the NimaGen sample-ID loci, collected in-stream from a small position table, and called genotypes written to the JSON.
- **Contamination estimate:** allele-balance based, at a fixed set of common SNP sites. Its design is still open.

## Outputs

Each run writes three files with the `--out` prefix. The JSON is the source of truth; the HTML and MultiQC files are rendered from it.

**`<prefix>.qc.json`** — a versioned schema, with fractions as 0–1 floats, not percentages:

```json
{
  "schema_version": "1.0",
  "tool": {"name": "warpqc", "version": "0.1.0", "command": "..."},
  "sample": "S01",
  "inputs": {"reference": "GRCh38.fa", "reference_md5": "...", "lanes": ["22XYZLT3:1", "22XYZLT3:2"]},
  "parameters": {"min_mapq": 20, "min_bq": 20, "optical_distance": 2500, "bin_size": 1000},
  "definitions": {"duplicates.fragments": "fragment-vs-fragment only"},
  "alignment": {...},
  "insert_size": {"median": 0, "mad": 0, "mean": 0, "sd": 0, "histogram": [...]},
  "duplicates": {"rate": 0, "optical_rate": 0, "cross_lane_fraction": 0, "est_library_size": 0},
  "coverage": {"filtered": {"mean": 0, "pct_ge": {"5": 0, "10": 0}}, "raw": {...}, "per_contig": {...}},
  "uniformity": {"pmad": 0, "pmad_gc": 0, "cv": 0, "evenness": 0, "fold80": 0, "gc_curve": [...]},
  "per_lane": {"22XYZLT3:1": {...}},
  "qc_status": {"overall": "pass", "checks": [{"metric": "coverage.filtered.mean", "value": 0, "status": "pass"}]}
}
```

**`<prefix>.qc.html`** — a single self-contained file with no external requests:

- A summary table at the top with pass/warn/fail badges from the thresholds.
- A plot for each of: insert size distribution, cumulative coverage (% bases ≥ X), depth histogram, per-contig mean depth, GC-bias curve, and the binned-coverage profile along the genome. Plus a per-lane table.
- Plots are inline SVG generated in Rust. Short axis labels and a legend only: no chart titles and no in-plot annotations.
- The full JSON is embedded in a `<script type="application/json">` tag, so the HTML alone is enough to re-extract every number.

**`<prefix>_mqc.json`** — MultiQC custom content: a general-stats row (mean coverage, % ≥ 30×, duplicate rate, median insert, PMAD) and line-graph sections for insert size and cumulative coverage. MultiQC over a run directory then gives the cross-sample view.

**Thresholds** — an optional `--thresholds qc.yaml` sets warn and fail bounds per JSON path:

```yaml
coverage.filtered.mean: {fail_below: 25, warn_below: 30}
coverage.filtered.pct_ge.20: {fail_below: 0.90, warn_below: 0.95}
duplicates.rate: {warn_above: 0.10, fail_above: 0.20}
uniformity.pmad_gc: {warn_above: 0.25, fail_above: 0.40}
```

The numbers above are placeholders. Real bounds should be set from the first 20–30 production samples.

**Exit codes**

- **0:** report written, whatever the QC status. In v0.1 there are no thresholds, so `qc_status.overall` is `"not_evaluated"`.
- **2:** input error: bad header, coordinate-sorted input, or a reference mismatch.
- **3:** I/O or resource error.

Failing QC never makes the job fail. The pipeline reads `qc_status.overall` instead.

## CLI and Snakemake integration

One Snakemake rule per sample: minibwa reads every lane's FASTQs and pipes SAM into WarpQC. Nothing is written but the reports.

**CLI**

```
warpqc run \
  --sample S01 \
  --reference GRCh38_no_alt.fa \
  --ref-cache GRCh38_no_alt.warpqc \
  --threads 16 \
  --thresholds qc.yaml \
  [--contigs chr1,...,chrM] [--exclude-bed lowmap.bed] [--bin-size 1000] \
  [--min-mapq 20] [--min-bq 20] [--optical-distance 2500] [--no-raw] \
  --out results/S01/S01 \
  -                        # SAM on stdin, or a SAM/BAM path

warpqc prepare-ref --reference GRCh38_no_alt.fa --out GRCh38_no_alt.warpqc
warpqc render --json S01.qc.json --out S01.qc.html   # re-render HTML without re-running
```

**Snakemake rule (sketch)**

```python
rule align_qc:
    input:
        r1=lambda w: sorted(lanes_r1(w.sample)),
        r2=lambda w: sorted(lanes_r2(w.sample)),
        idx="ref/GRCh38_no_alt.fa.mbw",
        cache="ref/GRCh38_no_alt.warpqc",
    output:
        json="qc/{sample}/{sample}.qc.json",
        html="qc/{sample}/{sample}.qc.html",
        mqc="qc/{sample}/{sample}_mqc.json",
    threads: 64
    resources: mem_mb=64000, runtime=180
    shell:
        "set -o pipefail; "
        "minibwa map -t {threads} ref/GRCh38_no_alt.fa "
        "  <(cat {input.r1}) <(cat {input.r2}) "
        "| warpqc run --sample {wildcards.sample} --reference ref/GRCh38_no_alt.fa "
        "    --ref-cache {input.cache} --threads 8 --out qc/{wildcards.sample}/{wildcards.sample} -"
```

**Notes for the rule**

- R1 and R2 lane lists must be sorted identically, so pairs stay in step. Concatenated `.fastq.gz` files are valid multi-member gzip, so `cat` is enough.
- minibwa doesn't handle ALT contigs, so the reference must be a no-ALT analysis set ([minibwa README](https://github.com/lh3/minibwa)). The same FASTA is passed to WarpQC.
- Most threads go to minibwa. WarpQC's 8 threads are for SAM parsing and accumulator updates, and it should keep up with minibwa's output rate.
- `set -o pipefail` makes the job fail if minibwa fails, instead of leaving a report built from a truncated stream. WarpQC also checks that the template count matches the FASTQ read count when `--expect-reads` is given.
- The memory budget is about 25 GB for WarpQC (see Architecture) plus minibwa's index. The figures in the sketch are a starting point; record the real values from the first runs.

## Rust implementation notes

The hot path is parsing SAM text and updating three shared structures: two depth arrays and the duplicate table. Everything else is per-thread histograms merged at the end.

**Crates**

| Need | Crate | Why |
| --- | --- | --- |
| SAM text parsing | hand-rolled on `memchr` | Zero-copy field splitting; a general SAM library would allocate per record |
| BAM input (test path), FASTA + `.fai` | `noodles` (bam, bgzf, fasta) | Pure Rust, no htslib build dependency |
| CLI | `clap` (derive) |  |
| Threads and channels | `crossbeam-channel`, `std::thread` | Fixed pipeline of reader → workers; `rayon` only for the final sweep |
| Duplicate table | `hashbrown` + `parking_lot::Mutex`, 256 shards | Shard picked by signature hash; low contention |
| JSON, YAML | `serde`, `serde_json`, `serde_yaml` |  |
| HTML | `askama` templates + hand-written SVG | No JS charting library, nothing fetched at view time |
| Allocator | `mimalloc` | Many short-lived per-batch buffers |

**Module layout**

```
src/
  main.rs            CLI: run, prepare-ref, render
  input/sam.rs       block reader: splits stdin into ~4 MB blocks ending on a template boundary
  input/bam.rs       BAM reader for tests and file input
  input/group.rs     template grouping, order validation
  record.rs          borrowed light record: flag, tid, pos, mapq, cigar, qual, tlen, NM, read-name fields
  reference.rs       non-N bitmask, contig set, per-bin GC; .warpqc cache read/write
  metrics/align.rs   flagstat, MAPQ, clipping, mismatch, chimeras
  metrics/insert.rs  insert-size histogram and stats
  metrics/dups.rs    signature packing, sharded table, fragment table
  metrics/optical.rs single-linkage clustering per group
  metrics/coverage.rs depth arrays, overlap handling, exclusion counters
  metrics/bins.rs    read-start bin counts (v0.2)
  summary/hist.rs    htsjdk-compatible histogram statistics
  summary/*.rs       final sweep: depth histogram, PMAD, GC curve, evenness, library size
  pipeline.rs        reader thread, workers, merge, report assembly
  report/{json,html,svg,multiqc,thresholds}.rs
```

**Concurrency model**

- One reader thread cuts stdin into blocks that never split a template, and sends them to N worker threads.
- Each worker handles whole templates. It updates its own thread-local histograms (alignment, insert size, bins, per lane).
- To decide duplicate status, the worker inserts the pair's signature into the duplicate-table shard. The insert returns "first seen" or "duplicate".
- Only first-seen pairs add to the filtered depth array and the insert-size histogram. The raw array is updated for every pair.
- Depth arrays are `[AtomicU16]`, saturating at 65,535 through a relaxed compare-exchange loop (`fetch_update` with `checked_add`). A check followed by a `fetch_add` can wrap to 0 when threads race. Reads spread across the genome, so contention is negligible.
- At end of stream, the thread-local accumulators are merged, then `rayon` sweeps the arrays per contig for the histogram, bins and GC.

**Duplicate table encoding**

- **Key:** packed into a `u128`. Each end gets 16 bits of contig index (decoy-containing references have thousands of contigs), 32 bits of unclipped 5′ position and 1 strand bit; library takes 8 bits. That is 106 bits, with no hashing collisions possible.
- **Value:** the first member's (lane u8, flags u8, tile u16, x u32, y u32): 12 bytes. The flags hold the R1-orientation bit for optical detection and "has a location". The key is stored as `[u32; 4]`, so a table slot is 28 bytes plus one control byte. No count is stored; counts come from the per-thread counters and the side list.
- **Further members:** the second and later members go into an append-only `(key, member)` list per shard (28 bytes each). At end of stream that list is sorted by key, and each group is assembled from it plus the first member in the table.
- **Optical clustering:** runs at end of stream over those groups only.
- **Shard choice:** bits 40–47 of the hash, which hashbrown uses neither for the bucket index nor for its 7-bit tag.

**Memory budget at 30× (estimate, to be measured)**

| Structure | Size |
| --- | --- |
| Filtered depth array (u16 × 3.1 Gb) | 6.2 GB |
| Raw depth array (`--no-raw` drops it) | 6.2 GB |
| Duplicate table, about 310 M pairs × 29 B per slot ÷ load factor (0.44–0.87) | \~10–20 GB |
| Side list of duplicate members, 28 B each (\~31 M at 10% duplicates) | \~0.9 GB |
| Non-N bitmask + bins + GC | \~0.5 GB |
| **Total** | **\~24–34 GB** |

The duplicate table grows linearly with depth, to roughly 20–40 GB at 60×. hashbrown's power-of-two bucket counts set where in the range a sample lands; `--expect-reads` reserves capacity up front, which avoids the transient copy during growth but not the load-factor spread. These are estimates to replace with measured RSS.

**Performance targets**

- Parse at least 300 MB/s of SAM on 8 worker threads, so WarpQC is never the slower end of the pipe.
- Total wall time should be minibwa's runtime plus under 2 minutes for the final sweep and rendering.
- Peak RSS stays within the budget above.

**Build and deploy:** a static `x86_64-unknown-linux-musl` release binary. Ship it in an Apptainer image together with a pinned minibwa, pinned by digest, the same way Theseus does.

## Testing and validation

Correctness is shown in three layers: synthetic unit fixtures, a small golden dataset, and full-genome concordance against Picard, mosdepth and DRAGEN. The acceptance bounds below are proposals to confirm before the build.

**Unit fixtures (hand-written SAM, a few records each)**

- [x] Overlapping mates: an overlap is counted once; a base failing BQ in one mate is still counted from the other.
- [x] Duplicates where the reverse read carries trailing soft clips, so the alignment ends differ but the unclipped 5′ ends match. They must group together.
- [x] Mate order invariance: the same pair with R1 and R2 swapped gives the same signature (property test with `proptest`).
- [x] Optical chain: A–B and B–C within the distance but A–C not, giving one cluster of 3 and 2 optical duplicates. A second lane carries a copy, which must count as a cross-lane duplicate.
- [x] Fragment duplicates; supplementary and secondary records excluded from coverage and duplicates.
- [x] Depth saturation at 65,535; N-masked positions excluded from the denominator.
- [ ] (v0.2) PMAD: no difference across a contig boundary; zero-count bins dropped; a hand-computed vector matches exactly.
- [x] Input errors: `SO:coordinate` header, a reappearing read name, `@SQ` mismatch against the reference, unparseable read names (warning, run continues).

**Golden dataset**

- HG002 reads from two lanes, restricted to chr20 (about 2 GB of FASTQ), aligned with the pinned minibwa.
- The expected JSON is committed; CI fails on any change outside stated float tolerances.

**Full-genome concordance (one PCR-free GIAB sample, then 10 production samples)**

The reference run takes the same minibwa output, merges it, coordinate-sorts it and runs the classic tools on it. The PMAD comparison uses DRAGEN 4.5 on the same FASTQ.

Settings the reference run needs so that it measures the same things:
- **Read groups:** one per lane (minibwa runs once for all lanes, so add RGs by lane). Picard compares tiles only within a read group, and tile numbers repeat across lanes.
- **Picard version:** ≥ 3.5.0 (htsjdk ≥ 5.0.0), for the pair-orientation rule (see Insert size).
- **MarkDuplicates first:** CollectInsertSizeMetrics and CollectWgsMetrics must run on its output.
- **Restrict to analysis contigs:** CollectWgsMetrics with `INTERVALS` set to the analysis contigs, otherwise its territory includes decoys and unplaced contigs.
- **Match WarpQC's read selection:** `COUNT_UNPAIRED=true`. Supplementary reads still differ (Picard keeps them), and the gap should be measured.
- **mosdepth for raw depth:** `-F 2820 -Q 0`, default mode, not `--fast-mode`.

| WarpQC metric | Compared against | Proposed acceptance |
| --- | --- | --- |
| Duplicate rate | Picard MarkDuplicates PERCENT\_DUPLICATION | within ±0.2 percentage points |
| Optical duplicates | Picard (OPTICAL\_DUPLICATE\_PIXEL\_DISTANCE=2500) | within ±5% relative |
| Estimated library size | Picard ESTIMATED\_LIBRARY\_SIZE | within ±5% relative |
| Filtered mean coverage | Picard CollectWgsMetrics MEAN\_COVERAGE | within ±1% (Picard caps at 250×) |
| % bases ≥ 5…100× | Picard PCT\_xX | within ±0.5 percentage points |
| Raw mean coverage | mosdepth (MAPQ 0, duplicates kept) | within ±0.5% |
| Insert median; mean and SD | Picard CollectInsertSizeMetrics | median exact; mean and SD within ±1% |
| PMAD (GC-corrected) | DRAGEN `cnv_metrics` PMAD, 10 samples | Spearman ρ ≥ 0.9, and the linear fit reported |

**Robustness and performance**

- [ ] Kill minibwa mid-stream: the job must fail, with no report written. If the stream stops cleanly between records, WarpQC sees a normal EOF, so this is guaranteed by `set -o pipefail` plus Snakemake deleting the outputs of failed jobs, not by WarpQC alone. WarpQC's share: write to a temporary file and rename at the end; exit 2 on a partial last line or an incomplete last template; and exit 2 on an `--expect-reads` mismatch. The rule should always pass `--expect-reads`, taken from the bcl-convert demultiplexing stats.
- [ ] A 30× sample on a worker node: record wall time and peak RSS for minibwa alone and for minibwa piped into WarpQC. Use them to set the rule's SLURM defaults.
- [ ] A 60× sample, to check that duplicate-table memory scales as estimated.

## Roadmap and open items

v0.1 is the scope of the first Claude Code handoff. v0.2 and v0.3 start only after the gate before them passes.

1. **v0.1 Core engine** (this handoff): streaming reader and grouping; alignment stats and insert size; duplicates incl. optical; depth arrays; JSON output.
    - **Gate 1:** unit fixtures and golden tests pass; Picard/mosdepth concordance within the bounds above.
2. **v0.2 Report and PMAD:** bins, PMAD raw and GC-corrected; HTML report and MultiQC output; thresholds YAML; per-lane tables.
    - **Gate 2:** PMAD calibrated against DRAGEN; SLURM defaults recorded from real runs.
3. **v0.3 Extras:** sex check, mtDNA copy number, NimaGen SNP fingerprint, contamination estimate.

**Open items**

- [ ] Confirm the proposed acceptance bounds in Testing and validation.
- [ ] Choose the reference build: GRCh38 no-ALT analysis set, with or without the hs38d1 decoys. minibwa needs no-ALT.
- [ ] Choose the `--exclude-bed` for PMAD bins (ENCODE blacklist, a mappability track, or both).
- [ ] Decide whether fragment duplicates should also match pair ends (full Picard parity) in v0.2.
- [ ] Measure the filtered-coverage gap from excluding supplementary records (Picard CollectWgsMetrics keeps them) on the golden and concordance runs; decide whether to keep the exclusion.
- [ ] Golden dataset: HG002 Illumina reads for chr20 and a pinned minibwa are not on the analysis machine yet; the expected JSON can't be committed until they are.
- [ ] Set real QC thresholds from the first 20–30 production samples.
- [ ] Pick a contamination method for v0.3.
- [x] Name settled on 25 October 2026: WarpQC (`warpqc` on the command line).
