//! Reference FASTA (via `.fai`): contig lengths, N runs and SAM-style M5
//! checksums, the analysis contig set, and the `.warpqc` cache.

use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use md5::{Digest, Md5};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::error::{IoContext, Result, input};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Contig {
    pub name: String,
    pub len: u64,
    /// MD5 of the upper-cased sequence, as in the SAM `@SQ M5` field.
    pub m5: String,
    /// Half-open 0-based runs of N/n.
    pub n_runs: Vec<(u64, u64)>,
}

impl Contig {
    pub fn non_n(&self) -> u64 {
        self.len - self.n_runs.iter().map(|(a, b)| b - a).sum::<u64>()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct Reference {
    pub contigs: Vec<Contig>,
    /// MD5 of the `.fai` file, used to check that a cache belongs to the
    /// FASTA given on the command line.
    pub fai_md5: String,
}

const CACHE_FORMAT: &str = "warpqc-ref";
const CACHE_VERSION: u32 = 1;

#[derive(Serialize, Deserialize)]
struct CacheFile {
    format: String,
    version: u32,
    reference: Reference,
}

pub(crate) fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

struct FaiEntry {
    name: String,
    len: u64,
    offset: u64,
    line_bases: u64,
    line_width: u64,
}

fn fai_path(fasta: &Path) -> PathBuf {
    let mut p = fasta.as_os_str().to_owned();
    p.push(".fai");
    PathBuf::from(p)
}

fn read_fai(fasta: &Path) -> Result<(Vec<FaiEntry>, String)> {
    let path = fai_path(fasta);
    let mut bytes = Vec::new();
    File::open(&path)
        .and_then(|mut f| f.read_to_end(&mut bytes))
        .map_err(|e| {
            input(format!(
                "cannot read FASTA index {} ({e}); run `samtools faidx` on the reference",
                path.display()
            ))
        })?;
    let mut entries = Vec::new();
    for (i, line) in bytes.split(|&c| c == b'\n').enumerate() {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&[u8]> = line.split(|&c| c == b'\t').collect();
        let num =
            |k: usize| -> Option<u64> { std::str::from_utf8(f.get(k)?).ok()?.trim().parse().ok() };
        let bad = || input(format!("{}: malformed line {}", path.display(), i + 1));
        entries.push(FaiEntry {
            name: String::from_utf8_lossy(f[0]).into_owned(),
            len: num(1).ok_or_else(bad)?,
            offset: num(2).ok_or_else(bad)?,
            line_bases: num(3).ok_or_else(bad)?,
            line_width: num(4).ok_or_else(bad)?,
        });
    }
    if entries.is_empty() {
        return Err(input(format!("{} lists no contigs", path.display())));
    }
    Ok((entries, hex(&Md5::digest(&bytes))))
}

fn scan_contig(fasta: &Path, e: &FaiEntry) -> Result<Contig> {
    let mut f = File::open(fasta).ctx(format!("opening {}", fasta.display()))?;
    f.seek(SeekFrom::Start(e.offset))
        .ctx(format!("seeking in {}", fasta.display()))?;
    let mut r = BufReader::with_capacity(1 << 20, f);
    let mut md5 = Md5::new();
    let mut n_runs = Vec::new();
    let mut run_start: Option<u64> = None;
    let mut pos = 0u64;
    let mut buf = Vec::with_capacity(e.line_width as usize + 2);
    while pos < e.len {
        buf.clear();
        let got = r
            .read_until(b'\n', &mut buf)
            .ctx(format!("reading {}", fasta.display()))?;
        if got == 0 {
            return Err(input(format!(
                "{}: contig {} is shorter than its .fai length",
                fasta.display(),
                e.name
            )));
        }
        while matches!(buf.last(), Some(b'\n' | b'\r')) {
            buf.pop();
        }
        if buf.first() == Some(&b'>') || buf.len() as u64 > e.line_bases {
            return Err(input(format!(
                "{}: .fai does not match the FASTA at contig {}; re-run `samtools faidx`",
                fasta.display(),
                e.name
            )));
        }
        for c in buf.iter_mut() {
            c.make_ascii_uppercase();
            if *c == b'N' {
                run_start.get_or_insert(pos);
            } else if let Some(s) = run_start.take() {
                n_runs.push((s, pos));
            }
            pos += 1;
        }
        md5.update(&buf);
    }
    if pos != e.len {
        return Err(input(format!(
            "{}: contig {} is longer than its .fai length",
            fasta.display(),
            e.name
        )));
    }
    if let Some(s) = run_start {
        n_runs.push((s, pos));
    }
    Ok(Contig {
        name: e.name.clone(),
        len: e.len,
        m5: hex(&md5.finalize()),
        n_runs,
    })
}

impl Reference {
    pub fn from_fasta(fasta: &Path) -> Result<Reference> {
        let (entries, fai_md5) = read_fai(fasta)?;
        let contigs = entries
            .par_iter()
            .map(|e| scan_contig(fasta, e))
            .collect::<Result<Vec<_>>>()?;
        Ok(Reference { contigs, fai_md5 })
    }

    /// Loads a cache written by `prepare-ref` and checks it against the
    /// FASTA's `.fai`.
    pub fn from_cache(cache: &Path, fasta: &Path) -> Result<Reference> {
        let f = File::open(cache).ctx(format!("opening {}", cache.display()))?;
        let c: CacheFile = serde_json::from_reader(BufReader::new(f)).map_err(|e| {
            input(format!(
                "{} is not a WarpQC reference cache: {e}",
                cache.display()
            ))
        })?;
        if c.format != CACHE_FORMAT || c.version != CACHE_VERSION {
            return Err(input(format!(
                "{} has cache format {} v{}; expected {CACHE_FORMAT} v{CACHE_VERSION}; re-run prepare-ref",
                cache.display(),
                c.format,
                c.version
            )));
        }
        let (_, fai_md5) = read_fai(fasta)?;
        if fai_md5 != c.reference.fai_md5 {
            return Err(input(format!(
                "{} was built from a different reference than {}; re-run prepare-ref",
                cache.display(),
                fasta.display()
            )));
        }
        Ok(c.reference)
    }

    pub fn load(fasta: &Path, cache: Option<&Path>) -> Result<Reference> {
        match cache {
            Some(c) => Reference::from_cache(c, fasta),
            None => Reference::from_fasta(fasta),
        }
    }

    pub fn write_cache(&self, path: &Path) -> Result<()> {
        let c = CacheFile {
            format: CACHE_FORMAT.into(),
            version: CACHE_VERSION,
            reference: self.clone(),
        };
        let json = serde_json::to_vec(&c).expect("cache serialises");
        crate::report::write_atomic(path, &json)
    }

    /// Digest of the ordered (name, length, M5) list; independent of FASTA
    /// line width and of the `.fai`.
    pub fn digest(&self) -> String {
        let mut h = Md5::new();
        for c in &self.contigs {
            h.update(format!("{}\t{}\t{}\n", c.name, c.len, c.m5).as_bytes());
        }
        hex(&h.finalize())
    }

    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.contigs.iter().position(|c| c.name == name)
    }
}

/// The contigs coverage is computed over, in reference order.
#[derive(Debug, Clone)]
pub struct AnalysisSet {
    /// Reference indices.
    pub contigs: Vec<usize>,
    /// Default contigs that the reference does not have.
    pub missing_defaults: Vec<String>,
}

fn strip_chr(name: &str) -> &str {
    name.strip_prefix("chr").unwrap_or(name)
}

pub fn is_autosome(name: &str) -> bool {
    matches!(strip_chr(name).parse::<u32>(), Ok(1..=22))
}

impl AnalysisSet {
    /// `spec` is a comma-separated list or a file with one contig per line;
    /// `None` selects chr1-22, X, Y, M (or 1-22, X, Y, MT).
    pub fn resolve(reference: &Reference, spec: Option<&str>) -> Result<AnalysisSet> {
        let mut missing_defaults = Vec::new();
        let names: Vec<String> = match spec {
            Some(s) if Path::new(s).is_file() => std::fs::read_to_string(s)
                .ctx(format!("reading {s}"))?
                .lines()
                .map(str::trim)
                .filter(|l| !l.is_empty() && !l.starts_with('#'))
                .map(String::from)
                .collect(),
            Some(s) => s
                .split(',')
                .map(str::trim)
                .filter(|l| !l.is_empty())
                .map(String::from)
                .collect(),
            None => {
                let chr: Vec<String> = (1..=22)
                    .map(|i| format!("chr{i}"))
                    .chain(["chrX", "chrY", "chrM"].map(String::from))
                    .collect();
                let bare: Vec<String> = (1..=22)
                    .map(|i| i.to_string())
                    .chain(["X", "Y", "MT"].map(String::from))
                    .collect();
                let present =
                    |v: &[String]| v.iter().filter(|n| reference.index_of(n).is_some()).count();
                let pick = if present(&chr) >= present(&bare) {
                    chr
                } else {
                    bare
                };
                let (have, miss): (Vec<String>, Vec<String>) = pick
                    .into_iter()
                    .partition(|n| reference.index_of(n).is_some());
                missing_defaults = miss;
                have
            }
        };
        let mut contigs = Vec::new();
        for n in &names {
            let i = reference
                .index_of(n)
                .ok_or_else(|| input(format!("analysis contig {n} is not in the reference")))?;
            if !contigs.contains(&i) {
                contigs.push(i);
            }
        }
        if contigs.is_empty() {
            return Err(input("no analysis contigs found in the reference"));
        }
        contigs.sort_unstable();
        Ok(AnalysisSet {
            contigs,
            missing_defaults,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mini() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/ref/mini.fa")
    }

    #[test]
    fn scans_mini_reference() {
        let r = Reference::from_fasta(&mini()).unwrap();
        let names: Vec<_> = r.contigs.iter().map(|c| c.name.as_str()).collect();
        assert_eq!(names, ["chr1", "chr2", "chrM", "chrUn_decoy"]);
        assert_eq!(r.contigs[0].len, 400);
        assert_eq!(r.contigs[0].n_runs, vec![(200, 220)]);
        assert_eq!(r.contigs[0].non_n(), 380);
        assert!(r.contigs[1].n_runs.is_empty());
        // Expected values from `samtools dict`.
        assert_eq!(r.contigs[0].m5, "36cc0a69c14906c47035ad94be13b787");
        assert_eq!(r.contigs[3].m5, "0bbd37a5c06775d5c2ed7c8f7a60d921");
    }

    #[test]
    fn cache_round_trip_and_mismatch() {
        let r = Reference::from_fasta(&mini()).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let cache = dir.path().join("mini.warpqc");
        r.write_cache(&cache).unwrap();
        assert_eq!(Reference::from_cache(&cache, &mini()).unwrap(), r);

        // A different FASTA (different .fai) must be rejected.
        let other = dir.path().join("other.fa");
        std::fs::write(&other, ">x\nACGT\n").unwrap();
        std::fs::write(dir.path().join("other.fa.fai"), "x\t4\t3\t4\t5\n").unwrap();
        let e = Reference::from_cache(&cache, &other).unwrap_err();
        assert_eq!(e.exit_code(), 2);
    }

    #[test]
    fn analysis_set_resolution() {
        let r = Reference::from_fasta(&mini()).unwrap();
        let d = AnalysisSet::resolve(&r, None).unwrap();
        assert_eq!(d.contigs, vec![0, 1, 2]);
        assert!(d.missing_defaults.contains(&"chrX".to_string()));
        let l = AnalysisSet::resolve(&r, Some("chr2,chr1")).unwrap();
        assert_eq!(l.contigs, vec![0, 1]);
        assert!(AnalysisSet::resolve(&r, Some("chr1,chr9")).is_err());
        assert!(
            is_autosome("chr7")
                && is_autosome("22")
                && !is_autosome("chrX")
                && !is_autosome("chr23")
        );
    }
}
