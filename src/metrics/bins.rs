//! Read-start counts in fixed-size bins on the autosomal analysis contigs,
//! the input to PMAD and the GC-bias curve.
//!
//! A read counts once, in the bin holding its leftmost aligned base, if it is
//! primary, not QC-fail, first-seen (not a duplicate) and MAPQ >= min_mapq.
//! Only full bins exist; a bin is excluded if it contains an N or overlaps
//! the `--exclude-bed`.

use std::path::Path;

use crate::error::{IoContext, Result, input};
use crate::metrics::Layout;
use crate::reference::{NO_GC, Reference, is_autosome};

pub struct BinContig {
    pub tid: u32,
    pub name: String,
    /// Index of the contig's first bin in the flat arrays.
    pub first: usize,
    pub n: usize,
}

pub struct BinLayout {
    pub size: u32,
    pub contigs: Vec<BinContig>,
    /// Per tid: index into `contigs`, or `u32::MAX`.
    by_tid: Vec<u32>,
    /// GC% per bin, `NO_GC` where the bin has an N.
    pub gc: Vec<u8>,
    /// Bins overlapping the exclude BED.
    pub bed_excluded: Vec<bool>,
}

impl BinLayout {
    pub fn new(layout: &Layout, reference: &Reference, tid_to_ref: &[usize]) -> BinLayout {
        let mut by_tid = vec![u32::MAX; layout.names.len()];
        let mut contigs = Vec::new();
        let mut gc = Vec::new();
        for &tid in &layout.analysis {
            let name = &layout.names[tid as usize];
            if !is_autosome(name) {
                continue;
            }
            let c = &reference.contigs[tid_to_ref[tid as usize]];
            by_tid[tid as usize] = contigs.len() as u32;
            contigs.push(BinContig {
                tid,
                name: name.clone(),
                first: gc.len(),
                n: c.gc.len(),
            });
            gc.extend_from_slice(&c.gc);
        }
        let n = gc.len();
        BinLayout {
            size: reference.bin_size,
            contigs,
            by_tid,
            gc,
            bed_excluded: vec![false; n],
        }
    }

    pub fn len(&self) -> usize {
        self.gc.len()
    }

    pub fn is_empty(&self) -> bool {
        self.gc.is_empty()
    }

    /// Flat bin index for a 0-based position, if it falls in a full bin of
    /// a binned contig.
    #[inline]
    pub fn index(&self, tid: i32, pos: i64) -> Option<usize> {
        if tid < 0 || pos < 0 {
            return None;
        }
        let k = *self.by_tid.get(tid as usize)?;
        if k == u32::MAX {
            return None;
        }
        let c = &self.contigs[k as usize];
        let b = (pos as u64 / self.size as u64) as usize;
        (b < c.n).then_some(c.first + b)
    }

    pub fn usable(&self, i: usize) -> bool {
        self.gc[i] != NO_GC && !self.bed_excluded[i]
    }

    /// Marks bins overlapping the BED intervals (0-based, half-open). Returns
    /// a warning naming contigs that are not in the reference.
    pub fn exclude_bed(&mut self, path: &Path, layout: &Layout) -> Result<Option<String>> {
        let text = std::fs::read_to_string(path).ctx(format!("reading {}", path.display()))?;
        let mut unknown: Vec<String> = Vec::new();
        for (i, line) in text.lines().enumerate() {
            let line = line.trim_end();
            if line.is_empty()
                || line.starts_with('#')
                || line.starts_with("track")
                || line.starts_with("browser")
            {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            let bad = || {
                input(format!(
                    "{}: malformed BED line {}: {line}",
                    path.display(),
                    i + 1
                ))
            };
            if f.len() < 3 {
                return Err(bad());
            }
            let (s, e): (u64, u64) = (
                f[1].parse().map_err(|_| bad())?,
                f[2].parse().map_err(|_| bad())?,
            );
            if e < s {
                return Err(bad());
            }
            let Some(tid) = layout.names.iter().position(|n| n == f[0]) else {
                if !unknown.iter().any(|u| u == f[0]) {
                    unknown.push(f[0].to_string());
                }
                continue;
            };
            let k = self.by_tid[tid];
            if k == u32::MAX || e == s {
                continue;
            }
            let c = &self.contigs[k as usize];
            let size = self.size as u64;
            let last = ((e - 1) / size).min(c.n as u64);
            for b in s / size..=last {
                if (b as usize) < c.n {
                    self.bed_excluded[c.first + b as usize] = true;
                }
            }
        }
        Ok((!unknown.is_empty()).then(|| {
            format!(
                "{}: contigs not in the reference were ignored: {}",
                path.display(),
                unknown.join(", ")
            )
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::header::Header;
    use crate::reference::AnalysisSet;
    use std::path::PathBuf;

    fn setup(bin: u32) -> (Layout, BinLayout) {
        let r = Reference::from_fasta(
            &PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/ref/mini.fa"),
            bin,
        )
        .unwrap();
        let h = Header::parse(
            "@SQ\tSN:chr1\tLN:400\n@SQ\tSN:chr2\tLN:300\n@SQ\tSN:chrM\tLN:100\n@SQ\tSN:chrUn_decoy\tLN:100\n",
        )
        .unwrap();
        let map = h.check(&r).unwrap();
        let set = AnalysisSet::resolve(&r, Some("chr1,chr2,chrM")).unwrap();
        let layout = Layout::new(&r, &set, &map).unwrap();
        let bins = BinLayout::new(&layout, &r, &map);
        (layout, bins)
    }

    #[test]
    fn autosomes_only_and_full_bins_only() {
        let (_, b) = setup(30);
        // chr1 400 -> 13 full bins, chr2 300 -> 10; chrM is not binned.
        assert_eq!(
            b.contigs
                .iter()
                .map(|c| (c.name.as_str(), c.n))
                .collect::<Vec<_>>(),
            [("chr1", 13), ("chr2", 10)]
        );
        assert_eq!(b.index(0, 0), Some(0));
        assert_eq!(b.index(0, 389), Some(12));
        assert_eq!(b.index(0, 390), None); // partial last bin
        assert_eq!(b.index(1, 0), Some(13));
        assert_eq!(b.index(2, 5), None); // chrM
        // chr1 bins 6 (180-210) and 7 (210-240) touch the N run 200-220.
        assert!(!b.usable(6) && !b.usable(7) && b.usable(5) && b.usable(8));
    }

    #[test]
    fn bed_exclusion() {
        let (layout, mut b) = setup(20);
        let bed = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/sam/bins.bed");
        let warn = b.exclude_bed(&bed, &layout).unwrap().unwrap();
        assert!(warn.contains("chrZ"), "{warn}");
        let excluded: Vec<usize> = (0..b.len()).filter(|&i| b.bed_excluded[i]).collect();
        // chr1 60-80 is bin 3; chr2 0-1 is chr2 bin 0 = flat 20.
        assert_eq!(excluded, vec![3, 20]);
    }
}
