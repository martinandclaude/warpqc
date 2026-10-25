//! Light alignment record shared by the SAM and BAM readers.
//!
//! An `Aln` is owned but reused: the worker keeps a pool of them and parses
//! into it, so the hot path does not allocate once the buffers have grown.

use crate::error::{Result, input};

pub mod flag {
    pub const PAIRED: u16 = 0x1;
    pub const PROPER_PAIR: u16 = 0x2;
    pub const UNMAPPED: u16 = 0x4;
    pub const MATE_UNMAPPED: u16 = 0x8;
    pub const REVERSE: u16 = 0x10;
    pub const MATE_REVERSE: u16 = 0x20;
    pub const READ1: u16 = 0x40;
    pub const READ2: u16 = 0x80;
    pub const SECONDARY: u16 = 0x100;
    pub const QC_FAIL: u16 = 0x200;
    pub const DUPLICATE: u16 = 0x400;
    pub const SUPPLEMENTARY: u16 = 0x800;
}

/// CIGAR operation codes, in BAM order (`MIDNSHP=X`).
pub mod op {
    pub const M: u8 = 0;
    pub const I: u8 = 1;
    pub const D: u8 = 2;
    pub const N: u8 = 3;
    pub const S: u8 = 4;
    pub const H: u8 = 5;
    pub const P: u8 = 6;
    pub const EQ: u8 = 7;
    pub const X: u8 = 8;

    pub fn from_sam(c: u8) -> Option<u8> {
        Some(match c {
            b'M' => M,
            b'I' => I,
            b'D' => D,
            b'N' => N,
            b'S' => S,
            b'H' => H,
            b'P' => P,
            b'=' => EQ,
            b'X' => X,
            _ => return None,
        })
    }

    /// M, = and X: the ops that put a read base on a reference position.
    #[inline]
    pub fn is_match(o: u8) -> bool {
        o == M || o == EQ || o == X
    }

    #[inline]
    pub fn consumes_ref(o: u8) -> bool {
        matches!(o, M | D | N | EQ | X)
    }

    #[inline]
    pub fn consumes_query(o: u8) -> bool {
        matches!(o, M | I | S | EQ | X)
    }

    #[inline]
    pub fn is_clip(o: u8) -> bool {
        o == S || o == H
    }
}

#[derive(Debug, Clone, Default)]
pub struct Aln {
    pub flag: u16,
    /// Index into the header's @SQ list, or -1.
    pub tid: i32,
    /// 0-based leftmost aligned position, or -1.
    pub pos: i64,
    pub mapq: u8,
    /// (op, length) pairs.
    pub cigar: Vec<(u8, u32)>,
    pub mate_tid: i32,
    pub mate_pos: i64,
    pub tlen: i64,
    /// Phred base qualities. Bases that are N in SEQ get quality 0, so a base
    /// quality filter also drops no-calls. Empty when QUAL is `*`.
    pub qual: Vec<u8>,
    /// NM tag, if present.
    pub nm: Option<u32>,
}

impl Aln {
    #[inline]
    pub fn has(&self, f: u16) -> bool {
        self.flag & f != 0
    }

    #[inline]
    pub fn is_mapped(&self) -> bool {
        !self.has(flag::UNMAPPED) && self.tid >= 0
    }

    #[inline]
    pub fn is_primary(&self) -> bool {
        !self.has(flag::SECONDARY | flag::SUPPLEMENTARY)
    }

    #[inline]
    pub fn is_reverse(&self) -> bool {
        self.has(flag::REVERSE)
    }

    /// Reference bases spanned by the alignment.
    pub fn ref_len(&self) -> i64 {
        self.cigar
            .iter()
            .filter(|(o, _)| op::consumes_ref(*o))
            .map(|&(_, l)| l as i64)
            .sum()
    }

    /// 0-based exclusive end of the alignment.
    pub fn end(&self) -> i64 {
        self.pos + self.ref_len()
    }

    /// Soft and hard clips before the first non-clip operation.
    pub fn leading_clip(&self) -> i64 {
        self.cigar
            .iter()
            .take_while(|(o, _)| op::is_clip(*o))
            .map(|&(_, l)| l as i64)
            .sum()
    }

    /// Soft and hard clips after the last non-clip operation.
    pub fn trailing_clip(&self) -> i64 {
        self.cigar
            .iter()
            .rev()
            .take_while(|(o, _)| op::is_clip(*o))
            .map(|&(_, l)| l as i64)
            .sum()
    }

    /// 0-based unclipped 5' position: POS minus the leading clips for a
    /// forward read, the last aligned base plus the trailing clips for a
    /// reverse read. May be negative or beyond the contig end.
    pub fn unclipped_5p(&self) -> i64 {
        if self.is_reverse() {
            self.end() - 1 + self.trailing_clip()
        } else {
            self.pos - self.leading_clip()
        }
    }

    /// Calls `f(ref_pos, query_pos, len)` for every M/=/X block.
    pub fn for_each_match_block(&self, mut f: impl FnMut(i64, usize, u32)) {
        let mut r = self.pos;
        let mut q = 0usize;
        for &(o, l) in &self.cigar {
            if op::is_match(o) {
                f(r, q, l);
            }
            if op::consumes_ref(o) {
                r += l as i64;
            }
            if op::consumes_query(o) {
                q += l as usize;
            }
        }
    }

    /// Base quality at query position `q`, 0 when QUAL is absent.
    #[inline]
    pub fn bq(&self, q: usize) -> u8 {
        self.qual.get(q).copied().unwrap_or(0)
    }
}

/// Maps RNAME to tid, caching the last hit (records arrive in runs on the
/// same contig far more often than not).
#[derive(Clone)]
pub struct TidMap {
    map: std::collections::HashMap<Vec<u8>, i32>,
    last: Option<(Vec<u8>, i32)>,
}

impl TidMap {
    pub fn new(names: &[String]) -> Self {
        let map = names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.as_bytes().to_vec(), i as i32))
            .collect();
        TidMap { map, last: None }
    }

    pub fn get(&mut self, name: &[u8]) -> Option<i32> {
        if name == b"*" {
            return Some(-1);
        }
        if let Some((n, t)) = &self.last
            && n.as_slice() == name
        {
            return Some(*t);
        }
        let t = *self.map.get(name)?;
        self.last = Some((name.to_vec(), t));
        Some(t)
    }
}

#[inline]
fn parse_uint(s: &[u8]) -> Option<u64> {
    if s.is_empty() || s.len() > 19 {
        return None;
    }
    let mut v: u64 = 0;
    for &c in s {
        let d = c.wrapping_sub(b'0');
        if d > 9 {
            return None;
        }
        v = v * 10 + d as u64;
    }
    Some(v)
}

#[inline]
fn parse_int(s: &[u8]) -> Option<i64> {
    match s.first() {
        Some(b'-') => parse_uint(&s[1..])
            .and_then(|v| i64::try_from(v).ok())
            .map(|v| -v),
        Some(b'+') => parse_uint(&s[1..]).and_then(|v| i64::try_from(v).ok()),
        _ => parse_uint(s).and_then(|v| i64::try_from(v).ok()),
    }
}

pub(crate) fn parse_u32(s: &[u8]) -> Option<u32> {
    parse_uint(s).and_then(|v| u32::try_from(v).ok())
}

fn bad(line: &[u8], what: &str) -> crate::error::Error {
    let shown = String::from_utf8_lossy(&line[..line.len().min(120)]);
    input(format!("malformed SAM record ({what}): {shown}"))
}

/// Parses one SAM record line (no trailing newline) into `aln` and returns
/// its QNAME.
pub fn parse_sam<'a>(line: &'a [u8], tids: &mut TidMap, aln: &mut Aln) -> Result<&'a [u8]> {
    let mut fields = line.split(|&c| c == b'\t');
    let mut next = |what: &str| fields.next().ok_or_else(|| bad(line, what));
    let qname = next("QNAME")?;
    let flag = next("FLAG")?;
    let rname = next("RNAME")?;
    let pos = next("POS")?;
    let mapq = next("MAPQ")?;
    let cigar = next("CIGAR")?;
    let rnext = next("RNEXT")?;
    let pnext = next("PNEXT")?;
    let tlen = next("TLEN")?;
    let seq = next("SEQ")?;
    let qual = next("QUAL")?;

    aln.flag = parse_uint(flag)
        .and_then(|v| u16::try_from(v).ok())
        .ok_or_else(|| bad(line, "FLAG"))?;
    aln.tid = tids
        .get(rname)
        .ok_or_else(|| bad(line, "RNAME not in @SQ header"))?;
    aln.pos = parse_int(pos).ok_or_else(|| bad(line, "POS"))? - 1;
    aln.mapq = parse_uint(mapq)
        .map(|v| v.min(255) as u8)
        .ok_or_else(|| bad(line, "MAPQ"))?;
    aln.mate_tid = if rnext == b"=" {
        aln.tid
    } else {
        tids.get(rnext)
            .ok_or_else(|| bad(line, "RNEXT not in @SQ header"))?
    };
    aln.mate_pos = parse_int(pnext).ok_or_else(|| bad(line, "PNEXT"))? - 1;
    aln.tlen = parse_int(tlen).ok_or_else(|| bad(line, "TLEN"))?;

    aln.cigar.clear();
    if cigar != b"*" {
        let mut n: u32 = 0;
        let mut have_digit = false;
        for &c in cigar {
            if c.is_ascii_digit() {
                n = n
                    .checked_mul(10)
                    .and_then(|v| v.checked_add((c - b'0') as u32))
                    .ok_or_else(|| bad(line, "CIGAR"))?;
                have_digit = true;
            } else {
                let o = op::from_sam(c).ok_or_else(|| bad(line, "CIGAR"))?;
                if !have_digit {
                    return Err(bad(line, "CIGAR"));
                }
                aln.cigar.push((o, n));
                n = 0;
                have_digit = false;
            }
        }
        if have_digit {
            return Err(bad(line, "CIGAR"));
        }
    }

    aln.qual.clear();
    if qual != b"*" {
        if seq != b"*" && seq.len() != qual.len() {
            return Err(bad(line, "SEQ and QUAL lengths differ"));
        }
        aln.qual.extend(qual.iter().map(|&q| q.saturating_sub(33)));
        if seq != b"*" {
            for (q, &b) in aln.qual.iter_mut().zip(seq) {
                if b == b'N' || b == b'n' {
                    *q = 0;
                }
            }
        }
    }

    aln.nm = None;
    for tag in fields {
        if tag.len() > 5 && &tag[..5] == b"NM:i:" {
            aln.nm = parse_u32(&tag[5..]);
            break;
        }
    }
    Ok(qname)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tids() -> TidMap {
        TidMap::new(&["chr1".to_string(), "chr2".to_string()])
    }

    fn parse(line: &str) -> Aln {
        let mut a = Aln::default();
        parse_sam(line.as_bytes(), &mut tids(), &mut a).unwrap();
        a
    }

    #[test]
    fn parses_fields() {
        let a = parse(
            "r1\t99\tchr1\t101\t60\t3S5M2I3M1D4M\t=\t201\t150\tNNAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIII\tAS:i:5\tNM:i:4",
        );
        assert_eq!(a.flag, 99);
        assert_eq!(a.tid, 0);
        assert_eq!(a.pos, 100);
        assert_eq!(a.mate_tid, 0);
        assert_eq!(a.mate_pos, 200);
        assert_eq!(a.tlen, 150);
        assert_eq!(a.nm, Some(4));
        assert_eq!(a.ref_len(), 13);
        assert_eq!(a.end(), 113);
        // N bases in SEQ get quality 0.
        assert_eq!(&a.qual[..3], &[0, 0, 40]);
    }

    #[test]
    fn rnext_other_contig_and_star() {
        let a = parse("r\t65\tchr1\t5\t0\t4M\tchr2\t9\t0\tACGT\t*");
        assert_eq!(a.mate_tid, 1);
        assert!(a.qual.is_empty());
        let u = parse("r\t77\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII");
        assert_eq!(u.tid, -1);
        assert_eq!(u.pos, -1);
        assert!(u.cigar.is_empty());
    }

    #[test]
    fn rejects_malformed() {
        let mut a = Aln::default();
        let mut t = tids();
        for l in [
            "r\t99\tchr1\t1\t60\t4Q\t=\t1\t0\tACGT\tIIII",
            "r\t99\tchr1\t1\t60\tM4\t=\t1\t0\tACGT\tIIII",
            "r\t99\tchr9\t1\t60\t4M\t=\t1\t0\tACGT\tIIII",
            "r\t99\tchr1\t1\t60\t4M\t=\t1\t0\tACGT\tIII",
            "r\t99\tchr1\t1",
            // 19 digits: fits u64 but not i64 (review finding: used to wrap).
            "r\t99\tchr1\t1\t60\t4M\t=\t1\t-9999999999999999999\tACGT\tIIII",
        ] {
            assert!(parse_sam(l.as_bytes(), &mut t, &mut a).is_err(), "{l}");
        }
    }

    #[test]
    fn unclipped_5p_forward_and_reverse() {
        // Forward: POS 101 (0-based 100), 2H3S leading clips -> 95.
        let f = parse("r\t99\tchr1\t101\t60\t2H3S10M\t=\t1\t0\tAAAAAAAAAAAAA\t*");
        assert_eq!(f.unclipped_5p(), 95);
        // Reverse: aligned 100..110 (last base 109) + 4S1H trailing -> 114.
        let r = parse("r\t147\tchr1\t101\t60\t10M4S1H\t=\t1\t0\tAAAAAAAAAAAAAA\t*");
        assert_eq!(r.unclipped_5p(), 114);
        // Reverse with a deletion in the middle: ref span 12.
        let d = parse("r\t147\tchr1\t101\t60\t5M2D5M\t=\t1\t0\tAAAAAAAAAA\t*");
        assert_eq!(d.unclipped_5p(), 111);
        // Leading clip past the contig start goes negative.
        let n = parse("r\t99\tchr1\t2\t60\t5S5M\t=\t1\t0\tAAAAAAAAAA\t*");
        assert_eq!(n.unclipped_5p(), -4);
    }

    #[test]
    fn match_blocks_skip_insertions_and_deletions() {
        let a = parse("r\t99\tchr1\t11\t60\t2S3M1I2M2D3=1X\t=\t1\t0\tAAAAAAAAAAAA\t*");
        let mut blocks = vec![];
        a.for_each_match_block(|r, q, l| blocks.push((r, q, l)));
        assert_eq!(
            blocks,
            vec![(10, 2, 3), (13, 6, 2), (17, 8, 3), (20, 11, 1)]
        );
    }
}
