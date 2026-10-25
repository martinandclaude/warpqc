//! Block reader: cuts the record stream into blocks of roughly `block_size`
//! bytes that always end on a template boundary, and checks aligner order on
//! the way.
//!
//! Order checks, all done here because this is the only sequential stage:
//! - a read name that reappears after its group closed (recent-names window);
//! - templates missing a primary mate (coordinate-sorted input without an
//!   `SO:coordinate` header shows up as nearly every template incomplete);
//! - a truncated stream: a final line without a newline, or a final template
//!   missing a mate.

use std::collections::{HashSet, VecDeque};
use std::hash::BuildHasher;
use std::io::BufRead;

use crate::error::{IoContext, Result, input};
use crate::record::{flag, parse_u32};

/// Source of newline-terminated SAM record lines.
pub trait LineSource: Send {
    /// SAM header text (all `@` lines).
    fn header(&self) -> &str;
    /// Appends the next record line, including its `\n`, to `buf`. Returns
    /// false at end of input. A final line without a newline is appended
    /// as-is so that the caller can report the truncation.
    fn next_line(&mut self, buf: &mut Vec<u8>) -> Result<bool>;
}

pub struct SamSource<R> {
    inner: R,
    header: String,
    /// First record line, read while looking for the end of the header.
    pending: Vec<u8>,
}

impl<R: BufRead + Send> SamSource<R> {
    pub fn new(mut inner: R) -> Result<Self> {
        let mut header = String::new();
        let mut line = Vec::new();
        loop {
            line.clear();
            if inner.read_until(b'\n', &mut line).ctx("reading SAM header")? == 0 {
                break;
            }
            if line.first() != Some(&b'@') {
                break;
            }
            header.push_str(&String::from_utf8_lossy(&line));
        }
        Ok(SamSource {
            inner,
            header,
            pending: line,
        })
    }
}

impl<R: BufRead + Send> LineSource for SamSource<R> {
    fn header(&self) -> &str {
        &self.header
    }

    fn next_line(&mut self, buf: &mut Vec<u8>) -> Result<bool> {
        if !self.pending.is_empty() {
            buf.append(&mut self.pending);
            return Ok(true);
        }
        Ok(self.inner.read_until(b'\n', buf).ctx("reading SAM input")? > 0)
    }
}

/// A run of complete templates, as SAM text with one record per line.
pub struct Block {
    pub data: Vec<u8>,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct OrderStats {
    pub templates: u64,
    pub incomplete: u64,
}

/// Window of recently closed read names checked for reappearance.
const RECENT_NAMES: usize = 1 << 16;
/// Incomplete templates are checked every this many templates and at EOF.
const CHECK_EVERY: u64 = 10_000;
/// More incomplete templates than this fraction means the input is not in
/// aligner order.
pub const MAX_INCOMPLETE_FRACTION: f64 = 0.01;

#[derive(Default)]
struct GroupState {
    paired: bool,
    primary_r1: bool,
    primary_r2: bool,
    primary_single: bool,
}

impl GroupState {
    fn complete(&self) -> bool {
        if self.paired {
            self.primary_r1 && self.primary_r2
        } else {
            self.primary_single
        }
    }
}

pub struct Blocker<S> {
    src: S,
    block_size: usize,
    /// Start of the line that was read but belongs to the next block.
    carry: Vec<u8>,
    cur_name: Vec<u8>,
    group: GroupState,
    recent: VecDeque<u64>,
    recent_set: HashSet<u64>,
    hasher: hashbrown::DefaultHashBuilder,
    pub stats: OrderStats,
    done: bool,
}

fn qname_and_flag(line: &[u8]) -> Result<(&[u8], u16)> {
    let shown = || String::from_utf8_lossy(&line[..line.len().min(120)]).into_owned();
    let t1 = memchr::memchr(b'\t', line).ok_or_else(|| input(format!("malformed SAM record: {}", shown())))?;
    let rest = &line[t1 + 1..];
    let t2 = memchr::memchr(b'\t', rest).ok_or_else(|| input(format!("malformed SAM record: {}", shown())))?;
    let f = parse_u32(&rest[..t2])
        .and_then(|v| u16::try_from(v).ok())
        .ok_or_else(|| input(format!("malformed SAM FLAG: {}", shown())))?;
    Ok((&line[..t1], f))
}

impl<S: LineSource> Blocker<S> {
    pub fn new(src: S, block_size: usize) -> Self {
        Blocker {
            src,
            block_size: block_size.max(1),
            carry: Vec::new(),
            cur_name: Vec::new(),
            group: GroupState::default(),
            recent: VecDeque::with_capacity(RECENT_NAMES),
            recent_set: HashSet::with_capacity(RECENT_NAMES),
            hasher: Default::default(),
            stats: OrderStats::default(),
            done: false,
        }
    }

    pub fn header(&self) -> &str {
        self.src.header()
    }

    fn close_group(&mut self) -> Result<()> {
        if self.cur_name.is_empty() {
            return Ok(());
        }
        self.stats.templates += 1;
        if !self.group.complete() {
            self.stats.incomplete += 1;
        }
        if self.stats.templates % CHECK_EVERY == 0 {
            self.check_incomplete()?;
        }
        let h = self.hasher.hash_one(&self.cur_name);
        if self.recent.len() == RECENT_NAMES {
            let old = self.recent.pop_front().unwrap();
            self.recent_set.remove(&old);
        }
        self.recent.push_back(h);
        self.recent_set.insert(h);
        self.group = GroupState::default();
        Ok(())
    }

    fn check_incomplete(&self) -> Result<()> {
        let s = self.stats;
        if s.templates > 0 && s.incomplete as f64 > MAX_INCOMPLETE_FRACTION * s.templates as f64 {
            return Err(input(format!(
                "{} of {} templates are missing a primary mate; the input is not in aligner \
                 (read-name grouped) order, e.g. coordinate-sorted without an SO:coordinate header",
                s.incomplete, s.templates
            )));
        }
        Ok(())
    }

    /// Updates grouping state for a line. Returns true if the line starts a
    /// new template.
    fn observe(&mut self, line: &[u8]) -> Result<bool> {
        if line.last() != Some(&b'\n') {
            return Err(input("input ends in the middle of a SAM record; the stream was truncated"));
        }
        let (name, f) = qname_and_flag(line)?;
        let new = name != self.cur_name.as_slice();
        if new {
            self.close_group()?;
            if self.recent_set.contains(&self.hasher.hash_one(name)) {
                return Err(input(format!(
                    "read name {} reappears after its records were closed; the input is not in \
                     aligner (read-name grouped) order",
                    String::from_utf8_lossy(name)
                )));
            }
            self.cur_name.clear();
            self.cur_name.extend_from_slice(name);
        }
        if f & (flag::SECONDARY | flag::SUPPLEMENTARY) == 0 {
            if f & flag::PAIRED != 0 {
                self.group.paired = true;
                if f & flag::READ1 != 0 {
                    self.group.primary_r1 = true;
                }
                if f & flag::READ2 != 0 {
                    self.group.primary_r2 = true;
                }
            } else {
                self.group.primary_single = true;
            }
        }
        Ok(new)
    }

    fn finish(&mut self) -> Result<()> {
        let last_incomplete = !self.cur_name.is_empty() && !self.group.complete();
        if last_incomplete && self.stats.incomplete == 0 {
            return Err(input(format!(
                "the last template ({}) is missing a primary mate; the stream was probably truncated",
                String::from_utf8_lossy(&self.cur_name)
            )));
        }
        self.close_group()?;
        self.check_incomplete()
    }

    /// The next block, or `None` at end of input.
    pub fn next_block(&mut self) -> Result<Option<Block>> {
        if self.done {
            return Ok(None);
        }
        let mut data = std::mem::take(&mut self.carry);
        data.reserve(self.block_size + 1024);
        loop {
            let start = data.len();
            if !self.src.next_line(&mut data)? {
                self.finish()?;
                self.done = true;
                return Ok(if data.is_empty() { None } else { Some(Block { data }) });
            }
            let new_template = self.observe(&data[start..])?;
            if new_template && start >= self.block_size {
                self.carry = data[start..].to_vec();
                data.truncate(start);
                return Ok(Some(Block { data }));
            }
        }
    }
}

/// Splits a block into lines without the trailing newline.
pub fn lines(block: &[u8]) -> impl Iterator<Item = &[u8]> {
    block
        .split(|&c| c == b'\n')
        .filter(|l| !l.is_empty())
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn blocker(text: &str, block_size: usize) -> Blocker<SamSource<&[u8]>> {
        Blocker::new(SamSource::new(text.as_bytes()).unwrap(), block_size)
    }

    fn collect(text: &str, block_size: usize) -> Result<Vec<Vec<u8>>> {
        let mut b = blocker(text, block_size);
        let mut out = vec![];
        while let Some(blk) = b.next_block()? {
            out.push(blk.data);
        }
        Ok(out)
    }

    fn pair(name: &str) -> String {
        format!(
            "{name}\t99\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\tIIII\n{name}\t147\tchr1\t5\t60\t4M\t=\t1\t-8\tACGT\tIIII\n"
        )
    }

    #[test]
    fn header_is_split_off() {
        let text = format!("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n{}", pair("a"));
        let b = blocker(&text, 1000);
        assert_eq!(b.header(), "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n");
        assert_eq!(collect(&text, 1000).unwrap().concat(), pair("a").into_bytes());
    }

    #[test]
    fn reappearing_name_is_an_error() {
        let text = format!("{}{}{}", pair("a"), pair("b"), pair("a"));
        let e = collect(&text, 1).unwrap_err();
        assert_eq!(e.exit_code(), 2);
        assert!(e.to_string().contains("reappears"), "{e}");
    }

    #[test]
    fn incomplete_templates_are_an_error() {
        // Coordinate order without the header tag: mates are not adjacent.
        let text = "a\t99\tchr1\t1\t60\t4M\t=\t50\t53\tACGT\tIIII\n\
                    b\t99\tchr1\t10\t60\t4M\t=\t60\t54\tACGT\tIIII\n\
                    a\t147\tchr1\t50\t60\t4M\t=\t1\t-53\tACGT\tIIII\n\
                    b\t147\tchr1\t60\t60\t4M\t=\t10\t-54\tACGT\tIIII\n";
        let e = collect(text, 1).unwrap_err();
        assert!(e.to_string().contains("not in aligner") || e.to_string().contains("reappears"), "{e}");
        // Every template is incomplete even when no name reappears.
        let text2 = "a\t99\tchr1\t1\t60\t4M\t=\t50\t53\tACGT\tIIII\n\
                     b\t147\tchr1\t60\t60\t4M\t=\t10\t-54\tACGT\tIIII\n";
        let e2 = collect(text2, 1).unwrap_err();
        assert!(e2.to_string().contains("missing a primary mate"), "{e2}");
    }

    #[test]
    fn truncation_is_an_error() {
        let full = format!("{}{}", pair("a"), pair("b"));
        // Cut in the middle of the final line.
        let cut = &full[..full.len() - 5];
        assert!(collect(cut, 10).unwrap_err().to_string().contains("truncated"));
        // Cut after R1 of the final template.
        let last_line = pair("b").lines().nth(1).unwrap().len() + 1;
        let cut2 = &full[..full.len() - last_line];
        assert!(cut2.ends_with('\n'));
        let e = collect(cut2, 10).unwrap_err().to_string();
        assert!(e.contains("last template"), "{e}");
    }

    #[test]
    fn unpaired_and_secondary_records() {
        let text = "s\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
                    s\t256\tchr1\t9\t0\t4M\t*\t0\t0\t*\t*\n";
        let mut b = blocker(text, 1);
        while b.next_block().unwrap().is_some() {}
        assert_eq!(b.stats.templates, 1);
        assert_eq!(b.stats.incomplete, 0);
    }

    proptest! {
        /// Blocks reassemble to the input and never split a template,
        /// whatever the block size.
        #[test]
        fn blocks_respect_template_boundaries(
            sizes in proptest::collection::vec(1usize..4, 1..40),
            block_size in 1usize..400,
        ) {
            let mut text = String::new();
            for (i, extra) in sizes.iter().enumerate() {
                text.push_str(&pair(&format!("t{i}")));
                for k in 1..*extra {
                    text.push_str(&format!("t{i}\t2145\tchr1\t{k}\t0\t2H2M\t=\t1\t0\tAC\tII\n"));
                }
            }
            let blocks = collect(&text, block_size).unwrap();
            prop_assert_eq!(blocks.concat(), text.as_bytes().to_vec());
            let mut seen = std::collections::HashSet::new();
            for blk in &blocks {
                let names: std::collections::BTreeSet<Vec<u8>> = lines(blk)
                    .map(|l| l.split(|&c| c == b'\t').next().unwrap().to_vec())
                    .collect();
                for n in names {
                    prop_assert!(seen.insert(n), "template split across blocks");
                }
            }
        }
    }
}
