//! Groups the records of a block into templates.

use crate::error::{Result, input};
use crate::input::sam::lines;
use crate::record::{Aln, TidMap, flag, parse_sam};

/// All records of one read name.
pub struct Template<'a> {
    pub qname: &'a [u8],
    pub alns: &'a [Aln],
    /// Primary first-of-pair record.
    pub r1: Option<usize>,
    /// Primary second-of-pair record.
    pub r2: Option<usize>,
    /// Primary record of an unpaired read.
    pub single: Option<usize>,
}

impl Template<'_> {
    pub fn get(&self, i: Option<usize>) -> Option<&Aln> {
        i.map(|i| &self.alns[i])
    }

    /// The primary records (R1, R2 or the single read) that are present.
    pub fn primaries(&self) -> impl Iterator<Item = &Aln> {
        [self.r1, self.r2, self.single]
            .into_iter()
            .flatten()
            .map(|i| &self.alns[i])
    }
}

/// Reusable record storage for one worker.
#[derive(Default)]
pub struct Grouper {
    pool: Vec<Aln>,
}

impl Grouper {
    fn slot(&mut self, n: usize) -> &mut Aln {
        if self.pool.len() <= n {
            self.pool.resize_with(n + 1, Aln::default);
        }
        &mut self.pool[n]
    }

    /// Parses `block` and calls `f` once per template.
    pub fn for_each_template(
        &mut self,
        block: &[u8],
        tids: &mut TidMap,
        mut f: impl FnMut(&Template) -> Result<()>,
    ) -> Result<()> {
        let mut n = 0usize;
        let mut cur: &[u8] = &[];
        for line in lines(block) {
            let mut tmp = std::mem::take(self.slot(n));
            let parsed = parse_sam(line, tids, &mut tmp);
            self.pool[n] = tmp;
            let qname = parsed?;
            if qname != cur && n > 0 {
                // The record just parsed starts the next template; emit the
                // finished one from slots 0..n, then move it to slot 0.
                self.emit(cur, n, &mut f)?;
                self.pool.swap(0, n);
                n = 0;
            }
            cur = qname;
            n += 1;
        }
        if n > 0 {
            self.emit(cur, n, &mut f)?;
        }
        Ok(())
    }

    fn emit(
        &self,
        qname: &[u8],
        n: usize,
        f: &mut impl FnMut(&Template) -> Result<()>,
    ) -> Result<()> {
        let alns = &self.pool[..n];
        let mut t = Template {
            qname,
            alns,
            r1: None,
            r2: None,
            single: None,
        };
        for (i, a) in alns.iter().enumerate() {
            if !a.is_primary() {
                continue;
            }
            let slot = if !a.has(flag::PAIRED) {
                &mut t.single
            } else if a.has(flag::READ1) && !a.has(flag::READ2) {
                &mut t.r1
            } else if a.has(flag::READ2) && !a.has(flag::READ1) {
                &mut t.r2
            } else {
                return Err(input(format!(
                    "paired record {} has neither or both of the READ1/READ2 flags",
                    String::from_utf8_lossy(qname)
                )));
            };
            if slot.replace(i).is_some() {
                return Err(input(format!(
                    "read name {} has more than one primary record for the same mate",
                    String::from_utf8_lossy(qname)
                )));
            }
        }
        if t.single.is_some() && (t.r1.is_some() || t.r2.is_some()) {
            return Err(input(format!(
                "read name {} mixes paired and unpaired primary records",
                String::from_utf8_lossy(qname)
            )));
        }
        f(&t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn groups_and_identifies_primaries() {
        let block = b"a\t99\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\tIIII\n\
a\t2147\tchr1\t20\t30\t2H2M\t=\t5\t0\tAC\tII\n\
a\t147\tchr1\t5\t60\t4M\t=\t1\t-8\tACGT\tIIII\n\
b\t0\tchr1\t9\t60\t4M\t*\t0\t0\tACGT\tIIII\n\
b\t256\tchr1\t30\t0\t4M\t*\t0\t0\t*\t*\n\
c\t77\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n\
c\t141\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n";
        let mut tids = TidMap::new(&["chr1".into()]);
        let mut g = Grouper::default();
        let mut seen = vec![];
        g.for_each_template(block, &mut tids, |t| {
            seen.push((
                String::from_utf8_lossy(t.qname).into_owned(),
                t.alns.len(),
                t.r1,
                t.r2,
                t.single,
                t.get(t.r1).map(|a| a.pos),
            ));
            Ok(())
        })
        .unwrap();
        assert_eq!(
            seen,
            vec![
                ("a".into(), 3, Some(0), Some(2), None, Some(0)),
                ("b".into(), 2, None, None, Some(0), None),
                ("c".into(), 2, Some(0), Some(1), None, Some(-1)),
            ]
        );
    }

    #[test]
    fn rejects_two_primary_r1() {
        let block = b"a\t99\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\tIIII\n\
a\t99\tchr1\t5\t60\t4M\t=\t1\t-8\tACGT\tIIII\n";
        let mut tids = TidMap::new(&["chr1".into()]);
        let e = Grouper::default()
            .for_each_template(block, &mut tids, |_| Ok(()))
            .unwrap_err();
        assert_eq!(e.exit_code(), 2);
    }
}
