//! SAM header parsing and the checks against the reference.

use crate::error::{Result, input};
use crate::reference::Reference;

#[derive(Debug, Clone, Default)]
pub struct Header {
    /// @SQ name and length, in header order (this order defines the tid).
    pub sq: Vec<(String, u64)>,
    pub sort_order: Option<String>,
    pub group_order: Option<String>,
}

fn tag<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    line.split('\t')
        .skip(1)
        .find_map(|f| f.strip_prefix(key)?.strip_prefix(':'))
}

impl Header {
    pub fn parse(text: &str) -> Result<Header> {
        let mut h = Header::default();
        for line in text.lines() {
            if line.starts_with("@HD") {
                h.sort_order = tag(line, "SO").map(String::from);
                h.group_order = tag(line, "GO").map(String::from);
            } else if line.starts_with("@SQ") {
                let name =
                    tag(line, "SN").ok_or_else(|| input(format!("@SQ line without SN: {line}")))?;
                let len = tag(line, "LN")
                    .and_then(|l| l.parse().ok())
                    .ok_or_else(|| input(format!("@SQ line without a valid LN: {line}")))?;
                h.sq.push((name.to_string(), len));
            }
        }
        Ok(h)
    }

    pub fn names(&self) -> Vec<String> {
        self.sq.iter().map(|(n, _)| n.clone()).collect()
    }

    /// Rejects coordinate-sorted input and checks that @SQ matches the
    /// reference by name and length. Returns header tid -> reference index.
    pub fn check(&self, reference: &Reference) -> Result<Vec<usize>> {
        if self.sort_order.as_deref() == Some("coordinate") {
            return Err(input(
                "input is coordinate-sorted (@HD SO:coordinate); WarpQC needs aligner order, \
                 with all records of a read name adjacent. Pipe the aligner output in directly",
            ));
        }
        if self.sq.is_empty() {
            return Err(input("SAM header has no @SQ lines"));
        }
        let mut map = Vec::with_capacity(self.sq.len());
        for (name, len) in &self.sq {
            let i = reference.index_of(name).ok_or_else(|| {
                input(format!("@SQ {name} is not in the reference; was the input aligned to a different FASTA?"))
            })?;
            let rlen = reference.contigs[i].len;
            if rlen != *len {
                return Err(input(format!(
                    "@SQ {name} has length {len} but the reference has {rlen}"
                )));
            }
            map.push(i);
        }
        if let Some(c) = reference
            .contigs
            .iter()
            .find(|c| !self.sq.iter().any(|(n, _)| n == &c.name))
        {
            return Err(input(format!(
                "reference contig {} is missing from the SAM header; was the input aligned to a different FASTA?",
                c.name
            )));
        }
        Ok(map)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn mini() -> Reference {
        Reference::from_fasta(
            &PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/ref/mini.fa"),
        )
        .unwrap()
    }

    const SQ: &str = "@SQ\tSN:chr1\tLN:400\n@SQ\tSN:chr2\tLN:300\n@SQ\tSN:chrM\tLN:100\n@SQ\tSN:chrUn_decoy\tLN:100\n";

    #[test]
    fn accepts_matching_header_in_any_order() {
        let h = Header::parse(&format!(
            "@HD\tVN:1.6\tSO:unsorted\tGO:query\n{SQ}@PG\tID:minibwa\n"
        ))
        .unwrap();
        assert_eq!(h.check(&mini()).unwrap(), vec![0, 1, 2, 3]);
        let swapped = "@SQ\tSN:chr2\tLN:300\n@SQ\tSN:chr1\tLN:400\n@SQ\tSN:chrM\tLN:100\n@SQ\tSN:chrUn_decoy\tLN:100\n";
        assert_eq!(
            Header::parse(swapped).unwrap().check(&mini()).unwrap(),
            vec![1, 0, 2, 3]
        );
        // queryname order keeps templates together, so it is accepted.
        let qn = Header::parse(&format!("@HD\tVN:1.6\tSO:queryname\n{SQ}")).unwrap();
        assert!(qn.check(&mini()).is_ok());
    }

    #[test]
    fn rejects_coordinate_sorted() {
        let h = Header::parse(&format!("@HD\tVN:1.6\tSO:coordinate\n{SQ}")).unwrap();
        let e = h.check(&mini()).unwrap_err();
        assert_eq!(e.exit_code(), 2);
        assert!(e.to_string().contains("coordinate-sorted"));
    }

    #[test]
    fn rejects_sq_mismatches() {
        let r = mini();
        for bad in [
            SQ.replace("LN:400", "LN:401"),
            SQ.replace("chrUn_decoy", "chrUn_other"),
            SQ.replace("@SQ\tSN:chrM\tLN:100\n", ""),
            String::new(),
        ] {
            assert_eq!(
                Header::parse(&bad)
                    .unwrap()
                    .check(&r)
                    .unwrap_err()
                    .exit_code(),
                2
            );
        }
    }
}
