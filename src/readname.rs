//! Illumina read-name parsing: `instrument:run:flowcell:lane:tile:x:y`,
//! optionally followed by an eighth UMI field.

use crate::record::parse_u32;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReadName<'a> {
    /// `flowcell:lane`, borrowed from the name.
    pub lane_key: &'a [u8],
    pub tile: u32,
    pub x: u32,
    pub y: u32,
}

pub fn parse_illumina(qname: &[u8]) -> Option<ReadName<'_>> {
    // Some pipelines keep the /1 and /2 mate suffix.
    let name = match qname {
        [rest @ .., b'/', b'1' | b'2'] => rest,
        _ => qname,
    };
    let mut colons = [0usize; 8];
    let mut n = 0;
    for (i, &c) in name.iter().enumerate() {
        if c == b':' {
            if n == colons.len() {
                return None;
            }
            colons[n] = i;
            n += 1;
        }
    }
    if n != 6 && n != 7 {
        return None;
    }
    let field = |k: usize| -> &[u8] {
        let start = if k == 0 { 0 } else { colons[k - 1] + 1 };
        let end = if k < n { colons[k] } else { name.len() };
        &name[start..end]
    };
    if field(2).is_empty() {
        return None;
    }
    parse_u32(field(3))?;
    Some(ReadName {
        lane_key: &name[colons[1] + 1..colons[3]],
        tile: parse_u32(field(4))?,
        x: parse_u32(field(5))?,
        y: parse_u32(field(6))?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seven_fields() {
        let r = parse_illumina(b"LH00123:45:22XYZLT3:2:1101:12345:6789").unwrap();
        assert_eq!(r.lane_key, b"22XYZLT3:2");
        assert_eq!((r.tile, r.x, r.y), (1101, 12345, 6789));
    }

    #[test]
    fn eight_fields_with_umi_and_mate_suffix() {
        let r = parse_illumina(b"LH00123:45:FC1:1:2678:10:20:ACGTACGT").unwrap();
        assert_eq!(r.lane_key, b"FC1:1");
        assert_eq!((r.tile, r.x, r.y), (2678, 10, 20));
        let s = parse_illumina(b"LH00123:45:FC1:1:2678:10:20/2").unwrap();
        assert_eq!(s.y, 20);
    }

    #[test]
    fn rejects_other_shapes() {
        for n in [
            &b"SRR123456.1"[..],
            b"HWUSI-EAS100R:6:73:941:1973",
            b"a:b:FC:x:1:2:3",
            b"a:b:FC:1:1:2:y",
            b"a:b::1:1:2:3",
            b"a:b:c:1:1:2:3:u:v",
        ] {
            assert!(parse_illumina(n).is_none(), "{}", String::from_utf8_lossy(n));
        }
    }
}
