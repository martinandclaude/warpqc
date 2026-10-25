//! End-of-stream analysis of duplicate groups: optical clustering, per-lane
//! attribution and cross-lane groups.
//!
//! Per-lane duplicates: a group of n members has n - 1 duplicates but no
//! natural keeper, so each member counts (n - 1)/n in its own lane. This is
//! the expected attribution under a uniformly random keeper; it sums to the
//! global count and does not favour any lane. It is kept as integer member
//! counts per group size, so the result is exact and independent of order.
//!
//! Optical clustering follows Picard MarkDuplicates: members are split by
//! R1 orientation (`orientationForOpticalDuplicates`), compared only within
//! the same lane (Picard: read group) and tile, linked when |dx| <= d and
//! |dy| <= d, and each single-linkage cluster of size k gives k - 1 optical
//! duplicates. Groups larger than Picard's `MAX_OPTICAL_DUPLICATE_SET_SIZE`
//! get no optical calls, as in Picard.

use std::collections::BTreeMap;

use crate::metrics::dups::{M_LOCATED, M_R1_REVERSE, Member};

pub const MAX_OPTICAL_SET: usize = 300_000;

#[derive(Debug, Clone, Default)]
pub struct DupAnalysis {
    pub groups: u64,
    pub dup_pairs: u64,
    pub optical: u64,
    pub cross_lane_dup_pairs: u64,
    /// Per lane: group size n -> members of that lane in groups of size n.
    pub lane_members: Vec<BTreeMap<u32, u64>>,
    pub lane_optical: Vec<u64>,
}

impl DupAnalysis {
    pub fn new(lanes: usize) -> Self {
        DupAnalysis {
            lane_members: vec![BTreeMap::new(); lanes],
            lane_optical: vec![0; lanes],
            ..Default::default()
        }
    }

    pub fn merge(&mut self, o: &DupAnalysis) {
        self.groups += o.groups;
        self.dup_pairs += o.dup_pairs;
        self.optical += o.optical;
        self.cross_lane_dup_pairs += o.cross_lane_dup_pairs;
        for (a, b) in self.lane_members.iter_mut().zip(&o.lane_members) {
            for (&n, &c) in b {
                *a.entry(n).or_default() += c;
            }
        }
        for (a, b) in self.lane_optical.iter_mut().zip(&o.lane_optical) {
            *a += b;
        }
    }

    /// Fractional duplicate pairs attributed to `lane`.
    pub fn lane_dups(&self, lane: usize) -> f64 {
        self.lane_members[lane]
            .iter()
            .map(|(&n, &c)| c as f64 * (n - 1) as f64 / n as f64)
            .sum()
    }
}

fn find(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

pub fn analyze_group(members: &[Member], distance: u32, acc: &mut DupAnalysis) {
    let n = members.len();
    if n < 2 {
        return;
    }
    acc.groups += 1;
    acc.dup_pairs += n as u64 - 1;
    for m in members {
        *acc.lane_members[m.lane as usize]
            .entry(n as u32)
            .or_default() += 1;
    }
    if members.iter().any(|m| m.lane != members[0].lane) {
        acc.cross_lane_dup_pairs += n as u64 - 1;
    }
    if n > MAX_OPTICAL_SET {
        return;
    }

    let mut located: Vec<Member> = members
        .iter()
        .copied()
        .filter(|m| m.flags & M_LOCATED != 0)
        .collect();
    located.sort_unstable_by_key(|m| (m.flags & M_R1_REVERSE, m.lane, m.tile, m.x, m.y));
    let bucket = |m: &Member| (m.flags & M_R1_REVERSE, m.lane, m.tile);
    let mut start = 0;
    while start < located.len() {
        let b = bucket(&located[start]);
        let end = start
            + located[start..]
                .iter()
                .take_while(|m| bucket(m) == b)
                .count();
        let ms = &located[start..end];
        if ms.len() > 1 {
            let mut parent: Vec<usize> = (0..ms.len()).collect();
            for i in 0..ms.len() {
                for j in i + 1..ms.len() {
                    if ms[j].x - ms[i].x > distance {
                        break;
                    }
                    if ms[i].y.abs_diff(ms[j].y) <= distance {
                        let (a, c) = (find(&mut parent, i), find(&mut parent, j));
                        parent[a] = c;
                    }
                }
            }
            let mut size = vec![0u64; ms.len()];
            for i in 0..ms.len() {
                let r = find(&mut parent, i);
                size[r] += 1;
            }
            let opt: u64 = size.iter().filter(|&&k| k > 1).map(|k| k - 1).sum();
            acc.optical += opt;
            acc.lane_optical[b.1 as usize] += opt;
        }
        start = end;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(lane: u8, tile: u16, x: u32, y: u32, r1_rev: bool) -> Member {
        Member {
            lane,
            tile,
            x,
            y,
            flags: M_LOCATED | if r1_rev { M_R1_REVERSE } else { 0 },
        }
    }

    #[test]
    fn chain_cluster_box_metric() {
        // Mirrors tests/fixtures/sam/optical_chain.sam group 1.
        let g = [
            m(0, 1101, 1000, 1000, false), // A
            m(0, 1101, 3000, 3000, false), // B: box-close to A, Euclidean 2828
            m(0, 1101, 5000, 5000, false), // C: box-close to B, far from A
            m(0, 1101, 1500, 1000, true),  // G: other orientation class
            m(1, 1101, 1000, 1000, false), // D: other lane
        ];
        let mut acc = DupAnalysis::new(2);
        analyze_group(&g, 2500, &mut acc);
        assert_eq!(acc.dup_pairs, 4);
        assert_eq!(acc.optical, 2);
        assert_eq!(acc.lane_optical, vec![2, 0]);
        assert_eq!(acc.cross_lane_dup_pairs, 4);
        // 4 members in lane 0 and 1 in lane 1, each counting 4/5.
        assert!((acc.lane_dups(0) - 3.2).abs() < 1e-12);
        assert!((acc.lane_dups(1) - 0.8).abs() < 1e-12);
    }

    #[test]
    fn per_lane_attribution_is_symmetric() {
        // Review finding: a fixed keeper rule favoured the lowest-named lane.
        let g = [m(0, 1101, 1, 1, false), m(1, 1101, 9, 9, false)];
        let mut acc = DupAnalysis::new(2);
        analyze_group(&g, 2500, &mut acc);
        assert_eq!((acc.lane_dups(0), acc.lane_dups(1)), (0.5, 0.5));
    }

    #[test]
    fn unlocated_members_are_not_optical() {
        let mut a = m(0, 1101, 1, 1, false);
        let mut b = a;
        a.flags = 0;
        b.flags = 0;
        let mut acc = DupAnalysis::new(1);
        analyze_group(&[a, b], 2500, &mut acc);
        assert_eq!((acc.dup_pairs, acc.optical), (1, 0));
    }

    #[test]
    fn distance_boundary_is_inclusive() {
        let g = [
            m(0, 1, 0, 0, false),
            m(0, 1, 2500, 2500, false),
            m(0, 1, 5001, 2500, false),
        ];
        let mut acc = DupAnalysis::new(1);
        analyze_group(&g, 2500, &mut acc);
        assert_eq!(acc.optical, 1);
    }
}
