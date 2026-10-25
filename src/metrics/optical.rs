//! End-of-stream analysis of duplicate groups: optical clustering, per-lane
//! attribution and cross-lane groups.
//!
//! Optical clustering follows Picard MarkDuplicates: members are split by
//! R1 orientation (`orientationForOpticalDuplicates`), compared only within
//! the same lane (Picard: read group) and tile, linked when |dx| <= d and
//! |dy| <= d, and each single-linkage cluster of size k gives k - 1 optical
//! duplicates. Groups larger than Picard's `MAX_OPTICAL_DUPLICATE_SET_SIZE`
//! get no optical calls, as in Picard.

use crate::metrics::dups::{M_LOCATED, M_R1_REVERSE, Member};

pub const MAX_OPTICAL_SET: usize = 300_000;

#[derive(Debug, Clone, Default)]
pub struct DupAnalysis {
    pub groups: u64,
    pub dup_pairs: u64,
    pub optical: u64,
    pub cross_lane_dup_pairs: u64,
    pub lane_dups: Vec<u64>,
    pub lane_optical: Vec<u64>,
}

impl DupAnalysis {
    pub fn new(lanes: usize) -> Self {
        DupAnalysis {
            lane_dups: vec![0; lanes],
            lane_optical: vec![0; lanes],
            ..Default::default()
        }
    }

    pub fn merge(&mut self, o: &DupAnalysis) {
        self.groups += o.groups;
        self.dup_pairs += o.dup_pairs;
        self.optical += o.optical;
        self.cross_lane_dup_pairs += o.cross_lane_dup_pairs;
        for (a, b) in self.lane_dups.iter_mut().zip(&o.lane_dups) {
            *a += b;
        }
        for (a, b) in self.lane_optical.iter_mut().zip(&o.lane_optical) {
            *a += b;
        }
    }
}

fn find(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

/// `rank[lane]` orders lanes by name so that the keeper, and with it the
/// per-lane attribution, does not depend on the order lanes were first seen.
pub fn analyze_group(members: &[Member], rank: &[u32], distance: u32, acc: &mut DupAnalysis) {
    let n = members.len();
    if n < 2 {
        return;
    }
    acc.groups += 1;
    acc.dup_pairs += n as u64 - 1;
    let order = |m: &Member| (rank[m.lane as usize], m.tile, m.x, m.y, m.flags);
    let keeper = (0..n).min_by_key(|&i| order(&members[i])).unwrap();
    for (i, m) in members.iter().enumerate() {
        if i != keeper {
            acc.lane_dups[m.lane as usize] += 1;
        }
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
        analyze_group(&g, &[0, 1], 2500, &mut acc);
        assert_eq!(acc.dup_pairs, 4);
        assert_eq!(acc.optical, 2);
        assert_eq!(acc.lane_optical, vec![2, 0]);
        assert_eq!(acc.cross_lane_dup_pairs, 4);
        assert_eq!(acc.lane_dups, vec![3, 1]);
    }

    #[test]
    fn keeper_follows_lane_name_rank_not_lane_id() {
        // Lane id 1 sorts first by name: its member is the keeper.
        let g = [m(0, 1101, 1, 1, false), m(1, 1101, 9, 9, false)];
        let mut acc = DupAnalysis::new(2);
        analyze_group(&g, &[1, 0], 2500, &mut acc);
        assert_eq!(acc.lane_dups, vec![1, 0]);
    }

    #[test]
    fn unlocated_members_are_not_optical() {
        let mut a = m(0, 1101, 1, 1, false);
        let mut b = a;
        a.flags = 0;
        b.flags = 0;
        let mut acc = DupAnalysis::new(1);
        analyze_group(&[a, b], &[0], 2500, &mut acc);
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
        analyze_group(&g, &[0], 2500, &mut acc);
        assert_eq!(acc.optical, 1);
    }
}
