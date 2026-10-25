//! Duplicate signatures and the sharded first-seen tables.
//!
//! A pair key packs library + both ends' (contig, unclipped 5', strand) into
//! 106 bits; it is exact, not a hash. Ends are ordered by the packed value,
//! i.e. by (contig, position, strand), which is Picard's set criterion:
//! library, both ends' contig and coordinate, and the orientation byte (with
//! a same-position RF forced to FR, which the strand tiebreak does here).

use std::hash::BuildHasher;

use hashbrown::HashTable;
use parking_lot::Mutex;

pub const N_SHARDS: usize = 256;

/// One end: 16-bit contig, 32-bit biased unclipped 5' position, strand bit.
pub fn pack_end(tid: u16, pos5: i64, reverse: bool) -> u64 {
    let biased = (pos5.clamp(i32::MIN as i64, i32::MAX as i64) + (1i64 << 31)) as u64;
    ((tid as u64) << 33) | (biased << 1) | reverse as u64
}

#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct PairKey([u32; 4]);

impl PairKey {
    /// Canonical key for two ends, in either order.
    pub fn new(library: u8, a: u64, b: u64) -> PairKey {
        let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
        let v = ((library as u128) << 98) | ((lo as u128) << 49) | hi as u128;
        PairKey([
            v as u32,
            (v >> 32) as u32,
            (v >> 64) as u32,
            (v >> 96) as u32,
        ])
    }
}

pub fn frag_key(library: u8, end: u64) -> u64 {
    ((library as u64) << 49) | end
}

/// Member flag: R1 is the reverse-strand read (Picard's
/// `orientationForOpticalDuplicates` RF class).
pub const M_R1_REVERSE: u8 = 1;
/// Member flag: the read name gave a tile and x/y.
pub const M_LOCATED: u8 = 2;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
#[repr(C)]
pub struct Member {
    pub lane: u8,
    pub flags: u8,
    pub tile: u16,
    pub x: u32,
    pub y: u32,
}

#[derive(Clone, Copy, Debug)]
#[repr(C)]
pub struct Entry {
    pub key: PairKey,
    pub member: Member,
}

pub struct PairShard {
    pub table: HashTable<Entry>,
    /// Second and later members of a group, unsorted.
    pub extra: Vec<Entry>,
}

pub struct DupTable {
    pairs: Box<[Mutex<PairShard>]>,
    frags: Box<[Mutex<HashTable<u64>>]>,
    hasher: hashbrown::DefaultHashBuilder,
}

#[inline]
fn shard_of(h: u64) -> usize {
    // hashbrown uses the low bits for the bucket and the top 7 for the
    // control tag; take the shard from bits that are neither.
    ((h >> 40) as usize) % N_SHARDS
}

impl DupTable {
    /// `expected_pairs` pre-sizes the tables so they do not double (and
    /// transiently hold two copies) while filling.
    pub fn new(expected_pairs: Option<u64>) -> DupTable {
        let per = expected_pairs.map_or(0, |n| (n as usize).div_ceil(N_SHARDS));
        let pairs = (0..N_SHARDS)
            .map(|_| {
                Mutex::new(PairShard {
                    table: HashTable::with_capacity(per),
                    extra: Vec::new(),
                })
            })
            .collect();
        let frags = (0..N_SHARDS)
            .map(|_| Mutex::new(HashTable::with_capacity(per / 50)))
            .collect();
        DupTable {
            pairs,
            frags,
            hasher: Default::default(),
        }
    }

    /// Records a pair; returns true if its signature is new.
    pub fn insert_pair(&self, key: PairKey, member: Member) -> bool {
        let h = self.hasher.hash_one(key);
        let mut s = self.pairs[shard_of(h)].lock();
        if s.table.find(h, |e| e.key == key).is_some() {
            s.extra.push(Entry { key, member });
            false
        } else {
            let hasher = &self.hasher;
            s.table
                .insert_unique(h, Entry { key, member }, |e| hasher.hash_one(e.key));
            true
        }
    }

    /// Records a fragment; returns true if its signature is new.
    pub fn insert_frag(&self, key: u64) -> bool {
        let h = self.hasher.hash_one(key);
        let mut s = self.frags[shard_of(h)].lock();
        if s.find(h, |&k| k == key).is_some() {
            false
        } else {
            let hasher = &self.hasher;
            s.insert_unique(h, key, |&k| hasher.hash_one(k));
            true
        }
    }

    /// Calls `f` with the members of every group of size >= 2, first-seen
    /// member first. Shards are visited in parallel.
    pub fn for_each_group<T: Send>(
        &mut self,
        init: impl Fn() -> T + Sync + Send,
        f: impl Fn(&mut T, &[Member]) + Sync + Send,
    ) -> Vec<T> {
        use rayon::prelude::*;
        let hasher = &self.hasher;
        self.pairs
            .par_iter_mut()
            .map(|shard| {
                let s = shard.get_mut();
                let mut acc = init();
                s.extra.sort_unstable_by_key(|e| e.key);
                let mut members = Vec::new();
                let mut i = 0;
                while i < s.extra.len() {
                    let key = s.extra[i].key;
                    let first = s
                        .table
                        .find(hasher.hash_one(key), |e| e.key == key)
                        .expect("extra member without a first member")
                        .member;
                    members.clear();
                    members.push(first);
                    while i < s.extra.len() && s.extra[i].key == key {
                        members.push(s.extra[i].member);
                        i += 1;
                    }
                    f(&mut acc, &members);
                }
                acc
            })
            .collect()
    }

    pub fn distinct_pairs(&self) -> u64 {
        self.pairs.iter().map(|s| s.lock().table.len() as u64).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn first_seen_then_duplicate() {
        let t = DupTable::new(Some(10));
        let k = PairKey::new(0, pack_end(0, 10, false), pack_end(0, 109, true));
        let m = Member::default();
        assert!(t.insert_pair(k, m));
        assert!(!t.insert_pair(k, m));
        assert!(t.insert_pair(
            PairKey::new(0, pack_end(0, 10, false), pack_end(0, 111, true)),
            m
        ));
        assert!(t.insert_frag(frag_key(0, pack_end(1, 50, false))));
        assert!(!t.insert_frag(frag_key(0, pack_end(1, 50, false))));
        assert!(t.insert_frag(frag_key(0, pack_end(1, 50, true))));
        assert_eq!(t.distinct_pairs(), 2);
    }

    #[test]
    fn groups_list_first_member_first() {
        let mut t = DupTable::new(None);
        let k = PairKey::new(0, pack_end(2, 5, false), pack_end(2, 90, true));
        for x in [7, 3, 9] {
            t.insert_pair(
                k,
                Member {
                    x,
                    ..Default::default()
                },
            );
        }
        let got = t.for_each_group(Vec::new, |acc, m| {
            acc.push(m.iter().map(|m| m.x).collect::<Vec<_>>())
        });
        let groups: Vec<_> = got.into_iter().flatten().collect();
        assert_eq!(groups, vec![vec![7, 3, 9]]);
    }

    #[test]
    fn negative_and_large_positions_stay_ordered() {
        assert!(pack_end(0, -5, false) < pack_end(0, 0, false));
        assert!(pack_end(0, i32::MAX as i64 - 1, true) < pack_end(1, -100, false));
    }

    proptest! {
        /// Swapping R1 and R2 gives the same key, ties included.
        #[test]
        fn mate_order_invariance(
            t1 in 0u16..4, p1 in -200i64..200, r1: bool,
            t2 in 0u16..4, p2 in -200i64..200, r2: bool,
        ) {
            let a = pack_end(t1, p1, r1);
            let b = pack_end(t2, p2, r2);
            prop_assert_eq!(PairKey::new(0, a, b), PairKey::new(0, b, a));
            // Same position, opposite strands: the key does not depend on
            // which mate is forward.
            let f = pack_end(t1, p1, false);
            let r = pack_end(t1, p1, true);
            prop_assert_eq!(PairKey::new(0, f, r), PairKey::new(0, r, f));
        }

        /// Different signatures never share a key.
        #[test]
        fn keys_are_exact(
            a in (0u16..3, -50i64..50, any::<bool>()),
            b in (0u16..3, -50i64..50, any::<bool>()),
            c in (0u16..3, -50i64..50, any::<bool>()),
            d in (0u16..3, -50i64..50, any::<bool>()),
        ) {
            let e = |x: (u16, i64, bool)| pack_end(x.0, x.1, x.2);
            let same = {
                let mut s1 = [e(a), e(b)]; s1.sort();
                let mut s2 = [e(c), e(d)]; s2.sort();
                s1 == s2
            };
            prop_assert_eq!(PairKey::new(0, e(a), e(b)) == PairKey::new(0, e(c), e(d)), same);
        }
    }
}
