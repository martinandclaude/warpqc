//! Statistics over integer-keyed count histograms (`h[v]` = number of
//! observations of value `v`), defined exactly as htsjdk's `Histogram` so
//! that Picard-comparable numbers come out identical.

pub fn count(h: &[u64]) -> u64 {
    h.iter().sum()
}

pub fn mean(h: &[u64]) -> f64 {
    let n = count(h);
    if n == 0 {
        return 0.0;
    }
    let s: f64 = h
        .iter()
        .enumerate()
        .map(|(v, &c)| v as f64 * c as f64)
        .sum();
    s / n as f64
}

/// Sample standard deviation (n - 1 denominator), as htsjdk.
pub fn sd(h: &[u64]) -> f64 {
    let n = count(h);
    if n < 2 {
        return 0.0;
    }
    let m = mean(h);
    let ss: f64 = h
        .iter()
        .enumerate()
        .map(|(v, &c)| c as f64 * (v as f64 - m).powi(2))
        .sum();
    (ss / (n - 1) as f64).sqrt()
}

/// htsjdk `getMedian` over (value, count) pairs sorted by value: the mean of
/// the n/2-th and (n/2+1)-th values for even n, the ceil(n/2)-th for odd n
/// (1-based ranks).
fn median_sorted(pairs: impl Iterator<Item = (f64, u64)>, n: u64) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let (lo_rank, hi_rank) = if n.is_multiple_of(2) {
        (n / 2, n / 2 + 1)
    } else {
        (n.div_ceil(2), n.div_ceil(2))
    };
    let (mut lo, mut hi) = (None, None);
    let mut total = 0u64;
    for (v, c) in pairs {
        if c == 0 {
            continue;
        }
        total += c;
        if lo.is_none() && total >= lo_rank {
            lo = Some(v);
        }
        if hi.is_none() && total >= hi_rank {
            hi = Some(v);
            break;
        }
    }
    (lo.unwrap() + hi.unwrap()) / 2.0
}

pub fn median(h: &[u64]) -> f64 {
    median_sorted(h.iter().enumerate().map(|(v, &c)| (v as f64, c)), count(h))
}

/// Median of |v - median|, unscaled (htsjdk `getMedianAbsoluteDeviation`).
pub fn mad(h: &[u64]) -> f64 {
    let n = count(h);
    if n == 0 {
        return 0.0;
    }
    let m = median(h);
    let mut dev: Vec<(f64, u64)> = h
        .iter()
        .enumerate()
        .filter(|(_, c)| **c > 0)
        .map(|(v, &c)| ((v as f64 - m).abs(), c))
        .collect();
    dev.sort_by(|a, b| a.0.total_cmp(&b.0));
    median_sorted(dev.into_iter(), n)
}

/// Smallest value with the largest count.
pub fn mode(h: &[u64]) -> u64 {
    let mut best = 0usize;
    for (v, &c) in h.iter().enumerate() {
        if c > h[best] {
            best = v;
        }
    }
    best as u64
}

/// First value whose cumulative fraction reaches `p` (htsjdk
/// `getPercentile`; no interpolation).
pub fn percentile(h: &[u64], p: f64) -> u64 {
    let n = count(h) as f64;
    let mut so_far = 0u64;
    for (v, &c) in h.iter().enumerate() {
        so_far += c;
        if c > 0 && so_far as f64 / n >= p {
            return v as u64;
        }
    }
    0
}

/// Fraction of observations with value >= `t`.
pub fn fraction_ge(h: &[u64], t: usize) -> f64 {
    let n = count(h);
    if n == 0 {
        return 0.0;
    }
    h.iter().skip(t).sum::<u64>() as f64 / n as f64
}

/// Drops values above `width` (htsjdk `trimByWidth`).
pub fn trim_by_width(h: &[u64], width: usize) -> &[u64] {
    &h[..(width + 1).min(h.len())]
}

/// Highest value with a non-zero count, plus one (length of the dense part).
pub fn used_len(h: &[u64]) -> usize {
    h.iter().rposition(|&c| c > 0).map_or(0, |i| i + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn h(values: &[usize]) -> Vec<u64> {
        let mut out = vec![0u64; values.iter().max().unwrap() + 1];
        for &v in values {
            out[v] += 1;
        }
        out
    }

    #[test]
    fn hand_computed_odd() {
        // [1,1,2,2,2,3,10]
        let x = h(&[1, 1, 2, 2, 2, 3, 10]);
        assert_eq!(count(&x), 7);
        assert_eq!(median(&x), 2.0);
        // |x-2| = [0,0,0,1,1,1,8] -> median 1
        assert_eq!(mad(&x), 1.0);
        assert_eq!(mean(&x), 3.0);
        // squared deviations sum 60, / 6
        assert!((sd(&x) - 10f64.sqrt()).abs() < 1e-12);
        assert_eq!(mode(&x), 2);
        assert_eq!(percentile(&x, 0.2), 1);
        assert_eq!(percentile(&x, 0.5), 2);
        assert_eq!(percentile(&x, 0.9), 10);
        assert!((fraction_ge(&x, 3) - 2.0 / 7.0).abs() < 1e-12);
        assert_eq!(trim_by_width(&x, 3), &[0, 2, 3, 1]);
        assert_eq!(used_len(&x), 11);
    }

    #[test]
    fn hand_computed_even() {
        // [1,2,3,4]: median 2.5; |x-2.5| = [0.5,0.5,1.5,1.5] -> MAD 1.0
        let x = h(&[1, 2, 3, 4]);
        assert_eq!(median(&x), 2.5);
        assert_eq!(mad(&x), 1.0);
        // ties for mode go to the smallest value
        assert_eq!(mode(&x), 1);
    }

    #[test]
    fn empty() {
        let x = vec![0u64; 5];
        assert_eq!(
            (median(&x), mad(&x), mean(&x), sd(&x)),
            (0.0, 0.0, 0.0, 0.0)
        );
    }
}
