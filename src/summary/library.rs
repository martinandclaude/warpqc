//! Picard's Lander-Waterman library size estimate
//! (`DuplicationMetrics.estimateLibrarySize`), reproduced step for step.

fn f(x: f64, c: f64, n: f64) -> f64 {
    c / x - 1.0 + (-n / x).exp()
}

/// `read_pairs` = pairs examined minus optical duplicate pairs;
/// `unique_pairs` = pairs examined minus duplicate pairs. Returns `None`
/// where Picard leaves the field empty.
pub fn estimate_library_size(read_pairs: u64, unique_pairs: u64) -> Option<u64> {
    let dups = read_pairs.checked_sub(unique_pairs)?;
    if read_pairs == 0 || dups == 0 {
        return None;
    }
    let (c, n) = (unique_pairs as f64, read_pairs as f64);
    let (mut m, mut big_m) = (1.0f64, 100.0f64);
    if f(m * c, c, n) < 0.0 {
        return None;
    }
    while f(big_m * c, c, n) > 0.0 {
        big_m *= 10.0;
    }
    for _ in 0..40 {
        let r = (m + big_m) / 2.0;
        let u = f(r * c, c, n);
        if u == 0.0 {
            break;
        } else if u > 0.0 {
            m = r;
        } else {
            big_m = r;
        }
    }
    Some((c * (m + big_m) / 2.0) as u64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_reference_implementation() {
        // Values from a line-by-line Python port of the Picard method.
        assert_eq!(estimate_library_size(1_000_000, 900_000), Some(4_660_793));
        assert_eq!(estimate_library_size(5, 2), Some(2));
        assert_eq!(estimate_library_size(10, 10), None);
        assert_eq!(estimate_library_size(0, 0), None);
    }
}
