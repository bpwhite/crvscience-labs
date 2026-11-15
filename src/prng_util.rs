//! Shared zero-dependency deterministic PRNG (`SplitMix64`, Knuth MMIX `DeterministicLcg64`),
//! 5-entity XML/SVG text escaping, and clade-separation permutation Z-score engine.

/// Deterministic SplitMix64 mixing function for index-addressable Monte Carlo seeds.
#[inline]
pub fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

/// Maps a seed through `splitmix64` to a uniform floating-point variate in `[0.0, 1.0)`.
#[inline]
pub fn rand_u01(seed: u64) -> f64 {
    (splitmix64(seed) >> 11) as f64 * (1.0 / ((1u64 << 53) as f64))
}

/// Zero-dependency Knuth MMIX 64-bit linear congruential generator (`LCG64`).
#[derive(Debug, Clone, Copy)]
pub struct DeterministicLcg64 {
    state: u64,
}

impl DeterministicLcg64 {
    /// Creates a new `DeterministicLcg64` seeded with `seed`.
    pub const fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Advances the internal state using MMIX LCG constants and returns the new 64-bit value.
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    /// Samples an unsigned integer in `[0, upper_exclusive)` using the upper 31 bits.
    #[inline]
    pub fn gen_range_usize(&mut self, upper_exclusive: usize) -> usize {
        if upper_exclusive <= 1 {
            return 0;
        }
        ((self.next_u64() >> 33) as usize) % upper_exclusive
    }

    /// In-place deterministic Fisher-Yates shuffle.
    pub fn shuffle<T>(&mut self, slice: &mut [T]) {
        for k in (1..slice.len()).rev() {
            let r = self.gen_range_usize(k + 1);
            slice.swap(k, r);
        }
    }
}

/// Escapes the 5 standard XML/SVG text entities (`&`, `<`, `>`, `"`, `'`).
pub fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

/// Computes the two-sample between-vs-within clade distance separation permutation Z-score
/// (`delta = mean_between - mean_within`) using `DeterministicLcg64` Fisher-Yates shuffles.
pub fn permutation_between_within_z<T: Clone + PartialEq>(
    labels: &[T],
    pair_dists: &[(usize, usize, f64)],
    obs_delta: f64,
    b_perms: usize,
    seed: u64,
) -> f64 {
    let perms = b_perms.max(1);
    let mut perm_deltas = Vec::with_capacity(perms);
    let mut rng = DeterministicLcg64::new(seed);
    for _ in 0..perms {
        let mut perm_sec = labels.to_vec();
        rng.shuffle(&mut perm_sec);
        let mut w_sum = 0.0;
        let mut w_cnt = 0usize;
        let mut b_sum = 0.0;
        let mut b_cnt = 0usize;
        for &(i, j, d) in pair_dists {
            if perm_sec[i] == perm_sec[j] {
                w_sum += d;
                w_cnt += 1;
            } else {
                b_sum += d;
                b_cnt += 1;
            }
        }
        let pd = (b_sum / (b_cnt.max(1) as f64)) - (w_sum / (w_cnt.max(1) as f64));
        perm_deltas.push(pd);
    }
    let perm_mean = perm_deltas.iter().sum::<f64>() / (perms as f64);
    let perm_var = perm_deltas
        .iter()
        .map(|&x| (x - perm_mean).powi(2))
        .sum::<f64>()
        / (perms as f64);
    let perm_sd = perm_var.sqrt().max(1e-6);
    (obs_delta - perm_mean) / perm_sd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splitmix64_and_rand_u01_are_deterministic() {
        let a = splitmix64(42);
        let b = splitmix64(42);
        assert_eq!(a, b);
        let u = rand_u01(42);
        assert!((0.0..1.0).contains(&u));
        assert_eq!(u, rand_u01(42));
    }

    #[test]
    fn deterministic_lcg64_matches_mmix_recurrence() {
        let mut rng = DeterministicLcg64::new(0x5DEE_CE66_Du64);
        let first = rng.next_u64();
        let expected = 0x5DEE_CE66_Du64
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        assert_eq!(first, expected);
    }

    #[test]
    fn xml_escape_handles_five_entities() {
        assert_eq!(
            xml_escape("<tag attr=\"x\" alt='y'>a & b</tag>"),
            "&lt;tag attr=&quot;x&quot; alt=&apos;y&apos;&gt;a &amp; b&lt;/tag&gt;"
        );
    }

    #[test]
    fn permutation_z_detects_strong_clade_separation() {
        let labels = vec![
            "A".to_string(),
            "A".to_string(),
            "B".to_string(),
            "B".to_string(),
        ];
        let pairs = vec![
            (0, 1, 0.02),
            (2, 3, 0.03),
            (0, 2, 0.45),
            (0, 3, 0.48),
            (1, 2, 0.44),
            (1, 3, 0.47),
        ];
        let obs_delta = ((0.45 + 0.48 + 0.44 + 0.47) / 4.0) - ((0.02 + 0.03) / 2.0);
        let z = permutation_between_within_z(&labels, &pairs, obs_delta, 100, 0x5DEE_CE66_Du64);
        assert!(z > 1.5, "expected positive clade separation Z, got {z}");
    }
}
