//! Locus region coherence (`region.rs`).
//!
//! GenBank records for one marker often cover different, barely overlapping stretches (e.g. the
//! 5' COI barcode region vs. 3' COI fragments, or 18S/28S partial domains). Mixing them in one
//! locus column puts non-homologous fragments side by side and destroys phylogenetic signal.
//!
//! Each candidate is mapped onto a long anchor sequence by the dominant k-mer diagonal (offset
//! between matching k-mer positions, tolerant to indels via binning), on both strands. The locus
//! "core" is the anchor stretch covered by at least `core_support` of mapped candidates; a
//! candidate must cover `min_core_coverage` of that core to stay in the locus.

use std::collections::HashMap;

/// Mapped span of a candidate on anchor coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AnchorSpan {
    pub start: usize,
    pub end: usize,
    /// Number of k-mer hits on the winning diagonal (evidence strength).
    pub hits: usize,
}

const K: usize = 10;
const OFFSET_BIN: i64 = 24;
const MIN_HITS: usize = 4;

fn encode(b: u8) -> Option<u64> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

fn kmer_codes(seq: &[u8]) -> Vec<(usize, u64)> {
    let mut out = Vec::with_capacity(seq.len());
    let mut code = 0u64;
    let mut valid = 0usize;
    let mask = (1u64 << (2 * K)) - 1;
    for (i, &b) in seq.iter().enumerate() {
        match encode(b) {
            Some(v) => {
                code = ((code << 2) | v) & mask;
                valid += 1;
                if valid >= K {
                    out.push((i + 1 - K, code));
                }
            }
            None => valid = 0,
        }
    }
    out
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' | b'U' => b'A',
            _ => b'N',
        })
        .collect()
}

/// k-mer index of an anchor sequence.
pub struct AnchorIndex {
    pub len: usize,
    positions: HashMap<u64, Vec<usize>>,
}

impl AnchorIndex {
    pub fn new(anchor: &[u8]) -> Self {
        let mut positions: HashMap<u64, Vec<usize>> = HashMap::new();
        for (pos, code) in kmer_codes(anchor) {
            positions.entry(code).or_default().push(pos);
        }
        Self { len: anchor.len(), positions }
    }

    /// Distinct k-mers / total k-mers. Real marker sequences are close to 1; a repetitive anchor
    /// (tandem repeats, periodic synthetic data) cannot place fragments unambiguously.
    pub fn complexity(&self) -> f64 {
        let total: usize = self.positions.values().map(|v| v.len()).sum();
        if total == 0 {
            0.0
        } else {
            self.positions.len() as f64 / total as f64
        }
    }

    fn map_strand(&self, seq: &[u8]) -> Option<AnchorSpan> {
        let mut bins: HashMap<i64, (usize, i64)> = HashMap::new(); // bin -> (hits, offset sum)
        for (qpos, code) in kmer_codes(seq) {
            if let Some(aps) = self.positions.get(&code) {
                if aps.len() > 8 {
                    continue; // low-complexity / repeated k-mer: uninformative
                }
                for &apos in aps {
                    let off = apos as i64 - qpos as i64;
                    let e = bins.entry(off.div_euclid(OFFSET_BIN)).or_insert((0, 0));
                    e.0 += 1;
                    e.1 += off;
                }
            }
        }
        // Merge each bin with its neighbours so a diagonal straddling a bin edge is not split
        let (best_bin, _) = bins
            .keys()
            .map(|&b| {
                let h: usize = (b - 1..=b + 1).filter_map(|x| bins.get(&x)).map(|v| v.0).sum();
                (b, h)
            })
            .max_by_key(|&(b, h)| (h, -b))?;
        let (hits, off_sum): (usize, i64) = (best_bin - 1..=best_bin + 1)
            .filter_map(|x| bins.get(&x))
            .fold((0, 0), |acc, v| (acc.0 + v.0, acc.1 + v.1));
        if hits < MIN_HITS {
            return None;
        }
        let offset = off_sum / hits as i64;
        let start = offset.max(0) as usize;
        let end = ((offset + seq.len() as i64).max(0) as usize).min(self.len);
        if end <= start {
            return None;
        }
        Some(AnchorSpan { start, end, hits })
    }

    /// Maps a candidate (either strand) onto anchor coordinates.
    pub fn map(&self, seq: &[u8]) -> Option<AnchorSpan> {
        let fwd = self.map_strand(seq);
        let rev = self.map_strand(&revcomp(seq));
        match (fwd, rev) {
            (Some(a), Some(b)) => Some(if b.hits > a.hits { b } else { a }),
            (a, b) => a.or(b),
        }
    }
}

/// Anchor stretch covered by at least `core_support` (fraction) of mapped spans.
pub fn core_window(spans: &[AnchorSpan], anchor_len: usize, core_support: f64) -> Option<(usize, usize)> {
    if spans.is_empty() || anchor_len == 0 {
        return None;
    }
    let mut depth = vec![0i64; anchor_len + 1];
    for s in spans {
        depth[s.start] += 1;
        depth[s.end.min(anchor_len)] -= 1;
    }
    let need = ((spans.len() as f64) * core_support).ceil().max(1.0) as i64;
    let (mut cur, mut best, mut run_start) = (0i64, None::<(usize, usize)>, None::<usize>);
    for (i, d) in depth.iter().enumerate().take(anchor_len) {
        cur += d;
        if cur >= need {
            run_start.get_or_insert(i);
        } else if let Some(s) = run_start.take() {
            if best.map(|(a, b)| i - s > b - a).unwrap_or(true) {
                best = Some((s, i));
            }
        }
    }
    if let Some(s) = run_start {
        if best.map(|(a, b)| anchor_len - s > b - a).unwrap_or(true) {
            best = Some((s, anchor_len));
        }
    }
    best
}

/// Fraction of the core window a span covers.
pub fn core_coverage(span: &AnchorSpan, core: (usize, usize)) -> f64 {
    let (cs, ce) = core;
    let ov = span.end.min(ce).saturating_sub(span.start.max(cs));
    ov as f64 / (ce - cs).max(1) as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pseudo(len: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..len)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn maps_fragments_and_finds_core() {
        let anchor = pseudo(1500, 7);
        let idx = AnchorIndex::new(&anchor);
        // 5' fragments (0..650) dominate; one 3' fragment (850..1500) does not overlap them
        let mut spans = Vec::new();
        for (s, e) in [(0, 650), (10, 660), (5, 640), (0, 700), (850, 1500)] {
            let mut frag = anchor[s..e].to_vec();
            frag[100] = b'A'; // a few substitutions
            frag[200] = b'C';
            spans.push(idx.map(&frag).expect("fragment maps"));
        }
        assert!(spans[0].start < 20 && (630..=670).contains(&spans[0].end));
        assert!(spans[4].start >= 830, "3' fragment maps to the 3' end: {:?}", spans[4]);
        let core = core_window(&spans, idx.len, 0.5).unwrap();
        assert!(core.0 <= 20 && core.1 >= 600, "core is the 5' region: {core:?}");
        assert!(core_coverage(&spans[0], core) > 0.9);
        assert!(core_coverage(&spans[4], core) < 0.05);
        // reverse-complemented fragment still maps
        let rc = revcomp(&anchor[0..650]);
        assert!(idx.map(&rc).map(|s| s.start < 20).unwrap_or(false));
        // unrelated sequence does not map
        assert!(idx.map(&pseudo(600, 99)).is_none());
    }
}
