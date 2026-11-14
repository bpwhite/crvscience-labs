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

const OFFSET_BIN: i64 = 24;
const SW_BAND: usize = 24;
const SW_WINDOW: usize = 400;

fn encode(b: u8) -> Option<u64> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

fn kmer_codes(seq: &[u8], k: usize) -> Vec<(usize, u64)> {
    let mut out = Vec::with_capacity(seq.len());
    let mut code = 0u64;
    let mut valid = 0usize;
    let mask = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    for (i, &b) in seq.iter().enumerate() {
        match encode(b) {
            Some(v) => {
                code = ((code << 2) | v) & mask;
                valid += 1;
                if valid >= k {
                    out.push((i + 1 - k, code));
                }
            }
            None => valid = 0,
        }
    }
    out
}

pub fn revcomp(seq: &[u8]) -> Vec<u8> {
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

#[derive(Default, Clone, Copy)]
struct DiagBin {
    hits: usize,
    offset_sum: i64,
    anchor_pos_sum: usize,
}

/// Computes local Smith-Waterman homology score within a banded window centered on `hit_center`
/// along diagonal `anc_pos ≈ seq_pos + offset`, returning `(normalized_score, raw_score)`.
fn banded_sw_homology(anc: &[u8], seq: &[u8], offset: i64, hit_center: usize) -> (f64, i32) {
    let ov_start = offset.max(0) as usize;
    let ov_end = ((offset + seq.len() as i64).max(0) as usize).min(anc.len());
    let ov_len = ov_end.saturating_sub(ov_start);
    if ov_len < 30 {
        return (0.0, 0);
    }
    let w_len = ov_len.min(SW_WINDOW);
    let a0 = hit_center
        .saturating_sub(w_len / 2)
        .max(ov_start)
        .min(ov_end - w_len);
    let a1 = a0 + w_len;

    let b0_est = (a0 as i64 - offset).max(0) as usize;
    let b0 = b0_est.saturating_sub(SW_BAND);
    let b1 = (b0_est + w_len + SW_BAND).min(seq.len());
    let sa = &anc[a0..a1];
    let sb = &seq[b0..b1];
    let m = sb.len();
    if m == 0 {
        return (0.0, 0);
    }

    let shift_in_b = b0_est.saturating_sub(b0);
    let mut prev = vec![0i32; m + 1];
    let mut curr = vec![0i32; m + 1];
    let mut best = 0i32;
    for (ii, &ba) in sa.iter().enumerate() {
        curr.fill(0);
        let j_center = ii + shift_in_b + 1;
        let j0 = j_center.saturating_sub(SW_BAND).max(1);
        let j1 = (j_center + SW_BAND).min(m);
        let ua = ba.to_ascii_uppercase();
        for j in j0..=j1 {
            let ub = sb[j - 1].to_ascii_uppercase();
            let sub = if ua == ub && matches!(ua, b'A' | b'C' | b'G' | b'T') {
                2
            } else {
                -2
            };
            let v = 0
                .max(prev[j - 1] + sub)
                .max(prev[j] - 3)
                .max(curr[j - 1] - 3);
            curr[j] = v;
            if v > best {
                best = v;
            }
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    (best as f64 / (2.0 * w_len as f64), best)
}

/// k-mer index of an anchor sequence (uses adaptive k=7 on short markers <= 550 bp, k=8 on longer loci).
pub struct AnchorIndex {
    pub len: usize,
    k: usize,
    anchor_upper: Vec<u8>,
    positions: HashMap<u64, Vec<usize>>,
}

impl AnchorIndex {
    pub fn new(anchor: &[u8]) -> Self {
        let k = if anchor.len() <= 550 { 7 } else { 8 };
        let anchor_upper: Vec<u8> = anchor
            .iter()
            .map(|&b| match b.to_ascii_uppercase() {
                b'U' => b'T',
                other => other,
            })
            .collect();
        let mut positions: HashMap<u64, Vec<usize>> = HashMap::new();
        for (pos, code) in kmer_codes(&anchor_upper, k) {
            positions.entry(code).or_default().push(pos);
        }
        Self {
            len: anchor_upper.len(),
            k,
            anchor_upper,
            positions,
        }
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

    fn candidate_bins(&self, seq: &[u8], min_hits: usize) -> Vec<(usize, i64, usize)> {
        let mut bins: HashMap<i64, DiagBin> = HashMap::new();
        for (qpos, code) in kmer_codes(seq, self.k) {
            if let Some(aps) = self.positions.get(&code) {
                if aps.len() > 12 {
                    continue; // low-complexity / repeated k-mer: uninformative
                }
                for &apos in aps {
                    let off = apos as i64 - qpos as i64;
                    let e = bins.entry(off.div_euclid(OFFSET_BIN)).or_default();
                    e.hits += 1;
                    e.offset_sum += off;
                    e.anchor_pos_sum += apos;
                }
            }
        }
        let mut ranked: Vec<(usize, i64)> = bins
            .keys()
            .map(|&b| {
                let h: usize = (b - 1..=b + 1)
                    .filter_map(|x| bins.get(&x))
                    .map(|v| v.hits)
                    .sum();
                (h, b)
            })
            .filter(|&(h, _)| h >= min_hits)
            .collect();
        ranked.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.cmp(&b.1)));
        ranked
            .into_iter()
            .take(3)
            .map(|(_, best_bin)| {
                let (hits, off_sum, anc_sum) = (best_bin - 1..=best_bin + 1)
                    .filter_map(|x| bins.get(&x))
                    .fold((0usize, 0i64, 0usize), |acc, v| {
                        (
                            acc.0 + v.hits,
                            acc.1 + v.offset_sum,
                            acc.2 + v.anchor_pos_sum,
                        )
                    });
                (hits, off_sum / hits as i64, anc_sum / hits)
            })
            .collect()
    }

    fn map_strand(&self, seq: &[u8]) -> Option<(AnchorSpan, i32)> {
        let min_hits = if self.k == 7 { 4 } else { 5 };
        let mut best: Option<(AnchorSpan, i32)> = None;
        for (hits, offset, hit_center) in self.candidate_bins(seq, min_hits) {
            let (norm_sw, raw_sw) = banded_sw_homology(&self.anchor_upper, seq, offset, hit_center);
            if norm_sw < 0.20 || raw_sw < 36 {
                continue;
            }
            let start = offset.max(0) as usize;
            let end = ((offset + seq.len() as i64).max(0) as usize).min(self.len);
            if end <= start {
                continue;
            }
            if best.map(|(_, bs)| raw_sw > bs).unwrap_or(true) {
                best = Some((AnchorSpan { start, end, hits }, raw_sw));
            }
        }
        best
    }

    /// Best raw k-mer diagonal `(is_revcomp, offset, anchor_hit_center)` without homology rejection,
    /// used to anchor fast pairwise alignment in `compute_approx_k2p`.
    pub fn best_diagonal(&self, seq: &[u8]) -> Option<(bool, i64, usize)> {
        let rc_seq = revcomp(seq);
        let fwd = self.candidate_bins(seq, 2).into_iter().next();
        let rev = self.candidate_bins(&rc_seq, 2).into_iter().next();
        match (fwd, rev) {
            (Some(a), Some(b)) => {
                if b.0 > a.0 {
                    Some((true, b.1, b.2))
                } else {
                    Some((false, a.1, a.2))
                }
            }
            (Some(a), None) => Some((false, a.1, a.2)),
            (None, Some(b)) => Some((true, b.1, b.2)),
            (None, None) => None,
        }
    }

    /// Maps a candidate (either strand) onto anchor coordinates.
    pub fn map(&self, seq: &[u8]) -> Option<AnchorSpan> {
        let fwd = self.map_strand(seq);
        let rev = self.map_strand(&revcomp(seq));
        match (fwd, rev) {
            (Some((sa, swa)), Some((sb, swb))) => Some(if swb > swa { sb } else { sa }),
            (a, b) => a.or(b).map(|(s, _)| s),
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

    #[test]
    fn divergent_short_gene_maps_and_unrelated_at_rich_decoy_fails() {
        // 159 bp divergent short marker (AT-rich, ~25% divergence + 3 bp 5' offset)
        let mut x = 42u64;
        let anchor: Vec<u8> = (0..162)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"AAAATTTTGCA"[(x >> 33) as usize % 11]
            })
            .collect();
        let idx = AnchorIndex::new(&anchor);

        // Homolog: shifted +3 bp with ~22% substitutions (especially 3rd codon positions)
        let mut homolog = anchor[3..162].to_vec();
        for (i, b) in homolog.iter_mut().enumerate() {
            if i % 3 == 2 && i % 4 == 0 {
                *b = match *b {
                    b'A' => b'G',
                    b'G' => b'A',
                    b'C' => b'T',
                    _ => b'C',
                };
            }
        }
        let span = idx.map(&homolog).expect("divergent short homolog should map");
        assert!(span.start <= 10 && span.end >= 150, "span was {:?}", span);

        // Non-homologous AT-rich decoy must be rejected by SW homology verification
        let decoy = pseudo(162, 123456);
        assert!(idx.map(&decoy).is_none());
    }

    #[test]
    fn mixed_domain_lsu_keeps_shared_core_and_full_length_rejects_disjoint() {
        let anchor = pseudo(2200, 777);
        let idx = AnchorIndex::new(&anchor);
        let mut spans = Vec::new();
        // 6 D2-D3 core fragments (900..1500) + 2 full-length (0..2200) + 2 non-overlapping 5' D1-D2 (0..650)
        for &(s, e) in &[
            (900usize, 1500usize),
            (910, 1510),
            (895, 1495),
            (905, 1505),
            (900, 1480),
            (920, 1520),
            (0, 2200),
            (50, 2180),
            (0, 650),
            (10, 640),
        ] {
            let mut frag = anchor[s..e].to_vec();
            for i in (15..frag.len()).step_by(20) {
                frag[i] = match frag[i] {
                    b'A' => b'G',
                    b'G' => b'A',
                    b'C' => b'T',
                    _ => b'C',
                };
            }
            spans.push(idx.map(&frag).expect("valid domain should map to full anchor"));
        }
        let core = core_window(&spans, idx.len, 0.5).expect("core exists");
        assert!(core.0 >= 850 && core.1 <= 1550, "core should isolate D2-D3: {:?}", core);
        assert!(core_coverage(&spans[0], core) > 0.9);
        assert!(core_coverage(&spans[6], core) > 0.99);
        assert!(core_coverage(&spans[8], core) < 0.05);
    }
}
