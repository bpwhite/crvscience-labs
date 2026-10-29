//! Tier 4: Complexity Theory & Information-Theoretic Genomics (`bio-suite complexity`)
//!
//! High-throughput, multi-threaded (`rayon`) implementations of:
//! - Sliding-window k-mer Shannon block entropy (H_1, H_2, H_3) & positional codon entropy
//! - Jensen-Shannon Divergence (JSD) against background base composition
//! - Algorithmic (Kolmogorov) complexity approximation via Lempel-Ziv 1976 (LZ76) parsing
//! - Chaos Game Representation (CGR) & Purine-Pyrimidine DNA Walks with Hurst Exponent (H)
//!   and Detrended Fluctuation Analysis (DFA alpha)
//! - Non-B DNA secondary structure propensity (G4Hunter G-quadruplex & Z-DNA CpG/YR alternation)

use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceComplexityMetrics {
    pub id: String,
    pub taxon: String,
    pub length_bp: usize,
    pub gc_percent: f64,
    pub gc3_percent: f64,
    pub cpg_obs_exp: f64,
    pub shannon_h1_bits: f64,
    pub shannon_h2_bits: f64,
    pub shannon_h3_bits: f64,
    pub codon_pos1_entropy: f64,
    pub codon_pos2_entropy: f64,
    pub codon_pos3_entropy: f64,
    pub lz76_phrases: usize,
    pub lz76_normalized: f64,
    pub dfa_alpha: f64,
    pub hurst_exponent: f64,
    pub g4hunter_max_score: f64,
    pub g4_motif_count: usize,
    pub zdna_propensity: f64,
    pub composite_complexity_index: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SlidingWindowRecord {
    pub window_id: String,
    pub start: usize,
    pub end: usize,
    pub gc_percent: f64,
    pub cpg_obs_exp: f64,
    pub shannon_h3_bits: f64,
    pub js_divergence_bg: f64,
    pub lz76_normalized: f64,
    pub dfa_alpha: f64,
    pub g4hunter_score: f64,
    pub zdna_score: f64,
    pub structural_complexity_z: f64,
}

/// Analyze a batch of sequences in parallel across all CPU threads using `rayon`.
pub fn analyze_sequences_parallel(
    sequences: &[(String, String, String)], // (id, taxon, seq)
) -> Vec<SequenceComplexityMetrics> {
    sequences
        .par_iter()
        .map(|(id, taxon, seq)| compute_sequence_complexity(id, taxon, seq))
        .collect()
}

pub fn compute_sequence_complexity(
    id: &str,
    taxon: &str,
    raw_seq: &str,
) -> SequenceComplexityMetrics {
    let seq: Vec<u8> = raw_seq
        .bytes()
        .filter(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
        .collect();
    let n = seq.len().max(1);

    // Base & CpG counts
    let mut a_cnt = 0usize;
    let mut c_cnt = 0usize;
    let mut g_cnt = 0usize;
    let mut t_cnt = 0usize;
    let mut cpg_cnt = 0usize;

    for i in 0..seq.len() {
        match seq[i] {
            b'A' => a_cnt += 1,
            b'C' => c_cnt += 1,
            b'G' => g_cnt += 1,
            b'T' => t_cnt += 1,
            _ => {}
        }
        if i + 1 < seq.len() && seq[i] == b'C' && seq[i + 1] == b'G' {
            cpg_cnt += 1;
        }
    }

    let gc_percent = ((c_cnt + g_cnt) as f64 / n as f64) * 100.0;
    let expected_cpg = (c_cnt as f64 * g_cnt as f64) / (n as f64).max(1.0);
    let cpg_obs_exp = if expected_cpg > 0.0 {
        cpg_cnt as f64 / expected_cpg
    } else {
        0.0
    };

    // Codon position GC3 & positional entropies
    let (codon_pos1_entropy, codon_pos2_entropy, codon_pos3_entropy, gc3_percent) =
        compute_codon_positional_metrics(&seq);

    let shannon_h1_bits = kmer_shannon_entropy(&seq, 1);
    let shannon_h2_bits = kmer_shannon_entropy(&seq, 2);
    let shannon_h3_bits = kmer_shannon_entropy(&seq, 3);

    // Lempel-Ziv 1976 (LZ76) algorithmic complexity
    let (lz76_phrases, lz76_normalized) = lempel_ziv_76(&seq);

    // Purine-Pyrimidine DNA walk -> DFA scaling exponent alpha & Hurst exponent H
    let dfa_alpha = compute_dna_walk_dfa(&seq);
    let hurst_exponent = compute_dna_walk_hurst(&seq);

    // Non-B DNA structural motifs (G4Hunter & Z-DNA alternating YR propensity)
    let (g4hunter_max_score, g4_motif_count) = scan_g4hunter(&seq, 25, 1.2);
    let zdna_propensity = compute_zdna_propensity(&seq);

    // Composite Information-Structural Complexity Index (0..100 scale)
    let entropy_efficiency = (shannon_h3_bits / 6.0).clamp(0.0, 1.0);
    let lz_component = lz76_normalized.clamp(0.0, 1.2) / 1.2;
    let dfa_deviation = (dfa_alpha - 0.5).abs().clamp(0.0, 0.5) * 2.0;
    let non_b_component = ((g4hunter_max_score / 2.5) + zdna_propensity * 2.0).clamp(0.0, 1.0);
    let composite_complexity_index = (entropy_efficiency * 35.0
        + lz_component * 35.0
        + dfa_deviation * 15.0
        + non_b_component * 15.0)
        * 100.0
        / 100.0;

    let _ = (a_cnt, t_cnt);

    SequenceComplexityMetrics {
        id: id.to_string(),
        taxon: taxon.to_string(),
        length_bp: seq.len(),
        gc_percent: round4(gc_percent),
        gc3_percent: round4(gc3_percent),
        cpg_obs_exp: round4(cpg_obs_exp),
        shannon_h1_bits: round4(shannon_h1_bits),
        shannon_h2_bits: round4(shannon_h2_bits),
        shannon_h3_bits: round4(shannon_h3_bits),
        codon_pos1_entropy: round4(codon_pos1_entropy),
        codon_pos2_entropy: round4(codon_pos2_entropy),
        codon_pos3_entropy: round4(codon_pos3_entropy),
        lz76_phrases,
        lz76_normalized: round4(lz76_normalized),
        dfa_alpha: round4(dfa_alpha),
        hurst_exponent: round4(hurst_exponent),
        g4hunter_max_score: round4(g4hunter_max_score),
        g4_motif_count,
        zdna_propensity: round4(zdna_propensity),
        composite_complexity_index: round4(composite_complexity_index),
    }
}

/// Multi-threaded sliding-window structural & algorithmic complexity scan across a chromosome region.
pub fn scan_sliding_windows_parallel(
    raw_seq: &str,
    window_size: usize,
    step_size: usize,
) -> Vec<SlidingWindowRecord> {
    let seq: Vec<u8> = raw_seq
        .bytes()
        .map(|b| b.to_ascii_uppercase())
        .filter(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
        .collect();

    if seq.len() < window_size {
        return Vec::new();
    }

    let bg_probs = base_distribution(&seq);
    let starts: Vec<usize> = (0..=seq.len() - window_size).step_by(step_size).collect();

    let mut raw_windows: Vec<SlidingWindowRecord> = starts
        .par_iter()
        .enumerate()
        .map(|(idx, &start)| {
            let end = start + window_size;
            let slice = &seq[start..end];

            let c_cnt = slice.iter().filter(|&&b| b == b'C').count();
            let g_cnt = slice.iter().filter(|&&b| b == b'G').count();
            let mut cpg_cnt = 0usize;
            for i in 0..slice.len().saturating_sub(1) {
                if slice[i] == b'C' && slice[i + 1] == b'G' {
                    cpg_cnt += 1;
                }
            }
            let gc_percent = ((c_cnt + g_cnt) as f64 / window_size as f64) * 100.0;
            let exp_cpg = (c_cnt as f64 * g_cnt as f64) / (window_size as f64);
            let cpg_obs_exp = if exp_cpg > 0.0 {
                cpg_cnt as f64 / exp_cpg
            } else {
                0.0
            };

            let shannon_h3 = kmer_shannon_entropy(slice, 3);
            let win_probs = base_distribution(slice);
            let js_div = jensen_shannon_divergence(&win_probs, &bg_probs);
            let (_, lz_norm) = lempel_ziv_76(slice);
            let dfa = compute_dna_walk_dfa(slice);
            let (g4_max, _) = scan_g4hunter(slice, 25, 1.2);
            let zdna = compute_zdna_propensity(slice);

            // Raw structural score combining compositional divergence, Non-B DNA propensity,
            // CpG density, and algorithmic compressibility deviation
            let raw_score = (js_div * 15.0)
                + (g4_max * 1.4)
                + (zdna * 4.5)
                + (cpg_obs_exp * 1.2)
                + ((1.0 - lz_norm).abs() * 3.5)
                + ((dfa - 0.5).abs() * 2.5);

            SlidingWindowRecord {
                window_id: format!("W{:03}", idx + 1),
                start,
                end,
                gc_percent: round4(gc_percent),
                cpg_obs_exp: round4(cpg_obs_exp),
                shannon_h3_bits: round4(shannon_h3),
                js_divergence_bg: round4(js_div),
                lz76_normalized: round4(lz_norm),
                dfa_alpha: round4(dfa),
                g4hunter_score: round4(g4_max),
                zdna_score: round4(zdna),
                structural_complexity_z: raw_score, // normalized to Z-score below
            }
        })
        .collect();

    // Standardize `structural_complexity_z` across all windows
    let n = raw_windows.len() as f64;
    let mean = raw_windows
        .iter()
        .map(|w| w.structural_complexity_z)
        .sum::<f64>()
        / n.max(1.0);
    let var = raw_windows
        .iter()
        .map(|w| (w.structural_complexity_z - mean).powi(2))
        .sum::<f64>()
        / n.max(1.0);
    let sd = var.sqrt().max(1e-9);

    for w in &mut raw_windows {
        w.structural_complexity_z = round4((w.structural_complexity_z - mean) / sd);
    }

    raw_windows
}

pub fn kmer_shannon_entropy(seq: &[u8], k: usize) -> f64 {
    if seq.len() < k || k == 0 {
        return 0.0;
    }
    let total = (seq.len() - k + 1) as f64;
    let mut counts: HashMap<&[u8], usize> = HashMap::new();
    for win in seq.windows(k) {
        *counts.entry(win).or_insert(0) += 1;
    }
    let mut h = 0.0;
    for &c in counts.values() {
        let p = c as f64 / total;
        if p > 0.0 {
            h -= p * p.log2();
        }
    }
    h
}

fn compute_codon_positional_metrics(seq: &[u8]) -> (f64, f64, f64, f64) {
    let mut pos_seqs = [Vec::new(), Vec::new(), Vec::new()];
    for (idx, &b) in seq.iter().enumerate() {
        pos_seqs[idx % 3].push(b);
    }
    let h1 = kmer_shannon_entropy(&pos_seqs[0], 1);
    let h2 = kmer_shannon_entropy(&pos_seqs[1], 1);
    let h3 = kmer_shannon_entropy(&pos_seqs[2], 1);
    let gc3_cnt = pos_seqs[2]
        .iter()
        .filter(|&&b| b == b'G' || b == b'C')
        .count();
    let gc3_pct = if !pos_seqs[2].is_empty() {
        (gc3_cnt as f64 / pos_seqs[2].len() as f64) * 100.0
    } else {
        0.0
    };
    (h1, h2, h3, gc3_pct)
}

/// Lempel-Ziv 1976 exhaustive complexity $c(n)$ normalized by $n / \log_4(n)$.
pub fn lempel_ziv_76(seq: &[u8]) -> (usize, f64) {
    let n = seq.len();
    if n == 0 {
        return (0, 0.0);
    }
    let mut phrases = 1usize;
    let mut i = 0usize;
    let mut k = 1usize;
    let mut l = 1usize;
    let mut k_max = 1usize;

    while u_idx(i, k, l) < n {
        if seq[i + k - 1] == seq[l + k - 1] {
            k += 1;
            if l + k > n {
                phrases += 1;
                break;
            }
        } else {
            if k > k_max {
                k_max = k;
            }
            i += 1;
            if i == l {
                phrases += 1;
                l += k_max;
                if l >= n {
                    break;
                }
                i = 0;
                k = 1;
                k_max = 1;
            } else {
                k = 1;
            }
        }
    }

    let nf = n as f64;
    let log4_n = (nf.ln() / 4.0f64.ln()).max(1.0);
    let expected = nf / log4_n;
    (phrases, phrases as f64 / expected)
}

#[inline]
fn u_idx(i: usize, k: usize, l: usize) -> usize {
    let _ = i;
    l + k - 1
}

/// Detrended Fluctuation Analysis (DFA) scaling exponent $\alpha$ on Purine (+1) / Pyrimidine (-1) DNA walk.
pub fn compute_dna_walk_dfa(seq: &[u8]) -> f64 {
    let n = seq.len();
    if n < 64 {
        return 0.5;
    }
    let mut x = Vec::with_capacity(n);
    let mut sum = 0.0;
    for &b in seq {
        let step = if b == b'A' || b == b'G' { 1.0 } else { -1.0 };
        sum += step;
        x.push(step);
    }
    let mean = sum / (n as f64);
    let mut y = Vec::with_capacity(n);
    let mut acc = 0.0;
    for &v in &x {
        acc += v - mean;
        y.push(acc);
    }

    let box_sizes = [8usize, 16, 32, 64, 128, 256];
    let mut log_s = Vec::new();
    let mut log_f = Vec::new();

    for &s in &box_sizes {
        if s * 2 > n {
            break;
        }
        let num_boxes = n / s;
        let mut f2_sum = 0.0;
        for b_idx in 0..num_boxes {
            let slice = &y[b_idx * s..(b_idx + 1) * s];
            // Linear least-squares trend inside box
            let sf = s as f64;
            let mean_t = (sf - 1.0) / 2.0;
            let mean_y = slice.iter().sum::<f64>() / sf;
            let mut num = 0.0;
            let mut den = 0.0;
            for (t, &yv) in slice.iter().enumerate() {
                let dt = (t as f64) - mean_t;
                num += dt * (yv - mean_y);
                den += dt * dt;
            }
            let slope = if den > 0.0 { num / den } else { 0.0 };
            let intercept = mean_y - slope * mean_t;
            for (t, &yv) in slice.iter().enumerate() {
                let trend = intercept + slope * (t as f64);
                let resid = yv - trend;
                f2_sum += resid * resid;
            }
        }
        let rms = (f2_sum / ((num_boxes * s) as f64)).sqrt();
        if rms > 0.0 {
            log_s.push((s as f64).ln());
            log_f.push(rms.ln());
        }
    }

    linear_slope(&log_s, &log_f).unwrap_or(0.5)
}

/// Rescaled-Range (R/S) Hurst Exponent $H$ on Purine/Pyrimidine DNA walk.
pub fn compute_dna_walk_hurst(seq: &[u8]) -> f64 {
    let n = seq.len();
    if n < 32 {
        return 0.5;
    }
    let window_sizes = [16usize, 32, 64, 128, 256, 512];
    let mut log_n = Vec::new();
    let mut log_rs = Vec::new();

    for &w in &window_sizes {
        if w > n {
            break;
        }
        let chunks = n / w;
        let mut rs_sum = 0.0;
        let mut valid_chunks = 0usize;
        for c in 0..chunks {
            let sub = &seq[c * w..(c + 1) * w];
            let vals: Vec<f64> = sub
                .iter()
                .map(|&b| if b == b'A' || b == b'G' { 1.0 } else { -1.0 })
                .collect();
            let mean = vals.iter().sum::<f64>() / (w as f64);
            let mut cum = 0.0;
            let mut min_c = 0.0f64;
            let mut max_c = 0.0f64;
            let mut var_sum = 0.0;
            for &v in &vals {
                let d = v - mean;
                cum += d;
                min_c = min_c.min(cum);
                max_c = max_c.max(cum);
                var_sum += d * d;
            }
            let r = max_c - min_c;
            let s = (var_sum / (w as f64)).sqrt();
            if s > 1e-9 && r > 0.0 {
                rs_sum += r / s;
                valid_chunks += 1;
            }
        }
        if valid_chunks > 0 {
            log_n.push((w as f64).ln());
            log_rs.push((rs_sum / (valid_chunks as f64)).ln());
        }
    }
    linear_slope(&log_n, &log_rs).unwrap_or(0.5)
}

/// G4Hunter sliding-window G-Quadruplex propensity scorer (Bedrat et al., Nucleic Acids Res 2016).
pub fn scan_g4hunter(seq: &[u8], window: usize, threshold: f64) -> (f64, usize) {
    if seq.len() < window {
        return (0.0, 0);
    }
    let mut per_base_score = vec![0.0f64; seq.len()];
    let mut i = 0;
    while i < seq.len() {
        let b = seq[i];
        if b == b'G' || b == b'C' {
            let mut run = 1usize;
            while i + run < seq.len() && seq[i + run] == b {
                run += 1;
            }
            let level = (run.min(4)) as f64 * if b == b'G' { 1.0 } else { -1.0 };
            for j in i..i + run {
                per_base_score[j] = level;
            }
            i += run;
        } else {
            i += 1;
        }
    }

    let mut max_abs = 0.0f64;
    let mut hits = 0usize;
    let mut cur_sum: f64 = per_base_score[..window].iter().sum();
    for start in 0..=seq.len() - window {
        if start > 0 {
            cur_sum += per_base_score[start + window - 1] - per_base_score[start - 1];
        }
        let score = (cur_sum / (window as f64)).abs();
        if score > max_abs {
            max_abs = score;
        }
        if score >= threshold {
            hits += 1;
        }
    }
    (max_abs, hits)
}

/// Z-DNA alternating purine-pyrimidine (CpG / TpG / CpA dinucleotide step) propensity score.
pub fn compute_zdna_propensity(seq: &[u8]) -> f64 {
    if seq.len() < 12 {
        return 0.0;
    }
    let mut alternating_runs = 0usize;
    let mut cur_run = 0usize;
    for i in 0..seq.len() - 1 {
        let is_pur1 = seq[i] == b'A' || seq[i] == b'G';
        let is_pur2 = seq[i + 1] == b'A' || seq[i + 1] == b'G';
        let is_gc_pair =
            (seq[i] == b'C' && seq[i + 1] == b'G') || (seq[i] == b'G' && seq[i + 1] == b'C');
        if is_pur1 != is_pur2 {
            cur_run += if is_gc_pair { 2 } else { 1 };
            if cur_run >= 8 {
                alternating_runs += 1;
            }
        } else {
            cur_run = 0;
        }
    }
    (alternating_runs as f64) / ((seq.len() as f64) / 100.0).max(1.0)
}

fn base_distribution(seq: &[u8]) -> [f64; 4] {
    let mut cnt = [0usize; 4];
    for &b in seq {
        match b {
            b'A' => cnt[0] += 1,
            b'C' => cnt[1] += 1,
            b'G' => cnt[2] += 1,
            b'T' => cnt[3] += 1,
            _ => {}
        }
    }
    let total = (cnt.iter().sum::<usize>() as f64).max(1.0);
    [
        cnt[0] as f64 / total,
        cnt[1] as f64 / total,
        cnt[2] as f64 / total,
        cnt[3] as f64 / total,
    ]
}

fn jensen_shannon_divergence(p: &[f64; 4], q: &[f64; 4]) -> f64 {
    let mut m = [0.0; 4];
    for i in 0..4 {
        m[i] = 0.5 * (p[i] + q[i]);
    }
    0.5 * kl_div(p, &m) + 0.5 * kl_div(q, &m)
}

fn kl_div(p: &[f64; 4], m: &[f64; 4]) -> f64 {
    let mut sum = 0.0;
    for i in 0..4 {
        if p[i] > 0.0 && m[i] > 0.0 {
            sum += p[i] * (p[i] / m[i]).log2();
        }
    }
    sum
}

fn linear_slope(xs: &[f64], ys: &[f64]) -> Option<f64> {
    if xs.len() < 2 || xs.len() != ys.len() {
        return None;
    }
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let mut num = 0.0;
    let mut den = 0.0;
    for (&x, &y) in xs.iter().zip(ys.iter()) {
        num += (x - mx) * (y - my);
        den += (x - mx) * (x - mx);
    }
    if den.abs() > 1e-12 {
        Some(num / den)
    } else {
        None
    }
}

#[inline]
pub fn round4(val: f64) -> f64 {
    (val * 10000.0).round() / 10000.0
}
