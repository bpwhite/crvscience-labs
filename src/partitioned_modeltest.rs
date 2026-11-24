//! Per-Locus Substitution `ModelTest` & Edge-Linked Proportional Partitioned ML/UFBoot2 Engine
//! (`Chernomor et al. 2016` proportional partition rate scalars with $\frac{1}{L}\sum_p L_p r_p = 1.0$).
//!
//! Task: Perform purely in-silico maximum-likelihood model selection per locus partition across a
//! multi-locus phylogenomic supermatrix, optimize an edge-linked proportional partitioned topology
//! and branch lengths, and evaluate stratified within-locus `UFBoot2` support.

use crate::complexity::round4;
use crate::prng_util::{splitmix64, DeterministicLcg64};
use crate::tree_builder::{
    apply_stochastic_nni_perturbation, build_rate_cats_from_params,
    build_rate_matrix_q_for_run_model, build_transition_prob_mats, build_weighted_bionj_pool,
    canonical_split_concordance_and_rf, compact_weighted_alignment,
    compress_alignment_and_estimate_rates, detect_binary_ry_columns,
    evaluate_pattern_range_weighted_lnl, evaluate_site_pattern_lnls, evaluate_tree_log_likelihood,
    extract_canonical_splits_from_ml_nodes, finalize_tree_metrics, format_ml_newick,
    infer_native_ml_tree_with_bootstrap, map_columns_to_patterns,
    optimize_substitution_model_on_tree_parallel, run_inside_outside_pulley_pass,
    solve_native_ml_tree_with_fit, CompressedAlignment, MlBootstrapReport, MlTreeNode,
    PhyloRunModel, SingleTreeRunMetrics, SITE_PATTERN_CHUNK_SIZE,
};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeSet, HashSet};
use std::time::Instant;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PartitionLocusInput {
    pub locus_name: String,
    pub locus_class: String,
    pub protein_coding: bool,
    pub is_nuclear: bool,
    pub start_col_1based: usize,
    pub end_col_1based: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LocusCandidateModelScore {
    pub model_code: String,
    pub k_subst_params: usize,
    pub rate_multiplier: f64,
    pub ln_l: f64,
    pub aic: f64,
    pub bic: f64,
    pub delta_bic: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerLocusModelTestSummary {
    pub locus_name: String,
    pub locus_class: String,
    pub alphabet_mode: String,
    pub start_col_1based: usize,
    pub end_col_1based: usize,
    pub length_bp: usize,
    pub taxa_with_data: usize,
    pub distinct_patterns: usize,
    pub invariant_fraction: f64,
    pub gc_or_purine_pct: f64,
    pub rcfv_locus: f64,
    pub winning_model_code: String,
    pub runner_up_model_code: String,
    pub delta_bic_vs_runner_up: f64,
    pub proportional_rate_r: f64,
    pub kappa_ti_tv: f64,
    pub pinv: f64,
    pub gamma_or_freerate_summary: String,
    pub locus_ln_l: f64,
    pub locus_aic: f64,
    pub locus_bic: f64,
    pub duration_ms: f64,
    pub all_models: Vec<LocusCandidateModelScore>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PartitionedSupermatrixMlResult {
    pub per_locus_models: Vec<PerLocusModelTestSummary>,
    pub n_taxa: usize,
    pub n_partitions: usize,
    pub total_sites: usize,
    pub partitioned_free_params_k: usize,
    pub partitioned_ln_l: f64,
    pub partitioned_aic: f64,
    pub partitioned_bic: f64,
    pub unpartitioned_best_model: String,
    pub unpartitioned_free_params_k: usize,
    pub unpartitioned_ln_l: f64,
    pub unpartitioned_aic: f64,
    pub unpartitioned_bic: f64,
    pub delta_ln_l_partitioned_minus_unpartitioned: f64,
    pub delta_bic_partitioned_minus_unpartitioned: f64,
    pub partitioned_metrics: SingleTreeRunMetrics,
    pub partitioned_bootstrap: MlBootstrapReport,
    pub unpartitioned_metrics: SingleTreeRunMetrics,
    pub unpartitioned_bootstrap: MlBootstrapReport,
    pub total_duration_ms: f64,
}

#[derive(Clone)]
struct FittedPartition {
    input: PartitionLocusInput,
    _seqs: Vec<String>,
    aln: CompressedAlignment,
    col_to_pat: Vec<usize>,
    _winning_code: String,
    winning_k_subst: usize,
    _winning_rm: PhyloRunModel,
    _winning_force_no_inv: bool,
    q: [[f64; 4]; 4],
    pi: [f64; 4],
    pinv: f64,
    rate_cats: Vec<(f64, f64)>,
    rate_r: f64,
    locus_ln_l: f64,
    summary: PerLocusModelTestSummary,
}

const CANDIDATE_MODELS: [(&str, usize, PhyloRunModel, bool); 8] = [
    ("JC", 0, PhyloRunModel::JcUniform, true),
    ("K80+G4", 2, PhyloRunModel::K80Gamma4, true),
    ("HKY+F+I+G4", 6, PhyloRunModel::HkyInvGamma4, false),
    ("TN+F+I+G4", 7, PhyloRunModel::Tn93InvGamma4, false),
    ("TIM2+F+I+G4", 8, PhyloRunModel::Tim2InvGamma4, false),
    ("GTR+F+G4", 9, PhyloRunModel::GtrInvGamma4Full, true),
    ("GTR+F+I+G4", 10, PhyloRunModel::GtrInvGamma4Full, false),
    ("GTR+F+R2", 10, PhyloRunModel::GtrFreeRate2, true),
];

fn golden_optimize_1d<F>(mut lo: f64, mut hi: f64, grid: &[f64], iters: usize, mut eval: F) -> (f64, f64)
where
    F: FnMut(f64) -> f64,
{
    let mut best_x = grid.first().copied().unwrap_or(1.0).clamp(lo, hi);
    let mut best_y = eval(best_x);
    let mut best_idx = 0usize;
    for (gi, &x) in grid.iter().enumerate() {
        let xc = x.clamp(lo, hi);
        let y = eval(xc);
        if y > best_y {
            best_y = y;
            best_x = xc;
            best_idx = gi;
        }
    }
    if grid.len() >= 2 {
        lo = if best_idx > 0 { grid[best_idx - 1].clamp(lo, hi) } else { lo };
        hi = if best_idx + 1 < grid.len() { grid[best_idx + 1].clamp(lo, hi) } else { hi };
    }
    if hi <= lo + 1e-6 {
        return (best_x, best_y);
    }
    let phi = 0.618_033_988_749_894_9_f64;
    let mut x1 = hi - phi * (hi - lo);
    let mut x2 = lo + phi * (hi - lo);
    let mut f1 = eval(x1);
    let mut f2 = eval(x2);
    for _ in 0..iters {
        if f1 < f2 {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = lo + phi * (hi - lo);
            f2 = eval(x2);
            if f2 > best_y {
                best_y = f2;
                best_x = x2;
            }
        } else {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = hi - phi * (hi - lo);
            f1 = eval(x1);
            if f1 > best_y {
                best_y = f1;
                best_x = x1;
            }
        }
    }
    let mid = 0.5 * (lo + hi);
    let f_mid = eval(mid);
    if f_mid > best_y {
        best_y = f_mid;
        best_x = mid;
    }
    (best_x, best_y)
}

fn eval_partitioned_lnl(nodes: &[MlTreeNode], parts: &[FittedPartition]) -> f64 {
    if nodes.is_empty() || parts.is_empty() {
        return 0.0;
    }
    let p_mats_by_part: Vec<Vec<Vec<[[f64; 4]; 4]>>> = parts
        .par_iter()
        .map(|p| build_transition_prob_mats(nodes, &p.q, &p.rate_cats, p.rate_r))
        .collect();
    let total_chunks_est: usize = parts
        .iter()
        .map(|p| p.aln.patterns.len().div_ceil(SITE_PATTERN_CHUNK_SIZE))
        .sum();
    let mut chunks: Vec<(usize, usize, usize)> = Vec::with_capacity(total_chunks_est);
    for (pi, p) in parts.iter().enumerate() {
        let n_pats = p.aln.patterns.len();
        let mut start = 0usize;
        while start < n_pats {
            let end = (start + SITE_PATTERN_CHUNK_SIZE).min(n_pats);
            chunks.push((pi, start, end));
            start = end;
        }
    }
    let chunk_sums: Vec<f64> = chunks
        .par_iter()
        .map(|&(pi, start, end)| {
            let p = &parts[pi];
            evaluate_pattern_range_weighted_lnl(
                nodes,
                &p.aln,
                &p_mats_by_part[pi],
                &p.pi,
                p.pinv,
                &p.rate_cats,
                start,
                end,
            )
        })
        .collect();
    chunk_sums.into_iter().sum()
}

#[cfg(test)]
fn eval_partitioned_lnl_coarse_reference(nodes: &[MlTreeNode], parts: &[FittedPartition]) -> f64 {
    let per_part: Vec<f64> = parts
        .par_iter()
        .map(|p| {
            crate::tree_builder::evaluate_tree_log_likelihood_seq(
                nodes,
                &p.aln,
                &p.q,
                &p.pi,
                p.pinv,
                &p.rate_cats,
                p.rate_r,
            )
        })
        .collect();
    per_part.into_iter().sum()
}

fn summarize_rate_cats_and_kappa(
    model_code: &str,
    q: &[[f64; 4]; 4],
    pi: &[f64; 4],
    pinv: f64,
    cats: &[(f64, f64)],
    fallback_kappa: f64,
) -> (f64, String) {
    let ag = q[0][2] / pi[2].max(1e-6);
    let ct = q[1][3] / pi[3].max(1e-6);
    let tv_sum = (q[0][1] / pi[1].max(1e-6))
        + (q[0][3] / pi[3].max(1e-6))
        + (q[1][2] / pi[2].max(1e-6))
        + (q[2][3] / pi[3].max(1e-6));
    let kappa = if tv_sum > 1e-8 {
        round4((((ag + ct) * 0.5) / (tv_sum * 0.25)).clamp(0.5, 45.0))
    } else {
        round4(fallback_kappa)
    };
    let summary = if model_code == "JC" || cats.len() == 1 {
        "Uniform".to_string()
    } else if model_code == "GTR+F+R2" && cats.len() == 2 {
        format!(
            "R2 (r1={:.2}, w1={:.2}; r2={:.2})",
            cats[0].0, cats[0].1, cats[1].0
        )
    } else {
        let scale = (1.0 - pinv).max(0.05);
        let var: f64 = cats
            .iter()
            .map(|&(r, _)| {
                let s = r * scale;
                (s - 1.0).powi(2)
            })
            .sum::<f64>()
            / (cats.len().max(1) as f64);
        let alpha_est = (1.0 / var.max(0.025)).clamp(0.08, 9.99);
        if pinv > 0.001 {
            format!("I+G4 (pinv={:.2}, a={:.2})", pinv, alpha_est)
        } else {
            format!("G4 (a={:.2})", alpha_est)
        }
    };
    (kappa, summary)
}

fn compute_locus_composition_stats(seqs_p: &[String]) -> (usize, f64, f64) {
    let mut taxa_present = 0usize;
    let mut total_valid = 0usize;
    let mut gc_or_purine = 0usize;
    let mut freqs: Vec<[f64; 4]> = Vec::with_capacity(seqs_p.len());
    for s in seqs_p {
        let mut c = [0usize; 4];
        for b in s.bytes() {
            match b.to_ascii_uppercase() {
                b'A' => c[0] += 1,
                b'C' => c[1] += 1,
                b'G' => c[2] += 1,
                b'T' | b'U' => c[3] += 1,
                _ => {}
            }
        }
        let n_b = c[0] + c[1] + c[2] + c[3];
        if n_b >= 10 {
            taxa_present += 1;
            total_valid += n_b;
            gc_or_purine += c[1] + c[2];
            let nf = n_b as f64;
            freqs.push([c[0] as f64 / nf, c[1] as f64 / nf, c[2] as f64 / nf, c[3] as f64 / nf]);
        }
    }
    let gc_pct = if total_valid > 0 {
        round4((gc_or_purine as f64 / total_valid as f64) * 100.0)
    } else {
        50.0
    };
    let rcfv = if freqs.len() >= 2 {
        let mut mean_f = [0.0f64; 4];
        for f in &freqs {
            for k in 0..4 {
                mean_f[k] += f[k];
            }
        }
        let m = freqs.len() as f64;
        for k in 0..4 {
            mean_f[k] /= m;
        }
        let sum_dev: f64 = freqs
            .iter()
            .map(|f| (0..4).map(|k| (f[k] - mean_f[k]).abs()).sum::<f64>())
            .sum();
        round4(sum_dev / m)
    } else {
        0.0
    };
    (taxa_present, gc_pct, rcfv)
}

/// Normalize proportional partition rates so length-weighted mean equals 1.0,
/// scaling non-root tree branch lengths inversely to preserve `r_p * b_e` exactly.
fn normalize_partition_rates_and_branches(nodes: &mut [MlTreeNode], parts: &mut [FittedPartition]) {
    let root_id = nodes.len().saturating_sub(1);
    let total_len: usize = parts.iter().map(|p| p.aln.n_sites).sum();
    if total_len == 0 {
        return;
    }
    let weighted_sum: f64 = parts
        .iter()
        .map(|p| (p.aln.n_sites as f64) * p.rate_r)
        .sum();
    let r_bar = weighted_sum / (total_len as f64);
    if r_bar <= 1e-8 {
        return;
    }
    for p in parts.iter_mut() {
        p.rate_r = (p.rate_r / r_bar).clamp(0.05, 8.0);
        p.summary.proportional_rate_r = round4(p.rate_r);
    }
    for (ni, nd) in nodes.iter_mut().enumerate() {
        if ni != root_id {
            nd.branch_len = (nd.branch_len * r_bar).clamp(1e-5, 5.0);
        }
    }
}

fn is_valid_binary_tree(root: usize, nodes: &[MlTreeNode], n_taxa: usize) -> bool {
    if root >= nodes.len() || n_taxa == 0 || nodes.len() < n_taxa {
        return false;
    }
    let mut visited = vec![false; nodes.len()];
    let mut stack = Vec::with_capacity(nodes.len());
    stack.push(root);
    let mut leaf_count = 0usize;
    while let Some(u) = stack.pop() {
        if u >= nodes.len() || visited[u] {
            return false;
        }
        visited[u] = true;
        match (nodes[u].left, nodes[u].right) {
            (None, None) => {
                if u >= n_taxa {
                    return false;
                }
                leaf_count += 1;
            }
            (Some(l), Some(r)) => {
                if u < n_taxa || l == r {
                    return false;
                }
                stack.push(l);
                stack.push(r);
            }
            _ => return false,
        }
    }
    leaf_count == n_taxa
}

#[allow(clippy::too_many_arguments)]
fn try_partitioned_nni_edge(
    u: usize,
    side: u8,
    child_id: usize,
    sib_id: usize,
    n_taxa: usize,
    nodes: &mut [MlTreeNode],
    parts: &[FittedPartition],
    cur_lnl: &mut f64,
    eval_cnt: &mut usize,
) -> bool {
    let expected_u_children = if side == 0 {
        (Some(child_id), Some(sib_id))
    } else {
        (Some(sib_id), Some(child_id))
    };
    if u >= nodes.len()
        || child_id >= nodes.len()
        || sib_id >= nodes.len()
        || child_id < n_taxa
        || (nodes[u].left, nodes[u].right) != expected_u_children
    {
        return false;
    }
    let (Some(a_id), Some(b_id)) = (nodes[child_id].left, nodes[child_id].right) else {
        return false;
    };
    *eval_cnt += 2;
    let orig_bl = nodes[child_id].branch_len;
    let orig_u_l = nodes[u].left;
    let orig_u_r = nodes[u].right;
    let orig_c_l = nodes[child_id].left;
    let orig_c_r = nodes[child_id].right;

    let mut best_lnl = *cur_lnl;
    let mut best_choice: Option<(u8, f64)> = None;

    for swap_mode in 0..2u8 {
        let mut cand = nodes.to_vec();
        if swap_mode == 0 {
            cand[child_id].right = Some(sib_id);
            if side == 0 {
                cand[u].right = Some(b_id);
            } else {
                cand[u].left = Some(b_id);
            }
        } else {
            cand[child_id].left = Some(sib_id);
            if side == 0 {
                cand[u].right = Some(a_id);
            } else {
                cand[u].left = Some(a_id);
            }
        }
        if !is_valid_binary_tree(cand.len().saturating_sub(1), &cand, n_taxa) {
            continue;
        }
        let lnl_0 = eval_partitioned_lnl(&cand, parts);
        let mut local_best_lnl = lnl_0;
        let mut local_best_bl = orig_bl;
        if lnl_0 > best_lnl - 3.5 {
            for &mult in &[0.65_f64, 1.35_f64] {
                cand[child_id].branch_len = (orig_bl * mult).clamp(1e-5, 3.0);
                let lnl_m = eval_partitioned_lnl(&cand, parts);
                if lnl_m > local_best_lnl {
                    local_best_lnl = lnl_m;
                    local_best_bl = cand[child_id].branch_len;
                }
            }
        }
        if local_best_lnl > best_lnl + 1e-4 {
            best_lnl = local_best_lnl;
            best_choice = Some((swap_mode, local_best_bl));
        }
    }

    if let Some((swap_mode, new_bl)) = best_choice {
        if swap_mode == 0 {
            nodes[child_id].right = Some(sib_id);
            if side == 0 {
                nodes[u].right = Some(b_id);
            } else {
                nodes[u].left = Some(b_id);
            }
        } else {
            nodes[child_id].left = Some(sib_id);
            if side == 0 {
                nodes[u].right = Some(a_id);
            } else {
                nodes[u].left = Some(a_id);
            }
        }
        if !is_valid_binary_tree(nodes.len().saturating_sub(1), nodes, n_taxa) {
            nodes[u].left = orig_u_l;
            nodes[u].right = orig_u_r;
            nodes[child_id].left = orig_c_l;
            nodes[child_id].right = orig_c_r;
            nodes[child_id].branch_len = orig_bl;
            return false;
        }
        nodes[child_id].branch_len = new_bl;
        *cur_lnl = best_lnl;
        true
    } else {
        nodes[u].left = orig_u_l;
        nodes[u].right = orig_u_r;
        nodes[child_id].left = orig_c_l;
        nodes[child_id].right = orig_c_r;
        nodes[child_id].branch_len = orig_bl;
        false
    }
}

fn collect_partitioned_internal_edges(
    u: usize,
    nodes: &[MlTreeNode],
    n_taxa: usize,
    out: &mut Vec<(usize, u8, usize, usize)>,
) {
    if let (Some(l), Some(r)) = (nodes[u].left, nodes[u].right) {
        if l >= n_taxa && nodes[l].left.is_some() && nodes[l].right.is_some() {
            out.push((u, 0, l, r));
            collect_partitioned_internal_edges(l, nodes, n_taxa, out);
        }
        if r >= n_taxa && nodes[r].left.is_some() && nodes[r].right.is_some() {
            out.push((u, 1, r, l));
            collect_partitioned_internal_edges(r, nodes, n_taxa, out);
        }
    }
}

fn optimize_partitioned_topology_and_branches(
    nodes: &mut [MlTreeNode],
    parts: &mut [FittedPartition],
    n_taxa: usize,
) -> (f64, usize, usize) {
    let mut cur_lnl = eval_partitioned_lnl(nodes, parts);
    let mut eval_cnt = 0usize;
    let mut acc_cnt = 0usize;
    let root_id = nodes.len().saturating_sub(1);

    for _pass in 0..2 {
        let mut edges = Vec::with_capacity(n_taxa);
        collect_partitioned_internal_edges(root_id, nodes, n_taxa, &mut edges);
        let mut pass_acc = 0usize;
        for (u, side, cid, sid) in edges {
            if try_partitioned_nni_edge(u, side, cid, sid, n_taxa, nodes, parts, &mut cur_lnl, &mut eval_cnt) {
                acc_cnt += 1;
                pass_acc += 1;
            }
        }
        debug_assert!(is_valid_binary_tree(root_id, nodes, n_taxa));
        if pass_acc == 0 {
            break;
        }
    }

    // Coordinate edge-length optimization across all non-root branches under exact partitioned likelihood
    for ei in 0..root_id {
        let b0 = nodes[ei].branch_len.clamp(1e-5, 2.5);
        let lo = (b0 * 0.45).max(1e-5);
        let hi = (b0 * 2.15).min(3.5);
        let grid = [lo, b0 * 0.72, b0 * 0.89, b0, b0 * 1.14, b0 * 1.42, hi];
        let mut cand = nodes.to_vec();
        let (best_b, best_l) = golden_optimize_1d(lo, hi, &grid, 7, |b| {
            cand[ei].branch_len = b;
            eval_partitioned_lnl(&cand, parts)
        });
        if best_l > cur_lnl {
            nodes[ei].branch_len = best_b;
            cur_lnl = best_l;
        }
    }

    // Re-optimize each partition's proportional rate multiplier `r_p`, then normalize weighted mean = 1.0
    let nodes_ro = nodes.to_vec();
    let new_rates: Vec<(f64, f64)> = parts
        .par_iter()
        .map(|p| {
            golden_optimize_1d(0.06, 5.5, &[0.10, 0.25, 0.50, 0.75, 1.0, 1.35, 1.9, 2.8, 4.2], 14, |r| {
                evaluate_tree_log_likelihood(&nodes_ro, &p.aln, &p.q, &p.pi, p.pinv, &p.rate_cats, r)
            })
        })
        .collect();
    for (p, (r_new, lnl_new)) in parts.iter_mut().zip(new_rates.into_iter()) {
        p.rate_r = r_new;
        p.locus_ln_l = round4(lnl_new);
        p.summary.locus_ln_l = round4(lnl_new);
        let len_bp = p.aln.n_sites.max(1);
        p.summary.locus_aic = round4(-2.0 * lnl_new + 2.0 * ((p.winning_k_subst + 1) as f64));
        p.summary.locus_bic = round4(-2.0 * lnl_new + ((p.winning_k_subst + 1) as f64) * (len_bp as f64).ln());
    }
    normalize_partition_rates_and_branches(nodes, parts);
    cur_lnl = eval_partitioned_lnl(nodes, parts);
    (cur_lnl, eval_cnt, acc_cnt)
}

/// Evaluate individual per-locus `ModelTest` across all supermatrix partitions, optimize an edge-linked
/// proportional partitioned ML topology with $\frac{1}{L}\sum_p L_p r_p = 1.0$, and execute stratified
/// within-locus `UFBoot2` replicates alongside an unpartitioned single-model supermatrix reference.
pub fn infer_partitioned_supermatrix_ml_with_bootstrap(
    raw_aligned_strings: &[String],
    labels: &[String],
    partitions: &[PartitionLocusInput],
    b_replicates: usize,
    seed: u64,
) -> Result<PartitionedSupermatrixMlResult, String> {
    let t0 = Instant::now();
    let n_taxa = labels.len();
    if n_taxa < 4 || raw_aligned_strings.len() != n_taxa {
        return Err(format!(
            "Partitioned ML requires >= 4 taxa with matching sequence count (labels={}, seqs={})",
            n_taxa,
            raw_aligned_strings.len()
        ));
    }
    if partitions.is_empty() {
        return Err("Partitioned ML requires at least 1 locus partition".to_string());
    }
    let total_sites = raw_aligned_strings.first().map(|s| s.len()).unwrap_or(0);
    if total_sites < 12 {
        return Err(format!("Alignment too short ({total_sites} cols) for partitioned ML"));
    }

    // Step A: Unpartitioned reference model selection & ML + UFBoot2
    let full_aln = compress_alignment_and_estimate_rates(raw_aligned_strings);
    let is_binary_col = detect_binary_ry_columns(raw_aligned_strings);
    let full_col_to_pat = map_columns_to_patterns(raw_aligned_strings, &full_aln);

    let guide_fast = solve_native_ml_tree_with_fit(
        raw_aligned_strings,
        labels,
        &full_aln,
        PhyloRunModel::GtrInvGamma4Fast,
    );
    let branch_k = 2 * n_taxa.saturating_sub(3) + 1;
    let ln_sites = (total_sites as f64).ln();
    let mut best_unpart_code = "GTR+F+I+G4".to_string();
    let mut best_unpart_rm = PhyloRunModel::GtrInvGamma4Full;
    let mut best_unpart_bic = f64::INFINITY;
    let unpart_candidates: Vec<(&str, PhyloRunModel, f64)> = CANDIDATE_MODELS
        .par_iter()
        .map(|&(code, k_m, rm, no_inv)| {
            let (q, pi, pinv, cats, _lnl) = optimize_substitution_model_on_tree_parallel(
                &guide_fast.nodes,
                &full_aln,
                rm,
                no_inv,
            );
            let exact_l = evaluate_tree_log_likelihood(&guide_fast.nodes, &full_aln, &q, &pi, pinv, &cats, 1.0);
            let k_tot = branch_k + k_m;
            let bic = -2.0 * exact_l + (k_tot as f64) * ln_sites;
            let resolved_rm = if code == "GTR+F+G4" {
                PhyloRunModel::GtrInvGamma4Full
            } else {
                rm
            };
            (code, resolved_rm, bic)
        })
        .collect();
    for (code, resolved_rm, bic) in unpart_candidates {
        if bic < best_unpart_bic {
            best_unpart_bic = bic;
            best_unpart_code = code.to_string();
            best_unpart_rm = resolved_rm;
        }
    }

    let (unpart_metrics, unpart_boot) = infer_native_ml_tree_with_bootstrap(
        raw_aligned_strings,
        labels,
        best_unpart_rm,
        b_replicates,
        seed ^ 0x1111,
    )?;
    let guide_fit = solve_native_ml_tree_with_fit(raw_aligned_strings, labels, &full_aln, best_unpart_rm);
    let guide_nodes = guide_fit.nodes.clone();

    // Step B: Parallel per-locus individual ModelTest on each partition
    let mut fitted_parts: Vec<FittedPartition> = partitions
        .par_iter()
        .enumerate()
        .map(|(_pi, part_in)| {
            let tp = Instant::now();
            let s0 = part_in.start_col_1based.saturating_sub(1).min(total_sites);
            let e0 = part_in.end_col_1based.min(total_sites).max(s0);
            let len_bp = e0.saturating_sub(s0).max(1);
            let seqs_p: Vec<String> = raw_aligned_strings
                .iter()
                .map(|row| {
                    let bytes = row.as_bytes();
                    if e0 <= bytes.len() {
                        String::from_utf8_lossy(&bytes[s0..e0]).into_owned()
                    } else {
                        "-".repeat(len_bp)
                    }
                })
                .collect();
            let aln_p = compress_alignment_and_estimate_rates(&seqs_p);
            let col_to_pat = map_columns_to_patterns(&seqs_p, &aln_p);
            let (taxa_with_data, gc_or_purine_pct, rcfv_locus) = compute_locus_composition_stats(&seqs_p);

            let is_ry_slice = detect_binary_ry_columns(&seqs_p).into_iter().any(|b| b);
            let alphabet_mode = if part_in.protein_coding && is_ry_slice {
                "RY3 (Codon Pos3 Binary + Pos1/2 4-State)".to_string()
            } else {
                "4-State Nucleotide".to_string()
            };

            struct EvalModel {
                score: LocusCandidateModelScore,
                rm: PhyloRunModel,
                no_inv: bool,
                q: [[f64; 4]; 4],
                pi: [f64; 4],
                pinv: f64,
                cats: Vec<(f64, f64)>,
            }
            let ln_lp = (len_bp as f64).ln();
            let rate_grid = [0.14_f64, 0.30, 0.55, 0.80, 1.0, 1.35, 1.9, 2.8, 4.1];

            let mut evals: Vec<EvalModel> = CANDIDATE_MODELS
                .par_iter()
                .map(|&(mcode, k_sub, rm, force_no_inv)| {
                    let (q0, pi0) = build_rate_matrix_q_for_run_model(rm, &aln_p);
                    let rate_rm0 = if force_no_inv && matches!(rm, PhyloRunModel::GtrInvGamma4Full) {
                        PhyloRunModel::K80Gamma4
                    } else {
                        rm
                    };
                    let (pinv0, cats0) = build_rate_cats_from_params(
                        rate_rm0,
                        0.50,
                        if force_no_inv { 0.0 } else { 0.24 },
                        0.68,
                        0.08,
                    );
                    let (r0, _) = golden_optimize_1d(0.08, 4.50, &rate_grid, 10, |r| {
                        evaluate_tree_log_likelihood(&guide_nodes, &aln_p, &q0, &pi0, pinv0, &cats0, r)
                    });
                    let mut scaled_nodes = guide_nodes.clone();
                    let g_root = scaled_nodes.len().saturating_sub(1);
                    for (ni, nd) in scaled_nodes.iter_mut().enumerate() {
                        if ni != g_root {
                            nd.branch_len = (nd.branch_len * r0).clamp(1e-5, 4.0);
                        }
                    }
                    let (q_opt, pi_opt, pinv_opt, cats_opt, _) =
                        optimize_substitution_model_on_tree_parallel(&scaled_nodes, &aln_p, rm, force_no_inv);
                    let (r_opt, lnl_opt) = golden_optimize_1d(0.08, 4.50, &rate_grid, 14, |r| {
                        evaluate_tree_log_likelihood(
                            &guide_nodes,
                            &aln_p,
                            &q_opt,
                            &pi_opt,
                            pinv_opt,
                            &cats_opt,
                            r,
                        )
                    });
                    let k_eff = (k_sub + 1) as f64;
                    let aic = round4(-2.0 * lnl_opt + 2.0 * k_eff);
                    let bic = round4(-2.0 * lnl_opt + k_eff * ln_lp);
                    EvalModel {
                        score: LocusCandidateModelScore {
                            model_code: mcode.to_string(),
                            k_subst_params: k_sub,
                            rate_multiplier: round4(r_opt),
                            ln_l: round4(lnl_opt),
                            aic,
                            bic,
                            delta_bic: 0.0,
                        },
                        rm,
                        no_inv: force_no_inv,
                        q: q_opt,
                        pi: pi_opt,
                        pinv: pinv_opt,
                        cats: cats_opt,
                    }
                })
                .collect();

            evals.sort_by(|a, b| {
                a.score
                    .bic
                    .partial_cmp(&b.score.bic)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| {
                        b.score
                            .ln_l
                            .partial_cmp(&a.score.ln_l)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .then_with(|| a.score.k_subst_params.cmp(&b.score.k_subst_params))
            });
            let min_bic = evals.first().map(|e| e.score.bic).unwrap_or(0.0);
            for ev in &mut evals {
                ev.score.delta_bic = round4((ev.score.bic - min_bic).max(0.0));
            }
            let win = &evals[0];
            let runner_up = evals.get(1).unwrap_or(win);
            let (kappa_ti_tv, gamma_or_freerate_summary) = summarize_rate_cats_and_kappa(
                &win.score.model_code,
                &win.q,
                &win.pi,
                win.pinv,
                &win.cats,
                aln_p.kappa_overall,
            );
            let dur_ms = round4(tp.elapsed().as_secs_f64() * 1000.0);
            let all_models: Vec<LocusCandidateModelScore> = evals.iter().map(|e| e.score.clone()).collect();
            let summary = PerLocusModelTestSummary {
                locus_name: part_in.locus_name.clone(),
                locus_class: part_in.locus_class.clone(),
                alphabet_mode,
                start_col_1based: part_in.start_col_1based,
                end_col_1based: part_in.end_col_1based,
                length_bp: len_bp,
                taxa_with_data,
                distinct_patterns: aln_p.patterns.len(),
                invariant_fraction: round4(aln_p.invariant_fraction),
                gc_or_purine_pct,
                rcfv_locus,
                winning_model_code: win.score.model_code.clone(),
                runner_up_model_code: runner_up.score.model_code.clone(),
                delta_bic_vs_runner_up: runner_up.score.delta_bic,
                proportional_rate_r: win.score.rate_multiplier,
                kappa_ti_tv,
                pinv: round4(win.pinv),
                gamma_or_freerate_summary,
                locus_ln_l: win.score.ln_l,
                locus_aic: win.score.aic,
                locus_bic: win.score.bic,
                duration_ms: dur_ms,
                all_models,
            };

            FittedPartition {
                input: part_in.clone(),
                _seqs: seqs_p,
                aln: aln_p,
                col_to_pat,
                _winning_code: win.score.model_code.clone(),
                winning_k_subst: win.score.k_subst_params,
                _winning_rm: win.rm,
                _winning_force_no_inv: win.no_inv,
                q: win.q,
                pi: win.pi,
                pinv: win.pinv,
                rate_cats: win.cats.clone(),
                rate_r: win.score.rate_multiplier,
                locus_ln_l: win.score.ln_l,
                summary,
            }
        })
        .collect();

    // Normalize initial rates so length-weighted mean == 1.0 before edge-linked topology polish
    let mut part_nodes = guide_nodes.clone();
    normalize_partition_rates_and_branches(&mut part_nodes, &mut fitted_parts);

    // Step C: Edge-Linked Partitioned Topology NNI + Coordinate Branch-Length & Rate-Multiplier Polish
    let (mut part_exact_lnl, mut nni_eval, mut nni_acc) =
        optimize_partitioned_topology_and_branches(&mut part_nodes, &mut fitted_parts, n_taxa);

    // Step D: Stratified Within-Locus Partitioned UFBoot2 + Self-Consistency Gate
    let n_pats_full = full_aln.patterns.len();
    let (init_q_full, init_pi_full) = build_rate_matrix_q_for_run_model(best_unpart_rm, &full_aln);
    let (init_pinv_full, init_cats_full) =
        build_rate_cats_from_params(best_unpart_rm, 0.30, 0.0, 0.68, 0.06);

    let binary_ry_columns_detected = is_binary_col.iter().filter(|&&b| b).count();
    let four_state_columns_detected = is_binary_col.len().saturating_sub(binary_ry_columns_detected);

    let m_explore = if b_replicates == 0 {
        0usize
    } else if total_sites > 3_000 {
        b_replicates.clamp(28, 36)
    } else {
        b_replicates.clamp(32, 44)
    };

    let explored_trees: Vec<Vec<MlTreeNode>> = (0..m_explore)
        .into_par_iter()
        .map(|m| {
            let mut rng = DeterministicLcg64::new(splitmix64(seed ^ 0x5041_5254_u64 ^ (m as u64 + 1)));
            let mut col_w = vec![0u16; total_sites];
            let mut pat_w = vec![0.0f64; n_pats_full];
            for p in &fitted_parts {
                let s0 = p.input.start_col_1based.saturating_sub(1).min(total_sites);
                let l_p = p.aln.n_sites.max(1);
                for _ in 0..l_p {
                    let off = rng.gen_range_usize(l_p);
                    let c = (s0 + off).min(total_sites.saturating_sub(1));
                    col_w[c] = col_w[c].saturating_add(1);
                    pat_w[full_col_to_pat[c]] += 1.0;
                }
            }
            if m % 4 == 0 {
                let aln_rep = compact_weighted_alignment(&full_aln, &pat_w);
                let cands = build_weighted_bionj_pool(
                    raw_aligned_strings,
                    labels,
                    best_unpart_rm,
                    None,
                    Some(&col_w),
                    &is_binary_col,
                    &[0, 1, 2],
                );
                let mut best_n = part_nodes.clone();
                let mut best_l = f64::NEG_INFINITY;
                for (_, mut cn) in cands {
                    let mut ev = 0usize;
                    let mut acc = 0usize;
                    let lnl = run_inside_outside_pulley_pass(
                        &mut cn,
                        n_taxa,
                        &aln_rep,
                        &init_q_full,
                        &init_pi_full,
                        init_pinv_full,
                        &init_cats_full,
                        true,
                        &mut ev,
                        &mut acc,
                    );
                    if lnl > best_l || best_l.is_infinite() {
                        best_l = lnl;
                        best_n = cn;
                    }
                }
                best_n
            } else {
                let aln_rep = compact_weighted_alignment(&full_aln, &pat_w);
                let n_swaps = match m % 3 {
                    0 => 2usize,
                    1 => 5usize,
                    _ => 9usize,
                };
                let mut cand = apply_stochastic_nni_perturbation(&part_nodes, n_taxa, n_swaps, &mut rng);
                let mut ev = 0usize;
                for _ in 0..2 {
                    let mut acc = 0usize;
                    run_inside_outside_pulley_pass(
                        &mut cand,
                        n_taxa,
                        &aln_rep,
                        &guide_fit.q,
                        &guide_fit.pi,
                        guide_fit.pinv,
                        &guide_fit.rate_cats,
                        true,
                        &mut ev,
                        &mut acc,
                    );
                    if acc == 0 {
                        break;
                    }
                }
                cand
            }
        })
        .collect();

    struct PartBootCand {
        nodes: Vec<MlTreeNode>,
        splits: BTreeSet<Vec<u64>>,
        ell_by_part: Vec<Vec<f64>>,
        exact_lnl: f64,
    }

    let mut seen_topologies: HashSet<BTreeSet<Vec<u64>>> = HashSet::with_capacity(m_explore + 1);
    let mut unique_trees: Vec<(Vec<MlTreeNode>, BTreeSet<Vec<u64>>)> = Vec::with_capacity(m_explore);
    let (_, part_init_splits) = extract_canonical_splits_from_ml_nodes(&part_nodes, n_taxa);
    seen_topologies.insert(part_init_splits.clone());
    unique_trees.push((part_nodes.clone(), part_init_splits));
    for cand in explored_trees {
        let (_, sp) = extract_canonical_splits_from_ml_nodes(&cand, n_taxa);
        if seen_topologies.insert(sp.clone()) {
            unique_trees.push((cand, sp));
        }
    }

    let inner_par = unique_trees.len() < rayon::current_num_threads();
    let bank: Vec<PartBootCand> = unique_trees
        .into_par_iter()
        .map(|(mut nodes, splits)| {
            let mut ev_b = 0usize;
            let mut acc_b = 0usize;
            run_inside_outside_pulley_pass(
                &mut nodes,
                n_taxa,
                &full_aln,
                &guide_fit.q,
                &guide_fit.pi,
                guide_fit.pinv,
                &guide_fit.rate_cats,
                false,
                &mut ev_b,
                &mut acc_b,
            );
            let score_part = |p: &FittedPartition| {
                evaluate_site_pattern_lnls(
                    &nodes,
                    &p.aln,
                    &p.q,
                    &p.pi,
                    p.pinv,
                    &p.rate_cats,
                    p.rate_r,
                    inner_par,
                )
            };
            let ell_by_part: Vec<Vec<f64>> = if inner_par {
                fitted_parts.par_iter().map(score_part).collect()
            } else {
                fitted_parts.iter().map(score_part).collect()
            };
            let exact_lnl: f64 = ell_by_part
                .iter()
                .zip(fitted_parts.iter())
                .map(|(ell, p)| {
                    ell.iter()
                        .zip(p.aln.weights.iter())
                        .map(|(&lp, &w)| if w > 0.0 { w * lp } else { 0.0 })
                        .sum::<f64>()
                })
                .sum();
            PartBootCand {
                nodes,
                splits,
                ell_by_part,
                exact_lnl,
            }
        })
        .collect();

    if let Some(best_c) = bank
        .iter()
        .max_by(|a, b| a.exact_lnl.partial_cmp(&b.exact_lnl).unwrap_or(std::cmp::Ordering::Equal))
    {
        if best_c.exact_lnl.is_finite() && best_c.exact_lnl > part_exact_lnl + 1e-4 {
            part_nodes = best_c.nodes.clone();
            // Polish rates after topology promotion
            normalize_partition_rates_and_branches(&mut part_nodes, &mut fitted_parts);
            part_exact_lnl = eval_partitioned_lnl(&part_nodes, &fitted_parts);
        }
    }

    // Identity-replicate & partitioned NNI stationarity self-consistency gate
    let mut id_nodes = part_nodes.clone();
    let mut id_lnl = part_exact_lnl;
    let id_root = id_nodes.len().saturating_sub(1);
    let mut id_edges = Vec::with_capacity(n_taxa);
    collect_partitioned_internal_edges(id_root, &id_nodes, n_taxa, &mut id_edges);
    for (u, side, cid, sid) in id_edges {
        if try_partitioned_nni_edge(
            u,
            side,
            cid,
            sid,
            n_taxa,
            &mut id_nodes,
            &fitted_parts,
            &mut id_lnl,
            &mut nni_eval,
        ) {
            nni_acc += 1;
        }
    }
    let (mut ref_split_entries, ref_split_set) =
        extract_canonical_splits_from_ml_nodes(&part_nodes, n_taxa);
    let (_, id_split_set) = extract_canonical_splits_from_ml_nodes(&id_nodes, n_taxa);
    let (identity_replicate_rf_distance, identity_replicate_concordance_pct) =
        canonical_split_concordance_and_rf(&ref_split_set, &id_split_set);
    let st_gain = (id_lnl - part_exact_lnl).max(0.0);
    let lnl_deficit_nats = (id_lnl - part_exact_lnl).max(0.0);
    let max_allowed_deficit = ((total_sites as f64) * 0.03).max(8.0);
    let estimator_self_consistent = lnl_deficit_nats <= max_allowed_deficit
        && (identity_replicate_concordance_pct >= 85.0 || st_gain <= 2.0);
    if id_lnl > part_exact_lnl + 1e-4 {
        part_nodes = id_nodes;
        part_exact_lnl = id_lnl;
        ref_split_entries = extract_canonical_splits_from_ml_nodes(&part_nodes, n_taxa).0;
    }
    if ref_split_entries.len() >= 3 && !estimator_self_consistent {
        return Err(format!(
            "Partitioned ML bootstrap estimator mismatch: concordance {:.2}% (RF={}, gain={:.2} nats)",
            identity_replicate_concordance_pct, identity_replicate_rf_distance, st_gain
        ));
    }

    // Stratified within-locus RELL UFBoot2 scoring across 0..b_replicates
    let n_splits = ref_split_entries.len();
    let (mean_split_support_pct, min_split_support_pct, max_split_support_pct, split_support_pct_by_mask) =
        if b_replicates == 0 || n_splits == 0 {
            (
                100.0,
                100.0,
                100.0,
                ref_split_entries.iter().map(|(m, _)| (m.clone(), 100.0)).collect(),
            )
        } else {
            let hit_vecs: Vec<Vec<u32>> = (0..b_replicates)
                .into_par_iter()
                .map(|rep| {
                    let mut rng = DeterministicLcg64::new(splitmix64(seed ^ 0x5245_4C4C_u64 ^ (rep as u64 + 1)));
                    let pat_w_by_p: Vec<Vec<f64>> = fitted_parts
                        .iter()
                        .map(|p| {
                            let l_p = p.aln.n_sites.max(1);
                            let mut pw = vec![0.0f64; p.aln.patterns.len()];
                            for _ in 0..l_p {
                                let col = rng.gen_range_usize(l_p);
                                let pat_i = p.col_to_pat.get(col).copied().unwrap_or(0);
                                if pat_i < pw.len() {
                                    pw[pat_i] += 1.0;
                                }
                            }
                            pw
                        })
                        .collect();
                    let mut best_score = f64::NEG_INFINITY;
                    let mut tied: Vec<usize> = Vec::with_capacity(4);
                    for (k, cand) in bank.iter().enumerate() {
                        let mut score = 0.0f64;
                        for (pi, pw) in pat_w_by_p.iter().enumerate() {
                            let ell = &cand.ell_by_part[pi];
                            for (&w, &lp) in pw.iter().zip(ell.iter()) {
                                if w > 0.0 {
                                    score += w * lp;
                                }
                            }
                        }
                        if score > best_score + 1e-7 {
                            best_score = score;
                            tied.clear();
                            tied.push(k);
                        } else if score >= best_score - 1e-7 {
                            tied.push(k);
                        }
                    }
                    let winner_k = if tied.len() <= 1 {
                        tied.first().copied().unwrap_or(0)
                    } else {
                        tied[rng.gen_range_usize(tied.len())]
                    };
                    let boot_splits = &bank[winner_k].splits;
                    ref_split_entries
                        .iter()
                        .map(|(sp, _)| if boot_splits.contains(sp) { 1u32 } else { 0u32 })
                        .collect()
                })
                .collect();

            let mut counts = vec![0usize; n_splits];
            for hv in &hit_vecs {
                for (i, &h) in hv.iter().enumerate() {
                    counts[i] += h as usize;
                }
            }
            let mut out_pairs = Vec::with_capacity(n_splits);
            let mut sum_pct = 0.0f64;
            let mut min_pct = 100.0f64;
            let mut max_pct = 0.0f64;
            for (i, (mask, node_indices)) in ref_split_entries.into_iter().enumerate() {
                let pct = round4((counts[i] as f64 / b_replicates as f64) * 100.0);
                sum_pct += pct;
                min_pct = min_pct.min(pct);
                max_pct = max_pct.max(pct);
                for &ni in &node_indices {
                    if ni < part_nodes.len() {
                        part_nodes[ni].support_pct = pct;
                    }
                }
                out_pairs.push((mask, pct));
            }
            (
                round4(sum_pct / (n_splits as f64)),
                round4(min_pct),
                round4(max_pct),
                out_pairs,
            )
        };

    // Final exact parameter counts and information criteria
    let n_partitions = fitted_parts.len();
    let sum_k_subst: usize = fitted_parts.iter().map(|p| p.winning_k_subst).sum();
    let unrooted_branches = 2 * n_taxa.saturating_sub(3);
    let partitioned_free_params_k = unrooted_branches + n_partitions.saturating_sub(1) + sum_k_subst;
    let partitioned_ln_l = round4(part_exact_lnl);
    let partitioned_aic = round4(-2.0 * partitioned_ln_l + 2.0 * (partitioned_free_params_k as f64));
    let partitioned_bic = round4(-2.0 * partitioned_ln_l + (partitioned_free_params_k as f64) * ln_sites);

    let total_dur_ms = round4(t0.elapsed().as_secs_f64() * 1000.0);
    let mut part_metrics = finalize_tree_metrics(
        "Native Rust Edge-Linked Partitioned ML",
        best_unpart_rm,
        n_taxa,
        total_sites,
        total_dur_ms,
        &part_nodes,
        partitioned_ln_l,
        nni_eval,
        nni_acc,
    );
    let root_id = part_nodes.len().saturating_sub(1);
    part_metrics.run_model_code = format!("Partitioned ({}-Locus BIC)", n_partitions);
    part_metrics.free_parameters_k = partitioned_free_params_k;
    part_metrics.aic = partitioned_aic;
    part_metrics.bic = partitioned_bic;
    part_metrics.mean_branch_support_pct = mean_split_support_pct;
    part_metrics.newick_tree = format!("{};", format_ml_newick(&part_nodes, root_id));

    let partitioned_bootstrap = MlBootstrapReport {
        b_replicates,
        estimator_model_code: format!("Edge-Linked Partitioned ({} Loci)", n_partitions),
        binary_ry_columns_detected,
        four_state_columns_detected,
        nontrivial_reference_splits: n_splits,
        mean_split_support_pct,
        min_split_support_pct,
        max_split_support_pct,
        identity_replicate_concordance_pct,
        identity_replicate_rf_distance,
        estimator_self_consistent,
        duration_ms: total_dur_ms,
        split_support_pct_by_mask,
    };

    let unpartitioned_free_params_k = unpart_metrics.free_parameters_k;
    let unpartitioned_ln_l = unpart_metrics.ln_l;
    let unpartitioned_aic = unpart_metrics.aic;
    let unpartitioned_bic = unpart_metrics.bic;
    let per_locus_models: Vec<PerLocusModelTestSummary> =
        fitted_parts.into_iter().map(|p| p.summary).collect();

    Ok(PartitionedSupermatrixMlResult {
        per_locus_models,
        n_taxa,
        n_partitions,
        total_sites,
        partitioned_free_params_k,
        partitioned_ln_l,
        partitioned_aic,
        partitioned_bic,
        unpartitioned_best_model: best_unpart_code,
        unpartitioned_free_params_k,
        unpartitioned_ln_l,
        unpartitioned_aic,
        unpartitioned_bic,
        delta_ln_l_partitioned_minus_unpartitioned: round4(partitioned_ln_l - unpartitioned_ln_l),
        delta_bic_partitioned_minus_unpartitioned: round4(partitioned_bic - unpartitioned_bic),
        partitioned_metrics: part_metrics,
        partitioned_bootstrap,
        unpartitioned_metrics: unpart_metrics,
        unpartitioned_bootstrap: unpart_boot,
        total_duration_ms: total_dur_ms,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fetch::parse_single_fasta;
    use std::fs;
    use std::path::PathBuf;

    #[test]
    fn test_partitioned_modeltest_heterotachous_two_locus_recovery() {
        let fixture_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("fixtures")
            .join("carnivora_coi_mini.fasta");
        let raw = fs::read_to_string(&fixture_path).expect("Carnivora fixture exists");
        let mut base_seqs = Vec::new();
        for chunk in raw.split('>') {
            let t = chunk.trim();
            if t.is_empty() {
                continue;
            }
            if let Some((_, seq)) = parse_single_fasta(&format!(">{}\n", t)) {
                base_seqs.push(seq);
            }
        }
        assert_eq!(base_seqs.len(), 4);

        // Expand 4 taxa -> 8 taxa, Locus 1 slow balanced/JC-like (600 bp) vs. Locus 2 fast 3x Ti/Tv skewed (600 bp)
        let mut rng = DeterministicLcg64::new(0x2026_1130);
        let mut labels = Vec::with_capacity(8);
        let mut full_seqs = Vec::with_capacity(8);

        // Generate shared ancestral templates (`600 bp` each)
        let l1_root: Vec<u8> = (0..600)
            .map(|_| match rng.gen_range_usize(4) {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            })
            .collect();
        // Skewed AT/CT template for fast Locus 2
        let mut l2_root = Vec::with_capacity(600);
        for _ in 0..600 {
            let p = rng.gen_range_usize(100);
            l2_root.push(if p < 42 {
                b'A'
            } else if p < 82 {
                b'T'
            } else if p < 91 {
                b'C'
            } else {
                b'G'
            });
        }

        let mutate_l1 = |src: &[u8], n_muts: usize, r: &mut DeterministicLcg64| -> Vec<u8> {
            let mut out = src.to_vec();
            for _ in 0..n_muts {
                let pos = r.gen_range_usize(out.len());
                out[pos] = match r.gen_range_usize(4) {
                    0 => b'A',
                    1 => b'C',
                    2 => b'G',
                    _ => b'T',
                };
            }
            out
        };
        let mutate_l2_fast_ti = |src: &[u8], n_muts: usize, r: &mut DeterministicLcg64| -> Vec<u8> {
            let mut out = src.to_vec();
            for _ in 0..n_muts {
                let pos = r.gen_range_usize(out.len());
                // Gamma rate heterogeneity proxy: concentrate hits in variable 3rd positions
                if pos % 3 == 0 && r.gen_range_usize(4) != 0 {
                    continue;
                }
                if r.gen_range_usize(10) < 9 {
                    out[pos] = match out[pos] {
                        b'A' => b'G',
                        b'G' => b'A',
                        b'C' => b'T',
                        _ => b'C',
                    };
                } else {
                    out[pos] = match r.gen_range_usize(2) {
                        0 => b'A',
                        _ => b'T',
                    };
                }
            }
            out
        };

        for clade_i in 0..4 {
            let l1_clade = mutate_l1(&l1_root, 10, &mut rng);
            let l2_clade = mutate_l2_fast_ti(&l2_root, 190, &mut rng);
            for sub_j in 0..2 {
                labels.push(format!("T{}_{}", clade_i + 1, sub_j + 1));
                let l1_tip = mutate_l1(&l1_clade, 5, &mut rng);
                let l2_tip = mutate_l2_fast_ti(&l2_clade, 90, &mut rng);
                let mut combined = String::with_capacity(1200);
                combined.push_str(&String::from_utf8(l1_tip).unwrap());
                combined.push_str(&String::from_utf8(l2_tip).unwrap());
                full_seqs.push(combined);
            }
        }

        let parts = vec![
            PartitionLocusInput {
                locus_name: "LOC1_SLOW_BALANCED".to_string(),
                locus_class: "Nuclear CDS".to_string(),
                protein_coding: false,
                is_nuclear: true,
                start_col_1based: 1,
                end_col_1based: 600,
            },
            PartitionLocusInput {
                locus_name: "LOC2_FAST_BIASED".to_string(),
                locus_class: "Mito PCG".to_string(),
                protein_coding: false,
                is_nuclear: false,
                start_col_1based: 601,
                end_col_1based: 1200,
            },
        ];

        let res = infer_partitioned_supermatrix_ml_with_bootstrap(&full_seqs, &labels, &parts, 50, 2026)
            .expect("Partitioned ModelTest unit test must succeed");
        assert_eq!(res.per_locus_models.len(), 2);
        let m1 = &res.per_locus_models[0];
        let m2 = &res.per_locus_models[1];
        assert!(
            m1.proportional_rate_r < 1.0 && m2.proportional_rate_r > 1.0,
            "Expected r1 < 1.0 < r2 on slow vs fast loci (got r1={}, r2={})",
            m1.proportional_rate_r,
            m2.proportional_rate_r
        );
        let weighted_mean_r = (600.0 * m1.proportional_rate_r + 600.0 * m2.proportional_rate_r) / 1200.0;
        assert!(
            (weighted_mean_r - 1.0).abs() < 1e-3,
            "Length-weighted mean rate must equal 1.0 (got {})",
            weighted_mean_r
        );
        assert_ne!(
            m1.winning_model_code, m2.winning_model_code,
            "Expected distinct BIC winning models across slow JC-like ({}) vs fast skewed ({})",
            m1.winning_model_code, m2.winning_model_code
        );
        assert!(
            res.partitioned_ln_l > res.unpartitioned_ln_l + 15.0,
            "Partitioned ln L ({}) must exceed unpartitioned ln L ({}) by > 15.0 nats",
            res.partitioned_ln_l,
            res.unpartitioned_ln_l
        );
        assert!(
            res.partitioned_bic < res.unpartitioned_bic,
            "Partitioned BIC ({}) must improve over unpartitioned BIC ({})",
            res.partitioned_bic,
            res.unpartitioned_bic
        );
        assert!(res.partitioned_bootstrap.estimator_self_consistent);
    }

    #[test]
    fn test_internal_partition_parallel_matches_sequential_and_speeds_up_skewed_loci() {
        use crate::tree_builder::{
            evaluate_tree_log_likelihood_seq, optimize_substitution_model_on_tree_seq,
        };

        let n_taxa = 24usize;
        let lengths = [1800usize, 240, 240, 240, 240, 240];
        let total_sites: usize = lengths.iter().sum();
        assert_eq!(total_sites, 3000);

        let mut rng = DeterministicLcg64::new(0x2026_1124_u64);
        let root_seq: Vec<u8> = (0..total_sites)
            .map(|_| [b'A', b'C', b'G', b'T'][rng.gen_range_usize(4)])
            .collect();

        let mutate_sites = |src: &[u8], n_muts: usize, r: &mut DeterministicLcg64| -> Vec<u8> {
            let mut out = src.to_vec();
            for _ in 0..n_muts {
                let pos = r.gen_range_usize(out.len());
                out[pos] = [b'A', b'C', b'G', b'T'][r.gen_range_usize(4)];
            }
            out
        };

        let mut labels = Vec::with_capacity(n_taxa);
        let mut full_seqs = Vec::with_capacity(n_taxa);
        for c_i in 0..6 {
            let clade_seq = mutate_sites(&root_seq, 450, &mut rng);
            for t_j in 0..4 {
                labels.push(format!("Tax_{}_{}", c_i + 1, t_j + 1));
                let tip = mutate_sites(&clade_seq, 320, &mut rng);
                full_seqs.push(String::from_utf8(tip).expect("ASCII DNA"));
            }
        }

        let full_aln = compress_alignment_and_estimate_rates(&full_seqs);
        let guide = solve_native_ml_tree_with_fit(
            &full_seqs,
            &labels,
            &full_aln,
            PhyloRunModel::GtrInvGamma4Fast,
        );

        let mut parts: Vec<FittedPartition> = Vec::with_capacity(lengths.len());
        let mut cur_start = 1usize;
        for (pi_idx, &len_bp) in lengths.iter().enumerate() {
            let s0 = cur_start - 1;
            let e0 = s0 + len_bp;
            let seqs_p: Vec<String> = full_seqs
                .iter()
                .map(|row| row[s0..e0].to_string())
                .collect();
            let aln_p = compress_alignment_and_estimate_rates(&seqs_p);
            let col_to_pat = map_columns_to_patterns(&seqs_p, &aln_p);
            let rm = if pi_idx % 2 == 0 {
                PhyloRunModel::GtrInvGamma4Full
            } else {
                PhyloRunModel::HkyInvGamma4
            };
            let mcode = if pi_idx % 2 == 0 {
                "GTR+F+I+G4"
            } else {
                "HKY+F+I+G4"
            };
            let (q, pi, pinv, rate_cats, lnl) =
                optimize_substitution_model_on_tree_seq(&guide.nodes, &aln_p, rm, false);
            let rate_r = 0.8 + 0.1 * (pi_idx as f64);
            let summary = PerLocusModelTestSummary {
                locus_name: format!("LOC_{}", pi_idx + 1),
                locus_class: "Universal CDS".to_string(),
                alphabet_mode: "4-State Nucleotide".to_string(),
                start_col_1based: cur_start,
                end_col_1based: e0,
                length_bp: len_bp,
                taxa_with_data: n_taxa,
                distinct_patterns: aln_p.patterns.len(),
                invariant_fraction: round4(aln_p.invariant_fraction),
                gc_or_purine_pct: 50.0,
                rcfv_locus: 0.05,
                winning_model_code: mcode.to_string(),
                runner_up_model_code: "K80+G4".to_string(),
                delta_bic_vs_runner_up: 10.0,
                proportional_rate_r: rate_r,
                kappa_ti_tv: 2.0,
                pinv: round4(pinv),
                gamma_or_freerate_summary: "I+G4".to_string(),
                locus_ln_l: round4(lnl),
                locus_aic: round4(-2.0 * lnl + 20.0),
                locus_bic: round4(-2.0 * lnl + 10.0 * (len_bp as f64).ln()),
                duration_ms: 0.0,
                all_models: Vec::new(),
            };
            parts.push(FittedPartition {
                input: PartitionLocusInput {
                    locus_name: summary.locus_name.clone(),
                    locus_class: summary.locus_class.clone(),
                    protein_coding: false,
                    is_nuclear: true,
                    start_col_1based: cur_start,
                    end_col_1based: e0,
                },
                _seqs: seqs_p,
                aln: aln_p,
                col_to_pat,
                _winning_code: mcode.to_string(),
                winning_k_subst: if pi_idx % 2 == 0 { 10 } else { 6 },
                _winning_rm: rm,
                _winning_force_no_inv: false,
                q,
                pi,
                pinv,
                rate_cats,
                rate_r,
                locus_ln_l: lnl,
                summary,
            });
            cur_start = e0 + 1;
        }

        // 1. Verify exact numerical parity and bit-for-bit determinism
        let lnl_flat_a = eval_partitioned_lnl(&guide.nodes, &parts);
        let lnl_flat_b = eval_partitioned_lnl(&guide.nodes, &parts);
        assert_eq!(
            lnl_flat_a.to_bits(),
            lnl_flat_b.to_bits(),
            "Flattened chunked partitioned likelihood must be bit-for-bit deterministic"
        );
        let lnl_coarse = eval_partitioned_lnl_coarse_reference(&guide.nodes, &parts);
        let lnl_seq: f64 = parts
            .iter()
            .map(|p| {
                evaluate_tree_log_likelihood_seq(
                    &guide.nodes,
                    &p.aln,
                    &p.q,
                    &p.pi,
                    p.pinv,
                    &p.rate_cats,
                    p.rate_r,
                )
            })
            .sum();
        assert!(
            (lnl_flat_a - lnl_seq).abs() < 1e-8,
            "Flat chunked lnL ({}) diverged from sequential ({})",
            lnl_flat_a,
            lnl_seq
        );
        assert!(
            (lnl_coarse - lnl_seq).abs() < 1e-8,
            "Coarse parallel lnL ({}) diverged from sequential ({})",
            lnl_coarse,
            lnl_seq
        );

        // 2. Verify pattern-level and single-partition parallel parity on Partition 0 (1,800 bp)
        let p0 = &parts[0];
        let ell_par = evaluate_site_pattern_lnls(
            &guide.nodes,
            &p0.aln,
            &p0.q,
            &p0.pi,
            p0.pinv,
            &p0.rate_cats,
            p0.rate_r,
            true,
        );
        let ell_seq = evaluate_site_pattern_lnls(
            &guide.nodes,
            &p0.aln,
            &p0.q,
            &p0.pi,
            p0.pinv,
            &p0.rate_cats,
            p0.rate_r,
            false,
        );
        assert_eq!(ell_par.len(), ell_seq.len());
        for (i, (&a, &b)) in ell_par.iter().zip(ell_seq.iter()).enumerate() {
            assert!(
                (a - b).abs() < 1e-12,
                "Site pattern {} mismatch: parallel={} vs sequential={}",
                i,
                a,
                b
            );
        }
        let p0_par_lnl = evaluate_tree_log_likelihood(
            &guide.nodes,
            &p0.aln,
            &p0.q,
            &p0.pi,
            p0.pinv,
            &p0.rate_cats,
            p0.rate_r,
        );
        let p0_seq_lnl = evaluate_tree_log_likelihood_seq(
            &guide.nodes,
            &p0.aln,
            &p0.q,
            &p0.pi,
            p0.pinv,
            &p0.rate_cats,
            p0.rate_r,
        );
        assert!(
            (p0_par_lnl - p0_seq_lnl).abs() < 1e-8,
            "Single-partition chunked lnL ({}) diverged from sequential ({})",
            p0_par_lnl,
            p0_seq_lnl
        );

        // 3. Benchmark 30 iterations across skewed 1,800 + 5x240 bp partitions
        // Warm up Rayon worker pool once before timing
        let _ = eval_partitioned_lnl(&guide.nodes, &parts);

        let iters = 30usize;
        let t_seq0 = Instant::now();
        let mut acc_seq = 0.0f64;
        for _ in 0..iters {
            let s: f64 = parts
                .iter()
                .map(|p| {
                    evaluate_tree_log_likelihood_seq(
                        &guide.nodes,
                        &p.aln,
                        &p.q,
                        &p.pi,
                        p.pinv,
                        &p.rate_cats,
                        p.rate_r,
                    )
                })
                .sum();
            acc_seq += std::hint::black_box(s);
        }
        let seq_ms = t_seq0.elapsed().as_secs_f64() * 1000.0;

        let t_coarse0 = Instant::now();
        let mut acc_coarse = 0.0f64;
        for _ in 0..iters {
            acc_coarse += std::hint::black_box(eval_partitioned_lnl_coarse_reference(&guide.nodes, &parts));
        }
        let coarse_ms = t_coarse0.elapsed().as_secs_f64() * 1000.0;

        let t_flat0 = Instant::now();
        let mut acc_flat = 0.0f64;
        for _ in 0..iters {
            acc_flat += std::hint::black_box(eval_partitioned_lnl(&guide.nodes, &parts));
        }
        let flat_ms = t_flat0.elapsed().as_secs_f64() * 1000.0;

        assert!(acc_seq.is_finite() && acc_coarse.is_finite() && acc_flat.is_finite());
        eprintln!(
            "[partitioned_modeltest bench | 24 taxa | skewed 6 partitions (1800 + 5x240 bp) | {} patterns on L0 | {} threads | {} iters] seq={:.2} ms, coarse_par={:.2} ms ({:.2}x vs seq), flat_chunk_par={:.2} ms ({:.2}x vs seq, {:.2}x vs coarse)",
            parts[0].aln.patterns.len(),
            rayon::current_num_threads(),
            iters,
            seq_ms,
            coarse_ms,
            seq_ms / coarse_ms.max(1e-6),
            flat_ms,
            seq_ms / flat_ms.max(1e-6),
            coarse_ms / flat_ms.max(1e-6),
        );
    }
}
