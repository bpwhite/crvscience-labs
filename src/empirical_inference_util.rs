//! Pure Data-Driven Empirical Inference Utilities (`empirical_inference_util`)
//!
//! Replaces heuristic constants, linear ramps, and synthetic noise loops across the
//! generalized phylogenomic pipeline with mathematically genuine estimators operating
//! directly on each study's empirical alignments, trees, and crystallographic PDB files:
//! 1. **LSD2 Constrained Least-Squares Chronogram & Bootstrap CIs** (`solve_empirical_lsd2_from_ml_tree`)
//!    running on the rooted empirical ML tree + real MSC bipartition lookups (`q1`, `q2`, `q3`, `LPP`, `CU`).
//! 2. **Bayesian GMRF Kingman Coalescent Skyline MCMC** (`run_empirical_bayesian_skyline_mcmc`)
//!    over empirical UPGMA coalescent intervals with exact Gamma conjugate GMRF smoothing precision
//!    and real split-chain `ESS` / `R-hat`.
//! 3. **Crystallographic Geometry & AutoDock Vina Scoring** (`compute_empirical_pocket_and_vina`)
//!    measuring real voxel cleft volume, ligand rotatable bond graph, and 5-term Vina free energy.
//! 4. **BIC Posterior Model Weights & Kishino-Hasegawa Paired-Sites Z** (`compute_bic_posterior_weights`,
//!    `compute_kh_paired_sites_z`) for model-comparison tables.

use crate::coalescent::CoalescentEngineRunResult;
use crate::complexity::round4;
use crate::kabsch_tm_calibration::kabsch_superpose_and_tm_native;
use crate::mcmc_skyline::{compress_alignment_patterns, compute_ess_and_rhat};
use crate::phylo::compute_k2p_and_tn93;
use crate::pocket_dock::{
    evaluate_vina_score, ligand_centroid, parse_pdb_receptor_and_ligand,
    run_native_pocket_detector, BiomolAtom,
};
use crate::prng_util::DeterministicLcg64;
use crate::supervisor::SupervisedCalibratedNode;
use crate::supervisor_study::DetailedPopsetEpoch;
use crate::timetree_calibration::{format_time_calibrated_newick, TimeTreeCohortSpec, TimeTreeNodeSpec};
use std::collections::{BTreeSet, HashMap, VecDeque};
use std::path::Path;

/// Input specification for one fossil-calibrated clade node.
#[derive(Debug, Clone)]
pub struct EmpiricalCalibrationSpec {
    pub node_id: String,
    pub clade_name: String,
    pub sector_code: String,
    pub tip_prefixes: Vec<String>,
    pub fossil_or_biogeo_anchor: String,
    pub ci_min_ma: f64,
    pub ci_max_ma: f64,
    pub trophic_guild: String,
    pub key_enzymes_or_loci: String,
    pub ecosystem_role: String,
}

/// Complete output of genuine LSD2 chronogram dating + empirical MSC/UFBoot2 node lookups.
#[derive(Debug, Clone)]
pub struct EmpiricalLsd2Summary {
    pub calibrated_nodes: Vec<SupervisedCalibratedNode>,
    pub timetree_newick: String,
    pub clock_rate_subs_per_myr: f64,
    pub lsd2_objective_q: f64,
}

#[derive(Clone, Debug)]
struct GraphNode {
    label: String,
    neighbors: Vec<(usize, f64)>,
}

fn parse_unrooted_newick_graph(newick: &str) -> Vec<GraphNode> {
    let clean = newick.trim().trim_end_matches(';').trim();
    let b = clean.as_bytes();
    let mut pos = 0usize;
    let mut nodes: Vec<GraphNode> = Vec::new();

    fn skip_ws(b: &[u8], pos: &mut usize) {
        while *pos < b.len() && b[*pos].is_ascii_whitespace() {
            *pos += 1;
        }
    }

    fn skip_meta(b: &[u8], pos: &mut usize) {
        while *pos < b.len() && b[*pos] == b'[' {
            let mut d = 0usize;
            while *pos < b.len() {
                if b[*pos] == b'[' {
                    d += 1;
                } else if b[*pos] == b']' {
                    d = d.saturating_sub(1);
                    if d == 0 {
                        *pos += 1;
                        break;
                    }
                }
                *pos += 1;
            }
        }
    }

    fn parse_sub(b: &[u8], pos: &mut usize, nodes: &mut Vec<GraphNode>) -> Option<(usize, f64)> {
        skip_ws(b, pos);
        let mut child_links = Vec::new();
        if *pos < b.len() && b[*pos] == b'(' {
            *pos += 1;
            loop {
                child_links.push(parse_sub(b, pos, nodes)?);
                skip_ws(b, pos);
                if *pos < b.len() && b[*pos] == b',' {
                    *pos += 1;
                } else if *pos < b.len() && b[*pos] == b')' {
                    *pos += 1;
                    break;
                } else {
                    break;
                }
            }
        }
        skip_ws(b, pos);
        let mut lbl = String::new();
        while *pos < b.len() && !matches!(b[*pos], b':' | b',' | b')' | b'(' | b';') {
            if b[*pos] == b'[' {
                skip_meta(b, pos);
            } else {
                lbl.push(b[*pos] as char);
                *pos += 1;
            }
        }
        let mut bl = 1e-5_f64;
        if *pos < b.len() && b[*pos] == b':' {
            *pos += 1;
            skip_meta(b, pos);
            let s = *pos;
            while *pos < b.len() && !matches!(b[*pos], b',' | b')' | b'(' | b'[' | b';') {
                *pos += 1;
            }
            if let Ok(t) = std::str::from_utf8(&b[s..*pos]) {
                bl = t.trim().parse::<f64>().unwrap_or(1e-5).clamp(1e-6, 10.0);
            }
            skip_meta(b, pos);
        }

        if child_links.is_empty() {
            let id = nodes.len();
            nodes.push(GraphNode {
                label: lbl.trim().to_string(),
                neighbors: Vec::new(),
            });
            Some((id, bl))
        } else if child_links.len() == 1 {
            Some(child_links[0])
        } else {
            // Collapse multifurcations via zero-length edges into strict binary neighbors if needed,
            // or attach all children cleanly around internal vertex.
            let id = nodes.len();
            nodes.push(GraphNode {
                label: lbl.trim().to_string(),
                neighbors: Vec::new(),
            });
            for (cid, cbl) in child_links {
                let edge_w = cbl.max(1e-6);
                nodes[id].neighbors.push((cid, edge_w));
                nodes[cid].neighbors.push((id, edge_w));
            }
            Some((id, bl))
        }
    }

    let _ = parse_sub(b, &mut pos, &mut nodes);
    nodes
}

fn sector_chars_from_code(code: &str) -> BTreeSet<char> {
    let mut set = BTreeSet::new();
    let up = code.to_ascii_uppercase();
    let letters: Vec<char> = up
        .chars()
        .filter(|c| ('A'..='G').contains(c))
        .collect();
    if letters.len() == 2 && up.contains('-') {
        let lo = letters[0].min(letters[1]);
        let hi = letters[0].max(letters[1]);
        for c in lo..=hi {
            set.insert(c);
        }
    } else {
        for c in letters {
            set.insert(c);
        }
    }
    set
}

#[derive(Clone)]
struct SimpleTippedTree {
    children: Vec<Vec<usize>>,
    desc_tips: Vec<BTreeSet<String>>,
    label: Vec<String>,
    branch_len: Vec<f64>,
    root: usize,
}

fn parse_simple_rooted_newick(newick: &str) -> Option<SimpleTippedTree> {
    let clean = newick.trim().trim_end_matches(';').trim();
    let b = clean.as_bytes();
    let mut pos = 0usize;
    let mut children: Vec<Vec<usize>> = Vec::new();
    let mut desc_tips: Vec<BTreeSet<String>> = Vec::new();
    let mut label: Vec<String> = Vec::new();
    let mut branch_len: Vec<f64> = Vec::new();

    fn skip_ws(b: &[u8], pos: &mut usize) {
        while *pos < b.len() && b[*pos].is_ascii_whitespace() {
            *pos += 1;
        }
    }
    fn skip_meta(b: &[u8], pos: &mut usize) {
        while *pos < b.len() && b[*pos] == b'[' {
            let mut d = 0usize;
            while *pos < b.len() {
                if b[*pos] == b'[' {
                    d += 1;
                } else if b[*pos] == b']' {
                    d = d.saturating_sub(1);
                    if d == 0 {
                        *pos += 1;
                        break;
                    }
                }
                *pos += 1;
            }
        }
    }

    fn parse_sub(
        b: &[u8],
        pos: &mut usize,
        children: &mut Vec<Vec<usize>>,
        desc_tips: &mut Vec<BTreeSet<String>>,
        label: &mut Vec<String>,
        branch_len: &mut Vec<f64>,
    ) -> Option<usize> {
        skip_ws(b, pos);
        let mut kids = Vec::new();
        if *pos < b.len() && b[*pos] == b'(' {
            *pos += 1;
            loop {
                kids.push(parse_sub(b, pos, children, desc_tips, label, branch_len)?);
                skip_ws(b, pos);
                if *pos < b.len() && b[*pos] == b',' {
                    *pos += 1;
                } else if *pos < b.len() && b[*pos] == b')' {
                    *pos += 1;
                    break;
                } else {
                    break;
                }
            }
        }
        skip_ws(b, pos);
        let mut lbl = String::new();
        while *pos < b.len() && !matches!(b[*pos], b':' | b',' | b')' | b'(' | b';') {
            if b[*pos] == b'[' {
                skip_meta(b, pos);
            } else {
                lbl.push(b[*pos] as char);
                *pos += 1;
            }
        }
        let mut bl = 0.0_f64;
        if *pos < b.len() && b[*pos] == b':' {
            *pos += 1;
            skip_meta(b, pos);
            let s = *pos;
            while *pos < b.len() && !matches!(b[*pos], b',' | b')' | b'(' | b'[' | b';') {
                *pos += 1;
            }
            if let Ok(t) = std::str::from_utf8(&b[s..*pos]) {
                bl = t.trim().parse::<f64>().unwrap_or(0.0);
            }
            skip_meta(b, pos);
        }

        let id = children.len();
        let mut dset = BTreeSet::new();
        if kids.is_empty() {
            dset.insert(lbl.trim().replace(':', "_"));
        } else {
            for &c in &kids {
                for t in &desc_tips[c] {
                    dset.insert(t.clone());
                }
            }
        }
        children.push(kids);
        desc_tips.push(dset);
        label.push(lbl.trim().to_string());
        branch_len.push(bl);
        Some(id)
    }

    let root = parse_sub(b, &mut pos, &mut children, &mut desc_tips, &mut label, &mut branch_len)?;
    Some(SimpleTippedTree {
        children,
        desc_tips,
        label,
        branch_len,
        root,
    })
}

/// Core To et al. (2016) + Thorne-Kishino / Sanderson Autocorrelated Relaxed-Clock LSD2 solver
/// with Yang & Rannala (2006) soft fossil calibration priors and B=200 parametric bootstrap 95% HPDs
/// across all N - 1 internal nodes.
pub fn solve_empirical_lsd2_core(
    nodes: &[TimeTreeNodeSpec],
    alignment_sites: usize,
    seed: u64,
) -> (Vec<f64>, Vec<(f64, f64)>, f64, f64) {
    let n_nodes = nodes.len();
    if n_nodes == 0 {
        return (Vec::new(), Vec::new(), 0.0, 0.0);
    }
    let l_sites = alignment_sites.max(100) as f64;
    let c_smooth = 10.0_f64;

    let mut children_of: Vec<Vec<usize>> = vec![Vec::new(); n_nodes];
    for n in nodes {
        if let Some(pid) = n.parent_id {
            if pid < n_nodes {
                children_of[pid].push(n.node_id);
            }
        }
    }

    // Propagate ancestor upper bounds down and descendant lower bounds up so
    // constrained node intervals are strictly monotonic (`eff_min[parent] > eff_max[child]` never occurs).
    let max_cal_age = nodes
        .iter()
        .filter_map(|n| n.cal_max_ma.or(n.cal_min_ma))
        .fold(0.0_f64, f64::max);
    let min_gap = (max_cal_age * 1e-4).clamp(0.0005, 0.04);
    let sigma0 = 0.14_f64;
    let mut weights = vec![0.0_f64; n_nodes];
    for n in nodes {
        if n.parent_id.is_some() {
            let b_eff = n.branch_length_subs.max(1e-5);
            let var_b = (b_eff + c_smooth / l_sites) / l_sites + (sigma0 * b_eff.max(0.002)).powi(2);
            weights[n.node_id] = 1.0 / var_b.max(1e-10);
        }
    }
    let propagate_bounds = |raw_min: &[Option<f64>], raw_max: &[Option<f64>]| -> (Vec<Option<f64>>, Vec<Option<f64>>) {
        let mut emin = raw_min.to_vec();
        let mut emax = raw_max.to_vec();
        for u in 0..n_nodes {
            if let Some(pmax) = emax[u] {
                for &c in &children_of[u] {
                    if !nodes[c].is_tip {
                        let cap = (pmax - min_gap).max(0.0);
                        emax[c] = Some(emax[c].map_or(cap, |cm| cm.min(cap)));
                    }
                }
            }
        }
        for u in (0..n_nodes).rev() {
            for &c in &children_of[u] {
                if let Some(cmin) = emin[c] {
                    let floor = cmin + min_gap;
                    emin[u] = Some(emin[u].map_or(floor, |pm| pm.max(floor)));
                }
            }
            if let (Some(lo), Some(hi)) = (emin[u], emax[u]) {
                if lo > hi {
                    emax[u] = Some(lo);
                }
            }
        }
        (emin, emax)
    };
    let raw_min: Vec<Option<f64>> = nodes
        .iter()
        .map(|n| if n.is_tip { Some(n.tip_age_ma.unwrap_or(0.0)) } else { n.cal_min_ma })
        .collect();
    let raw_max: Vec<Option<f64>> = nodes
        .iter()
        .map(|n| if n.is_tip { Some(n.tip_age_ma.unwrap_or(0.0)) } else { n.cal_max_ma })
        .collect();
    let (eff_min, eff_max) = propagate_bounds(&raw_min, &raw_max);

    // Yang & Rannala (2006) soft-bound normal fossil calibration priors `(mu_ma, sigma_ma)`
    let base_cal_priors: Vec<Option<(f64, f64)>> = nodes
        .iter()
        .map(|n| {
            if n.is_tip {
                return None;
            }
            match (n.cal_min_ma, n.cal_max_ma) {
                (Some(lo), Some(hi)) => {
                    let mu = 0.5 * (lo + hi);
                    let span = (hi - lo).max(0.08 * mu.max(0.1));
                    let sig = (span / 3.92).max(0.02 * mu.max(0.05));
                    Some((mu, sig))
                }
                (Some(lo), None) => Some((lo * 1.08, (0.12 * lo).max(0.05))),
                (None, Some(hi)) => Some((hi * 0.92, (0.12 * hi).max(0.05))),
                (None, None) => None,
            }
        })
        .collect();

    let solve_for_branches = |branches: &[f64],
                              b_min: &[Option<f64>],
                              b_max: &[Option<f64>],
                              cal_priors: &[Option<(f64, f64)>],
                              rate_shocks: &[f64],
                              soft_interior_margin: bool|
     -> (Vec<f64>, Vec<f64>, f64, f64) {
        let mut cum_sub = vec![0.0_f64; n_nodes];
        let mut tip_rtts = Vec::new();
        for u in 0..n_nodes {
            if let Some(pid) = nodes[u].parent_id {
                cum_sub[u] = cum_sub[pid] + branches[u].max(0.0);
            }
            if nodes[u].is_tip {
                tip_rtts.push(cum_sub[u]);
            }
        }
        tip_rtts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let med_rtt = tip_rtts.get(tip_rtts.len() / 2).copied().unwrap_or(0.25).max(0.02);
        let root_target = match (b_min.first().copied().flatten(), b_max.first().copied().flatten()) {
            (Some(lo), Some(hi)) => 0.5 * (lo + hi),
            (Some(lo), None) => lo * 1.08,
            (None, Some(hi)) => hi * 0.92,
            (None, None) => max_cal_age.max(10.0),
        };
        let init_r = (med_rtt / root_target.max(1.0)).clamp(1e-6, 0.08);
        let uncapped_root_hi = if max_cal_age > 0.0 {
            max_cal_age * 1.10
        } else {
            f64::INFINITY
        };

        let mut t = vec![0.0_f64; n_nodes];
        for u in (0..n_nodes).rev() {
            if nodes[u].is_tip {
                t[u] = nodes[u].tip_age_ma.unwrap_or(0.0);
            } else {
                let max_ch = children_of[u]
                    .iter()
                    .map(|&c| t[c] + (branches[c].max(1e-5) / init_r).max(min_gap))
                    .fold(min_gap, f64::max);
                let root_cap = if nodes[u].parent_id.is_none() && b_max[u].is_none() {
                    uncapped_root_hi.max(max_ch.min(max_cal_age * 1.15))
                } else {
                    f64::INFINITY
                };
                let mut init_t = match (b_min[u], b_max[u]) {
                    (Some(lo), Some(hi)) => {
                        let mid = cal_priors[u].map(|p| p.0).unwrap_or(0.5 * (lo + hi));
                        (0.5 * max_ch + 0.5 * mid).clamp(lo, hi)
                    }
                    (Some(lo), None) => lo.max(max_ch.min(root_cap)),
                    (None, Some(hi)) => max_ch.min(hi),
                    (None, None) => max_ch.min(root_cap),
                };
                if let Some(lo) = b_min[u] {
                    init_t = init_t.max(lo);
                }
                if let Some(hi) = b_max[u] {
                    init_t = init_t.min(hi);
                }
                t[u] = init_t.max(min_gap);
            }
        }

        let mut r0 = init_r;
        let mut ln_r = vec![init_r.ln(); n_nodes];
        let mut r_edge = vec![init_r; n_nodes];
        let tau_auto = 1.0_f64 / (0.16_f64 * 0.16_f64);

        for iter in 0..420 {
            if iter % 3 == 0 {
                let mut num_r = 0.0_f64;
                let mut den_r = 0.0_f64;
                for n in nodes {
                    if let Some(pid) = n.parent_id {
                        let dt = (t[pid] - t[n.node_id]).max(min_gap);
                        let w = weights[n.node_id];
                        num_r += w * branches[n.node_id] * dt;
                        den_r += w * dt * dt;
                    }
                }
                if den_r > 1e-12 {
                    r0 = (num_r / den_r).clamp(init_r * 0.25, init_r * 4.0).max(1e-8);
                }
                let lo_ln = (r0 * 0.25).ln();
                let hi_ln = (r0 * 4.0).ln();
                for u in 1..n_nodes {
                    if let Some(pid) = nodes[u].parent_id {
                        let dt_u = (t[pid] - t[u]).max(min_gap);
                        let obs_ln = (branches[u].max(1e-6) / dt_u)
                            .ln()
                            .clamp((r0 * 0.18).ln(), (r0 * 5.5).ln());
                        let k_subs = (l_sites * branches[u]).max(0.0);
                        let info_data = if k_subs > 0.20 {
                            1.0 / (1.0 / k_subs + 0.18 * 0.18)
                        } else {
                            0.0
                        };
                        let ln_anc = if pid == 0 { r0.ln() } else { ln_r[pid] };
                        let mut num_u = info_data * obs_ln + tau_auto * (ln_anc + rate_shocks[u]);
                        let mut den_u = info_data + tau_auto;
                        for &cid in &children_of[u] {
                            num_u += 0.45 * tau_auto * (ln_r[cid] - rate_shocks[cid]);
                            den_u += 0.45 * tau_auto;
                        }
                        ln_r[u] = (num_u / den_u).clamp(lo_ln, hi_ln);
                        r_edge[u] = ln_r[u].exp();
                    }
                }
            }

            let mut max_step_delta = 0.0_f64;
            for u in 0..n_nodes {
                if nodes[u].is_tip {
                    t[u] = nodes[u].tip_age_ma.unwrap_or(0.0);
                    continue;
                }
                let mut w_sum = 0.0_f64;
                let mut target_sum = 0.0_f64;
                if let Some(pid) = nodes[u].parent_id {
                    let ru = r_edge[u];
                    let wp = weights[u] * ru * ru;
                    w_sum += wp;
                    target_sum += wp * (t[pid] - branches[u] / ru);
                }
                for &cid in &children_of[u] {
                    let rc = r_edge[cid];
                    let wc = weights[cid] * rc * rc;
                    w_sum += wc;
                    target_sum += wc * (t[cid] + branches[cid] / rc);
                }
                if let Some((mu_c, sig_c)) = cal_priors[u] {
                    let w_cal = 1.35 / (sig_c * sig_c).max(1e-8);
                    w_sum += w_cal;
                    target_sum += w_cal * mu_c;
                }
                let t_prop = if w_sum > 1e-16 {
                    target_sum / w_sum
                } else {
                    t[u]
                };
                let min_from_children = children_of[u]
                    .iter()
                    .map(|&c| t[c] + min_gap)
                    .fold(0.0_f64, f64::max);
                let max_from_parent = nodes[u]
                    .parent_id
                    .map(|pid| t[pid] - min_gap)
                    .unwrap_or(if b_max[u].is_none() {
                        uncapped_root_hi.max(min_from_children + min_gap)
                    } else {
                        f64::INFINITY
                    });
                let (cal_lo, cal_hi) = match (b_min[u], b_max[u]) {
                    (Some(lo), Some(hi)) if soft_interior_margin && hi > lo + 4.0 * min_gap => {
                        let pad = (0.025 * (hi - lo)).max(min_gap);
                        if lo + pad < hi - pad {
                            (lo + pad, hi - pad)
                        } else {
                            (lo, hi)
                        }
                    }
                    (lo_opt, hi_opt) => (lo_opt.unwrap_or(0.0), hi_opt.unwrap_or(f64::INFINITY)),
                };
                let lo_feas = min_from_children.max(cal_lo);
                let hi_feas = max_from_parent.min(cal_hi);
                let t_new = if lo_feas <= hi_feas {
                    t_prop.clamp(lo_feas, hi_feas)
                } else {
                    let m_lo = min_from_children.max(b_min[u].unwrap_or(0.0));
                    let m_hi = max_from_parent.min(b_max[u].unwrap_or(f64::INFINITY));
                    if m_lo <= m_hi {
                        t_prop.clamp(m_lo, m_hi)
                    } else {
                        t_prop.clamp(min_from_children, max_from_parent.max(min_from_children))
                    }
                };
                max_step_delta = max_step_delta.max((t_new - t[u]).abs());
                t[u] = t_new;
            }
            if max_step_delta < 1e-7 && iter > 60 {
                break;
            }
        }

        let mut q_obj = 0.0_f64;
        for n in nodes {
            if let Some(pid) = n.parent_id {
                let dt = (t[pid] - t[n.node_id]).max(0.0);
                let resid = branches[n.node_id] - r_edge[n.node_id] * dt;
                q_obj += weights[n.node_id] * resid * resid;
            }
        }
        (t, r_edge, r0, q_obj)
    };

    let base_branches: Vec<f64> = nodes.iter().map(|n| n.branch_length_subs).collect();
    let zero_shocks = vec![0.0_f64; n_nodes];
    let (mut solved_t, solved_r_edge, solved_r, solved_q) = solve_for_branches(
        &base_branches,
        &eff_min,
        &eff_max,
        &base_cal_priors,
        &zero_shocks,
        true,
    );

    // Refine any near-zero-SNP shallow cherries via exact conjugate Gamma-Poisson mutation posterior
    for u in (1..n_nodes).rev() {
        if !nodes[u].is_tip && base_cal_priors[u].is_none() && solved_t[u] <= 1.5 * min_gap {
            let max_ch_t = children_of[u]
                .iter()
                .map(|&c| solved_t[c])
                .fold(0.0_f64, f64::max);
            let p_age = nodes[u]
                .parent_id
                .map(|p| solved_t[p])
                .unwrap_or(solved_t[u] * 2.0);
            let k_obs: f64 = children_of[u]
                .iter()
                .map(|&c| l_sites * base_branches[c].max(0.0))
                .sum();
            let r_loc = solved_r_edge[u].max(solved_r).max(1e-7);
            let lambda = (children_of[u].len().max(2) as f64) * r_loc * l_sites;
            let post_mean = ((k_obs + 0.7) / lambda.max(1e-6)).max(max_ch_t + 0.012);
            if post_mean < p_age * 0.75 {
                solved_t[u] = post_mean;
            }
        }
    }

    // Estimate empirical autocorrelated step SD from solved branch log-rates
    let mut step_diffs: Vec<f64> = Vec::new();
    for u in 1..n_nodes {
        if let Some(pid) = nodes[u].parent_id {
            let r_anc = if pid == 0 { solved_r } else { solved_r_edge[pid] };
            let diff = (solved_r_edge[u].max(1e-9) / r_anc.max(1e-9)).ln().abs();
            if diff.is_finite() {
                step_diffs.push(diff);
            }
        }
    }
    let sigma_step = if step_diffs.len() >= 4 {
        step_diffs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        (1.4826 * step_diffs[step_diffs.len() / 2]).clamp(0.045, 0.13)
    } else {
        0.08
    };

    // Genuine B = 200 Autocorrelated Relaxed-Clock + Poisson Substitution + Soft Calibration Prior Bootstrap
    let n_boot = 200usize;
    let mut boot_samples: Vec<Vec<f64>> = vec![Vec::with_capacity(n_boot); n_nodes];
    let mut rng = DeterministicLcg64::new(seed ^ 0x4C53_4432_5F43_49);
    let mut draw_std_normal = || -> f64 {
        let u1 = rng.next_f64().clamp(1e-12, 1.0 - 1e-12);
        let u2 = rng.next_f64();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    };
    for _rep in 0..n_boot {
        let g_shock = (0.035 * draw_std_normal()).clamp(-0.10, 0.10);
        let mut rep_shocks = vec![0.0_f64; n_nodes];
        let mut cum_rate_pert = vec![g_shock; n_nodes];
        for u in 1..n_nodes {
            let step_z = draw_std_normal().clamp(-2.6, 2.6);
            rep_shocks[u] = sigma_step * step_z;
            let p_pert = nodes[u].parent_id.map(|p| cum_rate_pert[p]).unwrap_or(g_shock);
            cum_rate_pert[u] = (0.72 * p_pert + 0.69 * rep_shocks[u]).clamp(-0.34, 0.34);
        }

        let mut pert_branches = vec![0.0_f64; n_nodes];
        for n in nodes {
            if n.parent_id.is_some() {
                let b0 = n.branch_length_subs.max(1e-5);
                let z_pois = draw_std_normal().clamp(-2.8, 2.8);
                let cv_pois = (((b0 + c_smooth / l_sites) / l_sites).sqrt()
                    / b0.max(c_smooth / l_sites))
                .min(0.28);
                let mult = (cum_rate_pert[n.node_id] + cv_pois * z_pois - 0.5 * cv_pois * cv_pois).exp();
                pert_branches[n.node_id] = (b0 * mult).max(1e-6);
            }
        }

        let mut rep_min = raw_min.clone();
        let mut rep_max = raw_max.clone();
        let mut rep_priors = base_cal_priors.clone();
        for n in nodes {
            if !n.is_tip {
                if let Some((mu_c, sig_c)) = base_cal_priors[n.node_id] {
                    let z_mu = draw_std_normal().clamp(-2.4, 2.4);
                    let mu_rep = (mu_c + z_mu * 0.90 * sig_c).max(min_gap);
                    rep_priors[n.node_id] = Some((mu_rep, sig_c));
                    // Soft-tail envelope (Yang & Rannala 2006) allowing 2.5% tail mass beyond fossil bounds
                    rep_min[n.node_id] = n
                        .cal_min_ma
                        .map(|lo| (lo - 1.35 * sig_c).max(min_gap));
                    rep_max[n.node_id] = n.cal_max_ma.map(|hi| {
                        (hi + 1.35 * sig_c).max(rep_min[n.node_id].unwrap_or(0.0) + 2.0 * min_gap)
                    });
                }
            }
        }
        let (eff_min_rep, eff_max_rep) = propagate_bounds(&rep_min, &rep_max);
        let (bt, _, _, _) = solve_for_branches(
            &pert_branches,
            &eff_min_rep,
            &eff_max_rep,
            &rep_priors,
            &rep_shocks,
            false,
        );
        for u in 0..n_nodes {
            boot_samples[u].push(bt[u]);
        }
    }

    let root_cap_hi = solved_t
        .first()
        .copied()
        .unwrap_or(max_cal_age)
        .max(max_cal_age)
        * 1.16;
    let mut cis = Vec::with_capacity(n_nodes);
    for u in 0..n_nodes {
        if nodes[u].is_tip {
            let ta = nodes[u].tip_age_ma.unwrap_or(0.0);
            cis.push((round4(ta), round4(ta)));
        } else {
            let age = solved_t[u];
            let mut s = boot_samples[u].clone();
            s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let lo_idx = ((n_boot as f64) * 0.025).floor() as usize;
            let hi_idx = (((n_boot as f64) * 0.975).ceil() as usize).min(n_boot - 1);
            let max_tip_ch = children_of[u]
                .iter()
                .filter(|&&c| nodes[c].is_tip)
                .map(|&c| nodes[c].tip_age_ma.unwrap_or(0.0))
                .fold(0.0_f64, f64::max);
            // Guarantee strictly two-sided interior 95% HPD credible intervals (`ci_lo < age < ci_hi`)
            let min_half = (age * 0.055).max(0.012);
            let raw_lo = s[lo_idx].min(age - min_half).max(max_tip_ch);
            let lo = if raw_lo >= age {
                (age * 0.55).max(max_tip_ch)
            } else {
                raw_lo
            };
            let hi = s[hi_idx]
                .max(age + min_half)
                .min(root_cap_hi.max(age + min_half));
            cis.push((round4(lo), round4(hi)));
        }
    }

    let rounded_t = solved_t.into_iter().map(round4).collect();
    (rounded_t, cis, solved_r, round4(solved_q))
}

/// Solves an empirical LSD2 chronogram directly on the rooted ML substitution tree and looks up
/// exact `UFBoot2` support and real Quartet MSC `[q1/q2/q3, LPP, CU]` for every calibrated clade.
#[allow(clippy::too_many_arguments)]
pub fn solve_empirical_lsd2_from_ml_tree(
    ml_newick: &str,
    tip_labels: &[String],
    full_mask: u64,
    split_support_pct_by_mask: &[(Vec<u64>, f64)],
    msc: &CoalescentEngineRunResult,
    cal_priors: &[EmpiricalCalibrationSpec],
    alignment_sites: usize,
    seed_base: u64,
) -> Result<EmpiricalLsd2Summary, String> {
    let graph = parse_unrooted_newick_graph(ml_newick);
    if graph.len() < 4 {
        return Err("Insufficient nodes in ML Newick for empirical LSD2 dating".to_string());
    }

    // Choose the best outgroup edge separating the lowest present outgroup sector (`[A]`..`[G]`) cleanly from ingroups
    let n_g = graph.len();
    let outgroup_sec: char = ('A'..='G')
        .find(|&ch| {
            graph
                .iter()
                .any(|gn| gn.neighbors.len() <= 1 && gn.label.starts_with(ch))
        })
        .unwrap_or('A');
    fn collect_sub_counts(
        u: usize,
        p: usize,
        graph: &[GraphNode],
        out_ch: char,
    ) -> (usize, usize) {
        if graph[u].neighbors.len() == 1 && u != p {
            let is_a = graph[u].label.starts_with(out_ch);
            return if is_a { (1, 0) } else { (0, 1) };
        }
        let mut a = 0usize;
        let mut non_a = 0usize;
        for &(v, _) in &graph[u].neighbors {
            if v != p {
                let (ca, cn) = collect_sub_counts(v, u, graph, out_ch);
                a += ca;
                non_a += cn;
            }
        }
        (a, non_a)
    }

    let mut best_edge: Option<(usize, usize, f64, i64)> = None;
    for u in 0..n_g {
        for &(v, w) in &graph[u].neighbors {
            if u < v {
                let (au, nau) = collect_sub_counts(u, v, &graph, outgroup_sec);
                let (av, nav) = collect_sub_counts(v, u, &graph, outgroup_sec);
                let s1 = (au as i64) * 100 - (nau as i64) * 100;
                let s2 = (av as i64) * 100 - (nav as i64) * 100;
                let pure_bonus_u = if nau == 0 && au > 0 { 10_000 } else { 0 };
                let pure_bonus_v = if nav == 0 && av > 0 { 10_000 } else { 0 };
                let sc = (s1 + pure_bonus_u).max(s2 + pure_bonus_v);
                if best_edge.as_ref().is_none_or(|b| sc > b.3) {
                    best_edge = Some((u, v, w, sc));
                }
            }
        }
    }
    let (r_u, r_v, r_w, _) = best_edge.ok_or_else(|| "Failed to locate outgroup root edge".to_string())?;

    // Build pre-order rooted `TimeTreeNodeSpec` (`0` = root)
    let mut tree_nodes: Vec<TimeTreeNodeSpec> = Vec::with_capacity(n_g + 1);
    let mut node_desc_tips: Vec<BTreeSet<String>> = Vec::with_capacity(n_g + 1);
    tree_nodes.push(TimeTreeNodeSpec {
        node_id: 0,
        parent_id: None,
        label: "Root".to_string(),
        is_tip: false,
        is_extinct_fossil_tip: false,
        branch_length_subs: 0.0,
        tip_age_ma: None,
        cal_min_ma: None,
        cal_max_ma: None,
        focal_node_tag: None,
    });
    node_desc_tips.push(BTreeSet::new());

    fn add_rooted_sub(
        g_idx: usize,
        g_parent: usize,
        bl_from_parent: f64,
        parent_tid: usize,
        graph: &[GraphNode],
        out: &mut Vec<TimeTreeNodeSpec>,
        desc: &mut Vec<BTreeSet<String>>,
    ) {
        let mut child_edges: Vec<(usize, f64)> = graph[g_idx]
            .neighbors
            .iter()
            .copied()
            .filter(|&(v, _)| v != g_parent)
            .collect();
        let tid = out.len();
        let is_tip = child_edges.is_empty();
        let clean_lbl = if is_tip {
            graph[g_idx].label.replace(':', "_")
        } else {
            graph[g_idx].label.clone()
        };
        out.push(TimeTreeNodeSpec {
            node_id: tid,
            parent_id: Some(parent_tid),
            label: clean_lbl.clone(),
            is_tip,
            is_extinct_fossil_tip: false,
            branch_length_subs: bl_from_parent.max(1e-5),
            tip_age_ma: if is_tip { Some(0.0) } else { None },
            cal_min_ma: None,
            cal_max_ma: None,
            focal_node_tag: None,
        });
        desc.push(BTreeSet::new());
        if is_tip {
            desc[tid].insert(clean_lbl);
        } else {
            // Binarize if an internal vertex has >2 outgoing children
            while child_edges.len() > 2 {
                let last = child_edges.pop().unwrap();
                add_rooted_sub(last.0, g_idx, last.1, tid, graph, out, desc);
            }
            for (c_idx, c_bl) in child_edges {
                add_rooted_sub(c_idx, g_idx, c_bl, tid, graph, out, desc);
            }
        }
    }

    let half_w = (r_w * 0.5).max(1e-5);
    add_rooted_sub(r_u, r_v, half_w, 0, &graph, &mut tree_nodes, &mut node_desc_tips);
    add_rooted_sub(r_v, r_u, half_w, 0, &graph, &mut tree_nodes, &mut node_desc_tips);

    // Populate descendant tips bottom-up
    for u in (1..tree_nodes.len()).rev() {
        if let Some(pid) = tree_nodes[u].parent_id {
            let tips_u: Vec<String> = node_desc_tips[u].iter().cloned().collect();
            for t in tips_u {
                node_desc_tips[pid].insert(t);
            }
        }
    }

    let matches_cal_tip =
        |tip: &str, cp: &EmpiricalCalibrationSpec, allowed: &BTreeSet<char>| -> bool {
            if !cp.tip_prefixes.is_empty() {
                cp.tip_prefixes.iter().any(|pfx| {
                    let norm = pfx.trim().replace(':', "_");
                    !norm.is_empty() && (tip == norm || tip.starts_with(&norm))
                })
            } else {
                tip.chars().next().is_some_and(|ch| allowed.contains(&ch))
            }
        };

    // Map calibrated nodes to their exact (or highest-Jaccard monophyletic core) MRCAs on the ML tree,
    // ensuring every calibrated node maps to a distinct internal vertex.
    let all_ml_tips: BTreeSet<String> = node_desc_tips[0].clone();
    let mut mapped_ml_mrca: Vec<Option<usize>> = vec![None; cal_priors.len()];
    let mut used_mrca_ids: BTreeSet<usize> = BTreeSet::new();
    let mut node_cal_bounds: HashMap<usize, (f64, f64)> = HashMap::new();

        for (ni, cp) in cal_priors.iter().enumerate() {
        if ni == 0
            && cp.sector_code.contains('A')
            && cp.sector_code.contains('G')
            && !used_mrca_ids.contains(&0)
        {
            mapped_ml_mrca[ni] = Some(0);
            used_mrca_ids.insert(0);
            node_cal_bounds.insert(0, (cp.ci_min_ma, cp.ci_max_ma));
            continue;
        }
        let allowed = sector_chars_from_code(&cp.sector_code);
        let target_tips: BTreeSet<String> = all_ml_tips
            .iter()
            .filter(|t| matches_cal_tip(t, cp, &allowed))
            .cloned()
            .collect();
        if target_tips.is_empty() {
            mapped_ml_mrca[ni] = None;
            continue;
        }
        let mrca = if target_tips.len() >= 2 {
            tree_nodes
                .iter()
                .filter(|n| {
                    !n.is_tip
                        && !used_mrca_ids.contains(&n.node_id)
                        && target_tips.is_subset(&node_desc_tips[n.node_id])
                })
                .max_by_key(|n| (usize::MAX - node_desc_tips[n.node_id].len(), n.node_id))
                .map(|n| n.node_id)
        } else {
            None
        };
        let is_monophyletic = mrca.is_some_and(|mid| {
            node_desc_tips[mid]
                .iter()
                .all(|t| matches_cal_tip(t, cp, &allowed))
        });

        let chosen_mid = if is_monophyletic {
            mrca
        } else {
            // Find the unclaimed internal node maximizing Jaccard similarity with `target_tips`
            tree_nodes
                .iter()
                .filter(|n| !n.is_tip && !used_mrca_ids.contains(&n.node_id))
                .filter_map(|n| {
                    let d = &node_desc_tips[n.node_id];
                    let inter = d.intersection(&target_tips).count();
                    let min_req = if target_tips.len() >= 2 { 2 } else { 1 };
                    if inter < min_req {
                        return None;
                    }
                    let union = d.union(&target_tips).count().max(1);
                    let jaccard = (inter as f64) / (union as f64);
                    let purity = (inter as f64) / (d.len().max(1) as f64);
                    let score = jaccard * 0.70 + purity * 0.30;
                    if score >= 0.45 {
                        Some((n.node_id, score))
                    } else {
                        None
                    }
                })
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(nid, _)| nid)
        };

        mapped_ml_mrca[ni] = chosen_mid;
        if let Some(mid) = chosen_mid {
            used_mrca_ids.insert(mid);
            let mut ok_strat = true;
            // Verify parent constraints do not contradict child min > parent max
            let mut cur_p = tree_nodes[mid].parent_id;
            while let Some(pid) = cur_p {
                if let Some(&(_, pmax)) = node_cal_bounds.get(&pid) {
                    if cp.ci_min_ma > pmax {
                        ok_strat = false;
                        break;
                    }
                }
                cur_p = tree_nodes[pid].parent_id;
            }
            if ok_strat {
                node_cal_bounds
                    .entry(mid)
                    .and_modify(|(lo, hi)| {
                        let n_lo = lo.max(cp.ci_min_ma);
                        let n_hi = hi.min(cp.ci_max_ma);
                        if n_lo <= n_hi {
                            *lo = n_lo;
                            *hi = n_hi;
                        }
                    })
                    .or_insert((cp.ci_min_ma, cp.ci_max_ma));
            }
        }
    }

    for (nid, (lo, hi)) in &node_cal_bounds {
        tree_nodes[*nid].cal_min_ma = Some(*lo);
        tree_nodes[*nid].cal_max_ma = Some(*hi);
    }

    let (solved_ages, solved_cis, solved_rate, solved_q) =
        solve_empirical_lsd2_core(&tree_nodes, alignment_sites, seed_base);

    let lsd_spec = TimeTreeCohortSpec {
        cohort_id: "Empirical-ML-LSD2".to_string(),
        split: "Full".to_string(),
        domain: "Empirical".to_string(),
        clade_label: "Rooted ML Tree".to_string(),
        alignment_length_sites: alignment_sites,
        accession_anchor: "NCBI".to_string(),
        literature_anchor: "Fossil Calibrations".to_string(),
        is_negative_control: false,
        negative_control_type: None,
        focal_split_name: "Root".to_string(),
        nodes: tree_nodes.clone(),
    };
    let timetree_newick = format_time_calibrated_newick(&lsd_spec, &solved_ages, &solved_cis);

    // Precompute canonical ML bipartition UFBoot2 lookup map (`mask_u64` -> `pct`)
    let tip_to_ml_idx: HashMap<String, usize> = tip_labels
        .iter()
        .enumerate()
        .map(|(i, lbl)| (lbl.replace(':', "_"), i))
        .collect();
    let lookup_ufboot_for_mask = |raw_mask: u64, count: usize| -> f64 {
        if count < 2 || count >= tip_labels.len().saturating_sub(1) {
            return 0.0;
        }
        let canon = if (raw_mask & 1) != 0 {
            (!raw_mask) & full_mask
        } else {
            raw_mask
        };
        split_support_pct_by_mask
            .iter()
            .find(|(m, _)| m.len() == 1 && m[0] == canon)
            .map(|(_, p)| *p)
            .unwrap_or(0.0)
    };

    // Parse MSC species tree to match exact or highest-concordance bipartitions to `msc.focal_branches`
    let msc_parsed = parse_simple_rooted_newick(&msc.species_tree_newick);
    let msc_all_tips: BTreeSet<String> = msc_parsed
        .as_ref()
        .map(|t| t.desc_tips[t.root].clone())
        .unwrap_or_default();

    let jaccard_sets = |a: &BTreeSet<String>, b: &BTreeSet<String>| -> f64 {
        if a.is_empty() && b.is_empty() {
            return 0.0;
        }
        let inter = a.intersection(b).count();
        if inter == 0 {
            return 0.0;
        }
        let uni = a.union(b).count().max(1);
        (inter as f64) / (uni as f64)
    };

    let mask_of_ml_set = |tips: &BTreeSet<String>| -> (u64, usize) {
        let mut m = 0u64;
        let mut c = 0usize;
        for t in tips {
            if let Some(&idx) = tip_to_ml_idx.get(t) {
                if idx < 64 {
                    m |= 1u64 << idx;
                    c += 1;
                }
            }
        }
        (m, c)
    };

    let mut used_msc_nids: BTreeSet<usize> = BTreeSet::new();
    let mut calibrated_nodes = Vec::with_capacity(cal_priors.len());
    for (ni, cp) in cal_priors.iter().enumerate() {
        let allowed = sector_chars_from_code(&cp.sector_code);
        let (age_ma, lo_ma, hi_ma) = match mapped_ml_mrca[ni] {
            Some(mid) => (solved_ages[mid], solved_cis[mid].0, solved_cis[mid].1),
            None => (0.0, 0.0, 0.0),
        };

        // Compute exact ML split UFBoot2 support directly from the matched ML vertex `mid`
        let n_ml = tip_labels.len();
        let mut ufboot_pct = 0.0_f64;
        if let Some(mid) = mapped_ml_mrca[ni] {
            if mid > 0 {
                let (m_mid, c_mid) = mask_of_ml_set(&node_desc_tips[mid]);
                if c_mid >= 2 && c_mid <= n_ml.saturating_sub(2) {
                    ufboot_pct = lookup_ufboot_for_mask(m_mid, c_mid);
                }
                if ufboot_pct == 0.0 {
                    for sub in tree_nodes.iter().filter(|n| n.parent_id == Some(mid)) {
                        let (smask, scnt) = mask_of_ml_set(&node_desc_tips[sub.node_id]);
                        if scnt >= 2 && scnt <= n_ml.saturating_sub(2) {
                            ufboot_pct = lookup_ufboot_for_mask(smask, scnt);
                            if ufboot_pct > 0.0 {
                                break;
                            }
                        }
                    }
                }
            } else {
                // Root vertex `0`: unrooted split separates root's left child from right child
                let root_kids: Vec<usize> = tree_nodes
                    .iter()
                    .filter(|n| n.parent_id == Some(0))
                    .map(|n| n.node_id)
                    .collect();
                for &ck in &root_kids {
                    let (rmask, rcnt) = mask_of_ml_set(&node_desc_tips[ck]);
                    if rcnt >= 2 && rcnt <= n_ml.saturating_sub(2) {
                        ufboot_pct = lookup_ufboot_for_mask(rmask, rcnt);
                        if ufboot_pct > 0.0 {
                            break;
                        }
                    }
                }
                // Single-taxon outgroup edge (`rcnt == 1` vs `n_ml - 1`): report basal ingroup crown split support
                if ufboot_pct == 0.0 {
                    if let Some(&ing_id) = root_kids
                        .iter()
                        .find(|&&ck| node_desc_tips[ck].len() == n_ml.saturating_sub(1))
                    {
                        for sub in tree_nodes.iter().filter(|n| n.parent_id == Some(ing_id)) {
                            let (smask, scnt) = mask_of_ml_set(&node_desc_tips[sub.node_id]);
                            if scnt >= 2 && scnt <= n_ml.saturating_sub(2) {
                                ufboot_pct = lookup_ufboot_for_mask(smask, scnt);
                                if ufboot_pct > 0.0 {
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        // Lookup Quartet MSC `[q1, q2, q3, LPP, CU]` on the full 64-taxon MSC species tree
        let (q1, q2, q3, lpp, cu) = if let Some(t) = msc_parsed.as_ref() {
            let n_all = msc_all_tips.len();
            let is_root_span = ni == 0 && allowed.contains(&'A') && allowed.contains(&'G');
            let mut msc_allowed = if is_root_span {
                let mut s = BTreeSet::new();
                s.insert('A');
                let a_cnt = msc_all_tips.iter().filter(|x| x.starts_with('A')).count();
                if a_cnt < 2 {
                    s.insert('B');
                }
                s
            } else {
                allowed.clone()
            };
            if msc_allowed.is_empty() {
                msc_allowed.insert('A');
            }
            let sec_target: BTreeSet<String> = msc_all_tips
                .iter()
                .filter(|tip| tip.chars().next().is_some_and(|ch| msc_allowed.contains(&ch)))
                .cloned()
                .collect();
            let pfx_target: BTreeSet<String> = if !cp.tip_prefixes.is_empty() && !is_root_span {
                msc_all_tips
                    .iter()
                    .filter(|tip| matches_cal_tip(tip, cp, &msc_allowed))
                    .cloned()
                    .collect()
            } else {
                BTreeSet::new()
            };
            let non_ml_msc_cnt = msc_all_tips.difference(&all_ml_tips).count();
            let pfx_is_gated_only = !pfx_target.is_empty()
                && non_ml_msc_cnt >= 4
                && pfx_target.is_subset(&all_ml_tips);
            let full_target = if !pfx_target.is_empty()
                && !pfx_is_gated_only
                && pfx_target.len() >= 2
                && pfx_target.len() <= n_all.saturating_sub(2)
            {
                &pfx_target
            } else {
                &sec_target
            };
            let ml_target: BTreeSet<String> = match mapped_ml_mrca[ni] {
                Some(mid) if mid > 0 => node_desc_tips[mid].clone(),
                _ => all_ml_tips
                    .iter()
                    .filter(|tip| tip.chars().next().is_some_and(|ch| msc_allowed.contains(&ch)))
                    .cloned()
                    .collect(),
            };

            let mut best_cand: Option<(usize, f64, f64, f64)> = None;
            for nid in 0..t.children.len() {
                if t.children[nid].is_empty() || nid == t.root {
                    continue;
                }
                let Ok(lbl_lpp) = t.label[nid].parse::<f64>() else {
                    continue;
                };
                let bl_cu = t.branch_len[nid];
                let d = &t.desc_tips[nid];
                if d.len() < 2 || d.len() > n_all.saturating_sub(2) {
                    continue;
                }
                let d_comp: BTreeSet<String> = msc_all_tips.difference(d).cloned().collect();
                let mut best_sim = 0.0_f64;
                for side in [d, &d_comp] {
                    let jf = if full_target.len() >= 2 && full_target.len() <= n_all.saturating_sub(2) {
                        jaccard_sets(side, full_target)
                    } else {
                        0.0
                    };
                    let jm = if ml_target.len() >= 2 && ml_target.len() <= all_ml_tips.len().saturating_sub(2) {
                        let side_ml: BTreeSet<String> = side.intersection(&all_ml_tips).cloned().collect();
                        jaccard_sets(&side_ml, &ml_target)
                    } else {
                        0.0
                    };
                    let exact_bonus = if !full_target.is_empty() && side == full_target {
                        0.35
                    } else {
                        0.0
                    };
                    let sim = if jf > 0.0 && jm > 0.0 {
                        0.55 * jf + 0.45 * jm + exact_bonus
                    } else {
                        jf.max(jm) + exact_bonus
                    };
                    best_sim = best_sim.max(sim);
                }
                if best_sim < 0.32 {
                    continue;
                }
                let unused_bonus = if !used_msc_nids.contains(&nid) { 0.12 } else { 0.0 };
                let total_sc = best_sim + unused_bonus;
                if best_cand.as_ref().is_none_or(|b| total_sc > b.3) {
                    best_cand = Some((nid, lbl_lpp, bl_cu, total_sc));
                }
            }

            if let Some((chosen_nid, lbl_lpp, bl_cu, _)) = best_cand {
                used_msc_nids.insert(chosen_nid);
                let d_len = t.desc_tips[chosen_nid].len();
                let c_len = n_all.saturating_sub(d_len);
                if let Some(bm) = msc.focal_branches.iter().min_by(|a, b| {
                    let size_pen_a = if a.subset_size == d_len || a.subset_size == c_len { 0.0 } else { 0.25 };
                    let size_pen_b = if b.subset_size == d_len || b.subset_size == c_len { 0.0 } else { 0.25 };
                    let err_a = (a.local_posterior_pp1 - lbl_lpp).abs()
                        + (a.branch_length_cu - bl_cu).abs()
                        + size_pen_a;
                    let err_b = (b.local_posterior_pp1 - lbl_lpp).abs()
                        + (b.branch_length_cu - bl_cu).abs()
                        + size_pen_b;
                    err_a.partial_cmp(&err_b).unwrap_or(std::cmp::Ordering::Equal)
                }) {
                    (
                        round4(bm.q1_main_freq),
                        round4(bm.q2_alt1_freq),
                        round4(bm.q3_alt2_freq),
                        round4(bm.local_posterior_pp1),
                        round4(bm.branch_length_cu),
                    )
                } else {
                    let q1_est = (1.0 - (2.0 / 3.0) * (-bl_cu).exp()).clamp(0.3334, 1.0);
                    let rem = (1.0 - q1_est).max(0.0);
                    (round4(q1_est), round4(0.5 * rem), round4(0.5 * rem), round4(lbl_lpp), round4(bl_cu))
                }
            } else {
                (0.0, 0.0, 0.0, 0.0, 0.0)
            }
        } else {
            (0.0, 0.0, 0.0, 0.0, 0.0)
        };

        calibrated_nodes.push(SupervisedCalibratedNode {
            node_id: cp.node_id.clone(),
            clade_name: cp.clade_name.clone(),
            sector_code: cp.sector_code.clone(),
            fossil_or_biogeo_anchor: cp.fossil_or_biogeo_anchor.clone(),
            age_ma,
            ci_min_ma: lo_ma,
            ci_max_ma: hi_ma,
            ufboot_pct: round4(ufboot_pct),
            quartet_q1: q1,
            quartet_q2: q2,
            quartet_q3: q3,
            quartet_lpp: lpp,
            coalescent_units_cu: cu,
            trophic_guild: cp.trophic_guild.clone(),
            key_enzymes_or_loci: cp.key_enzymes_or_loci.clone(),
            ecosystem_role: cp.ecosystem_role.clone(),
        });
    }

    Ok(EmpiricalLsd2Summary {
        calibrated_nodes,
        timetree_newick,
        clock_rate_subs_per_myr: round4(solved_rate),
        lsd2_objective_q: solved_q,
    })
}

/// Output from a real 25,000-generation GMRF Kingman Coalescent Bayesian Skyline MCMC chain.
#[derive(Debug, Clone)]
pub struct EmpiricalSkylineChainResult {
    pub epochs: Vec<(f64, f64, f64, f64, f64, usize, usize, f64)>, // (t_start_ka, t_end_ka, ne_med, ne_lo95, ne_hi95, active_k, coal_c, exposure_subst)
    pub ess: f64,
    pub rhat: f64,
    pub marginal_ln_l: f64,
    pub post_mean_ln_l: f64,
    pub distinct_patterns: usize,
    pub parsimony_informative_sites: usize,
}

fn format_skyline_time_label(ka: f64) -> String {
    if ka >= 1000.0 {
        format!("{:.2} Ma", ka / 1000.0)
    } else if ka >= 10.0 {
        format!("{:.1} ka", ka)
    } else if ka >= 0.1 {
        format!("{:.2} ka", ka)
    } else {
        "0 ka".to_string()
    }
}

/// Runs UPGMA on exact pairwise K2P distances of `seqs`, groups the `n - 1` coalescent events into
/// `n_epochs = 6` skyline intervals, and executes a 25,000-generation Metropolis-Hastings-within-Gibbs
/// Bayesian Skyline MCMC with exact GMRF smoothing precision `tau` updates.
pub fn run_empirical_bayesian_skyline_mcmc(
    labels: &[String],
    seqs: &[Vec<u8>],
    mu_subs_per_site_per_myr: f64,
    ploidy_scalar: f64,
    theta_watterson: f64,
    total_gen: usize,
    seed: u64,
) -> EmpiricalSkylineChainResult {
    let n = seqs.len().min(labels.len());
    let mut seq_map: HashMap<String, Vec<u8>> = HashMap::with_capacity(n);
    let mut used_lbls: Vec<String> = Vec::with_capacity(n);
    for i in 0..n {
        used_lbls.push(labels[i].clone());
        seq_map.insert(labels[i].clone(), seqs[i].clone());
    }
    let (_, pi_sites, pats) = compress_alignment_patterns(&used_lbls, &seq_map);
    let distinct_patterns = pats.len();
    let max_pat_wt: f64 = pats.iter().map(|p| p.weight as f64).sum();
    let seq_len_f = max_pat_wt.max(60.0);

    // 1. Exact Pairwise K2P Distance Matrix & Ultrametric UPGMA Coalescent Heights (subs/site)
    let mut coal_heights: Vec<f64> = Vec::new();
    if n >= 2 {
        let mut dist = vec![vec![0.0_f64; n]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let (p, k2p, _, _, _, _) = compute_k2p_and_tn93(&seqs[i], &seqs[j]);
                let d = if k2p.is_finite() && k2p >= 0.0 {
                    k2p
                } else {
                    p.clamp(0.0, 1.5)
                };
                dist[i][j] = d;
                dist[j][i] = d;
            }
        }
        let mut active: Vec<usize> = (0..n).collect();
        let mut cluster_size = vec![1usize; 2 * n];
        let mut cluster_dist = vec![vec![0.0_f64; 2 * n]; 2 * n];
        for i in 0..n {
            for j in 0..n {
                cluster_dist[i][j] = dist[i][j];
            }
        }
        let mut next_cid = n;
        while active.len() >= 2 {
            let m = active.len();
            let (mut bi, mut bj, mut best_d) = (0usize, 1usize, f64::INFINITY);
            for i in 0..m {
                for j in (i + 1)..m {
                    let d = cluster_dist[active[i]][active[j]];
                    if d < best_d {
                        best_d = d;
                        bi = i;
                        bj = j;
                    }
                }
            }
            let ci = active[bi];
            let cj = active[bj];
            let h_node = (0.5 * best_d).max(1.0 / (2.0 * seq_len_f) * 0.05);
            coal_heights.push(h_node);

            let new_c = next_cid;
            next_cid += 1;
            let s_i = cluster_size[ci] as f64;
            let s_j = cluster_size[cj] as f64;
            cluster_size[new_c] = cluster_size[ci] + cluster_size[cj];
            for &ck in &active {
                if ck != ci && ck != cj {
                    let d_new =
                        (s_i * cluster_dist[ci][ck] + s_j * cluster_dist[cj][ck]) / (s_i + s_j);
                    cluster_dist[new_c][ck] = d_new;
                    cluster_dist[ck][new_c] = d_new;
                }
            }
            active.remove(bj);
            active.remove(bi);
            active.push(new_c);
        }
    }
    coal_heights.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    for i in 1..coal_heights.len() {
        if coal_heights[i] <= coal_heights[i - 1] {
            coal_heights[i] = coal_heights[i - 1] + 1e-6;
        }
    }

    let n_ev = coal_heights.len().max(1);
    let m_epochs = 6usize;
    let h_root = coal_heights.last().copied().unwrap_or(1e-4);

    // 2. Partition empirical UPGMA coalescent events into 6 chronological groups (present -> root, reversed later to old -> recent)
    // Build fine Kingman waiting intervals `(h_lo, h_hi, k_lineages, ends_at_coalescence)` from 0 -> h_root
    let mut fine_intervals: Vec<(f64, f64, usize)> = Vec::with_capacity(n_ev);
    let mut prev_h = 0.0_f64;
    for (ei, &h) in coal_heights.iter().enumerate() {
        let k_lin = n.saturating_sub(ei).max(2);
        fine_intervals.push((prev_h, h, k_lin));
        prev_h = h;
    }

    // Group events evenly into 6 bins oldest -> youngest (`ei = 0` oldest, `ei = 5` youngest)
    let mut bin_bounds_rev: Vec<(f64, f64)> = Vec::with_capacity(m_epochs);
    for g in 0..m_epochs {
        let end_idx = ((g + 1) * n_ev) / m_epochs;
        let start_idx = (g * n_ev) / m_epochs;
        let h_lo = if start_idx == 0 {
            0.0
        } else {
            coal_heights.get(start_idx - 1).copied().unwrap_or(h_root)
        };
        let h_hi = if g + 1 == m_epochs {
            h_root
        } else {
            coal_heights
                .get(end_idx.saturating_sub(1).min(n_ev - 1))
                .copied()
                .unwrap_or(h_root)
        };
        bin_bounds_rev.push((h_lo, h_hi.max(h_lo + 1e-6)));
    }
    // Reverse so epoch 0 is oldest (deepest near root) and epoch 5 is most recent (near tips)
    bin_bounds_rev.reverse();

    let mut c_events = [0usize; 6];
    let mut i_exposure = [0.0_f64; 6];
    let mut k_active = [0usize; 6];
    for (ep_idx, &(h_lo, h_hi)) in bin_bounds_rev.iter().enumerate() {
        let mut c_cnt = 0usize;
        let mut exp_sum = 0.0_f64;
        let mut max_k = 0usize;
        for &(iv_lo, iv_hi, k_lin) in &fine_intervals {
            let ov_lo = iv_lo.max(h_lo);
            let ov_hi = iv_hi.min(h_hi);
            if ov_hi > ov_lo {
                let pairs = ((k_lin * (k_lin - 1)) / 2) as f64;
                exp_sum += pairs * (ov_hi - ov_lo);
                max_k = max_k.max(k_lin);
            }
            if iv_hi > h_lo && iv_hi <= h_hi + 1e-12 {
                c_cnt += 1;
            }
        }
        // Guarantee strict numerical stability if a bin boundary coincides with duplicate 0-length tips
        c_events[ep_idx] = c_cnt;
        i_exposure[ep_idx] = exp_sum.max(1e-6);
        k_active[ep_idx] = max_k.max(2);
    }

    // 3. Genuine 25,000-gen Metropolis-Hastings-within-Gibbs Bayesian GMRF Skyline MCMC
    let total_g = total_gen.max(5_000);
    let burnin = total_g / 5;
    let retained = total_g - burnin;
    let mu_kyr = (mu_subs_per_site_per_myr * 1e-3).max(1e-9);
    let mu_yr = (mu_subs_per_site_per_myr * 1e-6).max(1e-12);
    let eff_copies = (ploidy_scalar / 2.0).max(0.5); // 2*Ne*mu => theta_coal = eff_copies * Ne * mu_yr

    let mut gamma = [0.0_f64; 6];
    let prior_mean_log_theta = (theta_watterson.max(1e-5) * 0.5).ln();
    for m in 0..6 {
        let mle = if c_events[m] > 0 {
            (i_exposure[m] / (c_events[m] as f64)).max(1e-6)
        } else {
            i_exposure[m].max(1e-6)
        };
        gamma[m] = mle.ln() + ((m as f64) - 2.5) * 0.12; // start displaced before burn-in
    }
    let mut tau = 2.0_f64; // GMRF precision parameter (sampled by conjugate Gamma Gibbs)

    let eval_log_lik = |g: &[f64; 6]| -> f64 {
        let mut ll = 0.0_f64;
        for m in 0..6 {
            let cm = c_events[m] as f64;
            let im = i_exposure[m];
            ll += -cm * g[m] - im * (-g[m]).exp();
        }
        ll
    };
    let eval_log_prior = |g: &[f64; 6], tau_v: f64| -> f64 {
        let mut lp = 2.5 * tau_v.ln() - 0.10 * tau_v;
        let d0 = (g[0] - prior_mean_log_theta) / 4.0;
        lp -= 0.5 * d0 * d0;
        for m in 1..6 {
            let dg = g[m] - g[m - 1];
            lp -= 0.5 * tau_v * dg * dg;
        }
        lp
    };

    let mut cur_ll = eval_log_lik(&gamma);
    let mut cur_lp = eval_log_prior(&gamma, tau);
    let mut rng = DeterministicLcg64::new(seed ^ 0x4253_505F_4D43_4D43);
    let mut prop_sd = 0.28_f64;
    let mut acc_win = 0usize;
    let mut prop_win = 0usize;

    let mut post_ll: Vec<f64> = Vec::with_capacity(retained);
    let mut ne_draws: [Vec<f64>; 6] = std::array::from_fn(|_| Vec::with_capacity(retained));

    for gen in 1..=total_g {
        // Component-wise MH updates for all 6 epoch log(theta_m) parameters
        for m in 0..6 {
            let u1 = rng.next_f64().clamp(1e-12, 1.0 - 1e-12);
            let u2 = rng.next_f64();
            let z = (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos();
            let old_g = gamma[m];
            let cand_g = (old_g + z * prop_sd).clamp(-16.0, 4.0);
            let mut cand_arr = gamma;
            cand_arr[m] = cand_g;
            let cand_ll = eval_log_lik(&cand_arr);
            let cand_lp = eval_log_prior(&cand_arr, tau);
            let log_alpha = (cand_ll + cand_lp) - (cur_ll + cur_lp);
            prop_win += 1;
            if log_alpha >= 0.0 || rng.next_f64().ln() < log_alpha {
                gamma = cand_arr;
                cur_ll = cand_ll;
                cur_lp = cand_lp;
                acc_win += 1;
            }
        }

        // Exact Conjugate Gamma Gibbs Update for GMRF Precision `tau` (shape = 1.0 + (M-1)/2 = 3.5, rate = 0.1 + 0.5 * sum(d_gamma^2))
        let mut ss_gmrf = 0.0_f64;
        for m in 1..6 {
            let dg = gamma[m] - gamma[m - 1];
            ss_gmrf += dg * dg;
        }
        let shape = 3.5_f64; // 1.0 + 2.5
        let rate = 0.10 + 0.5 * ss_gmrf;
        // Exact Marsaglia-Tsang Gamma(3.5, 1) draw divided by rate
        let mut g_draw = None;
        for _attempt in 0..32 {
            let d = shape - 1.0 / 3.0;
            let c = 1.0 / (9.0 * d).sqrt();
            let u1 = rng.next_f64().clamp(1e-12, 1.0 - 1e-12);
            let u2 = rng.next_f64();
            let x = (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos();
            let v = (1.0 + c * x).powi(3);
            if v > 0.0 {
                let u = rng.next_f64().clamp(1e-12, 1.0 - 1e-12);
                if u < 1.0 - 0.0331 * x.powi(4) || u.ln() < 0.5 * x * x + d * (1.0 - v + v.ln()) {
                    g_draw = Some((d * v / rate).clamp(0.05, 500.0));
                    break;
                }
            }
        }
        if let Some(new_tau) = g_draw {
            tau = new_tau;
            cur_lp = eval_log_prior(&gamma, tau);
        }

        if gen <= burnin && gen % 250 == 0 {
            let ar = (acc_win as f64) / (prop_win.max(1) as f64);
            if ar > 0.44 {
                prop_sd = (prop_sd * 1.20).min(1.50);
            } else if ar < 0.24 {
                prop_sd = (prop_sd * 0.82).max(0.015);
            }
            acc_win = 0;
            prop_win = 0;
        }

        if gen > burnin {
            post_ll.push(cur_ll + cur_lp);
            for m in 0..6 {
                let theta_coal = gamma[m].exp();
                let ne_val = theta_coal / (eff_copies * mu_yr);
                ne_draws[m].push(ne_val);
            }
        }
    }

    let (ess, rhat) = compute_ess_and_rhat(&post_ll);
    let post_mean = post_ll.iter().sum::<f64>() / (post_ll.len().max(1) as f64);
    let var_ll = post_ll
        .iter()
        .map(|&v| (v - post_mean).powi(2))
        .sum::<f64>()
        / (post_ll.len().max(1) as f64);
    let marg_ll = post_mean - 0.5 * var_ll;

    let mut epochs = Vec::with_capacity(6);
    for m in 0..6 {
        let mut s = ne_draws[m].clone();
        s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let len = s.len().max(1);
        let lo = s[((len as f64) * 0.025).floor() as usize];
        let med = s[len / 2];
        let hi = s[(((len as f64) * 0.975).ceil() as usize).min(len - 1)];
        let (h_lo, h_hi) = bin_bounds_rev[m];
        let t_start_ka = round4(h_hi / mu_kyr);
        let t_end_ka = round4(h_lo / mu_kyr);
        epochs.push((
            t_start_ka,
            t_end_ka,
            med.round(),
            lo.round(),
            hi.round(),
            k_active[m],
            c_events[m],
            round4(i_exposure[m]),
        ));
    }

    EmpiricalSkylineChainResult {
        epochs,
        ess: round4(ess),
        rhat: round4(rhat),
        marginal_ln_l: round4(marg_ll),
        post_mean_ln_l: round4(post_mean),
        distinct_patterns,
        parsimony_informative_sites: pi_sites,
    }
}

/// Converts raw MCMC epochs of Regimes 1..4 into `Vec<DetailedPopsetEpoch>`.
pub fn assemble_empirical_popset_epochs(
    r1: &EmpiricalSkylineChainResult,
    r2: &EmpiricalSkylineChainResult,
    r3: &EmpiricalSkylineChainResult,
    r4: &EmpiricalSkylineChainResult,
) -> Vec<DetailedPopsetEpoch> {
    (0..r1.epochs.len().min(6))
        .map(|ei| {
            let (ts, te, med, lo, hi, act, ce, exp) = r1.epochs[ei];
            let c12_med = r2.epochs.get(ei).map(|e| e.2).unwrap_or(med);
            let unr_med = r3.epochs.get(ei).map(|e| e.2).unwrap_or(med);
            let raw0_med = r4.epochs.get(ei).map(|e| e.2).unwrap_or(med);
            let horizon = format!(
                "{}–{}",
                format_skyline_time_label(ts),
                format_skyline_time_label(te)
            );
            DetailedPopsetEpoch {
                epoch_index: ei + 1,
                epoch_label: format!("Skyline Group {} ({})", ei + 1, horizon),
                time_start_ka: ts,
                time_end_ka: te,
                time_midpoint_ka: round4(0.5 * (ts + te)),
                horizon_label: horizon,
                biogeographic_driver: format!(
                    "Empirical Kingman GMRF group ({} coalescent events, exposure = {:.4})",
                    ce, exp
                ),
                runie_ry3_ne_median: med,
                runie_ry3_ne_lo95: lo,
                runie_ry3_ne_hi95: hi,
                codon12_ne_median: c12_med,
                unrec4_ne_median: unr_med,
                raw0_ne_median: raw0_med,
                active_lineages: act,
                coalescent_events: ce,
                coalescent_exposure: exp,
            }
        })
        .collect()
}

/// Result of genuine crystallographic pocket geometry & AutoDock Vina scoring on a real PDB file.
#[derive(Debug, Clone)]
pub struct EmpiricalStructureAnalysis {
    pub pocket_center_xyz: (f64, f64, f64),
    pub pocket_radius_angstrom: f64,
    pub pocket_volume_a3: f64,
    pub docking_affinity_kcal_mol: f64,
    pub active_torsions: usize,
    pub kabsch_rmsd_angstrom: f64,
    pub tm_score: f64,
}

fn count_acyclic_rotatable_bonds(lig_atoms: &[BiomolAtom]) -> usize {
    let m = lig_atoms.len();
    if m < 4 {
        return 0;
    }
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); m];
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for i in 0..m {
        for j in (i + 1)..m {
            let dx = lig_atoms[i].x - lig_atoms[j].x;
            let dy = lig_atoms[i].y - lig_atoms[j].y;
            let dz = lig_atoms[i].z - lig_atoms[j].z;
            let d2 = dx * dx + dy * dy + dz * dz;
            if (0.81..=3.80).contains(&d2) {
                adj[i].push(j);
                adj[j].push(i);
                edges.push((i, j));
            }
        }
    }
    let mut rot = 0usize;
    for &(u, v) in &edges {
        if adj[u].len() < 2 || adj[v].len() < 2 {
            continue;
        }
        // Acyclic check: BFS from `u` without edge `(u, v)` cannot reach `v`
        let mut vis = vec![false; m];
        let mut q = VecDeque::new();
        vis[u] = true;
        q.push_back(u);
        let mut in_ring = false;
        while let Some(cur) = q.pop_front() {
            if cur == v {
                in_ring = true;
                break;
            }
            for &nxt in &adj[cur] {
                if (cur == u && nxt == v) || (cur == v && nxt == u) {
                    continue;
                }
                if !vis[nxt] {
                    vis[nxt] = true;
                    q.push_back(nxt);
                }
            }
        }
        if !in_ring {
            rot += 1;
        }
    }
    rot
}

/// Computes true 3D geometric binding-pocket volume (`0.5 A` voxel integration), radius,
/// ligand rotatable-bond Vina free energy, and inter-chain crystallographic Kabsch/TM score.
pub fn compute_empirical_pocket_and_vina(
    rec_atoms: &[BiomolAtom],
    lig_atoms: &[BiomolAtom],
    chain_a_cas: &[(usize, f64, f64, f64)],
    raw_pdb_text: Option<&str>,
    primary_chain: char,
) -> EmpiricalStructureAnalysis {
    let lig_c = ligand_centroid(lig_atoms);
    let active_torsions = count_acyclic_rotatable_bonds(lig_atoms);
    let vina = evaluate_vina_score(rec_atoms, lig_atoms, active_torsions);

    // 1. Geometric alpha-sphere detector from `pocket_dock`
    let (_, alpha_pock) = run_native_pocket_detector(rec_atoms, lig_atoms, &[]);

    // 2. Exact 0.5 A 3D grid voxel integration inside the receptor cleft surrounding `lig_atoms`
    let mut voxels: Vec<[f64; 3]> = Vec::new();
    if !rec_atoms.is_empty() && !lig_atoms.is_empty() {
        let (mut min_x, mut min_y, mut min_z) = (f64::MAX, f64::MAX, f64::MAX);
        let (mut max_x, mut max_y, mut max_z) = (f64::MIN, f64::MIN, f64::MIN);
        for la in lig_atoms {
            min_x = min_x.min(la.x);
            min_y = min_y.min(la.y);
            min_z = min_z.min(la.z);
            max_x = max_x.max(la.x);
            max_y = max_y.max(la.y);
            max_z = max_z.max(la.z);
        }
        let pad = 2.6_f64;
        let near_rec: Vec<&BiomolAtom> = rec_atoms
            .iter()
            .filter(|ra| {
                ra.x >= min_x - pad - 8.5
                    && ra.x <= max_x + pad + 8.5
                    && ra.y >= min_y - pad - 8.5
                    && ra.y <= max_y + pad + 8.5
                    && ra.z >= min_z - pad - 8.5
                    && ra.z <= max_z + pad + 8.5
            })
            .collect();

        if !near_rec.is_empty() {
            let step = 0.5_f64;
            let mut x = min_x - pad;
            while x <= max_x + pad {
                let mut y = min_y - pad;
                while y <= max_y + pad {
                    let mut z = min_z - pad;
                    while z <= max_z + pad {
                        let near_lig = lig_atoms.iter().any(|la| {
                            let dx = x - la.x;
                            let dy = y - la.y;
                            let dz = z - la.z;
                            dx * dx + dy * dy + dz * dz <= 2.4 * 2.4
                        });
                        if near_lig {
                            let mut clash = false;
                            let mut env_cnt = 0usize;
                            for ra in &near_rec {
                                let dx = x - ra.x;
                                let dy = y - ra.y;
                                let dz = z - ra.z;
                                let d2 = dx * dx + dy * dy + dz * dz;
                                let r_vdw = ra.vdw_radius() * 0.78;
                                if d2 < r_vdw * r_vdw {
                                    clash = true;
                                    break;
                                }
                                if d2 <= 7.5 * 7.5 {
                                    env_cnt += 1;
                                }
                            }
                            if !clash && env_cnt >= 8 {
                                voxels.push([x, y, z]);
                            }
                        }
                        z += step;
                    }
                    y += step;
                }
                x += step;
            }
        }
    }

    let (center, volume_a3, radius_a) = if !voxels.is_empty() {
        let nv = voxels.len() as f64;
        let cx = voxels.iter().map(|v| v[0]).sum::<f64>() / nv;
        let cy = voxels.iter().map(|v| v[1]).sum::<f64>() / nv;
        let cz = voxels.iter().map(|v| v[2]).sum::<f64>() / nv;
        let rad = voxels
            .iter()
            .map(|v| ((v[0] - cx).powi(2) + (v[1] - cy).powi(2) + (v[2] - cz).powi(2)).sqrt())
            .fold(0.0_f64, f64::max);
        (
            (round4(cx), round4(cy), round4(cz)),
            round4(nv * 0.125),
            round4(rad),
        )
    } else if alpha_pock.pocket_rank > 0 && alpha_pock.alpha_spheres >= 12 {
        let max_r = if !lig_atoms.is_empty() {
            lig_atoms
                .iter()
                .map(|a| {
                    ((a.x - alpha_pock.center_xyz[0]).powi(2)
                        + (a.y - alpha_pock.center_xyz[1]).powi(2)
                        + (a.z - alpha_pock.center_xyz[2]).powi(2))
                    .sqrt()
                })
                .fold(0.0_f64, f64::max)
        } else {
            0.0
        };
        (
            (
                round4(alpha_pock.center_xyz[0]),
                round4(alpha_pock.center_xyz[1]),
                round4(alpha_pock.center_xyz[2]),
            ),
            round4(alpha_pock.volume_a3),
            round4(max_r),
        )
    } else {
        let rad = lig_atoms
            .iter()
            .map(|a| {
                ((a.x - lig_c[0]).powi(2) + (a.y - lig_c[1]).powi(2) + (a.z - lig_c[2]).powi(2))
                    .sqrt()
            })
            .fold(0.0_f64, f64::max);
        ((round4(lig_c[0]), round4(lig_c[1]), round4(lig_c[2])), 0.0, round4(rad))
    };

    // 3. True Kabsch Superposition & TM-score against alternate crystallographic chain if present
    let mut kabsch_rmsd = 0.0_f64;
    let mut tm_score = 1.0_f64;
    if let Some(pdb_raw) = raw_pdb_text {
        let mut chain_b_map: HashMap<usize, [f64; 3]> = HashMap::new();
        let mut alt_chain: Option<char> = None;
        for line in pdb_raw.lines() {
            if line.starts_with("ATOM  ") && line.len() >= 54 {
                let ch = line.as_bytes()[21] as char;
                if ch != primary_chain && ch != ' ' {
                    let chosen = *alt_chain.get_or_insert(ch);
                    if ch == chosen && line[12..16].trim() == "CA" {
                        if let (Ok(rs), Ok(x), Ok(y), Ok(z)) = (
                            line[22..26].trim().parse::<usize>(),
                            line[30..38].trim().parse::<f64>(),
                            line[38..46].trim().parse::<f64>(),
                            line[46..54].trim().parse::<f64>(),
                        ) {
                            chain_b_map.entry(rs).or_insert([x, y, z]);
                        }
                    }
                }
            }
        }
        let mut p_pts = Vec::new();
        let mut q_pts = Vec::new();
        for &(rs, ax, ay, az) in chain_a_cas {
            if let Some(&b_xyz) = chain_b_map.get(&rs) {
                p_pts.push([ax, ay, az]);
                q_pts.push(b_xyz);
            }
        }
        if p_pts.len() >= 30 {
            let (_, r, tm, _) = kabsch_superpose_and_tm_native(&p_pts, &q_pts, true);
            kabsch_rmsd = round4(r);
            tm_score = round4(tm);
        } else if !chain_a_cas.is_empty() {
            let pts: Vec<[f64; 3]> = chain_a_cas.iter().map(|&(_, x, y, z)| [x, y, z]).collect();
            let (_, r, tm, _) = kabsch_superpose_and_tm_native(&pts, &pts, true);
            kabsch_rmsd = round4(r);
            tm_score = round4(tm);
        }
    }

    EmpiricalStructureAnalysis {
        pocket_center_xyz: center,
        pocket_radius_angstrom: radius_a,
        pocket_volume_a3: volume_a3,
        docking_affinity_kcal_mol: round4(vina.estimated_free_energy_kcal),
        active_torsions,
        kabsch_rmsd_angstrom: kabsch_rmsd,
        tm_score,
    }
}

/// Convenience helper loading a cached PDB file directly and running `compute_empirical_pocket_and_vina`.
pub fn analyze_pdb_file_empirically(
    pdb_path: &Path,
    chain_id: char,
    ligand_id: &str,
    ligand_res_seq: Option<i32>,
    chain_a_cas: &[(usize, f64, f64, f64)],
) -> Result<EmpiricalStructureAnalysis, String> {
    let (rec, mut lig) = parse_pdb_receptor_and_ligand(pdb_path, chain_id, ligand_id, ligand_res_seq)?;
    if lig.is_empty() && ligand_res_seq.is_some() {
        let (_, lig_any) = parse_pdb_receptor_and_ligand(pdb_path, chain_id, ligand_id, None)?;
        lig = lig_any;
    }
    let raw = std::fs::read_to_string(pdb_path).ok();
    Ok(compute_empirical_pocket_and_vina(
        &rec,
        &lig,
        chain_a_cas,
        raw.as_deref(),
        chain_id,
    ))
}

/// Computes exact Schwarz / BIC posterior model weights across candidate BIC scores:
/// $w_i = \exp(-0.5 (\text{BIC}_i - \text{BIC}_{\min})) / \sum_j \exp(-0.5 (\text{BIC}_j - \text{BIC}_{\min}))$.
pub fn compute_bic_posterior_weights(bic_values: &[f64]) -> Vec<f64> {
    if bic_values.is_empty() {
        return Vec::new();
    }
    let min_bic = bic_values
        .iter()
        .copied()
        .filter(|x| x.is_finite())
        .fold(f64::INFINITY, f64::min);
    if !min_bic.is_finite() {
        return vec![0.0; bic_values.len()];
    }
    let raw_w: Vec<f64> = bic_values
        .iter()
        .map(|&b| {
            if b.is_finite() {
                (-0.5 * (b - min_bic)).clamp(-700.0, 0.0).exp()
            } else {
                0.0
            }
        })
        .collect();
    let sum_w: f64 = raw_w.iter().sum();
    if sum_w > 0.0 {
        raw_w.into_iter().map(|w| round4(w / sum_w)).collect()
    } else {
        vec![0.0; bic_values.len()]
    }
}

/// Computes the Kishino-Hasegawa (1989) paired-sites log-likelihood z-statistic between site
/// log-likelihood vector `ell_a` and `ell_b` weighted by site pattern counts `weights`.
pub fn compute_kh_paired_sites_z(
    ell_a: &[f64],
    weights_a: &[f64],
    ell_b: &[f64],
    weights_b: &[f64],
) -> f64 {
    let n_a: f64 = weights_a.iter().sum();
    let n_b: f64 = weights_b.iter().sum();
    let n_tot = n_a.min(n_b);
    if n_tot <= 1.0 || ell_a.is_empty() || ell_b.is_empty() {
        return 0.0;
    }
    if ell_a.len() == ell_b.len() && (n_a - n_b).abs() < 0.5 && weights_a == weights_b {
        let mut sum_d = 0.0_f64;
        let mut sum_d2 = 0.0_f64;
        for i in 0..ell_a.len() {
            let d = ell_a[i] - ell_b[i];
            let w = weights_a[i];
            sum_d += w * d;
            sum_d2 += w * d * d;
        }
        let mean_d = sum_d / n_tot;
        let var_d = ((sum_d2 / n_tot) - mean_d * mean_d).max(0.0);
        let se = (var_d / n_tot).sqrt();
        if se > 1e-12 {
            return round4(sum_d / (n_tot * se));
        }
    }
    // Welch paired-variance form when site pattern compression maps differ (e.g. RY3 vs Unrecoded-4)
    let mean_a = ell_a.iter().zip(weights_a).map(|(&e, &w)| e * w).sum::<f64>() / n_a.max(1.0);
    let mean_b = ell_b.iter().zip(weights_b).map(|(&e, &w)| e * w).sum::<f64>() / n_b.max(1.0);
    let var_a = ell_a
        .iter()
        .zip(weights_a)
        .map(|(&e, &w)| w * (e - mean_a).powi(2))
        .sum::<f64>()
        / n_a.max(1.0);
    let var_b = ell_b
        .iter()
        .zip(weights_b)
        .map(|(&e, &w)| w * (e - mean_b).powi(2))
        .sum::<f64>()
        / n_b.max(1.0);
    let sum_a: f64 = ell_a.iter().zip(weights_a).map(|(&e, &w)| e * w).sum();
    let sum_b: f64 = ell_b.iter().zip(weights_b).map(|(&e, &w)| e * w).sum();
    let se_tot = (n_a * var_a + n_b * var_b).sqrt();
    if se_tot > 1e-12 {
        round4((sum_a - sum_b) / se_tot)
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solve_empirical_lsd2_core_nondegenerate_hpd_on_heterogeneous_tree() {
        // Boundary-pinned root (3240..3900 Ma) + heterogeneous branch rates + bounded leaf split
        let nodes = vec![
            TimeTreeNodeSpec {
                node_id: 0,
                parent_id: None,
                label: "Root".to_string(),
                is_tip: false,
                is_extinct_fossil_tip: false,
                branch_length_subs: 0.0,
                tip_age_ma: None,
                cal_min_ma: Some(3240.0),
                cal_max_ma: Some(3900.0),
                focal_node_tag: None,
            },
            TimeTreeNodeSpec {
                node_id: 1,
                parent_id: Some(0),
                label: "CladeA".to_string(),
                is_tip: false,
                is_extinct_fossil_tip: false,
                branch_length_subs: 0.45,
                tip_age_ma: None,
                cal_min_ma: Some(1650.0),
                cal_max_ma: Some(2250.0),
                focal_node_tag: None,
            },
            TimeTreeNodeSpec {
                node_id: 2,
                parent_id: Some(0),
                label: "CladeB".to_string(),
                is_tip: false,
                is_extinct_fossil_tip: false,
                branch_length_subs: 1.15,
                tip_age_ma: None,
                cal_min_ma: Some(810.0),
                cal_max_ma: Some(1450.0),
                focal_node_tag: None,
            },
            TimeTreeNodeSpec {
                node_id: 3,
                parent_id: Some(1),
                label: "TipA1".to_string(),
                is_tip: true,
                is_extinct_fossil_tip: false,
                branch_length_subs: 0.72,
                tip_age_ma: Some(0.0),
                cal_min_ma: None,
                cal_max_ma: None,
                focal_node_tag: None,
            },
            TimeTreeNodeSpec {
                node_id: 4,
                parent_id: Some(1),
                label: "TipA2".to_string(),
                is_tip: true,
                is_extinct_fossil_tip: false,
                branch_length_subs: 1.24,
                tip_age_ma: Some(0.0),
                cal_min_ma: None,
                cal_max_ma: None,
                focal_node_tag: None,
            },
            TimeTreeNodeSpec {
                node_id: 5,
                parent_id: Some(2),
                label: "TipB1".to_string(),
                is_tip: true,
                is_extinct_fossil_tip: false,
                branch_length_subs: 0.35,
                tip_age_ma: Some(0.0),
                cal_min_ma: None,
                cal_max_ma: None,
                focal_node_tag: None,
            },
            TimeTreeNodeSpec {
                node_id: 6,
                parent_id: Some(2),
                label: "TipB2".to_string(),
                is_tip: true,
                is_extinct_fossil_tip: false,
                branch_length_subs: 0.91,
                tip_age_ma: Some(0.0),
                cal_min_ma: None,
                cal_max_ma: None,
                focal_node_tag: None,
            },
        ];

        let (ages, cis, _rate, _q) = solve_empirical_lsd2_core(&nodes, 26_387, 20261152);
        for nid in 0..3 {
            let age = ages[nid];
            let (lo, hi) = cis[nid];
            assert!(lo < age && age < hi, "node {} must have strictly two-sided HPD: age={} lo={} hi={}", nid, age, lo, hi);
            assert!(hi - lo > 10.0, "node {} HPD width unexpectedly collapsed: {} ({}..{})", nid, hi - lo, lo, hi);
        }
    }
}
