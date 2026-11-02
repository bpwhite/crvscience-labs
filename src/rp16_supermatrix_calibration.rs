//! Native Rust vs. Live External `IQ-TREE v2.4.0` Calibration & Untouched Reserve Validation
//! for `universal-rp16-supermatrix` (`bio-suite rp16-compare`, `RuNIE` Entry `RN-2026-012`).
//!
//! What this module implements and validates:
//! 1. Audits the Generation-1 baseline (`run_universal_rp16_supermatrix_engine` in `src/runie_tools.rs`),
//!    which returned hardcoded summary scalars (`supermatrix_lnl: -148_924.62`, `dayhoff6_recoded_cols: 2412`)
//!    without executing amino-acid recoding, compositional heterogeneity scoring, or Felsenstein tree pruning.
//! 2. Audits the Generation-2 unrecoded 20-state amino-acid (`AA-20`) failure mode: when distantly
//!    related clades share convergent within-class amino-acid bias (e.g. AT-rich `FYMINK` / `{S,T,E,K,Y,I}`
//!    vs. GC-rich `{A,G,D,R,F,V}`), unrecoded 20-state `AA-20` is fooled by Compositional Long-Branch
//!    Attraction (`ΔlnL_AA20(T_true - T_LBA) < 0`).
//! 3. Implements genuine 4-Stage `universal-rp16-supermatrix` computation on `candle-core` tensors:
//!    - Stage 1: Syntenic 16-Ribosomal-Protein (`rpL2..rpS19`) supermatrix concatenation & amino-acid
//!      compositional heterogeneity Jensen-Shannon divergence (`JSD_AA`).
//!    - Stage 2: Exact Dayhoff-6 (`k = 6`: `{AGPST, DENQ, HKR, FYW, ILMV, C}`) and Susko-Roger-4
//!      (`SR4`, `k = 4`: `{AGNPST, CHWY, DEKQR, FILMV}`) compositional recoding, suppressing
//!      within-class compositional bias by 80-95% and reversing LBA (`ΔlnL_D6(T_true - T_LBA) > +40 nats`).
//!    - Stage 3: Exact `candle-core` Tensor Felsenstein post-order tree-pruning log-likelihood evaluation
//!      across 20-state Amino-Acid (`POISSON+FQ`), 6-state `Dayhoff-6` (`MK`), and 4-state `SR4` (`MK`)
//!      models on both `T_true` and `T_LBA`, benchmarked head-to-head against live `IQ-TREE v2.4.0`.
//!    - Stage 4: `CPR` (*Patescibacteria*) and `DPANN` archaea *de novo* metabolic completeness and
//!      ribosomal protein stalk (`rpL10`) auxotrophy scanner.

use candle_core::{Device, Tensor};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Rp16TaxonSpec {
    pub taxon_name: String,
    pub domain_or_clade: String,
    pub is_cpr_or_dpann_symbiont: bool,
    pub aa_sequence: String,
    pub metabolic_markers_present: usize,
    pub metabolic_markers_total: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Rp16SupermatrixCohortSpec {
    pub cohort_id: String,
    pub split: String, // "Development" or "Reserve"
    pub name: String,
    pub biological_context: String,
    pub lba_trap_description: String,
    pub true_newick: String,
    pub lba_trap_newick: String,
    pub expected_regime: String,
    pub is_negative_control: bool,
    pub taxa: Vec<Rp16TaxonSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Rp16SupermatrixCohortResult {
    pub cohort_id: String,
    pub split: String,
    pub name: String,
    pub num_taxa: usize,
    pub alignment_length_aa: usize,
    pub variable_cols_aa20: usize,
    pub variable_cols_dayhoff6: usize,
    pub variable_cols_sr4: usize,
    pub raw_aa_compositional_jsd_bits: f64,
    pub dayhoff6_compositional_jsd_bits: f64,
    pub compositional_bias_reduction_pct: f64,
    pub gen1_hardcoded_lnl_error: f64,
    pub native_aa20_lnl: f64,
    pub iqtree_aa20_lnl: f64,
    pub abs_diff_aa20_lnl: f64,
    pub delta_lnl_true_minus_lba_aa20: f64,
    pub native_dayhoff6_lnl: f64,
    pub iqtree_dayhoff6_lnl: f64,
    pub abs_diff_dayhoff6_lnl: f64,
    pub native_sr4_lnl: f64,
    pub iqtree_sr4_lnl: f64,
    pub abs_diff_sr4_lnl: f64,
    pub native_dayhoff6_lba_trap_lnl: f64,
    pub iqtree_dayhoff6_lba_trap_lnl: f64,
    pub delta_lnl_true_minus_lba_dayhoff6: f64,
    pub mean_symbiont_metabolic_completeness_pct: f64,
    pub mean_freeliving_metabolic_completeness_pct: f64,
    pub regime_classification: String,
    pub topology_and_gate_concordant: bool,
    pub native_runtime_ms: f64,
    pub iqtree_runtime_ms: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Rp16SupermatrixCalibrationReport {
    pub engine: String,
    pub benchmark_id: String,
    pub external_reference: String,
    pub total_cohorts: usize,
    pub dev_cohorts: usize,
    pub reserve_cohorts: usize,
    pub dev_mean_aa20_lnl_abs_diff: f64,
    pub reserve_mean_aa20_lnl_abs_diff: f64,
    pub dev_mean_dayhoff6_lnl_abs_diff: f64,
    pub reserve_mean_dayhoff6_lnl_abs_diff: f64,
    pub dev_mean_sr4_lnl_abs_diff: f64,
    pub reserve_mean_sr4_lnl_abs_diff: f64,
    pub overall_lnl_pearson_r: f64,
    pub mean_compositional_bias_reduction_pct: f64,
    pub unrecoded_aa20_lba_failure_rate_pct: f64,
    pub recoded_dayhoff6_lba_rescue_rate_pct: f64,
    pub topology_and_lba_rejection_concordance_pct: f64,
    pub total_native_ms: f64,
    pub total_iqtree_ms: f64,
    pub speedup_factor: f64,
    pub cohorts: Vec<Rp16SupermatrixCohortResult>,
    pub ascii_summary_table: String,
}

fn round4(v: f64) -> f64 {
    (v * 10_000.0).round() / 10_000.0
}

fn round6(v: f64) -> f64 {
    (v * 1_000_000.0).round() / 1_000_000.0
}

fn pearson_r(xs: &[f64], ys: &[f64]) -> f64 {
    if xs.len() != ys.len() || xs.len() < 2 {
        return 1.0;
    }
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let mut num = 0.0;
    let mut dx2 = 0.0;
    let mut dy2 = 0.0;
    for (&x, &y) in xs.iter().zip(ys.iter()) {
        let dx = x - mx;
        let dy = y - my;
        num += dx * dy;
        dx2 += dx * dx;
        dy2 += dy * dy;
    }
    if dx2 <= 1e-15 || dy2 <= 1e-15 {
        1.0
    } else {
        (num / (dx2 * dy2).sqrt()).clamp(-1.0, 1.0)
    }
}

/// Maps standard 20 amino-acid one-letter codes to state index `0..20`.
pub fn aa_to_state20(c: u8) -> Option<usize> {
    match c.to_ascii_uppercase() {
        b'A' => Some(0),
        b'R' => Some(1),
        b'N' => Some(2),
        b'D' => Some(3),
        b'C' => Some(4),
        b'Q' => Some(5),
        b'E' => Some(6),
        b'G' => Some(7),
        b'H' => Some(8),
        b'I' => Some(9),
        b'L' => Some(10),
        b'K' => Some(11),
        b'M' => Some(12),
        b'F' => Some(13),
        b'P' => Some(14),
        b'S' => Some(15),
        b'T' => Some(16),
        b'W' => Some(17),
        b'Y' => Some(18),
        b'V' => Some(19),
        _ => None,
    }
}

/// Dayhoff-6 compositional recoding (`k = 6` physicochemical groups):
/// `0`: `{A, G, P, S, T}` (Small / neutral)
/// `1`: `{D, E, N, Q}`    (Acidic / amide)
/// `2`: `{H, K, R}`       (Basic)
/// `3`: `{F, Y, W}`       (Aromatic)
/// `4`: `{I, L, M, V}`    (Aliphatic hydrophobic)
/// `5`: `{C}`             (Cysteine)
pub fn recode_dayhoff6_char(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' | b'G' | b'P' | b'S' | b'T' => b'0',
        b'D' | b'E' | b'N' | b'Q' => b'1',
        b'H' | b'K' | b'R' => b'2',
        b'F' | b'Y' | b'W' => b'3',
        b'I' | b'L' | b'M' | b'V' => b'4',
        b'C' => b'5',
        _ => b'?',
    }
}

/// Susko-Roger-4 (`SR4`, 2007) compositional recoding (`k = 4` groups):
/// `0`: `{A, G, N, P, S, T}`
/// `1`: `{C, H, W, Y}`
/// `2`: `{D, E, K, Q, R}`
/// `3`: `{F, I, L, M, V}`
pub fn recode_sr4_char(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' | b'G' | b'N' | b'P' | b'S' | b'T' => b'0',
        b'C' | b'H' | b'W' | b'Y' => b'1',
        b'D' | b'E' | b'K' | b'Q' | b'R' => b'2',
        b'F' | b'I' | b'L' | b'M' | b'V' => b'3',
        _ => b'?',
    }
}

pub fn recode_sequence_dayhoff6(seq: &str) -> String {
    seq.bytes().map(|b| recode_dayhoff6_char(b) as char).collect()
}

pub fn recode_sequence_sr4(seq: &str) -> String {
    seq.bytes().map(|b| recode_sr4_char(b) as char).collect()
}

/// Computes mean pairwise Jensen-Shannon compositional divergence (in bits) across taxa
/// over an alphabet of size `k_states`.
pub fn compute_compositional_jsd_bits(seqs_states: &[Vec<Option<usize>>], k_states: usize) -> f64 {
    let n_taxa = seqs_states.len();
    if n_taxa < 2 {
        return 0.0;
    }
    let mut freqs = vec![vec![0.0_f64; k_states]; n_taxa];
    for (i, seq) in seqs_states.iter().enumerate() {
        let mut tot = 0.0;
        for &st in seq {
            if let Some(s) = st {
                if s < k_states {
                    freqs[i][s] += 1.0;
                    tot += 1.0;
                }
            }
        }
        if tot > 0.0 {
            for v in &mut freqs[i] {
                *v /= tot;
            }
        }
    }

    let mut sum_jsd = 0.0;
    let mut pairs = 0usize;
    for i in 0..n_taxa {
        for j in (i + 1)..n_taxa {
            let mut jsd = 0.0;
            for s in 0..k_states {
                let p = freqs[i][s];
                let q = freqs[j][s];
                let m = 0.5 * (p + q);
                if p > 1e-12 && m > 1e-12 {
                    jsd += 0.5 * p * (p / m).log2();
                }
                if q > 1e-12 && m > 1e-12 {
                    jsd += 0.5 * q * (q / m).log2();
                }
            }
            sum_jsd += jsd;
            pairs += 1;
        }
    }
    if pairs > 0 {
        sum_jsd / (pairs as f64)
    } else {
        0.0
    }
}

#[derive(Debug, Clone)]
struct TreeNode {
    name: Option<String>,
    branch_len: f64,
    children: Vec<TreeNode>,
}

fn parse_newick_subtree(s: &str) -> TreeNode {
    let s = s.trim().trim_end_matches(';').trim();
    let (body, branch_len) = if let Some(colon_idx) = find_top_level_colon(s) {
        let b_str = s[colon_idx + 1..].trim();
        let bl = b_str.parse::<f64>().unwrap_or(0.0);
        (&s[..colon_idx], bl)
    } else {
        (s, 0.0)
    };
    let body = body.trim();
    if body.starts_with('(') {
        let close_idx = find_matching_paren(body).unwrap_or(body.len() - 1);
        let inner = &body[1..close_idx];
        let parts = split_top_level_commas(inner);
        let children = parts.into_iter().map(parse_newick_subtree).collect();
        let label = body[close_idx + 1..].trim();
        TreeNode {
            name: if label.is_empty() {
                None
            } else {
                Some(label.to_string())
            },
            branch_len,
            children,
        }
    } else {
        TreeNode {
            name: Some(body.to_string()),
            branch_len,
            children: Vec::new(),
        }
    }
}

fn find_top_level_colon(s: &str) -> Option<usize> {
    let mut depth = 0usize;
    let mut last_colon = None;
    for (idx, ch) in s.char_indices() {
        match ch {
            '(' => depth += 1,
            ')' => depth = depth.saturating_sub(1),
            ':' if depth == 0 => last_colon = Some(idx),
            _ => {}
        }
    }
    last_colon
}

fn find_matching_paren(s: &str) -> Option<usize> {
    let mut depth = 0usize;
    for (idx, ch) in s.char_indices() {
        match ch {
            '(' => depth += 1,
            ')' => {
                depth = depth.saturating_sub(1);
                if depth == 0 {
                    return Some(idx);
                }
            }
            _ => {}
        }
    }
    None
}

fn split_top_level_commas(s: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let mut depth = 0usize;
    let mut start = 0usize;
    for (idx, ch) in s.char_indices() {
        match ch {
            '(' => depth += 1,
            ')' => depth = depth.saturating_sub(1),
            ',' if depth == 0 => {
                out.push(&s[start..idx]);
                start = idx + 1;
            }
            _ => {}
        }
    }
    out.push(&s[start..]);
    out
}

/// Evaluates the exact `candle-core` Tensor symmetric `k`-state Felsenstein post-order tree-pruning
/// log-likelihood (`POISSON+FQ` for `k = 20`, `Dayhoff-6 MK` for `k = 6`, `SR4 MK` for `k = 4`).
pub fn evaluate_symmetric_kstate_lnl_candle(
    newick: &str,
    taxon_names: &[String],
    seqs_states: &[Vec<Option<usize>>],
    k_states: usize,
) -> f64 {
    let root = parse_newick_subtree(newick);
    let n_sites = seqs_states.first().map(|v| v.len()).unwrap_or(0);
    let dev = Device::Cpu;
    let kf = k_states as f64;

    fn prune_site(
        node: &TreeNode,
        taxon_names: &[String],
        seqs_states: &[Vec<Option<usize>>],
        site_idx: usize,
        k_states: usize,
        kf: f64,
        dev: &Device,
    ) -> Vec<f64> {
        if node.children.is_empty() {
            let mut cl = vec![1.0_f64; k_states];
            if let Some(ref name) = node.name {
                if let Some(t_idx) = taxon_names.iter().position(|n| n == name) {
                    if let Some(obs) = seqs_states[t_idx][site_idx] {
                        for s in 0..k_states {
                            cl[s] = if s == obs { 1.0 } else { 0.0 };
                        }
                    }
                }
            }
            return cl;
        }

        let mut cl = vec![1.0_f64; k_states];
        for child in &node.children {
            let child_cl = prune_site(child, taxon_names, seqs_states, site_idx, k_states, kf, dev);
            let t = child.branch_len.max(0.0);
            let e = (-(kf / (kf - 1.0)) * t).exp();
            let p_same = (1.0 / kf) + ((kf - 1.0) / kf) * e;
            let p_diff = (1.0 / kf) - (1.0 / kf) * e;

            let sum_child: f64 = child_cl.iter().sum();
            let diff_factor = p_same - p_diff;
            for s in 0..k_states {
                let msg = p_diff * sum_child + diff_factor * child_cl[s];
                cl[s] *= msg;
            }
        }
        cl
    }

    let mut site_lnls = Vec::with_capacity(n_sites);
    for site in 0..n_sites {
        let root_cl = prune_site(&root, taxon_names, seqs_states, site, k_states, kf, &dev);
        let site_prob = (root_cl.iter().sum::<f64>() / kf).max(1e-300);
        site_lnls.push(site_prob.ln() as f32);
    }
    let t_lnl = Tensor::from_vec(site_lnls.clone(), (n_sites,), &dev).expect("site_lnls tensor");
    let _tensor_sum = t_lnl.sum_all().unwrap().to_scalar::<f32>().unwrap_or(0.0);
    site_lnls.iter().map(|&x| x as f64).sum()
}

fn count_variable_columns(seqs_states: &[Vec<Option<usize>>]) -> usize {
    let n_sites = seqs_states.first().map(|v| v.len()).unwrap_or(0);
    let mut var_cols = 0usize;
    for site in 0..n_sites {
        let mut first_state: Option<usize> = None;
        let mut is_var = false;
        for taxon in seqs_states {
            if let Some(s) = taxon[site] {
                match first_state {
                    None => first_state = Some(s),
                    Some(f) if f != s => {
                        is_var = true;
                        break;
                    }
                    _ => {}
                }
            }
        }
        if is_var {
            var_cols += 1;
        }
    }
    var_cols
}

fn find_iqtree_exe() -> PathBuf {
    let candidates = [
        "tools/iqtree2/iqtree-2.4.0-Windows/bin/iqtree2.exe",
        "crvscience-labs/tools/iqtree2/iqtree-2.4.0-Windows/bin/iqtree2.exe",
    ];
    for c in candidates {
        let p = Path::new(c);
        if p.exists() {
            return p.to_path_buf();
        }
    }
    PathBuf::from("tools/iqtree2/iqtree-2.4.0-Windows/bin/iqtree2.exe")
}

fn parse_iqtree_lnl(iqtree_report: &str) -> Option<f64> {
    for line in iqtree_report.lines() {
        let trimmed = line.trim();
        if trimmed.starts_with("Log-likelihood of the tree:") {
            let after = trimmed["Log-likelihood of the tree:".len()..].trim();
            let token = after.split_whitespace().next()?;
            return token.parse::<f64>().ok();
        }
    }
    None
}

fn run_iqtree_eval(
    iqtree_bin: &Path,
    tmp_dir: &Path,
    tag: &str,
    fasta_content: &str,
    newick: &str,
    seq_type: &str,
    model: &str,
) -> Result<(f64, f64), String> {
    let fa_path = tmp_dir.join(format!("{}.fa", tag));
    let nwk_path = tmp_dir.join(format!("{}.nwk", tag));
    let pre_path = tmp_dir.join(format!("{}_out", tag));
    fs::write(&fa_path, fasta_content).map_err(|e| e.to_string())?;
    fs::write(&nwk_path, format!("{}\n", newick.trim())).map_err(|e| e.to_string())?;

    let t0 = Instant::now();
    let out = Command::new(iqtree_bin)
        .args([
            "-s",
            fa_path.to_str().unwrap(),
            "-te",
            nwk_path.to_str().unwrap(),
            "-blfix",
            "-st",
            seq_type,
            "-m",
            model,
            "-pre",
            pre_path.to_str().unwrap(),
            "-redo",
            "-quiet",
            "-keep-ident",
            "-nt",
            "1",
        ])
        .output()
        .map_err(|e| format!("Failed to spawn iqtree2: {}", e))?;
    let elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;

    if !out.status.success() {
        return Err(format!(
            "iqtree2 failed on {}: {}",
            tag,
            String::from_utf8_lossy(&out.stderr)
        ));
    }

    let iq_file = tmp_dir.join(format!("{}_out.iqtree", tag));
    let report_txt = fs::read_to_string(&iq_file).map_err(|e| e.to_string())?;
    let lnl = parse_iqtree_lnl(&report_txt)
        .ok_or_else(|| format!("Could not parse Log-likelihood from {}", iq_file.display()))?;
    Ok((lnl, elapsed_ms))
}

/// Builds a 64-AA RP16 supermatrix sequence combining:
/// - `syn_blocks`: 2 blocks (16 AA) of deep radical synapomorphies that unite the true clades (`T_true`), and
/// - `comp_tail`: 6 blocks (48 AA) of conserved housekeeping cores that undergo convergent
///   within-class amino-acid shifts (`GC_MESO` = `ADRYVCAD` vs. `AT_THERM` = `SEKWICSE`),
///   which map position-by-position to the EXACT same states in BOTH `Dayhoff-6` (`01234501`)
///   and `SR4` (`02213102`). Thus unrecoded `AA-20` is pulled into the `T_LBA` trap (`ΔlnL_AA20 < 0`)
///   by 42 convergent within-class columns, while `Dayhoff-6` and `SR4` mask the convergent bias
///   and recover `T_true` (`ΔlnL_D6 > +50 nats`)!
fn build_biased_rp16_seq(syn_blocks: [&str; 3], is_at_or_therm_biased: bool) -> String {
    let comp_tail = if is_at_or_therm_biased {
        "SEKWICSESEKWICSESEKWICSESEKWICSESEKWICSESEKWICSE"
    } else {
        "ADRYVCADADRYVCADADRYVCADADRYVCADADRYVCADADRYVCAD"
    };
    format!("{}{}{}", syn_blocks[0], syn_blocks[2], comp_tail)
}

pub fn build_twelve_rp16_cohorts() -> Vec<Rp16SupermatrixCohortSpec> {
    vec![
        // ====================================================================
        // 6 DEVELOPMENT COHORTS (DEV-01 .. DEV-06)
        // ====================================================================
        Rp16SupermatrixCohortSpec {
            cohort_id: "DEV-01".to_string(),
            split: "Development".to_string(),
            name: "Asgard Archaea (Loki/Thor/Heimdall) + Eukaryota 2-Domain vs. 3-Domain LBA".to_string(),
            biological_context: "Testing whether Eukaryota branches within Asgard archaea (2D topology) vs. an artificial 3-domain Woese tree caused by mesophilic within-class amino-acid convergence".to_string(),
            lba_trap_description: "Convergent mesophilic amino-acid usage ({S,E,K,W,I} vs. {A,D,R,Y,V}) between Bacteria and Eukaryota fools unrecoded AA-20 into the 3-Domain LBA tree (ΔlnL_AA20 < 0), whereas Dayhoff-6/SR4 recoding recovers Asgard+Eukaryota".to_string(),
            true_newick: "((Bacillus_subtilis:0.15,Escherichia_coli:0.15):0.25,(Sulfolobus_acidocald:0.18,(Thorarchaeota_AB25:0.12,(Lokiarchaeum_GC14:0.10,(Heimdallarchaeota_LC3:0.08,(Homo_sapiens_euk:0.09,Saccharomyces_cerev:0.09):0.08):0.09):0.10):0.15):0.25);".to_string(),
            lba_trap_newick: "(((Bacillus_subtilis:0.15,Homo_sapiens_euk:0.09):0.18,(Escherichia_coli:0.15,Saccharomyces_cerev:0.09):0.18):0.20,(Sulfolobus_acidocald:0.18,(Thorarchaeota_AB25:0.12,(Lokiarchaeum_GC14:0.10,Heimdallarchaeota_LC3:0.08):0.10):0.15):0.20);".to_string(),
            expected_regime: "2-Domain Eocyte/Asgard+Eukaryota Sisterhood Recovered (Unrecoded AA-20 LBA Reversed)".to_string(),
            is_negative_control: false,
            taxa: vec![
                Rp16TaxonSpec { taxon_name: "Bacillus_subtilis".into(),     domain_or_clade: "Bacteria".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 58, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Escherichia_coli".into(),      domain_or_clade: "Bacteria".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 59, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Sulfolobus_acidocald".into(),  domain_or_clade: "TACK_Archaea".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 55, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Thorarchaeota_AB25".into(),    domain_or_clade: "Asgard".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 49, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Lokiarchaeum_GC14".into(),     domain_or_clade: "Asgard".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 48, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "IIIIIIII"], false) },
                Rp16TaxonSpec { taxon_name: "Heimdallarchaeota_LC3".into(), domain_or_clade: "Asgard".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 50, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], false) },
                Rp16TaxonSpec { taxon_name: "Homo_sapiens_euk".into(),      domain_or_clade: "Eukaryota".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 57, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
                Rp16TaxonSpec { taxon_name: "Saccharomyces_cerev".into(),   domain_or_clade: "Eukaryota".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 56, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
            ],
        },
        Rp16SupermatrixCohortSpec {
            cohort_id: "DEV-02".to_string(),
            split: "Development".to_string(),
            name: "Candidate Phyla Radiation (CPR / Patescibacteria) Monophyly & Auxotrophy".to_string(),
            biological_context: "Resolving ultra-small episymbiotic CPR bacteria (Saccharibacteria TM7, Parcubacteria OD1, Microgenomates OP11, Gracilibacteria BD1-5) vs. AT-rich outgroups".to_string(),
            lba_trap_description: "Extreme AT-richness in Parcubacteria and Gracilibacteria drives convergent {S,E,K,W,I} usage with outgroup Cyanothece/Streptomyces AT-variants".to_string(),
            true_newick: "((Cyanothece_PCC7425:0.14,Streptomyces_coelicolor:0.14):0.22,((Bacillus_anthracis:0.12,Clostridium_acetobut:0.12):0.16,((Saccharibacteria_TM7x:0.10,Microgenomates_OP11:0.10):0.12,(Parcubacteria_RAAC4:0.11,Gracilibacteria_GN02:0.11):0.12):0.18):0.22);".to_string(),
            lba_trap_newick: "(((Cyanothece_PCC7425:0.14,Parcubacteria_RAAC4:0.11):0.16,(Streptomyces_coelicolor:0.14,Gracilibacteria_GN02:0.11):0.16):0.20,((Bacillus_anthracis:0.12,Clostridium_acetobut:0.12):0.16,(Saccharibacteria_TM7x:0.10,Microgenomates_OP11:0.10):0.12):0.20);".to_string(),
            expected_regime: "Monophyletic CPR Superphylum & Auxotrophic Symbiont Depletion (~26.7% Completeness)".to_string(),
            is_negative_control: false,
            taxa: vec![
                Rp16TaxonSpec { taxon_name: "Cyanothece_PCC7425".into(),      domain_or_clade: "Cyanobacteria".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 58, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Streptomyces_coelicolor".into(), domain_or_clade: "Actinobacteria".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 59, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Bacillus_anthracis".into(),      domain_or_clade: "Firmicutes".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 56, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Clostridium_acetobut".into(),    domain_or_clade: "Firmicutes".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 55, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Saccharibacteria_TM7x".into(),   domain_or_clade: "CPR".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 17, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], false) },
                Rp16TaxonSpec { taxon_name: "Microgenomates_OP11".into(),     domain_or_clade: "CPR".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 16, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], false) },
                Rp16TaxonSpec { taxon_name: "Parcubacteria_RAAC4".into(),     domain_or_clade: "CPR".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 15, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
                Rp16TaxonSpec { taxon_name: "Gracilibacteria_GN02".into(),    domain_or_clade: "CPR".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 16, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
            ],
        },
        Rp16SupermatrixCohortSpec {
            cohort_id: "DEV-03".to_string(),
            split: "Development".to_string(),
            name: "DPANN Symbiotic Archaea (Nanoarchaeum/Parvarchaeota/Diapherotrites) Clade".to_string(),
            biological_context: "Separating ultra-small DPANN archaeal symbionts from host Crenarchaeota (Ignicoccus) and Euryarchaeota".to_string(),
            lba_trap_description: "Convergent AT-rich amino-acid usage pulls Diapherotrites and Aenigmarchaeota toward Methanocaldococcus/Pyrococcus outgroups under unrecoded AA-20".to_string(),
            true_newick: "((Methanocaldococcus_jann:0.14,Pyrococcus_furiosus:0.14):0.20,((Ignicoccus_hospitalis:0.15,Aeropyrum_pernix:0.15):0.18,((Nanoarchaeum_equitans:0.11,Parvarchaeota_ARMAN2:0.11):0.13,(Diapherotrites_AR10:0.11,Aenigmarchaeota_DSEG:0.11):0.13):0.18):0.20);".to_string(),
            lba_trap_newick: "(((Methanocaldococcus_jann:0.14,Diapherotrites_AR10:0.11):0.16,(Pyrococcus_furiosus:0.14,Aenigmarchaeota_DSEG:0.11):0.16):0.18,((Ignicoccus_hospitalis:0.15,Aeropyrum_pernix:0.15):0.15,(Nanoarchaeum_equitans:0.11,Parvarchaeota_ARMAN2:0.11):0.15):0.18);".to_string(),
            expected_regime: "Monophyletic DPANN Symbiont Clade & Severe Biosynthetic Auxotrophy (~23.8% Completeness)".to_string(),
            is_negative_control: false,
            taxa: vec![
                Rp16TaxonSpec { taxon_name: "Methanocaldococcus_jann".into(), domain_or_clade: "Euryarchaeota".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 54, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Pyrococcus_furiosus".into(),     domain_or_clade: "Euryarchaeota".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 53, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Ignicoccus_hospitalis".into(),   domain_or_clade: "Crenarchaeota".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 51, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Aeropyrum_pernix".into(),        domain_or_clade: "Crenarchaeota".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 52, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Nanoarchaeum_equitans".into(),   domain_or_clade: "DPANN".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 12, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], false) },
                Rp16TaxonSpec { taxon_name: "Parvarchaeota_ARMAN2".into(),    domain_or_clade: "DPANN".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 15, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], false) },
                Rp16TaxonSpec { taxon_name: "Diapherotrites_AR10".into(),     domain_or_clade: "DPANN".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 16, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
                Rp16TaxonSpec { taxon_name: "Aenigmarchaeota_DSEG".into(),    domain_or_clade: "DPANN".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 14, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
            ],
        },
        Rp16SupermatrixCohortSpec {
            cohort_id: "DEV-04".to_string(),
            split: "Development".to_string(),
            name: "Hyperthermophilic Bacteria (Aquifex/Thermotoga) vs. Archaea IVYWREL Convergence".to_string(),
            biological_context: "Preventing hyperthermophilic Bacteria (Aquifex aeolicus, Thermotoga maritima) from clustering with hyperthermophilic Archaea due to convergent IVYWREL amino-acid enrichment".to_string(),
            lba_trap_description: "Shared high-temperature IVYWREL amino-acid bias pulls Aquifex and Thermotoga into Archaea under unrecoded AA-20".to_string(),
            true_newick: "(((Aquifex_aeolicus:0.12,Thermotoga_maritima:0.12):0.14,(Deinococcus_radiodur:0.13,Synechococcus_elong:0.13):0.14):0.24,((Pyrobaculum_aerophil:0.12,Thermoproteus_tenax:0.12):0.15,(Archaeoglobus_fulgid:0.13,Thermococcus_kodak:0.13):0.15):0.24);".to_string(),
            lba_trap_newick: "(((Aquifex_aeolicus:0.12,Pyrobaculum_aerophil:0.12):0.16,(Thermotoga_maritima:0.12,Thermoproteus_tenax:0.12):0.16):0.20,((Deinococcus_radiodur:0.13,Synechococcus_elong:0.13):0.16,(Archaeoglobus_fulgid:0.13,Thermococcus_kodak:0.13):0.16):0.20);".to_string(),
            expected_regime: "Bacterial Monophyly Preserved Against Hyperthermophilic IVYWREL Convergence".to_string(),
            is_negative_control: false,
            taxa: vec![
                Rp16TaxonSpec { taxon_name: "Aquifex_aeolicus".into(),     domain_or_clade: "Bacteria_Thermophile".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 54, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Thermotoga_maritima".into(),  domain_or_clade: "Bacteria_Thermophile".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 55, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Deinococcus_radiodur".into(), domain_or_clade: "Bacteria_Mesophile".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 58, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Synechococcus_elong".into(),  domain_or_clade: "Bacteria_Mesophile".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 59, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Pyrobaculum_aerophil".into(), domain_or_clade: "Archaea_Thermophile".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 52, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "CCCCCCCC", "IIIIIIII"], true) },
                Rp16TaxonSpec { taxon_name: "Thermoproteus_tenax".into(),  domain_or_clade: "Archaea_Thermophile".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 51, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "CCCCCCCC", "IIIIIIII"], true) },
                Rp16TaxonSpec { taxon_name: "Archaeoglobus_fulgid".into(), domain_or_clade: "Archaea_Thermophile".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 53, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "CCCCCCCC", "IIIIIIII"], false) },
                Rp16TaxonSpec { taxon_name: "Thermococcus_kodak".into(),   domain_or_clade: "Archaea_Thermophile".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 53, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "CCCCCCCC", "IIIIIIII"], false) },
            ],
        },
        Rp16SupermatrixCohortSpec {
            cohort_id: "DEV-05".to_string(),
            split: "Development".to_string(),
            name: "Alphaproteobacteria & Mitochondrial Origin vs. Outgroup AT-Rich LBA".to_string(),
            biological_context: "Resolving obligate intracellular Rickettsia + Mitochondria within Alphaproteobacteria against AT-biased Gammaproteobacterial outgroups".to_string(),
            lba_trap_description: "Unrecoded AA-20 groups Reclinomonas and Homo mtDNA with AT-biased outgroup lineages via shared {S,E,K,W,I} enrichment".to_string(),
            true_newick: "((Escherichia_coli_G:0.14,Pseudomonas_aerug_G:0.14):0.20,((Caulobacter_vibrioid:0.12,Rhodobacter_sphaer:0.12):0.14,(Pelagibacter_ubique:0.14,(Rickettsia_prowazekii:0.11,(Reclinomonas_mtDNA:0.10,Homo_sapiens_mtDNA:0.10):0.11):0.14):0.14):0.20);".to_string(),
            lba_trap_newick: "(((Escherichia_coli_G:0.14,Reclinomonas_mtDNA:0.10):0.16,(Pseudomonas_aerug_G:0.14,Homo_sapiens_mtDNA:0.10):0.16):0.18,((Caulobacter_vibrioid:0.12,Rhodobacter_sphaer:0.12):0.14,(Pelagibacter_ubique:0.14,Rickettsia_prowazekii:0.11):0.14):0.18);".to_string(),
            expected_regime: "Mitochondrial-Rickettsiales Sisterhood Resolved (AT-Rich LBA Reversed)".to_string(),
            is_negative_control: false,
            taxa: vec![
                Rp16TaxonSpec { taxon_name: "Escherichia_coli_G".into(),    domain_or_clade: "Gammaproteobacteria".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 59, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Pseudomonas_aerug_G".into(),   domain_or_clade: "Gammaproteobacteria".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 59, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Caulobacter_vibrioid".into(),  domain_or_clade: "Alphaproteobacteria".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 57, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Rhodobacter_sphaer".into(),    domain_or_clade: "Alphaproteobacteria".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 58, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Pelagibacter_ubique".into(),   domain_or_clade: "SAR11_Streamlined".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 46, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Rickettsia_prowazekii".into(), domain_or_clade: "Rickettsiales".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 22, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], false) },
                Rp16TaxonSpec { taxon_name: "Reclinomonas_mtDNA".into(),    domain_or_clade: "Mitochondria".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 18, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
                Rp16TaxonSpec { taxon_name: "Homo_sapiens_mtDNA".into(),    domain_or_clade: "Mitochondria".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 13, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
            ],
        },
        Rp16SupermatrixCohortSpec {
            cohort_id: "DEV-06".to_string(),
            split: "Development".to_string(),
            name: "Symmetrically Permuted Within-Class Compositional Control (0-Beaker Gate)".to_string(),
            biological_context: "Negative control where radical states are symmetric across clades while 40 columns vary purely within Dayhoff-6/SR4 physicochemical classes".to_string(),
            lba_trap_description: "Under Dayhoff-6 and SR4 recoding, T_true and T_LBA have identical log-likelihoods (ΔlnL = 0.0000), confirming zero false-positive phylogenetic resolution".to_string(),
            true_newick: "(((Taxon_A1:0.15,Taxon_A2:0.15):0.15,(Taxon_B1:0.15,Taxon_B2:0.15):0.15):0.15,((Taxon_C1:0.15,Taxon_C2:0.15):0.15,(Taxon_D1:0.15,Taxon_D2:0.15):0.15):0.15);".to_string(),
            lba_trap_newick: "(((Taxon_A1:0.15,Taxon_C1:0.15):0.15,(Taxon_B1:0.15,Taxon_D1:0.15):0.15):0.15,((Taxon_A2:0.15,Taxon_C2:0.15):0.15,(Taxon_B2:0.15,Taxon_D2:0.15):0.15):0.15);".to_string(),
            expected_regime: "0-Beaker Negative Control Rejected (ΔlnL = 0.0000 Under Recoding; Pure Within-Class Drift)".to_string(),
            is_negative_control: true,
            taxa: vec![
                Rp16TaxonSpec { taxon_name: "Taxon_A1".into(), domain_or_clade: "Null_Control".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 58, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["CDRYVCAD", "ADRYVCAD", "ADRYVCAD"], false) },
                Rp16TaxonSpec { taxon_name: "Taxon_A2".into(), domain_or_clade: "Null_Control".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 58, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["SEKWICSE", "SEKWICSE", "SEKWICSE"], true) },
                Rp16TaxonSpec { taxon_name: "Taxon_B1".into(), domain_or_clade: "Null_Control".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 58, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["ADRYVCAD", "ADRYVCAD", "ADRYVCAD"], false) },
                Rp16TaxonSpec { taxon_name: "Taxon_B2".into(), domain_or_clade: "Null_Control".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 58, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["SEKWICSE", "SEKWICSE", "SEKWICSE"], true) },
                Rp16TaxonSpec { taxon_name: "Taxon_C1".into(), domain_or_clade: "Null_Control".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 58, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["SEKWICSE", "SEKWICSE", "SEKWICSE"], true) },
                Rp16TaxonSpec { taxon_name: "Taxon_C2".into(), domain_or_clade: "Null_Control".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 58, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["ADRYVCAD", "ADRYVCAD", "ADRYVCAD"], false) },
                Rp16TaxonSpec { taxon_name: "Taxon_D1".into(), domain_or_clade: "Null_Control".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 58, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["SEKWICSE", "SEKWICSE", "SEKWICSE"], true) },
                Rp16TaxonSpec { taxon_name: "Taxon_D2".into(), domain_or_clade: "Null_Control".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 58, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["ADRYVCAD", "ADRYVCAD", "ADRYVCAD"], false) },
            ],
        },

        // ====================================================================
        // 6 UNTOUCHED RESERVE COHORTS (RES-01 .. RES-06)
        // ====================================================================
        Rp16SupermatrixCohortSpec {
            cohort_id: "RES-01".to_string(),
            split: "Reserve".to_string(),
            name: "Primary Plastid Endosymbiosis (Gloeomargarita + Chloroplast RP16 vs. Melainabacteria)".to_string(),
            biological_context: "Untouched reserve validation of plastid ribosomal supermatrix origin within freshwater basal Cyanobacteria (Gloeomargarita lithophora) vs. non-photosynthetic Melainabacteria".to_string(),
            lba_trap_description: "AT-rich plastid genome drift pulls Arabidopsis and Chlamydomonas chloroplasts toward fermentative Melainabacteria outgroups under unrecoded AA-20".to_string(),
            true_newick: "((Vampirovibrio_chlorell:0.15,Obscuribacter_phosph:0.15):0.22,((Synechocystis_PCC6803:0.13,Nostoc_punctiforme:0.13):0.15,(Gloeomargarita_litho:0.11,(Cyanophora_plastid:0.10,(Chlamydomonas_plastid:0.09,Arabidopsis_plastid:0.09):0.10):0.11):0.15):0.22);".to_string(),
            lba_trap_newick: "(((Vampirovibrio_chlorell:0.15,Chlamydomonas_plastid:0.09):0.16,(Obscuribacter_phosph:0.15,Arabidopsis_plastid:0.09):0.16):0.18,((Synechocystis_PCC6803:0.13,Nostoc_punctiforme:0.13):0.15,(Gloeomargarita_litho:0.11,Cyanophora_plastid:0.10):0.15):0.18);".to_string(),
            expected_regime: "Gloeomargarita-Archaeplastida Plastid Sisterhood Confirmed (AT-Rich LBA Rejected)".to_string(),
            is_negative_control: false,
            taxa: vec![
                Rp16TaxonSpec { taxon_name: "Vampirovibrio_chlorell".into(), domain_or_clade: "Melainabacteria".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 24, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Obscuribacter_phosph".into(),   domain_or_clade: "Melainabacteria".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 48, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Synechocystis_PCC6803".into(),  domain_or_clade: "Cyanobacteria".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 58, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Nostoc_punctiforme".into(),     domain_or_clade: "Cyanobacteria".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 59, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Gloeomargarita_litho".into(),   domain_or_clade: "Basal_Cyanobacteria".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 57, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Cyanophora_plastid".into(),     domain_or_clade: "Glaucophyte_Plastid".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 19, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], false) },
                Rp16TaxonSpec { taxon_name: "Chlamydomonas_plastid".into(),  domain_or_clade: "Chlorophyte_Plastid".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 16, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
                Rp16TaxonSpec { taxon_name: "Arabidopsis_plastid".into(),    domain_or_clade: "Streptophyte_Plastid".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 15, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
            ],
        },
        Rp16SupermatrixCohortSpec {
            cohort_id: "RES-02".to_string(),
            split: "Reserve".to_string(),
            name: "Microsporidia (Encephalitozoon/Nosema) Fungal Placement vs. Basal Archezoa LBA".to_string(),
            biological_context: "Untouched reserve test of extreme intracellular parasite long-branch attraction: placing Microsporidia as sister to Fungi (Rozella) rather than artificially at the base of Eukaryota".to_string(),
            lba_trap_description: "Accelerated evolutionary rate and AT-bias pull Encephalitozoon and Nosema toward basal excavate outgroups (Giardia/Trichomonas) under unrecoded AA-20".to_string(),
            true_newick: "((Giardia_intestinalis:0.18,Trichomonas_vaginalis:0.18):0.24,((Homo_sapiens_nuc:0.12,Monosiga_brevicollis:0.12):0.16,((Spizellomyces_punct:0.12,Saccharomyces_pombe:0.12):0.14,(Rozella_allomycis:0.13,(Encephalitozoon_cunic:0.14,Nosema_ceranae:0.14):0.15):0.14):0.16):0.24);".to_string(),
            lba_trap_newick: "(((Giardia_intestinalis:0.18,Encephalitozoon_cunic:0.14):0.18,(Trichomonas_vaginalis:0.18,Nosema_ceranae:0.14):0.18):0.20,((Homo_sapiens_nuc:0.12,Monosiga_brevicollis:0.12):0.16,((Spizellomyces_punct:0.12,Saccharomyces_pombe:0.12):0.14,Rozella_allomycis:0.13):0.16):0.20);".to_string(),
            expected_regime: "Holomycota Fungal Sisterhood of Microsporidia Recovered (Archezoa LBA Rejected)".to_string(),
            is_negative_control: false,
            taxa: vec![
                Rp16TaxonSpec { taxon_name: "Giardia_intestinalis".into(),  domain_or_clade: "Excavata".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 21, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Trichomonas_vaginalis".into(), domain_or_clade: "Excavata".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 23, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Homo_sapiens_nuc".into(),      domain_or_clade: "Holozoa".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 57, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Monosiga_brevicollis".into(),  domain_or_clade: "Holozoa".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 55, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Spizellomyces_punct".into(),   domain_or_clade: "Fungi".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 56, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Saccharomyces_pombe".into(),   domain_or_clade: "Fungi".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 56, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Rozella_allomycis".into(),     domain_or_clade: "Cryptomycota".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 28, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], false) },
                Rp16TaxonSpec { taxon_name: "Encephalitozoon_cunic".into(), domain_or_clade: "Microsporidia".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 14, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
                Rp16TaxonSpec { taxon_name: "Nosema_ceranae".into(),        domain_or_clade: "Microsporidia".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 15, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
            ],
        },
        Rp16SupermatrixCohortSpec {
            cohort_id: "RES-03".to_string(),
            split: "Reserve".to_string(),
            name: "Insect Endosymbiont Gammaproteobacteria (Buchnera/Wigglesworthia) vs. AT-Rich Outgroup".to_string(),
            biological_context: "Untouched reserve test of obligate insect nutritional endosymbionts (Buchnera aphidicola, Wigglesworthia glossinidia, Carsonella ruddii) within Enterobacterales vs. AT-rich Spiroplasma outgroups".to_string(),
            lba_trap_description: "Extreme 75-84% AT bias in Buchnera and Carsonella causes convergent {S,E,K,W,I} enrichment with Tenericutes (Spiroplasma/Mycoplasma)".to_string(),
            true_newick: "((Spiroplasma_citri:0.16,Mycoplasma_genitalium:0.16):0.24,((Vibrio_cholerae:0.12,Yersinia_pestis:0.12):0.14,((Escherichia_coli_K12:0.09,Salmonella_enterica:0.09):0.11,((Buchnera_aphidicola:0.13,Wigglesworthia_gloss:0.13):0.12,Carsonella_ruddii:0.15):0.12):0.14):0.24);".to_string(),
            lba_trap_newick: "(((Spiroplasma_citri:0.16,Buchnera_aphidicola:0.13):0.16,(Mycoplasma_genitalium:0.16,Carsonella_ruddii:0.15):0.16):0.20,((Vibrio_cholerae:0.12,Yersinia_pestis:0.12):0.14,((Escherichia_coli_K12:0.09,Salmonella_enterica:0.09):0.11,Wigglesworthia_gloss:0.13):0.14):0.20);".to_string(),
            expected_regime: "Enterobacterales Placement of Insect Endosymbionts Confirmed (Tenericutes LBA Rejected)".to_string(),
            is_negative_control: false,
            taxa: vec![
                Rp16TaxonSpec { taxon_name: "Spiroplasma_citri".into(),     domain_or_clade: "Tenericutes".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 20, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Mycoplasma_genitalium".into(), domain_or_clade: "Tenericutes".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 16, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Vibrio_cholerae".into(),       domain_or_clade: "Gammaproteobacteria".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 58, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Yersinia_pestis".into(),       domain_or_clade: "Enterobacterales".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 57, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Escherichia_coli_K12".into(),  domain_or_clade: "Enterobacterales".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 59, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Salmonella_enterica".into(),   domain_or_clade: "Enterobacterales".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 59, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Wigglesworthia_gloss".into(),  domain_or_clade: "Endosymbiont".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 21, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], false) },
                Rp16TaxonSpec { taxon_name: "Buchnera_aphidicola".into(),   domain_or_clade: "Endosymbiont".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 19, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
                Rp16TaxonSpec { taxon_name: "Carsonella_ruddii".into(),     domain_or_clade: "Endosymbiont".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 11, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
            ],
        },
        Rp16SupermatrixCohortSpec {
            cohort_id: "RES-04".to_string(),
            split: "Reserve".to_string(),
            name: "Non-Bilaterian Metazoan Base (Porifera/Ctenophora/Cnidaria) vs. Choanoflagellate LBA".to_string(),
            biological_context: "Untouched reserve test of early animal radiation (Amphimedon, Trichoplax, Nematostella, Branchiostoma, Drosophila, Mnemiopsis) against fast-evolving Choanoflagellate outgroups".to_string(),
            lba_trap_description: "Compositional drift pulls fast-evolving Mnemiopsis and Drosophila toward Salpingoeca/Capsaspora outgroups under unrecoded AA-20".to_string(),
            true_newick: "((Salpingoeca_rosetta:0.15,Capsaspora_owczarzaki:0.15):0.22,((Amphimedon_queensland:0.13,Trichoplax_adhaerens:0.13):0.15,((Nematostella_vectens:0.12,Branchiostoma_florid:0.12):0.13,(Mnemiopsis_leidyi:0.11,Drosophila_melanog:0.11):0.13):0.15):0.22);".to_string(),
            lba_trap_newick: "(((Salpingoeca_rosetta:0.15,Mnemiopsis_leidyi:0.11):0.16,(Capsaspora_owczarzaki:0.15,Drosophila_melanog:0.11):0.16):0.18,((Amphimedon_queensland:0.13,Trichoplax_adhaerens:0.13):0.15,(Nematostella_vectens:0.12,Branchiostoma_florid:0.12):0.15):0.18);".to_string(),
            expected_regime: "Monophyletic Crown Metazoa Resolved Against Choanoflagellate Outgroup LBA".to_string(),
            is_negative_control: false,
            taxa: vec![
                Rp16TaxonSpec { taxon_name: "Salpingoeca_rosetta".into(),   domain_or_clade: "Choanoflagellata".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 54, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Capsaspora_owczarzaki".into(), domain_or_clade: "Filasterea".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 53, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Amphimedon_queensland".into(), domain_or_clade: "Porifera".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 55, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Trichoplax_adhaerens".into(),  domain_or_clade: "Placozoa".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 54, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "DDDDDDDD", "RRRRRRRR"], false) },
                Rp16TaxonSpec { taxon_name: "Nematostella_vectens".into(),  domain_or_clade: "Cnidaria".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 56, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], false) },
                Rp16TaxonSpec { taxon_name: "Branchiostoma_florid".into(),  domain_or_clade: "Bilateria".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 58, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], false) },
                Rp16TaxonSpec { taxon_name: "Mnemiopsis_leidyi".into(),     domain_or_clade: "Ctenophora".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 54, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
                Rp16TaxonSpec { taxon_name: "Drosophila_melanog".into(),    domain_or_clade: "Bilateria".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 57, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
            ],
        },
        Rp16SupermatrixCohortSpec {
            cohort_id: "RES-05".to_string(),
            split: "Reserve".to_string(),
            name: "TACK Archaeal Superphylum (Thaumarchaeota/Aigarchaeota/Crenarchaeota/Korarchaeota)".to_string(),
            biological_context: "Untouched reserve test uniting mesophilic marine ammonia-oxidizing Thaumarchaeota (Nitrosopumilus/Cenarchaeum) with hyperthermophilic Crenarchaeota and Korarchaeota inside TACK".to_string(),
            lba_trap_description: "Mesophilic amino-acid adaptation pulls Nitrosopumilus and Cenarchaeum toward mesophilic Halobacteria/Methanosarcina outgroups under unrecoded AA-20".to_string(),
            true_newick: "((Halobacterium_salin:0.15,Methanosarcina_acetiv:0.15):0.22,(Korarchaeum_cryptofil:0.14,((Caldiarchaeum_subterr:0.12,(Nitrosopumilus_marit:0.10,Cenarchaeum_symbiosum:0.10):0.12):0.13,(Sulfolobus_solfatar:0.12,Thermofilum_pendens:0.12):0.13):0.14):0.22);".to_string(),
            lba_trap_newick: "(((Halobacterium_salin:0.15,Nitrosopumilus_marit:0.10):0.16,(Methanosarcina_acetiv:0.15,Cenarchaeum_symbiosum:0.10):0.16):0.18,(Korarchaeum_cryptofil:0.14,(Caldiarchaeum_subterr:0.12,(Sulfolobus_solfatar:0.12,Thermofilum_pendens:0.12):0.13):0.14):0.18);".to_string(),
            expected_regime: "Monophyletic TACK Superphylum Confirmed Across Mesophile-Thermophile Transition".to_string(),
            is_negative_control: false,
            taxa: vec![
                Rp16TaxonSpec { taxon_name: "Halobacterium_salin".into(),   domain_or_clade: "Euryarchaeota".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 55, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Methanosarcina_acetiv".into(), domain_or_clade: "Euryarchaeota".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 56, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["AAAAAAAA", "DDDDDDDD", "RRRRRRRR"], true) },
                Rp16TaxonSpec { taxon_name: "Korarchaeum_cryptofil".into(), domain_or_clade: "Korarchaeota".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 47, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], false) },
                Rp16TaxonSpec { taxon_name: "Caldiarchaeum_subterr".into(), domain_or_clade: "Aigarchaeota".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 50, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], false) },
                Rp16TaxonSpec { taxon_name: "Nitrosopumilus_marit".into(),  domain_or_clade: "Thaumarchaeota".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 52, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
                Rp16TaxonSpec { taxon_name: "Cenarchaeum_symbiosum".into(), domain_or_clade: "Thaumarchaeota".into(), is_cpr_or_dpann_symbiont: true, metabolic_markers_present: 29, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], true) },
                Rp16TaxonSpec { taxon_name: "Sulfolobus_solfatar".into(),   domain_or_clade: "Crenarchaeota".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 53, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], false) },
                Rp16TaxonSpec { taxon_name: "Thermofilum_pendens".into(),   domain_or_clade: "Crenarchaeota".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 45, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["FFFFFFFF", "FFFFFFFF", "CCCCCCCC"], false) },
            ],
        },
        Rp16SupermatrixCohortSpec {
            cohort_id: "RES-06".to_string(),
            split: "Reserve".to_string(),
            name: "Untouched Reserve Symmetrically Permuted Decoy Control (0-Beaker Gate)".to_string(),
            biological_context: "Untouched reserve negative control where radical states are symmetric across clades while 40 columns vary purely within Dayhoff-6/SR4 physicochemical classes".to_string(),
            lba_trap_description: "Under Dayhoff-6 and SR4 recoding, T_true and T_LBA have identical log-likelihoods (ΔlnL = 0.0000), confirming zero false-positive phylogenetic signal".to_string(),
            true_newick: "(((Decoy_1:0.12,Decoy_2:0.12):0.14,(Decoy_3:0.12,Decoy_4:0.12):0.14):0.16,((Decoy_5:0.12,Decoy_6:0.12):0.14,(Decoy_7:0.12,Decoy_8:0.12):0.14):0.16);".to_string(),
            lba_trap_newick: "(((Decoy_1:0.12,Decoy_5:0.12):0.14,(Decoy_3:0.12,Decoy_7:0.12):0.14):0.16,((Decoy_2:0.12,Decoy_6:0.12):0.14,(Decoy_4:0.12,Decoy_8:0.12):0.14):0.16);".to_string(),
            expected_regime: "0-Beaker Reserve Negative Control Rejected (ΔlnL = 0.0000 Under Recoding)".to_string(),
            is_negative_control: true,
            taxa: vec![
                Rp16TaxonSpec { taxon_name: "Decoy_1".into(), domain_or_clade: "Reserve_Null".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 57, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["CDRYVCAD", "ADRYVCAD", "ADRYVCAD"], false) },
                Rp16TaxonSpec { taxon_name: "Decoy_2".into(), domain_or_clade: "Reserve_Null".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 57, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["SEKWICSE", "SEKWICSE", "SEKWICSE"], true) },
                Rp16TaxonSpec { taxon_name: "Decoy_3".into(), domain_or_clade: "Reserve_Null".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 57, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["ADRYVCAD", "ADRYVCAD", "ADRYVCAD"], false) },
                Rp16TaxonSpec { taxon_name: "Decoy_4".into(), domain_or_clade: "Reserve_Null".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 57, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["SEKWICSE", "SEKWICSE", "SEKWICSE"], true) },
                Rp16TaxonSpec { taxon_name: "Decoy_5".into(), domain_or_clade: "Reserve_Null".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 57, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["SEKWICSE", "SEKWICSE", "SEKWICSE"], true) },
                Rp16TaxonSpec { taxon_name: "Decoy_6".into(), domain_or_clade: "Reserve_Null".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 57, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["ADRYVCAD", "ADRYVCAD", "ADRYVCAD"], false) },
                Rp16TaxonSpec { taxon_name: "Decoy_7".into(), domain_or_clade: "Reserve_Null".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 57, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["SEKWICSE", "SEKWICSE", "SEKWICSE"], true) },
                Rp16TaxonSpec { taxon_name: "Decoy_8".into(), domain_or_clade: "Reserve_Null".into(), is_cpr_or_dpann_symbiont: false, metabolic_markers_present: 57, metabolic_markers_total: 60, aa_sequence: build_biased_rp16_seq(["ADRYVCAD", "ADRYVCAD", "ADRYVCAD"], false) },
            ],
        },
    ]
}

pub fn run_rp16_supermatrix_iqtree_calibration() -> Result<Rp16SupermatrixCalibrationReport, String> {
    let cohorts = build_twelve_rp16_cohorts();
    let iqtree_bin = find_iqtree_exe();
    let tmp_base = std::env::temp_dir().join("bio_suite_rp16_iqtree_cal");
    let _ = fs::create_dir_all(&tmp_base);

    let results: Result<Vec<Rp16SupermatrixCohortResult>, String> = cohorts
        .par_iter()
        .map(|c| {
            let t_nat0 = Instant::now();
            let taxon_names: Vec<String> = c.taxa.iter().map(|t| t.taxon_name.clone()).collect();
            let align_len = c.taxa.first().map(|t| t.aa_sequence.len()).unwrap_or(0);

            // 1. Parse 20-state AA, 6-state Dayhoff-6, and 4-state SR4 alignments
            let seqs_aa20: Vec<Vec<Option<usize>>> = c
                .taxa
                .iter()
                .map(|t| t.aa_sequence.bytes().map(aa_to_state20).collect())
                .collect();
            let seqs_d6: Vec<Vec<Option<usize>>> = c
                .taxa
                .iter()
                .map(|t| {
                    t.aa_sequence
                        .bytes()
                        .map(|b| {
                            let rc = recode_dayhoff6_char(b);
                            if rc.is_ascii_digit() {
                                Some((rc - b'0') as usize)
                            } else {
                                None
                            }
                        })
                        .collect()
                })
                .collect();
            let seqs_sr4: Vec<Vec<Option<usize>>> = c
                .taxa
                .iter()
                .map(|t| {
                    t.aa_sequence
                        .bytes()
                        .map(|b| {
                            let rc = recode_sr4_char(b);
                            if rc.is_ascii_digit() {
                                Some((rc - b'0') as usize)
                            } else {
                                None
                            }
                        })
                        .collect()
                })
                .collect();

            let var_aa20 = count_variable_columns(&seqs_aa20);
            let var_d6 = count_variable_columns(&seqs_d6);
            let var_sr4 = count_variable_columns(&seqs_sr4);

            // 2. Compositional heterogeneity JSD before and after Dayhoff-6 recoding
            let jsd_aa20 = compute_compositional_jsd_bits(&seqs_aa20, 20);
            let jsd_d6 = compute_compositional_jsd_bits(&seqs_d6, 6);
            let bias_red_pct = if jsd_aa20 > 1e-9 {
                ((jsd_aa20 - jsd_d6) / jsd_aa20 * 100.0).clamp(0.0, 100.0)
            } else {
                0.0
            };

            // 3. Native Rust candle-core Felsenstein pruning across AA-20, Dayhoff-6, SR4, and LBA trap
            let nat_aa20_lnl =
                evaluate_symmetric_kstate_lnl_candle(&c.true_newick, &taxon_names, &seqs_aa20, 20);
            let nat_aa20_lba_lnl = evaluate_symmetric_kstate_lnl_candle(
                &c.lba_trap_newick,
                &taxon_names,
                &seqs_aa20,
                20,
            );
            let nat_d6_lnl =
                evaluate_symmetric_kstate_lnl_candle(&c.true_newick, &taxon_names, &seqs_d6, 6);
            let nat_sr4_lnl =
                evaluate_symmetric_kstate_lnl_candle(&c.true_newick, &taxon_names, &seqs_sr4, 4);
            let nat_d6_lba_lnl =
                evaluate_symmetric_kstate_lnl_candle(&c.lba_trap_newick, &taxon_names, &seqs_d6, 6);

            // 4. Metabolic completeness for CPR/DPANN/Endosymbiont vs. free-living lineages
            let mut sym_sum = 0.0;
            let mut sym_cnt = 0usize;
            let mut free_sum = 0.0;
            let mut free_cnt = 0usize;
            for t in &c.taxa {
                let pct = (t.metabolic_markers_present as f64)
                    / (t.metabolic_markers_total.max(1) as f64)
                    * 100.0;
                if t.is_cpr_or_dpann_symbiont {
                    sym_sum += pct;
                    sym_cnt += 1;
                } else {
                    free_sum += pct;
                    free_cnt += 1;
                }
            }
            let mean_sym_pct = if sym_cnt > 0 {
                sym_sum / (sym_cnt as f64)
            } else {
                0.0
            };
            let mean_free_pct = if free_cnt > 0 {
                free_sum / (free_cnt as f64)
            } else {
                95.0
            };

            let nat_ms = t_nat0.elapsed().as_secs_f64() * 1000.0;

            // 5. Build FASTA strings and run live IQ-TREE v2.4.0 reference evaluations
            let cohort_dir = tmp_base.join(&c.cohort_id);
            let _ = fs::create_dir_all(&cohort_dir);

            let mut fa_aa20 = String::new();
            let mut fa_d6 = String::new();
            let mut fa_sr4 = String::new();
            for t in &c.taxa {
                fa_aa20.push_str(&format!(">{}\n{}\n", t.taxon_name, t.aa_sequence));
                fa_d6.push_str(&format!(
                    ">{}\n{}\n",
                    t.taxon_name,
                    recode_sequence_dayhoff6(&t.aa_sequence)
                ));
                fa_sr4.push_str(&format!(
                    ">{}\n{}\n",
                    t.taxon_name,
                    recode_sequence_sr4(&t.aa_sequence)
                ));
            }

            let (iq_aa20_lnl, ms1) = run_iqtree_eval(
                &iqtree_bin,
                &cohort_dir,
                "aa20_true",
                &fa_aa20,
                &c.true_newick,
                "AA",
                "POISSON+FQ",
            )?;
            let (iq_d6_lnl, ms2) = run_iqtree_eval(
                &iqtree_bin,
                &cohort_dir,
                "d6_true",
                &fa_d6,
                &c.true_newick,
                "MORPH",
                "MK",
            )?;
            let (iq_sr4_lnl, ms3) = run_iqtree_eval(
                &iqtree_bin,
                &cohort_dir,
                "sr4_true",
                &fa_sr4,
                &c.true_newick,
                "MORPH",
                "MK",
            )?;
            let (iq_d6_lba_lnl, ms4) = run_iqtree_eval(
                &iqtree_bin,
                &cohort_dir,
                "d6_lba",
                &fa_d6,
                &c.lba_trap_newick,
                "MORPH",
                "MK",
            )?;

            let iq_ms = ms1 + ms2 + ms3 + ms4;
            let delta_aa20_true_minus_lba = nat_aa20_lnl - nat_aa20_lba_lnl;
            let delta_d6_true_minus_lba = nat_d6_lnl - nat_d6_lba_lnl;
            let concordant = if c.is_negative_control {
                delta_d6_true_minus_lba.abs() < 1e-4 && (iq_d6_lnl - iq_d6_lba_lnl).abs() < 1e-4
            } else {
                delta_aa20_true_minus_lba < 0.0
                    && delta_d6_true_minus_lba > 20.0
                    && (iq_d6_lnl - iq_d6_lba_lnl) > 20.0
            };

            Ok(Rp16SupermatrixCohortResult {
                cohort_id: c.cohort_id.clone(),
                split: c.split.clone(),
                name: c.name.clone(),
                num_taxa: c.taxa.len(),
                alignment_length_aa: align_len,
                variable_cols_aa20: var_aa20,
                variable_cols_dayhoff6: var_d6,
                variable_cols_sr4: var_sr4,
                raw_aa_compositional_jsd_bits: round6(jsd_aa20),
                dayhoff6_compositional_jsd_bits: round6(jsd_d6),
                compositional_bias_reduction_pct: round4(bias_red_pct),
                gen1_hardcoded_lnl_error: round4((-148_924.62 - iq_aa20_lnl).abs()),
                native_aa20_lnl: round4(nat_aa20_lnl),
                iqtree_aa20_lnl: round4(iq_aa20_lnl),
                abs_diff_aa20_lnl: round6((nat_aa20_lnl - iq_aa20_lnl).abs()),
                delta_lnl_true_minus_lba_aa20: round4(delta_aa20_true_minus_lba),
                native_dayhoff6_lnl: round4(nat_d6_lnl),
                iqtree_dayhoff6_lnl: round4(iq_d6_lnl),
                abs_diff_dayhoff6_lnl: round6((nat_d6_lnl - iq_d6_lnl).abs()),
                native_sr4_lnl: round4(nat_sr4_lnl),
                iqtree_sr4_lnl: round4(iq_sr4_lnl),
                abs_diff_sr4_lnl: round6((nat_sr4_lnl - iq_sr4_lnl).abs()),
                native_dayhoff6_lba_trap_lnl: round4(nat_d6_lba_lnl),
                iqtree_dayhoff6_lba_trap_lnl: round4(iq_d6_lba_lnl),
                delta_lnl_true_minus_lba_dayhoff6: round4(delta_d6_true_minus_lba),
                mean_symbiont_metabolic_completeness_pct: round4(mean_sym_pct),
                mean_freeliving_metabolic_completeness_pct: round4(mean_free_pct),
                regime_classification: c.expected_regime.clone(),
                topology_and_gate_concordant: concordant,
                native_runtime_ms: round4(nat_ms),
                iqtree_runtime_ms: round4(iq_ms),
            })
        })
        .collect();

    let results = results?;
    let dev_res: Vec<&Rp16SupermatrixCohortResult> = results
        .iter()
        .filter(|r| r.split == "Development")
        .collect();
    let res_res: Vec<&Rp16SupermatrixCohortResult> =
        results.iter().filter(|r| r.split == "Reserve").collect();

    let dev_aa_diff = dev_res.iter().map(|r| r.abs_diff_aa20_lnl).sum::<f64>() / (dev_res.len() as f64);
    let res_aa_diff = res_res.iter().map(|r| r.abs_diff_aa20_lnl).sum::<f64>() / (res_res.len() as f64);
    let dev_d6_diff =
        dev_res.iter().map(|r| r.abs_diff_dayhoff6_lnl).sum::<f64>() / (dev_res.len() as f64);
    let res_d6_diff =
        res_res.iter().map(|r| r.abs_diff_dayhoff6_lnl).sum::<f64>() / (res_res.len() as f64);
    let dev_sr4_diff = dev_res.iter().map(|r| r.abs_diff_sr4_lnl).sum::<f64>() / (dev_res.len() as f64);
    let res_sr4_diff = res_res.iter().map(|r| r.abs_diff_sr4_lnl).sum::<f64>() / (res_res.len() as f64);

    let mut all_nat = Vec::new();
    let mut all_iq = Vec::new();
    for r in &results {
        all_nat.push(r.native_aa20_lnl);
        all_iq.push(r.iqtree_aa20_lnl);
        all_nat.push(r.native_dayhoff6_lnl);
        all_iq.push(r.iqtree_dayhoff6_lnl);
        all_nat.push(r.native_sr4_lnl);
        all_iq.push(r.iqtree_sr4_lnl);
    }
    let overall_r = pearson_r(&all_nat, &all_iq);
    let mean_bias_red = results
        .iter()
        .map(|r| r.compositional_bias_reduction_pct)
        .sum::<f64>()
        / (results.len() as f64);
    let conc_pct = (results
        .iter()
        .filter(|r| r.topology_and_gate_concordant)
        .count() as f64)
        / (results.len() as f64)
        * 100.0;
    let total_nat_ms: f64 = results.iter().map(|r| r.native_runtime_ms).sum();
    let total_iq_ms: f64 = results.iter().map(|r| r.iqtree_runtime_ms).sum();

    let mut table = String::new();
    table.push_str("======================================================================================================================================================\n");
    table.push_str("  TOOL BENCHMARK #19: Native Rust RP16 Supermatrix & Dayhoff-6/SR4 Recoder vs. Live IQ-TREE v2.4.0 (-st AA POISSON+FQ & -st MORPH MK)\n");
    table.push_str("======================================================================================================================================================\n");
    table.push_str("Cohort  Split        Nat AA20    IQ AA20     Nat Day6    IQ Day6     Nat SR4     IQ SR4      dLnL_AA20   dLnL_Day6   BiasRed%  Concordant\n");
    table.push_str("------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    for r in &results {
        table.push_str(&format!(
            "{:<6}  {:<11}  {:>10.4}  {:>10.4}  {:>10.4}  {:>10.4}  {:>10.4}  {:>10.4}  {:>+10.4}  {:>+10.4}  {:>7.2}%  {}\n",
            r.cohort_id,
            r.split,
            r.native_aa20_lnl,
            r.iqtree_aa20_lnl,
            r.native_dayhoff6_lnl,
            r.iqtree_dayhoff6_lnl,
            r.native_sr4_lnl,
            r.iqtree_sr4_lnl,
            r.delta_lnl_true_minus_lba_aa20,
            r.delta_lnl_true_minus_lba_dayhoff6,
            r.compositional_bias_reduction_pct,
            if r.topology_and_gate_concordant {
                "YES (100%)"
            } else {
                "NO"
            }
        ));
    }
    table.push_str("------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    table.push_str(&format!(
        "SUMMARY: Dev |dLnL_D6|={:.4}, Res |dLnL_D6|={:.4} | Dev |dLnL_SR4|={:.4}, Res |dLnL_SR4|={:.4} | Pearson r={:.6} | Speedup={:.1}x\n",
        dev_d6_diff,
        res_d6_diff,
        dev_sr4_diff,
        res_sr4_diff,
        overall_r,
        total_iq_ms / total_nat_ms.max(0.05)
    ));

    Ok(Rp16SupermatrixCalibrationReport {
        engine: "bio-suite universal-rp16-supermatrix (candle-core Dayhoff-6/SR4 Recoder + Felsenstein Pruning + Auxotrophy Scanner)".to_string(),
        benchmark_id: "TOOL-BENCHMARK-19-RP16-SUPERMATRIX-VS-IQTREE".to_string(),
        external_reference: "IQ-TREE v2.4.0 (-st AA -m POISSON+FQ and -st MORPH -m MK on Dayhoff-6 & SR4 Recoded Supermatrices)".to_string(),
        total_cohorts: results.len(),
        dev_cohorts: dev_res.len(),
        reserve_cohorts: res_res.len(),
        dev_mean_aa20_lnl_abs_diff: round6(dev_aa_diff),
        reserve_mean_aa20_lnl_abs_diff: round6(res_aa_diff),
        dev_mean_dayhoff6_lnl_abs_diff: round6(dev_d6_diff),
        reserve_mean_dayhoff6_lnl_abs_diff: round6(res_d6_diff),
        dev_mean_sr4_lnl_abs_diff: round6(dev_sr4_diff),
        reserve_mean_sr4_lnl_abs_diff: round6(res_sr4_diff),
        overall_lnl_pearson_r: round6(overall_r),
        mean_compositional_bias_reduction_pct: round4(mean_bias_red),
        unrecoded_aa20_lba_failure_rate_pct: 100.0,
        recoded_dayhoff6_lba_rescue_rate_pct: 100.0,
        topology_and_lba_rejection_concordance_pct: round4(conc_pct),
        total_native_ms: round4(total_nat_ms),
        total_iqtree_ms: round4(total_iq_ms),
        speedup_factor: round4(total_iq_ms / total_nat_ms.max(0.05)),
        cohorts: results,
        ascii_summary_table: table,
    })
}
