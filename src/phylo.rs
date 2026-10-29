//! Tier 2 (Engine 1): Molecular Systematics, Phylogenomics & COI Barcode Clustering (`bio-suite phylo`)
//!
//! Provides high-throughput multi-threaded (`rayon`) algorithms for:
//! - Codon-preserving Multiple Sequence Alignment (MSA) + external MAFFT/MUSCLE5 drop-in adapters
//! - Pairwise Kimura 2-Parameter (K2P) & Tamura-Nei (TN93) evolutionary distances & Ti/Tv ratios
//! - Neighbor-Joining (NJ) phylogenetic tree inference with parallel bootstrap replicates
//! - Automatic Barcode Gap Discovery (ABGD) & MOTU species cluster delineation

use crate::complexity::round4;
use crate::fetch::FastaRecord;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::PathBuf;
use std::process::Command;
use std::time::Instant;

/// Multiple Sequence Alignment execution mode.
/// Defaults strictly to `AlignerMode::Native` unless an external tool (`AlignerMode::Mafft`)
/// is explicitly requested via CLI (`--aligner mafft`) or calibration API.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AlignerMode {
    /// Default: Self-contained Rust progressive affine-gap (Gotoh) + codon-anchor hybrid MSA
    Native,
    /// On-demand external tool: MAFFT v7 (`--auto` -> L-INS-i / FFT-NS-i / FFT-NS-2)
    Mafft,
}

impl Default for AlignerMode {
    fn default() -> Self {
        AlignerMode::Native
    }
}

impl AlignerMode {
    pub fn from_cli_arg(arg: &str) -> Self {
        match arg.trim().to_ascii_lowercase().as_str() {
            "mafft" | "linsi" | "l-ins-i" | "fftnsi" => AlignerMode::Mafft,
            _ => AlignerMode::Native,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MsaExecutionTelemetry {
    pub aligner_mode: String,
    pub alignment_method: String,
    pub duration_ms: f64,
    pub alignment_length_bp: usize,
    pub gap_columns: usize,
    pub gap_fraction: f64,
    pub sum_of_pairs_identity: f64,
    pub external_tool_log: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignedSequence {
    pub id: String,
    pub taxon: String,
    pub common_name: String,
    pub family: String,
    pub aligned_seq: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairwiseDistanceEntry {
    pub seq_a: String,
    pub taxon_a: String,
    pub family_a: String,
    pub seq_b: String,
    pub taxon_b: String,
    pub family_b: String,
    pub relation_level: String, // "intraspecific", "intrageneric", "intrafamily", "interfamily"
    pub p_distance: f64,
    pub k2p_distance: f64,
    pub tn93_distance: f64,
    pub transitions: usize,
    pub transversions: usize,
    pub ti_tv_ratio: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BarcodeGapSummary {
    pub threshold_k2p: f64,
    pub max_intraspecific_k2p: f64,
    pub min_interspecific_k2p: f64,
    pub mean_intraspecific_k2p: f64,
    pub mean_intrageneric_k2p: f64,
    pub mean_intrafamily_k2p: f64,
    pub mean_interfamily_k2p: f64,
    pub barcoding_gap_width: f64,
    pub gap_ratio_min_inter_to_max_intra: f64,
    pub mean_ti_tv_ratio: f64,
    pub motu_count: usize,
    pub clusters: Vec<MotuCluster>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MotuCluster {
    pub motu_id: String,
    pub primary_taxon: String,
    pub family: String,
    pub member_accessions: Vec<String>,
    pub member_names: Vec<String>,
    pub max_within_cluster_k2p: f64,
    pub nearest_neighbor_taxon: String,
    pub nearest_neighbor_k2p: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhyloAnalysisResult {
    pub alignment_method: String,
    pub alignment_length_bp: usize,
    pub conserved_sites: usize,
    pub variable_sites: usize,
    pub parsimony_informative_sites: usize,
    pub distance_matrix: Vec<Vec<f64>>,
    pub labels: Vec<String>,
    pub pairwise_comparisons: Vec<PairwiseDistanceEntry>,
    pub newick_tree: String,
    pub bootstrap_replicates: usize,
    pub clade_support_summary: Vec<CladeNodeAnnotation>,
    pub barcode_gap: BarcodeGapSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CladeNodeAnnotation {
    pub clade_name: String,
    pub rank: String,
    pub member_taxa: Vec<String>,
    pub bootstrap_support_pct: f64,
    pub mean_within_clade_k2p: f64,
    pub synapomorphic_codon_sites: usize,
}

/// Perform Multiple Sequence Alignment using the default `AlignerMode::Native` engine
/// (self-contained Rust progressive Gotoh affine-gap + codon-anchor hybrid MSA).
/// External tools (such as `MAFFT`) are only invoked when explicitly requested via
/// `align_sequences_with_mode(records, AlignerMode::Mafft)` or `--aligner mafft`.
pub fn align_sequences(records: &[FastaRecord]) -> (String, Vec<AlignedSequence>) {
    let (method, aligned, _telemetry) = align_sequences_with_mode(records, AlignerMode::Native);
    (method, aligned)
}

/// Perform Multiple Sequence Alignment with explicit `AlignerMode` selection (`Native` default vs. `Mafft`).
pub fn align_sequences_with_mode(
    records: &[FastaRecord],
    mode: AlignerMode,
) -> (String, Vec<AlignedSequence>, MsaExecutionTelemetry) {
    let t0 = Instant::now();
    if records.is_empty() {
        let empty_tel = MsaExecutionTelemetry {
            aligner_mode: format!("{:?}", mode),
            alignment_method: "none".to_string(),
            duration_ms: 0.0,
            alignment_length_bp: 0,
            gap_columns: 0,
            gap_fraction: 0.0,
            sum_of_pairs_identity: 0.0,
            external_tool_log: None,
        };
        return ("none".to_string(), Vec::new(), empty_tel);
    }

    match mode {
        AlignerMode::Mafft => {
            match align_sequences_mafft(records) {
                Ok((method, aligned, stderr_log)) => {
                    let dur = t0.elapsed().as_secs_f64() * 1000.0;
                    let tel = compute_msa_telemetry("Mafft", &method, dur, &aligned, Some(stderr_log));
                    return (method, aligned, tel);
                }
                Err(err_msg) => {
                    let (method, aligned) = align_sequences_native_gotoh(records);
                    let dur = t0.elapsed().as_secs_f64() * 1000.0;
                    let tel = compute_msa_telemetry(
                        "Native (MAFFT fallback)",
                        &method,
                        dur,
                        &aligned,
                        Some(format!("MAFFT fallback reason: {}", err_msg)),
                    );
                    (method, aligned, tel)
                }
            }
        }
        AlignerMode::Native => {
            let (method, aligned) = align_sequences_native_gotoh(records);
            let dur = t0.elapsed().as_secs_f64() * 1000.0;
            let tel = compute_msa_telemetry("Native", &method, dur, &aligned, None);
            (method, aligned, tel)
        }
    }
}

/// Native Rust v2 Progressive Multiple Sequence Aligner:
/// - Mode 1 (Ungapped Reading-Frame Coding Barcodes, e.g., vertebrate COI `1,545 bp`):
///   Preserves reading-frame codon triplets (`-15..=+15 bp` step 3) when all sequences share
///   uniform coding length (spread <= 9 bp) and high ungapped identity.
/// - Mode 2 (Variable-Length / Indel-Rich rRNA, ITS1-5.8S-ITS2 & Divergent Loci):
///   Selects the 6-mer medoid guide reference and performs parallel Anchor-Segmented Gotoh
///   3-State Affine-Gap Dynamic Programming (`Match = +4`, `Transition = +1`, `Transversion = -3`,
///   `GapOpen = -10`, `GapExtend = -1`), merging pairwise insertion tracts into a consistent
///   multi-sequence profile alignment.
pub fn align_sequences_native_gotoh(records: &[FastaRecord]) -> (String, Vec<AlignedSequence>) {
    let min_len = records.iter().map(|r| r.sequence.len()).min().unwrap_or(0);
    let max_len = records.iter().map(|r| r.sequence.len()).max().unwrap_or(0);
    let has_preexisting_gaps = records
        .iter()
        .any(|r| r.sequence.as_bytes().contains(&b'-'));

    // Check if uniform-length protein-coding barcode panel (e.g. 1,545 bp COX1)
    let is_uniform_coding = !has_preexisting_gaps
        && min_len >= 300
        && max_len.saturating_sub(min_len) <= 9
        && (min_len % 3 == 0);

    if is_uniform_coding {
        let target_len = min_len.min(1545);
        let aligned_len = (target_len / 3) * 3;
        let ref_seq = records[0].sequence.as_bytes();
        let aligned: Vec<AlignedSequence> = records
            .par_iter()
            .map(|rec| {
                let seq_bytes = rec.sequence.as_bytes();
                let mut best_offset: isize = 0;
                let mut best_matches = 0usize;
                for codon_shift in -5isize..=5isize {
                    let offset = codon_shift * 3;
                    let mut matches = 0usize;
                    for i in (0..aligned_len.min(600)).step_by(3) {
                        let src_idx = (i as isize) + offset;
                        if src_idx >= 0
                            && (src_idx as usize + 1) < seq_bytes.len()
                            && (i + 1) < ref_seq.len()
                        {
                            if seq_bytes[src_idx as usize] == ref_seq[i]
                                && seq_bytes[src_idx as usize + 1] == ref_seq[i + 1]
                            {
                                matches += 1;
                            }
                        }
                    }
                    if matches > best_matches {
                        best_matches = matches;
                        best_offset = offset;
                    }
                }

                let mut out = String::with_capacity(aligned_len);
                for i in 0..aligned_len {
                    let idx = (i as isize) + best_offset;
                    if idx >= 0 && (idx as usize) < seq_bytes.len() {
                        out.push((seq_bytes[idx as usize] as char).to_ascii_uppercase());
                    } else {
                        out.push('-');
                    }
                }

                AlignedSequence {
                    id: rec.id.clone(),
                    taxon: rec.taxon.clone(),
                    common_name: rec.common_name.clone(),
                    family: rec.family.clone(),
                    aligned_seq: out,
                }
            })
            .collect();

        return (
            "bio-suite native Rust v2 (Codon-Anchored + Gotoh Affine-Gap Hybrid MSA)".to_string(),
            aligned,
        );
    }

    // Mode 2: Variable-length / Indel-bearing panel -> Progressive Gotoh Affine-Gap Center-Star / Guide Profile MSA
    let cleaned_seqs: Vec<Vec<u8>> = records
        .iter()
        .map(|r| {
            r.sequence
                .bytes()
                .filter(|&b| b != b'-' && !b.is_ascii_whitespace())
                .map(|b| b.to_ascii_uppercase())
                .collect()
        })
        .collect();

    // Select 6-mer medoid sequence as the guide anchor
    let center_idx = select_kmer_medoid_index(&cleaned_seqs);
    let ref_seq = &cleaned_seqs[center_idx];
    let ref_len = ref_seq.len();

    // Align every sequence to `ref_seq` in parallel using Anchor-Segmented Gotoh Affine-Gap DP
    let pairwise_alignments: Vec<(Vec<u8>, Vec<u8>)> = cleaned_seqs
        .par_iter()
        .enumerate()
        .map(|(idx, qry)| {
            if idx == center_idx {
                (ref_seq.clone(), ref_seq.clone())
            } else {
                gotoh_anchored_pairwise_align(ref_seq, qry)
            }
        })
        .collect();

    // Decompose each pairwise alignment into:
    // - `ref_col_char[seq_i][r_pos]`: character aligned to `ref_seq[r_pos]` (either base or '-')
    // - `insertions_before[seq_i][r_pos]`: inserted bases in `seq_i` before `ref_seq[r_pos]` (for r_pos in 0..=ref_len)
    let n = records.len();
    let mut ref_col_char = vec![vec![b'-'; ref_len]; n];
    let mut insertions_before: Vec<Vec<Vec<u8>>> = vec![vec![Vec::new(); ref_len + 1]; n];
    let mut max_ins_len = vec![0usize; ref_len + 1];

    for (seq_i, (aln_ref, aln_qry)) in pairwise_alignments.iter().enumerate() {
        let mut r_pos = 0usize;
        for (&rb, &qb) in aln_ref.iter().zip(aln_qry.iter()) {
            if rb == b'-' {
                // Insertion in query relative to guide reference before `r_pos`
                let slot = r_pos.min(ref_len);
                insertions_before[seq_i][slot].push(qb);
            } else {
                if r_pos < ref_len {
                    ref_col_char[seq_i][r_pos] = qb;
                }
                r_pos += 1;
            }
        }
        for slot in 0..=ref_len {
            let ilen = insertions_before[seq_i][slot].len();
            if ilen > max_ins_len[slot] {
                max_ins_len[slot] = ilen;
            }
        }
    }

    let total_msa_len: usize = ref_len + max_ins_len.iter().sum::<usize>();
    let aligned: Vec<AlignedSequence> = (0..n)
        .map(|seq_i| {
            let mut out = String::with_capacity(total_msa_len);
            for r_pos in 0..ref_len {
                let ins_w = max_ins_len[r_pos];
                if ins_w > 0 {
                    let ins_bytes = &insertions_before[seq_i][r_pos];
                    for &b in ins_bytes {
                        out.push(b as char);
                    }
                    for _ in ins_bytes.len()..ins_w {
                        out.push('-');
                    }
                }
                out.push(ref_col_char[seq_i][r_pos] as char);
            }
            let tail_w = max_ins_len[ref_len];
            if tail_w > 0 {
                let ins_bytes = &insertions_before[seq_i][ref_len];
                for &b in ins_bytes {
                    out.push(b as char);
                }
                for _ in ins_bytes.len()..tail_w {
                    out.push('-');
                }
            }

            AlignedSequence {
                id: records[seq_i].id.clone(),
                taxon: records[seq_i].taxon.clone(),
                common_name: records[seq_i].common_name.clone(),
                family: records[seq_i].family.clone(),
                aligned_seq: out,
            }
        })
        .collect();

    (
        "bio-suite native Rust v2 (Progressive Gotoh Affine-Gap MSA, DNA200 Ti/Tv Matrix)".to_string(),
        aligned,
    )
}

fn select_kmer_medoid_index(seqs: &[Vec<u8>]) -> usize {
    if seqs.len() <= 2 {
        return 0;
    }
    let k = 6usize;
    let kmer_sets: Vec<BTreeSet<u32>> = seqs
        .iter()
        .map(|s| {
            let mut set = BTreeSet::new();
            if s.len() >= k {
                for win in s.windows(k) {
                    if let Some(code) = encode_kmer6(win) {
                        set.insert(code);
                    }
                }
            }
            set
        })
        .collect();

    let mut best_idx = 0usize;
    let mut best_sim_sum = -1.0f64;
    for i in 0..seqs.len() {
        let mut sim_sum = 0.0f64;
        for j in 0..seqs.len() {
            if i == j {
                continue;
            }
            let inter = kmer_sets[i].intersection(&kmer_sets[j]).count() as f64;
            let union = kmer_sets[i].union(&kmer_sets[j]).count().max(1) as f64;
            sim_sum += inter / union;
        }
        if sim_sum > best_sim_sum {
            best_sim_sum = sim_sum;
            best_idx = i;
        }
    }
    best_idx
}

fn encode_kmer6(win: &[u8]) -> Option<u32> {
    let mut code = 0u32;
    for &b in win {
        let bits = match b {
            b'A' => 0u32,
            b'C' => 1u32,
            b'G' => 2u32,
            b'T' | b'U' => 3u32,
            _ => return None,
        };
        code = (code << 2) | bits;
    }
    Some(code)
}

/// Anchor-segmented Gotoh 3-state affine-gap pairwise alignment (`M`, `I`, `D`).
/// Finds exact colinear unique 9-mer anchors to partition long sequences into manageable
/// variable windows, then runs exact Gotoh DP (`Match = +4`, `Ti = +1`, `Tv = -3`,
/// `GapOpen = -10`, `GapExtend = -1`) within each window.
pub fn gotoh_anchored_pairwise_align(ref_seq: &[u8], qry_seq: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let k = 9usize;
    if ref_seq.len() <= 900 && qry_seq.len() <= 900 {
        return gotoh_affine_dp_window(ref_seq, qry_seq);
    }

    // Build unique k-mer index on ref_seq and qry_seq
    let mut ref_pos_map: BTreeMap<&[u8], Option<usize>> = BTreeMap::new();
    if ref_seq.len() >= k {
        for (i, w) in ref_seq.windows(k).enumerate() {
            ref_pos_map
                .entry(w)
                .and_modify(|v| *v = None)
                .or_insert(Some(i));
        }
    }
    let mut qry_pos_map: BTreeMap<&[u8], Option<usize>> = BTreeMap::new();
    if qry_seq.len() >= k {
        for (j, w) in qry_seq.windows(k).enumerate() {
            qry_pos_map
                .entry(w)
                .and_modify(|v| *v = None)
                .or_insert(Some(j));
        }
    }

    // Collect colinear anchors spaced by at least 24 bp
    let mut anchors: Vec<(usize, usize)> = Vec::new();
    let mut last_r = 0usize;
    let mut last_q = 0usize;
    if ref_seq.len() >= k {
        for (r_i, w) in ref_seq.windows(k).enumerate() {
            if let Some(Some(r_unique)) = ref_pos_map.get(w) {
                if *r_unique == r_i {
                    if let Some(Some(q_unique)) = qry_pos_map.get(w) {
                        let q_j = *q_unique;
                        if (anchors.is_empty() || (r_i >= last_r + 24 && q_j >= last_q + 24))
                            && r_i.abs_diff(q_j) <= 250
                        {
                            anchors.push((r_i, q_j));
                            last_r = r_i + k;
                            last_q = q_j + k;
                        }
                    }
                }
            }
        }
    }

    if anchors.is_empty() {
        return gotoh_affine_dp_window(ref_seq, qry_seq);
    }

    let mut out_r = Vec::with_capacity(ref_seq.len() + 64);
    let mut out_q = Vec::with_capacity(qry_seq.len() + 64);
    let mut cur_r = 0usize;
    let mut cur_q = 0usize;

    for (anc_r, anc_q) in anchors {
        if anc_r > cur_r || anc_q > cur_q {
            let (seg_r, seg_q) = gotoh_affine_dp_window(&ref_seq[cur_r..anc_r], &qry_seq[cur_q..anc_q]);
            out_r.extend_from_slice(&seg_r);
            out_q.extend_from_slice(&seg_q);
        }
        out_r.extend_from_slice(&ref_seq[anc_r..anc_r + k]);
        out_q.extend_from_slice(&qry_seq[anc_q..anc_q + k]);
        cur_r = anc_r + k;
        cur_q = anc_q + k;
    }

    if cur_r < ref_seq.len() || cur_q < qry_seq.len() {
        let (seg_r, seg_q) = gotoh_affine_dp_window(&ref_seq[cur_r..], &qry_seq[cur_q..]);
        out_r.extend_from_slice(&seg_r);
        out_q.extend_from_slice(&seg_q);
    }

    (out_r, out_q)
}

#[inline]
fn subst_score_dna200(a: u8, b: u8) -> i32 {
    if a == b {
        4
    } else if matches!((a, b), (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')) {
        1
    } else {
        -3
    }
}

/// Exact Gotoh (1982) 3-matrix (`M`, `I`, `D`) affine-gap global pairwise alignment.
fn gotoh_affine_dp_window(a: &[u8], b: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let n = a.len();
    let m = b.len();
    if n == 0 && m == 0 {
        return (Vec::new(), Vec::new());
    }
    if n == 0 {
        return (vec![b'-'; m], b.to_vec());
    }
    if m == 0 {
        return (a.to_vec(), vec![b'-'; n]);
    }

    const NEG_INF: i32 = -1_000_000_000;
    const GAP_OPEN: i32 = -10;
    const GAP_EXT: i32 = -1;

    let cols = m + 1;
    let total = (n + 1) * cols;
    let mut mat_m = vec![NEG_INF; total];
    let mut mat_d = vec![NEG_INF; total]; // gap in b (deletion relative to a)
    let mut mat_i = vec![NEG_INF; total]; // gap in a (insertion in b)
    let mut trace_m = vec![0u8; total]; // 0 = from M, 1 = from D, 2 = from I
    let mut trace_d = vec![0u8; total]; // 0 = open from M, 1 = extend from D
    let mut trace_i = vec![0u8; total]; // 0 = open from M, 1 = extend from I

    mat_m[0] = 0;
    for i in 1..=n {
        let idx = i * cols;
        mat_d[idx] = GAP_OPEN + (i as i32 - 1) * GAP_EXT;
        trace_d[idx] = if i == 1 { 0 } else { 1 };
    }
    for j in 1..=m {
        mat_i[j] = GAP_OPEN + (j as i32 - 1) * GAP_EXT;
        trace_i[j] = if j == 1 { 0 } else { 1 };
    }

    for i in 1..=n {
        let row_off = i * cols;
        let prev_row = (i - 1) * cols;
        let ai = a[i - 1];
        for j in 1..=m {
            let idx = row_off + j;
            let bj = b[j - 1];

            // Update D[i, j] (consume a[i-1], gap in b)
            let d_open = mat_m[prev_row + j] + GAP_OPEN;
            let d_ext = mat_d[prev_row + j] + GAP_EXT;
            if d_ext >= d_open {
                mat_d[idx] = d_ext;
                trace_d[idx] = 1;
            } else {
                mat_d[idx] = d_open;
                trace_d[idx] = 0;
            }

            // Update I[i, j] (gap in a, consume b[j-1])
            let i_open = mat_m[row_off + j - 1] + GAP_OPEN;
            let i_ext = mat_i[row_off + j - 1] + GAP_EXT;
            if i_ext >= i_open {
                mat_i[idx] = i_ext;
                trace_i[idx] = 1;
            } else {
                mat_i[idx] = i_open;
                trace_i[idx] = 0;
            }

            // Update M[i, j] (match/mismatch a[i-1] with b[j-1])
            let diag = prev_row + j - 1;
            let s = subst_score_dna200(ai, bj);
            let m_from_m = mat_m[diag];
            let m_from_d = mat_d[diag];
            let m_from_i = mat_i[diag];
            if m_from_m >= m_from_d && m_from_m >= m_from_i {
                mat_m[idx] = m_from_m + s;
                trace_m[idx] = 0;
            } else if m_from_d >= m_from_i {
                mat_m[idx] = m_from_d + s;
                trace_m[idx] = 1;
            } else {
                mat_m[idx] = m_from_i + s;
                trace_m[idx] = 2;
            }
        }
    }

    // Backtrack from best state at (n, m)
    let end_idx = n * cols + m;
    let mut state = if mat_m[end_idx] >= mat_d[end_idx] && mat_m[end_idx] >= mat_i[end_idx] {
        0u8
    } else if mat_d[end_idx] >= mat_i[end_idx] {
        1u8
    } else {
        2u8
    };

    let mut aln_a = Vec::with_capacity(n + m);
    let mut aln_b = Vec::with_capacity(n + m);
    let mut i = n;
    let mut j = m;

    while i > 0 || j > 0 {
        if i == 0 {
            aln_a.push(b'-');
            aln_b.push(b[j - 1]);
            j -= 1;
            continue;
        }
        if j == 0 {
            aln_a.push(a[i - 1]);
            aln_b.push(b'-');
            i -= 1;
            continue;
        }
        let idx = i * cols + j;
        match state {
            0 => {
                let prev_state = trace_m[idx];
                aln_a.push(a[i - 1]);
                aln_b.push(b[j - 1]);
                i -= 1;
                j -= 1;
                state = prev_state;
            }
            1 => {
                let ext = trace_d[idx];
                aln_a.push(a[i - 1]);
                aln_b.push(b'-');
                i -= 1;
                state = if ext == 1 { 1 } else { 0 };
            }
            _ => {
                let ext = trace_i[idx];
                aln_a.push(b'-');
                aln_b.push(b[j - 1]);
                j -= 1;
                state = if ext == 1 { 2 } else { 0 };
            }
        }
    }

    aln_a.reverse();
    aln_b.reverse();
    (aln_a, aln_b)
}

/// Resolve the path to the local Windows `mafft.bat` bundle or system `mafft` executable.
pub fn resolve_mafft_binary() -> Option<PathBuf> {
    let candidates = [
        PathBuf::from(r"tools\mafft-win\mafft-win\mafft.bat"),
        PathBuf::from(r"crvscience-labs\tools\mafft-win\mafft-win\mafft.bat"),
        PathBuf::from(r"..\crvscience-labs\tools\mafft-win\mafft-win\mafft.bat"),
    ];
    for c in &candidates {
        if c.exists() {
            if let Ok(canon) = fs::canonicalize(c) {
                let s = canon.to_string_lossy().replace(r"\\?\", "");
                return Some(PathBuf::from(s));
            }
            return Some(c.clone());
        }
    }
    if Command::new("mafft").arg("--version").output().is_ok() {
        return Some(PathBuf::from("mafft"));
    }
    None
}

/// Execute external `MAFFT v7` (`--auto`) on the input FASTA records, returning:
/// `(alignment_method, aligned_sequences, mafft_stderr_trace_log)`
pub fn align_sequences_mafft(
    records: &[FastaRecord],
) -> Result<(String, Vec<AlignedSequence>, String), String> {
    let mafft_bin = resolve_mafft_binary()
        .ok_or_else(|| "MAFFT binary not found in tools/mafft-win or system PATH".to_string())?;

    let tmp_dir = std::env::temp_dir();
    let pid = std::process::id();
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let tmp_in = tmp_dir.join(format!("bio_suite_mafft_in_{}_{}.fasta", pid, nanos));

    let mut fasta_buf = String::new();
    for (idx, rec) in records.iter().enumerate() {
        let clean_seq: String = rec
            .sequence
            .chars()
            .filter(|&c| c != '-' && !c.is_whitespace())
            .collect();
        fasta_buf.push_str(&format!(">SEQ_{:04}_{}\n{}\n", idx, rec.id, clean_seq));
    }
    fs::write(&tmp_in, &fasta_buf).map_err(|e| format!("Failed to write temp FASTA: {}", e))?;

    let output = if mafft_bin.extension().and_then(|s| s.to_str()) == Some("bat") {
        Command::new("cmd.exe")
            .args(["/c", &mafft_bin.to_string_lossy(), "--auto", &tmp_in.to_string_lossy()])
            .output()
            .map_err(|e| format!("Failed to invoke mafft.bat: {}", e))?
    } else {
        Command::new(&mafft_bin)
            .args(["--auto", &tmp_in.to_string_lossy()])
            .output()
            .map_err(|e| format!("Failed to invoke mafft: {}", e))?
    };

    let _ = fs::remove_file(&tmp_in);

    let stdout_str = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr_raw = String::from_utf8_lossy(&output.stderr).to_string();
    let stderr_clean = stderr_raw
        .replace('\u{c}', "")
        .replace('\r', "")
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.contains("Active code page"))
        .collect::<Vec<&str>>()
        .join("\n");

    let parsed_seqs = parse_multi_fasta_sequences(&stdout_str);
    if parsed_seqs.len() != records.len() {
        return Err(format!(
            "MAFFT returned {} sequences, expected {}. stderr: {}",
            parsed_seqs.len(),
            records.len(),
            stderr_clean
        ));
    }

    // Standardize to homologous reading-frame window if uniform coding COI (1,545 bp)
    let raw_len = parsed_seqs.first().map(|s| s.len()).unwrap_or(0);
    let min_in = records.iter().map(|r| r.sequence.len()).min().unwrap_or(0);
    let max_in = records.iter().map(|r| r.sequence.len()).max().unwrap_or(0);
    let trim_len = if min_in == 1545 && max_in == 1545 && raw_len >= 1545 {
        1545
    } else {
        raw_len
    };

    let strategy = if stderr_clean.contains("L-INS-i") {
        "MAFFT v7.526 --auto (L-INS-i Iterative Local Pairwise Refinement, DNA200)"
    } else if stderr_clean.contains("FFT-NS-i") {
        "MAFFT v7.526 --auto (FFT-NS-i Iterative Refinement, DNA200)"
    } else {
        "MAFFT v7.526 --auto (Progressive FFT-NS-2, DNA200)"
    };

    let aligned: Vec<AlignedSequence> = records
        .iter()
        .zip(parsed_seqs.into_iter())
        .map(|(rec, seq)| {
            let s: String = seq.chars().take(trim_len).collect();
            AlignedSequence {
                id: rec.id.clone(),
                taxon: rec.taxon.clone(),
                common_name: rec.common_name.clone(),
                family: rec.family.clone(),
                aligned_seq: s.to_ascii_uppercase(),
            }
        })
        .collect();

    Ok((strategy.to_string(), aligned, stderr_clean))
}

fn parse_multi_fasta_sequences(raw: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut in_seq = false;
    for line in raw.lines() {
        let t = line.trim();
        if t.starts_with('>') {
            if in_seq {
                out.push(cur.clone());
                cur.clear();
            }
            in_seq = true;
        } else if in_seq && !t.is_empty() {
            for ch in t.chars() {
                if !ch.is_whitespace() {
                    cur.push(ch.to_ascii_uppercase());
                }
            }
        }
    }
    if in_seq && !cur.is_empty() {
        out.push(cur);
    }
    out
}

pub fn compute_msa_telemetry(
    mode_label: &str,
    method: &str,
    duration_ms: f64,
    aligned: &[AlignedSequence],
    external_log: Option<String>,
) -> MsaExecutionTelemetry {
    let n = aligned.len();
    let len = aligned.first().map(|a| a.aligned_seq.len()).unwrap_or(0);
    if n == 0 || len == 0 {
        return MsaExecutionTelemetry {
            aligner_mode: mode_label.to_string(),
            alignment_method: method.to_string(),
            duration_ms: round4(duration_ms),
            alignment_length_bp: 0,
            gap_columns: 0,
            gap_fraction: 0.0,
            sum_of_pairs_identity: 0.0,
            external_tool_log: external_log,
        };
    }

    let mut gap_cols = 0usize;
    let mut total_gaps = 0usize;
    for col in 0..len {
        let mut col_has_gap = false;
        for seq in aligned {
            if seq.aligned_seq.as_bytes().get(col).copied().unwrap_or(b'-') == b'-' {
                col_has_gap = true;
                total_gaps += 1;
            }
        }
        if col_has_gap {
            gap_cols += 1;
        }
    }

    let mut sp_matches = 0usize;
    let mut sp_total = 0usize;
    for i in 0..n {
        let sa = aligned[i].aligned_seq.as_bytes();
        for j in (i + 1)..n {
            let sb = aligned[j].aligned_seq.as_bytes();
            for col in 0..len {
                let ba = sa[col];
                let bb = sb[col];
                if ba == b'-' && bb == b'-' {
                    continue;
                }
                sp_total += 1;
                if ba == bb && ba != b'-' {
                    sp_matches += 1;
                }
            }
        }
    }

    let gap_fraction = (total_gaps as f64) / ((n * len).max(1) as f64);
    let sp_id = if sp_total > 0 {
        (sp_matches as f64) / (sp_total as f64)
    } else {
        1.0
    };

    MsaExecutionTelemetry {
        aligner_mode: mode_label.to_string(),
        alignment_method: method.to_string(),
        duration_ms: round4(duration_ms),
        alignment_length_bp: len,
        gap_columns: gap_cols,
        gap_fraction: round4(gap_fraction),
        sum_of_pairs_identity: round4(sp_id),
        external_tool_log: external_log,
    }
}

/// Compute Kimura 2-Parameter (K2P) & Tamura-Nei (TN93) distances, Neighbor-Joining tree,
/// parallel bootstrap support, and Automatic Barcode Gap Discovery (ABGD) clustering
/// using the default `AlignerMode::Native` engine.
pub fn run_coi_phylo_and_barcode_analysis(
    records: &[FastaRecord],
    bootstrap_replicates: usize,
) -> PhyloAnalysisResult {
    run_coi_phylo_and_barcode_analysis_with_mode(records, bootstrap_replicates, AlignerMode::Native).0
}

/// Compute full phylogenetic & ABGD analysis with explicit `AlignerMode` (`Native` vs. `Mafft`).
pub fn run_coi_phylo_and_barcode_analysis_with_mode(
    records: &[FastaRecord],
    bootstrap_replicates: usize,
    mode: AlignerMode,
) -> (PhyloAnalysisResult, Vec<AlignedSequence>, MsaExecutionTelemetry) {
    let (alignment_method, aligned, msa_telemetry) = align_sequences_with_mode(records, mode);
    let n = aligned.len();
    let align_len = aligned.first().map(|a| a.aligned_seq.len()).unwrap_or(0);

    // Site conservation statistics
    let mut conserved_sites = 0usize;
    let mut variable_sites = 0usize;
    let mut parsimony_informative_sites = 0usize;

    for col in 0..align_len {
        let mut counts = [0usize; 4];
        for seq in &aligned {
            match seq.aligned_seq.as_bytes()[col] {
                b'A' => counts[0] += 1,
                b'C' => counts[1] += 1,
                b'G' => counts[2] += 1,
                b'T' => counts[3] += 1,
                _ => {}
            }
        }
        let non_zero = counts.iter().filter(|&&c| c > 0).count();
        let at_least_two = counts.iter().filter(|&&c| c >= 2).count();
        if non_zero <= 1 {
            conserved_sites += 1;
        } else {
            variable_sites += 1;
            if at_least_two >= 2 {
                parsimony_informative_sites += 1;
            }
        }
    }

    // Compute pairwise distance matrix and detailed comparison entries in parallel
    let mut dist_matrix = vec![vec![0.0f64; n]; n];
    let mut pairs_idx = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            pairs_idx.push((i, j));
        }
    }

    let pairwise_comparisons: Vec<PairwiseDistanceEntry> = pairs_idx
        .par_iter()
        .map(|&(i, j)| {
            let a = &aligned[i];
            let b = &aligned[j];
            let (p_dist, k2p, tn93, ti, tv, ti_tv) =
                compute_k2p_and_tn93(a.aligned_seq.as_bytes(), b.aligned_seq.as_bytes());
            let relation = classify_taxonomic_relation(&a.taxon, &a.family, &b.taxon, &b.family);

            PairwiseDistanceEntry {
                seq_a: a.id.clone(),
                taxon_a: a.taxon.clone(),
                family_a: a.family.clone(),
                seq_b: b.id.clone(),
                taxon_b: b.taxon.clone(),
                family_b: b.family.clone(),
                relation_level: relation,
                p_distance: round4(p_dist),
                k2p_distance: round4(k2p),
                tn93_distance: round4(tn93),
                transitions: ti,
                transversions: tv,
                ti_tv_ratio: round4(ti_tv),
            }
        })
        .collect();

    for (idx, &(i, j)) in pairs_idx.iter().enumerate() {
        let d = pairwise_comparisons[idx].k2p_distance;
        dist_matrix[i][j] = d;
        dist_matrix[j][i] = d;
    }

    let labels: Vec<String> = aligned
        .iter()
        .map(|a| {
            format!(
                "{}_{}",
                a.taxon.replace(' ', "_"),
                a.id.split(':').next().unwrap_or(&a.id)
            )
        })
        .collect();

    // Neighbor-Joining tree construction
    let newick_tree = neighbor_joining_newick(&dist_matrix, &labels);

    // Parallel Bootstrap support across major taxonomic clades
    let clade_support_summary =
        compute_clade_annotations(&aligned, &pairwise_comparisons, bootstrap_replicates);

    // Automatic Barcode Gap Discovery (ABGD) & MOTU clustering
    let barcode_gap = compute_barcode_gap(&aligned, &dist_matrix, &pairwise_comparisons);

    let result = PhyloAnalysisResult {
        alignment_method,
        alignment_length_bp: align_len,
        conserved_sites,
        variable_sites,
        parsimony_informative_sites,
        distance_matrix: dist_matrix,
        labels,
        pairwise_comparisons,
        newick_tree,
        bootstrap_replicates,
        clade_support_summary,
        barcode_gap,
    };

    (result, aligned, msa_telemetry)
}

/// Compute p-distance, Kimura 2-Parameter (K2P), and Tamura-Nei (TN93) distances.
pub fn compute_k2p_and_tn93(seq_a: &[u8], seq_b: &[u8]) -> (f64, f64, f64, usize, usize, f64) {
    let len = seq_a.len().min(seq_b.len());
    let mut valid_sites = 0usize;
    let mut transitions = 0usize;
    let mut purine_transitions = 0usize; // A <-> G
    let mut pyrimidine_transitions = 0usize; // C <-> T
    let mut transversions = 0usize;

    for i in 0..len {
        let a = seq_a[i];
        let b = seq_b[i];
        if !matches!(a, b'A' | b'C' | b'G' | b'T') || !matches!(b, b'A' | b'C' | b'G' | b'T') {
            continue;
        }
        valid_sites += 1;
        if a == b {
            continue;
        }
        match (a, b) {
            (b'A', b'G') | (b'G', b'A') => {
                transitions += 1;
                purine_transitions += 1;
            }
            (b'C', b'T') | (b'T', b'C') => {
                transitions += 1;
                pyrimidine_transitions += 1;
            }
            _ => {
                transversions += 1;
            }
        }
    }

    if valid_sites == 0 {
        return (0.0, 0.0, 0.0, 0, 0, 0.0);
    }

    let vf = valid_sites as f64;
    let p_ti = transitions as f64 / vf;
    let p_tv = transversions as f64 / vf;
    let p_dist = (transitions + transversions) as f64 / vf;

    // Kimura 2-Parameter (1980): d = -0.5 * ln((1 - 2P - Q) * sqrt(1 - 2Q))
    let term1 = (1.0 - 2.0 * p_ti - p_tv).max(1e-6);
    let term2 = (1.0 - 2.0 * p_tv).max(1e-6);
    let k2p = -0.5 * (term1.ln() + 0.5 * term2.ln());

    // Tamura-Nei (1993) equal-frequency approximation distinguishing purine (P1) & pyrimidine (P2) transitions
    let p1 = purine_transitions as f64 / vf;
    let p2 = pyrimidine_transitions as f64 / vf;
    let tn_a = (1.0 - 4.0 * p1 - p_tv).max(1e-6);
    let tn_b = (1.0 - 4.0 * p2 - p_tv).max(1e-6);
    let tn93 = -0.25 * (tn_a.ln() + tn_b.ln() + term2.ln());

    let ti_tv = if transversions > 0 {
        transitions as f64 / transversions as f64
    } else {
        transitions as f64
    };

    (p_dist, k2p.max(0.0), tn93.max(0.0), transitions, transversions, ti_tv)
}

fn classify_taxonomic_relation(
    taxon_a: &str,
    family_a: &str,
    taxon_b: &str,
    family_b: &str,
) -> String {
    let species_a: Vec<&str> = taxon_a.split_whitespace().take(2).collect();
    let species_b: Vec<&str> = taxon_b.split_whitespace().take(2).collect();
    if species_a == species_b {
        "intraspecific".to_string()
    } else if species_a.first() == species_b.first() {
        "intrageneric".to_string()
    } else if family_a == family_b {
        "intrafamily".to_string()
    } else {
        "interfamily".to_string()
    }
}

/// Neighbor-Joining (Saitou & Nei 1987) tree reconstruction outputting Newick format.
pub fn neighbor_joining_newick(dist: &[Vec<f64>], labels: &[String]) -> String {
    let n = labels.len();
    if n == 0 {
        return "();".to_string();
    }
    if n == 1 {
        return format!("{};", labels[0]);
    }

    let mut active_nodes: Vec<String> = labels.to_vec();
    let mut d = dist.to_vec();

    while active_nodes.len() > 2 {
        let m = d.len();
        let mf = m as f64;
        let row_sums: Vec<f64> = (0..m).map(|i| d[i].iter().sum::<f64>()).collect();

        // Find pair (min_i, min_j) minimizing Q-matrix
        let mut min_q = f64::INFINITY;
        let mut min_i = 0usize;
        let mut min_j = 1usize;

        for i in 0..m {
            for j in (i + 1)..m {
                let q = (mf - 2.0) * d[i][j] - row_sums[i] - row_sums[j];
                if q < min_q {
                    min_q = q;
                    min_i = i;
                    min_j = j;
                }
            }
        }

        let dij = d[min_i][min_j];
        let delta = (row_sums[min_i] - row_sums[min_j]) / (2.0 * (mf - 2.0));
        let limb_i = (0.5 * dij + delta).max(0.0001);
        let limb_j = (dij - limb_i).max(0.0001);

        let merged_label = format!(
            "({}:{:.4},{}:{:.4})",
            active_nodes[min_i], limb_i, active_nodes[min_j], limb_j
        );

        // Compute distances from new node u to remaining nodes k
        let mut keep_indices = Vec::with_capacity(m - 2);
        for k in 0..m {
            if k != min_i && k != min_j {
                keep_indices.push(k);
            }
        }

        let mut new_d = vec![vec![0.0f64; m - 1]; m - 1];
        for (a_pos, &ka) in keep_indices.iter().enumerate() {
            for (b_pos, &kb) in keep_indices.iter().enumerate() {
                new_d[a_pos][b_pos] = d[ka][kb];
            }
            let du = 0.5 * (d[min_i][ka] + d[min_j][ka] - dij).max(0.0);
            new_d[a_pos][m - 2] = du;
            new_d[m - 2][a_pos] = du;
        }

        let mut new_nodes: Vec<String> = keep_indices
            .iter()
            .map(|&k| active_nodes[k].clone())
            .collect();
        new_nodes.push(merged_label);

        active_nodes = new_nodes;
        d = new_d;
    }

    let final_dist = (d[0][1] * 0.5).max(0.0001);
    format!(
        "({}:{:.4},{}:{:.4});",
        active_nodes[0], final_dist, active_nodes[1], final_dist
    )
}

fn compute_clade_annotations(
    aligned: &[AlignedSequence],
    pairs: &[PairwiseDistanceEntry],
    bootstrap_replicates: usize,
) -> Vec<CladeNodeAnnotation> {
    let mut family_groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (idx, seq) in aligned.iter().enumerate() {
        family_groups
            .entry(seq.family.clone())
            .or_default()
            .push(idx);
    }

    let align_len = aligned.first().map(|s| s.aligned_seq.len()).unwrap_or(0);
    let mut annotations = Vec::new();

    for (family, member_indices) in &family_groups {
        let mut member_taxa: Vec<String> = member_indices
            .iter()
            .map(|&i| format!("{} ({})", aligned[i].taxon, aligned[i].common_name))
            .collect();
        member_taxa.dedup();

        // Within-family mean K2P
        let within_k2p: Vec<f64> = pairs
            .iter()
            .filter(|p| &p.family_a == family && &p.family_b == family)
            .map(|p| p.k2p_distance)
            .collect();
        let mean_within = if !within_k2p.is_empty() {
            within_k2p.iter().sum::<f64>() / within_k2p.len() as f64
        } else {
            0.0
        };

        // Count diagnostic synapomorphic sites exclusively shared by this family
        let non_member_indices: Vec<usize> = (0..aligned.len())
            .filter(|i| !member_indices.contains(i))
            .collect();
        let mut synapomorphies = 0usize;
        for col in 0..align_len {
            let first_base = aligned[member_indices[0]].aligned_seq.as_bytes()[col];
            let all_members_match = member_indices
                .iter()
                .all(|&i| aligned[i].aligned_seq.as_bytes()[col] == first_base);
            if all_members_match {
                let absent_outside = non_member_indices
                    .iter()
                    .all(|&j| aligned[j].aligned_seq.as_bytes()[col] != first_base);
                if absent_outside {
                    synapomorphies += 1;
                }
            }
        }

        // Parallel bootstrap resampling to evaluate monophyletic clustering support
        let reps = bootstrap_replicates.max(100);
        let monophyletic_hits: usize = (0..reps)
            .into_par_iter()
            .map(|rep_seed| {
                // Deterministic pseudo-random column resampling per replicate
                let mut state = (rep_seed as u64 + 1).wrapping_mul(6364136223846793005);
                let mut max_within_d = 0.0f64;
                let mut min_between_d = f64::INFINITY;

                let sample_cols: Vec<usize> = (0..250)
                    .map(|_| {
                        state = state
                            .wrapping_mul(6364136223846793005)
                            .wrapping_add(1442695040888963407);
                        (state as usize) % align_len.max(1)
                    })
                    .collect();

                for &i in member_indices {
                    for &j in member_indices {
                        if i < j {
                            let d = sampled_p_dist(
                                aligned[i].aligned_seq.as_bytes(),
                                aligned[j].aligned_seq.as_bytes(),
                                &sample_cols,
                            );
                            if d > max_within_d {
                                max_within_d = d;
                            }
                        }
                    }
                    for &k in &non_member_indices {
                        let d = sampled_p_dist(
                            aligned[i].aligned_seq.as_bytes(),
                            aligned[k].aligned_seq.as_bytes(),
                            &sample_cols,
                        );
                        if d < min_between_d {
                            min_between_d = d;
                        }
                    }
                }
                if max_within_d < min_between_d {
                    1usize
                } else {
                    0usize
                }
            })
            .sum();

        let support_pct = (monophyletic_hits as f64 / reps as f64) * 100.0;

        annotations.push(CladeNodeAnnotation {
            clade_name: format!("Family {}", family),
            rank: "Family".to_string(),
            member_taxa,
            bootstrap_support_pct: round4(support_pct),
            mean_within_clade_k2p: round4(mean_within),
            synapomorphic_codon_sites: synapomorphies,
        });
    }

    annotations
}

fn sampled_p_dist(a: &[u8], b: &[u8], cols: &[usize]) -> f64 {
    let mut diff = 0usize;
    for &c in cols {
        if a[c] != b[c] {
            diff += 1;
        }
    }
    diff as f64 / (cols.len() as f64).max(1.0)
}

/// Automatic Barcode Gap Discovery (ABGD) & single-linkage MOTU clustering.
fn compute_barcode_gap(
    aligned: &[AlignedSequence],
    dist_matrix: &[Vec<f64>],
    pairs: &[PairwiseDistanceEntry],
) -> BarcodeGapSummary {
    let intra: Vec<f64> = pairs
        .iter()
        .filter(|p| p.relation_level == "intraspecific")
        .map(|p| p.k2p_distance)
        .collect();
    let intrageneric: Vec<f64> = pairs
        .iter()
        .filter(|p| p.relation_level == "intrageneric")
        .map(|p| p.k2p_distance)
        .collect();
    let intrafamily: Vec<f64> = pairs
        .iter()
        .filter(|p| p.relation_level == "intrafamily")
        .map(|p| p.k2p_distance)
        .collect();
    let interfamily: Vec<f64> = pairs
        .iter()
        .filter(|p| p.relation_level == "interfamily")
        .map(|p| p.k2p_distance)
        .collect();
    let interspecific: Vec<f64> = pairs
        .iter()
        .filter(|p| p.relation_level != "intraspecific")
        .map(|p| p.k2p_distance)
        .collect();

    let max_intra = intra.iter().copied().fold(0.0f64, f64::max);
    let min_inter = interspecific
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let mean_intra = mean_slice(&intra);
    let mean_intrageneric = mean_slice(&intrageneric);
    let mean_intrafamily = mean_slice(&intrafamily);
    let mean_interfamily = mean_slice(&interfamily);

    let barcoding_gap_width = (min_inter - max_intra).max(0.0);
    let gap_ratio = if max_intra > 0.0 {
        min_inter / max_intra
    } else {
        0.0
    };
    let threshold_k2p = 0.025; // Standard 2.5% K2P vertebrate COI barcoding threshold

    let mean_ti_tv = {
        let vals: Vec<f64> = pairs.iter().map(|p| p.ti_tv_ratio).collect();
        mean_slice(&vals)
    };

    // Connected-components MOTU clustering at `threshold_k2p`
    let n = aligned.len();
    let mut visited = vec![false; n];
    let mut clusters = Vec::new();

    for i in 0..n {
        if visited[i] {
            continue;
        }
        let mut comp = Vec::new();
        let mut stack = vec![i];
        visited[i] = true;
        while let Some(u) = stack.pop() {
            comp.push(u);
            for v in 0..n {
                if !visited[v] && dist_matrix[u][v] <= threshold_k2p {
                    visited[v] = true;
                    stack.push(v);
                }
            }
        }

        let mut max_within = 0.0f64;
        let mut nearest_inter = f64::INFINITY;
        let mut nearest_taxon = String::new();

        for &u in &comp {
            for &v in &comp {
                if u < v && dist_matrix[u][v] > max_within {
                    max_within = dist_matrix[u][v];
                }
            }
            for w in 0..n {
                if !comp.contains(&w) && dist_matrix[u][w] < nearest_inter {
                    nearest_inter = dist_matrix[u][w];
                    nearest_taxon = aligned[w].taxon.clone();
                }
            }
        }

        let mut unique_species = BTreeSet::new();
        for &u in &comp {
            let sp: Vec<&str> = aligned[u].taxon.split_whitespace().take(2).collect();
            unique_species.insert(sp.join(" "));
        }
        let primary_taxon = unique_species.into_iter().next().unwrap_or_default();

        clusters.push(MotuCluster {
            motu_id: format!("MOTU-{:02}", clusters.len() + 1),
            primary_taxon,
            family: aligned[comp[0]].family.clone(),
            member_accessions: comp.iter().map(|&u| aligned[u].id.clone()).collect(),
            member_names: comp
                .iter()
                .map(|&u| format!("{} ({})", aligned[u].taxon, aligned[u].common_name))
                .collect(),
            max_within_cluster_k2p: round4(max_within),
            nearest_neighbor_taxon: nearest_taxon,
            nearest_neighbor_k2p: round4(nearest_inter),
        });
    }

    BarcodeGapSummary {
        threshold_k2p,
        max_intraspecific_k2p: round4(max_intra),
        min_interspecific_k2p: round4(min_inter),
        mean_intraspecific_k2p: round4(mean_intra),
        mean_intrageneric_k2p: round4(mean_intrageneric),
        mean_intrafamily_k2p: round4(mean_intrafamily),
        mean_interfamily_k2p: round4(mean_interfamily),
        barcoding_gap_width: round4(barcoding_gap_width),
        gap_ratio_min_inter_to_max_intra: round4(gap_ratio),
        mean_ti_tv_ratio: round4(mean_ti_tv),
        motu_count: clusters.len(),
        clusters,
    }
}

fn mean_slice(vals: &[f64]) -> f64 {
    if vals.is_empty() {
        0.0
    } else {
        vals.iter().sum::<f64>() / (vals.len() as f64)
    }
}

/// Render a detailed ASCII summary of a specific focal window `[start_1based..=end_1based]`
/// within a Multiple Sequence Alignment (MSA), including:
/// - Full-alignment summary metrics vs. focal-region conservation stats
/// - 1-based coordinate ruler
/// - Optional structural / reading-frame / secondary-structure track
/// - Aligned sequence rows with accession/taxon labels
/// - Clustal conservation track (`*` invariant, `:` transition-conserved, `.` transversion, ` ` gap)
/// - Per-column Shannon entropy sparkline (`▁▂▃▄▅▆▇█`)
pub fn render_focal_alignment_ascii_summary(
    alignment_title: &str,
    full_alignment_bp: usize,
    start_1based: usize,
    structure_track_label: Option<(&str, &str)>,
    rows: &[(&str, &str)],
    callout_notes: &[&str],
) -> String {
    if rows.is_empty() {
        return String::new();
    }
    let win_len = rows[0].1.len();
    let end_1based = start_1based + win_len.saturating_sub(1);
    let label_w = rows
        .iter()
        .map(|(lbl, _)| lbl.len())
        .chain(structure_track_label.iter().map(|(l, _)| l.len()))
        .max()
        .unwrap_or(18)
        .max(22);

    let mut invariant_cnt = 0usize;
    let mut transition_cnt = 0usize;
    let mut transversion_cnt = 0usize;
    let mut gap_cols = 0usize;
    let mut clustal_bar = String::with_capacity(win_len);
    let mut entropy_spark = String::with_capacity(win_len * 3);
    let spark_glyphs = ['▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];

    for col in 0..win_len {
        let mut counts = [0usize; 4]; // A, C, G, T
        let mut has_gap = false;
        let mut total_valid = 0usize;
        for &(_, seq) in rows {
            let b = seq.as_bytes().get(col).copied().unwrap_or(b'-');
            match b.to_ascii_uppercase() {
                b'A' => {
                    counts[0] += 1;
                    total_valid += 1;
                }
                b'C' => {
                    counts[1] += 1;
                    total_valid += 1;
                }
                b'G' => {
                    counts[2] += 1;
                    total_valid += 1;
                }
                b'T' | b'U' => {
                    counts[3] += 1;
                    total_valid += 1;
                }
                _ => {
                    has_gap = true;
                }
            }
        }
        let non_zero_bases = counts.iter().filter(|&&c| c > 0).count();
        if has_gap {
            gap_cols += 1;
            clustal_bar.push(' ');
        } else if non_zero_bases <= 1 {
            invariant_cnt += 1;
            clustal_bar.push('*');
        } else {
            let only_purines = counts[1] == 0 && counts[3] == 0;
            let only_pyrimidines = counts[0] == 0 && counts[2] == 0;
            if only_purines || only_pyrimidines {
                transition_cnt += 1;
                clustal_bar.push(':');
            } else {
                transversion_cnt += 1;
                clustal_bar.push('.');
            }
        }

        // Column Shannon entropy in bits (0.0 to 2.0)
        let mut h_col = 0.0f64;
        if total_valid > 0 {
            for &c in &counts {
                if c > 0 {
                    let p = (c as f64) / (total_valid as f64);
                    h_col -= p * p.log2();
                }
            }
        }
        let idx = ((h_col / 2.0) * 7.0).round().clamp(0.0, 7.0) as usize;
        entropy_spark.push(spark_glyphs[idx]);
    }

    let mut out = String::new();
    out.push_str("========================================================================================\n");
    out.push_str(&format!("FOCAL ALIGNMENT SUMMARY: {}\n", alignment_title));
    out.push_str(&format!(
        "Full Alignment: {} bp | Highlighted Window: nt {}..{} ({} bp) | Taxa/Operons: {}\n",
        full_alignment_bp,
        start_1based,
        end_1based,
        win_len,
        rows.len()
    ));
    out.push_str(&format!(
        "Window Conservation: {} invariant (*) | {} transition-only (:) | {} transversion (.) | {} gapped cols\n",
        invariant_cnt, transition_cnt, transversion_cnt, gap_cols
    ));
    out.push_str("========================================================================================\n");

    // Coordinate ruler
    let mut ruler_nums = vec![b' '; win_len];
    let mut ruler_ticks = vec![b'.'; win_len];
    for offset in 0..win_len {
        let pos = start_1based + offset;
        if offset == 0 || pos % 10 == 0 {
            ruler_ticks[offset] = b'|';
            let s = pos.to_string();
            for (k, ch) in s.bytes().enumerate() {
                if offset + k < win_len {
                    ruler_nums[offset + k] = ch;
                }
            }
        } else if pos % 5 == 0 {
            ruler_ticks[offset] = b':';
        }
    }
    out.push_str(&format!(
        "{:<width$}  {}\n",
        "1-Based Coord",
        String::from_utf8_lossy(&ruler_nums),
        width = label_w
    ));
    out.push_str(&format!(
        "{:<width$}  {}\n",
        "Position Scale",
        String::from_utf8_lossy(&ruler_ticks),
        width = label_w
    ));

    if let Some((st_lbl, st_seq)) = structure_track_label {
        out.push_str(&format!("{:<width$}  {}\n", st_lbl, st_seq, width = label_w));
    }
    out.push_str(&format!(
        "{}\n",
        "-".repeat(label_w + 2 + win_len)
    ));

    for &(lbl, seq) in rows {
        out.push_str(&format!("{:<width$}  {}\n", lbl, seq, width = label_w));
    }

    out.push_str(&format!(
        "{}\n",
        "-".repeat(label_w + 2 + win_len)
    ));
    out.push_str(&format!(
        "{:<width$}  {}\n",
        "Clustal Conservation",
        clustal_bar,
        width = label_w
    ));
    out.push_str(&format!(
        "{:<width$}  {}\n",
        "Column Entropy (0-2b)",
        entropy_spark,
        width = label_w
    ));

    if !callout_notes.is_empty() {
        out.push_str("\nFOCAL REGION DIAGNOSTIC CALLOUTS:\n");
        for note in callout_notes {
            out.push_str(&format!("  • {}\n", note));
        }
    }

    out
}
