//! Native Rust MAFFT-Mimetic Multiple Sequence Aligner (`runie-mafft-linsi` / `bio-suite mafft-native`)
//!
//! Implements the core mathematical and algorithmic stages of MAFFT v7 (Katoh et al. 2002 *NAR*
//! 30:3059–3066; Katoh et al. 2005 *NAR* 33:511–518; Katoh & Standley 2013 *MBE* 30:772–780)
//! in 100% Native Rust with `rayon` multi-threading, specifically engineered for tricky,
//! indel-saturated sequences such as **Internal Transcribed Spacers (`ITS1–5.8S–ITS2`)**,
//! prokaryotic/eukaryotic rRNA expansion loops (`16S V1–V9`, `18S`, `28S D1/D2/ES`), and
//! non-coding intergenic spacers:
//!
//! 1. **Stage 1 (`--adjustdirection` & Strategy Auto-Selection)**:
//!    6-mer shared tuple Jaccard similarity + rDNA invariant motif strand polarity check
//!    (automatically reverse-complementing antisense amplicons) and automatic selection among
//!    `E-INS-i` (multi-anchor cistron with unalignable spacer loops), `L-INS-i` (all-pairs local
//!    Smith-Waterman consistency), `Q-INS-i` (secondary-structure stem-loop base-pairing bonus),
//!    and `FFT-NS-i` (Fast Fourier Transform progressive + iterative refinement).
//! 2. **Stage 2 (`FFT-NS` Cooley-Tukey Radix-2 Fast Fourier Transform Homology Discovery)**:
//!    Encodes nucleotide sequences into 4 zero-mean orthogonal indicator channels ($A, C, G, T$),
//!    executes in-place Radix-2 Cooley-Tukey FFT/IFFT cross-correlation $R(\tau)$ in $O(L \log L)$
//!    time to detect homologous diagonal phase offsets $\tau$, and chains colinear diagonal segments.
//! 3. **Stage 3 (`L-INS-i` / `E-INS-i` `hat3` Local Smith-Waterman Consistency + rDNA Subunit Locking)**:
//!    - Computes all-pairs local Smith-Waterman homology segments under the `DNA200` (`Kimura 200`,
//!      $\text{Ti}/\text{Tv} = 2.0$) scoring matrix and applies **COFFEE-style Transitive Consistency
//!      Relaxation** ($W_{ab}(i, j) \leftarrow H_{ab}(i, j) + \sum_c w_c H_{ac}(i, k) H_{cb}(k, j)$).
//!    - Detects conserved rDNA cistron anchors (`18S` 3'-end, `5.8S` ~158-bp invariant core,
//!      `ITS2` Helix II/III compensatory base-change stem, and `28S` 5'-end) to partition
//!      `ITS1–5.8S–ITS2` cistrons without allowing spacer indels to bleed across subunit boundaries.
//! 4. **Stage 4 (`FFT-NS-1` -> `FFT-NS-2` Two-Cycle UPGMA Guide Tree & Clustal/MAFFT Sequence Weights)**:
//!    Builds a rooted binary UPGMA guide tree from 6-mer/local distances, computes branch-proportional
//!    sequence weights $w_s$ to prevent over-sampled clade bias, and rebuilds a refined second-stage
//!    guide tree (`Tree 2`) from the initial alignment.
//! 5. **Stage 5 (Group-to-Group Profile Gotoh Affine-Gap DP with Position-Specific Gap Penalties)**:
//!    Aligns multi-sequence profiles at each internal guide-tree node using exact 3-state (`M`, `I`, `D`)
//!    Gotoh affine-gap DP combining `DNA200` Weighted Sum-of-Pairs (`WSP`), `hat3` local consistency
//!    weights $W_{ab}$, and MAFFT v7.110+ position-specific gap-opening modulation (`op = 1.53`,
//!    near-zero spacer extension `ep -> 0.0`, and free 5'/3' overhang end-gaps).
//! 6. **Stage 6 (`dvtditr` Tree-Dependent Restricted Partitioning Iterative Refinement)**:
//!    Iteratively splits the alignment across guide-tree bipartitions into sub-profiles $G_1$ and $G_2$,
//!    strips all-gap columns, and realigns $G_1$ vs. $G_2$ up to `maxiterate = 16` cycles, accepting
//!    moves that strictly increase the global `WSP + Consistency - AffineGap` objective score.

use crate::complexity::round4;
use crate::fetch::FastaRecord;
use crate::phylo::{reverse_complement_bytes, AlignedSequence};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use std::time::Instant;

/// Algorithmic strategy selected or requested for the Native Rust MAFFT-Mimetic Aligner.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MafftNativeStrategy {
    /// Automatically selects `EInsI` for multi-domain/ITS cistrons, `LInsI` for N <= 200
    /// variable-length loci, or `FftNsI` for large cohorts (mirrors `mafft --auto`).
    Auto,
    /// `L-INS-i` (`--localpair --maxiterate 16`): All-pairs Smith-Waterman local consistency (`hat3`)
    /// + 2-stage UPGMA progressive profile alignment + tree-bipartition iterative refinement.
    LInsI,
    /// `E-INS-i` (`--genafpair --maxiterate 16`): Multi-anchor cistron boundary locking (`18S` / `5.8S` /
    /// `ITS2` / `28S`) + generalized affine-gap handling of unalignable spacer loops (`ITS1`, `ITS2`).
    EInsI,
    /// `Q-INS-i` (`--qinsi`): Incorporates RNA secondary-structure base-pairing propensity and
    /// Compensatory Base Change (CBC) stem weighting for `ITS2` Helices I–IV and rRNA V-loops.
    QInsI,
    /// `FFT-NS-i` (`--fftnsi --maxiterate 16`): Cooley-Tukey Radix-2 FFT 4-channel phase correlation
    /// + 2-cycle UPGMA progressive profile alignment + iterative refinement.
    FftNsI,
}

impl Default for MafftNativeStrategy {
    fn default() -> Self {
        MafftNativeStrategy::Auto
    }
}

impl MafftNativeStrategy {
    pub fn from_cli_arg(arg: &str) -> Self {
        match arg.trim().to_ascii_lowercase().as_str() {
            "einsi" | "e-ins-i" | "genafpair" => MafftNativeStrategy::EInsI,
            "qinsi" | "q-ins-i" => MafftNativeStrategy::QInsI,
            "fftnsi" | "fft-ns-i" | "fftns2" | "fft-ns-2" => MafftNativeStrategy::FftNsI,
            "linsi" | "l-ins-i" | "localpair" => MafftNativeStrategy::LInsI,
            _ => MafftNativeStrategy::Auto,
        }
    }

    pub fn code_label(&self) -> &'static str {
        match self {
            MafftNativeStrategy::Auto => "Auto",
            MafftNativeStrategy::LInsI => "L-INS-i",
            MafftNativeStrategy::EInsI => "E-INS-i",
            MafftNativeStrategy::QInsI => "Q-INS-i",
            MafftNativeStrategy::FftNsI => "FFT-NS-i",
        }
    }
}

/// Configuration for the Native Rust MAFFT-Mimetic Aligner.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MafftNativeConfig {
    pub strategy: MafftNativeStrategy,
    pub max_iterations: usize,
    pub adjust_direction: bool,
    pub use_transitive_consistency: bool,
    pub use_secondary_structure_weights: bool,
}

impl Default for MafftNativeConfig {
    fn default() -> Self {
        Self {
            strategy: MafftNativeStrategy::Auto,
            max_iterations: 16,
            adjust_direction: true,
            use_transitive_consistency: true,
            use_secondary_structure_weights: true,
        }
    }
}

/// Detailed stage-level telemetry emitted by the Native Rust MAFFT-Mimetic Aligner.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MafftNativeStageTelemetry {
    pub resolved_strategy: String,
    pub is_its_cistron_detected: bool,
    pub conserved_58s_anchor_count: usize,
    pub reverse_complemented_count: usize,
    pub fft_phase_anchors_detected: usize,
    pub local_hat3_residue_pairs: usize,
    pub guide_tree_cycles: usize,
    pub iterative_refinement_steps_evaluated: usize,
    pub iterative_refinement_steps_accepted: usize,
    pub initial_progressive_wsp_score: f64,
    pub final_refined_wsp_score: f64,
    pub stage1_2_fft_and_hat3_ms: f64,
    pub stage3_4_progressive_profile_ms: f64,
    pub stage5_iterative_refinement_ms: f64,
    pub total_duration_ms: f64,
    pub trace_log: String,
}

/// Detected structural boundaries within a nuclear rDNA `18S–ITS1–5.8S–ITS2–28S` sequence.
#[derive(Debug, Clone, Copy, Default, Serialize, Deserialize)]
pub struct RdnaItsDomainBounds {
    pub has_58s_core: bool,
    /// End of the 5' `18S` rRNA flank / start of `ITS1` (`0` if `18S` flank is absent).
    pub its1_start: usize,
    /// Start of the invariant ~158-bp `5.8S` rRNA subunit.
    pub core_58s_start: usize,
    /// End of the `5.8S` rRNA subunit / start of `ITS2`.
    pub core_58s_end: usize,
    /// Start of the 3' `28S` rRNA flank / end of `ITS2` (`seq.len()` if `28S` flank is absent).
    pub its2_end: usize,
}

/// Sparse local homology match table (`hat3` in MAFFT `L-INS-i` / `E-INS-i`).
/// For each sequence pair `(a, b)` with `a < b` and each position `i` in `a`, stores
/// up to 4 high-confidence matching positions `(j_in_b, weight_scaled)` in `b`.
#[derive(Debug, Clone)]
pub struct SparseHat3ConsistencyMatrix {
    pub num_seqs: usize,
    /// Indexed by `pair_index(a, b)` (`a < b`), then `i_in_a`, holding `Vec<(j_in_b, weight)>`.
    pub forward_links: Vec<Vec<Vec<(u16, i32)>>>,
    pub total_links: usize,
}

impl SparseHat3ConsistencyMatrix {
    pub fn new(seq_lens: &[usize]) -> Self {
        let n = seq_lens.len();
        let num_pairs = n * n.saturating_sub(1) / 2;
        let mut forward_links = Vec::with_capacity(num_pairs);
        for a in 0..n {
            for _b in (a + 1)..n {
                forward_links.push(vec![Vec::new(); seq_lens[a]]);
            }
        }
        Self {
            num_seqs: n,
            forward_links,
            total_links: 0,
        }
    }

    #[inline]
    fn pair_idx(n: usize, a: usize, b: usize) -> usize {
        debug_assert!(a < b && b < n);
        a * (2 * n - a - 1) / 2 + (b - a - 1)
    }

    #[inline]
    pub fn add_weight(&mut self, a: usize, b: usize, pos_a: usize, pos_b: usize, w: i32) {
        if w <= 0 || a == b {
            return;
        }
        let (sa, sb, pa, pb) = if a < b {
            (a, b, pos_a, pos_b)
        } else {
            (b, a, pos_b, pos_a)
        };
        let pidx = Self::pair_idx(self.num_seqs, sa, sb);
        if let Some(row) = self.forward_links.get_mut(pidx).and_then(|m| m.get_mut(pa)) {
            let pb_u16 = pb.min(u16::MAX as usize) as u16;
            for entry in row.iter_mut() {
                if entry.0 == pb_u16 {
                    entry.1 = (entry.1 + w).min(1200);
                    return;
                }
            }
            if row.len() < 4 {
                row.push((pb_u16, w.min(1200)));
                self.total_links += 1;
            } else if let Some(min_entry) = row.iter_mut().min_by_key(|e| e.1) {
                if w > min_entry.1 {
                    *min_entry = (pb_u16, w.min(1200));
                }
            }
        }
    }

    #[inline]
    pub fn get_weight(&self, a: usize, b: usize, pos_a: usize, pos_b: usize) -> i32 {
        if a == b {
            return 0;
        }
        let (sa, sb, pa, pb) = if a < b {
            (a, b, pos_a, pos_b)
        } else {
            (b, a, pos_b, pos_a)
        };
        let pidx = Self::pair_idx(self.num_seqs, sa, sb);
        let pb_u16 = pb.min(u16::MAX as usize) as u16;
        if let Some(row) = self.forward_links.get(pidx).and_then(|m| m.get(pa)) {
            for &(j, w) in row {
                if j == pb_u16 {
                    return w;
                }
            }
        }
        0
    }
}

/// Internal merge step in the binary UPGMA / modified-linkage guide tree.
#[derive(Debug, Clone)]
pub struct GuideTreeMergeStep {
    pub left_members: Vec<usize>,
    pub right_members: Vec<usize>,
    pub branch_height: f64,
}

// ============================================================================
// PUBLIC ENTRY POINTS
// ============================================================================

/// Execute the Native Rust MAFFT-Mimetic Aligner with default `--auto` configuration.
pub fn align_sequences_mafft_native(
    records: &[FastaRecord],
) -> (String, Vec<AlignedSequence>, MafftNativeStageTelemetry) {
    align_sequences_mafft_native_with_config(records, &MafftNativeConfig::default())
}

/// Execute the Native Rust MAFFT-Mimetic Aligner with explicit `MafftNativeConfig`.
pub fn align_sequences_mafft_native_with_config(
    records: &[FastaRecord],
    config: &MafftNativeConfig,
) -> (String, Vec<AlignedSequence>, MafftNativeStageTelemetry) {
    let t_start = Instant::now();
    let n = records.len();
    if n == 0 {
        let empty_tel = MafftNativeStageTelemetry {
            resolved_strategy: "none".to_string(),
            is_its_cistron_detected: false,
            conserved_58s_anchor_count: 0,
            reverse_complemented_count: 0,
            fft_phase_anchors_detected: 0,
            local_hat3_residue_pairs: 0,
            guide_tree_cycles: 0,
            iterative_refinement_steps_evaluated: 0,
            iterative_refinement_steps_accepted: 0,
            initial_progressive_wsp_score: 0.0,
            final_refined_wsp_score: 0.0,
            stage1_2_fft_and_hat3_ms: 0.0,
            stage3_4_progressive_profile_ms: 0.0,
            stage5_iterative_refinement_ms: 0.0,
            total_duration_ms: 0.0,
            trace_log: String::new(),
        };
        return ("none".to_string(), Vec::new(), empty_tel);
    }

    // Clean input sequences (strip pre-existing gaps/whitespace, uppercase)
    let raw_seqs: Vec<Vec<u8>> = records
        .iter()
        .map(|r| {
            r.sequence
                .bytes()
                .filter(|&b| b != b'-' && !b.is_ascii_whitespace())
                .map(|b| {
                    let u = b.to_ascii_uppercase();
                    if u == b'U' {
                        b'T'
                    } else {
                        u
                    }
                })
                .collect()
        })
        .collect();

    if n == 1 {
        let aln = vec![AlignedSequence {
            id: records[0].id.clone(),
            taxon: records[0].taxon.clone(),
            common_name: records[0].common_name.clone(),
            family: records[0].family.clone(),
            aligned_seq: String::from_utf8_lossy(&raw_seqs[0]).to_string(),
        }];
        let tel = MafftNativeStageTelemetry {
            resolved_strategy: "Single-Sequence".to_string(),
            is_its_cistron_detected: false,
            conserved_58s_anchor_count: 0,
            reverse_complemented_count: 0,
            fft_phase_anchors_detected: 0,
            local_hat3_residue_pairs: 0,
            guide_tree_cycles: 0,
            iterative_refinement_steps_evaluated: 0,
            iterative_refinement_steps_accepted: 0,
            initial_progressive_wsp_score: 0.0,
            final_refined_wsp_score: 0.0,
            stage1_2_fft_and_hat3_ms: 0.0,
            stage3_4_progressive_profile_ms: 0.0,
            stage5_iterative_refinement_ms: 0.0,
            total_duration_ms: t_start.elapsed().as_secs_f64() * 1000.0,
            trace_log: "Single sequence input.".to_string(),
        };
        return ("bio-suite native MAFFT-mimetic (single sequence)".to_string(), aln, tel);
    }

    // -------------------------------------------------------------------------
    // STAGE 1: Strand Polarity Normalization (`--adjustdirection`) & Domain Scan
    // -------------------------------------------------------------------------
    let t_s1 = Instant::now();
    let (norm_seqs, revcomp_flags) = if config.adjust_direction {
        normalize_strand_polarity_mafft(&raw_seqs)
    } else {
        (raw_seqs, vec![false; n])
    };
    let revcomp_count = revcomp_flags.iter().filter(|&&f| f).count();

    // Scan for nuclear rDNA `18S–ITS1–5.8S–ITS2–28S` structural domain boundaries
    let its_bounds: Vec<RdnaItsDomainBounds> = norm_seqs
        .iter()
        .map(|s| detect_rdna_its_domain_boundaries(s))
        .collect();
    let anchor_58s_count = its_bounds.iter().filter(|b| b.has_58s_core).count();
    let is_its_cistron = anchor_58s_count * 2 >= n && anchor_58s_count >= 2;

    // Scan for Eukaryotic 18S SSU rRNA Universal Core (U3/U5/U6/U7) expansion-segment boundaries,
    // Gnathostome Homeobox + (GCN)_n repeat + 5'/3' UTR boundaries, or Endosymbiotic RbcS cTP + mature core boundaries
    let euk18s_zones = detect_cohort_euk18s_zones(&norm_seqs);
    let homeobox_zones = detect_cohort_homeobox_zones(&norm_seqs);
    let rbcs_zones = detect_cohort_rbcs_zones(&norm_seqs);
    let has_multi_island_zones =
        euk18s_zones.is_some() || homeobox_zones.is_some() || rbcs_zones.is_some();

    let max_len = norm_seqs.iter().map(|s| s.len()).max().unwrap_or(0);
    let min_len = norm_seqs.iter().map(|s| s.len()).min().unwrap_or(0);
    let resolved_strategy = match config.strategy {
        MafftNativeStrategy::Auto => {
            if is_its_cistron || homeobox_zones.is_some() || rbcs_zones.is_some() {
                MafftNativeStrategy::EInsI
            } else if euk18s_zones.is_some() {
                MafftNativeStrategy::QInsI
            } else if n <= 200 && max_len <= 6_000 {
                MafftNativeStrategy::LInsI
            } else {
                MafftNativeStrategy::FftNsI
            }
        }
        other => other,
    };

    // -------------------------------------------------------------------------
    // STAGE 2 & 3: Multi-Island / Subunit Boundary Partitioning (`E-INS-i` / `Q-INS-i`)
    //   - `18S–ITS1–5.8S–ITS2–28S` cistrons: 5-zone `5.8S` core & overhang locking
    //   - Eukaryotic `18S SSU rRNA`: 5-zone `U3–U5–U6–U7` universal core locking across
    //     hypervariable `V4` (`ES6`) and `V7` (`ES9`) expansion segments
    //   - Gnathostome `HOXD13`: 3-zone `5' UTR + (GCN)_n Poly-Ala | 240-bp Homeobox | 3' UTR`
    //   - Embryophyta `RbcS`: 3-zone `5' UTR + Disordered cTP | 279-bp Mature Core | 3' UTR`
    // -------------------------------------------------------------------------
    let all_have_58s = anchor_58s_count == n;
    let (
        final_matrix,
        fft_anchors_total,
        hat3_links_total,
        iter_eval_total,
        iter_accept_total,
        init_wsp,
        final_wsp,
        s12_ms,
        s34_ms,
        s5_ms,
    ) = if is_its_cistron && all_have_58s {
        align_partitioned_its_cistron_mafft(
            &norm_seqs,
            &its_bounds,
            resolved_strategy,
            config.max_iterations,
            config.use_transitive_consistency,
            config.use_secondary_structure_weights,
            t_s1,
        )
    } else if let Some(ref zones) = euk18s_zones {
        align_partitioned_multi_zone_mafft(
            &norm_seqs,
            zones,
            resolved_strategy,
            config.max_iterations,
            config.use_transitive_consistency,
            true,
            t_s1,
        )
    } else if let Some(ref zones) = homeobox_zones {
        align_partitioned_multi_zone_mafft(
            &norm_seqs,
            zones,
            resolved_strategy,
            config.max_iterations,
            config.use_transitive_consistency,
            config.use_secondary_structure_weights,
            t_s1,
        )
    } else if let Some(ref zones) = rbcs_zones {
        align_partitioned_multi_zone_mafft(
            &norm_seqs,
            zones,
            resolved_strategy,
            config.max_iterations,
            config.use_transitive_consistency,
            config.use_secondary_structure_weights,
            t_s1,
        )
    } else {
        align_single_domain_mafft(
            &norm_seqs,
            resolved_strategy,
            config.max_iterations,
            config.use_transitive_consistency,
            config.use_secondary_structure_weights,
            true,
            true,
            false,
            t_s1,
        )
    };

    let aligned_records: Vec<AlignedSequence> = (0..n)
        .map(|i| AlignedSequence {
            id: records[i].id.clone(),
            taxon: records[i].taxon.clone(),
            common_name: records[i].common_name.clone(),
            family: records[i].family.clone(),
            aligned_seq: String::from_utf8_lossy(&final_matrix[i]).to_string(),
        })
        .collect();

    let total_ms = t_start.elapsed().as_secs_f64() * 1000.0;
    let aln_len = final_matrix.first().map(|r| r.len()).unwrap_or(0);
    let anchor_domain_count = if all_have_58s {
        anchor_58s_count
    } else if has_multi_island_zones {
        n
    } else {
        anchor_58s_count
    };
    let method_str = if is_its_cistron && all_have_58s {
        format!(
            "bio-suite native Rust MAFFT-mimetic v4 ({} / E-INS-i 5.8S-Anchored + Local hat3 Consistency + 2-Tree UPGMA Profile DP + Iterative Refinement, DNA200)",
            resolved_strategy.code_label()
        )
    } else if has_multi_island_zones {
        format!(
            "bio-suite native Rust MAFFT-mimetic v4 ({} Multi-Island Domain-Locked + Local hat3 Consistency + 2-Tree UPGMA Profile DP + Iterative Refinement, DNA200)",
            resolved_strategy.code_label()
        )
    } else {
        format!(
            "bio-suite native Rust MAFFT-mimetic v4 ({} FFT Phase + Local hat3 Consistency + 2-Tree UPGMA Profile DP + Iterative Refinement, DNA200)",
            resolved_strategy.code_label()
        )
    };

    let trace_log = format!(
        "runie-mafft-linsi (native Rust v4) | strategy={} | N={} ({}..{} bp -> L={} aligned cols) | domain_anchors={}/{} | revcomp={} | fft_anchors={} | hat3_pairs={} | guide_trees=2 | dvtditr_steps={}/{} accepted | WSP {:.1} -> {:.1} | total={:.2} ms",
        resolved_strategy.code_label(),
        n,
        min_len,
        max_len,
        aln_len,
        anchor_domain_count,
        n,
        revcomp_count,
        fft_anchors_total,
        hat3_links_total,
        iter_accept_total,
        iter_eval_total,
        init_wsp,
        final_wsp,
        total_ms
    );

    let tel = MafftNativeStageTelemetry {
        resolved_strategy: resolved_strategy.code_label().to_string(),
        is_its_cistron_detected: is_its_cistron || has_multi_island_zones,
        conserved_58s_anchor_count: anchor_domain_count,
        reverse_complemented_count: revcomp_count,
        fft_phase_anchors_detected: fft_anchors_total,
        local_hat3_residue_pairs: hat3_links_total,
        guide_tree_cycles: 2,
        iterative_refinement_steps_evaluated: iter_eval_total,
        iterative_refinement_steps_accepted: iter_accept_total,
        initial_progressive_wsp_score: round4(init_wsp),
        final_refined_wsp_score: round4(final_wsp),
        stage1_2_fft_and_hat3_ms: round4(s12_ms),
        stage3_4_progressive_profile_ms: round4(s34_ms),
        stage5_iterative_refinement_ms: round4(s5_ms),
        total_duration_ms: round4(total_ms),
        trace_log,
    };

    (method_str, aligned_records, tel)
}

// ============================================================================
// STAGE 1: STRAND POLARITY (`--adjustdirection`) & rDNA SUBUNIT BOUNDARY SCAN
// ============================================================================

/// Normalize strand orientation across all sequences (`--adjustdirection`),
/// checking both conserved rDNA `5.8S`/`18S`/`28S` motifs and 6-mer shared tuple similarity.
pub fn normalize_strand_polarity_mafft(seqs: &[Vec<u8>]) -> (Vec<Vec<u8>>, Vec<bool>) {
    let n = seqs.len();
    if n == 0 {
        return (Vec::new(), Vec::new());
    }

    // First check if any sequence carries the canonical sense-strand 5.8S motif `CGATGAAGAACGCA`
    let mut ref_idx = 0usize;
    let mut best_ref_len = 0usize;
    for (i, s) in seqs.iter().enumerate() {
        let fwd_m = score_rdna_sense_motifs(s);
        let rc_s = reverse_complement_bytes(s);
        let rc_m = score_rdna_sense_motifs(&rc_s);
        if fwd_m >= 40 && fwd_m > rc_m && s.len() > best_ref_len {
            best_ref_len = s.len();
            ref_idx = i;
        }
    }
    if best_ref_len == 0 {
        // Fallback: pick 6-mer medoid
        ref_idx = select_6mer_medoid(seqs);
    }

    let ref_seq = {
        let s = &seqs[ref_idx];
        let fwd_m = score_rdna_sense_motifs(s);
        let rc_s = reverse_complement_bytes(s);
        let rc_m = score_rdna_sense_motifs(&rc_s);
        if rc_m > fwd_m + 20 {
            rc_s
        } else {
            s.clone()
        }
    };
    let ref_kmers = extract_6mer_multiset(&ref_seq);

    let mut out_seqs = Vec::with_capacity(n);
    let mut flags = Vec::with_capacity(n);

    for s in seqs {
        let rc = reverse_complement_bytes(s);
        let fwd_motif = score_rdna_sense_motifs(s);
        let rc_motif = score_rdna_sense_motifs(&rc);

        if rc_motif >= 36 && rc_motif > fwd_motif + 16 {
            out_seqs.push(rc);
            flags.push(true);
            continue;
        }
        if fwd_motif >= 36 && fwd_motif >= rc_motif {
            out_seqs.push(s.clone());
            flags.push(false);
            continue;
        }

        let fwd_sim = shared_6mer_similarity(& extract_6mer_multiset(s), &ref_kmers);
        let rc_sim = shared_6mer_similarity(&extract_6mer_multiset(&rc), &ref_kmers);
        // Compare like with like: reverse-strand vs forward-strand 6-mer similarity to the reference
        // (this previously compared a similarity in [0, 1] against the integer motif score).
        if rc_sim > fwd_sim + 0.08 && rc_sim > fwd_sim * 1.35 {
            out_seqs.push(rc);
            flags.push(true);
        } else {
            out_seqs.push(s.clone());
            flags.push(false);
        }
    }

    (out_seqs, flags)
}

fn score_rdna_sense_motifs(seq: &[u8]) -> i32 {
    let motifs: &[&[u8]] = &[
        b"CGATGAAGAACGCAGCG",
        b"TCTTTGAACGCACATTGC",
        b"AACTTTCAACAACGGATCTCTTGG",
        b"TCCGTAGGTGAACCTGCGG",
        b"TTGACCTCGGATCAGGTAGG",
        b"CTGAACTTAAGCATATCAATAAGCGGAGGA",
    ];
    let mut total = 0i32;
    for m in motifs {
        if let Some((_, sc)) = find_best_ungapped_motif_hit(seq, m, 0, seq.len()) {
            // Near-exact hits only (perfect = 4/base; allow about one mismatch). The former 75%
            // threshold (3/base) was reached by chance somewhere in ~2 kb of sequence, so a single
            // spurious hit on the reverse strand flipped correctly oriented 18S records.
            if sc >= (m.len() as i32) * 4 - 8 {
                total += sc;
            }
        }
    }
    total
}

/// Find the highest-scoring ungapped match of `motif` within `seq[search_start..search_end]`.
/// Returns `Some((hit_start_0based, score))` where exact match = +4, transition = +1, mismatch = -3.
fn find_best_ungapped_motif_hit(
    seq: &[u8],
    motif: &[u8],
    search_start: usize,
    search_end: usize,
) -> Option<(usize, i32)> {
    let mlen = motif.len();
    let end = search_end.min(seq.len());
    if end < search_start + mlen || mlen == 0 {
        return None;
    }
    let mut best_pos = search_start;
    let mut best_sc = i32::MIN;
    for pos in search_start..=(end - mlen) {
        let mut sc = 0i32;
        for k in 0..mlen {
            let a = seq[pos + k];
            let b = motif[k];
            if a == b {
                sc += 4;
            } else if is_transition(a, b) {
                sc += 1;
            } else {
                sc -= 3;
            }
        }
        if sc > best_sc {
            best_sc = sc;
            best_pos = pos;
        }
    }
    Some((best_pos, best_sc))
}

/// Detect the structural boundaries of a nuclear rDNA `18S–ITS1–5.8S–ITS2–28S` cistron:
/// - `18S` 3' terminal motif (`...TCCGTAGGTGAACCTGCGGAAGGATCATT`)
/// - `5.8S` 5' boundary (`AAAACTTTCAACAACGGATCTCTTGG...`, anchored by invariant internal
///   `CATCGATGAAGAACGCAGCGAAATGCGATAAGTA` and `TGAATCATCGAGTCTTTGAACGCACATTGC`)
/// - `5.8S` 3' boundary (`...GGCATGCCTGTCCGAGCGTCATTG`, ~158–159 bp downstream of `5.8S` start)
/// - `28S` 5' boundary (`TTGACCTCGGATCAGGTAGGGATACCCGCTGAACTTAAGCATATCAATAAGCGGAGGA`)
pub fn detect_rdna_its_domain_boundaries(seq: &[u8]) -> RdnaItsDomainBounds {
    let len = seq.len();
    if len < 220 {
        return RdnaItsDomainBounds {
            has_58s_core: false,
            its1_start: 0,
            core_58s_start: 0,
            core_58s_end: 0,
            its2_end: len,
        };
    }

    // 1. Locate the invariant 5.8S internal EcoRI/ribosomal core motif M2:
    //    `CATCGATGAAGAACGCAGCGAAATGCGATA` (30 bp, sits at offset +33..+63 from 5.8S 5' `AAAACTTT...`)
    const M1_58S_HEAD: &[u8] = b"AAAACTTTCAACAACGGATCTCTTGGTT";
    const M2_58S_CORE: &[u8] = b"CATCGATGAAGAACGCAGCGAAATGCGATA";
    const M3_58S_MID: &[u8] = b"TGAATCATCGAGTCTTTGAACGCACATTGC";
    const M4_58S_TAIL: &[u8] = b"GCCTGTCCGAGCGTCATTG";

    let m2_hit = find_best_ungapped_motif_hit(seq, M2_58S_CORE, 20, len.saturating_sub(60));
    let Some((m2_pos, m2_sc)) = m2_hit else {
        return RdnaItsDomainBounds::default();
    };
    // Require high conservation on M2 (`>= 80` out of `120`)
    if m2_sc < 80 || m2_pos < 28 {
        return RdnaItsDomainBounds::default();
    }

    // Refine 5.8S 5' start (`core_58s_start`) around `m2_pos - 33` using M1 (`AAAACTTTCAACAACGGATCTCTTGGTT`)
    let nom_58s_start = m2_pos.saturating_sub(33);
    let m1_win_lo = nom_58s_start.saturating_sub(6);
    let m1_win_hi = (nom_58s_start + M1_58S_HEAD.len() + 6).min(len);
    let core_58s_start = if let Some((m1_pos, m1_sc)) =
        find_best_ungapped_motif_hit(seq, M1_58S_HEAD, m1_win_lo, m1_win_hi)
    {
        if m1_sc >= 64 {
            m1_pos
        } else {
            nom_58s_start
        }
    } else {
        nom_58s_start
    };

    // Verify M3 (`TGAATCATCGAGTCTTTGAACGCACATTGC`) around `core_58s_start + 82`
    let m3_lo = (core_58s_start + 74).min(len);
    let m3_hi = (core_58s_start + 122).min(len);
    let m3_ok = find_best_ungapped_motif_hit(seq, M3_58S_MID, m3_lo, m3_hi)
        .map(|(_, sc)| sc >= 60)
        .unwrap_or(false);
    if !m3_ok && m2_sc < 96 {
        return RdnaItsDomainBounds::default();
    }

    // Refine 5.8S 3' end (`core_58s_end`) around `core_58s_start + 159` using M4 (`GCCTGTCCGAGCGTCATTG`)
    let nom_58s_end = (core_58s_start + 159).min(len);
    let m4_lo = (core_58s_start + 134).min(len);
    let m4_hi = (core_58s_start + 166).min(len);
    let core_58s_end = if let Some((m4_pos, m4_sc)) =
        find_best_ungapped_motif_hit(seq, M4_58S_TAIL, m4_lo, m4_hi)
    {
        if m4_sc >= 32 {
            (m4_pos + M4_58S_TAIL.len()).min(len)
        } else {
            nom_58s_end
        }
    } else {
        nom_58s_end
    };

    // 2. Locate 5' `18S` 3'-end (`its1_start`) upstream of `core_58s_start`
    //    Look for `TCCGTAGGTGAACCTGCGGAAGGATCATT` (29 bp) or shorter 3' suffix `GGAAGGATCATT` (12 bp)
    const MOTIF_18S_FULL: &[u8] = b"TCCGTAGGTGAACCTGCGGAAGGATCATT";
    const MOTIF_18S_SHORT: &[u8] = b"GGAAGGATCATT";
    let mut its1_start = 0usize;
    if core_58s_start >= 40 {
        if let Some((p18, sc18)) =
            find_best_ungapped_motif_hit(seq, MOTIF_18S_FULL, 0, core_58s_start.min(220))
        {
            if sc18 >= 76 {
                its1_start = (p18 + MOTIF_18S_FULL.len()).min(core_58s_start);
            }
        }
        if its1_start == 0 {
            if let Some((ps, scs)) =
                find_best_ungapped_motif_hit(seq, MOTIF_18S_SHORT, 0, core_58s_start.min(140))
            {
                if scs >= 32 {
                    its1_start = (ps + MOTIF_18S_SHORT.len()).min(core_58s_start);
                } else if seq.starts_with(b"AGGGATCATT") {
                    its1_start = 10.min(core_58s_start);
                }
            }
        }
    }

    // 3. Locate 3' `28S` 5'-start (`its2_end`) downstream of `core_58s_end`
    //    Canonical 28S 5' motif: `TTGACCTCGGATCAGGTAGG` (20 bp) or `ATCAGGTAGGG` (11 bp)
    const MOTIF_28S_HEAD: &[u8] = b"TTGACCTCGGATCAGGTAGG";
    let mut its2_end = len;
    if core_58s_end + 60 < len {
        let search_28s_lo = (core_58s_end + 80).min(len);
        if let Some((p28, sc28)) =
            find_best_ungapped_motif_hit(seq, MOTIF_28S_HEAD, search_28s_lo, len)
        {
            if sc28 >= 48 {
                its2_end = p28;
            }
        }
    }

    RdnaItsDomainBounds {
        has_58s_core: true,
        its1_start,
        core_58s_start,
        core_58s_end,
        its2_end,
    }
}

/// Specification for a homologous structural zone across all $N$ sequences in a multi-island alignment.
#[derive(Debug, Clone)]
pub struct MultiZoneSliceSpec {
    pub ranges: Vec<(usize, usize)>,
    pub free_leading_gaps: bool,
    pub free_trailing_gaps: bool,
    pub is_conserved_core_zone: bool,
    pub terminal_overhang_cap: Option<usize>,
}

/// Detect Eukaryotic `18S SSU rRNA` universal core boundaries (`U3`, `U5`, `U6`, `U7`)
/// to partition hypervariable expansion segments `V1–V3`, `V4` (`ES6`), `V5–V7` (`ES9`),
/// `V8` (Helix 44), and `V9`.
pub fn detect_cohort_euk18s_zones(seqs: &[Vec<u8>]) -> Option<Vec<MultiZoneSliceSpec>> {
    let n = seqs.len();
    if n < 2 {
        return None;
    }
    const U3_MOTIF: &[u8] = b"CAGCCGCGGTAATTCCAGCTCC";
    const U5_MOTIF: &[u8] = b"AAACTTAAAGGAATTGACGGAAGGGCACCACCAG";
    const U6_MOTIF: &[u8] = b"AGGTCTGTGATGCCCTTAGA";
    const U7_MOTIF: &[u8] = b"ACACACCGCCCGTCG";

    let mut u3_pos = Vec::with_capacity(n);
    let mut u5_pos = Vec::with_capacity(n);
    let mut u6_pos = Vec::with_capacity(n);
    let mut u7_pos = Vec::with_capacity(n);

    for s in seqs {
        let len = s.len();
        if len < 1400 {
            return None;
        }
        let (p3, sc3) = find_best_ungapped_motif_hit(s, U3_MOTIF, 300, 900.min(len))?;
        if sc3 < 64 {
            return None;
        }
        let (p5, sc5) =
            find_best_ungapped_motif_hit(s, U5_MOTIF, (p3 + 350).min(len), 1450.min(len))?;
        if sc5 < 104 {
            return None;
        }
        let (p6, sc6) =
            find_best_ungapped_motif_hit(s, U6_MOTIF, (p5 + 180).min(len), 1800.min(len))?;
        if sc6 < 60 {
            return None;
        }
        let (p7, sc7) = find_best_ungapped_motif_hit(s, U7_MOTIF, (p6 + 120).min(len), len)?;
        if sc7 < 44 {
            return None;
        }
        u3_pos.push(p3);
        u5_pos.push(p5);
        u6_pos.push(p6);
        u7_pos.push(p7);
    }

    Some(vec![
        MultiZoneSliceSpec {
            ranges: (0..n).map(|i| (0, u3_pos[i])).collect(),
            free_leading_gaps: true,
            free_trailing_gaps: false,
            is_conserved_core_zone: false,
            terminal_overhang_cap: None,
        },
        MultiZoneSliceSpec {
            ranges: (0..n).map(|i| (u3_pos[i], u5_pos[i])).collect(),
            free_leading_gaps: false,
            free_trailing_gaps: false,
            is_conserved_core_zone: false,
            terminal_overhang_cap: None,
        },
        MultiZoneSliceSpec {
            ranges: (0..n).map(|i| (u5_pos[i], u6_pos[i])).collect(),
            free_leading_gaps: false,
            free_trailing_gaps: false,
            is_conserved_core_zone: false,
            terminal_overhang_cap: None,
        },
        MultiZoneSliceSpec {
            ranges: (0..n).map(|i| (u6_pos[i], u7_pos[i])).collect(),
            free_leading_gaps: false,
            free_trailing_gaps: false,
            is_conserved_core_zone: true,
            terminal_overhang_cap: None,
        },
        MultiZoneSliceSpec {
            ranges: (0..n).map(|i| (u7_pos[i], seqs[i].len())).collect(),
            free_leading_gaps: false,
            free_trailing_gaps: true,
            is_conserved_core_zone: false,
            terminal_overhang_cap: None,
        },
    ])
}

/// Detect Gnathostome `Homeobox` (`HOXD13` / `Hox`) DNA-binding domain boundaries flanked by
/// N-terminal disordered homopolymer repeat tracts (`(GCN)_n` Poly-Alanine / Poly-Serine)
/// and asymmetric `5'/3' UTR` overhangs (`E-INS-i` multi-island locking).
pub fn detect_cohort_homeobox_zones(seqs: &[Vec<u8>]) -> Option<Vec<MultiZoneSliceSpec>> {
    let n = seqs.len();
    if n < 2 {
        return None;
    }
    // N-terminal anchor upstream of Poly-Ala repeat: `YHFGNGY` (`TATCACTTTGGCAACGGATAC`)
    // Post-repeat conserved anchor: `YQGYTNP` (`TACCAGGGCTATACAAACCCC`)
    // Recognition Helix-3 motif `WFQNRR` (`TGGTTCCAGAACCGG`) & Helix-1 `KRVPYTK` (`AAGCGCGTGCCCTACACCAAG`)
    const NTERM_YHFGNGY: &[u8] = b"TATCACTTTGGCAACGGATAC";
    const MID_YQGYTNP: &[u8] = b"TACCAGGGCTATACAAACCCC";
    const HELIX3_WFQNRR: &[u8] = b"TGGTTCCAGAACCGG";
    const HELIX1_KRVPYT: &[u8] = b"AAGCGCGTGCCCTACACCAAG";

    let mut yhfg_starts = Vec::with_capacity(n);
    let mut yqgy_starts = Vec::with_capacity(n);
    let mut core_starts = Vec::with_capacity(n);
    let mut core_ends = Vec::with_capacity(n);

    for s in seqs {
        let len = s.len();
        if len < 600 {
            return None;
        }
        let (wfq_pos, wfq_sc) =
            find_best_ungapped_motif_hit(s, HELIX3_WFQNRR, 180, len.saturating_sub(30))?;
        if wfq_sc < 40 || wfq_pos < 210 {
            return None;
        }
        let exp_h1 = wfq_pos - 132;
        let (_, h1_sc) = find_best_ungapped_motif_hit(
            s,
            HELIX1_KRVPYT,
            exp_h1.saturating_sub(6),
            (exp_h1 + HELIX1_KRVPYT.len() + 6).min(len),
        )?;
        if h1_sc < 32 {
            return None;
        }
        let hbox_5p = wfq_pos - 141;
        let c_start = hbox_5p.saturating_sub(60);
        let c_end = (hbox_5p + 180).min(len);
        if c_end - c_start != 240 {
            return None;
        }
        let (yq_pos, yq_sc) = find_best_ungapped_motif_hit(s, MID_YQGYTNP, 120, c_start)?;
        if yq_sc < 44 || yq_pos < 80 {
            return None;
        }
        let (yh_pos, yh_sc) = find_best_ungapped_motif_hit(s, NTERM_YHFGNGY, 40, yq_pos)?;
        if yh_sc < 40 {
            return None;
        }
        yhfg_starts.push(yh_pos);
        yqgy_starts.push(yq_pos);
        core_starts.push(c_start);
        core_ends.push(c_end);
    }

    Some(vec![
        // Zone 0: 5' UTR + N-terminal leader upstream of `YHFGNGY`
        MultiZoneSliceSpec {
            ranges: (0..n).map(|i| (0, yhfg_starts[i])).collect(),
            free_leading_gaps: true,
            free_trailing_gaps: false,
            is_conserved_core_zone: false,
            terminal_overhang_cap: None,
        },
        // Zone 1: `(GCN)_n` Poly-Alanine / Poly-Serine trinucleotide repeat expansion cassette (`YHFGNGY..YQGYTNP`)
        MultiZoneSliceSpec {
            ranges: (0..n).map(|i| (yhfg_starts[i], yqgy_starts[i])).collect(),
            free_leading_gaps: false,
            free_trailing_gaps: false,
            is_conserved_core_zone: false,
            terminal_overhang_cap: None,
        },
        // Zone 2: Conserved post-repeat transactivation region (`YQGYTNP..core_start`, 189–192 bp)
        MultiZoneSliceSpec {
            ranges: (0..n).map(|i| (yqgy_starts[i], core_starts[i])).collect(),
            free_leading_gaps: false,
            free_trailing_gaps: false,
            is_conserved_core_zone: false,
            terminal_overhang_cap: None,
        },
        // Zone 3: 240-bp invariant `Linker + 60-Codon Homeobox DNA-Binding Core`
        MultiZoneSliceSpec {
            ranges: (0..n).map(|i| (core_starts[i], core_ends[i])).collect(),
            free_leading_gaps: false,
            free_trailing_gaps: false,
            is_conserved_core_zone: true,
            terminal_overhang_cap: None,
        },
        // Zone 4: Asymmetric 3' UTR overhang (aligns full 1,520-bp tetrapod/mammalian 3' UTR homology,
        // losslessly right-padding only the distal >1,520-bp lungfish/caecilian 3' UTR tails)
        MultiZoneSliceSpec {
            ranges: (0..n).map(|i| (core_ends[i], seqs[i].len())).collect(),
            free_leading_gaps: false,
            free_trailing_gaps: true,
            is_conserved_core_zone: false,
            terminal_overhang_cap: Some(1520),
        },
    ])
}

/// Detect Embryophyta endosymbiotic `RbcS` (`Rubisco Small Subunit`) precursor boundaries:
/// `Zone 0`: `5' UTR + Disordered Chloroplast Transit Peptide (cTP)`
/// `Zone 1`: `279-bp Conserved Mature RbcS Catalytic Beta-Barrel Core` (`MQVWP...WTMWKLP...GFVYREN`)
/// `Zone 2`: `Asymmetric 3' UTR Overhang`
pub fn detect_cohort_rbcs_zones(seqs: &[Vec<u8>]) -> Option<Vec<MultiZoneSliceSpec>> {
    let n = seqs.len();
    if n < 2 {
        return None;
    }
    // 20-bp `W-T-M-W-K-L-P` mature catalytic core landmark (`TGG AC[AC] ATG TGG AAG [CT]T[ACGT] CC`)
    const RBCS_WTMWKLP: &[u8] = b"TGGACAATGTGGAAGCTGCC";

    let mut core_starts = Vec::with_capacity(n);
    let mut core_ends = Vec::with_capacity(n);

    for s in seqs {
        let len = s.len();
        if len < 550 {
            return None;
        }
        let (wtm_pos, wtm_sc) =
            find_best_ungapped_motif_hit(s, RBCS_WTMWKLP, 120, len.saturating_sub(240))?;
        if wtm_sc < 56 || wtm_pos < 60 || wtm_pos + 240 > len {
            return None;
        }
        let c_start = wtm_pos - 39;
        let c_end = wtm_pos + 240;
        core_starts.push(c_start);
        core_ends.push(c_end);
    }

    Some(vec![
        MultiZoneSliceSpec {
            ranges: (0..n).map(|i| (0, core_starts[i])).collect(),
            free_leading_gaps: true,
            free_trailing_gaps: false,
            is_conserved_core_zone: false,
            terminal_overhang_cap: None,
        },
        MultiZoneSliceSpec {
            ranges: (0..n).map(|i| (core_starts[i], core_ends[i])).collect(),
            free_leading_gaps: false,
            free_trailing_gaps: false,
            is_conserved_core_zone: true,
            terminal_overhang_cap: None,
        },
        MultiZoneSliceSpec {
            ranges: (0..n).map(|i| (core_ends[i], seqs[i].len())).collect(),
            free_leading_gaps: false,
            free_trailing_gaps: true,
            is_conserved_core_zone: false,
            terminal_overhang_cap: None,
        },
    ])
}

// ============================================================================
// STAGE 2: COOLEY-TUKEY RADIX-2 FAST FOURIER TRANSFORM (`FFT-NS`)
// ============================================================================

/// In-place Radix-2 Cooley-Tukey Complex Fast Fourier Transform (`invert = false` for FFT,
/// `invert = true` for IFFT). `re` and `im` lengths must be a power of two.
pub fn fft_radix2_in_place(re: &mut [f64], im: &mut [f64], invert: bool) {
    let n = re.len();
    debug_assert!(n.is_power_of_two() && re.len() == im.len());
    if n <= 1 {
        return;
    }

    // Bit-reversal permutation
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while (j & bit) != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            re.swap(i, j);
            im.swap(i, j);
        }
    }

    // Butterfly stages
    let mut len = 2usize;
    while len <= n {
        let half = len >> 1;
        let ang = 2.0 * PI / (len as f64) * if invert { 1.0 } else { -1.0 };
        let wlen_re = ang.cos();
        let wlen_im = ang.sin();

        for i in (0..n).step_by(len) {
            let mut w_re = 1.0f64;
            let mut w_im = 0.0f64;
            for k in 0..half {
                let u_re = re[i + k];
                let u_im = im[i + k];
                let v_re = re[i + k + half] * w_re - im[i + k + half] * w_im;
                let v_im = re[i + k + half] * w_im + im[i + k + half] * w_re;
                re[i + k] = u_re + v_re;
                im[i + k] = u_im + v_im;
                re[i + k + half] = u_re - v_re;
                im[i + k + half] = u_im - v_im;
                let next_w_re = w_re * wlen_re - w_im * wlen_im;
                w_im = w_re * wlen_im + w_im * wlen_re;
                w_re = next_w_re;
            }
        }
        len <<= 1;
    }

    if invert {
        let inv_n = 1.0 / (n as f64);
        for x in re.iter_mut() {
            *x *= inv_n;
        }
        for y in im.iter_mut() {
            *y *= inv_n;
        }
    }
}

/// Discover homologous diagonal segments between `seq_a` and `seq_b` using 4-channel
/// Cooley-Tukey FFT cross-correlation (`FFT-NS`, Katoh et al. 2002).
/// Returns colinear homologous anchor points `(pos_a, pos_b, window_len, weight)`.
pub fn find_fft_homologous_anchors(
    seq_a: &[u8],
    seq_b: &[u8],
) -> Vec<(usize, usize, usize, i32)> {
    let la = seq_a.len();
    let lb = seq_b.len();
    if la < 20 || lb < 20 {
        return Vec::new();
    }

    let n_fft = (la + lb).next_power_of_two();
    let mut sum_corr_re = vec![0.0f64; n_fft];
    let mut sum_corr_im = vec![0.0f64; n_fft];

    let bases = [b'A', b'C', b'G', b'T'];
    let mut fa_re = vec![0.0f64; n_fft];
    let mut fa_im = vec![0.0f64; n_fft];
    let mut fb_re = vec![0.0f64; n_fft];
    let mut fb_im = vec![0.0f64; n_fft];

    for &base in &bases {
        fa_re.fill(0.0);
        fa_im.fill(0.0);
        fb_re.fill(0.0);
        fb_im.fill(0.0);

        for (i, &b) in seq_a.iter().enumerate() {
            fa_re[i] = if b == base { 0.75 } else { -0.25 };
        }
        for (j, &b) in seq_b.iter().enumerate() {
            fb_re[j] = if b == base { 0.75 } else { -0.25 };
        }

        fft_radix2_in_place(&mut fa_re, &mut fa_im, false);
        fft_radix2_in_place(&mut fb_re, &mut fb_im, false);

        // Multiply FA(k) * conj(FB(k))
        for k in 0..n_fft {
            sum_corr_re[k] += fa_re[k] * fb_re[k] + fa_im[k] * fb_im[k];
            sum_corr_im[k] += fa_im[k] * fb_re[k] - fa_re[k] * fb_im[k];
        }
    }

    fft_radix2_in_place(&mut sum_corr_re, &mut sum_corr_im, true);

    // Collect lag scores for `lag = i_a as isize - j_b as isize`
    let mut lags: Vec<(isize, f64)> = Vec::with_capacity(la + lb);
    let mut mean = 0.0f64;
    let mut m2 = 0.0f64;
    let mut count = 0.0f64;

    for lag in -(lb as isize - 12)..=(la as isize - 12) {
        let idx = if lag >= 0 {
            lag as usize
        } else {
            (n_fft as isize + lag) as usize
        };
        let val = sum_corr_re[idx];
        count += 1.0;
        let d = val - mean;
        mean += d / count;
        m2 += d * (val - mean);
        lags.push((lag, val));
    }
    let sd = (m2 / count.max(1.0)).sqrt().max(1e-6);

    // Select top peak lags with Z >= 2.2
    lags.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    let mut candidate_segments: Vec<(usize, usize, usize, i32)> = Vec::new();
    let win = 16usize;

    for &(lag, val) in lags.iter().take(6) {
        if (val - mean) / sd < 2.2 {
            break;
        }
        let start_a = if lag > 0 { lag as usize } else { 0 };
        let start_b = if lag < 0 { (-lag) as usize } else { 0 };
        let diag_len = (la.saturating_sub(start_a)).min(lb.saturating_sub(start_b));
        if diag_len < win {
            continue;
        }

        let mut k = 0usize;
        while k + win <= diag_len {
            let mut sc = 0i32;
            let mut matches = 0usize;
            for off in 0..win {
                let ba = seq_a[start_a + k + off];
                let bb = seq_b[start_b + k + off];
                sc += dna200_pair_score(ba, bb);
                if ba == bb {
                    matches += 1;
                }
            }
            if matches >= 11 && sc >= 700 {
                let w = ((sc as usize) / win).clamp(80, 420) as i32;
                candidate_segments.push((start_a + k, start_b + k, win, w));
                k += win / 2;
            } else {
                k += 4;
            }
        }
    }

    if candidate_segments.is_empty() {
        return Vec::new();
    }

    // Colinear DP chaining of FFT segments
    candidate_segments.sort_by_key(|&(pa, pb, _, _)| (pa, pb));
    let m = candidate_segments.len();
    let mut dp: Vec<i64> = candidate_segments
        .iter()
        .map(|&(_, _, len, w)| (len as i64) * (w as i64))
        .collect();
    let mut prev = vec![None::<usize>; m];
    let mut best_end = 0usize;
    let mut best_sc = dp[0];

    for i in 1..m {
        let (ai, bi, len_i, wi) = candidate_segments[i];
        let base_i = (len_i as i64) * (wi as i64);
        for u in 0..i {
            let (au, bu, len_u, _) = candidate_segments[u];
            if ai >= au + len_u / 2 && bi >= bu + len_u / 2 && ai > au && bi > bu {
                let drift = ((ai - au) as i64 - (bi - bu) as i64).abs();
                let cand = dp[u] + base_i - drift * 12;
                if cand > dp[i] {
                    dp[i] = cand;
                    prev[i] = Some(u);
                }
            }
        }
        if dp[i] > best_sc {
            best_sc = dp[i];
            best_end = i;
        }
    }

    let mut chain = Vec::new();
    let mut cur = Some(best_end);
    while let Some(idx) = cur {
        chain.push(candidate_segments[idx]);
        cur = prev[idx];
    }
    chain.reverse();
    chain
}

// ============================================================================
// STAGE 3: `L-INS-i` / `E-INS-i` LOCAL SMITH-WATERMAN & TRANSITIVE CONSISTENCY (`hat3`)
// ============================================================================

#[inline]
pub fn is_transition(a: u8, b: u8) -> bool {
    matches!(
        (a, b),
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')
    )
}

/// MAFFT `DNA200` (`Kimura 200`, Ti/Tv = 2.0) pairwise residue score scaled by 100.
#[inline]
pub fn dna200_pair_score(a: u8, b: u8) -> i32 {
    if a == b && (a == b'A' || a == b'C' || a == b'G' || a == b'T') {
        100
    } else if is_transition(a, b) {
        20
    } else if a == b'N' || b == b'N' {
        0
    } else {
        -75
    }
}

/// Compute local Smith-Waterman + multi-segment homologous tracts between `seq_a` and `seq_b`
/// (`tbfast-pair` `alg=L` / `alg=E` in MAFFT `L-INS-i` / `E-INS-i`).
/// Returns `(normalized_local_similarity, matched_residue_pairs_with_weights, fft_anchor_count)`.
fn compute_local_pair_hat3_matches(
    seq_a: &[u8],
    seq_b: &[u8],
    use_rna_stem_bonus: bool,
    stem_prop_a: &[i32],
    stem_prop_b: &[i32],
) -> (f64, Vec<(usize, usize, i32)>, usize) {
    let la = seq_a.len();
    let lb = seq_b.len();
    if la == 0 || lb == 0 {
        return (0.0, Vec::new(), 0);
    }

    // 1. FFT homologous anchors (`FFT-NS`)
    let fft_anchors = find_fft_homologous_anchors(seq_a, seq_b);
    let fft_count = fft_anchors.len();

    // 2. If sequences are very long (> 1,400 bp), use anchor-seeded banded local extraction;
    //    for ITS / rRNA amplicons (<= 1,400 bp), run exact Smith-Waterman + semi-global local tract extraction.
    let mut matches: Vec<(usize, usize, i32)> = Vec::new();
    for &(pa, pb, wlen, w) in &fft_anchors {
        for off in 0..wlen {
            if pa + off < la && pb + off < lb && seq_a[pa + off] == seq_b[pb + off] {
                matches.push((pa + off, pb + off, w));
            }
        }
    }

    if la <= 1_800 && lb <= 1_800 {
        let sw_tracts = extract_smith_waterman_local_tracts(
            seq_a,
            seq_b,
            use_rna_stem_bonus,
            stem_prop_a,
            stem_prop_b,
        );
        for (pa, pb, w) in sw_tracts {
            matches.push((pa, pb, w));
        }
    }

    // Compute normalized local similarity in [0, 1]
    let min_l = la.min(lb).max(1) as f64;
    let unique_covered = {
        let mut cov_a = vec![false; la];
        let mut id_cnt = 0usize;
        for &(pa, pb, _) in &matches {
            if pa < la && pb < lb && !cov_a[pa] && seq_a[pa] == seq_b[pb] {
                cov_a[pa] = true;
                id_cnt += 1;
            }
        }
        id_cnt as f64
    };
    let sim = (unique_covered / min_l).clamp(0.0, 1.0);

    (sim, matches, fft_count)
}

/// Extract high-confidence local homologous alignment tracts via Smith-Waterman affine-gap DP
/// (`alg=L` in MAFFT `tbfast-pair`: `op = -200`, `ep = -10`), filtering out short random coincidences
/// (`tract_len >= 8` and `identity >= 62%`) so non-homologous spacer loops receive zero `hat3` weight.
fn extract_smith_waterman_local_tracts(
    a: &[u8],
    b: &[u8],
    use_rna_stem_bonus: bool,
    stem_a: &[i32],
    stem_b: &[i32],
) -> Vec<(usize, usize, i32)> {
    let n = a.len();
    let m = b.len();
    if n < 8 || m < 8 {
        return Vec::new();
    }

    const GAP_OPEN: i32 = -220;
    const GAP_EXT: i32 = -12;

    let cols = m + 1;
    let total = (n + 1) * cols;
    let mut h = vec![0i32; total];
    let mut e = vec![-1_000_000i32; total];
    let mut f = vec![-1_000_000i32; total];
    let mut trace = vec![0u8; total]; // 0 = stop, 1 = diag, 2 = up (delete), 3 = left (insert)

    let mut peaks: Vec<(i32, usize, usize)> = Vec::new();

    for i in 1..=n {
        let row = i * cols;
        let prev = (i - 1) * cols;
        let ai = a[i - 1];
        let sa = if use_rna_stem_bonus {
            stem_a.get(i - 1).copied().unwrap_or(0)
        } else {
            0
        };

        for j in 1..=m {
            let idx = row + j;
            let bj = b[j - 1];

            let e_open = h[prev + j] + GAP_OPEN;
            let e_ext = e[prev + j] + GAP_EXT;
            let e_val = e_open.max(e_ext);
            e[idx] = e_val;

            let f_open = h[row + j - 1] + GAP_OPEN;
            let f_ext = f[row + j - 1] + GAP_EXT;
            let f_val = f_open.max(f_ext);
            f[idx] = f_val;

            let mut sub = dna200_pair_score(ai, bj);
            if use_rna_stem_bonus {
                let sb = stem_b.get(j - 1).copied().unwrap_or(0);
                if sa > 0 && sb > 0 && (ai == bj || is_transition(ai, bj)) {
                    sub += (sa.min(sb)) / 2;
                }
            }
            let d_val = h[prev + j - 1] + sub;

            let mut best = 0i32;
            let mut tr = 0u8;
            if d_val > best {
                best = d_val;
                tr = 1;
            }
            if e_val > best {
                best = e_val;
                tr = 2;
            }
            if f_val > best {
                best = f_val;
                tr = 3;
            }
            h[idx] = best;
            trace[idx] = tr;

            if best >= 680 && tr == 1 {
                peaks.push((best, i, j));
            }
        }
    }

    if peaks.is_empty() {
        return Vec::new();
    }

    // Sort local alignment peaks descending by score and extract non-overlapping contiguous runs
    peaks.sort_unstable_by(|x, y| y.0.cmp(&x.0));

    let mut used_a = vec![false; n];
    let mut used_b = vec![false; m];
    let mut tracts: Vec<(usize, usize, i32)> = Vec::new();

    for &(peak_sc, end_i, end_j) in peaks.iter().take(64) {
        if used_a[end_i - 1] || used_b[end_j - 1] {
            continue;
        }
        let mut i = end_i;
        let mut j = end_j;
        let mut path: Vec<(usize, usize)> = Vec::new();

        while i > 0 && j > 0 {
            let idx = i * cols + j;
            let tr = trace[idx];
            if tr == 0 || h[idx] <= 0 {
                break;
            }
            if used_a[i - 1] || used_b[j - 1] {
                break;
            }
            match tr {
                1 => {
                    path.push((i - 1, j - 1));
                    i -= 1;
                    j -= 1;
                }
                2 => {
                    i -= 1;
                }
                _ => {
                    j -= 1;
                }
            }
        }

        if path.len() < 8 {
            continue;
        }
        path.reverse();

        // Filter contiguous gap-free sub-runs of length >= 6 within the local path
        // so only genuine homologous blocks get stored in `hat3`
        let mut run_start = 0usize;
        while run_start < path.len() {
            let mut run_end = run_start + 1;
            while run_end < path.len()
                && path[run_end].0 == path[run_end - 1].0 + 1
                && path[run_end].1 == path[run_end - 1].1 + 1
            {
                run_end += 1;
            }
            let run_len = run_end - run_start;
            if run_len >= 6 {
                let mut exact = 0usize;
                let mut ti = 0usize;
                for &(pa, pb) in &path[run_start..run_end] {
                    if a[pa] == b[pb] {
                        exact += 1;
                    } else if is_transition(a[pa], b[pb]) {
                        ti += 1;
                    }
                }
                // Require strong local identity (`(exact * 2 + ti) >= run_len * 13 / 10`, i.e. >= 65% weighted id)
                if exact * 2 + ti >= (run_len * 13) / 10 {
                    let base_w = ((peak_sc / 8).clamp(110, 480)
                        + ((exact * 180) / run_len.max(1)) as i32)
                        .min(600);
                    for &(pa, pb) in &path[run_start..run_end] {
                        used_a[pa] = true;
                        used_b[pb] = true;
                        tracts.push((pa, pb, base_w));
                    }
                }
            }
            run_start = run_end;
        }
    }

    tracts
}

/// Compute local RNA stem-loop base-pairing propensity (`Q-INS-i` lightweight local hairpin scanner)
/// for each nucleotide position in `seq` (detecting inverted complementary stems of 4–8 bp
/// separated by 3–18 nt loops, characteristic of `ITS2` Helices I–IV and `16S` V-loops).
fn compute_local_stem_propensity(seq: &[u8]) -> Vec<i32> {
    let n = seq.len();
    let mut prop = vec![0i32; n];
    if n < 14 {
        return prop;
    }

    for i in 0..n.saturating_sub(11) {
        for loop_len in 4..=16 {
            let j_start = i + 4 + loop_len;
            if j_start + 4 > n {
                break;
            }
            let mut pairs = 0i32;
            let mut gc_pairs = 0i32;
            for k in 0..5 {
                if i + k < j_start && j_start + 4 >= k && j_start + 4 - k < n {
                    let b1 = seq[i + k];
                    let b2 = seq[j_start + 4 - k];
                    match (b1, b2) {
                        (b'G', b'C') | (b'C', b'G') => {
                            pairs += 1;
                            gc_pairs += 1;
                        }
                        (b'A', b'T') | (b'T', b'A') | (b'G', b'T') | (b'T', b'G') => {
                            pairs += 1;
                        }
                        _ => {}
                    }
                }
            }
            if pairs >= 4 && gc_pairs >= 2 {
                let bonus = pairs * 6 + gc_pairs * 4;
                for k in 0..5 {
                    prop[i + k] = prop[i + k].max(bonus);
                    if j_start + 4 >= k && j_start + 4 - k < n {
                        prop[j_start + 4 - k] = prop[j_start + 4 - k].max(bonus);
                    }
                }
            }
        }
    }
    prop
}

/// Build the all-pairs `hat3` local consistency matrix with COFFEE-style transitive consistency
/// relaxation (`W_ab(i, j) <- H_ab(i, j) + sum_c w_c H_ac(i, k) H_cb(k, j)`), plus the initial
/// pairwise evolutionary distance matrix $D_{ab}$.
fn build_hat3_and_distance_matrix(
    seqs: &[Vec<u8>],
    use_transitive: bool,
    use_rna_stem: bool,
) -> (SparseHat3ConsistencyMatrix, Vec<Vec<f64>>, usize) {
    let n = seqs.len();
    let lens: Vec<usize> = seqs.iter().map(|s| s.len()).collect();
    let mut hat3 = SparseHat3ConsistencyMatrix::new(&lens);
    let mut dist_mat = vec![vec![0.0f64; n]; n];

    let stem_props: Vec<Vec<i32>> = if use_rna_stem {
        seqs.par_iter()
            .map(|s| compute_local_stem_propensity(s))
            .collect()
    } else {
        vec![Vec::new(); n]
    };

    let kmer_sets: Vec<Vec<u16>> = seqs.iter().map(|s| extract_6mer_multiset(s)).collect();

    // Parallelize all-pairs `(a, b)` local Smith-Waterman + FFT anchor computation
    let pair_coords: Vec<(usize, usize)> = (0..n)
        .flat_map(|a| ((a + 1)..n).map(move |b| (a, b)))
        .collect();

    let pair_results: Vec<(usize, usize, f64, Vec<(usize, usize, i32)>, usize)> = pair_coords
        .par_iter()
        .map(|&(a, b)| {
            let (local_sim, tracts, fft_cnt) = compute_local_pair_hat3_matches(
                &seqs[a],
                &seqs[b],
                use_rna_stem,
                &stem_props[a],
                &stem_props[b],
            );
            let kmer_sim = shared_6mer_similarity(&kmer_sets[a], &kmer_sets[b]);
            let combined_sim = (0.65 * local_sim + 0.35 * kmer_sim).clamp(0.0, 1.0);
            let d = (1.0 - combined_sim).clamp(0.0001, 0.9999);
            (a, b, d, tracts, fft_cnt)
        })
        .collect();

    let mut total_fft_anchors = 0usize;
    for (a, b, d, tracts, fft_cnt) in pair_results {
        dist_mat[a][b] = d;
        dist_mat[b][a] = d;
        total_fft_anchors += fft_cnt;
        for (pa, pb, w) in tracts {
            hat3.add_weight(a, b, pa, pb, w);
        }
    }

    // COFFEE / MAFFT L-INS-i Transitive Consistency Relaxation across intermediate sequences `c`:
    // If `a[i]` aligns locally to `c[k]` AND `c[k]` aligns locally to `b[j]`, reinforce `(a[i], b[j])`.
    if use_transitive && n >= 3 && n <= 64 {
        let base_hat3 = hat3.clone();
        for a in 0..n {
            for b in (a + 1)..n {
                for c in 0..n {
                    if c == a || c == b {
                        continue;
                    }
                    let w_ac_cb =
                        ((1.0 - dist_mat[a][c] * 0.5) * (1.0 - dist_mat[c][b] * 0.5)).clamp(0.2, 1.0);
                    for pos_a in 0..lens[a] {
                        let links_ac = get_links_between(&base_hat3, a, c, pos_a);
                        if links_ac.is_empty() {
                            continue;
                        }
                        for (pos_c, w_ac) in links_ac {
                            let links_cb = get_links_between(&base_hat3, c, b, pos_c);
                            for (pos_b, w_cb) in links_cb {
                                let trans_w = (((w_ac.min(w_cb) as f64) * 0.35 * w_ac_cb)
                                    / ((n - 2).max(1) as f64))
                                    as i32;
                                if trans_w >= 12 {
                                    hat3.add_weight(a, b, pos_a, pos_b, trans_w);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    (hat3, dist_mat, total_fft_anchors)
}

fn get_links_between(
    hat3: &SparseHat3ConsistencyMatrix,
    src: usize,
    dst: usize,
    pos_src: usize,
) -> Vec<(usize, i32)> {
    let n = hat3.num_seqs;
    if src < dst {
        let pidx = SparseHat3ConsistencyMatrix::pair_idx(n, src, dst);
        if let Some(row) = hat3.forward_links.get(pidx).and_then(|m| m.get(pos_src)) {
            return row.iter().map(|&(j, w)| (j as usize, w)).collect();
        }
        Vec::new()
    } else {
        // Scan inverse links for `(dst, src)` where `j == pos_src`
        let pidx = SparseHat3ConsistencyMatrix::pair_idx(n, dst, src);
        let mut out = Vec::new();
        if let Some(mat) = hat3.forward_links.get(pidx) {
            let target_u16 = pos_src.min(u16::MAX as usize) as u16;
            for (pos_dst, row) in mat.iter().enumerate() {
                for &(j, w) in row {
                    if j == target_u16 {
                        out.push((pos_dst, w));
                    }
                }
            }
        }
        out
    }
}

// ============================================================================
// STAGE 4: TWO-CYCLE UPGMA / MODIFIED-LINKAGE GUIDE TREE & SEQUENCE WEIGHTS
// ============================================================================

/// Construct a rooted binary UPGMA / MAFFT modified-linkage (`sueff = 0.10`) guide tree
/// and compute Clustal/MAFFT branch-proportional sequence weights `w[0..n]`.
pub fn build_mafft_guide_tree_and_weights(
    dist_matrix: &[Vec<f64>],
) -> (Vec<GuideTreeMergeStep>, Vec<f64>) {
    let n = dist_matrix.len();
    if n <= 1 {
        return (Vec::new(), vec![1.0; n]);
    }

    let mut active_clusters: Vec<(Vec<usize>, f64)> =
        (0..n).map(|i| (vec![i], 0.0f64)).collect();
    let mut cluster_dists = dist_matrix.to_vec();
    let mut steps = Vec::with_capacity(n - 1);
    let mut seq_weights = vec![0.0f64; n];

    // MAFFT `sueff_global = 0.10`: blends 90% UPGMA average linkage + 10% minimum linkage
    const SUEFF: f64 = 0.10;

    while active_clusters.len() > 1 {
        let m = active_clusters.len();
        let mut best_u = 0usize;
        let mut best_v = 1usize;
        let mut best_d = f64::MAX;

        for u in 0..m {
            for v in (u + 1)..m {
                if cluster_dists[u][v] < best_d {
                    best_d = cluster_dists[u][v];
                    best_u = u;
                    best_v = v;
                }
            }
        }

        let new_height = (best_d * 0.5).max(
            active_clusters[best_u]
                .1
                .max(active_clusters[best_v].1)
                + 0.0005,
        );
        let branch_u = (new_height - active_clusters[best_u].1).max(0.001);
        let branch_v = (new_height - active_clusters[best_v].1).max(0.001);

        let size_u = active_clusters[best_u].0.len() as f64;
        let size_v = active_clusters[best_v].0.len() as f64;
        for &idx in &active_clusters[best_u].0 {
            seq_weights[idx] += branch_u / size_u;
        }
        for &idx in &active_clusters[best_v].0 {
            seq_weights[idx] += branch_v / size_v;
        }

        let left_members = active_clusters[best_u].0.clone();
        let right_members = active_clusters[best_v].0.clone();
        steps.push(GuideTreeMergeStep {
            left_members: left_members.clone(),
            right_members: right_members.clone(),
            branch_height: new_height,
        });

        let mut merged_members = left_members;
        merged_members.extend(right_members);

        // Compute distances from merged cluster to all remaining clusters
        let mut new_row = Vec::with_capacity(m - 1);
        for k in 0..m {
            if k == best_u || k == best_v {
                continue;
            }
            let du = cluster_dists[best_u][k];
            let dv = cluster_dists[best_v][k];
            let avg_d = (du * size_u + dv * size_v) / (size_u + size_v);
            let min_d = du.min(dv);
            let blended = (1.0 - SUEFF) * avg_d + SUEFF * min_d;
            new_row.push(blended);
        }

        // Rebuild `active_clusters` and `cluster_dists`
        let mut next_clusters = Vec::with_capacity(m - 1);
        let mut kept_indices = Vec::with_capacity(m - 2);
        for k in 0..m {
            if k != best_u && k != best_v {
                next_clusters.push(active_clusters[k].clone());
                kept_indices.push(k);
            }
        }
        next_clusters.push((merged_members, new_height));

        let new_m = next_clusters.len();
        let mut next_dists = vec![vec![0.0f64; new_m]; new_m];
        for (i_new, &i_old) in kept_indices.iter().enumerate() {
            for (j_new, &j_old) in kept_indices.iter().enumerate() {
                next_dists[i_new][j_new] = cluster_dists[i_old][j_old];
            }
            next_dists[i_new][new_m - 1] = new_row[i_new];
            next_dists[new_m - 1][i_new] = new_row[i_new];
        }

        active_clusters = next_clusters;
        cluster_dists = next_dists;
    }

    // Normalize sequence weights so mean weight = 1.0
    let sum_w: f64 = seq_weights.iter().sum();
    if sum_w > 1e-9 {
        let scale = (n as f64) / sum_w;
        for w in seq_weights.iter_mut() {
            *w = (*w * scale).clamp(0.25, 2.50);
        }
    } else {
        seq_weights.fill(1.0);
    }

    (steps, seq_weights)
}

/// Recompute pairwise evolutionary distances from an aligned sequence matrix (`FFT-NS-2` tree refinement).
fn compute_aligned_distance_matrix(aligned_matrix: &[Vec<u8>]) -> Vec<Vec<f64>> {
    let n = aligned_matrix.len();
    let mut dist = vec![vec![0.0f64; n]; n];
    for a in 0..n {
        for b in (a + 1)..n {
            let sa = &aligned_matrix[a];
            let sb = &aligned_matrix[b];
            let len = sa.len().min(sb.len());
            let mut compared = 0usize;
            let mut diff_ti = 0usize;
            let mut diff_tv = 0usize;
            let mut internal_gaps = 0usize;
            for c in 0..len {
                let ba = sa[c];
                let bb = sb[c];
                if ba != b'-' && bb != b'-' {
                    compared += 1;
                    if ba != bb {
                        if is_transition(ba, bb) {
                            diff_ti += 1;
                        } else {
                            diff_tv += 1;
                        }
                    }
                } else if ba != b'-' || bb != b'-' {
                    internal_gaps += 1;
                }
            }
            let d = if compared >= 12 {
                let p_sub = (diff_ti as f64 + 1.5 * diff_tv as f64) / (compared as f64);
                let p_gap = (internal_gaps as f64) / ((compared + internal_gaps).max(1) as f64);
                (0.82 * p_sub + 0.18 * p_gap).clamp(0.0001, 0.9999)
            } else {
                0.85
            };
            dist[a][b] = d;
            dist[b][a] = d;
        }
    }
    dist
}

// ============================================================================
// STAGE 5: GROUP-TO-GROUP PROFILE GOTOH AFFINE-GAP DYNAMIC PROGRAMMING
// ============================================================================

/// Representation of a sub-alignment Profile $G$ over `member_indices` of the original sequences.
#[derive(Debug, Clone)]
struct AlignmentProfile {
    /// Indices of original sequences (`0..N`) belonging to this profile.
    members: Vec<usize>,
    /// Aligned rows (`members.len() x num_cols`).
    rows: Vec<Vec<u8>>,
    /// For each member `m_idx` and column `c`, the 0-based ungapped residue coordinate `Some(pos)`
    /// in original sequence `members[m_idx]`, or `None` if `rows[m_idx][c] == b'-'`.
    orig_pos: Vec<Vec<Option<u16>>>,
}

impl AlignmentProfile {
    fn from_single_sequence(seq_idx: usize, seq: &[u8]) -> Self {
        let orig: Vec<Option<u16>> = (0..seq.len())
            .map(|p| Some(p.min(u16::MAX as usize) as u16))
            .collect();
        Self {
            members: vec![seq_idx],
            rows: vec![seq.to_vec()],
            orig_pos: vec![orig],
        }
    }

    #[inline]
    fn num_cols(&self) -> usize {
        self.rows.first().map(|r| r.len()).unwrap_or(0)
    }
}

/// Align two multi-sequence profiles `prof_a` and `prof_b` via exact Gotoh 3-state (`M`, `I`, `D`)
/// affine-gap dynamic programming combining:
/// 1. Weighted Sum-of-Pairs (`WSP`) `DNA200` score,
/// 2. `L-INS-i` / `E-INS-i` `hat3` local residue consistency score $W_{ab}(i, j)$,
/// 3. MAFFT v7.110+ position-specific gap-open and gap-extension penalties (`op = 1.53`, `ep -> 0.0`),
/// 4. Semi-global free 5' leading and/or 3' trailing overhang end-gaps.
#[allow(clippy::too_many_arguments)]
fn align_profiles_gotoh_mafft(
    prof_a: &AlignmentProfile,
    prof_b: &AlignmentProfile,
    hat3: &SparseHat3ConsistencyMatrix,
    seq_weights: &[f64],
    free_leading_gaps: bool,
    free_trailing_gaps: bool,
    is_conserved_core_zone: bool,
) -> AlignmentProfile {
    let la = prof_a.num_cols();
    let lb = prof_b.num_cols();
    let na = prof_a.members.len();
    let nb = prof_b.members.len();

    if la == 0 && lb == 0 {
        let mut members = prof_a.members.clone();
        members.extend_from_slice(&prof_b.members);
        return AlignmentProfile {
            members,
            rows: vec![Vec::new(); na + nb],
            orig_pos: vec![Vec::new(); na + nb],
        };
    }
    if la == 0 {
        let mut members = prof_a.members.clone();
        members.extend_from_slice(&prof_b.members);
        let mut rows = vec![vec![b'-'; lb]; na];
        rows.extend(prof_b.rows.clone());
        let mut orig_pos = vec![vec![None; lb]; na];
        orig_pos.extend(prof_b.orig_pos.clone());
        return AlignmentProfile {
            members,
            rows,
            orig_pos,
        };
    }
    if lb == 0 {
        let mut members = prof_a.members.clone();
        members.extend_from_slice(&prof_b.members);
        let mut rows = prof_a.rows.clone();
        rows.extend(vec![vec![b'-'; la]; nb]);
        let mut orig_pos = prof_a.orig_pos.clone();
        orig_pos.extend(vec![vec![None; la]; nb]);
        return AlignmentProfile {
            members,
            rows,
            orig_pos,
        };
    }

    // Precompute normalized member weights within Profile A and Profile B
    let wa: Vec<f64> = prof_a.members.iter().map(|&m| seq_weights[m]).collect();
    let wb: Vec<f64> = prof_b.members.iter().map(|&m| seq_weights[m]).collect();
    let sum_wa: f64 = wa.iter().sum::<f64>().max(1e-9);
    let sum_wb: f64 = wb.iter().sum::<f64>().max(1e-9);
    let norm_wa: Vec<f64> = wa.iter().map(|&w| w / sum_wa).collect();
    let norm_wb: Vec<f64> = wb.iter().map(|&w| w / sum_wb).collect();

    // Precompute per-column nucleotide frequencies `[A, C, G, T]` and gap fractions `gap_a[c]`, `gap_b[d]`
    let mut freq_a = vec![[0.0f64; 4]; la];
    let mut gap_a = vec![0.0f64; la];
    for c in 0..la {
        let mut g = 0.0f64;
        for (mi, row) in prof_a.rows.iter().enumerate() {
            let w = norm_wa[mi];
            match row[c] {
                b'A' => freq_a[c][0] += w,
                b'C' => freq_a[c][1] += w,
                b'G' => freq_a[c][2] += w,
                b'T' => freq_a[c][3] += w,
                b'-' => g += w,
                _ => {}
            }
        }
        gap_a[c] = g.clamp(0.0, 1.0);
    }

    let mut freq_b = vec![[0.0f64; 4]; lb];
    let mut gap_b = vec![0.0f64; lb];
    for d in 0..lb {
        let mut g = 0.0f64;
        for (mj, row) in prof_b.rows.iter().enumerate() {
            let w = norm_wb[mj];
            match row[d] {
                b'A' => freq_b[d][0] += w,
                b'C' => freq_b[d][1] += w,
                b'G' => freq_b[d][2] += w,
                b'T' => freq_b[d][3] += w,
                b'-' => g += w,
                _ => {}
            }
        }
        gap_b[d] = g.clamp(0.0, 1.0);
    }

    // Position-specific gap-opening penalties (`op = 1.53` in MAFFT `tbfast.exe`, scaled to integer DP units):
    // In conserved rDNA core subunits (`5.8S`), gap open is boosted (`-640`) so internal gaps are prevented;
    // in spacer zones (`ITS1`, `ITS2`, `16S V-loops`), gap open is `-410` modulated by `(1 - 0.55 * gap_frac)`
    // so gaps cleanly coalesce at existing gap boundaries, with near-zero gap extension (`-6` per column,
    // matching MAFFT's `-0.00` to `-0.02` extension penalty).
    let base_gap_open: f64 = if is_conserved_core_zone {
        -680.0
    } else {
        -460.0
    };
    let gap_ext: i32 = if is_conserved_core_zone { -28 } else { -12 };

    let open_pen_a: Vec<i32> = (0..la)
        .map(|c| {
            let factor = (1.0 - 0.55 * gap_a[c]).clamp(0.38, 1.0);
            (base_gap_open * factor).round() as i32
        })
        .collect();
    let open_pen_b: Vec<i32> = (0..lb)
        .map(|d| {
            let factor = (1.0 - 0.55 * gap_b[d]).clamp(0.38, 1.0);
            (base_gap_open * factor).round() as i32
        })
        .collect();

    // Precompute sparse consistency bonus matrix `cons_bonus[c * lb + d]` from `hat3`
    let mut cons_bonus = vec![0i32; la * lb];
    for (mi, &seq_a_id) in prof_a.members.iter().enumerate() {
        let wa_m = norm_wa[mi];
        for (mj, &seq_b_id) in prof_b.members.iter().enumerate() {
            let pair_w = wa_m * norm_wb[mj];
            let (sa, sb, flip) = if seq_a_id < seq_b_id {
                (seq_a_id, seq_b_id, false)
            } else {
                (seq_b_id, seq_a_id, true)
            };
            let pidx = SparseHat3ConsistencyMatrix::pair_idx(hat3.num_seqs, sa, sb);
            let Some(pair_table) = hat3.forward_links.get(pidx) else {
                continue;
            };

            if !flip {
                // Build quick lookup from `pos_b -> col_d` in `prof_b`
                let max_pos_b = prof_b.orig_pos[mj]
                    .iter()
                    .filter_map(|&x| x)
                    .max()
                    .unwrap_or(0) as usize;
                let mut pos_b_to_col = vec![u16::MAX; max_pos_b + 1];
                for (d, &op) in prof_b.orig_pos[mj].iter().enumerate() {
                    if let Some(pb) = op {
                        pos_b_to_col[pb as usize] = d as u16;
                    }
                }
                for (c, &op_a) in prof_a.orig_pos[mi].iter().enumerate() {
                    if let Some(pa) = op_a {
                        if let Some(links) = pair_table.get(pa as usize) {
                            for &(pb_u16, w) in links {
                                if let Some(&d_u16) = pos_b_to_col.get(pb_u16 as usize) {
                                    if d_u16 != u16::MAX {
                                        let d = d_u16 as usize;
                                        cons_bonus[c * lb + d] +=
                                            ((w as f64) * pair_w * 0.85).round() as i32;
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                // `sa == seq_b_id < sb == seq_a_id`
                let max_pos_a = prof_a.orig_pos[mi]
                    .iter()
                    .filter_map(|&x| x)
                    .max()
                    .unwrap_or(0) as usize;
                let mut pos_a_to_col = vec![u16::MAX; max_pos_a + 1];
                for (c, &op) in prof_a.orig_pos[mi].iter().enumerate() {
                    if let Some(pa) = op {
                        pos_a_to_col[pa as usize] = c as u16;
                    }
                }
                for (d, &op_b) in prof_b.orig_pos[mj].iter().enumerate() {
                    if let Some(pb) = op_b {
                        if let Some(links) = pair_table.get(pb as usize) {
                            for &(pa_u16, w) in links {
                                if let Some(&c_u16) = pos_a_to_col.get(pa_u16 as usize) {
                                    if c_u16 != u16::MAX {
                                        let c = c_u16 as usize;
                                        cons_bonus[c * lb + d] +=
                                            ((w as f64) * pair_w * 0.85).round() as i32;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Gotoh 3-Matrix (`M`, `D`, `I`) Dynamic Programming over `(la + 1) x (lb + 1)`
    const NEG_INF: i32 = -1_000_000_000;
    const TERM_GAP_EXT: i32 = -1; // Near-zero terminal overhang tie-breaker (`-0.01` per column)

    let cols = lb + 1;
    let total = (la + 1) * cols;
    let mut mat_m = vec![NEG_INF; total];
    let mut mat_d = vec![NEG_INF; total]; // gap in Profile B (consume col of A)
    let mut mat_i = vec![NEG_INF; total]; // gap in Profile A (consume col of B)
    let mut tr_m = vec![0u8; total]; // 0 = from M, 1 = from D, 2 = from I
    let mut tr_d = vec![0u8; total]; // 0 = open from M, 1 = extend from D
    let mut tr_i = vec![0u8; total]; // 0 = open from M, 1 = extend from I

    mat_m[0] = 0;
    for i in 1..=la {
        let idx = i * cols;
        if free_leading_gaps {
            mat_d[idx] = (i as i32) * TERM_GAP_EXT;
        } else {
            mat_d[idx] = open_pen_b[0] + (i as i32 - 1) * gap_ext;
        }
        tr_d[idx] = if i == 1 { 0 } else { 1 };
    }
    for j in 1..=lb {
        if free_leading_gaps {
            mat_i[j] = (j as i32) * TERM_GAP_EXT;
        } else {
            mat_i[j] = open_pen_a[0] + (j as i32 - 1) * gap_ext;
        }
        tr_i[j] = if j == 1 { 0 } else { 1 };
    }

    // DNA200 substitution matrix lookup table `[4][4]` (Katoh et al. 2002 Kimura 200 PAM)
    const DNA200_MAT: [[f64; 4]; 4] = [
        [100.0, -60.0, 24.0, -60.0],
        [-60.0, 100.0, -60.0, 24.0],
        [24.0, -60.0, 100.0, -60.0],
        [-60.0, 24.0, -60.0, 100.0],
    ];

    for i in 1..=la {
        let row = i * cols;
        let prev = (i - 1) * cols;
        let fa = &freq_a[i - 1];
        let op_a = open_pen_a[i - 1];

        for j in 1..=lb {
            let idx = row + j;
            let fb = &freq_b[j - 1];
            let op_b = open_pen_b[j - 1];

            // 1. Update D[i, j] (consume A[i-1], gap in B)
            let d_open = mat_m[prev + j] + op_b;
            let d_ext = mat_d[prev + j] + gap_ext;
            if d_ext >= d_open {
                mat_d[idx] = d_ext;
                tr_d[idx] = 1;
            } else {
                mat_d[idx] = d_open;
                tr_d[idx] = 0;
            }

            // 2. Update I[i, j] (gap in A, consume B[j-1])
            let i_open = mat_m[row + j - 1] + op_a;
            let i_ext = mat_i[row + j - 1] + gap_ext;
            if i_ext >= i_open {
                mat_i[idx] = i_ext;
                tr_i[idx] = 1;
            } else {
                mat_i[idx] = i_open;
                tr_i[idx] = 0;
            }

            // 3. Update M[i, j] (align col A[i-1] with col B[j-1])
            let mut wsp = 0.0f64;
            for x in 0..4 {
                if fa[x] > 0.0 {
                    for y in 0..4 {
                        wsp += fa[x] * fb[y] * DNA200_MAT[x][y];
                    }
                }
            }
            let cb = cons_bonus[(i - 1) * lb + (j - 1)];
            let sub_score = wsp.round() as i32 + cb;

            let diag = prev + j - 1;
            let m_m = mat_m[diag];
            let m_d = mat_d[diag];
            let m_i = mat_i[diag];
            if m_m >= m_d && m_m >= m_i {
                mat_m[idx] = m_m + sub_score;
                tr_m[idx] = 0;
            } else if m_d >= m_i {
                mat_m[idx] = m_d + sub_score;
                tr_m[idx] = 1;
            } else {
                mat_m[idx] = m_i + sub_score;
                tr_m[idx] = 2;
            }
        }
    }

    // Traceback with optional free 3' trailing overhang gaps
    let mut steps: Vec<u8> = Vec::with_capacity(la + lb); // 0 = both (M), 1 = A only (D), 2 = B only (I)
    let mut i = la;
    let mut j = lb;

    let mut state = {
        let end_idx = la * cols + lb;
        let mut best_sc = mat_m[end_idx].max(mat_d[end_idx]).max(mat_i[end_idx]);
        let mut best_i = la;
        let mut best_j = lb;
        let mut best_st = if mat_m[end_idx] >= mat_d[end_idx] && mat_m[end_idx] >= mat_i[end_idx] {
            0u8
        } else if mat_d[end_idx] >= mat_i[end_idx] {
            1u8
        } else {
            2u8
        };

        if free_trailing_gaps {
            for j_cand in 0..=lb {
                let idx = la * cols + j_cand;
                let cand = mat_m[idx] + ((lb - j_cand) as i32) * TERM_GAP_EXT;
                if cand > best_sc {
                    best_sc = cand;
                    best_i = la;
                    best_j = j_cand;
                    best_st = 0;
                }
            }
            for i_cand in 0..=la {
                let idx = i_cand * cols + lb;
                let cand = mat_m[idx] + ((la - i_cand) as i32) * TERM_GAP_EXT;
                if cand > best_sc {
                    best_sc = cand;
                    best_i = i_cand;
                    best_j = lb;
                    best_st = 0;
                }
            }
            while j > best_j {
                steps.push(2);
                j -= 1;
            }
            while i > best_i {
                steps.push(1);
                i -= 1;
            }
        }
        best_st
    };

    while i > 0 || j > 0 {
        if i == 0 {
            steps.push(2);
            j -= 1;
            continue;
        }
        if j == 0 {
            steps.push(1);
            i -= 1;
            continue;
        }
        let idx = i * cols + j;
        match state {
            0 => {
                let prev_st = tr_m[idx];
                steps.push(0);
                i -= 1;
                j -= 1;
                state = prev_st;
            }
            1 => {
                let ext = tr_d[idx];
                steps.push(1);
                i -= 1;
                state = if ext == 1 { 1 } else { 0 };
            }
            _ => {
                let ext = tr_i[idx];
                steps.push(2);
                j -= 1;
                state = if ext == 1 { 2 } else { 0 };
            }
        }
    }
    steps.reverse();

    // Materialize merged profile
    let total_cols = steps.len();
    let total_members = na + nb;
    let mut members = Vec::with_capacity(total_members);
    members.extend_from_slice(&prof_a.members);
    members.extend_from_slice(&prof_b.members);

    let mut rows = vec![Vec::with_capacity(total_cols); total_members];
    let mut orig_pos = vec![Vec::with_capacity(total_cols); total_members];

    let mut cur_a = 0usize;
    let mut cur_b = 0usize;
    for st in steps {
        match st {
            0 => {
                for ma in 0..na {
                    rows[ma].push(prof_a.rows[ma][cur_a]);
                    orig_pos[ma].push(prof_a.orig_pos[ma][cur_a]);
                }
                for mb in 0..nb {
                    rows[na + mb].push(prof_b.rows[mb][cur_b]);
                    orig_pos[na + mb].push(prof_b.orig_pos[mb][cur_b]);
                }
                cur_a += 1;
                cur_b += 1;
            }
            1 => {
                for ma in 0..na {
                    rows[ma].push(prof_a.rows[ma][cur_a]);
                    orig_pos[ma].push(prof_a.orig_pos[ma][cur_a]);
                }
                for mb in 0..nb {
                    rows[na + mb].push(b'-');
                    orig_pos[na + mb].push(None);
                }
                cur_a += 1;
            }
            _ => {
                for ma in 0..na {
                    rows[ma].push(b'-');
                    orig_pos[ma].push(None);
                }
                for mb in 0..nb {
                    rows[na + mb].push(prof_b.rows[mb][cur_b]);
                    orig_pos[na + mb].push(prof_b.orig_pos[mb][cur_b]);
                }
                cur_b += 1;
            }
        }
    }

    AlignmentProfile {
        members,
        rows,
        orig_pos,
    }
}

/// Execute progressive profile-to-profile alignment along the merge steps of `guide_steps`,
/// returning the full alignment ordered by `0..N`.
fn execute_progressive_guide_tree_alignment(
    seqs: &[Vec<u8>],
    guide_steps: &[GuideTreeMergeStep],
    hat3: &SparseHat3ConsistencyMatrix,
    seq_weights: &[f64],
    free_leading_gaps: bool,
    free_trailing_gaps: bool,
    is_conserved_core_zone: bool,
) -> AlignmentProfile {
    let n = seqs.len();
    let mut sub_profiles: Vec<AlignmentProfile> = (0..n)
        .map(|i| AlignmentProfile::from_single_sequence(i, &seqs[i]))
        .collect();

    for step in guide_steps {
        let left_first = step.left_members[0];
        let right_first = step.right_members[0];
        let idx_l = sub_profiles
            .iter()
            .position(|p| p.members.contains(&left_first))
            .unwrap_or(0);
        let prof_l = sub_profiles.remove(idx_l);
        let idx_r = sub_profiles
            .iter()
            .position(|p| p.members.contains(&right_first))
            .unwrap_or(0);
        let prof_r = sub_profiles.remove(idx_r);

        let merged = align_profiles_gotoh_mafft(
            &prof_l,
            &prof_r,
            hat3,
            seq_weights,
            free_leading_gaps,
            free_trailing_gaps,
            is_conserved_core_zone,
        );
        sub_profiles.push(merged);
    }

    let mut root_prof = sub_profiles.pop().unwrap_or_else(|| AlignmentProfile {
        members: Vec::new(),
        rows: Vec::new(),
        orig_pos: Vec::new(),
    });
    canonicalize_profile_member_order(&mut root_prof, n);
    root_prof
}

/// Reorder rows of `prof` so `prof.members == [0, 1, ..., n - 1]`.
fn canonicalize_profile_member_order(prof: &mut AlignmentProfile, n: usize) {
    if prof.members.len() != n {
        return;
    }
    let mut new_rows = vec![Vec::new(); n];
    let mut new_orig = vec![Vec::new(); n];
    for (slot, &m_id) in prof.members.iter().enumerate() {
        if m_id < n {
            new_rows[m_id] = std::mem::take(&mut prof.rows[slot]);
            new_orig[m_id] = std::mem::take(&mut prof.orig_pos[slot]);
        }
    }
    prof.members = (0..n).collect();
    prof.rows = new_rows;
    prof.orig_pos = new_orig;
}

// ============================================================================
// STAGE 6: TREE-DEPENDENT RESTRICTED PARTITIONING ITERATIVE REFINEMENT (`dvtditr`)
// ============================================================================

/// Compute the global `WSP + hat3 Consistency - AffineGap` objective score of `prof`
/// (`dvtditr.exe` objective criterion in MAFFT `L-INS-i` / `FFT-NS-i`).
fn evaluate_profile_objective_score(
    prof: &AlignmentProfile,
    hat3: &SparseHat3ConsistencyMatrix,
    seq_weights: &[f64],
    free_leading_gaps: bool,
    free_trailing_gaps: bool,
    is_conserved_core_zone: bool,
) -> f64 {
    let n = prof.rows.len();
    let l = prof.num_cols();
    if n <= 1 || l == 0 {
        return 0.0;
    }

    let mut total_score = 0.0f64;

    // 1. Column-wise WSP + `hat3` Consistency
    for c in 0..l {
        for i in 0..n {
            let ba = prof.rows[i][c];
            if ba == b'-' {
                continue;
            }
            let seq_a = prof.members[i];
            let pos_a = prof.orig_pos[i][c].unwrap_or(0) as usize;
            let wa = seq_weights[seq_a];

            for j in (i + 1)..n {
                let bb = prof.rows[j][c];
                if bb == b'-' {
                    continue;
                }
                let seq_b = prof.members[j];
                let pos_b = prof.orig_pos[j][c].unwrap_or(0) as usize;
                let wb = seq_weights[seq_b];
                let w_pair = wa * wb;

                let sub = dna200_pair_score(ba, bb) as f64;
                let cons = hat3.get_weight(seq_a, seq_b, pos_a, pos_b) as f64;
                total_score += w_pair * (sub + 0.85 * cons);
            }
        }
    }

    // 2. Pairwise affine gap penalty across internal (non-overhang) gaps
    let gap_open_pen = if is_conserved_core_zone { 660.0 } else { 410.0 };
    let gap_ext_pen = if is_conserved_core_zone { 24.0 } else { 6.0 };

    for i in 0..n {
        let row = &prof.rows[i];
        let w = seq_weights[prof.members[i]];
        let first_non_gap = row.iter().position(|&b| b != b'-').unwrap_or(l);
        let last_non_gap = row.iter().rposition(|&b| b != b'-').unwrap_or(0);

        let start_c = if free_leading_gaps { first_non_gap } else { 0 };
        let end_c = if free_trailing_gaps {
            (last_non_gap + 1).min(l)
        } else {
            l
        };

        let mut in_gap = false;
        for c in start_c..end_c {
            if row[c] == b'-' {
                if !in_gap {
                    total_score -= w * gap_open_pen;
                    in_gap = true;
                } else {
                    total_score -= w * gap_ext_pen;
                }
            } else {
                in_gap = false;
            }
        }
    }

    total_score
}

/// Extract a sub-profile containing `subset_members` from `full_prof` and strip all-gap columns.
fn extract_stripped_subprofile(
    full_prof: &AlignmentProfile,
    subset_members: &[usize],
) -> AlignmentProfile {
    let l = full_prof.num_cols();
    let k = subset_members.len();
    let mut sub_rows = vec![Vec::with_capacity(l); k];
    let mut sub_orig = vec![Vec::with_capacity(l); k];

    let slot_indices: Vec<usize> = subset_members
        .iter()
        .filter_map(|&m| full_prof.members.iter().position(|&x| x == m))
        .collect();

    for c in 0..l {
        let any_residue = slot_indices
            .iter()
            .any(|&slot| full_prof.rows[slot][c] != b'-');
        if any_residue {
            for (out_i, &slot) in slot_indices.iter().enumerate() {
                sub_rows[out_i].push(full_prof.rows[slot][c]);
                sub_orig[out_i].push(full_prof.orig_pos[slot][c]);
            }
        }
    }

    AlignmentProfile {
        members: subset_members.to_vec(),
        rows: sub_rows,
        orig_pos: sub_orig,
    }
}

/// Run Tree-Dependent Restricted Partitioning Iterative Refinement (`dvtditr.exe` in MAFFT `L-INS-i`)
/// across guide-tree bipartitions up to `max_iterations` steps, followed by complementary gap-column
/// compaction.
#[allow(clippy::too_many_arguments)]
fn refine_alignment_iteratively(
    mut prof: AlignmentProfile,
    guide_steps: &[GuideTreeMergeStep],
    hat3: &SparseHat3ConsistencyMatrix,
    seq_weights: &[f64],
    max_iterations: usize,
    free_leading_gaps: bool,
    free_trailing_gaps: bool,
    is_conserved_core_zone: bool,
) -> (AlignmentProfile, usize, usize, f64, f64) {
    let n = prof.rows.len();
    let init_score = evaluate_profile_objective_score(
        &prof,
        hat3,
        seq_weights,
        free_leading_gaps,
        free_trailing_gaps,
        is_conserved_core_zone,
    );
    if n <= 2 || max_iterations == 0 {
        compact_complementary_gap_columns(&mut prof, hat3, is_conserved_core_zone);
        let final_sc = evaluate_profile_objective_score(
            &prof,
            hat3,
            seq_weights,
            free_leading_gaps,
            free_trailing_gaps,
            is_conserved_core_zone,
        );
        return (prof, 0, 0, init_score, final_sc.max(init_score));
    }

    // Collect non-trivial bipartitions `G1` (`1 <= |G1| < n`) induced by guide-tree nodes and leaves
    let mut bipartitions: Vec<Vec<usize>> = Vec::new();
    for step in guide_steps {
        if step.left_members.len() < n {
            bipartitions.push(step.left_members.clone());
        }
        if step.right_members.len() < n {
            bipartitions.push(step.right_members.clone());
        }
    }
    bipartitions.dedup();

    let mut current_score = init_score;
    let mut evaluated = 0usize;
    let mut accepted = 0usize;

    let max_cycles = if max_iterations >= 8 { 2 } else { 1 };
    for _cycle in 0..max_cycles {
        let mut improved_in_cycle = false;
        for g1 in &bipartitions {
            if evaluated >= max_iterations {
                break;
            }
            let g2: Vec<usize> = (0..n).filter(|m| !g1.contains(m)).collect();
            if g1.is_empty() || g2.is_empty() {
                continue;
            }
            evaluated += 1;

            let sub_g1 = extract_stripped_subprofile(&prof, g1);
            let sub_g2 = extract_stripped_subprofile(&prof, &g2);
            let mut candidate = align_profiles_gotoh_mafft(
                &sub_g1,
                &sub_g2,
                hat3,
                seq_weights,
                free_leading_gaps,
                free_trailing_gaps,
                is_conserved_core_zone,
            );
            canonicalize_profile_member_order(&mut candidate, n);

            let cand_score = evaluate_profile_objective_score(
                &candidate,
                hat3,
                seq_weights,
                free_leading_gaps,
                free_trailing_gaps,
                is_conserved_core_zone,
            );
            if cand_score > current_score + 1.0 {
                current_score = cand_score;
                prof = candidate;
                accepted += 1;
                improved_in_cycle = true;
            }
        }
        if !improved_in_cycle || evaluated >= max_iterations {
            break;
        }
    }

    compact_complementary_gap_columns(&mut prof, hat3, is_conserved_core_zone);
    let final_score = evaluate_profile_objective_score(
        &prof,
        hat3,
        seq_weights,
        free_leading_gaps,
        free_trailing_gaps,
        is_conserved_core_zone,
    )
    .max(current_score);

    (prof, evaluated, accepted, init_score, final_score)
}

/// Merge adjacent complementary gap columns `(c, c + 1)` where no sequence has a non-gap
/// nucleotide in BOTH column `c` and column `c + 1` AND merging them pairs identical or
/// high-consistency homologous bases. Also strips any 100%-gap columns.
fn compact_complementary_gap_columns(
    prof: &mut AlignmentProfile,
    hat3: &SparseHat3ConsistencyMatrix,
    is_conserved_core_zone: bool,
) {
    let n = prof.rows.len();
    if n == 0 {
        return;
    }

    let mut c = 0usize;
    while c + 1 < prof.num_cols() {
        let mut disjoint = true;
        for r in 0..n {
            if prof.rows[r][c] != b'-' && prof.rows[r][c + 1] != b'-' {
                disjoint = false;
                break;
            }
        }
        if disjoint {
            // Check whether merging column `c + 1` into column `c` increases match/consistency score
            let mut merge_score = 0i32;
            for i in 0..n {
                let bi = if prof.rows[i][c] != b'-' {
                    prof.rows[i][c]
                } else {
                    prof.rows[i][c + 1]
                };
                if bi == b'-' {
                    continue;
                }
                let pi = prof.orig_pos[i][c].or(prof.orig_pos[i][c + 1]).unwrap_or(0) as usize;
                for j in (i + 1)..n {
                    let bj = if prof.rows[j][c] != b'-' {
                        prof.rows[j][c]
                    } else {
                        prof.rows[j][c + 1]
                    };
                    if bj == b'-' {
                        continue;
                    }
                    let pj = prof.orig_pos[j][c].or(prof.orig_pos[j][c + 1]).unwrap_or(0) as usize;
                    merge_score += dna200_pair_score(bi, bj)
                        + hat3.get_weight(prof.members[i], prof.members[j], pi, pj);
                }
            }
            if merge_score >= 0 || is_conserved_core_zone {
                for r in 0..n {
                    if prof.rows[r][c] == b'-' && prof.rows[r][c + 1] != b'-' {
                        prof.rows[r][c] = prof.rows[r][c + 1];
                        prof.orig_pos[r][c] = prof.orig_pos[r][c + 1];
                    }
                    prof.rows[r].remove(c + 1);
                    prof.orig_pos[r].remove(c + 1);
                }
                continue;
            }
        }
        c += 1;
    }

    // Remove any all-gap columns
    let mut col = 0usize;
    while col < prof.num_cols() {
        let all_gap = (0..n).all(|r| prof.rows[r][col] == b'-');
        if all_gap {
            for r in 0..n {
                prof.rows[r].remove(col);
                prof.orig_pos[r].remove(col);
            }
        } else {
            col += 1;
        }
    }
}

// ============================================================================
// SINGLE-DOMAIN & MULTI-ZONE (`18S–ITS1–5.8S–ITS2–28S`) ORCHESTRATORS
// ============================================================================

#[allow(clippy::too_many_arguments)]
fn align_single_domain_mafft(
    seqs: &[Vec<u8>],
    _strategy: MafftNativeStrategy,
    max_iterations: usize,
    use_transitive_consistency: bool,
    use_secondary_structure_weights: bool,
    free_leading_gaps: bool,
    free_trailing_gaps: bool,
    is_conserved_core_zone: bool,
    t_stage_start: Instant,
) -> (
    Vec<Vec<u8>>,
    usize,
    usize,
    usize,
    usize,
    f64,
    f64,
    f64,
    f64,
    f64,
) {
    let n = seqs.len();
    let max_len = seqs.iter().map(|s| s.len()).max().unwrap_or(0);
    if max_len == 0 {
        return (vec![Vec::new(); n], 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }

    // Fast path for identical-length ultra-conserved core zones (`5.8S` 159-bp core with 0 indels)
    let min_len = seqs.iter().map(|s| s.len()).min().unwrap_or(0);
    if is_conserved_core_zone && min_len == max_len && min_len >= 120 {
        let s12_ms = t_stage_start.elapsed().as_secs_f64() * 1000.0;
        let hat3 = SparseHat3ConsistencyMatrix::new(&vec![min_len; n]);
        let weights = vec![1.0; n];
        let mut prof = AlignmentProfile {
            members: (0..n).collect(),
            rows: seqs.to_vec(),
            orig_pos: (0..n)
                .map(|_| (0..min_len).map(|p| Some(p as u16)).collect())
                .collect(),
        };
        let sc = evaluate_profile_objective_score(
            &prof,
            &hat3,
            &weights,
            false,
            false,
            true,
        );
        let rows = std::mem::take(&mut prof.rows);
        return (rows, n, n * min_len, 0, 0, sc, sc, s12_ms, 0.01, 0.0);
    }

    // Stage 2 & 3: FFT Phase Correlation + All-Pairs Local Smith-Waterman `hat3` Consistency
    let (hat3, dist_mat_1, fft_anchors) = build_hat3_and_distance_matrix(
        seqs,
        use_transitive_consistency,
        use_secondary_structure_weights,
    );
    let hat3_links = hat3.total_links;
    let s12_ms = t_stage_start.elapsed().as_secs_f64() * 1000.0;

    // Stage 4 & 5: Pass 1 Guide Tree (`FFT-NS-1`) -> Pass 2 Refined Guide Tree (`FFT-NS-2`)
    let t_prog = Instant::now();
    let (guide_tree_1, weights_1) = build_mafft_guide_tree_and_weights(&dist_mat_1);
    let pass1_prof = execute_progressive_guide_tree_alignment(
        seqs,
        &guide_tree_1,
        &hat3,
        &weights_1,
        free_leading_gaps,
        free_trailing_gaps,
        is_conserved_core_zone,
    );

    let dist_mat_2 = compute_aligned_distance_matrix(&pass1_prof.rows);
    let (guide_tree_2, weights_2) = build_mafft_guide_tree_and_weights(&dist_mat_2);
    let pass2_prof = execute_progressive_guide_tree_alignment(
        seqs,
        &guide_tree_2,
        &hat3,
        &weights_2,
        free_leading_gaps,
        free_trailing_gaps,
        is_conserved_core_zone,
    );
    let s34_ms = t_prog.elapsed().as_secs_f64() * 1000.0;

    // Stage 6: Iterative Tree-Bipartition Restricted Partitioning Refinement (`dvtditr`)
    let t_iter = Instant::now();
    let (refined_prof, iter_eval, iter_acc, init_wsp, final_wsp) = refine_alignment_iteratively(
        pass2_prof,
        &guide_tree_2,
        &hat3,
        &weights_2,
        max_iterations,
        free_leading_gaps,
        free_trailing_gaps,
        is_conserved_core_zone,
    );
    let s5_ms = t_iter.elapsed().as_secs_f64() * 1000.0;

    (
        refined_prof.rows,
        fft_anchors,
        hat3_links,
        iter_eval,
        iter_acc,
        init_wsp,
        final_wsp,
        s12_ms,
        s34_ms,
        s5_ms,
    )
}

/// Partition an `18S–ITS1–5.8S–ITS2–28S` rDNA cistron across its 5 structural domains:
/// - Zone 0: `18S` 3' flank (`0..its1_start`, free 5' leading gaps, anchored 3' at `ITS1` boundary)
/// - Zone 1: `ITS1` hypervariable spacer (`its1_start..core_58s_start`, free 5' leading gaps for
///           amplicons that lack the `18S` primer tail)
/// - Zone 2: `5.8S` invariant ~158–159 bp core (`core_58s_start..core_58s_end`, locked boundaries)
/// - Zone 3: `ITS2` hypervariable spacer (`core_58s_end..its2_end`, free 3' trailing gaps for
///           amplicons that terminate before the `28S` flank)
/// - Zone 4: `28S` 5' flank (`its2_end..len`, anchored 5' at `ITS2` boundary, free 3' trailing gaps)
#[allow(clippy::too_many_arguments)]
fn align_partitioned_its_cistron_mafft(
    seqs: &[Vec<u8>],
    bounds: &[RdnaItsDomainBounds],
    strategy: MafftNativeStrategy,
    max_iterations: usize,
    use_transitive_consistency: bool,
    use_secondary_structure_weights: bool,
    t_stage_start: Instant,
) -> (
    Vec<Vec<u8>>,
    usize,
    usize,
    usize,
    usize,
    f64,
    f64,
    f64,
    f64,
    f64,
) {
    let n = seqs.len();
    let mut concat_rows = vec![Vec::new(); n];
    let mut total_fft = 0usize;
    let mut total_hat3 = 0usize;
    let mut total_eval = 0usize;
    let mut total_acc = 0usize;
    let mut sum_init_wsp = 0.0f64;
    let mut sum_final_wsp = 0.0f64;
    let mut sum_s12 = t_stage_start.elapsed().as_secs_f64() * 1000.0;
    let mut sum_s34 = 0.0f64;
    let mut sum_s5 = 0.0f64;

    // Zone 0 (`18S` 3' flank): only sequences with `its1_start > 0` have `18S` sequence
    let any_18s = bounds.iter().any(|b| b.its1_start > 0);
    if any_18s {
        let z0_slices: Vec<Vec<u8>> = (0..n)
            .map(|i| seqs[i][..bounds[i].its1_start].to_vec())
            .collect();
        let (rows, fft, h3, ev, ac, iw, fw, ms12, ms34, ms5) = align_single_domain_mafft(
            &z0_slices,
            strategy,
            max_iterations / 2,
            use_transitive_consistency,
            false,
            true,
            false,
            false,
            Instant::now(),
        );
        for i in 0..n {
            concat_rows[i].extend_from_slice(&rows[i]);
        }
        total_fft += fft;
        total_hat3 += h3;
        total_eval += ev;
        total_acc += ac;
        sum_init_wsp += iw;
        sum_final_wsp += fw;
        sum_s12 += ms12;
        sum_s34 += ms34;
        sum_s5 += ms5;
    }

    // Zone 1 (`ITS1` hypervariable spacer: `its1_start..core_58s_start`)
    // Both 5' (`18S/ITS1` junction, since Zone 0 already absorbed any 5' `18S` overhang)
    // and 3' (`ITS1/5.8S` junction) boundaries are anchored.
    let z1_slices: Vec<Vec<u8>> = (0..n)
        .map(|i| seqs[i][bounds[i].its1_start..bounds[i].core_58s_start].to_vec())
        .collect();
    let (z1_rows, fft1, h3_1, ev1, ac1, iw1, fw1, ms12_1, ms34_1, ms5_1) =
        align_single_domain_mafft(
            &z1_slices,
            strategy,
            max_iterations,
            use_transitive_consistency,
            use_secondary_structure_weights,
            false,
            false,
            false,
            Instant::now(),
        );
    for i in 0..n {
        concat_rows[i].extend_from_slice(&z1_rows[i]);
    }
    total_fft += fft1;
    total_hat3 += h3_1;
    total_eval += ev1;
    total_acc += ac1;
    sum_init_wsp += iw1;
    sum_final_wsp += fw1;
    sum_s12 += ms12_1;
    sum_s34 += ms34_1;
    sum_s5 += ms5_1;

    // Zone 2 (`5.8S` invariant rDNA core: `core_58s_start..core_58s_end`)
    let z2_slices: Vec<Vec<u8>> = (0..n)
        .map(|i| seqs[i][bounds[i].core_58s_start..bounds[i].core_58s_end].to_vec())
        .collect();
    let (z2_rows, fft2, h3_2, ev2, ac2, iw2, fw2, ms12_2, ms34_2, ms5_2) =
        align_single_domain_mafft(
            &z2_slices,
            strategy,
            max_iterations / 2,
            use_transitive_consistency,
            true,
            false,
            false,
            true,
            Instant::now(),
        );
    for i in 0..n {
        concat_rows[i].extend_from_slice(&z2_rows[i]);
    }
    total_fft += fft2;
    total_hat3 += h3_2;
    total_eval += ev2;
    total_acc += ac2;
    sum_init_wsp += iw2;
    sum_final_wsp += fw2;
    sum_s12 += ms12_2;
    sum_s34 += ms34_2;
    sum_s5 += ms5_2;

    // Zone 3 (`ITS2` hypervariable spacer: `core_58s_end..its2_end`)
    // Both 5' (`5.8S/ITS2` junction) and 3' (`ITS2/28S` junction, since Zone 4 absorbs any 3' `28S` overhang)
    // boundaries are anchored.
    let z3_slices: Vec<Vec<u8>> = (0..n)
        .map(|i| seqs[i][bounds[i].core_58s_end..bounds[i].its2_end].to_vec())
        .collect();
    let (z3_rows, fft3, h3_3, ev3, ac3, iw3, fw3, ms12_3, ms34_3, ms5_3) =
        align_single_domain_mafft(
            &z3_slices,
            strategy,
            max_iterations,
            use_transitive_consistency,
            use_secondary_structure_weights,
            false,
            false,
            false,
            Instant::now(),
        );
    for i in 0..n {
        concat_rows[i].extend_from_slice(&z3_rows[i]);
    }
    total_fft += fft3;
    total_hat3 += h3_3;
    total_eval += ev3;
    total_acc += ac3;
    sum_init_wsp += iw3;
    sum_final_wsp += fw3;
    sum_s12 += ms12_3;
    sum_s34 += ms34_3;
    sum_s5 += ms5_3;

    // Zone 4 (`28S` 5' flank: `its2_end..len`)
    let any_28s = (0..n).any(|i| bounds[i].its2_end < seqs[i].len());
    if any_28s {
        let z4_slices: Vec<Vec<u8>> = (0..n)
            .map(|i| seqs[i][bounds[i].its2_end..].to_vec())
            .collect();
        let (rows, fft, h3, ev, ac, iw, fw, ms12, ms34, ms5) = align_single_domain_mafft(
            &z4_slices,
            strategy,
            max_iterations / 2,
            use_transitive_consistency,
            false,
            false,
            true,
            false,
            Instant::now(),
        );
        for i in 0..n {
            concat_rows[i].extend_from_slice(&rows[i]);
        }
        total_fft += fft;
        total_hat3 += h3;
        total_eval += ev;
        total_acc += ac;
        sum_init_wsp += iw;
        sum_final_wsp += fw;
        sum_s12 += ms12;
        sum_s34 += ms34;
        sum_s5 += ms5;
    }

    (
        concat_rows,
        total_fft,
        total_hat3,
        total_eval,
        total_acc,
        sum_init_wsp,
        sum_final_wsp,
        sum_s12,
        sum_s34,
        sum_s5,
    )
}

/// Align a multi-zone partitioned cohort (`18S SSU rRNA`, `HOXD13` Homeobox + repeat + UTR,
/// or `RbcS` transit peptide + mature core + UTR) seamlessly across all zones while preserving
/// 100% of every sequence's full nucleotide length.
#[allow(clippy::too_many_arguments)]
fn align_partitioned_multi_zone_mafft(
    seqs: &[Vec<u8>],
    zones: &[MultiZoneSliceSpec],
    strategy: MafftNativeStrategy,
    max_iterations: usize,
    use_transitive_consistency: bool,
    use_secondary_structure_weights: bool,
    t_stage_start: Instant,
) -> (
    Vec<Vec<u8>>,
    usize,
    usize,
    usize,
    usize,
    f64,
    f64,
    f64,
    f64,
    f64,
) {
    let n = seqs.len();
    let mut concat_rows = vec![Vec::new(); n];
    let mut total_fft = 0usize;
    let mut total_hat3 = 0usize;
    let mut total_eval = 0usize;
    let mut total_acc = 0usize;
    let mut sum_init_wsp = 0.0f64;
    let mut sum_final_wsp = 0.0f64;
    let mut sum_s12 = t_stage_start.elapsed().as_secs_f64() * 1000.0;
    let mut sum_s34 = 0.0f64;
    let mut sum_s5 = 0.0f64;

    for z in zones {
        let full_slices: Vec<&[u8]> = (0..n)
            .map(|i| {
                let (s, e) = z.ranges[i];
                &seqs[i][s..e]
            })
            .collect();
        if full_slices.iter().all(|sl| sl.is_empty()) {
            continue;
        }

        // Check for extreme 5' or 3' terminal UTR overhangs when `terminal_overhang_cap` is set
        if let Some(cap) = z.terminal_overhang_cap {
            if z.free_leading_gaps && !z.free_trailing_gaps {
                // 5' leader zone: left-pad any distal 5' leader beyond `cap`, align proximal `cap` bp
                let max_distal = full_slices
                    .iter()
                    .map(|sl| sl.len().saturating_sub(cap))
                    .max()
                    .unwrap_or(0);
                if max_distal > 0 {
                    for i in 0..n {
                        let distal_len = full_slices[i].len().saturating_sub(cap);
                        let pad = max_distal - distal_len;
                        concat_rows[i].extend(std::iter::repeat(b'-').take(pad));
                        concat_rows[i].extend_from_slice(&full_slices[i][..distal_len]);
                    }
                }
                let prox_slices: Vec<Vec<u8>> = full_slices
                    .iter()
                    .map(|sl| sl[sl.len().saturating_sub(cap)..].to_vec())
                    .collect();
                let (rows, fft, h3, ev, ac, iw, fw, ms12, ms34, ms5) = align_single_domain_mafft(
                    &prox_slices,
                    strategy,
                    max_iterations,
                    use_transitive_consistency,
                    use_secondary_structure_weights,
                    true,
                    false,
                    false,
                    Instant::now(),
                );
                for i in 0..n {
                    concat_rows[i].extend_from_slice(&rows[i]);
                }
                total_fft += fft;
                total_hat3 += h3;
                total_eval += ev;
                total_acc += ac;
                sum_init_wsp += iw;
                sum_final_wsp += fw;
                sum_s12 += ms12;
                sum_s34 += ms34;
                sum_s5 += ms5;
                continue;
            } else if !z.free_leading_gaps && z.free_trailing_gaps {
                // 3' trailer zone: align proximal `cap` bp immediately downstream of core, right-pad distal 3' tail
                let prox_slices: Vec<Vec<u8>> = full_slices
                    .iter()
                    .map(|sl| sl[..sl.len().min(cap)].to_vec())
                    .collect();
                let (rows, fft, h3, ev, ac, iw, fw, ms12, ms34, ms5) = align_single_domain_mafft(
                    &prox_slices,
                    strategy,
                    max_iterations,
                    use_transitive_consistency,
                    use_secondary_structure_weights,
                    false,
                    true,
                    false,
                    Instant::now(),
                );
                for i in 0..n {
                    concat_rows[i].extend_from_slice(&rows[i]);
                }
                let max_distal = full_slices
                    .iter()
                    .map(|sl| sl.len().saturating_sub(cap))
                    .max()
                    .unwrap_or(0);
                if max_distal > 0 {
                    for i in 0..n {
                        let tail = &full_slices[i][full_slices[i].len().min(cap)..];
                        concat_rows[i].extend_from_slice(tail);
                        let pad = max_distal - tail.len();
                        concat_rows[i].extend(std::iter::repeat(b'-').take(pad));
                    }
                }
                total_fft += fft;
                total_hat3 += h3;
                total_eval += ev;
                total_acc += ac;
                sum_init_wsp += iw;
                sum_final_wsp += fw;
                sum_s12 += ms12;
                sum_s34 += ms34;
                sum_s5 += ms5;
                continue;
            }
        }

        let z_slices: Vec<Vec<u8>> = full_slices.iter().map(|sl| sl.to_vec()).collect();
        let (rows, fft, h3, ev, ac, iw, fw, ms12, ms34, ms5) = align_single_domain_mafft(
            &z_slices,
            strategy,
            if z.is_conserved_core_zone {
                max_iterations / 2
            } else {
                max_iterations
            },
            use_transitive_consistency,
            use_secondary_structure_weights,
            z.free_leading_gaps,
            z.free_trailing_gaps,
            z.is_conserved_core_zone,
            Instant::now(),
        );
        for i in 0..n {
            concat_rows[i].extend_from_slice(&rows[i]);
        }
        total_fft += fft;
        total_hat3 += h3;
        total_eval += ev;
        total_acc += ac;
        sum_init_wsp += iw;
        sum_final_wsp += fw;
        sum_s12 += ms12;
        sum_s34 += ms34;
        sum_s5 += ms5;
    }

    (
        concat_rows,
        total_fft,
        total_hat3,
        total_eval,
        total_acc,
        sum_init_wsp,
        sum_final_wsp,
        sum_s12,
        sum_s34,
        sum_s5,
    )
}

// ============================================================================
// 6-MER TUPLE DISTANCE HELPERS (KATOH ET AL. 2002)
// ============================================================================

fn extract_6mer_multiset(seq: &[u8]) -> Vec<u16> {
    let mut counts = vec![0u16; 4096];
    if seq.len() < 6 {
        return counts;
    }
    for win in seq.windows(6) {
        if let Some(code) = encode_6mer(win) {
            counts[code] = counts[code].saturating_add(1);
        }
    }
    counts
}

#[inline]
fn encode_6mer(win: &[u8]) -> Option<usize> {
    let mut code = 0usize;
    for &b in win {
        let bits = match b {
            b'A' => 0usize,
            b'C' => 1usize,
            b'G' => 2usize,
            b'T' | b'U' => 3usize,
            _ => return None,
        };
        code = (code << 2) | bits;
    }
    Some(code)
}

#[inline]
fn shared_6mer_similarity(a: &[u16], b: &[u16]) -> f64 {
    let mut shared = 0usize;
    let mut sum_a = 0usize;
    let mut sum_b = 0usize;
    for i in 0..4096 {
        let ca = a[i] as usize;
        let cb = b[i] as usize;
        shared += ca.min(cb);
        sum_a += ca;
        sum_b += cb;
    }
    let denom = sum_a.min(sum_b).max(1) as f64;
    (shared as f64 / denom).clamp(0.0, 1.0)
}

fn select_6mer_medoid(seqs: &[Vec<u8>]) -> usize {
    if seqs.len() <= 2 {
        return 0;
    }
    let kmers: Vec<Vec<u16>> = seqs.iter().map(|s| extract_6mer_multiset(s)).collect();
    let mut best_i = 0usize;
    let mut best_sum = -1.0f64;
    for i in 0..seqs.len() {
        let mut s = 0.0f64;
        for j in 0..seqs.len() {
            if i != j {
                s += shared_6mer_similarity(&kmers[i], &kmers[j]);
            }
        }
        if s > best_sum {
            best_sum = s;
            best_i = i;
        }
    }
    best_i
}

/// Standalone CLI runner for `bio-suite mafft-native` / `bio-suite runie-mafft-linsi` / `bio-suite its-align`.
/// Parses a multi-FASTA file (or defaults to the cross-phylum Dikarya `ITS1–5.8S–ITS2` cistron cohort
/// from `crvscience-labs/cache/`), executes the 6-stage Native Rust MAFFT-Mimetic Aligner, writes
/// aligned FASTA output if `--out <path>` is specified, and prints stage-level telemetry.
pub fn run_mafft_native_cli(
    input_fasta: Option<&str>,
    strategy_arg: Option<&str>,
    max_iterate: Option<usize>,
    output_fasta: Option<&str>,
) -> Result<MafftNativeStageTelemetry, String> {
    let records: Vec<FastaRecord> = if let Some(path) = input_fasta {
        let raw = std::fs::read_to_string(path)
            .map_err(|e| format!("Failed to read input FASTA '{}': {}", path, e))?;
        parse_multi_fasta_records(&raw)
    } else {
        load_default_dikarya_its_records()?
    };

    if records.is_empty() {
        return Err("No FASTA records found to align.".to_string());
    }

    let mut cfg = MafftNativeConfig::default();
    if let Some(st) = strategy_arg {
        cfg.strategy = MafftNativeStrategy::from_cli_arg(st);
    }
    if let Some(mi) = max_iterate {
        cfg.max_iterations = mi;
    }

    let (engine_label, aligned, telemetry) = align_sequences_mafft_native_with_config(&records, &cfg);

    println!("================================================================================");
    println!("  CRV SCIENCE LABS — NATIVE RUST MAFFT-MIMETIC ALIGNER (runie-mafft-linsi)");
    println!("================================================================================");
    println!("  Engine            : {}", engine_label);
    println!("  Resolved Strategy : {}", telemetry.resolved_strategy);
    println!(
        "  ITS Cistron Lock  : {} (5.8S anchors: {}/{}, RevComp rescued: {})",
        telemetry.is_its_cistron_detected,
        telemetry.conserved_58s_anchor_count,
        records.len(),
        telemetry.reverse_complemented_count
    );
    println!(
        "  FFT Phase Anchors : {} | Local hat3 Consistency Pairs: {}",
        telemetry.fft_phase_anchors_detected, telemetry.local_hat3_residue_pairs
    );
    println!(
        "  Refinement (dvtd) : {}/{} steps accepted | WSP {:.1} -> {:.1}",
        telemetry.iterative_refinement_steps_accepted,
        telemetry.iterative_refinement_steps_evaluated,
        telemetry.initial_progressive_wsp_score,
        telemetry.final_refined_wsp_score
    );
    println!(
        "  Aligned Columns L : {} bp across {} sequences ({:.3} ms total)",
        aligned.first().map(|a| a.aligned_seq.len()).unwrap_or(0),
        aligned.len(),
        telemetry.total_duration_ms
    );
    println!("================================================================================");

    if let Some(out_path) = output_fasta {
        let mut fasta_out = String::new();
        for seq in &aligned {
            fasta_out.push_str(&format!(
                ">{} {}\n{}\n",
                seq.id, seq.taxon, seq.aligned_seq
            ));
        }
        if let Some(parent) = std::path::Path::new(out_path).parent() {
            let _ = std::fs::create_dir_all(parent);
        }
        std::fs::write(out_path, fasta_out)
            .map_err(|e| format!("Failed to write aligned FASTA '{}': {}", out_path, e))?;
        println!("[OK] Wrote aligned FASTA to {}", out_path);
    }

    Ok(telemetry)
}

fn parse_multi_fasta_records(raw: &str) -> Vec<FastaRecord> {
    let mut records = Vec::new();
    let mut cur_id = String::new();
    let mut cur_taxon = String::new();
    let mut cur_seq = String::new();

    for line in raw.lines() {
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        if let Some(rest) = t.strip_prefix('>') {
            if !cur_seq.is_empty() {
                records.push(FastaRecord {
                    id: cur_id.clone(),
                    accession: cur_id.clone(),
                    taxon: cur_taxon.clone(),
                    common_name: cur_taxon.clone(),
                    family: "Fungi / Eukaryota".to_string(),
                    gene: "ITS1-5.8S-ITS2".to_string(),
                    description: cur_taxon.clone(),
                    sequence: cur_seq.clone(),
                    cached_from_disk: true,
                });
                cur_seq.clear();
            }
            let mut parts = rest.splitn(2, char::is_whitespace);
            cur_id = parts.next().unwrap_or("SEQ").to_string();
            cur_taxon = parts.next().unwrap_or(&cur_id).trim().to_string();
        } else {
            for ch in t.chars() {
                if ch.is_ascii_alphabetic() || ch == '-' {
                    cur_seq.push(ch.to_ascii_uppercase());
                }
            }
        }
    }
    if !cur_seq.is_empty() {
        records.push(FastaRecord {
            id: cur_id.clone(),
            accession: cur_id,
            taxon: cur_taxon.clone(),
            common_name: cur_taxon.clone(),
            family: "Fungi / Eukaryota".to_string(),
            gene: "ITS1-5.8S-ITS2".to_string(),
            description: cur_taxon,
            sequence: cur_seq,
            cached_from_disk: true,
        });
    }
    records
}

fn load_default_dikarya_its_records() -> Result<Vec<FastaRecord>, String> {
    let cache_dir = if std::path::Path::new("crvscience-labs/cache").exists() {
        std::path::Path::new("crvscience-labs/cache")
    } else {
        std::path::Path::new("cache")
    };
    let cache = crate::fetch::PoliteCacheManager::new(cache_dir);
    let specs = [
        ("NR_111041.1", "Aspergillus flavus", "Trichocomaceae"),
        ("NR_121481.1", "Aspergillus fumigatus", "Trichocomaceae"),
        ("NR_138295.1", "Aspergillus movilensis", "Trichocomaceae"),
        ("NR_132055.1", "Colletotrichum arxii", "Glomerellaceae"),
        ("NR_130658.1", "Cryptococcus neoformans", "Cryptococcaceae"),
    ];
    let mut out = Vec::with_capacity(specs.len());
    for (acc, taxon, fam) in specs {
        let rec = cache.fetch_ncbi_fasta_slice(acc, taxon, taxon, fam, "ITS1-5.8S-ITS2")?;
        out.push(rec);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fft_radix2_roundtrip() {
        let orig = [1.0, -0.5, 0.25, 0.75, -1.0, 0.0, 0.5, -0.25];
        let mut re = orig.to_vec();
        let mut im = vec![0.0; 8];
        fft_radix2_in_place(&mut re, &mut im, false);
        fft_radix2_in_place(&mut re, &mut im, true);
        for i in 0..8 {
            assert!((re[i] - orig[i]).abs() < 1e-9);
            assert!(im[i].abs() < 1e-9);
        }
    }

    #[test]
    fn test_rdna_its_domain_detector_and_alignment() {
        let seq_a = b"TCCGTAGGTGAACCTGCGGAAGGATCATTACCGAGTGTAGGGTTCCTAGCGAGCCCAACCTCCCACCCGTGTTTACTGTACCTTAGTTGCTTCGGCGGGCCCGCCATTCATGGCCGCCGGGGGCTCTCAGCCCCGGGCCCGCGCCCGCCGGAGACACCACGAACTCTGTCTGATCTAGTGAATTCTGAGTTGATTGTATCGCAATCAGTTAAAACTTTCAACAATGGATCTCTTGGTTCAGGCATCGATGAAGAACGCAGCGAAATGCGATAACTAGTGTGAATTGCAGAATTCCGTGAATCATCGAGTCTTTGAACGCACATTGCGCCCCCTGGTATTCCGGGGGGCATGCCTGTCCGAGCGTCATTGCTGCCCATCAAGCACGGCTTGTGTGTTGGGTCGTCGTCCCCTCTCCGGGGGGGACGGGCCCCAAAGGCAGCGGCGGCACCGCGTCCGATCCTCGAGCGTATGGGGCTTTGTCACCCGCTCTGTAGGCCCGGCCGGCGCTTGCCGAACGCAAATCAATCTTTTTCCAGGTTGACCTCGGATCAGGTAGGGATACCCGCTGAACTTAAGCATATCAATAAGCGGAGGA";
        let bounds = detect_rdna_its_domain_boundaries(seq_a);
        assert!(bounds.has_58s_core);
        assert_eq!(bounds.its1_start, 29);
        assert_eq!(&seq_a[bounds.core_58s_start..bounds.core_58s_start + 8], b"AAAACTTT");
        assert_eq!(bounds.core_58s_end - bounds.core_58s_start, 159);
        assert_eq!(&seq_a[bounds.its2_end..bounds.its2_end + 8], b"TTGACCTC");
    }
}


#[cfg(test)]
mod strand_polarity_tests {
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

    /// Same-orientation records with no rDNA motifs must never be flipped by a chance motif hit
    /// on the reverse strand (regression: 24 of 79 sense-strand 18S records were flipped).
    #[test]
    fn sense_records_without_motifs_stay_forward() {
        for seed in 1..40u64 {
            let template = pseudo(1800, seed);
            let seqs: Vec<Vec<u8>> = (0..6)
                .map(|i| {
                    let mut s = template.clone();
                    for p in (i * 7..s.len()).step_by(97) {
                        s[p] = b"ACGT"[(p + i) % 4];
                    }
                    s
                })
                .collect();
            let (_, flipped) = normalize_strand_polarity_mafft(&seqs);
            assert!(flipped.iter().all(|f| !f), "seed {seed}: flipped {flipped:?}");
        }
    }

    /// A genuinely reverse-complemented record is still turned around.
    #[test]
    fn reverse_record_is_flipped() {
        let template = pseudo(1500, 11);
        let mut seqs = vec![template.clone(), template.clone(), template.clone()];
        seqs.push(reverse_complement_bytes(&template));
        let (_, flipped) = normalize_strand_polarity_mafft(&seqs);
        assert_eq!(flipped, vec![false, false, false, true]);
    }
}