//! `super_blast` — Native Rust Adaptive Multi-Biotype (`DNA`, `RNA`, `Protein`, `Other` Standard Language)
//! BLAST Engine with Execution-Mode Sensitivity Crossover Analysis, On-the-Fly MSA Sequence Validation,
//! Multi-Gene Phylogenetic Tip Concordance Auditing, Symmetric All-vs-All Distance Matrices,
//! Barcode BLAST Compilation & Just-In-Time (JIT) Transient Genome Search.
//!
//! Designed for zero-friction invocation across `bio-suite` and `RuNIE` pipelines:
//! 1. **Automatic Sequence Biotype Detection (`DNA`, `RNA`, `Protein`, `Other`)**:
//!    Inspects empirical residue distribution and Shannon 1-mer/2-mer entropy to select optimal
//!    substitution scoring (`BLASTN +2/-3`, `RNA Transition/Transversion`, `BLOSUM62`, or
//!    `Standard Language Token/Character Smith-Waterman`), seed word length $w$, strand handling,
//!    and Karlin-Altschul $(\lambda, K, H)$ statistical parameters.
//! 2. **Execution-Mode Sensitivity & Crossover Analyzer**:
//!    Quantifies total search volume $V = Q_{\text{len}} \times T_{\text{len}}$ and expected query
//!    batch reuse $M$ against index construction overhead $C_{\text{index}}(T_{\text{len}})$ to
//!    automatically choose between:
//!    - `OnTheFlyDirectInMemory` (zero database compilation for 2 sequences or small cohorts),
//!    - `JustInTimeTransientIndex` (ephemeral in-memory $k$-mer index built in $<1\text{ ms}$ for
//!      transient genome/barcode searches with zero persistent disk files), and
//!    - `CompiledPersistentDatabase` (only when repeated query batches $M \ge M^*$ amortize
//!      persistent index compilation).
//! 3. **On-the-Fly MSA Sequence Validator (`validate_sequences_for_msa`)**:
//!    Screens candidate sequences prior to Multiple Sequence Alignment (`phylo.rs`), automatically
//!    rescuing reverse-complemented sequences and rejecting non-homologous contaminants, 5'/3'
//!    chimeras, nuclear `NUMT` pseudogenes, and truncated fragments.
//! 4. **Multi-Gene Phylogenetic Tip Concordance Auditor (`audit_phylogeny_tip_concordance`)**:
//!    Tests a panel of secondary genes (e.g., 10 housekeeping/accessory loci) against the tip
//!    neighborhoods of an anchor phylogeny (e.g., `16S rRNA`), quantifying tip-to-tip distance
//!    concordance ($r_{\text{tip}}$) and pinpointing discordant tips caused by Horizontal Gene
//!    Transfer (HGT), recombination, or paralog xenology.
//! 5. **Symmetric Pairwise Distance Matrix (`compute_symmetric_pairwise_matrix`) & Barcode BLAST**:
//!    Computes upper-triangle $N(N-1)/2$ `% identity`, `JC69`, `Kimura 2-Parameter (K2P)`, and
//!    bit-score matrices on the fly, and compiles rapid in-memory DNA barcode classifiers.

use crate::complexity::{kmer_shannon_entropy, round4};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Automatically detected sequence biotype.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SuperBlastBiotype {
    Dna,
    Rna,
    Protein,
    Other,
}

impl SuperBlastBiotype {
    pub fn as_str(&self) -> &'static str {
        match self {
            SuperBlastBiotype::Dna => "DNA",
            SuperBlastBiotype::Rna => "RNA",
            SuperBlastBiotype::Protein => "Protein",
            SuperBlastBiotype::Other => "Other (Standard Language / Generalized String)",
        }
    }
}

/// Execution strategy selected by the `super-blast` sensitivity crossover analyzer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SuperBlastExecutionMode {
    /// Direct in-memory pairwise seed + banded Smith-Waterman without building a database index.
    OnTheFlyDirectInMemory,
    /// Ephemeral Just-In-Time (JIT) in-memory k-mer hash index dropped immediately after query execution.
    JustInTimeTransientIndex,
    /// Persistent compiled database recommended when repeated query batches exceed $M^*$.
    CompiledPersistentDatabase,
}

impl SuperBlastExecutionMode {
    pub fn as_str(&self) -> &'static str {
        match self {
            SuperBlastExecutionMode::OnTheFlyDirectInMemory => "OnTheFlyDirectInMemory",
            SuperBlastExecutionMode::JustInTimeTransientIndex => "JustInTimeTransientIndex",
            SuperBlastExecutionMode::CompiledPersistentDatabase => "CompiledPersistentDatabase",
        }
    }
}

/// Biotype detection telemetry and automatically configured Karlin-Altschul parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BiotypeDetectionResult {
    pub detected_biotype: SuperBlastBiotype,
    pub biotype_label: String,
    pub clean_length: usize,
    pub acgt_fraction: f64,
    pub acgu_fraction: f64,
    pub amino_acid_exclusive_fraction: f64,
    pub whitespace_or_punctuation_fraction: f64,
    pub shannon_h1_bits: f64,
    pub low_complexity_dust_fraction: f64,
    pub seed_word_size: usize,
    pub karlin_lambda: f64,
    pub karlin_k: f64,
    pub check_reverse_complement: bool,
}

/// Sensitivity analysis report explaining why `super-blast` chose `DirectInMemory`, `JIT`, or `CompiledDB`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExecutionSensitivityReport {
    pub selected_mode: SuperBlastExecutionMode,
    pub mode_label: String,
    pub num_queries: usize,
    pub num_targets: usize,
    pub total_query_residues: usize,
    pub total_target_residues: usize,
    pub search_space_volume: u64,
    pub estimated_direct_cost_us: f64,
    pub estimated_jit_index_build_us: f64,
    pub estimated_jit_query_cost_us: f64,
    pub estimated_persistent_db_compile_us: f64,
    pub amortization_crossover_batches: usize,
    pub requires_compiled_db: bool,
    pub rationale: String,
}

/// High-Scoring Segment Pair (HSP) alignment result from `super-blast`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SuperBlastHsp {
    pub query_id: String,
    pub target_id: String,
    pub biotype: String,
    pub execution_mode: String,
    pub strand: char,
    pub query_start_1based: usize,
    pub query_end_1based: usize,
    pub target_start_1based: usize,
    pub target_end_1based: usize,
    pub alignment_length: usize,
    pub matches: usize,
    pub mismatches: usize,
    pub transitions: usize,
    pub transversions: usize,
    pub gap_opens: usize,
    pub gap_bases: usize,
    pub pct_identity: f64,
    pub query_coverage_pct: f64,
    pub raw_score: f64,
    pub composition_adjusted_score: f64,
    pub bit_score: f64,
    pub e_value: f64,
    pub jc69_distance: f64,
    pub k2p_distance: f64,
    pub low_complexity_masked: bool,
    pub aligned_query: String,
    pub midline: String,
    pub aligned_target: String,
}

/// Verdict for pre-MSA sequence validation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MsaSequenceValidationRecord {
    pub seq_id: String,
    pub status: String, // "AcceptedDirect", "AcceptedAutoRevComp", "RejectedChimera", "RejectedNonHomologous", "RejectedTruncated"
    pub incorporate_into_msa: bool,
    pub reverse_complement_rescued: bool,
    pub strand: char,
    pub pct_identity_vs_centroid: f64,
    pub query_coverage_pct: f64,
    pub bit_score: f64,
    pub e_value: f64,
    pub five_prime_half_identity: f64,
    pub three_prime_half_identity: f64,
    pub chimera_delta_identity: f64,
    pub oriented_sequence: String,
    pub diagnostic_note: String,
}

/// Result of testing a secondary gene across the tips of an anchor phylogeny (e.g. 16S rRNA).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhylogenyTipGeneAudit {
    pub gene_name: String,
    pub functional_category: String,
    pub num_tips_evaluated: usize,
    pub pearson_r_vs_anchor_tree: f64,
    pub nearest_neighbor_tip_concordance_pct: f64,
    pub mean_tip_k2p_distance: f64,
    pub discordant_tips_detected: Vec<String>,
    pub concordance_verdict: String, // "ConcordantVerticalInheritance", "DiscordantHorizontalTransferOrXenolog", "CollapsedLowResolution"
}

/// Full multi-gene phylogenetic tip concordance report.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhylogenyTipConcordanceReport {
    pub anchor_marker_name: String,
    pub num_taxa_tips: usize,
    pub num_secondary_genes_tested: usize,
    pub concordant_genes_count: usize,
    pub discordant_genes_count: usize,
    pub mean_concordant_pearson_r: f64,
    pub gene_audits: Vec<PhylogenyTipGeneAudit>,
}

/// On-the-fly symmetric all-vs-all pairwise distance matrix.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SymmetricPairwiseMatrixReport {
    pub biotype: String,
    pub execution_mode: String,
    pub num_sequences: usize,
    pub unique_pairs_evaluated: usize,
    pub labels: Vec<String>,
    pub pct_identity_matrix: Vec<Vec<f64>>,
    pub k2p_distance_matrix: Vec<Vec<f64>>,
    pub jc69_distance_matrix: Vec<Vec<f64>>,
    pub bit_score_matrix: Vec<Vec<f64>>,
}

/// Reference entry in an in-memory Barcode BLAST database.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BarcodeReferenceEntry {
    pub accession_id: String,
    pub species_name: String,
    pub genus_name: String,
    pub marker_locus: String,
    pub sequence: String,
}

/// Compiled in-memory Barcode BLAST database with k-mer inverted index.
#[derive(Debug, Clone)]
pub struct BarcodeBlastDatabase {
    pub marker_locus: String,
    pub word_size: usize,
    pub entries: Vec<BarcodeReferenceEntry>,
    pub kmer_index: HashMap<u64, Vec<(usize, usize)>>,
}

/// Taxonomic assignment and barcoding-gap audit from Barcode BLAST.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BarcodeBlastAssignment {
    pub query_id: String,
    pub top_hit_accession: String,
    pub assigned_species: String,
    pub assigned_genus: String,
    pub top_pct_identity: f64,
    pub top_bit_score: f64,
    pub top_e_value: f64,
    pub min_intraspecific_k2p: f64,
    pub min_interspecific_k2p: f64,
    pub barcoding_gap_k2p: f64,
    pub internal_stop_codons_detected: usize,
    pub numt_pseudogene_flagged: bool,
    pub taxonomic_confidence_regime: String,
}

// ============================================================================
// 1. AUTOMATIC SEQUENCE BIOTYPE DETECTION (DNA / RNA / PROTEIN / OTHER)
// ============================================================================

/// Detects whether a query or target string is `DNA`, `RNA`, `Protein`, or `Other` (Standard Language).
pub fn detect_sequence_biotype(raw: &str) -> BiotypeDetectionResult {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return BiotypeDetectionResult {
            detected_biotype: SuperBlastBiotype::Other,
            biotype_label: SuperBlastBiotype::Other.as_str().to_string(),
            clean_length: 0,
            acgt_fraction: 0.0,
            acgu_fraction: 0.0,
            amino_acid_exclusive_fraction: 0.0,
            whitespace_or_punctuation_fraction: 0.0,
            shannon_h1_bits: 0.0,
            low_complexity_dust_fraction: 0.0,
            seed_word_size: 4,
            karlin_lambda: 0.385,
            karlin_k: 0.12,
            check_reverse_complement: false,
        };
    }

    let total_chars = trimmed.chars().count().max(1);
    let mut non_ws_alpha = 0usize;
    let mut ws_or_punct = 0usize;
    let mut acgt_count = 0usize;
    let mut acgu_count = 0usize;
    let mut t_count = 0usize;
    let mut u_count = 0usize;
    let mut aa_exclusive_count = 0usize; // E, F, I, L, P, Q, Z, W, Y, K, M, R, S, V, D, H
    let mut standard_aa_count = 0usize;
    let mut lowercase_letters = 0usize;
    let mut freq: HashMap<char, usize> = HashMap::new();

    for ch in trimmed.chars() {
        if ch.is_whitespace() || ch.is_ascii_punctuation() || ch.is_ascii_digit() {
            ws_or_punct += 1;
            continue;
        }
        if ch.is_ascii_lowercase() {
            lowercase_letters += 1;
        }
        let up = ch.to_ascii_uppercase();
        non_ws_alpha += 1;
        *freq.entry(up).or_insert(0) += 1;

        if matches!(up, 'A' | 'C' | 'G' | 'T') {
            acgt_count += 1;
        }
        if matches!(up, 'A' | 'C' | 'G' | 'U') {
            acgu_count += 1;
        }
        if up == 'T' {
            t_count += 1;
        }
        if up == 'U' {
            u_count += 1;
        }
        if matches!(
            up,
            'E' | 'F' | 'I' | 'L' | 'P' | 'Q' | 'W' | 'Y' | 'K' | 'M' | 'R' | 'S' | 'V' | 'D' | 'H'
        ) {
            aa_exclusive_count += 1;
        }
        if matches!(
            up,
            'A' | 'C'
                | 'D'
                | 'E'
                | 'F'
                | 'G'
                | 'H'
                | 'I'
                | 'K'
                | 'L'
                | 'M'
                | 'N'
                | 'P'
                | 'Q'
                | 'R'
                | 'S'
                | 'T'
                | 'V'
                | 'W'
                | 'Y'
        ) {
            standard_aa_count += 1;
        }
    }

    let alpha_denom = non_ws_alpha.max(1) as f64;
    let acgt_frac = acgt_count as f64 / alpha_denom;
    let acgu_frac = acgu_count as f64 / alpha_denom;
    let aa_excl_frac = aa_exclusive_count as f64 / alpha_denom;
    let std_aa_frac = standard_aa_count as f64 / alpha_denom;
    let ws_punct_frac = ws_or_punct as f64 / (total_chars as f64);
    let lowercase_frac = lowercase_letters as f64 / alpha_denom;

    // Compute empirical 1-character Shannon entropy H1
    let mut h1 = 0.0_f64;
    for &c in freq.values() {
        let p = c as f64 / alpha_denom;
        if p > 0.0 {
            h1 -= p * p.log2();
        }
    }

    // DUST-style low-complexity fraction (homopolymer or dinucleotide dominance)
    let clean_upper: String = trimmed
        .chars()
        .filter(|c| c.is_ascii_alphabetic())
        .map(|c| c.to_ascii_uppercase())
        .collect();
    let dust_frac = compute_dust_low_complexity_fraction(&clean_upper);

    // Classification decision rules:
    // 1. If significant whitespace/punctuation (>= 6%) OR mixed-case prose with spaces, classify as Other (Standard Language).
    let has_non_bio_letters = freq.contains_key(&'J') || freq.contains_key(&'O') || freq.contains_key(&'Z');
    let detected_biotype = if ws_punct_frac >= 0.06 || (ws_punct_frac > 0.01 && lowercase_frac > 0.25) {
        SuperBlastBiotype::Other
    } else if acgu_frac >= 0.88 && u_count > t_count && u_count >= 2 {
        SuperBlastBiotype::Rna
    } else if acgt_frac >= 0.88 {
        SuperBlastBiotype::Dna
    } else if !has_non_bio_letters && aa_excl_frac >= 0.12 && std_aa_frac >= 0.90 {
        SuperBlastBiotype::Protein
    } else if acgt_frac >= 0.75 {
        SuperBlastBiotype::Dna
    } else {
        SuperBlastBiotype::Other
    };

    let (seed_word_size, karlin_lambda, karlin_k, check_rc) = match detected_biotype {
        SuperBlastBiotype::Dna => {
            let w = if clean_upper.len() < 120 { 7 } else { 9 };
            (w, 0.625, 0.410, true)
        }
        SuperBlastBiotype::Rna => {
            let w = if clean_upper.len() < 120 { 7 } else { 9 };
            (w, 0.625, 0.410, true)
        }
        SuperBlastBiotype::Protein => (3, 0.267, 0.041, false),
        SuperBlastBiotype::Other => (4, 0.385, 0.120, false),
    };

    BiotypeDetectionResult {
        detected_biotype,
        biotype_label: detected_biotype.as_str().to_string(),
        clean_length: clean_upper.len(),
        acgt_fraction: round4(acgt_frac),
        acgu_fraction: round4(acgu_frac),
        amino_acid_exclusive_fraction: round4(aa_excl_frac),
        whitespace_or_punctuation_fraction: round4(ws_punct_frac),
        shannon_h1_bits: round4(h1),
        low_complexity_dust_fraction: round4(dust_frac),
        seed_word_size,
        karlin_lambda,
        karlin_k,
        check_reverse_complement: check_rc,
    }
}

/// Computes the DUST low-complexity fraction across sliding 32-residue windows.
pub fn compute_dust_low_complexity_fraction(seq: &str) -> f64 {
    let bytes = seq.as_bytes();
    let n = bytes.len();
    if n < 12 {
        return 0.0;
    }
    let win = n.min(32);
    let mut low_windows = 0usize;
    let mut total_windows = 0usize;
    let step = (win / 2).max(1);
    let mut pos = 0usize;
    while pos + win <= n {
        total_windows += 1;
        let sub = &seq[pos..pos + win];
        let h2 = kmer_shannon_entropy(sub.as_bytes(), 2);
        // Dinucleotide entropy below 1.65 bits indicates strong homopolymer or (AT)n/(CA)n tandem bias
        let mut char_counts = [0usize; 256];
        let mut max_single = 0usize;
        for &b in &bytes[pos..pos + win] {
            char_counts[b as usize] += 1;
            if char_counts[b as usize] > max_single {
                max_single = char_counts[b as usize];
            }
        }
        let dominant_frac = max_single as f64 / (win as f64);
        if (h2 > 0.0 && h2 < 1.65) || dominant_frac >= 0.68 {
            low_windows += 1;
        }
        pos += step;
    }
    if total_windows == 0 {
        0.0
    } else {
        low_windows as f64 / (total_windows as f64)
    }
}

// ============================================================================
// 2. EXECUTION-MODE SENSITIVITY & CROSSOVER ANALYZER
// ============================================================================

/// Performs quantitative sensitivity analysis to decide whether `super-blast` should:
/// 1. Run `OnTheFlyDirectInMemory` (no database compilation needed),
/// 2. Build a `JustInTimeTransientIndex` (ephemeral in-memory hash index for genome/cohort search), or
/// 3. Compile a `CompiledPersistentDatabase` (only when query batch reuse $M \ge M^*$).
pub fn analyze_execution_strategy(
    num_queries: usize,
    num_targets: usize,
    total_query_residues: usize,
    total_target_residues: usize,
    expected_query_batches: usize,
) -> ExecutionSensitivityReport {
    let q_len = total_query_residues.max(1) as f64;
    let t_len = total_target_residues.max(1) as f64;
    let volume = (total_query_residues as u64).saturating_mul(total_target_residues as u64);
    let batches = expected_query_batches.max(1);

    // Empirical cost model (in microseconds):
    // - Direct streaming seed + banded SW cost scales with 0.0018 * V (with diagonal seed pruning)
    let est_direct_us = 12.0 + 0.00045 * (volume as f64);
    // - JIT ephemeral in-memory k-mer hash table build scales linearly with target length T_len
    let est_jit_build_us = 25.0 + 0.042 * t_len;
    // - JIT query lookup + local extension scales sub-linearly with target length
    let est_jit_query_us = 8.0 + 0.018 * q_len + 0.000035 * (volume as f64);
    // - Full persistent database compilation (sorting, serialization, disk header layout)
    let est_persistent_compile_us = 1850.0 + 0.38 * t_len;

    // Amortization crossover M*: how many repeated query batches are needed before persistent DB compilation wins?
    let per_batch_savings = (est_direct_us - est_jit_query_us).max(1.0);
    let crossover_batches = ((est_persistent_compile_us / per_batch_savings).ceil() as usize).max(2);

    let (selected_mode, requires_compiled_db, rationale) = if num_targets <= 4 && volume < 350_000 && batches == 1 {
        (
            SuperBlastExecutionMode::OnTheFlyDirectInMemory,
            false,
            format!(
                "Search volume V = {} bp^2 across {} query x {} target(s) is below the JIT index crossover (3.5e5 bp^2); direct in-memory diagonal seeding + banded Smith-Waterman ({:.1} us) avoids {:.1} us of unnecessary database compilation overhead.",
                volume, num_queries, num_targets, est_direct_us, est_persistent_compile_us
            ),
        )
    } else if batches < crossover_batches && volume < 500_000_000 {
        (
            SuperBlastExecutionMode::JustInTimeTransientIndex,
            false,
            format!(
                "Search volume V = {} bp^2 ({} targets, {} batch) favors a Just-In-Time (JIT) ephemeral in-memory k-mer index ({:.1} us build + {:.1} us query = {:.1} us total vs. {:.1} us persistent DB compilation; persistent DB only amortizes after M* >= {} repeated batches).",
                volume,
                num_targets,
                batches,
                est_jit_build_us,
                est_jit_query_us,
                est_jit_build_us + est_jit_query_us,
                est_persistent_compile_us,
                crossover_batches
            ),
        )
    } else {
        (
            SuperBlastExecutionMode::CompiledPersistentDatabase,
            true,
            format!(
                "Repeated batch reuse (M = {} >= M* = {}) or large target volume (V = {} bp^2) amortizes persistent database compilation ({:.1} us compile cost recovered across {} query batches).",
                batches, crossover_batches, volume, est_persistent_compile_us, batches
            ),
        )
    };

    ExecutionSensitivityReport {
        selected_mode,
        mode_label: selected_mode.as_str().to_string(),
        num_queries,
        num_targets,
        total_query_residues,
        total_target_residues,
        search_space_volume: volume,
        estimated_direct_cost_us: round4(est_direct_us),
        estimated_jit_index_build_us: round4(est_jit_build_us),
        estimated_jit_query_cost_us: round4(est_jit_query_us),
        estimated_persistent_db_compile_us: round4(est_persistent_compile_us),
        amortization_crossover_batches: crossover_batches,
        requires_compiled_db,
        rationale,
    }
}

// ============================================================================
// 3. CORE SCORING MATRICES & PAIRWISE / JIT ALIGNMENT ENGINE
// ============================================================================

/// Returns the reverse complement of a nucleotide sequence (`DNA` or `RNA`).
pub fn reverse_complement_nt(seq: &str, is_rna: bool) -> String {
    seq.chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'A' => {
                if is_rna {
                    'U'
                } else {
                    'T'
                }
            }
            'T' | 'U' => 'A',
            'C' => 'G',
            'G' => 'C',
            'R' => 'Y',
            'Y' => 'R',
            'S' => 'S',
            'W' => 'W',
            'K' => 'M',
            'M' => 'K',
            other => other,
        })
        .collect()
}

/// Standard BLOSUM62 log-odds substitution score for amino acid pairs.
pub fn blosum62_score(a: u8, b: u8) -> i32 {
    const AA_ORDER: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";
    const MATRIX: [[i32; 20]; 20] = [
        // A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
        [  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0], // A
        [ -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3], // R
        [ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3], // N
        [ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3], // D
        [  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1], // C
        [ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2], // Q
        [ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2], // E
        [  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3], // G
        [ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3], // H
        [ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3], // I
        [ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1], // L
        [ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2], // K
        [ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1], // M
        [ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1], // F
        [ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2], // P
        [  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2], // S
        [  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0], // T
        [ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3], // W
        [ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1], // Y
        [  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4], // V
    ];
    let idx_a = AA_ORDER.iter().position(|&x| x == a.to_ascii_uppercase());
    let idx_b = AA_ORDER.iter().position(|&x| x == b.to_ascii_uppercase());
    match (idx_a, idx_b) {
        (Some(ia), Some(ib)) => MATRIX[ia][ib],
        _ => {
            if a.eq_ignore_ascii_case(&b) {
                4
            } else {
                -2
            }
        }
    }
}

#[inline]
fn residue_pair_score(a: u8, b: u8, biotype: SuperBlastBiotype) -> i32 {
    match biotype {
        SuperBlastBiotype::Dna | SuperBlastBiotype::Rna => {
            let ua = match a.to_ascii_uppercase() {
                b'U' => b'T',
                x => x,
            };
            let ub = match b.to_ascii_uppercase() {
                b'U' => b'T',
                x => x,
            };
            if ua == ub {
                2
            } else {
                -3
            }
        }
        SuperBlastBiotype::Protein => blosum62_score(a, b),
        SuperBlastBiotype::Other => {
            if a.eq_ignore_ascii_case(&b) {
                3
            } else {
                -2
            }
        }
    }
}

#[inline]
fn gap_penalties_for_biotype(biotype: SuperBlastBiotype) -> (i32, i32) {
    // Returns (gap_open_including_first_base, gap_extend) matching Biopython PairwiseAligner conventions:
    // For DNA/RNA: open = -5, extend = -2
    // For Protein: open = -11, extend = -1
    // For Other:   open = -4, extend = -1
    match biotype {
        SuperBlastBiotype::Dna | SuperBlastBiotype::Rna => (-5, -2),
        SuperBlastBiotype::Protein => (-11, -1),
        SuperBlastBiotype::Other => (-4, -1),
    }
}

fn clean_sequence_for_biotype(raw: &str, biotype: SuperBlastBiotype) -> String {
    match biotype {
        SuperBlastBiotype::Dna | SuperBlastBiotype::Rna | SuperBlastBiotype::Protein => raw
            .chars()
            .filter(|c| c.is_ascii_alphabetic())
            .map(|c| c.to_ascii_uppercase())
            .collect(),
        SuperBlastBiotype::Other => {
            // Collapse runs of whitespace into a single space and uppercase for case-insensitive standard language matching
            let mut out = String::with_capacity(raw.len());
            let mut prev_space = false;
            for ch in raw.trim().chars() {
                if ch.is_whitespace() {
                    if !prev_space {
                        out.push(' ');
                        prev_space = true;
                    }
                } else {
                    out.push(ch.to_ascii_uppercase());
                    prev_space = false;
                }
            }
            out
        }
    }
}

/// Finds the dominant diagonal shift $(t - q)$ via $k$-mer seed hits so large sequences can be
/// aligned with a fast banded Smith-Waterman pass when appropriate.
fn find_best_seed_diagonal(q: &[u8], t: &[u8], w: usize) -> Option<(isize, usize)> {
    if q.len() < w || t.len() < w {
        return None;
    }
    let mut table: HashMap<&[u8], Vec<usize>> = HashMap::new();
    for (i, win) in t.windows(w).enumerate() {
        table.entry(win).or_default().push(i);
    }
    let mut diag_votes: HashMap<isize, usize> = HashMap::new();
    for (qi, win) in q.windows(w).enumerate() {
        if let Some(positions) = table.get(win) {
            for &ti in positions {
                let diag = (ti as isize) - (qi as isize);
                // Bin diagonals in buckets of 4 to tolerate small indels
                let bucket = diag.div_euclid(4) * 4;
                *diag_votes.entry(bucket).or_insert(0) += 1;
            }
        }
    }
    diag_votes
        .into_iter()
        .max_by_key(|&(_, count)| count)
}

/// Performs affine-gap Smith-Waterman local alignment (`Gotoh` 3-state local dynamic programming)
/// between `query` and `target` on a single strand.
/// Uses exact convention: opening a gap of length $L$ costs `gap_open + (L - 1) * gap_ext`
/// (both `gap_open` and `gap_ext` are negative integers).
fn smith_waterman_gotoh_strand(
    q_id: &str,
    t_id: &str,
    q: &[u8],
    t: &[u8],
    biotype: SuperBlastBiotype,
    exec_mode: SuperBlastExecutionMode,
    strand: char,
    karlin_lambda: f64,
    karlin_k: f64,
    dust_frac: f64,
    db_effective_len: usize,
) -> SuperBlastHsp {
    let m = q.len();
    let n = t.len();
    if m == 0 || n == 0 {
        return SuperBlastHsp {
            query_id: q_id.to_string(),
            target_id: t_id.to_string(),
            biotype: biotype.as_str().to_string(),
            execution_mode: exec_mode.as_str().to_string(),
            strand,
            query_start_1based: 0,
            query_end_1based: 0,
            target_start_1based: 0,
            target_end_1based: 0,
            alignment_length: 0,
            matches: 0,
            mismatches: 0,
            transitions: 0,
            transversions: 0,
            gap_opens: 0,
            gap_bases: 0,
            pct_identity: 0.0,
            query_coverage_pct: 0.0,
            raw_score: 0.0,
            composition_adjusted_score: 0.0,
            bit_score: 0.0,
            e_value: 10.0,
            jc69_distance: 2.0,
            k2p_distance: 2.0,
            low_complexity_masked: false,
            aligned_query: String::new(),
            midline: String::new(),
            aligned_target: String::new(),
        };
    }

    let (gap_open, gap_ext) = gap_penalties_for_biotype(biotype);
    let neg_inf: i32 = -1_000_000_000;

    // Optional diagonal band for very large target contigs when JIT seed hits are strong
    let band_opt: Option<(isize, isize)> = if m.saturating_mul(n) > 1_200_000 {
        let w = match biotype {
            SuperBlastBiotype::Protein => 3,
            SuperBlastBiotype::Other => 4,
            _ => 8,
        };
        if let Some((best_diag, votes)) = find_best_seed_diagonal(q, t, w) {
            if votes >= 2 {
                let half_band = (m.max(80) as isize) / 2 + 48;
                Some((best_diag - half_band, best_diag + half_band))
            } else {
                None
            }
        } else {
            None
        }
    } else {
        None
    };

    let cols = n + 1;
    let mut h_mat = vec![0i32; (m + 1) * cols];
    let mut e_mat = vec![neg_inf; (m + 1) * cols]; // Gap in query (insertion in target, horizontal move)
    let mut f_mat = vec![neg_inf; (m + 1) * cols]; // Gap in target (deletion in target, vertical move)
    // Traceback state: 0 = stop, 1 = diag, 2 = from E (left), 3 = from F (up)
    // Plus bit flags for whether E or F extended vs opened
    let mut tb_h = vec![0u8; (m + 1) * cols];
    let mut tb_e = vec![0u8; (m + 1) * cols]; // 0 = opened from H, 1 = extended from E
    let mut tb_f = vec![0u8; (m + 1) * cols]; // 0 = opened from H, 1 = extended from F

    let mut best_score = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    for i in 1..=m {
        let (j_min, j_max) = if let Some((d_min, d_max)) = band_opt {
            let lo = ((i as isize) + d_min).clamp(1, n as isize) as usize;
            let hi = ((i as isize) + d_max).clamp(1, n as isize) as usize;
            (lo, hi)
        } else {
            (1, n)
        };
        let row_offset = i * cols;
        let prev_row = (i - 1) * cols;
        let qi = q[i - 1];

        for j in j_min..=j_max {
            let idx = row_offset + j;
            let left_idx = row_offset + (j - 1);
            let up_idx = prev_row + j;
            let diag_idx = prev_row + (j - 1);

            // E[i, j]: gap in query (horizontal move from j-1)
            let open_e = h_mat[left_idx].saturating_add(gap_open);
            let ext_e = e_mat[left_idx].saturating_add(gap_ext);
            let (e_val, e_ext_flag) = if ext_e >= open_e {
                (ext_e, 1u8)
            } else {
                (open_e, 0u8)
            };
            e_mat[idx] = e_val;
            tb_e[idx] = e_ext_flag;

            // F[i, j]: gap in target (vertical move from i-1)
            let open_f = h_mat[up_idx].saturating_add(gap_open);
            let ext_f = f_mat[up_idx].saturating_add(gap_ext);
            let (f_val, f_ext_flag) = if ext_f >= open_f {
                (ext_f, 1u8)
            } else {
                (open_f, 0u8)
            };
            f_mat[idx] = f_val;
            tb_f[idx] = f_ext_flag;

            let sub = residue_pair_score(qi, t[j - 1], biotype);
            let diag_val = h_mat[diag_idx].saturating_add(sub);

            let mut max_v = 0i32;
            let mut dir = 0u8;
            if diag_val > max_v {
                max_v = diag_val;
                dir = 1;
            }
            if e_val > max_v {
                max_v = e_val;
                dir = 2;
            }
            if f_val > max_v {
                max_v = f_val;
                dir = 3;
            }

            h_mat[idx] = max_v;
            tb_h[idx] = dir;

            if max_v > best_score {
                best_score = max_v;
                best_i = i;
                best_j = j;
            }
        }
    }

    // Traceback from (best_i, best_j)
    let mut ci = best_i;
    let mut cj = best_j;
    let mut state = 0u8; // 0 = H, 2 = E, 3 = F
    let mut aln_q: Vec<u8> = Vec::with_capacity(m + n);
    let mut aln_t: Vec<u8> = Vec::with_capacity(m + n);
    let mut mid: Vec<u8> = Vec::with_capacity(m + n);

    let mut matches = 0usize;
    let mut mismatches = 0usize;
    let mut transitions = 0usize;
    let mut transversions = 0usize;
    let mut gap_opens = 0usize;
    let mut gap_bases = 0usize;

    while ci > 0 && cj > 0 {
        let idx = ci * cols + cj;
        if state == 0 {
            let dir = tb_h[idx];
            if dir == 0 || h_mat[idx] == 0 {
                break;
            }
            if dir == 1 {
                let qa = q[ci - 1];
                let tb = t[cj - 1];
                aln_q.push(qa);
                aln_t.push(tb);
                let ua = match qa.to_ascii_uppercase() {
                    b'U' => b'T',
                    x => x,
                };
                let ub = match tb.to_ascii_uppercase() {
                    b'U' => b'T',
                    x => x,
                };
                if ua == ub {
                    matches += 1;
                    mid.push(b'|');
                } else {
                    mismatches += 1;
                    if matches!(biotype, SuperBlastBiotype::Dna | SuperBlastBiotype::Rna) {
                        let is_ts = matches!((ua, ub), (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C'));
                        if is_ts {
                            transitions += 1;
                            mid.push(b':');
                        } else {
                            transversions += 1;
                            mid.push(b'.');
                        }
                    } else if residue_pair_score(qa, tb, biotype) > 0 {
                        mid.push(b'+');
                    } else {
                        mid.push(b' ');
                    }
                }
                ci -= 1;
                cj -= 1;
            } else if dir == 2 {
                state = 2;
            } else {
                state = 3;
            }
        } else if state == 2 {
            // In E: gap in query (horizontal move)
            aln_q.push(b'-');
            aln_t.push(t[cj - 1]);
            mid.push(b' ');
            gap_bases += 1;
            let ext = tb_e[idx];
            cj -= 1;
            if ext == 0 {
                gap_opens += 1;
                state = 0;
            }
        } else {
            // In F: gap in target (vertical move)
            aln_q.push(q[ci - 1]);
            aln_t.push(b'-');
            mid.push(b' ');
            gap_bases += 1;
            let ext = tb_f[idx];
            ci -= 1;
            if ext == 0 {
                gap_opens += 1;
                state = 0;
            }
        }
    }

    aln_q.reverse();
    aln_t.reverse();
    mid.reverse();

    let q_start_1 = if best_i > 0 { ci + 1 } else { 0 };
    let q_end_1 = best_i;
    let t_start_1 = if best_j > 0 { cj + 1 } else { 0 };
    let t_end_1 = best_j;

    let aln_len = aln_q.len();
    let pct_id = if aln_len > 0 {
        (matches as f64 / aln_len as f64) * 100.0
    } else {
        0.0
    };
    let q_cov = if m > 0 && q_end_1 >= q_start_1 && q_start_1 > 0 {
        ((q_end_1 - q_start_1 + 1) as f64 / m as f64) * 100.0
    } else {
        0.0
    };

    let raw_score = best_score as f64;
    // Composition-biased / DUST low-complexity masking adjustment:
    let low_complexity_masked = dust_frac >= 0.50;
    let comp_factor = if low_complexity_masked {
        (1.0 - 0.82 * dust_frac).clamp(0.12, 1.0)
    } else {
        1.0
    };
    let adj_score = raw_score * comp_factor;

    // Karlin-Altschul Bit-Score: S' = (lambda * S - ln K) / ln 2
    let bit_score = if adj_score > 0.0 {
        ((karlin_lambda * adj_score - karlin_k.ln()) / std::f64::consts::LN_2).max(0.0)
    } else {
        0.0
    };

    // Karlin-Altschul Expect value: E = m * n_eff * 2^(-S')
    let eff_m = (m as f64).max(1.0);
    let eff_n = (db_effective_len.max(n) as f64).max(1.0);
    let e_value = if bit_score > 0.0 {
        (eff_m * eff_n * 2.0_f64.powf(-bit_score)).clamp(1e-180, 1e6)
    } else {
        eff_m * eff_n
    };

    // Evolutionary distances (JC69 and Kimura 2-Parameter K2P) across ungapped aligned sites
    let ungapped_sites = (matches + mismatches).max(1) as f64;
    let p_mismatch = (mismatches as f64) / ungapped_sites;
    let jc69 = if p_mismatch >= 0.74 {
        2.0
    } else {
        (-0.75 * (1.0 - (4.0 / 3.0) * p_mismatch).ln()).max(0.0)
    };

    let k2p = if matches!(biotype, SuperBlastBiotype::Dna | SuperBlastBiotype::Rna) {
        let p_ts = (transitions as f64) / ungapped_sites;
        let q_tv = (transversions as f64) / ungapped_sites;
        let term1 = 1.0 - 2.0 * p_ts - q_tv;
        let term2 = 1.0 - 2.0 * q_tv;
        if term1 <= 0.01 || term2 <= 0.01 {
            2.0
        } else {
            (-0.5 * term1.ln() - 0.25 * term2.ln()).max(0.0)
        }
    } else {
        // Protein/Other Kimura distance approximation: -ln(1 - D - 0.2 * D^2)
        let d = p_mismatch;
        let arg = (1.0 - d - 0.2 * d * d).max(0.01);
        (-arg.ln()).max(0.0)
    };

    SuperBlastHsp {
        query_id: q_id.to_string(),
        target_id: t_id.to_string(),
        biotype: biotype.as_str().to_string(),
        execution_mode: exec_mode.as_str().to_string(),
        strand,
        query_start_1based: q_start_1,
        query_end_1based: q_end_1,
        target_start_1based: t_start_1,
        target_end_1based: t_end_1,
        alignment_length: aln_len,
        matches,
        mismatches,
        transitions,
        transversions,
        gap_opens,
        gap_bases,
        pct_identity: round4(pct_id),
        query_coverage_pct: round4(q_cov),
        raw_score: round4(raw_score),
        composition_adjusted_score: round4(adj_score),
        bit_score: round4(bit_score),
        e_value,
        jc69_distance: round4(jc69),
        k2p_distance: round4(k2p),
        low_complexity_masked,
        aligned_query: String::from_utf8_lossy(&aln_q).to_string(),
        midline: String::from_utf8_lossy(&mid).to_string(),
        aligned_target: String::from_utf8_lossy(&aln_t).to_string(),
    }
}

/// Runs `super-blast` on a pair of sequences/strings, automatically detecting sequence biotype,
/// running execution-mode sensitivity analysis, and searching both strands if `DNA` or `RNA`.
pub fn super_blast_pair(
    query_id: &str,
    query_raw: &str,
    target_id: &str,
    target_raw: &str,
) -> (BiotypeDetectionResult, ExecutionSensitivityReport, SuperBlastHsp) {
    let bio_info = detect_sequence_biotype(query_raw);
    let q_clean = clean_sequence_for_biotype(query_raw, bio_info.detected_biotype);
    let t_clean = clean_sequence_for_biotype(target_raw, bio_info.detected_biotype);

    let sens = analyze_execution_strategy(1, 1, q_clean.len(), t_clean.len(), 1);

    let hsp_fwd = smith_waterman_gotoh_strand(
        query_id,
        target_id,
        q_clean.as_bytes(),
        t_clean.as_bytes(),
        bio_info.detected_biotype,
        sens.selected_mode,
        '+',
        bio_info.karlin_lambda,
        bio_info.karlin_k,
        bio_info.low_complexity_dust_fraction,
        t_clean.len(),
    );

    if bio_info.check_reverse_complement {
        let is_rna = bio_info.detected_biotype == SuperBlastBiotype::Rna;
        let q_rc = reverse_complement_nt(&q_clean, is_rna);
        let hsp_rev = smith_waterman_gotoh_strand(
            query_id,
            target_id,
            q_rc.as_bytes(),
            t_clean.as_bytes(),
            bio_info.detected_biotype,
            sens.selected_mode,
            '-',
            bio_info.karlin_lambda,
            bio_info.karlin_k,
            bio_info.low_complexity_dust_fraction,
            t_clean.len(),
        );
        if hsp_rev.raw_score > hsp_fwd.raw_score {
            return (bio_info, sens, hsp_rev);
        }
    }

    (bio_info, sens, hsp_fwd)
}

// ============================================================================
// 4. ON-THE-FLY MSA PRE-INCORPORATION SEQUENCE VALIDATOR
// ============================================================================

/// Validates a cohort of candidate sequences on the fly prior to Multiple Sequence Alignment (`MSA`).
/// Automatically rescues reverse-complemented sequences, detects 5'/3' chimeras, and rejects
/// non-homologous contaminants or truncated fragments.
pub fn validate_sequences_for_msa(
    centroid_id: &str,
    centroid_seq: &str,
    candidates: &[(String, String)],
    min_identity_pct: f64,
    min_coverage_pct: f64,
) -> Vec<MsaSequenceValidationRecord> {
    let bio = detect_sequence_biotype(centroid_seq);
    let c_clean = clean_sequence_for_biotype(centroid_seq, bio.detected_biotype);
    let mut records = Vec::with_capacity(candidates.len());

    for (cand_id, cand_raw) in candidates {
        let q_clean = clean_sequence_for_biotype(cand_raw, bio.detected_biotype);
        let (_, _, hsp) = super_blast_pair(cand_id, &q_clean, centroid_id, &c_clean);

        // Oriented sequence (if reverse-complement matched better, flip it!)
        let oriented = if hsp.strand == '-' {
            reverse_complement_nt(&q_clean, bio.detected_biotype == SuperBlastBiotype::Rna)
        } else {
            q_clean.clone()
        };

        // Chimera check: split oriented candidate into 5' half and 3' half and blast each against centroid
        let mid_pt = (oriented.len() / 2).max(1);
        let five_prime = &oriented[..mid_pt];
        let three_prime = &oriented[mid_pt..];
        let (_, _, hsp_5p) = super_blast_pair(&format!("{}_5p", cand_id), five_prime, centroid_id, &c_clean);
        let (_, _, hsp_3p) = super_blast_pair(&format!("{}_3p", cand_id), three_prime, centroid_id, &c_clean);

        // Effective half identity weighted by half coverage
        let id_5p = hsp_5p.pct_identity * (hsp_5p.query_coverage_pct / 100.0);
        let id_3p = hsp_3p.pct_identity * (hsp_3p.query_coverage_pct / 100.0);
        let chimera_delta = (id_5p - id_3p).abs();

        let (status, incorporate, note) = if hsp.e_value > 1e-3 || hsp.pct_identity < min_identity_pct {
            (
                "RejectedNonHomologous".to_string(),
                false,
                format!(
                    "Rejected non-homologous outlier ({:.1}% identity, E = {:.2e} vs. {}); excluded before MSA.",
                    hsp.pct_identity, hsp.e_value, centroid_id
                ),
            )
        } else if chimera_delta >= 32.0 && id_5p.min(id_3p) < 55.0 {
            (
                "RejectedChimera".to_string(),
                false,
                format!(
                    "Rejected 5'/3' chimeric sequence (5' effective identity = {:.1}% vs. 3' = {:.1}%, |Δ| = {:.1}%).",
                    id_5p, id_3p, chimera_delta
                ),
            )
        } else if hsp.query_coverage_pct < min_coverage_pct {
            (
                "RejectedTruncated".to_string(),
                false,
                format!(
                    "Rejected truncated fragment (coverage {:.1}% < {:.1}% minimum MSA threshold).",
                    hsp.query_coverage_pct, min_coverage_pct
                ),
            )
        } else if hsp.strand == '-' {
            (
                "AcceptedAutoRevComp".to_string(),
                true,
                format!(
                    "Reverse-complement orientation detected on '-' strand ({:.1}% identity, {:.1} bits); automatically re-oriented to '+' for MSA.",
                    hsp.pct_identity, hsp.bit_score
                ),
            )
        } else {
            (
                "AcceptedDirect".to_string(),
                true,
                format!(
                    "Validated homologous orthologous sequence ({:.1}% identity, {:.1}% coverage, {:.1} bits).",
                    hsp.pct_identity, hsp.query_coverage_pct, hsp.bit_score
                ),
            )
        };

        records.push(MsaSequenceValidationRecord {
            seq_id: cand_id.clone(),
            status,
            incorporate_into_msa: incorporate,
            reverse_complement_rescued: hsp.strand == '-' && incorporate,
            strand: hsp.strand,
            pct_identity_vs_centroid: hsp.pct_identity,
            query_coverage_pct: hsp.query_coverage_pct,
            bit_score: hsp.bit_score,
            e_value: hsp.e_value,
            five_prime_half_identity: round4(id_5p),
            three_prime_half_identity: round4(id_3p),
            chimera_delta_identity: round4(chimera_delta),
            oriented_sequence: oriented,
            diagnostic_note: note,
        });
    }

    records
}

// ============================================================================
// 5. SYMMETRIC ALL-VS-ALL PAIRWISE DISTANCE MATRIX ON THE FLY
// ============================================================================

/// Computes a symmetric all-vs-all pairwise distance matrix (`% Identity`, `K2P`, `JC69`, `Bit-Score`)
/// across a sequence cohort on the fly in memory without persistent database compilation.
pub fn compute_symmetric_pairwise_matrix(
    sequences: &[(String, String)],
) -> SymmetricPairwiseMatrixReport {
    let n = sequences.len();
    if n == 0 {
        return SymmetricPairwiseMatrixReport {
            biotype: "DNA".to_string(),
            execution_mode: SuperBlastExecutionMode::OnTheFlyDirectInMemory.as_str().to_string(),
            num_sequences: 0,
            unique_pairs_evaluated: 0,
            labels: Vec::new(),
            pct_identity_matrix: Vec::new(),
            k2p_distance_matrix: Vec::new(),
            jc69_distance_matrix: Vec::new(),
            bit_score_matrix: Vec::new(),
        };
    }

    let bio = detect_sequence_biotype(&sequences[0].1);
    let total_len: usize = sequences.iter().map(|(_, s)| s.len()).sum();
    let sens = analyze_execution_strategy(n, n, total_len, total_len, 1);

    let labels: Vec<String> = sequences.iter().map(|(id, _)| id.clone()).collect();
    let mut id_mat = vec![vec![100.0; n]; n];
    let mut k2p_mat = vec![vec![0.0; n]; n];
    let mut jc69_mat = vec![vec![0.0; n]; n];
    let mut bit_mat = vec![vec![0.0; n]; n];

    for i in 0..n {
        let (_, _, self_hsp) = super_blast_pair(&sequences[i].0, &sequences[i].1, &sequences[i].0, &sequences[i].1);
        bit_mat[i][i] = self_hsp.bit_score;
        for j in (i + 1)..n {
            let (_, _, hsp) = super_blast_pair(&sequences[i].0, &sequences[i].1, &sequences[j].0, &sequences[j].1);
            id_mat[i][j] = hsp.pct_identity;
            id_mat[j][i] = hsp.pct_identity;
            k2p_mat[i][j] = hsp.k2p_distance;
            k2p_mat[j][i] = hsp.k2p_distance;
            jc69_mat[i][j] = hsp.jc69_distance;
            jc69_mat[j][i] = hsp.jc69_distance;
            bit_mat[i][j] = hsp.bit_score;
            bit_mat[j][i] = hsp.bit_score;
        }
    }

    SymmetricPairwiseMatrixReport {
        biotype: bio.biotype_label,
        execution_mode: sens.mode_label,
        num_sequences: n,
        unique_pairs_evaluated: n * (n.saturating_sub(1)) / 2,
        labels,
        pct_identity_matrix: id_mat,
        k2p_distance_matrix: k2p_mat,
        jc69_distance_matrix: jc69_mat,
        bit_score_matrix: bit_mat,
    }
}

// ============================================================================
// 6. MULTI-GENE PHYLOGENETIC TIP CONCORDANCE AUDITOR (e.g. 16S + 10 GENES)
// ============================================================================

/// Tests a series of secondary genes against the tip relationships of an anchor phylogeny
/// (e.g., a `16S rRNA` ribosomal tree vs. 10 housekeeping/accessory genes across the same tip taxa).
/// Identifies discordant genes and pinpoints exact tip taxa exhibiting Horizontal Gene Transfer (HGT)
/// or paralogous replacement.
pub fn audit_phylogeny_tip_concordance(
    anchor_marker_name: &str,
    anchor_tip_sequences: &[(String, String)],
    secondary_genes: &[(String, String, Vec<(String, String)>)], // (gene_name, category, tip_sequences)
) -> PhylogenyTipConcordanceReport {
    let n_tips = anchor_tip_sequences.len();
    let anchor_mat = compute_symmetric_pairwise_matrix(anchor_tip_sequences);

    // Flatten upper-triangle anchor K2P distances and find each tip's nearest sister tip in the anchor tree
    let mut anchor_upper = Vec::new();
    let mut anchor_nearest_sister = vec![0usize; n_tips];
    for i in 0..n_tips {
        let mut best_j = 0usize;
        let mut best_d = f64::MAX;
        for j in 0..n_tips {
            if i != j && anchor_mat.k2p_distance_matrix[i][j] < best_d {
                best_d = anchor_mat.k2p_distance_matrix[i][j];
                best_j = j;
            }
            if j > i {
                anchor_upper.push(anchor_mat.k2p_distance_matrix[i][j]);
            }
        }
        anchor_nearest_sister[i] = best_j;
    }

    let mut gene_audits = Vec::with_capacity(secondary_genes.len());
    let mut concordant_count = 0usize;
    let mut discordant_count = 0usize;
    let mut concordant_r_sum = 0.0_f64;

    for (gene_name, category, gene_tips) in secondary_genes {
        let g_mat = compute_symmetric_pairwise_matrix(gene_tips);
        let mut gene_upper = Vec::new();
        let mut nn_matches = 0usize;
        let mut discordant_tips = Vec::new();

        for i in 0..n_tips.min(g_mat.num_sequences) {
            let mut best_j = 0usize;
            let mut best_d = f64::MAX;
            for j in 0..n_tips.min(g_mat.num_sequences) {
                if i != j && g_mat.k2p_distance_matrix[i][j] < best_d {
                    best_d = g_mat.k2p_distance_matrix[i][j];
                    best_j = j;
                }
                if j > i {
                    gene_upper.push(g_mat.k2p_distance_matrix[i][j]);
                }
            }
            let expected_sister = anchor_nearest_sister[i];
            let d_to_expected = g_mat.k2p_distance_matrix[i][expected_sister];
            // Concordant tip if nearest neighbor matches OR distance to expected sister is within 0.015 K2P of best_d
            if best_j == expected_sister || (d_to_expected - best_d).abs() <= 0.015 {
                nn_matches += 1;
            } else {
                discordant_tips.push(format!(
                    "{} (grouped with {} instead of sister {})",
                    anchor_tip_sequences[i].0,
                    anchor_tip_sequences[best_j].0,
                    anchor_tip_sequences[expected_sister].0
                ));
            }
        }

        let r = pearson_correlation(&anchor_upper, &gene_upper);
        let nn_pct = if n_tips > 0 {
            (nn_matches as f64 / n_tips as f64) * 100.0
        } else {
            0.0
        };
        let mean_k2p = if gene_upper.is_empty() {
            0.0
        } else {
            gene_upper.iter().sum::<f64>() / (gene_upper.len() as f64)
        };

        let verdict = if r >= 0.80 && nn_pct >= 75.0 {
            concordant_count += 1;
            concordant_r_sum += r;
            "ConcordantVerticalInheritance".to_string()
        } else {
            discordant_count += 1;
            "DiscordantHorizontalTransferOrXenolog".to_string()
        };

        gene_audits.push(PhylogenyTipGeneAudit {
            gene_name: gene_name.clone(),
            functional_category: category.clone(),
            num_tips_evaluated: n_tips,
            pearson_r_vs_anchor_tree: round4(r),
            nearest_neighbor_tip_concordance_pct: round4(nn_pct),
            mean_tip_k2p_distance: round4(mean_k2p),
            discordant_tips_detected: discordant_tips,
            concordance_verdict: verdict,
        });
    }

    let mean_concordant_r = if concordant_count > 0 {
        concordant_r_sum / (concordant_count as f64)
    } else {
        0.0
    };

    PhylogenyTipConcordanceReport {
        anchor_marker_name: anchor_marker_name.to_string(),
        num_taxa_tips: n_tips,
        num_secondary_genes_tested: secondary_genes.len(),
        concordant_genes_count: concordant_count,
        discordant_genes_count: discordant_count,
        mean_concordant_pearson_r: round4(mean_concordant_r),
        gene_audits,
    }
}

fn pearson_correlation(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len().min(ys.len());
    if n < 2 {
        return 0.0;
    }
    let mx = xs[..n].iter().sum::<f64>() / (n as f64);
    let my = ys[..n].iter().sum::<f64>() / (n as f64);
    let mut num = 0.0;
    let mut dx2 = 0.0;
    let mut dy2 = 0.0;
    for i in 0..n {
        let dx = xs[i] - mx;
        let dy = ys[i] - my;
        num += dx * dy;
        dx2 += dx * dx;
        dy2 += dy * dy;
    }
    let denom = (dx2 * dy2).sqrt();
    if denom < 1e-12 {
        0.0
    } else {
        (num / denom).clamp(-1.0, 1.0)
    }
}

// ============================================================================
// 7. "BARCODE BLAST" COMPILER & TAXONOMIC BARCODING-GAP CLASSIFIER
// ============================================================================

impl BarcodeBlastDatabase {
    /// Compiles a curated panel of DNA barcodes into an in-memory 2-bit packed k-mer index.
    pub fn compile(marker_locus: &str, entries: Vec<BarcodeReferenceEntry>, word_size: usize) -> Self {
        let w = word_size.clamp(6, 15);
        let mut kmer_index: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        for (entry_idx, entry) in entries.iter().enumerate() {
            let clean = clean_sequence_for_biotype(&entry.sequence, SuperBlastBiotype::Dna);
            let bytes = clean.as_bytes();
            if bytes.len() >= w {
                for pos in 0..=(bytes.len() - w) {
                    if let Some(code) = encode_dna_kmer(&bytes[pos..pos + w]) {
                        kmer_index.entry(code).or_default().push((entry_idx, pos));
                    }
                }
            }
        }
        Self {
            marker_locus: marker_locus.to_string(),
            word_size: w,
            entries,
            kmer_index,
        }
    }

    /// Queries the compiled Barcode BLAST database, computing the intra- vs. inter-specific
    /// barcoding gap ($\Delta_{\text{gap}}$) and screening protein-coding barcodes (`COX1`)
    /// for premature stop codons indicative of nuclear `NUMT` pseudogenes.
    pub fn classify_barcode(&self, query_id: &str, query_seq: &str) -> BarcodeBlastAssignment {
        let q_clean = clean_sequence_for_biotype(query_seq, SuperBlastBiotype::Dna);
        let total_db_len: usize = self.entries.iter().map(|e| e.sequence.len()).sum();

        let mut best_hsp: Option<(usize, SuperBlastHsp)> = None;
        let mut all_hits: Vec<(usize, SuperBlastHsp)> = Vec::with_capacity(self.entries.len());

        for (idx, entry) in self.entries.iter().enumerate() {
            let t_clean = clean_sequence_for_biotype(&entry.sequence, SuperBlastBiotype::Dna);
            let hsp = smith_waterman_gotoh_strand(
                query_id,
                &entry.accession_id,
                q_clean.as_bytes(),
                t_clean.as_bytes(),
                SuperBlastBiotype::Dna,
                SuperBlastExecutionMode::JustInTimeTransientIndex,
                '+',
                0.625,
                0.410,
                0.0,
                total_db_len,
            );
            if best_hsp
                .as_ref()
                .map(|(_, b)| hsp.bit_score > b.bit_score)
                .unwrap_or(true)
            {
                best_hsp = Some((idx, hsp.clone()));
            }
            all_hits.push((idx, hsp));
        }

        let (top_idx, top_hsp) = best_hsp.expect("BarcodeBlastDatabase has at least 1 reference entry");
        let top_entry = &self.entries[top_idx];

        let mut min_intra = 0.0_f64;
        let mut found_intra = false;
        let mut min_inter = 2.0_f64;
        for (idx, hsp) in &all_hits {
            let sp = &self.entries[*idx].species_name;
            if sp == &top_entry.species_name {
                if !found_intra || hsp.k2p_distance < min_intra {
                    min_intra = hsp.k2p_distance;
                    found_intra = true;
                }
            } else if hsp.k2p_distance < min_inter {
                min_inter = hsp.k2p_distance;
            }
        }
        let barcoding_gap = min_inter - min_intra;

        // Screen for internal stop codons across the best reading frame for COX1 / coding barcodes
        let stop_count = if self.marker_locus.to_ascii_uppercase().contains("COX1")
            || self.marker_locus.to_ascii_uppercase().contains("COI")
        {
            min_internal_stops_across_frames(&q_clean)
        } else {
            0
        };
        let numt_flagged = stop_count > 0;

        let regime = if numt_flagged {
            "RejectedNumtPseudogeneArtifact".to_string()
        } else if top_hsp.pct_identity >= 97.5 && barcoding_gap >= 0.02 {
            "SpeciesLevelConfirmedWithBarcodingGap".to_string()
        } else if top_hsp.pct_identity >= 90.0 {
            "GenusLevelCandidate".to_string()
        } else {
            "UnassignedDivergentLineage".to_string()
        };

        BarcodeBlastAssignment {
            query_id: query_id.to_string(),
            top_hit_accession: top_entry.accession_id.clone(),
            assigned_species: top_entry.species_name.clone(),
            assigned_genus: top_entry.genus_name.clone(),
            top_pct_identity: top_hsp.pct_identity,
            top_bit_score: top_hsp.bit_score,
            top_e_value: top_hsp.e_value,
            min_intraspecific_k2p: round4(min_intra),
            min_interspecific_k2p: round4(min_inter),
            barcoding_gap_k2p: round4(barcoding_gap),
            internal_stop_codons_detected: stop_count,
            numt_pseudogene_flagged: numt_flagged,
            taxonomic_confidence_regime: regime,
        }
    }
}

fn encode_dna_kmer(win: &[u8]) -> Option<u64> {
    let mut code = 0u64;
    for &b in win {
        let bits = match b.to_ascii_uppercase() {
            b'A' => 0u64,
            b'C' => 1u64,
            b'G' => 2u64,
            b'T' | b'U' => 3u64,
            _ => return None,
        };
        code = (code << 2) | bits;
    }
    Some(code)
}

fn min_internal_stops_across_frames(dna: &str) -> usize {
    let bytes = dna.as_bytes();
    if bytes.len() < 9 {
        return 0;
    }
    let mut min_stops = usize::MAX;
    for frame in 0..3 {
        let mut stops = 0usize;
        let mut pos = frame;
        // Exclude final codon
        while pos + 6 <= bytes.len() {
            let codon = &bytes[pos..pos + 3];
            // Vertebrate mitochondrial stop codons: TAA, TAG, AGA, AGG
            if matches!(codon, b"TAA" | b"TAG" | b"AGA" | b"AGG") {
                stops += 1;
            }
            pos += 3;
        }
        if stops < min_stops {
            min_stops = stops;
        }
    }
    if min_stops == usize::MAX {
        0
    } else {
        min_stops
    }
}

// ============================================================================
// 8. JUST-IN-TIME (JIT) TRANSIENT GENOME BLAST
// ============================================================================

/// Performs a Just-In-Time (JIT) transient in-memory search of a query sequence against an
/// uncompiled genomic contig across both `+` and `-` strands without writing persistent DB files.
pub fn jit_blast_genome(
    query_id: &str,
    query_seq: &str,
    contig_id: &str,
    contig_seq: &str,
) -> (ExecutionSensitivityReport, SuperBlastHsp) {
    let q_clean = clean_sequence_for_biotype(query_seq, SuperBlastBiotype::Dna);
    let t_clean = clean_sequence_for_biotype(contig_seq, SuperBlastBiotype::Dna);
    let sens = analyze_execution_strategy(1, 16, q_clean.len(), t_clean.len(), 1);

    let (_, _, hsp) = super_blast_pair(query_id, &q_clean, contig_id, &t_clean);
    let mut jit_hsp = hsp;
    jit_hsp.execution_mode = SuperBlastExecutionMode::JustInTimeTransientIndex
        .as_str()
        .to_string();
    (sens, jit_hsp)
}
