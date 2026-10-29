//! Core Module Fidelity & GitHub Search Engine Indexing Compliance Suite (`core_fidelity.rs`)
//!
//! Validates the mathematical, algorithmic, and biological fidelity of all four
//! `bio-suite` engines (`fetch`, `phylo`, `genome_target`, `complexity`) using compact,
//! deterministic test fixtures (`tests/fixtures/`), and enforces strict GitHub
//! repository hygiene and Blackbird Code Search indexing constraints:
//! - Hard GitHub file limit: No file >= 100 MB (104,857,600 bytes) anywhere in the repository
//! - GitHub Code Search limit: All tracked source/test/doc files < 384 KB (393,216 bytes)
//! - Valid UTF-8 text without NUL bytes in tracked source files
//! - Compact test footprint: `tests/` maintains <= 10 total files

use bio_suite::complexity::{
    analyze_sequences_parallel, compute_dna_walk_dfa, compute_dna_walk_hurst,
    compute_sequence_complexity, compute_zdna_propensity, kmer_shannon_entropy, lempel_ziv_76,
    scan_g4hunter, scan_sliding_windows_parallel,
};
use bio_suite::dag_pipelines::{
    analyze_idp_sticker_spacer_grammar, compute_3d_contact_laplacian_spectral_metrics,
    compute_apc_corrected_mi_and_potts_couplings, compute_mutational_spectrum_jsd_and_ks,
    evaluate_potts_hamiltonian_escape_barrier, run_chromatin_condensate_scaffolder_agent,
    run_coevo_epistatic_escape_agent, run_master_dag_chain_agent,
    run_somatic_clonal_synthetic_lethal_agent,
};
use bio_suite::fetch::{parse_single_fasta, FastaRecord, GenomicAnnotation, CARNIVORA_COI_PANEL};
use bio_suite::genome_target::{
    compute_wakefield_fine_mapping, intersect_high_complexity_windows, BIO_SUITE_TOOL_REGISTRY,
};
use bio_suite::phylo::{
    align_sequences, compute_k2p_and_tn93, neighbor_joining_newick,
    run_coi_phylo_and_barcode_analysis,
};
use bio_suite::pipelines::{
    analyze_gfa_pangenome_bubble, compute_codon_dn_ds_and_relax_k,
    compute_isoform_splicing_entropy, compute_kl_divergence_bits,
    compute_rin_betweenness_centrality, compute_rscu_cosine_distance, compute_rscu_vector,
    compute_tetranucleotide_spectrum, profile_amp_peptide_biophysics, run_bgc_amp_miner_agent,
    run_hgt_resistome_scanner_agent, run_pangenome_cryptic_target_agent,
    run_pseudogene_trait_loss_agent,
};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

const GITHUB_HARD_FILE_LIMIT_BYTES: u64 = 100 * 1024 * 1024; // 100 MB hard GitHub limit
const GITHUB_SEARCH_INDEX_LIMIT_BYTES: u64 = 384 * 1024; // 384 KB GitHub Blackbird Code Search limit

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Parse a multi-record FASTA fixture file into `Vec<FastaRecord>` using `parse_single_fasta`
/// on each entry so we test the exact FASTA & IUPAC normalization logic used in production.
fn load_multi_fasta_fixture(rel_path: &str) -> Vec<FastaRecord> {
    let full_path = repo_root().join(rel_path);
    let raw = fs::read_to_string(&full_path)
        .unwrap_or_else(|e| panic!("Failed to read fixture {}: {}", full_path.display(), e));

    let mut records = Vec::new();
    for chunk in raw.split('>') {
        let trimmed = chunk.trim();
        if trimmed.is_empty() {
            continue;
        }
        let single_fasta = format!(">{}\n", trimmed);
        let (header, sequence) = parse_single_fasta(&single_fasta)
            .unwrap_or_else(|| panic!("Failed to parse FASTA entry in {}", rel_path));

        // Header format: "ACCESSION Genus species ..."
        let words: Vec<&str> = header.split_whitespace().collect();
        let id = words.first().copied().unwrap_or("UNKNOWN").to_string();
        let taxon = if words.len() >= 3 {
            format!("{} {}", words[1], words[2])
        } else {
            "Unknown sp.".to_string()
        };
        let family = if taxon.starts_with("Ursus") {
            "Ursidae".to_string()
        } else if taxon.starts_with("Canis") || taxon.starts_with("Vulpes") {
            "Canidae".to_string()
        } else {
            "Hominidae".to_string()
        };

        records.push(FastaRecord {
            id: id.clone(),
            accession: id,
            taxon: taxon.clone(),
            common_name: taxon,
            family,
            gene: "COX1".to_string(),
            description: header,
            sequence,
            cached_from_disk: true,
        });
    }
    records
}

#[test]
fn test_fetch_fasta_parser_and_iupac_reading_frame_preservation() {
    // Verify that `parse_single_fasta` preserves exact reading-frame length when encountering
    // IUPAC degenerate bases (R, Y, W, S, K, M, B, D, H, V, N) by converting them to 'N'
    // instead of dropping characters and inducing a frameshift.
    let raw_with_iupac = ">NC_009971.1_test Ursus thibetanus COX1\n\
                          ATGTTCATTAATCGWCTGACTATYCTCAACTAACCACAAAGACATTGGAACACTATACCT\n\
                          ATTATTCGGCGCATGAGCTGGAGTCCTAGGCACAGCTCTAAGCCTCCTTATTCGA\n";
    let (header, seq) = parse_single_fasta(raw_with_iupac).expect("Should parse valid FASTA");
    assert_eq!(header, "NC_009971.1_test Ursus thibetanus COX1");
    assert_eq!(
        seq.len(),
        115,
        "Sequence length must preserve all 115 bases including IUPAC 'W' and 'Y' converted to 'N'"
    );
    assert!(
        seq.contains('N'),
        "IUPAC degenerate bases 'W' and 'Y' must be normalized to 'N' in-frame"
    );

    // Also test against the compact Carnivora COI fixture (4 records, 300 bp each = 100 codons)
    let records = load_multi_fasta_fixture("tests/fixtures/carnivora_coi_mini.fasta");
    assert_eq!(records.len(), 4, "Fixture must contain 4 Carnivora records");
    for rec in &records {
        assert_eq!(
            rec.sequence.len(),
            300,
            "Every Carnivora mini COI fixture sequence must be 300 bp (100 codons) for {}",
            rec.id
        );
        assert_eq!(rec.sequence.len() % 3, 0);
    }

    // Check Ursus arctos record #0 specifically has 'N' from IUPAC 'W', 'L', 'R', 'Y' in-frame
    assert!(
        records[0].sequence.contains('N'),
        "Ursus arctos fixture IUPAC degenerate bases must be preserved as 'N' in-frame"
    );

    // Verify curated 16-specimen panel metadata integrity
    assert_eq!(CARNIVORA_COI_PANEL.len(), 16);
}

#[test]
fn test_phylo_alignment_k2p_tn93_nj_tree_and_abgd_barcoding_gap() {
    let records = load_multi_fasta_fixture("tests/fixtures/carnivora_coi_mini.fasta");
    let (_, aligned) = align_sequences(&records);
    assert_eq!(aligned.len(), 4);
    for a in &aligned {
        assert_eq!(a.aligned_seq.len(), 300);
    }

    // Direct pairwise check between the two intraspecific Canis lupus isolates (records 2 and 3)
    let (p_intra, k2p_intra, tn93_intra, _ti_intra, _tv_intra, _) = compute_k2p_and_tn93(
        aligned[2].aligned_seq.as_bytes(),
        aligned[3].aligned_seq.as_bytes(),
    );
    assert!(
        p_intra < 0.01,
        "Intraspecific Canis lupus vs Canis lupus familiaris p-distance should be < 1%, got {}",
        p_intra
    );
    assert!(k2p_intra >= p_intra);
    assert!(tn93_intra >= p_intra);

    // Direct pairwise check between interfamily Ursus arctos (0) and Canis lupus (2)
    let (p_inter, k2p_inter, tn93_inter, ti_inter, tv_inter, _) = compute_k2p_and_tn93(
        aligned[0].aligned_seq.as_bytes(),
        aligned[2].aligned_seq.as_bytes(),
    );
    assert!(
        k2p_inter > 0.10,
        "Interfamily Ursidae vs Canidae K2P distance must exceed 10%, got {}",
        k2p_inter
    );
    assert!(k2p_inter > p_inter);
    assert!(tn93_inter > p_inter);
    assert!(ti_inter > tv_inter, "Mitochondrial COI must exhibit transition bias (Ti > Tv)");

    // Full phylogenetic & ABGD barcoding analysis on the 4-specimen fixture panel
    let phylo = run_coi_phylo_and_barcode_analysis(&records, 200);
    assert_eq!(phylo.alignment_length_bp, 300);
    assert!(phylo.conserved_sites > 220);
    assert!(phylo.variable_sites > 30);
    assert_eq!(phylo.pairwise_comparisons.len(), 6); // 4 * 3 / 2

    // Verify clear barcoding gap between intraspecific (< 1%) and interspecific pairs
    assert!(
        phylo.barcode_gap.max_intraspecific_k2p < 0.025,
        "Max intraspecific K2P ({}) must fall below the 2.5% threshold",
        phylo.barcode_gap.max_intraspecific_k2p
    );

    // Verify Newick tree syntax and presence of all 4 labels
    assert!(phylo.newick_tree.starts_with('('));
    assert!(phylo.newick_tree.ends_with(");"));
    for label in &phylo.labels {
        assert!(
            phylo.newick_tree.contains(label),
            "Newick tree must contain taxon label {}",
            label
        );
    }

    // Verify standalone Neighbor-Joining tree on a 3-taxon additive matrix
    let d_mat = vec![
        vec![0.0, 0.02, 0.10],
        vec![0.02, 0.0, 0.10],
        vec![0.10, 0.10, 0.0],
    ];
    let lbls = vec![
        "TaxonA".to_string(),
        "TaxonB".to_string(),
        "TaxonC".to_string(),
    ];
    let nwk = neighbor_joining_newick(&d_mat, &lbls);
    assert!(nwk.contains("TaxonA") && nwk.contains("TaxonB") && nwk.contains("TaxonC"));
}

#[test]
fn test_complexity_entropy_lz76_dfa_hurst_g4hunter_and_zdna_fidelity() {
    let motif_records = load_multi_fasta_fixture("tests/fixtures/chr1_agrn_motif_mini.fasta");
    assert_eq!(motif_records.len(), 1);
    let chr1_seq = &motif_records[0].sequence;
    assert_eq!(chr1_seq.len(), 1800);

    // Slice the 1,800 bp fixture into its three constituent 360 bp structural regimes:
    // - 0..360: GC-rich coding regime
    // - 360..720: Alternating purine-pyrimidine Z-DNA tract (TG/CA/CG repeats)
    // - 720..1200: G-quadruplex (G4) + Z-DNA hotspot
    let seq_coding = &chr1_seq[0..360];
    let seq_zdna_repeat = &chr1_seq[360..720];
    let seq_g4_hotspot = &chr1_seq[720..1200];

    let batch_inputs = vec![
        (
            "REGIME_CODING".to_string(),
            "Homo sapiens".to_string(),
            seq_coding.to_string(),
        ),
        (
            "REGIME_ZDNA_REPEAT".to_string(),
            "Homo sapiens".to_string(),
            seq_zdna_repeat.to_string(),
        ),
        (
            "REGIME_G4_HOTSPOT".to_string(),
            "Homo sapiens".to_string(),
            seq_g4_hotspot.to_string(),
        ),
    ];
    let metrics = analyze_sequences_parallel(&batch_inputs);
    assert_eq!(metrics.len(), 3);

    let m_coding = &metrics[0];
    let m_zdna = &metrics[1];
    let m_g4 = &metrics[2];

    // Also verify single-sequence call matches parallel batch call
    let single_check = compute_sequence_complexity("REGIME_G4_HOTSPOT", "Homo sapiens", seq_g4_hotspot);
    assert_eq!(single_check.g4hunter_max_score, m_g4.g4hunter_max_score);

    // 1. G4Hunter propensity must peak strongly in the G4 hotspot regime
    assert!(
        m_g4.g4hunter_max_score >= 2.0,
        "G4 hotspot regime must yield G4Hunter max score >= 2.0, got {}",
        m_g4.g4hunter_max_score
    );
    assert!(
        m_g4.g4_motif_count > 0,
        "G4 hotspot regime must have positive G4 motif hits"
    );

    // 2. Z-DNA alternating purine-pyrimidine propensity must peak in the Z-DNA repeat regime
    assert!(
        m_zdna.zdna_propensity > m_coding.zdna_propensity,
        "Alternating YR Z-DNA tract must have higher Z-DNA propensity ({}) than coding regime ({})",
        m_zdna.zdna_propensity,
        m_coding.zdna_propensity
    );

    // 3. Lempel-Ziv 1976 (LZ76) algorithmic complexity must be lower on simple dinucleotide repeats than coding sequence
    assert!(
        m_zdna.lz76_normalized < m_coding.lz76_normalized,
        "Simple dinucleotide repeat tract must have lower normalized LZ76 complexity ({}) than coding sequence ({})",
        m_zdna.lz76_normalized,
        m_coding.lz76_normalized
    );

    // 4. K-mer Shannon entropy bounds: 0 <= H_k <= 2*k bits, and finite scaling exponents
    for m in &metrics {
        assert!(m.shannon_h1_bits > 0.5 && m.shannon_h1_bits <= 2.0001);
        assert!(m.shannon_h2_bits > 0.8 && m.shannon_h2_bits <= 4.0001);
        assert!(m.shannon_h3_bits > 1.0 && m.shannon_h3_bits <= 6.0001);
        assert!(m.dfa_alpha.is_finite());
        assert!(m.hurst_exponent.is_finite());
    }
    assert!(m_coding.dfa_alpha > 0.1 && m_coding.dfa_alpha < 1.5);

    // 5. Direct unit checks on primitives
    let homopolymer = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    assert_eq!(kmer_shannon_entropy(homopolymer, 1), 0.0);
    let (lz_homo_phrases, lz_homo_norm) = lempel_ziv_76(homopolymer);
    assert!(lz_homo_phrases <= 3 && lz_homo_norm < 0.25);

    let g4_oligo = b"GGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGG";
    let (g4_score, g4_hits) = scan_g4hunter(g4_oligo, 25, 1.2);
    assert!(g4_score > 1.5 && g4_hits > 0);

    let zdna_oligo = b"CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG";
    assert!(compute_zdna_propensity(zdna_oligo) > 10.0);
    assert!(compute_dna_walk_dfa(chr1_seq.as_bytes()) > 0.1);
    assert!(compute_dna_walk_hurst(chr1_seq.as_bytes()) > 0.1);
}

#[test]
fn test_sliding_window_scan_annotation_intersection_and_wakefield_fine_mapping() {
    let motif_records = load_multi_fasta_fixture("tests/fixtures/chr1_agrn_motif_mini.fasta");
    let chr1_seq = &motif_records[0].sequence;

    let windows = scan_sliding_windows_parallel(chr1_seq, 300, 150);
    assert!(!windows.is_empty());

    let annotations = vec![
        GenomicAnnotation {
            feature_id: "GENE_AGRN_HOTSPOT".to_string(),
            feature_type: "CDS".to_string(),
            symbol: "AGRN".to_string(),
            start: 600,
            end: 1250,
            strand: '+',
            description: "AGRN 5' CpG/G4 exon".to_string(),
        },
        GenomicAnnotation {
            feature_id: "REPEAT_VNTR".to_string(),
            feature_type: "repeat_region".to_string(),
            symbol: "VNTR-1p36".to_string(),
            start: 350,
            end: 1150,
            strand: '+',
            description: "Subtelomeric tandem repeat".to_string(),
        },
    ];

    let regions = intersect_high_complexity_windows(&windows, &annotations, 0.5);
    assert!(
        !regions.is_empty(),
        "High-complexity hotspot window must exceed Z >= 0.5 and intersect annotations"
    );
    let top_region = &regions[0];
    assert!(
        top_region.overlapping_cds_count >= 1 || top_region.overlapping_repeat_count >= 1,
        "Top complexity region ({}..{}) must overlap AGRN CDS or VNTR repeat",
        top_region.start,
        top_region.end
    );
    assert!(top_region.narrative_summary.contains("Region-X01"));

    // Test Wakefield ABF fine-mapping (`bio-suite target`)
    let variants = vec![
        ("rs_lead_causal".to_string(), 650usize, 0.48, 0.06),
        ("rs_secondary".to_string(), 710usize, 0.21, 0.08),
        ("rs_null_noise".to_string(), 1400usize, 0.02, 0.09),
    ];
    let credible_set = compute_wakefield_fine_mapping(&variants, 0.04);
    assert_eq!(credible_set.len(), 3);
    assert_eq!(credible_set[0].rsid, "rs_lead_causal");
    assert!(
        credible_set[0].posterior_inclusion_prob > 0.90,
        "Lead variant with Z = 8.0 must have PIP > 0.90, got {}",
        credible_set[0].posterior_inclusion_prob
    );
    let pip_sum: f64 = credible_set
        .iter()
        .map(|v| v.posterior_inclusion_prob)
        .sum();
    assert!(
        (pip_sum - 1.0).abs() < 0.01,
        "Posterior inclusion probabilities across locus must sum to ~1.0, got {}",
        pip_sum
    );

    // Verify all 22 entries in `BIO_SUITE_TOOL_REGISTRY` are GPL-3.0 compatible
    assert_eq!(BIO_SUITE_TOOL_REGISTRY.len(), 22);
    for tool in BIO_SUITE_TOOL_REGISTRY {
        assert!(
            tool.gpl3_compatible,
            "Tool {} must be GPL-3.0 compatible",
            tool.tool_name
        );
        assert!(!tool.spdx_license.is_empty());
        assert!(!tool.native_rust_fallback.is_empty());
    }
}

#[test]
fn test_pipelines_3_to_6_hgt_bgc_pseudoloss_and_pancryptic_fidelity() {
    // 1. Pipeline 3 (`hgt-scan`): Tetranucleotide KL-divergence, RSCU distance & Genomic Island detection
    let gc_rich_host = b"CGCGCCGACGAGTTCCTCGCCGGCCTGGACGCCGTCGGCAAGCGCCTGGCCATCGAGGCCCTCAAGGACGGC";
    let at_rich_island = b"TTAATTAAATGAAATTTAATAATGCATATTATAAATTTGCAATTAATAGAAATTTAATTATGCAATTTAATG";
    let spec_host = compute_tetranucleotide_spectrum(gc_rich_host);
    let spec_island = compute_tetranucleotide_spectrum(at_rich_island);
    let kl_div = compute_kl_divergence_bits(&spec_island, &spec_host);
    assert!(
        kl_div > 0.02,
        "AT-rich foreign island vs GC-rich host must yield positive tetranucleotide KL-divergence, got {}",
        kl_div
    );
    let rscu_host = compute_rscu_vector(gc_rich_host);
    let rscu_island = compute_rscu_vector(at_rich_island);
    let rscu_dist = compute_rscu_cosine_distance(&rscu_island, &rscu_host);
    assert!(
        rscu_dist > 0.05,
        "Foreign island vs host must exhibit distinct RSCU codon usage bias, got {}",
        rscu_dist
    );

    let hgt_report = run_hgt_resistome_scanner_agent().expect("Pipeline 3 hgt-scan must succeed");
    assert_eq!(hgt_report.pipeline_id, "hgt-scan-v1");
    assert_eq!(hgt_report.genomic_islands.len(), 1);
    assert!(hgt_report.genomic_islands[0].max_composite_anomaly_z >= 2.0);
    assert!(!hgt_report.genomic_islands[0].resistome_cargo.is_empty());

    // 2. Pipeline 4 (`bgc-mine`): Amphipathic AMP biophysics, MIBiG dereplication & ARTS self-resistance
    let amp = profile_amp_peptide_biophysics(
        "TEST_AMP_01",
        "GIGKFLHSAKKFGKAFVGEIMNSKWRKVGRK",
        "Amphipathic alpha-helix",
    );
    assert_eq!(amp.length_aa, 31);
    assert!(amp.net_charge_ph7 >= 5.0);
    assert!(amp.hydrophobic_moment_mu_h > 0.25);
    assert!(amp.membrane_disruption_prob > 0.50);

    let bgc_report = run_bgc_amp_miner_agent().expect("Pipeline 4 bgc-mine must succeed");
    assert_eq!(bgc_report.pipeline_id, "bgc-mine-v1");
    assert_eq!(bgc_report.total_clusters_detected, 3);
    assert_eq!(bgc_report.known_clusters_dereplicated, 1);
    assert_eq!(bgc_report.novel_prioritized_clusters, 2);
    assert!(bgc_report.novel_product_dossiers[0].self_resistance_hook.is_some());

    // 3. Pipeline 5 (`pseudo-loss`): Codon dN/dS, premature stop detection & HyPhy RELAX k
    let intact_cds = "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCC";
    let pseudo_cds = "ATGTCCCTGTGGATGTGACTCCTGCCCCTGCTGGAGCTGTTGGCCCTCTGGTGACCTGACACAGCCGCAGCC";
    let (omega, relax_k, _dn, stops) =
        compute_codon_dn_ds_and_relax_k(intact_cds, pseudo_cds, 0.12);
    assert!(stops >= 1, "Pseudogene sequence must contain premature TGA stop codons");
    assert!(omega > 0.20 && relax_k < 1.0);

    let pseudo_report =
        run_pseudogene_trait_loss_agent().expect("Pipeline 5 pseudo-loss must succeed");
    assert_eq!(pseudo_report.pipeline_id, "pseudo-loss-v1");
    assert_eq!(pseudo_report.sequencing_errors_filtered, 1);
    assert_eq!(pseudo_report.validated_losses_and_targets, 3);

    // 4. Pipeline 6 (`pan-cryptic`): GFA walk entropy, Ollivier-Ricci curvature, H_iso & RIN betweenness
    let mini_gfa = "\
H\tVN:Z:1.0
S\tN1\tACGT
S\tN2\tCACA
S\tN3\tGAGA
S\tN4\tTGCA
L\tN1\t+\tN2\t+\t0M
L\tN1\t+\tN3\t+\t0M
L\tN2\t+\tN3\t+\t0M
L\tN2\t+\tN4\t+\t0M
L\tN3\t+\tN4\t+\t0M
P\tP1\tN1+,N2+,N4+\t*
P\tP2\tN1+,N3+,N4+\t*
P\tP3\tN1+,N2+,N3+,N4+\t*
";
    let bubble = analyze_gfa_pangenome_bubble("B1", "TEST_LOCUS", "chr1:1-100", mini_gfa);
    assert_eq!(bubble.gfa_node_count, 4);
    assert_eq!(bubble.gfa_edge_count, 5);
    assert_eq!(bubble.cyclomatic_number_beta1, 2);
    assert!(bubble.walk_shannon_entropy_bits > 0.8);

    let (h_wt, h_sv, jsd) =
        compute_isoform_splicing_entropy(&[0.85, 0.10, 0.05], &[0.30, 0.35, 0.35]);
    assert!(h_sv > h_wt && jsd > 0.10);

    let cb = compute_rin_betweenness_centrality(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
    assert_eq!(cb.len(), 5);
    assert!(cb[2] > cb[0], "Central bridge node in linear chain must have peak betweenness");

    let pan_report =
        run_pangenome_cryptic_target_agent().expect("Pipeline 6 pan-cryptic must succeed");
    assert_eq!(pan_report.pipeline_id, "pan-cryptic-v1");
    assert_eq!(pan_report.stage1_pangenome_bubbles.len(), 2);
    assert!(pan_report.stage4_cryptic_allosteric_pockets[0].disease_isoform_pocket_miner_prob > 0.85);
}

#[test]
fn test_pipelines_7_to_9_and_master_dag_chaining_fidelity() {
    // 1. Pipeline 7 (`co-evo`): APC-corrected Mutual Information & Potts Hamiltonian Escape Barrier
    // Columns (0, 1) co-vary strictly (H<->D, Q<->E, A<->N) while columns 2..5 vary independently
    let msa: Vec<&[u8]> = vec![
        b"HDKMWY", b"HDKLVY", b"HDRMWF", b"HDRLVF", b"QEKMWY", b"QEKLVY", b"QERMWF", b"QERLVF",
        b"ANKMWY", b"ANKLVY", b"ANRMWF", b"ANRLVF",
    ];
    let couplings = compute_apc_corrected_mi_and_potts_couplings(&msa);
    assert!(!couplings.is_empty());
    assert_eq!((couplings[0].0, couplings[0].1), (0, 1));
    assert!(couplings[0].2 > 1.4 && couplings[0].3 > 1.0);

    let (de1_ortho, de_net_ortho, _) =
        evaluate_potts_hamiltonian_escape_barrier(&[2.35, 0.45], &[(0, 1, 3.12)], 0.85);
    let (_, de_net_allo, _) =
        evaluate_potts_hamiltonian_escape_barrier(&[6.40, 4.85, 3.90], &[(0, 1, 0.45)], 0.85);
    assert!(
        de_net_ortho < 0.0 && de_net_allo > 10.0 && de1_ortho > 0.0,
        "Compensatory orthosteric path must be thermodynamically permissive while frustrated allosteric hub blocks escape"
    );

    let coevo_report = run_coevo_epistatic_escape_agent().expect("Pipeline 7 must succeed");
    assert_eq!(coevo_report.pipeline_id, "co-evo-v1");
    assert_eq!(coevo_report.stage4_dead_end_target_dossiers.len(), 2);

    // 2. Pipeline 8 (`chrom-condensate`): 3D Contact Graph Laplacian & IDP Sticker-Spacer Grammar
    let clique_4 = vec![
        vec![0.0, 10.0, 10.0, 10.0],
        vec![10.0, 0.0, 10.0, 10.0],
        vec![10.0, 10.0, 0.0, 10.0],
        vec![10.0, 10.0, 10.0, 0.0],
    ];
    let (vn_s, fiedler_l2, _) = compute_3d_contact_laplacian_spectral_metrics(&clique_4);
    assert!(
        vn_s > 1.4 && fiedler_l2 > 20.0,
        "Dense condensate hub must yield positive Von Neumann entropy ({}) and strong Fiedler connectivity ({})",
        vn_s,
        fiedler_l2
    );

    let (sticker_f, spacer_f, h_gram, csat) = analyze_idp_sticker_spacer_grammar(
        "SYSGYSSYGQSGYSSYSGYSSYQQNTGYSSYRGGRGGRGGYSSYQQSYSGYSSYGQSRGGFGGSRGGYSSYSGYSSYQQS",
    );
    assert!(sticker_f > 0.30 && spacer_f > 0.40 && h_gram > 0.8 && csat < 10.0);

    let cond_report =
        run_chromatin_condensate_scaffolder_agent().expect("Pipeline 8 must succeed");
    assert_eq!(cond_report.pipeline_id, "chrom-condensate-v1");
    assert_eq!(cond_report.stage4_condensate_therapeutics.len(), 2);

    // 3. Pipeline 9 (`clonal-synth`): 96-Channel Mutational JSD & KS Complexity + CRISPR Synthetic Lethality
    let neutral_96 = vec![1.0f64; 96];
    let mut apobec_96 = vec![0.2f64; 96];
    for ch in 48..56 {
        apobec_96[ch] = 25.0;
    }
    let (h_apobec, jsd_apobec, ks_apobec) =
        compute_mutational_spectrum_jsd_and_ks(&apobec_96, &neutral_96);
    assert!(
        h_apobec < 6.0 && jsd_apobec > 0.25 && ks_apobec > 0.30,
        "APOBEC kataegis burst must show significant JSD ({}) and KS divergence ({})",
        jsd_apobec,
        ks_apobec
    );

    let clonal_report =
        run_somatic_clonal_synthetic_lethal_agent().expect("Pipeline 9 must succeed");
    assert_eq!(clonal_report.pipeline_id, "clonal-synth-v1");
    assert_eq!(clonal_report.stage3_synthetic_lethal_pairs.len(), 2);

    // 4. Master DAG Chain (`dag-chain`): Verify all 9 pipelines interlock cleanly
    let dag_report = run_master_dag_chain_agent().expect("Master DAG must succeed");
    assert_eq!(dag_report.total_pipelines_executed, 9);
    assert_eq!(dag_report.dag_edges.len(), 3);
    assert_eq!(dag_report.tool_reuse_matrix.len(), 8);

    // 5. Multi-Pipeline Discovery Sweep (`sweep`): Verify 0★, 1★, 2★, and 3★ rubric gates (including rRNA 16S/18S/23S/ITS calibration) and disk telemetry
    let root = repo_root();
    let sweep_report = bio_suite::sweep::run_multi_pipeline_discovery_sweep(&root)
        .expect("Multi-pipeline sweep must succeed");
    assert_eq!(sweep_report.total_runs_executed, 13);
    assert_eq!(sweep_report.zero_star_excluded_count, 5);
    assert_eq!(sweep_report.one_star_baseline_count, 2);
    assert_eq!(sweep_report.two_star_validation_count, 2);
    assert_eq!(sweep_report.three_star_synthesis_count, 4);
    assert_eq!(sweep_report.rrna_cross_domain_markers.len(), 4);
    assert_eq!(sweep_report.rrna_23s_ptc_epistatic_pairs.len(), 3);
    assert!(sweep_report.disk_snapshots.len() >= 4);

    // 6. Haplotype Clustering & Adaptive Visualization (`haplo-cluster`):
    //    Verify Small-Dataset ASCII viewpoint (h <= 30) and Large-Dataset Super-Cluster + SVG/DOT (h > 30)
    let haplo_report = bio_suite::haplotype::run_haplotype_cluster_pipeline_agent()
        .expect("Haplotype clustering pipeline must succeed");
    assert_eq!(haplo_report.pipeline_id, "haplo-cluster-v1");
    assert_eq!(haplo_report.small_dataset_cohort.diversity_metrics.total_sequences_n, 24);
    assert_eq!(haplo_report.small_dataset_cohort.diversity_metrics.unique_haplotypes_h, 10);
    assert_eq!(haplo_report.small_dataset_cohort.clusters.len(), 3);
    assert!(haplo_report
        .small_dataset_cohort
        .visualizations
        .ascii_cluster_network_viewpoint
        .contains("CLUSTER-01"));
    assert!(haplo_report
        .small_dataset_cohort
        .visualizations
        .ascii_polymorphic_site_matrix
        .contains("POLYMORPHIC SEGREGATING SITES ALIGNMENT"));
    assert_eq!(haplo_report.large_dataset_cohort.diversity_metrics.total_sequences_n, 120);
    assert_eq!(haplo_report.large_dataset_cohort.diversity_metrics.unique_haplotypes_h, 40);
    assert_eq!(haplo_report.large_dataset_cohort.clusters.len(), 5);
    assert!(haplo_report
        .large_dataset_cohort
        .visualizations
        .svg_vector_network
        .starts_with("<svg"));
    assert!(haplo_report
        .large_dataset_cohort
        .visualizations
        .graphviz_dot_network
        .contains("graph HaplotypeMinimumSpanningNetwork"));
}

#[test]
fn test_github_search_indexing_and_file_size_rules() {
    let root = repo_root();

    // 1. Verify `.gitattributes` and `.gitignore` exist to configure GitHub Linguist & exclude runtime artifacts
    assert!(
        root.join(".gitattributes").exists(),
        ".gitattributes must exist for GitHub Linguist search indexing"
    );
    assert!(
        root.join(".gitignore").exists(),
        ".gitignore must exist to prevent tracking build or runtime cache blobs"
    );

    // 2. Keep only a small number of files in `tests/` (<= 10 files total)
    let mut test_files = Vec::new();
    collect_files_recursive(&root.join("tests"), &mut test_files);
    assert!(
        !test_files.is_empty() && test_files.len() <= 10,
        "tests/ directory must remain compact (1..=10 files), found {} files",
        test_files.len()
    );

    // 3. Walk all non-ignored repository directories (`src/`, `tests/`, root configs)
    //    and verify NO file exceeds the 100 MB GitHub limit or the 384 KB GitHub Code Search limit.
    let mut tracked_candidates = Vec::new();
    collect_repo_source_files(&root, &mut tracked_candidates);

    // Also include `git ls-files` if git is available
    if let Ok(output) = Command::new("git")
        .arg("ls-files")
        .current_dir(&root)
        .output()
    {
        if output.status.success() {
            let stdout = String::from_utf8_lossy(&output.stdout);
            for line in stdout.lines() {
                let p = root.join(line.trim());
                if p.exists() && !tracked_candidates.contains(&p) {
                    tracked_candidates.push(p);
                }
            }
        }
    }

    assert!(!tracked_candidates.is_empty());

    for file_path in &tracked_candidates {
        let meta = fs::metadata(file_path)
            .unwrap_or_else(|e| panic!("Failed to stat {}: {}", file_path.display(), e));
        let size = meta.len();

        assert!(
            size < GITHUB_HARD_FILE_LIMIT_BYTES,
            "GitHub Hard Limit Violation: {} is {} bytes (>= 100 MB)",
            file_path.display(),
            size
        );
        assert!(
            size < GITHUB_SEARCH_INDEX_LIMIT_BYTES,
            "GitHub Blackbird Code Search Indexing Violation: {} is {} bytes (>= 384 KB limit)",
            file_path.display(),
            size
        );

        // Verify valid UTF-8 for all tracked source/fixture/config files
        let content = fs::read_to_string(file_path).unwrap_or_else(|e| {
            panic!(
                "GitHub Search Indexing requires valid UTF-8 text in tracked file {}: {}",
                file_path.display(),
                e
            )
        });
        assert!(
            !content.contains('\0'),
            "Tracked file {} must not contain NUL bytes",
            file_path.display()
        );
    }
}

fn collect_files_recursive(dir: &Path, out: &mut Vec<PathBuf>) {
    if let Ok(entries) = fs::read_dir(dir) {
        for entry in entries.flatten() {
            let p = entry.path();
            if p.is_dir() {
                collect_files_recursive(&p, out);
            } else if p.is_file() {
                out.push(p);
            }
        }
    }
}

fn collect_repo_source_files(dir: &Path, out: &mut Vec<PathBuf>) {
    if let Ok(entries) = fs::read_dir(dir) {
        for entry in entries.flatten() {
            let p = entry.path();
            let name = entry.file_name().to_string_lossy().to_string();
            if p.is_dir() {
                if matches!(
                    name.as_str(),
                    ".git" | "target" | "cache" | ".cache" | "runs" | "node_modules"
                ) {
                    continue;
                }
                collect_repo_source_files(&p, out);
            } else if p.is_file() {
                out.push(p);
            }
        }
    }
}
