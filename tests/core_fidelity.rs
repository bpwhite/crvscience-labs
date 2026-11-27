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

    // Verify all 33 entries in `BIO_SUITE_TOOL_REGISTRY` are GPL-3.0 compatible
    assert_eq!(BIO_SUITE_TOOL_REGISTRY.len(), 33);
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

    // 5. Multi-Pipeline Discovery Sweep (`sweep`): Verify 0★, 1★, 2★, and 3★ rubric gates (including rRNA 16S/18S/23S/ITS calibration),
    //    empirical vs. procedural burn-in provenance separation, and disk telemetry
    let root = repo_root();
    let sweep_report = bio_suite::sweep::run_multi_pipeline_discovery_sweep(&root)
        .expect("Multi-pipeline sweep must succeed");
    assert_eq!(
        sweep_report.data_provenance,
        "EMPIRICAL DATA: NCBI REFSEQ/GENBANK CACHED"
    );
    assert!(sweep_report.publishable_standalone);
    assert_eq!(sweep_report.total_runs_executed, 32);
    assert_eq!(sweep_report.zero_star_excluded_count, 8);
    assert_eq!(sweep_report.one_star_baseline_count, 3);
    assert_eq!(sweep_report.two_star_validation_count, 4);
    assert_eq!(sweep_report.three_star_synthesis_count, 8);
    assert_eq!(sweep_report.four_star_emergent_count, 6);
    assert_eq!(sweep_report.five_star_significant_count, 3);
    assert_eq!(sweep_report.synthesis_5_eskape_chromosomes.len(), 8);
    assert_eq!(sweep_report.synthesis_6_eukaryote_prokaryote_rrna.len(), 15);
    assert_eq!(sweep_report.rrna_cross_domain_markers.len(), 4);
    assert_eq!(sweep_report.rrna_23s_ptc_epistatic_pairs.len(), 3);
    assert!(sweep_report.disk_snapshots.len() >= 4);

    // 5b. Procedural Burn-In Calibration (`burn-in`): Verify strict "CALIBRATION DATA: PROCEDURALLY GENERATED" stamp and non-publishable gate
    let burnin_report = bio_suite::sweep::run_procedural_burnin_calibration_sweep()
        .expect("Procedural burn-in calibration must succeed");
    assert_eq!(
        burnin_report.data_provenance,
        "CALIBRATION DATA: PROCEDURALLY GENERATED"
    );
    assert!(!burnin_report.publishable_standalone);
    assert_eq!(burnin_report.total_trials, 14);
    assert_eq!(burnin_report.zero_star_failures, 4);
    assert_eq!(burnin_report.four_plus_star_signals, 4);
    for t in &burnin_report.trials {
        assert_eq!(t.data_provenance, "CALIBRATION DATA: PROCEDURALLY GENERATED");
        assert!(!t.publishable_standalone);
    }

    // 6. Haplotype Clustering & Adaptive Visualization (`haplo-cluster`):
    //    Verify Small-Dataset ASCII viewpoint (h <= 30) and Large-Dataset Super-Cluster + SVG/DOT
    let haplo_report = bio_suite::haplotype::run_haplotype_cluster_pipeline_agent()
        .expect("Haplotype clustering pipeline must succeed");
    assert_eq!(haplo_report.pipeline_id, "haplo-cluster-v1");
    assert_eq!(haplo_report.small_dataset_cohort.diversity_metrics.total_sequences_n, 5);
    assert_eq!(haplo_report.small_dataset_cohort.diversity_metrics.unique_haplotypes_h, 5);
    assert!(haplo_report.small_dataset_cohort.clusters.len() >= 2);
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
    assert_eq!(haplo_report.large_dataset_cohort.diversity_metrics.total_sequences_n, 4);
    assert_eq!(haplo_report.large_dataset_cohort.diversity_metrics.unique_haplotypes_h, 4);
    assert_eq!(haplo_report.large_dataset_cohort.clusters.len(), 3);
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
        let is_preprint_pdf = file_path.extension().and_then(|e| e.to_str()) == Some("pdf")
            && file_path
                .parent()
                .and_then(|p| p.file_name())
                .and_then(|n| n.to_str())
                == Some("Preprints");
        if is_preprint_pdf {
            let bytes = fs::read(file_path)
                .unwrap_or_else(|e| panic!("Failed to read PDF {}: {}", file_path.display(), e));
            assert!(
                bytes.starts_with(b"%PDF-"),
                "Preprint archive PDF {} must begin with %PDF- magic header",
                file_path.display()
            );
            continue;
        }
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

#[test]
fn test_disease_genes_32_panel_atlas_and_pleiotropy_fidelity() {
    let atlas = bio_suite::disease_panels::build_disease_genes_atlas();
    assert_eq!(atlas.total_panels, 32, "Expected 32 curated disease panels");
    assert!(
        atlas.unique_gene_count >= 500,
        "Expected >= 500 unique human disease genes across 32 panels, got {}",
        atlas.unique_gene_count
    );
    assert!(
        atlas.total_panel_gene_assignments >= 850,
        "Expected >= 850 total gene-panel assignments, got {}",
        atlas.total_panel_gene_assignments
    );
    assert!(atlas.pleiotropic_gene_count_ge_2 >= 150);
    assert!(atlas.high_pleiotropy_hubs_ge_4 >= 10);
    assert!(!atlas.top_pleiotropic_hubs.is_empty());
    assert!(!atlas.top_panel_overlaps.is_empty());

    // Verify APOE tiered discovery & calibration sweep
    let cache = bio_suite::fetch::PoliteCacheManager::default_workspace_cache();
    let ws_root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let apoe_sweep = bio_suite::disease_panels::run_apoe_tiered_discovery_sweep(ws_root, &cache)
        .expect("APOE tiered discovery sweep should succeed");
    assert_eq!(apoe_sweep.target_gene, "APOE");
    assert_eq!(apoe_sweep.pleiotropic_panel_count, 5);
    assert_eq!(apoe_sweep.total_runs_executed, 10);
    assert_eq!(apoe_sweep.zero_star_burnin_count, 3);
    assert_eq!(apoe_sweep.zero_star_empirical_reject_count, 2);
    assert_eq!(apoe_sweep.unpublished_1_to_3_star_calibration_count, 3);
    assert_eq!(apoe_sweep.published_4_to_5_star_discovery_count, 2);
    assert_eq!(apoe_sweep.chr19q13_real_window_telemetry.len(), 6);
    assert_eq!(apoe_sweep.isoform_epistasis_matrix.len(), 6);
    // Ensure all 0★ and 1★-3★ runs are strictly marked unpublished
    for run in &apoe_sweep.tiered_runs {
        if run.star_rating <= 3 {
            assert!(
                !run.published_article && !run.publishable_standalone,
                "Run {} ({}★) must remain unpublished",
                run.run_id,
                run.star_rating
            );
        } else {
            assert!(
                run.published_article && run.publishable_standalone,
                "Run {} ({}★) must be marked published",
                run.run_id,
                run.star_rating
            );
        }
    }

    // Verify Avian 4-Locus (16S, 18S, 28S, ITS) Tiered Systematics Sweep
    let avian_sweep =
        bio_suite::avian_systematics::run_avian_tiered_systematics_sweep(ws_root, &cache)
            .expect("Avian 4-locus tiered systematics sweep should succeed");
    assert_eq!(avian_sweep.markers_evaluated.len(), 4);
    assert_eq!(avian_sweep.supermatrix_partition_scheme.len(), 6);
    assert_eq!(avian_sweep.total_runs_executed, 10);
    assert_eq!(avian_sweep.zero_star_burnin_count, 3);
    assert_eq!(avian_sweep.zero_star_empirical_reject_count, 2);
    assert_eq!(avian_sweep.unpublished_1_to_3_star_calibration_count, 3);
    assert_eq!(avian_sweep.published_4_to_5_star_discovery_count, 2);
    assert!(!avian_sweep.empirical_locus_telemetry.is_empty());
    assert!(avian_sweep
        .combined_partitioned_focal_alignment_ascii
        .contains("P1:18S"));
    assert!(avian_sweep
        .combined_partitioned_focal_alignment_ascii
        .contains("P4:16S-mtDNA"));
    for run in &avian_sweep.tiered_runs {
        if run.star_rating <= 3 {
            assert!(
                !run.published_article && !run.publishable_standalone,
                "Avian run {} ({}★) must remain unpublished",
                run.run_id,
                run.star_rating
            );
        } else {
            assert!(
                run.published_article && run.publishable_standalone,
                "Avian run {} ({}★) must be marked published",
                run.run_id,
                run.star_rating
            );
        }
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
                    ".git" | "target" | "cache" | ".cache" | "runs" | "node_modules" | "tools"
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

#[test]
fn test_native_aligner_default_policy_and_gotoh_affine_gap_fidelity() {
    use bio_suite::fetch::FastaRecord;
    use bio_suite::phylo::{align_sequences_with_mode, AlignerMode};

    // 1. Default policy MUST be AlignerMode::Native
    assert_eq!(
        AlignerMode::default(),
        AlignerMode::Native,
        "Default aligner mode must be Native unless an external tool is explicitly specified"
    );
    assert_eq!(AlignerMode::from_cli_arg("native"), AlignerMode::Native);
    assert_eq!(AlignerMode::from_cli_arg("mafft"), AlignerMode::Mafft);

    // 2. Test Gotoh 3-state affine-gap progressive MSA on sequences with internal indels
    let records = vec![
        FastaRecord {
            id: "SEQ_REF.1".to_string(),
            accession: "SEQ_REF.1".to_string(),
            taxon: "Taxon_A".to_string(),
            common_name: "Ref".to_string(),
            family: "FamA".to_string(),
            gene: "ITS".to_string(),
            description: "SEQ_REF.1".to_string(),
            sequence: "ATCGATCGATCGAAAACCCCGGGGTTTTATCGATCGATCG".to_string(),
            cached_from_disk: true,
        },
        FastaRecord {
            id: "SEQ_DEL.1".to_string(),
            accession: "SEQ_DEL.1".to_string(),
            taxon: "Taxon_B".to_string(),
            common_name: "Del4".to_string(),
            family: "FamA".to_string(),
            gene: "ITS".to_string(),
            description: "SEQ_DEL.1".to_string(),
            // 4-bp deletion of CCCC in the middle
            sequence: "ATCGATCGATCGAAAAGGGGTTTTATCGATCGATCG".to_string(),
            cached_from_disk: true,
        },
        FastaRecord {
            id: "SEQ_INS.1".to_string(),
            accession: "SEQ_INS.1".to_string(),
            taxon: "Taxon_C".to_string(),
            common_name: "Ins3".to_string(),
            family: "FamA".to_string(),
            gene: "ITS".to_string(),
            description: "SEQ_INS.1".to_string(),
            // 3-bp insertion of TTT after GGGG
            sequence: "ATCGATCGATCGAAAACCCCGGGGTTTTTTATCGATCGATCG".to_string(),
            cached_from_disk: true,
        },
    ];

    let (_method, aligned, telemetry) = align_sequences_with_mode(&records, AlignerMode::Native);
    assert_eq!(telemetry.aligner_mode, "Native");
    assert_eq!(aligned.len(), 3);
    assert_eq!(aligned[0].aligned_seq.len(), aligned[1].aligned_seq.len());
    assert_eq!(aligned[1].aligned_seq.len(), aligned[2].aligned_seq.len());
    // Single contiguous 4-bp gap should be opened in SEQ_DEL.1 ("----")
    assert!(
        aligned[1].aligned_seq.contains("----"),
        "Expected affine-gap aligner to group the 4-bp deletion into a contiguous '----' gap, got {}",
        aligned[1].aligned_seq
    );
    assert!(telemetry.sum_of_pairs_identity > 0.90);

    // 3. Verify ModelTestMode default policy == ModelTestMode::Native and 24-model evaluation
    use bio_suite::modeltest::{run_modeltest_with_mode, ModelTestMode};
    assert_eq!(
        ModelTestMode::default(),
        ModelTestMode::Native,
        "Default modeltest mode must be Native unless an external tool is explicitly specified"
    );
    assert_eq!(ModelTestMode::from_cli_arg("native"), ModelTestMode::Native);
    assert_eq!(
        ModelTestMode::from_cli_arg("jmodeltest"),
        ModelTestMode::JModelTest
    );
    let mt_res = run_modeltest_with_mode(&aligned, ModelTestMode::Native);
    assert_eq!(mt_res.engine_mode, "Native");
    assert_eq!(mt_res.models_evaluated.len(), 24);
    assert!(mt_res.best_bic_ln_l < 0.0);

    // 4. Verify TreeBuilderMode & CoalescentMode default policies == Native and Quartet MSC fidelity
    use bio_suite::coalescent::{
        compute_abba_baba_patterson_d, compute_coalescent_units_branch_length,
        compute_local_posterior_probability, run_native_quartet_coalescent, CoalescentMode,
    };
    use bio_suite::tree_builder::TreeBuilderMode;

    assert_eq!(
        TreeBuilderMode::default(),
        TreeBuilderMode::Native,
        "Default tree builder mode must be Native unless an external tool is explicitly specified"
    );
    assert_eq!(
        CoalescentMode::default(),
        CoalescentMode::Native,
        "Default coalescent mode must be Native unless --coalescent astral is explicitly specified"
    );
    assert_eq!(CoalescentMode::from_cli_arg("native"), CoalescentMode::Native);
    assert_eq!(CoalescentMode::from_cli_arg("astral"), CoalescentMode::Astral);

    let gts = vec![
        "(((TaxA:0.01,TaxB:0.01):0.03,TaxC:0.02):0.04,(TaxD:0.01,TaxE:0.01):0.04);".to_string(),
        "(((TaxA:0.01,TaxB:0.01):0.03,TaxC:0.02):0.04,(TaxD:0.01,TaxE:0.01):0.04);".to_string(),
        "(((TaxA:0.01,TaxB:0.01):0.03,TaxC:0.02):0.04,(TaxD:0.01,TaxE:0.01):0.04);".to_string(),
        "(((TaxA:0.01,TaxC:0.01):0.02,TaxB:0.02):0.04,(TaxD:0.01,TaxE:0.01):0.04);".to_string(),
    ];
    let (msc_res, _sp_tree, _taxa) =
        run_native_quartet_coalescent("TEST-COALESCENT-5TAXA", &gts)
            .expect("Native quartet coalescent must succeed on 5-taxon gene tree ensemble");
    assert_eq!(msc_res.taxa_count, 5);
    assert_eq!(msc_res.gene_tree_count, 4);
    assert_eq!(msc_res.total_quartets_evaluated, 20);
    assert!(msc_res.satisfied_quartet_score >= 18);
    assert!(msc_res.normalized_quartet_score >= 0.90);

    let lpp = compute_local_posterior_probability(65.0, 17.0, 18.0);
    assert!(lpp > 0.99, "Expected high LPP for q1=0.65 over 100 loci, got {}", lpp);
    let cu = compute_coalescent_units_branch_length(0.6503, 424.0);
    assert!(
        (cu - 0.645).abs() < 0.02,
        "Expected ~0.645 CU for Homininae q1=0.6503, got {}",
        cu
    );

    let abba = compute_abba_baba_patterson_d(
        "TEST-ABBA",
        "P1",
        b"ACGTACGTACGTACGT",
        "P2",
        b"ACGTTCGTACGTTCGT",
        "P3",
        b"ACGTTCGTACGTTCGT",
        "Out",
        b"ACGTACGTACGTACGT",
        "Test introgression P2<->P3",
    );
    assert_eq!(abba.abba_shared_p2_p3_sites, 2);
    assert_eq!(abba.baba_shared_p1_p3_sites, 0);
    assert!((abba.patterson_d_stat - 1.0).abs() < 1e-6);

    // 5. Verify SelectionMode default policy == SelectionMode::Native and 61-State / 60-State Goldman-Yang + BEB fidelity
    use bio_suite::codon_selection::{
        annotate_primate_lysozyme_domain, run_native_codon_selection, CodonModelSpace,
        GoldmanYangEigenSystem, SelectionMode,
    };
    assert_eq!(
        SelectionMode::default(),
        SelectionMode::Native,
        "Default selection mode must be Native unless --selection paml is explicitly specified"
    );
    assert_eq!(SelectionMode::from_cli_arg("native"), SelectionMode::Native);
    assert_eq!(SelectionMode::from_cli_arg("paml"), SelectionMode::Paml);
    assert_eq!(SelectionMode::from_cli_arg("codeml"), SelectionMode::Paml);

    let space_univ = CodonModelSpace::from_alignment(&[], 0);
    let space_mt = CodonModelSpace::from_alignment(&[], 1);
    assert_eq!(space_univ.sense_to_64.len(), 61);
    assert_eq!(space_mt.sense_to_64.len(), 60);

    let eig_test = GoldmanYangEigenSystem::build(&space_univ, 2.5, 0.3, None);
    let p_mat = eig_test.transition_matrix(0.10);
    assert_eq!(p_mat.len(), 61 * 61);
    let row0_sum: f64 = p_mat[0..61].iter().sum();
    assert!((row0_sum - 1.0).abs() < 1e-6);

    let mini_codon_seqs = vec![
        ("Tax1".to_string(), "AAGGTCTTTGAAAGGTGTGAGTTGGCCAGAACTCTG".to_string()),
        ("Tax2".to_string(), "AAGGTCTTTGAAAGGTGTGAGTTGGCCAGAACTCTG".to_string()),
        ("Tax3".to_string(), "AAGATCTTTGAAAGGTGTGAGTTGGCCAGAACTCTG".to_string()),
        ("Tax4".to_string(), "AAGATCTTTGAAAGGTGTGAGTTGGCCAGAACTCTG".to_string()),
    ];
    let mini_tree = "((Tax1:0.02,Tax2:0.02):0.05,(Tax3:0.02,Tax4:0.02) #1:0.08);";
    let sel_res = run_native_codon_selection(
        &mini_codon_seqs,
        mini_tree,
        0,
        true,
        &annotate_primate_lysozyme_domain,
    );
    assert_eq!(sel_res.sense_codon_states, 61);
    assert_eq!(sel_res.codon_sites, 12);
    assert!(sel_res.m0_lnl < 0.0);
    assert!(sel_res.m2a_lnl >= sel_res.m1a_lnl - 1e-4);

    // 6. Verify OrthologyMode default policy == OrthologyMode::Native and Smith-Waterman BLOSUM62 + RBH Orthogroup fidelity
    use bio_suite::orthology::{
        build_cohort4_euk_globin_apoe_paralogs, run_native_orthology_engine,
        smith_waterman_blosum62_local, OrthologyMode,
    };
    assert_eq!(
        OrthologyMode::default(),
        OrthologyMode::Native,
        "Default orthology mode must be Native unless --orthology diamond is explicitly specified"
    );
    assert_eq!(OrthologyMode::from_cli_arg("native"), OrthologyMode::Native);
    assert_eq!(OrthologyMode::from_cli_arg("diamond"), OrthologyMode::Diamond);

    let c4_prots = build_cohort4_euk_globin_apoe_paralogs();
    assert_eq!(c4_prots.len(), 20);
    let sw_self = smith_waterman_blosum62_local(&c4_prots[0], &c4_prots[0], 3997, 20)
        .expect("Self hit on HUMAN_HBA1 must succeed");
    assert_eq!(sw_self.raw_sw_score, 733);
    assert!((sw_self.bit_score - 286.0).abs() <= 1.0);

    let (orth_res, _hits) = run_native_orthology_engine(&c4_prots);
    assert_eq!(orth_res.protein_count, 20);
    assert_eq!(orth_res.inferred_orthogroups_count, 5);
    assert!(orth_res.in_paralog_pairs_count >= 2);

    // 7. Verify SyntenyMode default policy == SyntenyMode::Native and 4-Stage (w,k)-Minimizer + Concave-Gap Chaining + SV Breakpoint fidelity
    use bio_suite::synteny::{
        reverse_complement_dna, run_native_synteny_engine, sketch_canonical_minimizers,
        SyntenyMode,
    };
    assert_eq!(
        SyntenyMode::default(),
        SyntenyMode::Native,
        "Default synteny mode must be Native unless --synteny minimap2 is explicitly specified"
    );
    assert_eq!(SyntenyMode::from_cli_arg("native"), SyntenyMode::Native);
    assert_eq!(SyntenyMode::from_cli_arg("minimap2"), SyntenyMode::Minimap2);
    assert_eq!(SyntenyMode::from_cli_arg("asm20"), SyntenyMode::Minimap2);

    // Construct a deterministic 3,000 bp non-repetitive target sequence and a query with a 1,000 bp inversion in the middle
    let mut target_synth = String::with_capacity(3000);
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    for _ in 0..3000 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let b = match state & 3 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        };
        target_synth.push(b);
    }
    let seeds_t = sketch_canonical_minimizers(&target_synth, 5, 13);
    assert!(!seeds_t.is_empty(), "Minimizer sketch must yield canonical seeds");

    let mut query_synth = String::with_capacity(3000);
    query_synth.push_str(&target_synth[0..1000]);
    query_synth.push_str(&reverse_complement_dna(&target_synth[1000..2000]));
    query_synth.push_str(&target_synth[2000..3000]);

    let syn_res = run_native_synteny_engine(
        "SYNTH_TARGET_3KB",
        &target_synth,
        "SYNTH_QUERY_INV1KB",
        &query_synth,
        "test",
    );
    assert!(
        syn_res.primary_chains_count >= 3,
        "Expected 3 primary chains (SYN -> INV -> SYN), got {}",
        syn_res.primary_chains_count
    );
    assert!(
        syn_res.target_coverage_pct >= 95.0,
        "Expected >=95% target coverage across 3 kb synthetic inversion benchmark, got {}%",
        syn_res.target_coverage_pct
    );
    assert!(
        syn_res.sv_calls.iter().any(|sv| sv.sv_class == "INV" && sv.strand == '-'),
        "Expected SyRI split-chain caller to identify '-' strand INV block"
    );

    // 8. Verify GeneCallerMode default policy == GeneCallerMode::Native and 4-Stage ORF Graph + Shine-Dalgarno RBS scoring fidelity
    use bio_suite::gene_caller::{
        encode_prodigal_bitmaps, shine_dalgarno_exact, GeneCallerMode,
    };
    assert_eq!(
        GeneCallerMode::default(),
        GeneCallerMode::Native,
        "Default gene caller mode must be Native"
    );
    let rbs_test_dna = "AATTAAAGGAGGTTTTTATTATGAAACGTTTAATT";
    let (seq_bm, _rseq_bm, _useq_bm, slen_bm, _gc_bm) = encode_prodigal_bitmaps(rbs_test_dna);
    assert_eq!(slen_bm, rbs_test_dna.len());
    let mut rwt = [0.0f64; 28];
    for (idx, w) in rwt.iter_mut().enumerate() {
        *w = idx as f64;
    }
    let sd_bin = shine_dalgarno_exact(&seq_bm, 6, 20, &rwt);
    assert_eq!(
        sd_bin, 27,
        "Canonical 6-mer AGGAGG with 8 bp spacer must map to highest RBS bin 27, got {}",
        sd_bin
    );

    // 9. Verify ProfileHmmMode default policy == ProfileHmmMode::Native and 4-Stage Plan-7 Forward/Backward + Optimal Accuracy + Null2 fidelity
    use bio_suite::profile_hmm::{
        build_native_plan7_from_msa, score_sequence_plan7_native, ProfileHmmMode,
    };
    assert_eq!(
        ProfileHmmMode::default(),
        ProfileHmmMode::Native,
        "Default Profile HMM mode must be Native unless --hmm hmmer3 is explicitly specified"
    );
    assert_eq!(ProfileHmmMode::from_cli_arg("native"), ProfileHmmMode::Native);
    assert_eq!(ProfileHmmMode::from_cli_arg("hmmer3"), ProfileHmmMode::Hmmer3);

    let msa_bzip = [
        ("GCN4_YEAST", "DPAAKRARNTEAARRSRARKLQRMKQLEDKVEELLSKNYHLENEVARLKKLVGER"),
        ("JUN_HUMAN",  "ESQERIKAERKRMRNRIAASKCRKRKLERIARLEEKVKTLKAQNSELASTANMLR"),
        ("FOS_HUMAN",  "SPEEEEKRRIRRERNKMAAAKCRNRRRELTDTLQAETDQLEDEKSALQTEIANLL"),
    ];
    let hmm_bzip = build_native_plan7_from_msa(
        "TEST_bZIP",
        "PF00170",
        "Test bZIP Leucine Zipper",
        &msa_bzip,
        0.5,
    )
    .expect("Native MSA-to-Plan-7 builder must succeed");
    assert_eq!(hmm_bzip.m, 55);
    let hit_gcn4 = score_sequence_plan7_native(
        &hmm_bzip,
        "GCN4_TARGET",
        "Yeast GCN4",
        "MSEYQPSLFALNDPAAKRARNTEAARRSRARKLQRMKQLEDKVEELLSKNYHLENEVARLKKLVGER",
        10.0,
        1.0,
    );
    assert!(
        hit_gcn4.seq_bitscore > 30.0,
        "Expected strong Plan-7 Forward log-odds score (>30 bits) on true bZIP target, got {}",
        hit_gcn4.seq_bitscore
    );
    assert!(!hit_gcn4.domains.is_empty());
    assert!(hit_gcn4.domains[0].mean_post_acc > 0.70);

    // 10. Verify Native 4-Stage rGFA Pangenome Constructor, Onodera-Sadakane-Shibuya Superbubble Decomposer & Walk Genotyper fidelity
    use bio_suite::pangenome::{
        analyze_rgfa_pangenome, build_native_rgfa_pangenome, build_pangenome_benchmark_cohort,
    };
    let pan_cohort = build_pangenome_benchmark_cohort();
    assert_eq!(pan_cohort.len(), 10);
    assert_eq!(
        pan_cohort.iter().filter(|c| c.kingdom == "Prokaryota").count(),
        5
    );
    assert_eq!(
        pan_cohort.iter().filter(|c| c.kingdom == "Eukaryota").count(),
        5
    );
    let left_flank = "ACGTACGTAGCTAGCTAGCTGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC";
    let del_core = "TTGGAACCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAATT";
    let right_flank = "GGCCAATTGGCCAATTGGCCAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT";
    let h0 = format!("{}{}{}", left_flank, del_core, right_flank);
    let h1 = format!("{}{}", left_flank, right_flank);
    let haps = vec![("H0_REF".to_string(), h0), ("H1_DEL".to_string(), h1)];
    let (segs, lnks, gfa_txt) = build_native_rgfa_pangenome(&haps, 50);
    assert_eq!(segs.len(), 3);
    assert_eq!(lnks.len(), 3);
    assert!(gfa_txt.contains("S\ts1\t"));
    let (summary, bubs) = analyze_rgfa_pangenome(&pan_cohort[0], &haps, &segs, &lnks);
    assert_eq!(summary.num_superbubbles, 1);
    assert_eq!(bubs[0].num_alleles, 2);
    assert_eq!(bubs[0].haplotype_calls, vec![">s2:100:+", "*:0:+"]);
    assert!((bubs[0].bubble_walk_entropy_bits - 1.0).abs() < 1e-6);

    // 11. Verify Tool Pair 14 (`amr-compare` vs. NCBI AMRFinderPlus v4.2.7): fam.tsv hierarchy climb & HGT D_KL separation
    use bio_suite::amr_hgt::{
        climb_fam_hierarchy, compute_hgt_compositional_signature, load_fam_hierarchy,
    };
    let fam_path = repo_root().join("crvscience-labs/tools/amrfinder/data/2026-08-07.1/fam.tsv");
    if fam_path.exists() {
        let fam_map = load_fam_hierarchy(&fam_path).expect("fam.tsv parser must succeed");
        assert!(
            fam_map.len() >= 2200,
            "Expected >= 2,200 fam.tsv hierarchy nodes, got {}",
            fam_map.len()
        );
        // Verify allele node blaTEM-156 stays at blaTEM-156 at 100% identity/coverage
        let exact_node = climb_fam_hierarchy(&fam_map, "blaTEM-156", 100.0, 100.0, false)
            .expect("blaTEM-156 hierarchy lookup must succeed");
        assert_eq!(exact_node.node_id, "blaTEM-156");
    }
    let (kl_kpc, _, _, ctx_kpc, _) =
        compute_hgt_compositional_signature("MSIQHFRVALIPFFAAFCLPVFA", true, "blaKPC-2", "BETA-LACTAM");
    let (kl_pmr, _, _, ctx_pmr, _) =
        compute_hgt_compositional_signature("MHFLRRPISLRQRLILTIGAILL", true, "pmrB_C84R", "COLISTIN");
    assert!(
        kl_kpc > 0.35 && kl_pmr < 0.10 && ctx_kpc.starts_with("HGT_") && ctx_pmr == "CHROMOSOMAL_INTRINSIC_MUTATION",
        "Expected mobile blaKPC-2 island D_KL ({}) >> chromosomal pmrB_C84R D_KL ({})",
        kl_kpc,
        kl_pmr
    );
}

#[test]
fn test_bgc_rule_engine_architecture_guards_and_stockholm_parser() {
    use bio_suite::bgc_mine::{
        call_clusters_for_contig, classify_core_gene, parse_stockholm, BgcGene, DomainHit, RuleSet,
        BGC_FAMILIES,
    };
    use std::collections::HashMap;

    // 1. Stockholm parser: interleaved blocks are concatenated, '.' gaps become '-', GA is read.
    let sto = "# STOCKHOLM 1.0\n#=GF GA   21.50 21.50;\nseqA  AC.DE\nseqB  AC-DE\n\nseqA  FG.HI\nseqB  FGXHI\n//\n";
    let (rows, ga) = parse_stockholm(sto);
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0], ("seqA".to_string(), "AC-DEFG-HI".to_string()));
    assert_eq!(rows[1].1, "AC-DEFGXHI");
    assert!((ga - 21.5).abs() < 1e-6);
    assert_eq!(BGC_FAMILIES.len(), 34);

    let gene = |i: usize, accs: &[&str]| BgcGene {
        seq_id: "t".to_string(),
        gene_index: i,
        beg: 1 + i * 1500,
        end: 1200 + i * 1500,
        strand: '+',
        length_aa: 400,
        hits: accs
            .iter()
            .map(|a| DomainHit { accession: a.to_string(), bitscore: 50.0, i_evalue: 1e-10, pvalue: 1e-13, ali_from: 1, ali_to: 300 })
            .collect(),
        core_class: None,
        support_weight: 0.0,
        smoothed_support: 0.0,
    };
    let weights: HashMap<String, f32> = BGC_FAMILIES.iter().map(|(a, _, w)| (a.to_string(), *w)).collect();

    // 2. Architecture guards: a lone ketosynthase (fatty-acid synthase fabF-like) is not a core PKS gene,
    //    a ketosynthase fused to an acyltransferase is; an adenylation domain alone (acyl-CoA ligase) is not NRPS.
    assert_eq!(classify_core_gene(&gene(0, &["PF00109", "PF02801"]), RuleSet::V2), None);
    assert_eq!(classify_core_gene(&gene(0, &["PF00109", "PF00698"]), RuleSet::V2), Some("PKS-I"));
    assert_eq!(classify_core_gene(&gene(0, &["PF00501"]), RuleSet::V2), None);
    assert_eq!(classify_core_gene(&gene(0, &["PF00668", "PF00501", "PF00550"]), RuleSet::V2), Some("NRPS"));

    // 3. A fatty-acid-synthase-like operon (ACP + acyltransferase + reductase + one lone ketosynthase)
    //    must not produce a cluster under either rule set.
    let mut fas = vec![
        gene(0, &["PF00550"]),
        gene(1, &["PF00698"]),
        gene(2, &["PF08659"]),
        gene(3, &["PF00109", "PF02801"]),
    ];
    assert!(call_clusters_for_contig("t", &mut fas.clone(), &weights, RuleSet::V1).is_empty());
    assert!(call_clusters_for_contig("t", &mut fas, &weights, RuleSet::V2).is_empty());

    // 4. Stand-alone adenylation gene + neighbouring carrier-protein gene (enterobactin-like):
    //    invisible to the V1 rules, called as one NRPS-like cluster by V2.
    let build = || vec![
        gene(0, &[]),
        gene(1, &["PF00501"]),
        gene(2, &["PF00550"]),
        gene(3, &[]),
    ];
    assert!(call_clusters_for_contig("t", &mut build(), &weights, RuleSet::V1).is_empty());
    let v2 = call_clusters_for_contig("t", &mut build(), &weights, RuleSet::V2);
    assert_eq!(v2.len(), 1);
    assert!(v2[0].classes.contains(&"NRPS-like".to_string()));
}

#[test]
fn test_bgc_rule_set_v3_extensions_and_frozen_parameters() {
    use bio_suite::bgc_mine::{
        call_clusters_with_params, classify_core_gene, BgcGene, ClusterParams, DomainHit, RuleSet,
        BGC_FAMILIES, BGC_FAMILIES_EXTRA, BGC_FAMILIES_EXTRA2, FROZEN_V3_PARAMS,
    };
    use std::collections::HashMap;

    let gene = |i: usize, accs: &[&str]| BgcGene {
        seq_id: "t".to_string(),
        gene_index: i,
        beg: 1 + i * 1500,
        end: 1200 + i * 1500,
        strand: '+',
        length_aa: 400,
        hits: accs
            .iter()
            .map(|a| DomainHit { accession: a.to_string(), bitscore: 60.0, i_evalue: 1e-12, pvalue: 1e-15, ali_from: 1, ali_to: 300 })
            .collect(),
        core_class: None,
        support_weight: 0.0,
        smoothed_support: 0.0,
    };
    let mut weights: HashMap<String, f32> = BGC_FAMILIES.iter().map(|(a, _, w)| (a.to_string(), *w)).collect();
    weights.extend(BGC_FAMILIES_EXTRA.iter().map(|(a, _, w)| (a.to_string(), *w)));
    weights.extend(BGC_FAMILIES_EXTRA2.iter().map(|(a, _, w)| (a.to_string(), *w)));
    assert_eq!(BGC_FAMILIES_EXTRA.len(), 14);
    assert_eq!(BGC_FAMILIES_EXTRA2.len(), 8);
    assert_eq!(weights.len(), 56);

    // Published defaults are untouched by the parameterisation.
    let d = ClusterParams::default();
    assert_eq!((d.support_window, d.extend_gap, d.pad), (4, 5000, 2000));
    assert!((d.seed_min - 1.8).abs() < 1e-6);
    // Frozen calibration values.
    assert_eq!((FROZEN_V3_PARAMS.support_window, FROZEN_V3_PARAMS.extend_gap, FROZEN_V3_PARAMS.pad), (6, 5000, 1000));
    assert!((FROZEN_V3_PARAMS.seed_min - 1.2).abs() < 1e-6);

    // V3-only core definitions: bacterial terpene synthase, adenylation + terminal reductase,
    // fungal starter-unit transacylase fused to a ketosynthase.
    assert_eq!(classify_core_gene(&gene(0, &["PF19086"]), RuleSet::V2), None);
    assert_eq!(classify_core_gene(&gene(0, &["PF19086"]), RuleSet::V3), Some("Terpene"));
    assert_eq!(classify_core_gene(&gene(0, &["PF00501", "PF07993"]), RuleSet::V2), None);
    assert_eq!(classify_core_gene(&gene(0, &["PF00501", "PF07993"]), RuleSet::V3), Some("NRPS"));
    assert_eq!(classify_core_gene(&gene(0, &["PF00109", "PF16073"]), RuleSet::V2), None);
    assert_eq!(classify_core_gene(&gene(0, &["PF00109", "PF16073"]), RuleSet::V3), Some("PKS-I"));

    // V4-specific refinements: RamS/AmfS precursor (PF19402) is RiPP under V4, while lone YcaO (PF02624)
    // is classified as RiPP-weak under V4 so it requires neighbouring support (>= 1.8) instead of seeding alone.
    assert_eq!(classify_core_gene(&gene(0, &["PF19402"]), RuleSet::V3), None);
    assert_eq!(classify_core_gene(&gene(0, &["PF19402"]), RuleSet::V4), Some("RiPP"));
    assert_eq!(classify_core_gene(&gene(0, &["PF02624"]), RuleSet::V3), Some("RiPP"));
    assert_eq!(classify_core_gene(&gene(0, &["PF02624"]), RuleSet::V4), Some("RiPP-weak"));

    let lone_ycao = || vec![gene(0, &[]), gene(1, &["PF02624"]), gene(2, &["PF00067"]), gene(3, &[])];
    assert_eq!(call_clusters_with_params("t", &mut lone_ycao(), &weights, RuleSet::V3, &FROZEN_V3_PARAMS).len(), 1);
    assert!(call_clusters_with_params("t", &mut lone_ycao(), &weights, RuleSet::V4, &FROZEN_V3_PARAMS).is_empty());

    // Stand-alone bacterial terpene synthase (PF19086) seeds a cluster under V4 even without extra tailoring hits.
    let lone_terpene = || vec![gene(0, &[]), gene(1, &["PF19086"]), gene(2, &[])];
    assert!(call_clusters_with_params("t", &mut lone_terpene(), &weights, RuleSet::V3, &FROZEN_V3_PARAMS).is_empty());
    let v4_terp = call_clusters_with_params("t", &mut lone_terpene(), &weights, RuleSet::V4, &FROZEN_V3_PARAMS);
    assert_eq!(v4_terp.len(), 1);
    assert!(v4_terp[0].classes.contains(&"Terpene".to_string()));

    // Guards remain: a lone ketosynthase and a lone adenylation gene are still not core genes under V3/V4.
    assert_eq!(classify_core_gene(&gene(0, &["PF00109"]), RuleSet::V3), None);
    assert_eq!(classify_core_gene(&gene(0, &["PF00501"]), RuleSet::V3), None);
    assert_eq!(classify_core_gene(&gene(0, &["PF00109"]), RuleSet::V4), None);
    assert_eq!(classify_core_gene(&gene(0, &["PF00501"]), RuleSet::V4), None);

    // An adenylation gene immediately next to an MbtH partner gene is an NRPS-like cluster under V3 and V4,
    // whereas an adenylation gene separated by 2 intervening genes from a carrier protein (Treponema-like)
    // is rejected by V4's immediate-adjacency guard.
    let build = || vec![gene(0, &[]), gene(1, &["PF00501"]), gene(2, &["PF03621"]), gene(3, &[])];
    assert!(call_clusters_with_params("t", &mut build(), &weights, RuleSet::V2, &FROZEN_V3_PARAMS).is_empty());
    let v3 = call_clusters_with_params("t", &mut build(), &weights, RuleSet::V3, &FROZEN_V3_PARAMS);
    assert_eq!(v3.len(), 1);
    assert!(v3[0].classes.contains(&"NRPS-like".to_string()));
    let v4 = call_clusters_with_params("t", &mut build(), &weights, RuleSet::V4, &FROZEN_V3_PARAMS);
    assert_eq!(v4.len(), 1);
    assert!(v4[0].classes.contains(&"NRPS-like".to_string()));

    let distant_nrps_like = || vec![gene(0, &["PF00501"]), gene(1, &[]), gene(2, &[]), gene(3, &["PF00550"])];
    assert_eq!(call_clusters_with_params("t", &mut distant_nrps_like(), &weights, RuleSet::V3, &FROZEN_V3_PARAMS).len(), 1);
    assert!(call_clusters_with_params("t", &mut distant_nrps_like(), &weights, RuleSet::V4, &FROZEN_V3_PARAMS).is_empty());
}

#[test]
fn test_runie_pipeline_inference_engine_all_3_tiers_and_9_studies() {
    use bio_suite::fetch::PoliteCacheManager;
    use bio_suite::runie::{run_runie_master_pipeline_inference, RUNIE_TOOL_VOCABULARY};

    let cache = PoliteCacheManager::default_workspace_cache();
    let ws_root = Path::new(env!("CARGO_MANIFEST_DIR"));

    // 1. Verify open-ended tool vocabulary (38 calibrated states including mafft-native & bio-suite-supervisor)
    assert_eq!(RUNIE_TOOL_VOCABULARY.len(), 38);

    // 2. Execute full RuNIE GA-MCMC Step-Up Pipeline Synthesis & Multi-Study Validation
    let report = run_runie_master_pipeline_inference(&cache, ws_root)
        .expect("RuNIE master pipeline inference must succeed");

    assert_eq!(report.engine_acronym, "RuNIE");
    assert_eq!(report.total_tool_vocabulary_size, 38);
    assert!(report.autonomous_tools_synthesized_and_calibrated >= 6);
    assert!(report.total_scientific_studies_solved >= 9);
    assert_eq!(report.tier_validation_summaries.len(), 3);

    // 3. Verify all synthesized native Rust domain tools pass head-to-head calibration
    for cal in &report.synthesized_tool_calibrations {
        assert!(
            cal.calibration_gate_passed,
            "Synthesized tool {} must pass mandatory calibration gate",
            cal.tool_id
        );
        assert!(
            cal.pearson_r_vs_reference >= 0.990,
            "Synthesized tool {} Pearson r {} must be >= 0.990",
            cal.tool_id,
            cal.pearson_r_vs_reference
        );
        assert!(cal.negative_zero_beaker_controls >= 1);
        assert!(cal.signal_beakers_earned >= 2);
    }

    // 4. Verify all 3 Experimental Tiers (Tier 1: Exploratory, Tier 2: Validation, Tier 3: Scientific) pass
    for (idx, tier_sum) in report.tier_validation_summaries.iter().enumerate() {
        assert_eq!(tier_sum.tier, (idx + 1) as u8);
        assert!(tier_sum.validation_gate_passed);
        assert!(tier_sum.studies_evaluated >= 9);
        assert!(tier_sum.mean_posterior_probability >= 0.85);
        assert!(tier_sum.mean_log_bayes_factor_bf10 >= 50.0);
    }

    // 5. Verify all 9 canonical scientific study pipelines converge via 4 Step-Up Generations
    for study in &report.study_dossiers {
        assert_eq!(study.step_up_generations.len(), 4);
        assert_eq!(study.top_candidate_pipelines.len(), 4);
        assert!(
            study.winning_posterior_probability >= 0.88,
            "Study #{} ({}) winning MAP posterior {} must be >= 0.88",
            study.study_number,
            study.study_id,
            study.winning_posterior_probability
        );
        assert!(
            study.winning_log_bayes_factor_bf10 >= 80.0,
            "Study #{} ({}) winning ln(BF_10) {} must be >= 80.0",
            study.study_number,
            study.study_id,
            study.winning_log_bayes_factor_bf10
        );
        assert!(study.winning_signal_z_score >= 4.50);
        assert_eq!(study.signal_beakers, 5);
        assert!(!study.emergent_scientific_conclusion.is_empty());

        // Verify Gen 2 (K_min failure state) has lower log-posterior than Gen 1, Gen 3, and Gen 4
        let g1 = &study.step_up_generations[0];
        let g2 = &study.step_up_generations[1];
        let g3 = &study.step_up_generations[2];
        let g4 = &study.step_up_generations[3];
        assert!(!g1.kmin_knowledge_recovered);
        assert!(!g2.kmin_knowledge_recovered);
        assert!(g3.kmin_knowledge_recovered);
        assert!(g4.kmin_knowledge_recovered);
        assert!(g2.log_posterior_unnormalized < g1.log_posterior_unnormalized);
        assert!(g3.log_posterior_unnormalized > g1.log_posterior_unnormalized);
        assert!(g4.log_posterior_unnormalized > g3.log_posterior_unnormalized);
    }
}

#[test]
fn test_pseudo_loss_macse_and_paml_relax_calibration() {
    use bio_suite::pseudo_loss_calibration::{
        build_12_vertebrate_cohorts, run_naive_triplet_codon_aligner,
        run_native_frameshift_coding_aligner, FROZEN_FRAMESHIFT_PARAMS,
    };

    let labs_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let cohorts = build_12_vertebrate_cohorts(labs_dir)
        .expect("12 vertebrate RefSeq cohorts must load from cache/pseudo_loss/refseq_cds_anchors.fasta");
    assert_eq!(cohorts.len(), 12);

    // 1. Verify 6 Development + 6 Untouched Reserve split
    let dev_count = cohorts.iter().filter(|c| c.scope == "development").count();
    let res_count = cohorts.iter().filter(|c| c.scope == "reserve").count();
    assert_eq!((dev_count, res_count), (6, 6));

    // 2. Verify Generation-1 Naive triplet aligner fails on -1 bp frameshifts (GULO DEV-01)
    //    whereas the calibrated 3-frame coding aligner recovers all 4 frameshift sites (codon 24)
    //    and all 4 internal stop codons (codon 52).
    let gulo = &cohorts[0];
    let gen1 = run_naive_triplet_codon_aligner(&gulo.reliable_seqs, &gulo.lr_seqs);
    assert_eq!(gen1.total_frameshift_sites, 0);
    assert!(gen1.total_internal_stop_sites >= 8); // spurious downstream stops from frame shift

    let nat = run_native_frameshift_coding_aligner(
        &gulo.reliable_seqs,
        &gulo.lr_seqs,
        &FROZEN_FRAMESHIFT_PARAMS,
    );
    assert_eq!(nat.total_frameshift_sites, 4);
    assert_eq!(nat.total_internal_stop_sites, 4);
    assert!(nat.shared_across_sister_lineages);
    assert!(!nat.singleton_homopolymer_artifact);

    // 3. Verify singleton homopolymer artifact controls (DEV-06 ABCA4 6G and RES-06 HBB 5A)
    let abca4 = &cohorts[5];
    let nat_abca4 = run_native_frameshift_coding_aligner(
        &abca4.reliable_seqs,
        &abca4.lr_seqs,
        &FROZEN_FRAMESHIFT_PARAMS,
    );
    assert_eq!(nat_abca4.total_frameshift_sites, 1);
    assert!(!nat_abca4.shared_across_sister_lineages);
    assert!(nat_abca4.singleton_homopolymer_artifact);

    let hbb = &cohorts[11];
    let nat_hbb = run_native_frameshift_coding_aligner(
        &hbb.reliable_seqs,
        &hbb.lr_seqs,
        &FROZEN_FRAMESHIFT_PARAMS,
    );
    assert_eq!(nat_hbb.total_frameshift_sites, 1);
    assert!(!nat_hbb.shared_across_sister_lineages);
    assert!(nat_hbb.singleton_homopolymer_artifact);
}

#[test]
fn test_chrom_condensate_scipy_calibration() {
    use bio_suite::chrom_condensate_calibration::run_chrom_condensate_scipy_calibration;

    let report = run_chrom_condensate_scipy_calibration(None)
        .expect("3D Chromatin Laplacian & TAD Insulation calibration vs SciPy v1.18.1 must succeed");
    assert_eq!(report.total_cohorts, 12);
    assert_eq!(report.dev_cohorts, 6);
    assert_eq!(report.reserve_cohorts, 6);
    assert!(report.dev_mean_abs_err_lambda2 <= 1e-4);
    assert!(report.reserve_mean_abs_err_lambda2 <= 1e-4);
    assert!(report.dev_mean_abs_err_norm_lambda2 <= 1e-4);
    assert!(report.reserve_mean_abs_err_norm_lambda2 <= 1e-4);
    assert!(report.dev_mean_abs_err_vn_entropy <= 1e-4);
    assert!(report.reserve_mean_abs_err_vn_entropy <= 1e-4);
    assert!(report.dev_mean_abs_err_insulation <= 1e-4);
    assert!(report.reserve_mean_abs_err_insulation <= 1e-4);
    assert_eq!(report.tad_boundary_exact_match_rate_pct, 100.0);
    assert_eq!(report.regime_concordance_rate_pct, 100.0);
}

#[test]
fn test_clonal_synth_scipy_calibration() {
    use bio_suite::clonal_synth_calibration::run_clonal_synth_vs_scipy_calibration;
    use std::path::Path;

    let report = run_clonal_synth_vs_scipy_calibration(Path::new("."))
        .expect("Somatic Clonal 96-SBS JSD, NNLS & CRISPR Synthetic Lethality calibration vs SciPy v1.18.1 must succeed");
    assert_eq!(report.total_cohorts, 12);
    assert_eq!(report.dev_cohorts, 6);
    assert_eq!(report.reserve_cohorts, 6);
    assert!(report.overall_mae_h96_bits <= 1e-4);
    assert!(report.overall_mae_jsd_bits <= 1e-4);
    assert!(report.overall_mae_nnls_weight <= 1e-4);
    assert!(report.overall_mae_delta_chronos <= 1e-4);
    assert!(report.overall_mae_welch_t <= 1e-4);
    assert_eq!(report.overall_regime_concordance_pct, 100.0);
}

#[test]
fn test_finemap_coloc_scipy_calibration() {
    use bio_suite::finemap_coloc_calibration::run_finemap_coloc_scipy_calibration;

    let report = run_finemap_coloc_scipy_calibration();
    assert_eq!(report.total_cohorts, 12);
    assert_eq!(report.dev_cohorts, 6);
    assert_eq!(report.reserve_cohorts, 6);
    assert_eq!(report.total_variants_evaluated, 240);
    assert!(report.mae_lead_ln_abf <= 1e-4);
    assert!(report.mae_lead_pip_gwas <= 1e-4);
    assert!(report.mae_credible_entropy_bits <= 1e-4);
    assert!(report.mae_coloc_pp_h3 <= 1e-4);
    assert!(report.mae_coloc_pp_h4 <= 1e-4);
    assert_eq!(report.credible_set_concordance_rate, 1.0);
    assert_eq!(report.regime_concordance_rate, 1.0);
    // Verify negative controls DEV-06 (H0 null) and RES-06 (H3 linkage decoy)
    assert_eq!(report.cohorts[5].native_regime, "UnderpoweredNullH0Control");
    assert!(report.cohorts[5].native_pp_h0 >= 0.90);
    assert_eq!(report.cohorts[11].native_regime, "DistinctLinkageH3Regime");
    assert!(report.cohorts[11].native_pp_h3 >= 0.99);
}

#[test]
fn test_mr_twas_scipy_calibration() {
    use bio_suite::mr_twas_calibration::run_mr_twas_scipy_calibration;

    let report = run_mr_twas_scipy_calibration();
    assert_eq!(report.total_cohorts, 12);
    assert_eq!(report.dev_cohorts, 6);
    assert_eq!(report.reserve_cohorts, 6);
    assert_eq!(report.total_instruments_evaluated, 120);
    assert!(report.mae_ivw_beta <= 1e-4);
    assert!(report.mae_ivw_se <= 1e-4);
    assert!(report.mae_cochran_q <= 1e-4);
    assert!(report.mae_egger_slope <= 1e-4);
    assert!(report.mae_egger_intercept <= 1e-4);
    assert!(report.mae_weighted_median_beta <= 1e-4);
    assert_eq!(report.regime_concordance_rate, 1.0);
    // Verify negative controls DEV-06 (CRP null) and RES-06 (HDL horizontal pleiotropy)
    assert_eq!(report.cohorts[5].native_regime, "NullNonCausalInstrumentSet");
    assert_eq!(
        report.cohorts[11].native_regime,
        "HorizontalPleiotropyConfoundedReject"
    );
}

#[test]
fn test_splicing_crispr_scipy_calibration() {
    use bio_suite::splicing_crispr_calibration::run_splicing_crispr_scipy_calibration;

    let report = run_splicing_crispr_scipy_calibration()
        .expect("Alternative Splicing & BAGEL2 CRISPR Essentiality calibration vs. SciPy v1.18.1 must succeed");
    assert_eq!(report.total_cohorts, 12);
    assert_eq!(report.development_cohorts, 6);
    assert_eq!(report.untouched_reserve_cohorts, 6);
    assert_eq!(report.total_replicates_evaluated, 36);
    assert_eq!(report.total_guides_evaluated, 72);
    assert!(report.mean_abs_delta_psi_error <= 1e-4);
    assert!(report.mean_abs_delta_psi_z_error <= 1e-4);
    assert!(report.mean_abs_isoform_jsd_error <= 1e-4);
    assert!(report.mean_abs_bagel_ln_bf_error <= 1e-4);
    assert!(report.mean_abs_posterior_ess_error <= 1e-4);
    assert_eq!(report.regime_concordance_percent, 100.0);
    // Verify negative controls DEV-06 (GAPDH/AAVS1 null) and RES-06 (PVT1 8q24 CN=18 amplicon DSB decoy)
    assert_eq!(
        report.cohorts[5].native_regime,
        "NullHousekeepingIsoformNoiseControl"
    );
    assert_eq!(
        report.cohorts[11].native_regime,
        "OffTargetCopyNumberAmplificationArtifactReject"
    );
}

#[test]
fn test_kmer_complexity_scipy_calibration() {
    use bio_suite::kmer_complexity_calibration::run_kmer_complexity_scipy_calibration;

    let report = run_kmer_complexity_scipy_calibration()
        .expect("k-mer Shannon Entropy, LZ76 NCD & Fractal DNA Walk DFA calibration vs. SciPy v1.18.1 must succeed");
    assert_eq!(report.total_cohorts, 12);
    assert_eq!(report.development_cohorts, 6);
    assert_eq!(report.untouched_reserve_cohorts, 6);
    assert_eq!(report.total_base_pairs_evaluated, 18432);
    assert!(report.mean_abs_shannon_h3_error <= 1e-4);
    assert!(report.mean_abs_shannon_h6_error <= 1e-4);
    assert!(report.mean_abs_lz76_norm_error <= 1e-4);
    assert!(report.mean_abs_lz76_ncd_error <= 1e-4);
    assert!(report.mean_abs_dfa_alpha_error <= 1e-4);
    assert!(report.mean_abs_hurst_h_error <= 1e-4);
    assert_eq!(report.regime_concordance_percent, 100.0);
    // Verify negative controls DEV-06 (homopolymer artifact) and RES-06 (I.I.D. white noise)
    assert_eq!(
        report.cohorts[5].native_regime,
        "HomopolymerOrDinucleotideLowComplexityArtifactReject"
    );
    assert_eq!(
        report.cohorts[11].native_regime,
        "MemorylessIIDBernoulliWhiteNoiseControlReject"
    );
}

#[test]
fn test_super_blast_calibration() {
    use bio_suite::super_blast_calibration::run_super_blast_calibration;

    let report = run_super_blast_calibration()
        .expect("Super-BLAST Tier 1 calibration sweep vs. Biopython v1.88 & SciPy v1.18.1 must succeed");
    assert_eq!(report.total_cohorts, 12);
    assert_eq!(report.development_cohorts, 6);
    assert_eq!(report.untouched_reserve_cohorts, 6);
    assert_eq!(report.biotype_classification_accuracy_pct, 100.0);
    assert!(report.mean_abs_raw_score_error <= 1e-4);
    assert!(report.mean_abs_bit_score_error <= 1e-4);
    assert!(report.mean_abs_pct_identity_error <= 1e-4);
    assert!(report.mean_abs_k2p_distance_error <= 1e-4);
    assert_eq!(report.regime_concordance_percent, 100.0);
    // Verify MSA pre-validation rescue & rejections
    assert_eq!(report.msa_pre_validation_audit[0].status, "AcceptedDirect");
    assert_eq!(report.msa_pre_validation_audit[1].status, "AcceptedAutoRevComp");
    assert!(report.msa_pre_validation_audit[1].reverse_complement_rescued);
    assert!(!report.msa_pre_validation_audit[2].incorporate_into_msa);
    assert!(!report.msa_pre_validation_audit[3].incorporate_into_msa);
    // Verify 16S vs 10-gene tip concordance audit
    assert_eq!(report.phylogeny_16s_10gene_tip_audit.num_secondary_genes_tested, 10);
    assert_eq!(report.phylogeny_16s_10gene_tip_audit.concordant_genes_count, 8);
    assert_eq!(report.phylogeny_16s_10gene_tip_audit.discordant_genes_count, 2);
    // Verify negative controls DEV-06 (random shuffle) and RES-06 (low-complexity ATn repeat)
    assert_eq!(
        report.cohorts[5].native_regime,
        "NonHomologousRandomShuffleControlReject"
    );
    assert_eq!(
        report.cohorts[11].native_regime,
        "LowComplexityCompositionBiasedArtifactReject"
    );
}

#[test]
fn test_ml_tree_bootstrap_estimator_equivalence_and_source_guardrails() {
    use bio_suite::tree_builder::{infer_native_ml_tree_with_bootstrap, PhyloRunModel};

    fn collect_rs_files(dir: &Path, out: &mut Vec<PathBuf>) {
        if let Ok(entries) = fs::read_dir(dir) {
            for ent in entries.flatten() {
                let p = ent.path();
                if p.is_dir() {
                    collect_rs_files(&p, out);
                } else if p.extension().and_then(|e| e.to_str()) == Some("rs") {
                    out.push(p);
                }
            }
        }
    }

    // 1. Repository-Wide Static Architecture Guardrail:
    // Ban ad-hoc uniform-JC/NJ bootstrap helpers outside `tree_builder.rs` and `phylo.rs`
    let src_root = repo_root().join("src");
    let mut rs_files = Vec::new();
    collect_rs_files(&src_root, &mut rs_files);
    assert!(!rs_files.is_empty(), "Expected Rust source files under src/");

    let mut saw_acan_unified_call = false;
    let forbidden_patterns = [
        "fn fast_bionj_splits",
        "fn run_real_multinomial_bootstrap",
        "clamped_p / b_factor",
        "\"FastNJ Boot",
    ];
    for p in &rs_files {
        let fname = p.file_name().and_then(|s| s.to_str()).unwrap_or("");
        let content = fs::read_to_string(p)
            .unwrap_or_else(|e| panic!("Failed reading {}: {}", p.display(), e));
        if fname == "supervisor_acanthochitonidae.rs" {
            if content.contains("infer_native_ml_tree_with_bootstrap") {
                saw_acan_unified_call = true;
            }
        }
        if fname == "tree_builder.rs" || fname == "phylo.rs" {
            continue;
        }
        for pat in &forbidden_patterns {
            assert!(
                !content.contains(pat),
                "Architecture Guardrail Violation in {}: forbidden mismatched bootstrap pattern `{}` found; use `tree_builder::infer_native_ml_tree_with_bootstrap` instead",
                p.display(),
                pat
            );
        }
    }
    assert!(
        saw_acan_unified_call,
        "supervisor_acanthochitonidae.rs must invoke `infer_native_ml_tree_with_bootstrap`"
    );

    // 2. End-to-End API Guardrail Check on 8-taxon Carnivora COI RY3 Panel (>= 5 nontrivial splits)
    let base_recs = load_multi_fasta_fixture("tests/fixtures/carnivora_coi_mini.fasta");
    assert_eq!(base_recs.len(), 4);
    let flip_ti = |b: u8| -> u8 {
        match b.to_ascii_uppercase() {
            b'A' => b'G',
            b'G' => b'A',
            b'C' => b'T',
            b'T' => b'C',
            _ => b'A',
        }
    };
    let mut labels = Vec::with_capacity(8);
    let mut unrec_seqs = Vec::with_capacity(8);
    for (ti, rec) in base_recs.iter().enumerate() {
        let base_bytes = rec.sequence.as_bytes();
        // Create two divergent sister sub-lineages per base sequence with distinct shared + tip substitutions
        for sub_idx in 0..2 {
            let mut bvec = base_bytes.to_vec();
            let clade_stride = 19 + ti * 13;
            for k in (clade_stride..bvec.len()).step_by(43) {
                if k % 3 != 2 {
                    bvec[k] = flip_ti(bvec[k]);
                }
            }
            let tip_offset = 5 + ti * 7 + sub_idx * 17;
            for k in (tip_offset..bvec.len()).step_by(29) {
                bvec[k] = flip_ti(bvec[k]);
            }
            labels.push(format!("{}_{}", rec.id.replace('.', "_"), sub_idx + 1));
            unrec_seqs.push(String::from_utf8(bvec).unwrap());
        }
    }
    let ry3_seqs: Vec<String> = unrec_seqs
        .iter()
        .map(|s| {
            s.bytes()
                .enumerate()
                .map(|(i, b)| {
                    if i % 3 == 2 {
                        match b.to_ascii_uppercase() {
                            b'A' | b'G' => 'A',
                            b'C' | b'T' => 'C',
                            _ => '-',
                        }
                    } else if matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T') {
                        b.to_ascii_uppercase() as char
                    } else {
                        '-'
                    }
                })
                .collect()
        })
        .collect();

    let (metrics, boot) = infer_native_ml_tree_with_bootstrap(
        &ry3_seqs,
        &labels,
        PhyloRunModel::GtrFreeRate2,
        50,
        2026,
    )
    .expect("Model-matched ML bootstrap on 8-taxon Carnivora RY3 panel must succeed");

    assert_eq!(boot.binary_ry_columns_detected, 100);
    assert_eq!(boot.four_state_columns_detected, 200);
    assert_eq!(boot.nontrivial_reference_splits, 5);
    assert!(boot.estimator_self_consistent);
    assert!(
        boot.identity_replicate_concordance_pct >= 95.0,
        "Expected >= 95.0% stationarity concordance, got {:.2}%",
        boot.identity_replicate_concordance_pct
    );
    assert!(
        boot.mean_split_support_pct >= 75.0 && metrics.mean_branch_support_pct == boot.mean_split_support_pct,
        "Expected strong synchronized UFBoot2 support on Carnivora RY3 panel (got {:.2}%)",
        boot.mean_split_support_pct
    );
}

#[test]
fn test_seqnet_coding_ortholog_qc_and_post_msa_row_guard() {
    use bio_suite::fetch::FastaRecord;
    use bio_suite::supermatrix_audit_util::{
        align_and_trim_codon_core_with_qc, align_and_trim_frame0_codon_core_with_qc,
        ncbi_table_to_paml_icode,
    };
    use bio_suite_seqnet::coding_ortholog::{
        audit_aligned_codon_rows, audit_coding_ortholog_cohort, fetch_refseq_gene_orthologs,
        CodingOrthologCandidate, CodingOrthologConfig, ProteinMotifRule,
    };

    // 20-codon diverse Table-11 block (~60 nt per block); repeated 7x = 140 codons (420 bp)
    let block = "GCTGAACTGAAAGGTGTTATCGACCGTACTTCTAACTTCTACCAGATGCCGTGGCACGAC";
    let make_cds = |i_motif: &str, c_motif: &str, stop_at_35: bool| -> String {
        let mut codons: Vec<&str> = Vec::with_capacity(140);
        let unit = [
            "GCT", "GAA", "CTG", "AAA", "GGT", "GTT", "ATC", "GAC", "CGT", "ACT",
            "TCT", "AAC", "TTC", "TAC", "CAG", "ATG", "CCG", "TGG", "CAC", "GAT",
        ];
        for k in 0..140 {
            codons.push(unit[k % unit.len()]);
        }
        let _ = block;
        if stop_at_35 {
            codons[35] = "TGA";
        }
        // Place 5-aa I-helix at codon 70..75
        let i_c = match i_motif {
            "AGHET" => ["GCT", "GGT", "CAT", "GAA", "ACT"],
            _ => ["GCT", "GGT", "CTG", "CAG", "ACT"],
        };
        for (o, c) in i_c.iter().enumerate() {
            codons[70 + o] = c;
        }
        // Place 10-aa Cys-loop at codon 110..120
        let c_c = match c_motif {
            "FGNGQRACIG" => ["TTC", "GGT", "AAC", "GGT", "CAG", "CGT", "GCT", "TGC", "ATC", "GGT"],
            _ => ["TTC", "GGT", "AAC", "GGT", "CAG", "CGT", "GCT", "GCT", "ATC", "GGT"],
        };
        for (o, c) in c_c.iter().enumerate() {
            codons[110 + o] = c;
        }
        codons.join("")
    };

    // Repeat AAA (Lys) and CCT (Pro) so motifs remain embedded while amino acid k-mer profile diverges sharply
    let make_outlier_cds = || -> String {
        let mut codons = vec!["AAA"; 140];
        let i_c = ["GCT", "GGT", "CAT", "GAA", "ACT"];
        for (o, c) in i_c.iter().enumerate() {
            codons[70 + o] = c;
        }
        let c_c = ["TTC", "GGT", "AAC", "GGT", "CAG", "CGT", "GCT", "TGC", "ATC", "GGT"];
        for (o, c) in c_c.iter().enumerate() {
            codons[110 + o] = c;
        }
        codons.join("")
    };

    let cfg = CodingOrthologConfig {
        cache_key: "test_cyp102a_qc".to_string(),
        genetic_code_table: 11,
        preserve_frame0: true,
        min_bp: 300,
        max_bp: 1000,
        max_domain_bp: Some(420),
        min_aa_len: 120,
        max_internal_stops: 0,
        target_n: 12,
        pin_species_substr: Some("megaterium".to_string()),
        blocked_organisms: vec!["anthracis".to_string()],
        required_motifs: vec![
            ProteinMotifRule {
                name: "I-helix".to_string(),
                pos_min: 55,
                pos_max: 85,
                patterns: vec!["AGHET".to_string(), "SGHET".to_string()],
            },
            ProteinMotifRule {
                name: "Cys-loop".to_string(),
                pos_min: 100,
                pos_max: 130,
                patterns: vec!["F..G...C.G".to_string()],
            },
        ],
        min_anchor_aa_kmer_sim: 0.35,
        min_row_coverage: 0.45,
        min_consensus_aa_identity: 0.18,
        min_consensus_dayhoff6_sim: 0.35,
    };

    let valid_cds = make_cds("AGHET", "FGNGQRACIG", false);
    let cands = vec![
        CodingOrthologCandidate {
            accession: "NC_VALID_1".to_string(),
            organism: "Priestia megaterium".to_string(),
            title: "CYP102A1".to_string(),
            subname: String::new(),
            sequence_raw: valid_cds.clone(),
        },
        CodingOrthologCandidate {
            accession: "NC_VALID_2".to_string(),
            organism: "Bacillus subtilis".to_string(),
            title: "CYP102A2".to_string(),
            subname: String::new(),
            sequence_raw: valid_cds.clone(),
        },
        CodingOrthologCandidate {
            accession: "NC_VALID_3".to_string(),
            organism: "Bacillus licheniformis".to_string(),
            title: "CYP102A3".to_string(),
            subname: String::new(),
            sequence_raw: valid_cds.clone(),
        },
        CodingOrthologCandidate {
            accession: "NC_VALID_4".to_string(),
            organism: "Bacillus velezensis".to_string(),
            title: "CYP102A4".to_string(),
            subname: String::new(),
            sequence_raw: valid_cds.clone(),
        },
        CodingOrthologCandidate {
            accession: "NC_BLOCKED".to_string(),
            organism: "Bacillus anthracis".to_string(),
            title: "CYP102A_blocked".to_string(),
            subname: String::new(),
            sequence_raw: valid_cds.clone(),
        },
        CodingOrthologCandidate {
            accession: "NC_STOP_TGA".to_string(),
            organism: "Geobacillus stearothermophilus".to_string(),
            title: "truncated pseudogene".to_string(),
            subname: String::new(),
            sequence_raw: make_cds("AGHET", "FGNGQRACIG", true),
        },
        CodingOrthologCandidate {
            accession: "NC_NO_MOTIF".to_string(),
            organism: "Paenibacillus polymyxa".to_string(),
            title: "missing heme motif".to_string(),
            subname: String::new(),
            sequence_raw: make_cds("AGLQT", "FGNGQRAAIG", false),
        },
        CodingOrthologCandidate {
            accession: "NC_OUTLIER".to_string(),
            organism: "Streptomyces coelicolor".to_string(),
            title: "frameshifted non-homolog".to_string(),
            subname: String::new(),
            sequence_raw: make_outlier_cds(),
        },
    ];

    let rep = audit_coding_ortholog_cohort(&cands, &cfg);
    assert_eq!(rep.accepted.len(), 4, "Expected 4 valid orthologs accepted");
    assert_eq!(rep.quarantined.len(), 4, "Expected all 4 invalid orthologs quarantined");
    assert!(rep.quarantined.iter().any(|q| q.accession == "NC_BLOCKED" && q.reason == "BLOCKED_ORGANISM"));
    assert!(rep.quarantined.iter().any(|q| q.accession == "NC_STOP_TGA" && q.reason == "TRUNCATED_OR_EARLY_STOP"));
    assert!(rep.quarantined.iter().any(|q| q.accession == "NC_NO_MOTIF" && q.reason == "MISSING_ORTHOLOG_MOTIF"));
    assert!(rep.quarantined.iter().any(|q| q.accession == "NC_OUTLIER" && q.reason == "NON_HOMOLOGOUS_PROTEIN_OUTLIER"));

    // Post-MSA row guard check: all-gap and unrelated tryptophan/poly-K rows must be quarantined
    let coherent_aa = "AELKGVIDRTSNFYQMPWHD".repeat(6);
    let coherent_nt = "GCT".repeat(120);
    let rows = vec![
        ("OK1".to_string(), coherent_aa.clone(), coherent_nt.clone()),
        ("OK2".to_string(), coherent_aa.clone(), coherent_nt.clone()),
        ("OK3".to_string(), coherent_aa.clone(), coherent_nt.clone()),
        ("OK4".to_string(), coherent_aa.clone(), coherent_nt.clone()),
        ("GAP_ROW".to_string(), "-".repeat(120), "---".repeat(120)),
        ("GARBAGE_ROW".to_string(), "W".repeat(120), "TGG".repeat(120)),
    ];
    let (kept_idx, q_rows) = audit_aligned_codon_rows(&rows, 0.45, 0.18, 0.35);
    assert_eq!(kept_idx, vec![0, 1, 2, 3]);
    assert_eq!(q_rows.len(), 2);

    // End-to-end align_and_trim_frame0_codon_core_with_qc drops non-homologous frameshift row
    let recs = vec![
        FastaRecord { id: "R1".into(), accession: "R1".into(), taxon: "Sp1".into(), common_name: String::new(), family: "F".into(), gene: "G".into(), description: String::new(), sequence: valid_cds.clone(), cached_from_disk: true },
        FastaRecord { id: "R2".into(), accession: "R2".into(), taxon: "Sp2".into(), common_name: String::new(), family: "F".into(), gene: "G".into(), description: String::new(), sequence: valid_cds.clone(), cached_from_disk: true },
        FastaRecord { id: "R3".into(), accession: "R3".into(), taxon: "Sp3".into(), common_name: String::new(), family: "F".into(), gene: "G".into(), description: String::new(), sequence: valid_cds.clone(), cached_from_disk: true },
        FastaRecord { id: "R4".into(), accession: "R4".into(), taxon: "Sp4".into(), common_name: String::new(), family: "F".into(), gene: "G".into(), description: String::new(), sequence: valid_cds, cached_from_disk: true },
        FastaRecord { id: "BAD".into(), accession: "BAD".into(), taxon: "Bad".into(), common_name: String::new(), family: "F".into(), gene: "G".into(), description: String::new(), sequence: "TTT".repeat(140), cached_from_disk: true },
    ];
    let (nt_kept, _, bp, _, bad_q) = align_and_trim_frame0_codon_core_with_qc(&recs, 11, 0.70, 100);
    assert!(bp >= 300);
    assert!(nt_kept.len() == 4 && nt_kept.iter().all(|(id, _)| id != "BAD"));
    assert!(bad_q.iter().any(|q| q.id == "BAD"));

    // Verify existing cached RefSeq CYP102A panel passes 12/12 with 0 quarantined if cache present
    let val_cache = Path::new(env!("CARGO_MANIFEST_DIR")).join("cache/val_set_070_gram_positive");
    if val_cache.join("gp_fam04_refseq_cyp102a_v2.json").is_file() {
        let mut cyp_cfg = CodingOrthologConfig {
            cache_key: "gp_fam04_refseq_cyp102a_v2".to_string(),
            genetic_code_table: 11,
            preserve_frame0: true,
            min_bp: 2900,
            max_bp: 3500,
            max_domain_bp: Some(1395),
            min_aa_len: 420,
            max_internal_stops: 0,
            target_n: 12,
            pin_species_substr: Some("megaterium".to_string()),
            blocked_organisms: vec!["anthracis".to_string(), "cereus".to_string()],
            required_motifs: vec![
                ProteinMotifRule {
                    name: "I-helix".to_string(),
                    pos_min: 225,
                    pos_max: 315,
                    patterns: vec!["AGHET".to_string(), "SGHET".to_string(), "AGFET".to_string()],
                },
                ProteinMotifRule {
                    name: "Cys-loop".to_string(),
                    pos_min: 350,
                    pos_max: 435,
                    patterns: vec!["F..G...C.G".to_string()],
                },
            ],
            ..CodingOrthologConfig::default()
        };
        cyp_cfg.target_n = 12;
        let cached_rep = fetch_refseq_gene_orthologs(&[], &cyp_cfg, &val_cache);
        assert_eq!(cached_rep.accepted.len(), 12, "Cached gp_fam04_refseq_cyp102a_v2 must yield 12/12");
        assert!(cached_rep.quarantined.is_empty(), "Cached gp_fam04_refseq_cyp102a_v2 must have 0 quarantined");
    }

    // Verify NCBI-to-PAML icode mapping (Table 1/5/11 -> 0 Universal; Table 2 -> 1 Vertebrate Mito)
    assert_eq!(ncbi_table_to_paml_icode(1), 0);
    assert_eq!(ncbi_table_to_paml_icode(11), 0);
    assert_eq!(ncbi_table_to_paml_icode(5), 0);
    assert_eq!(ncbi_table_to_paml_icode(2), 1);

    // Generalized Table-5 codon QC drops non-COX1 >2400bp contigs and non-homologous frameshifts
    let mito_cox1_unit = "GCTGAATTAGGTCAACCAGGTTCATTATTAGGTGATGATCAAATTTATAATGTAATTGTT";
    let mito_cox1 = mito_cox1_unit.repeat(7);
    let mito_recs = vec![
        FastaRecord { id: "M1".into(), accession: "M1".into(), taxon: "Pest1".into(), common_name: String::new(), family: "F".into(), gene: "COI".into(), description: String::new(), sequence: mito_cox1.clone(), cached_from_disk: true },
        FastaRecord { id: "M2".into(), accession: "M2".into(), taxon: "Pest2".into(), common_name: String::new(), family: "F".into(), gene: "COI".into(), description: String::new(), sequence: mito_cox1.clone(), cached_from_disk: true },
        FastaRecord { id: "M3".into(), accession: "M3".into(), taxon: "Pest3".into(), common_name: String::new(), family: "F".into(), gene: "COI".into(), description: String::new(), sequence: mito_cox1.clone(), cached_from_disk: true },
        FastaRecord { id: "M4".into(), accession: "M4".into(), taxon: "Pest4".into(), common_name: String::new(), family: "F".into(), gene: "COI".into(), description: String::new(), sequence: mito_cox1, cached_from_disk: true },
        FastaRecord { id: "M_BAD_SHORT".into(), accession: "MBS".into(), taxon: "BadS".into(), common_name: String::new(), family: "F".into(), gene: "COI".into(), description: String::new(), sequence: "TTT".repeat(140), cached_from_disk: true },
        FastaRecord { id: "M_BAD_CONTIG".into(), accession: "MBC".into(), taxon: "BadC".into(), common_name: String::new(), family: "F".into(), gene: "COI".into(), description: String::new(), sequence: "CCCCAAAAATTTTGGGG".repeat(180), cached_from_disk: true },
    ];
    let (m_nt, _, m_bp, _, m_q) = align_and_trim_codon_core_with_qc(&mito_recs, 5, false, 0.55, 80, 0.45, 0.18, 0.35);
    assert!(m_bp >= 300);
    assert_eq!(m_nt.len(), 4);
    assert!(m_nt.iter().all(|(id, _)| id != "M_BAD_SHORT" && id != "M_BAD_CONTIG"));
    assert!(m_q.iter().any(|q| q.id == "M_BAD_SHORT"));
    assert!(m_q.iter().any(|q| q.id == "M_BAD_CONTIG"));
}

#[test]
fn test_declarative_study_config_json_and_no_new_per_study_rs_files() {
    use bio_suite::supermatrix_audit_util::DeclarativeStudyConfig;

    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    for (fname, exp_loci, exp_pdb, exp_gcode) in [
        ("eukaryotic_nucleus_preprint32.json", 18usize, "1BYU", 1u8),
        ("gram_positive_preprint31.json", 10usize, "1BU7", 11u8),
        ("us_ag_pests_preprint30.json", 16usize, "1TQN", 5u8),
        ("china_australia_pests_preprint29.json", 16usize, "1TQN", 5u8),
        ("acanthochitonidae_preprint28.json", 17usize, "1TQN", 5u8),
        ("nematoda_preprint27.json", 16usize, "1TQN", 5u8),
        ("isoptera_preprint26.json", 16usize, "1TQN", 5u8),
        ("nudibranch_preprint25.json", 16usize, "1TQN", 5u8),
        ("invasive_snails_preprint24.json", 16usize, "1TQN", 5u8),
        ("bemisia_preprint23.json", 16usize, "1TQN", 5u8),
        ("anopheles_preprint22.json", 17usize, "1TQN", 5u8),
        ("basal_metazoa_preprint05.json", 20usize, "1LI1", 1u8),
    ] {
        let cfg_path = manifest_dir.join("benchmarks/studies").join(fname);
        let cfg = DeclarativeStudyConfig::load(&cfg_path)
            .unwrap_or_else(|e| panic!("{} must load into DeclarativeStudyConfig: {}", fname, e));
        assert_eq!(cfg.species_slots.len(), 64, "{} species_slots != 64", fname);
        assert_eq!(cfg.locus_defs.len(), exp_loci, "{} locus_defs", fname);
        assert_eq!(cfg.functional_families.len(), 4, "{} functional_families != 4", fname);
        assert_eq!(cfg.skyline_cohorts.len(), 12, "{} skyline_cohorts != 12", fname);
        assert_eq!(cfg.skyline_epochs.len(), 6, "{} skyline_epochs != 6", fname);
        assert_eq!(cfg.haplotype_panels.len(), 2, "{} haplotype_panels != 2", fname);
        assert_eq!(cfg.calibrated_nodes.len(), 8, "{} calibrated_nodes != 8", fname);
        assert_eq!(cfg.abba_baba_contrasts.len(), 4, "{} abba_baba_contrasts != 4", fname);
        assert_eq!(cfg.pdb_template.pdb_id, exp_pdb, "{} pdb_id", fname);
        assert_eq!(cfg.genetic_code_table, exp_gcode, "{} genetic_code_table", fname);
    }

    let cfg_path = manifest_dir.join("benchmarks/studies/gram_positive_preprint31.json");
    let cfg = DeclarativeStudyConfig::load(&cfg_path)
        .expect("gram_positive_preprint31.json must load into DeclarativeStudyConfig");

    let json_roundtrip = serde_json::to_string(&cfg).expect("serialize DeclarativeStudyConfig");
    let cfg_rt: DeclarativeStudyConfig =
        serde_json::from_str(&json_roundtrip).expect("deserialize roundtrip DeclarativeStudyConfig");
    assert_eq!(cfg_rt.study_slug, cfg.study_slug);
    assert_eq!(cfg_rt.run_id, cfg.run_id);

    // Freeze the 17 legacy supervisor_*.rs files: all new studies must use JSON configs, never new .rs files.
    let mut supervisor_files: Vec<String> = fs::read_dir(manifest_dir.join("src"))
        .expect("read src dir")
        .filter_map(|e| e.ok())
        .filter_map(|e| e.file_name().into_string().ok())
        .filter(|name| name.starts_with("supervisor_") && name.ends_with(".rs"))
        .collect();
    supervisor_files.sort();
    let expected_frozen = vec![
        "supervisor_acanthochitonidae.rs",
        "supervisor_basal_metazoa.rs",
        "supervisor_basal_metazoa_figures.rs",
        "supervisor_basal_metazoa_manuscript.rs",
        "supervisor_china_aus_pests.rs",
        "supervisor_china_aus_pests_manuscript.rs",
        "supervisor_figures.rs",
        "supervisor_gram_positive.rs",
        "supervisor_gram_positive_manuscript.rs",
        "supervisor_isoptera.rs",
        "supervisor_manuscript.rs",
        "supervisor_nematoda.rs",
        "supervisor_nudibranchs.rs",
        "supervisor_snails.rs",
        "supervisor_study.rs",
        "supervisor_us_ag_pests.rs",
        "supervisor_us_ag_pests_manuscript.rs",
    ];
    assert_eq!(
        supervisor_files, expected_frozen,
        "Per-study Rust files are frozen at 17 legacy files; author new studies as JSON configs in benchmarks/studies/<study>.json"
    );
}

#[test]
fn test_bootstrap_control_no_rell_collapse_on_unaligned_noise() {
    use bio_suite::prng_util::DeterministicLcg64;
    use bio_suite::tree_builder::{infer_native_ml_tree_with_bootstrap_control, PhyloRunModel};

    let bases = [b'A', b'C', b'G', b'T'];
    let mut rng = DeterministicLcg64::new(0x5241_5730_5F54_4553);
    let mut labels = Vec::with_capacity(12);
    let mut seqs = Vec::with_capacity(12);
    for t in 0..12 {
        labels.push(format!("NOISE_{:02}", t + 1));
        let mut s = String::with_capacity(1800);
        for _ in 0..1800 {
            s.push(bases[rng.gen_range_usize(4)] as char);
        }
        seqs.push(s);
    }
    let (_, rep) = infer_native_ml_tree_with_bootstrap_control(
        &seqs,
        &labels,
        PhyloRunModel::HkyInvGamma4,
        40,
        20261122,
    )
    .expect("bootstrap control on unaligned noise must succeed");
    assert_eq!(rep.b_replicates, 40);
    assert!(
        rep.mean_split_support_pct < 55.0,
        "non-gating control must not suffer RELL collapse on unaligned noise (got {:.2}%)",
        rep.mean_split_support_pct
    );
}

mod empirical_anti_shortcut_tests {
    use bio_suite::coalescent::run_native_quartet_coalescent;
    use bio_suite::empirical_inference_util::{
        analyze_pdb_file_empirically, compute_bic_posterior_weights, compute_kh_paired_sites_z,
        run_empirical_bayesian_skyline_mcmc, solve_empirical_lsd2_from_ml_tree,
        EmpiricalCalibrationSpec,
    };
    use std::fs;
    use std::path::PathBuf;

    #[test]
    fn test_no_hardcoded_shortcuts_in_generalized_empirical_engines() {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let src_dir = root.join("src");

        // Guard Blackbird < 384 KB source file ceiling across all src/*.rs
        fn check_sizes(dir: &std::path::Path) {
            if let Ok(entries) = fs::read_dir(dir) {
                for ent in entries.flatten() {
                    let p = ent.path();
                    if p.is_dir() {
                        check_sizes(&p);
                    } else if p.extension().is_some_and(|e| e == "rs") {
                        let sz = fs::metadata(&p).map(|m| m.len()).unwrap_or(0);
                        assert!(
                            sz < 393_216,
                            "Rust source {} exceeds 384 KB Blackbird limit ({} bytes)",
                            p.display(),
                            sz
                        );
                    }
                }
            }
        }
        check_sizes(&src_dir);

        let audit_util = fs::read_to_string(src_dir.join("supermatrix_audit_util.rs")).unwrap();
        let emp_util = fs::read_to_string(src_dir.join("empirical_inference_util.rs")).unwrap();
        let study_full = fs::read_to_string(src_dir.join("supervisor_study.rs")).unwrap();
        let study_gen_end = study_full.find("evaluate_population_cohort_4regimes").unwrap_or(0) + 12000;
        let study_slice = &study_full[..study_gen_end.min(study_full.len())];

        let figs_full = fs::read_to_string(src_dir.join("supervisor_figures.rs")).unwrap();
        let figs_start = figs_full.find("fn collect_contiguous_sector_runs").unwrap();
        let figs_end = figs_full.find("fn generate_acanthochitonidae_locked_figures").unwrap_or(figs_full.len());
        let figs_slice = &figs_full[figs_start..figs_end];

        let ms_full = fs::read_to_string(src_dir.join("supervisor_manuscript.rs")).unwrap();
        let ms_start = ms_full.find("pub fn export_declarative_manuscript_and_staged_bundle").unwrap();
        let ms_slice = &ms_full[ms_start..];

        let banned_tokens = [
            "n.age_ma * 0.89",
            "n.age_ma * 1.11",
            "0.32 * (",
            "-9.42",
            "1094.2",
            "ACGTACGT",
            "append_rotated_acgt",
            "aa_volume_proxy",
            "build_unified_rooted_tree",
            "normalized_quartet_score - 0.14",
            ".clamp(0.80, 5.50)",
            "+4.86 SD",
            "+3.72 SD",
            "+2.18 SD",
            "-3.94 SD",
            "11,436,768 quartets",
            "\"112ka\"",
            "ESS > 480",
            "R-hat = 1.000",
            "N = 432",
            "pop_s * 1.06",
            "ml_s * 0.85",
            "unwrap_or(92.4)",
            "November 21, 2026",
            "&[0.15, 0.49, 1.03, 2.33]",
            "4.98,",
        ];

        for (fname, slice) in [
            ("supermatrix_audit_util.rs", audit_util.as_str()),
            ("empirical_inference_util.rs", emp_util.as_str()),
            ("supervisor_study.rs(generalized)", study_slice),
            ("supervisor_figures.rs(generalized)", figs_slice),
            ("supervisor_manuscript.rs(declarative)", ms_slice),
        ] {
            for tok in banned_tokens {
                assert!(
                    !slice.contains(tok),
                    "Banned hardcoded shortcut token `{}` reappeared inside {}",
                    tok,
                    fname
                );
            }
        }

        let mcmc_skyline_src = fs::read_to_string(src_dir.join("mcmc_skyline.rs")).unwrap();
        for banned_floor in [
            ".max(12_840)",
            ".max(9_420)",
            ".max(9_680)",
            ".max(9_820)",
            ".max(1_850_000)",
            ".max(1_420)",
            ".max(3_680)",
        ] {
            assert!(
                !mcmc_skyline_src.contains(banned_floor),
                "Banned artificial timing floor `{}` reappeared inside mcmc_skyline.rs",
                banned_floor
            );
        }

        let anopheles_src = fs::read_to_string(src_dir.join("anopheles_biosuite.rs")).unwrap();
        for banned_floor in [".max(940)", ".max(260)", ".max(620)"] {
            assert!(
                !anopheles_src.contains(banned_floor),
                "Banned artificial timing floor `{}` reappeared inside anopheles_biosuite.rs",
                banned_floor
            );
        }
    }

    #[test]
    fn test_empirical_algorithms_respond_to_data_perturbation() {
        // 1. LSD2 Sensitivity & Non-Monophyletic Zeroing
        let ml_a = "((A_T001:0.10,A_T002:0.10)100:0.15,((C_T017:0.06,C_T018:0.06)100:0.12,(D_T025:0.04,D_T026:0.04)100:0.08)100:0.15);";
        let ml_b = "((A_T001:0.10,A_T002:0.10)100:0.15,((C_T017:0.06,C_T018:0.06)100:0.12,(D_T025:0.02,D_T026:0.02)100:0.26)100:0.15);";
        let tips: Vec<String> = ["A_T001", "A_T002", "C_T017", "C_T018", "D_T025", "D_T026"]
            .into_iter()
            .map(str::to_string)
            .collect();
        let full_mask = (1u64 << 6) - 1;
        let split_map = vec![(vec![(1u64 << 4) | (1u64 << 5)], 97.5)];
        let gts = vec![ml_a.to_string(), ml_a.to_string()];
        let (msc_res, _, _) = run_native_quartet_coalescent("lsd2_pert", &gts).unwrap();

        let cal_specs = vec![
            EmpiricalCalibrationSpec {
                node_id: "Node I".to_string(),
                clade_name: "Root".to_string(),
                sector_code: "[A]-[D]".to_string(),
                tip_prefixes: Vec::new(),
                fossil_or_biogeo_anchor: "Root bound".to_string(),
                ci_min_ma: 1000.0,
                ci_max_ma: 1200.0,
                trophic_guild: String::new(),
                key_enzymes_or_loci: String::new(),
                ecosystem_role: String::new(),
            },
            EmpiricalCalibrationSpec {
                node_id: "Node D".to_string(),
                clade_name: "Ingroup D".to_string(),
                sector_code: "[D]".to_string(),
                tip_prefixes: Vec::new(),
                fossil_or_biogeo_anchor: "Loose wide prior".to_string(),
                ci_min_ma: 80.0,
                ci_max_ma: 950.0,
                trophic_guild: String::new(),
                key_enzymes_or_loci: String::new(),
                ecosystem_role: String::new(),
            },
            EmpiricalCalibrationSpec {
                node_id: "Node Absent".to_string(),
                clade_name: "Absent Sector E".to_string(),
                sector_code: "[E]".to_string(),
                tip_prefixes: Vec::new(),
                fossil_or_biogeo_anchor: "None".to_string(),
                ci_min_ma: 200.0,
                ci_max_ma: 400.0,
                trophic_guild: String::new(),
                key_enzymes_or_loci: String::new(),
                ecosystem_role: String::new(),
            },
        ];

        let lsd_a = solve_empirical_lsd2_from_ml_tree(
            ml_a,
            &tips,
            full_mask,
            &split_map,
            &msc_res,
            &cal_specs,
            1200,
            42,
        )
        .expect("LSD2 A");
        let lsd_b = solve_empirical_lsd2_from_ml_tree(
            ml_b,
            &tips,
            full_mask,
            &split_map,
            &msc_res,
            &cal_specs,
            1200,
            42,
        )
        .expect("LSD2 B");

        assert!(lsd_a.timetree_newick.contains("[&age_ma="));
        assert_eq!(lsd_a.calibrated_nodes[1].ufboot_pct, 97.5);
        assert!(lsd_a.calibrated_nodes[1].quartet_lpp > 0.0);
        // Absent clade must strictly report 0.0 instead of fabricating an age or support
        assert_eq!(lsd_a.calibrated_nodes[2].age_ma, 0.0);
        assert_eq!(lsd_a.calibrated_nodes[2].ci_min_ma, 0.0);
        assert_eq!(lsd_a.calibrated_nodes[2].ci_max_ma, 0.0);
        assert_eq!(lsd_a.calibrated_nodes[2].ufboot_pct, 0.0);
        assert_eq!(lsd_a.calibrated_nodes[2].quartet_lpp, 0.0);
        // Doubling/shifting branch lengths must noticeably shift solved node age
        let delta_age = (lsd_a.calibrated_nodes[1].age_ma - lsd_b.calibrated_nodes[1].age_ma).abs();
        assert!(
            delta_age > 1.0,
            "LSD2 age failed to respond to ML branch perturbation (A={}, B={})",
            lsd_a.calibrated_nodes[1].age_ma,
            lsd_b.calibrated_nodes[1].age_ma
        );

        // 2. Bayesian Skyline MCMC Sensitivity
        let base_seq: Vec<u8> = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
            .iter()
            .cycle()
            .take(180)
            .copied()
            .collect();
        let mut lbls = Vec::with_capacity(10);
        let mut star_seqs = Vec::with_capacity(10);
        let mut deep_seqs = Vec::with_capacity(10);
        for i in 0..10 {
            lbls.push(format!("ISO_{:02}", i + 1));
            let mut s1 = base_seq.clone();
            s1[i * 3] = match s1[i * 3] {
                b'A' => b'G',
                b'G' => b'A',
                b'C' => b'T',
                _ => b'C',
            };
            star_seqs.push(s1);

            let mut s2 = base_seq.clone();
            let offset = if i < 5 { 0 } else { 15 };
            for p in 0..14 {
                let idx = (offset + p * 7) % 180;
                s2[idx] = match s2[idx] {
                    b'A' => b'T',
                    b'T' => b'A',
                    b'C' => b'G',
                    _ => b'C',
                };
            }
            s2[i * 2] = if i % 2 == 0 { b'G' } else { b'C' };
            deep_seqs.push(s2);
        }
        let sky_star = run_empirical_bayesian_skyline_mcmc(&lbls, &star_seqs, 0.01, 1.0, 0.005, 4_000, 101);
        let sky_deep = run_empirical_bayesian_skyline_mcmc(&lbls, &deep_seqs, 0.01, 1.0, 0.04, 4_000, 101);
        assert!(
            (sky_star.marginal_ln_l - sky_deep.marginal_ln_l).abs() > 5.0,
            "Skyline lnL must change sharply between star and deep alignments ({} vs {})",
            sky_star.marginal_ln_l,
            sky_deep.marginal_ln_l
        );
        let max_ne_ratio = sky_star
            .epochs
            .iter()
            .zip(sky_deep.epochs.iter())
            .map(|(a, b)| (a.2.max(b.2)) / a.2.min(b.2).max(1.0))
            .fold(1.0_f64, f64::max);
        assert!(
            max_ne_ratio > 1.5,
            "Skyline Ne profile failed to scale with underlying sequence divergence (ratio={})",
            max_ne_ratio
        );

        // 3. 3D Pocket Volume + AutoDock Vina Sensitivity
        let cache_1byu = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("cache/pdb/1BYU.pdb");
        let tmp_dir = std::env::temp_dir();
        let orig_pdb_path = if cache_1byu.is_file() {
            cache_1byu.clone()
        } else {
            let p = tmp_dir.join("crv_test_pocket_orig.pdb");
            let mut pdb = String::new();
            let mut ser = 1usize;
            for &(x, y, z, res) in &[
                (10.0, 0.0, 0.0, "GLY"), (8.5, 4.0, 2.0, "ALA"), (4.0, 6.5, -2.5, "VAL"),
                (-3.5, 6.0, 3.0, "SER"), (-6.5, 2.0, -3.0, "THR"), (-6.0, -4.5, 2.0, "LYS"),
                (-2.0, -6.5, -2.5, "ASP"), (4.0, -6.0, 2.5, "ASN"), (7.5, -3.5, -2.0, "GLU"),
                (0.0, 0.0, 6.0, "PHE"), (0.0, 0.0, -6.0, "LEU"), (5.0, 5.0, 5.0, "ILE"),
            ] {
                pdb.push_str(&format!(
                    "ATOM  {:>5}  CA  {} A{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00 18.00           C\n",
                    ser, res, ser, x, y, z
                ));
                pdb.push_str(&format!(
                    "ATOM  {:>5}  N   {} A{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00 18.00           N\n",
                    ser + 100, res, ser, x + 0.8, y - 0.6, z + 0.4
                ));
                ser += 1;
            }
            pdb.push_str("HETATM  901  PG  GDP A 220       0.200   0.100   0.000  1.00 15.00           P\n");
            pdb.push_str("HETATM  902  O1B GDP A 220       1.600   0.400   0.200  1.00 15.00           O\n");
            pdb.push_str("HETATM  903  C1' GDP A 220      -1.200  -0.300   0.400  1.00 15.00           C\n");
            pdb.push_str("END\n");
            fs::write(&p, pdb).unwrap();
            p
        };
        let shifted_path = tmp_dir.join("crv_test_pocket_shifted15A.pdb");
        let raw_pdb = fs::read_to_string(&orig_pdb_path).unwrap();
        let mut shifted_pdb = String::with_capacity(raw_pdb.len());
        let mut cas = Vec::new();
        for line in raw_pdb.lines() {
            if line.starts_with("ATOM  ") && line.len() >= 54 && line.get(12..16).map(|s| s.trim()) == Some("CA") {
                if let (Ok(r), Ok(x), Ok(y), Ok(z)) = (
                    line[22..26].trim().parse::<usize>(),
                    line[30..38].trim().parse::<f64>(),
                    line[38..46].trim().parse::<f64>(),
                    line[46..54].trim().parse::<f64>(),
                ) {
                    cas.push((r, x, y, z));
                }
            }
            if line.starts_with("HETATM") && line.len() >= 54 && line.get(17..20).map(|s| s.trim()) == Some("GDP") {
                if let (Ok(x), Ok(y), Ok(z)) = (
                    line[30..38].trim().parse::<f64>(),
                    line[38..46].trim().parse::<f64>(),
                    line[46..54].trim().parse::<f64>(),
                ) {
                    shifted_pdb.push_str(&format!(
                        "{}{:>8.3}{:>8.3}{:>8.3}{}\n",
                        &line[..30],
                        x + 15.0,
                        y + 15.0,
                        z + 15.0,
                        &line[54..]
                    ));
                    continue;
                }
            }
            shifted_pdb.push_str(line);
            shifted_pdb.push('\n');
        }
        fs::write(&shifted_path, &shifted_pdb).unwrap();

        let emp_orig = analyze_pdb_file_empirically(&orig_pdb_path, 'A', "GDP", None, &cas).unwrap();
        let emp_shift = analyze_pdb_file_empirically(&shifted_path, 'A', "GDP", None, &cas).unwrap();
        let _ = fs::remove_file(&shifted_path);
        assert!(emp_orig.pocket_volume_a3 > 10.0, "Expected non-trivial voxel cavity (got {})", emp_orig.pocket_volume_a3);
        assert!(emp_orig.docking_affinity_kcal_mol < -1.0, "Expected bound Vina affinity < -1 kcal/mol (got {})", emp_orig.docking_affinity_kcal_mol);
        assert!(
            emp_shift.docking_affinity_kcal_mol > emp_orig.docking_affinity_kcal_mol + 0.5,
            "Displaced +15A ligand must weaken Vina affinity (orig={}, shifted={})",
            emp_orig.docking_affinity_kcal_mol,
            emp_shift.docking_affinity_kcal_mol
        );

        // 4. BIC Posterior Weights & KH Paired-Sites Z
        let w = compute_bic_posterior_weights(&[1000.0, 1004.0, 1040.0]);
        let w_sum: f64 = w.iter().sum();
        assert!((w_sum - 1.0).abs() < 1e-2 && w[0] > w[1] && w[1] > w[2]);
        let ell_a: Vec<f64> = (0..40).map(|i| -1.2 - (i % 3) as f64 * 0.05).collect();
        let ell_b: Vec<f64> = (0..40).map(|i| -2.0 - (i % 5) as f64 * 0.1).collect();
        let wts = vec![2.0; 40];
        let z_ab = compute_kh_paired_sites_z(&ell_a, &wts, &ell_b, &wts);
        let z_ba = compute_kh_paired_sites_z(&ell_b, &wts, &ell_a, &wts);
        assert!((z_ab + z_ba).abs() < 1e-6 && z_ab > 2.0, "KH Z must be antisymmetric and significant (got {} / {})", z_ab, z_ba);
    }

    #[test]
    fn test_dynamic_coverage_gate_relaxes_evenly_by_sector() {
        use bio_suite::supervisor_study::{
            passes_supermatrix_phylogram_coverage_gate_with_thresholds,
            resolve_dynamic_supermatrix_coverage_gate,
        };

        let total_cols = 10_000usize;
        let n_loci = 18usize;

        // 1. All 7 sectors already have >= 4 OTUs at >= 80% -> retains 0.80 strict gate
        let mut full_stats = Vec::new();
        for sec in 0..7 {
            for _ in 0..5 {
                full_stats.push((sec, 8_500usize, 15usize));
            }
        }
        let res_full = resolve_dynamic_supermatrix_coverage_gate(
            &full_stats, total_cols, n_loci, true, 0.80, 0.55, 4,
        );
        assert!((res_full.min_coverage_ratio - 0.80).abs() < 1e-9);
        assert_eq!(res_full.min_loci, 8);
        assert!(!res_full.is_dynamic_relaxed);

        // 2. Sectors [D] (3) and [G] (6) only reach >= 4 OTUs around 62% -> steps down cleanly to 0.60
        let mut sparse_stats = Vec::new();
        for sec in 0..7 {
            if sec == 3 || sec == 6 {
                sparse_stats.push((sec, 8_400, 14));
                for _ in 0..3 {
                    sparse_stats.push((sec, 6_200, 11));
                }
                sparse_stats.push((sec, 4_500, 8));
            } else {
                for _ in 0..5 {
                    sparse_stats.push((sec, 8_600, 16));
                }
            }
        }
        let res_relaxed = resolve_dynamic_supermatrix_coverage_gate(
            &sparse_stats, total_cols, n_loci, true, 0.80, 0.55, 4,
        );
        assert!(
            (res_relaxed.min_coverage_ratio - 0.60).abs() < 1e-9,
            "Expected relaxed 0.60 ratio to cover >=4 OTUs in [D] & [G], got {}",
            res_relaxed.min_coverage_ratio
        );
        assert!(res_relaxed.is_dynamic_relaxed);
        assert_eq!(res_relaxed.min_loci, 8);
        let retained_d = sparse_stats
            .iter()
            .filter(|&&(sec, vb, lp)| {
                sec == 3
                    && passes_supermatrix_phylogram_coverage_gate_with_thresholds(
                        vb,
                        total_cols,
                        lp,
                        res_relaxed.min_coverage_ratio,
                        res_relaxed.min_loci,
                    )
            })
            .count();
        assert_eq!(retained_d, 4);

        // 3. Respects floor_ratio when a sector is universally sparse below 0.55
        let mut ultra_sparse = sparse_stats.clone();
        for item in &mut ultra_sparse {
            if item.0 == 3 && item.1 == 6_200 {
                item.1 = 5_600;
            }
        }
        let res_floor = resolve_dynamic_supermatrix_coverage_gate(
            &ultra_sparse, total_cols, n_loci, true, 0.80, 0.55, 4,
        );
        assert!(
            (res_floor.min_coverage_ratio - 0.55).abs() < 1e-9,
            "Expected clamp to floor 0.55, got {}",
            res_floor.min_coverage_ratio
        );
    }

    #[test]
    fn test_lsd2_preserves_strict_parent_child_monotonicity_under_tight_upper_bounds() {
        use bio_suite::empirical_inference_util::solve_empirical_lsd2_core;
        use bio_suite::timetree_calibration::TimeTreeNodeSpec;

        // Node 0 (Root) -> {Node 1 (Parent capped at 2690 Ma), Node 3 (Tip)}
        // Node 1 -> {Node 2 (unconstrained internal child with very long descendant branches), Node 4 (Tip)}
        // Node 2 -> {Node 5 (Tip), Node 6 (Tip)}
        let nodes = vec![
            TimeTreeNodeSpec {
                node_id: 0,
                parent_id: None,
                label: "Root".to_string(),
                is_tip: false,
                is_extinct_fossil_tip: false,
                branch_length_subs: 0.0,
                tip_age_ma: None,
                cal_min_ma: Some(3000.0),
                cal_max_ma: Some(3800.0),
                focal_node_tag: None,
            },
            TimeTreeNodeSpec {
                node_id: 1,
                parent_id: Some(0),
                label: "ParentTightMax".to_string(),
                is_tip: false,
                is_extinct_fossil_tip: false,
                branch_length_subs: 0.05,
                tip_age_ma: None,
                cal_min_ma: Some(2200.0),
                cal_max_ma: Some(2690.0),
                focal_node_tag: None,
            },
            TimeTreeNodeSpec {
                node_id: 2,
                parent_id: Some(1),
                label: "LongBranchChild".to_string(),
                is_tip: false,
                is_extinct_fossil_tip: false,
                branch_length_subs: 0.01,
                tip_age_ma: None,
                cal_min_ma: None,
                cal_max_ma: None,
                focal_node_tag: None,
            },
            TimeTreeNodeSpec {
                node_id: 3,
                parent_id: Some(0),
                label: "OutTip".to_string(),
                is_tip: true,
                is_extinct_fossil_tip: false,
                branch_length_subs: 0.40,
                tip_age_ma: Some(0.0),
                cal_min_ma: None,
                cal_max_ma: None,
                focal_node_tag: None,
            },
            TimeTreeNodeSpec {
                node_id: 4,
                parent_id: Some(1),
                label: "SisterTip".to_string(),
                is_tip: true,
                is_extinct_fossil_tip: false,
                branch_length_subs: 0.30,
                tip_age_ma: Some(0.0),
                cal_min_ma: None,
                cal_max_ma: None,
                focal_node_tag: None,
            },
            TimeTreeNodeSpec {
                node_id: 5,
                parent_id: Some(2),
                label: "LongTipA".to_string(),
                is_tip: true,
                is_extinct_fossil_tip: false,
                branch_length_subs: 0.55,
                tip_age_ma: Some(0.0),
                cal_min_ma: None,
                cal_max_ma: None,
                focal_node_tag: None,
            },
            TimeTreeNodeSpec {
                node_id: 6,
                parent_id: Some(2),
                label: "LongTipB".to_string(),
                is_tip: true,
                is_extinct_fossil_tip: false,
                branch_length_subs: 0.55,
                tip_age_ma: Some(0.0),
                cal_min_ma: None,
                cal_max_ma: None,
                focal_node_tag: None,
            },
        ];

        let (solved_t, _cis, rate, _q) = solve_empirical_lsd2_core(&nodes, 5000, 42);
        assert!(rate > 0.0 && rate.is_finite());
        assert!(solved_t[1] <= 2690.0 + 1e-6, "Node 1 must respect cal_max 2690: {}", solved_t[1]);
        for node in &nodes {
            if let Some(pid) = node.parent_id {
                assert!(
                    solved_t[pid] > solved_t[node.node_id],
                    "Strict parent > child monotonicity violated: node {} ({}) vs child {} ({})",
                    pid,
                    solved_t[pid],
                    node.node_id,
                    solved_t[node.node_id]
                );
            }
        }
    }
}


