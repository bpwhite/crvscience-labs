//! `bio-suite` — CRV Science Labs High-Throughput Bioinformatics, Molecular Systematics
//! & Information-Theoretic Genome Complexity Engine
//!
//! Released under GNU General Public License v3.0 (GPL-3.0-or-later).
//! Repository: https://github.com/bpwhite/crvscience-labs.git

use bio_suite::agent::{
    run_chr1_complexity_agent, run_coi_barcode_cluster_agent_with_mode,
};
use bio_suite::amr_hgt::run_amr_hgt_calibration_benchmark;
use bio_suite::avian_systematics::run_avian_tiered_systematics_sweep;
use bio_suite::bgc_mine::run_bgc_calibration_benchmark;
use bio_suite::clonal_synth_calibration::run_clonal_synth_vs_scipy_calibration;
use bio_suite::coalescent::run_coalescent_native_vs_astral_calibration;
use bio_suite::codon_selection::run_codon_selection_native_vs_paml_calibration;
use bio_suite::dag_pipelines::{
    run_chromatin_condensate_scaffolder_agent, run_coevo_epistatic_escape_agent,
    run_master_dag_chain_agent, run_somatic_clonal_synthetic_lethal_agent,
};
use bio_suite::disease_panels::{build_disease_genes_atlas, run_apoe_tiered_discovery_sweep};
use bio_suite::fetch::PoliteCacheManager;
use bio_suite::gene_caller::run_gene_caller_native_vs_prodigal_calibration;
use bio_suite::genome_target::BIO_SUITE_TOOL_REGISTRY;
use bio_suite::haplotype::run_haplotype_cluster_pipeline_agent;
use bio_suite::modeltest::{
    run_factorial_msa_x_modeltest_calibration, run_modeltest_native_vs_jmodeltest_calibration,
    run_modeltest_with_mode, ModelTestMode,
};
use bio_suite::msa_calibration::run_msa_native_vs_mafft_calibration;
use bio_suite::na_structure_1d2d3d::run_na_1d2d3d_exploratory_calibration_sweep;
use bio_suite::nonb_dna::run_nonb_dna_vs_g4hunter_zhunt_calibration;
use bio_suite::orthology::run_orthology_native_vs_diamond_calibration;
use bio_suite::pangenome::run_pangenome_calibration_benchmark;
use bio_suite::phylo::{align_sequences_with_mode, AlignerMode};
use bio_suite::pipelines::{
    run_bgc_amp_miner_agent, run_hgt_resistome_scanner_agent, run_pangenome_cryptic_target_agent,
    run_pseudogene_trait_loss_agent,
};
use bio_suite::pocket_dock::run_pocket_dock_calibration_benchmark;
use bio_suite::profile_hmm::run_profile_hmm_vs_hmmer3_calibration;
use bio_suite::pseudo_loss_calibration::run_pseudo_loss_calibration;
use bio_suite::rna_fold::run_rna_fold_vs_viennarna_benchmark;
use bio_suite::runie::run_runie_master_pipeline_inference;
use bio_suite::sweep::{
    run_multi_pipeline_discovery_sweep, run_procedural_burnin_calibration_sweep,
};
use bio_suite::synteny::run_synteny_native_vs_minimap2_calibration;
use bio_suite::tree_builder::run_tree_builder_native_vs_iqtree_calibration;
use serde::Serialize;
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

fn emit_json_output<T: Serialize>(report: &T, duration_ms: u128, out_file: Option<&String>) {
    let json = serde_json::to_string_pretty(report).expect("Failed to serialize pipeline report");
    if let Some(path) = out_file {
        if let Some(parent) = PathBuf::from(path).parent() {
            let _ = fs::create_dir_all(parent);
        }
        fs::write(path, &json).expect("Failed to write output file");
        eprintln!(
            "[bio-suite] Completed in {} ms. Output written to {}",
            duration_ms, path
        );
    } else {
        println!("{}", json);
    }
}

fn dispatch_pipeline(
    workflow: &str,
    cache: &PoliteCacheManager,
    out_file: Option<&String>,
    aligner_mode: AlignerMode,
) {
    match workflow {
        "coi-barcode-cluster" => {
            eprintln!(
                "[bio-suite] Launching Pipeline 1: coi-barcode-cluster (500 Rayon bootstrap replicates, aligner={:?})...",
                aligner_mode
            );
            match run_coi_barcode_cluster_agent_with_mode(cache, 500, aligner_mode) {
                Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
                Err(e) => {
                    eprintln!("[bio-suite] Pipeline failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "modeltest" => {
            let modeltest_mode = env::args()
                .position(|a| a == "--modeltest")
                .and_then(|idx| env::args().nth(idx + 1))
                .map(|s| ModelTestMode::from_cli_arg(&s))
                .unwrap_or(ModelTestMode::Native);
            eprintln!(
                "[bio-suite] Launching Nucleotide Substitution Model Selection (engine={:?}, aligner={:?})...",
                modeltest_mode, aligner_mode
            );
            match cache.fetch_carnivora_coi_panel() {
                Ok(records) => {
                    let (_, aligned, _) = align_sequences_with_mode(&records, aligner_mode);
                    let res = run_modeltest_with_mode(&aligned, modeltest_mode);
                    eprintln!("\n{}", res.ascii_model_ranking_table);
                    emit_json_output(&res, res.duration_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] ModelTest failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "tree-compare" | "tree-calibration" | "iqtree-compare" => {
            eprintln!("[bio-suite] Launching 54-Taxon (106-Node) Jawed Vertebrate Tree Builder Calibration: Native ML Tree Builder (BIONJ+NNI) vs. IQ-TREE v2.4.0 (Using Native ModelTest Only)...");
            match run_tree_builder_native_vs_iqtree_calibration(cache) {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_summary_table);
                    eprintln!("\n{}", report.ascii_vertebrate_cladogram);
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Tree builder calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "coalescent-compare" | "astral-compare" | "coalescent-calibration" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching Tool Pair 4 Calibration: Native Rust Quartet Multispecies Coalescent (MSC + LPP + ABBA-BABA) vs. ASTRAL v5.7.8...");
            match run_coalescent_native_vs_astral_calibration(&ws_root, cache) {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_master_summary_table);
                    for c in &report.cohorts {
                        eprintln!("\n{}", c.ascii_cohort_comparison_table);
                    }
                    eprintln!("\n{}", report.ascii_mammalian_species_cladogram);
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Coalescent calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "selection-compare" | "codeml-compare" | "paml-compare" | "selection-calibration" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching Tool Pair 5 Calibration: Native Rust 61-State Goldman-Yang Codon CTMC + 4D Grid BEB vs. PAML v4.10.10 (codeml)...");
            match run_codon_selection_native_vs_paml_calibration(cache, &ws_root) {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_summary_table);
                    eprintln!("\n{}", report.ascii_site_posterior_table);
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Codon selection calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "orthology-compare" | "diamond-compare" | "orthology-calibration" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching Tool Pair 6 Calibration: Native Rust 4-Stage Smith-Waterman BLOSUM62, Emms-Kelly Length Normalization & RBH Orthology Engine vs. DIAMOND v2.2.8...");
            match run_orthology_native_vs_diamond_calibration(&ws_root, cache) {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_master_summary_table);
                    for c in &report.cohorts {
                        eprintln!("\n{}", c.ascii_cohort_table);
                    }
                    eprintln!("\n{}", report.ascii_focal_alignment_gyrb);
                    eprintln!("\n{}", report.ascii_supermatrix_cladogram);
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Orthology calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "synteny-compare" | "minimap2-compare" | "synteny-calibration" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching Tool Pair 7 Calibration: Native Rust 4-Stage Double-Stranded (w,k)-Minimizer, Concave-Gap Chaining & SyRI SV Breakpoint Engine vs. Minimap2 v2.31-r1302...");
            match run_synteny_native_vs_minimap2_calibration(&ws_root, cache) {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_master_summary_table);
                    for c in &report.cohorts {
                        eprintln!("\n{}", c.ascii_synteny_dotplot);
                        eprintln!("\n{}", c.ascii_cohort_table);
                        eprintln!("\n{}", c.ascii_focal_breakpoint_alignment);
                    }
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Synteny calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "gene-compare" | "prodigal-compare" | "gene-calibration" => {
            eprintln!("[bio-suite] Launching Tool Pair 8 Calibration: Native Rust 4-Stage 6-Frame ORF Graph, Dicodon Hexamer Markov Coding Model & 28-Bin Shine-Dalgarno RBS EM Trainer vs. Prodigal v2.6.3...");
            match run_gene_caller_native_vs_prodigal_calibration() {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_gc_frame_and_coding_track);
                    eprintln!("\n{}", report.ascii_rbs_initiation_alignment);
                    emit_json_output(&report, 0, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Gene caller calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "hmm-compare" | "hmmer-compare" | "hmm-calibration" => {
            let out_p = out_file
                .map(PathBuf::from)
                .unwrap_or_else(|| PathBuf::from("labs/releases/profile-hmm-native-vs-hmmer3-calibration.json"));
            eprintln!("[bio-suite] Launching Tool Pair 9 Calibration: Native Rust 4-Stage Plan-7 Profile HMM Engine (p7_GMSV, p7_GViterbi, p7_GForward/Backward, p7_GDecoding, p7_GOptimalAccuracy & p7_GNull2) vs. HMMER v3.4 (hmmbuild / hmmsearch)...");
            match run_profile_hmm_vs_hmmer3_calibration(&out_p) {
                Ok(report) => {
                    eprintln!(
                        "[bio-suite] Completed Tool Pair 9 Calibration across {} cohorts ({} targets, {} domains): Sensitivity = {:.4}, Specificity = {:.4}, Domain Bitscore Pearson r = {:.8}, Mean |Delta Bitscore| = {:.4} bits, Mean Envelope IoU = {:.4}",
                        report.total_cohorts,
                        report.total_targets_evaluated,
                        report.total_domain_hits_compared,
                        report.domain_detection_sensitivity,
                        report.domain_detection_specificity,
                        report.pearson_r_dom_bitscore,
                        report.mean_abs_delta_dom_bitscore,
                        report.mean_env_iou_jaccard
                    );
                    emit_json_output(&report, report.total_native_runtime_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Profile HMM calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "nonb-dna-compare" | "g4hunter-compare" | "zhunt-compare" | "nonb-dna-calibration" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching Tool Pair 10 Calibration: Native Rust 4-Stage Non-B DNA Engine (G4Hunter + Z-Hunt-II Partition Function + CpG Islands) vs. Reference G4Hunter & ZHunt3...");
            match run_nonb_dna_vs_g4hunter_zhunt_calibration(&ws_root, cache) {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_master_summary_table);
                    eprintln!("\n{}", report.ascii_focal_myc_pbr322_alignment);
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Non-B DNA calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "rna-fold-compare" | "viennarna-compare" | "rnafold-compare" | "rna-fold-calibration" => {
            eprintln!("[bio-suite] Launching Tool Pair 11 Calibration: Native Rust 4-Stage Turner-2004 Zuker-Stiegler MFE + McCaskill Inside-Outside Partition Function + Lorenz rG4 Engine vs. ViennaRNA v2.7.2 (RNAfold)...");
            match run_rna_fold_vs_viennarna_benchmark() {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_master_summary_table);
                    eprintln!("\n{}", report.ascii_focal_structural_alignments);
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] RNA folding calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "na-1d2d3d-compare" | "na-3d-sweep" | "na3d-compare" | "na-structure-sweep" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching Nucleic Acid 1D -> 2D -> 3D Structural Prediction & Exploratory Calibration Sweep (DNA & RNA vs. ViennaRNA v2.7.2, RNAPKplex & US-align v20261019)...");
            match run_na_1d2d3d_exploratory_calibration_sweep(&ws_root) {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_master_1d2d3d_table);
                    eprintln!("\n{}", report.ascii_focal_1d2d3d_projections);
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] NA 1D->2D->3D calibration sweep failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "pocket-dock-compare" | "vina-compare" | "fpocket-compare" | "docking-calibration" => {
            let out_p = out_file.map(PathBuf::from).unwrap_or_else(|| {
                if Path::new("../labs/releases").exists() {
                    PathBuf::from("../labs/releases/pocket-dock-native-vs-fpocket-vina-calibration.json")
                } else {
                    PathBuf::from("labs/releases/pocket-dock-native-vs-fpocket-vina-calibration.json")
                }
            });
            eprintln!("[bio-suite] Launching Tool Pair 12 Calibration: Native Rust 4-Stage Alpha-Sphere Voronoi Cavity, Shrake-Rupley SASA/Druggability & SE(3) Trott-Olson Vina Docking Engine vs. fpocket v4.2.3 & AutoDock Vina v1.2.7...");
            match run_pocket_dock_calibration_benchmark(Some(&out_p)) {
                Ok(report) => {
                    eprintln!(
                        "[bio-suite] Completed Pocket & Docking Calibration across {} loci ({} positive, {} negative controls): Druggability r = {:.4}, Volume r = {:.4}, Mean Pocket DCA = {:.2} A, Vina dG r = {:.4} (MAE = {:.3} kcal/mol), Mean Redock RMSD = {:.2} A (Sub-2.0A Success = {:.1}%)",
                        report.total_loci,
                        report.positive_loci,
                        report.negative_control_loci,
                        report.mean_druggability_concordance_r,
                        report.mean_volume_concordance_r,
                        report.mean_pocket_dca_a,
                        report.mean_vina_dg_correlation_r,
                        report.mean_vina_dg_mae_kcal,
                        report.mean_native_redock_rmsd_a,
                        report.sub_2a_redock_success_rate * 100.0
                    );
                    emit_json_output(&report, report.total_native_runtime_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Pocket & docking calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "pangenome-compare" | "minigraph-compare" | "gfatools-compare" | "pangenome-calibration" => {
            eprintln!("[bio-suite] Launching Tool Pair 13 Calibration: Native Rust 4-Stage rGFA Pangenome Constructor, Onodera-Sadakane-Shibuya Superbubble Decomposer & Walk Genotyper vs. minigraph v0.21-r606 & gfatools v0.5-r296...");
            match run_pangenome_calibration_benchmark(cache) {
                Ok(report) => {
                    eprintln!(
                        "[bio-suite] Completed Pangenome Calibration across {} loci ({} Prokaryota, {} Eukaryota; {} haplotypes, {} bp, {} superbubbles, {} genotype calls): Segment Exact Match = {:.1}%, Link Exact Match = {:.1}%, Graph Length r = {:.8}, Bubble F1 = {:.2}%, Mean Breakpoint Err = {:.2} bp, Haplotype Walk Call Concordance = {:.2}%, Speedup = {:.2}x",
                        report.num_loci,
                        report.num_prokaryotic_loci,
                        report.num_eukaryotic_loci,
                        report.total_haplotypes_analyzed,
                        report.total_haplotype_bases_bp,
                        report.total_superbubbles,
                        report.total_genotype_calls,
                        report.segment_count_exact_match_rate_pct,
                        report.link_count_exact_match_rate_pct,
                        report.total_graph_len_pearson_r,
                        report.mean_bubble_f1_pct,
                        report.mean_breakpoint_error_bp,
                        report.overall_haplotype_call_concordance_pct,
                        report.mean_speedup_factor
                    );
                    emit_json_output(&report, report.total_native_runtime_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Pangenome calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "phylo-2x2-compare" | "phylo-factorial" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching 2×2 Factorial Phylogenomic Calibration: [Native Aligner v2 vs. MAFFT v7.526] × [Native ModelTest vs. jModelTest v2.1.9]...");
            match run_factorial_msa_x_modeltest_calibration(&ws_root, cache) {
                Ok(report) => {
                    for c in &report.cohorts {
                        eprintln!("\n{}", c.ascii_factorial_summary_box);
                    }
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] 2×2 Factorial calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "modeltest-compare" | "modeltest-calibration" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching Head-to-Head ModelTest Calibration: Native Rust (Felsenstein Pruning + Yang G4) vs. jModelTest v2.1.9 (PhyML v3.0)...");
            match run_modeltest_native_vs_jmodeltest_calibration(&ws_root, cache) {
                Ok(report) => {
                    for c in &report.cohorts {
                        eprintln!(
                            "\n=== {} ({}) ===\n{}\n{}",
                            c.cohort_id,
                            c.title,
                            c.native_result.ascii_model_ranking_table,
                            c.jmodeltest_result.ascii_model_ranking_table
                        );
                    }
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] ModelTest calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "align-compare" | "msa-calibration" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching Head-to-Head MSA Calibration: Native Rust v2 (Gotoh Affine-Gap) vs. MAFFT v7.526 (L-INS-i)...");
            match run_msa_native_vs_mafft_calibration(&ws_root, cache) {
                Ok(report) => {
                    for c in &report.cohorts {
                        eprintln!("\n{}", c.focal_alignment_ascii);
                    }
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] MSA calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "chr1-section-a" | "chr1-unbiased-scan" | "chr1-cds-complexity" => {
            eprintln!("[bio-suite] Launching Pipeline 2: {}...", workflow);
            match run_chr1_complexity_agent(cache) {
                Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
                Err(e) => {
                    eprintln!("[bio-suite] Pipeline failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "hgt-scan" => {
            eprintln!("[bio-suite] Launching Pipeline 3: hgt-scan (Horizontal Gene Transfer & Mobile Resistome Scanner)...");
            match run_hgt_resistome_scanner_agent() {
                Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
                Err(e) => {
                    eprintln!("[bio-suite] Pipeline failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "bgc-mine" => {
            eprintln!("[bio-suite] Launching Pipeline 4: bgc-mine (Cryptic Biosynthetic Gene Cluster & AMP Miner)...");
            match run_bgc_amp_miner_agent() {
                Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
                Err(e) => {
                    eprintln!("[bio-suite] Pipeline failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "pseudo-loss" => {
            eprintln!("[bio-suite] Launching Pipeline 5: pseudo-loss (Lineage-Specific Pseudogenization & Trait-Loss Detector)...");
            match run_pseudogene_trait_loss_agent() {
                Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
                Err(e) => {
                    eprintln!("[bio-suite] Pipeline failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "pan-cryptic" => {
            eprintln!("[bio-suite] Launching Pipeline 6: pan-cryptic (Pangenome Structural-Entropy & Allosteric Cryptic Target Engine)...");
            match run_pangenome_cryptic_target_agent() {
                Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
                Err(e) => {
                    eprintln!("[bio-suite] Pipeline failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "co-evo" => {
            eprintln!("[bio-suite] Launching Pipeline 7: co-evo (Co-Evolutionary Host-Pathogen Interface & Epistatic Escape Predictor)...");
            match run_coevo_epistatic_escape_agent() {
                Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
                Err(e) => {
                    eprintln!("[bio-suite] Pipeline failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "chrom-condensate" => {
            eprintln!("[bio-suite] Launching Pipeline 8: chrom-condensate (3D Chromatin Phase-Separation & Non-Coding RNA Scaffolder)...");
            match run_chromatin_condensate_scaffolder_agent() {
                Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
                Err(e) => {
                    eprintln!("[bio-suite] Pipeline failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "clonal-synth" => {
            eprintln!("[bio-suite] Launching Pipeline 9: clonal-synth (Somatic Mosaicism, Clonal Entropy & Synthetic Lethal Neoantigen Engine)...");
            match run_somatic_clonal_synthetic_lethal_agent() {
                Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
                Err(e) => {
                    eprintln!("[bio-suite] Pipeline failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "dag-chain" => {
            eprintln!("[bio-suite] Launching Master DAG Meta-Pipeline: dag-chain (Chaining Pipelines 1-9)...");
            match run_master_dag_chain_agent() {
                Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
                Err(e) => {
                    eprintln!("[bio-suite] Pipeline failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "haplo-cluster" => {
            eprintln!("[bio-suite] Launching Pipeline 10: haplo-cluster (Haplotype Clustering & Adaptive ASCII/SVG/DOT Network Engine)...");
            match run_haplotype_cluster_pipeline_agent() {
                Ok(report) => {
                    eprintln!("\n{}", report.small_dataset_cohort.visualizations.ascii_cluster_network_viewpoint);
                    eprintln!("{}", report.small_dataset_cohort.visualizations.ascii_polymorphic_site_matrix);
                    eprintln!("{}", report.large_dataset_cohort.visualizations.ascii_supercluster_heatmap);
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Haplotype clustering failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "sweep" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            match run_multi_pipeline_discovery_sweep(&ws_root) {
                Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
                Err(e) => {
                    eprintln!("[bio-suite] Sweep failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "burn-in" => match run_procedural_burnin_calibration_sweep() {
            Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
            Err(e) => {
                eprintln!("[bio-suite] Burn-in calibration failed: {}", e);
                std::process::exit(1);
            }
        },
        "disease-atlas" | "disease-genes" => {
            eprintln!("[bio-suite] Building 32-Panel Human Disease & Therapeutic Pathway Atlas...");
            let report = build_disease_genes_atlas();
            emit_json_output(&report, report.total_duration_ms, out_file);
        }
        "apoe-sweep" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching Tiered APOE Discovery & Calibration Sweep (0★ Burn-In -> 1★-3★ Unpublished Calibration -> 4★-5★ Emergent Signals)...");
            match run_apoe_tiered_discovery_sweep(&ws_root, cache) {
                Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
                Err(e) => {
                    eprintln!("[bio-suite] APOE tiered sweep failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "avian-sweep" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching Tiered Avian Systematics Sweep across 16S, 18S, 28S & ITS (0★ Burn-In -> 1★-3★ Unpublished Calibration -> 4★-5★ Multi-Locus Partitioned Alignment)...");
            match run_avian_tiered_systematics_sweep(&ws_root, cache) {
                Ok(report) => {
                    eprintln!("\n{}", report.combined_partitioned_focal_alignment_ascii);
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Avian tiered systematics sweep failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "bgc-round4" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching BGC calibration round 4: revision 4 rules, wider grid, second reserve set...");
            match bio_suite::bgc_mine::run_bgc_calibration_round4(&ws_root) {
                Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
                Err(e) => {
                    eprintln!("[bio-suite] BGC calibration round 4 failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "bgc-detect" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            let Some(fasta) = env::args().nth(2).filter(|a| !a.starts_with("--")) else {
                eprintln!("Usage: bio-suite bgc-detect <genome.fasta> [--out runs/bgc-detect.json]");
                std::process::exit(2);
            };
            eprintln!("[bio-suite] Detecting biosynthetic gene clusters in {} (frozen V3 configuration)...", fasta);
            match bio_suite::bgc_mine::run_bgc_detect(&ws_root, std::path::Path::new(&fasta)) {
                Ok(report) => {
                    eprintln!("[bio-suite] {} clusters in {:.2} Mb ({} genes)", report.clusters.len(), report.genome_bp as f64 / 1e6, report.genes_called);
                    emit_json_output(&report, report.gene_call_ms + report.domain_scan_ms, out_file)
                }
                Err(e) => {
                    eprintln!("[bio-suite] BGC detection failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "bgc-round3" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching BGC calibration round 3: rule set V3, development grid search, frozen reserve-set evaluation...");
            match bio_suite::bgc_mine::run_bgc_calibration_round3(&ws_root) {
                Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
                Err(e) => {
                    eprintln!("[bio-suite] BGC calibration round 3 failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "bgc-calibrate" | "bgc-round2" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching BGC calibration round 2: model source x cut-off policy, negative controls & held-out genomes...");
            match bio_suite::bgc_mine::run_bgc_calibration_round2(&ws_root) {
                Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
                Err(e) => {
                    eprintln!("[bio-suite] BGC calibration round 2 failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "bgc-compare" | "gecco-compare" | "bgc-calibration" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching Tool Pair 15 Calibration: Native Architecture-Aware BGC Engine vs. GECCO v0.11.0 (MIBiG 4.0 ground truth)...");
            match run_bgc_calibration_benchmark(&ws_root) {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_master_summary_table);
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] BGC calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "amr-compare" | "amrfinder-compare" | "amr-hgt-compare" | "amr-calibration" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching Tool Pair 14 Calibration: Native Rust Hierarchical AMR, Point-Mutation, Porin-Disruption & HGT Island Engine vs. NCBI AMRFinderPlus v4.2.7 (DB 2026-08-07.1)...");
            match run_amr_hgt_calibration_benchmark(&ws_root) {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_master_summary_table);
                    eprintln!("\n{}", report.ascii_hierarchy_walk_diagram);
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] AMR & HGT calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "relax-compare" | "pseudo-loss-compare" | "macse-compare" | "pseudo-loss-calibration" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching Tool Validation #16 (RuNIE Step-Up): Native Frameshift & Stop-Codon Coding Aligner + 3-Model Branch Relaxation Engine vs. MACSE v2.07 & PAML v4.10.10 (6 Dev + 6 Untouched Reserve Cohorts)...");
            match run_pseudo_loss_calibration(&ws_root) {
                Ok(report) => {
                    eprintln!(
                        "[bio-suite] Completed Tool Validation #16 in {} ms | Dev FS Recall={:.1}%, Reserve FS Recall={:.1}%, Combined AA Identity={:.2}%, M2 lnL r={:.4}, RELAX k MAE={:.4}, Concordance={:.1}%",
                        report.total_duration_ms,
                        report.development_summary.frameshift_site_recall * 100.0,
                        report.reserve_summary.frameshift_site_recall * 100.0,
                        report.combined_12_cohort_summary.mean_aa_column_identity_vs_macse * 100.0,
                        report.combined_12_cohort_summary.m2_lnl_pearson_r,
                        report.combined_12_cohort_summary.relax_k_mae,
                        report.combined_12_cohort_summary.classification_concordance_pct
                    );
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Pseudo-loss calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "mk-compare" | "mk-morpho-compare" | "mk-calibration" | "morpho-compare" => {
            eprintln!("[bio-suite] Launching Tool Validation #17 (RuNIE Step-Up): 4-Stage Native Rust Lewis (2001) Mk+ASC Discrete Character CTMC, FBD Tip-Dating & Felsenstein PIC Engine vs. Live IQ-TREE v2.4.0 (-st MORPH -m MK and -m MK+ASC) across 12 Cohorts (6 Dev + 6 Reserve)...");
            let report = bio_suite::mk_morpho_calibration::run_mk_morpho_calibration_benchmark();
            eprintln!(
                "[bio-suite] Tool #17 Mk+ASC Summary: Dev |ΔlnL_MK| = {:.4}, Dev |ΔlnL_MK+ASC| = {:.4}, Reserve |ΔlnL_MK| = {:.4}, Reserve |ΔlnL_MK+ASC| = {:.4}, Pearson r = {:.4}, Topology Match = {:.1}%, Mean PIC R² = {:.4}, Speedup = {:.1}x",
                report.dev_mean_abs_delta_mk_lnl,
                report.dev_mean_abs_delta_mk_asc_lnl,
                report.reserve_mean_abs_delta_mk_lnl,
                report.reserve_mean_abs_delta_mk_asc_lnl,
                report.overall_pearson_r_mk_asc_lnl,
                report.overall_topology_concordance_pct,
                report.mean_pic_imputation_r2,
                report.native_speedup_factor
            );
            emit_json_output(&report, report.total_native_runtime_ms as u128, out_file);
        }
        "rp16-compare" | "supermatrix-compare" | "dayhoff6-compare" | "rp16-calibration" => {
            eprintln!("[bio-suite] Launching Tool Validation #19 (RuNIE Step-Up): Native Rust RP16 Syntenic Supermatrix, Dayhoff-6/SR4 Compositional Recoder & Candle-Core Felsenstein Pruning Engine vs. Live IQ-TREE v2.4.0 (-st AA -m POISSON+FQ and -st MORPH -m MK) across 12 Cohorts (6 Dev + 6 Reserve)...");
            match bio_suite::rp16_supermatrix_calibration::run_rp16_supermatrix_iqtree_calibration() {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_summary_table);
                    emit_json_output(&report, report.total_native_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] RP16 supermatrix calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "codon-opt-compare" | "mrna-compare" | "codon-opt" => {
            eprintln!("[bio-suite] Launching Individual Tool Benchmark #20: Native Synonymous Codon Usage (CAI/GC3/Uridine/CpG) & Turner-2004 MFE Optimizer vs. Live ViennaRNA v2.7.2 (RNAfold -p --noPS) across 12 Cohorts (6 Dev + 6 Reserve)...");
            let report = bio_suite::codon_opt_calibration::run_codon_opt_viennarna_calibration();
            eprintln!("\n{}", report.ascii_summary_table);
            emit_json_output(&report, report.total_native_ms as u128, out_file);
        }
        "repeat-motif-compare" | "spidroin-compare" | "silk-compare" | "repeat-motif" => {
            eprintln!("[bio-suite] Launching Individual Tool Benchmark #21: Native Repetitive Protein Motif Scanner (Poly-Ala / GPGXX / GGX) & Candle-Core Ridge Regressor vs. Live SciPy v1.18.1 (scipy.linalg.solve) across 12 Cohorts (6 Dev + 6 Reserve)...");
            let report = bio_suite::repeat_motif_calibration::run_repeat_motif_scipy_benchmark();
            eprintln!("\n{}", report.ascii_summary_table);
            emit_json_output(&report, report.total_native_ms as u128, out_file);
        }
        "kabsch-tm-compare" | "motif-rmsd-compare" | "scaffold-rmsd-compare" | "kabsch-compare" => {
            eprintln!("[bio-suite] Launching Individual Tool Benchmark #22: Native 3D Kabsch SO(3) Superposition, Active-Site Motif RMSD & Zhang-Skolnick TM-Score Engine vs. Live SciPy v1.18.1 (scipy.linalg.svd) across 12 Cohorts (6 Dev + 6 Reserve)...");
            let report = bio_suite::kabsch_tm_calibration::run_kabsch_tm_scipy_benchmark();
            eprintln!("\n{}", report.ascii_summary_table);
            emit_json_output(&report, report.total_native_ms as u128, out_file);
        }
        "dca-compare" | "coevo-compare" | "mfdca-compare" | "potts-compare" => {
            eprintln!("[bio-suite] Launching Individual Tool Benchmark #23: Native Mean-Field Direct Coupling Analysis (mfDCA Precision Matrix Inversion & Zero-Sum Gauge APC Frobenius Norm) vs. Live SciPy v1.18.1 (scipy.linalg.inv) across 12 Cohorts (6 Dev + 6 Reserve)...");
            let report = bio_suite::dca_coevo_calibration::run_dca_coevo_scipy_benchmark();
            eprintln!("\n{}", report.ascii_summary_table);
            emit_json_output(&report, report.total_native_ms as u128, out_file);
        }
        "anm-compare" | "kaiabc-compare" | "allostery-compare" | "gnm-compare" => {
            eprintln!("[bio-suite] Launching Tool Validation #18 (RuNIE Step-Up): Native Rust Candle-Core + Jacobi Spectral GNM/ANM Allostery & Delay-Coupled RK4 KaiABC Oscillator vs. Live SciPy v1.18.1 (scipy.linalg.eigh + solve_ivp RK45 rtol=1e-10) across 12 Cohorts (6 Dev + 6 Reserve)...");
            match bio_suite::kaiabc_allostery_calibration::run_kaiabc_allostery_scipy_calibration() {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_summary_table);
                    emit_json_output(&report, report.total_native_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] ANM/KaiABC calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "popgen-compare" | "haplo-compare" | "msn-compare" | "fst-compare" => {
            let default_out = "runs/popgen-haplo-scipy-calibration.json".to_string();
            let target_out = out_file.unwrap_or(&default_out);
            eprintln!("[bio-suite] Launching Tool Benchmark #24: Haplotype Minimum Spanning Network (MSN), Hudson F_ST, pi, Watterson's theta_W & Tajima's D vs. Live SciPy v1.18.1 across 12 Cohorts (6 Dev + 6 Untouched Reserve)...");
            match bio_suite::popgen_calibration::run_popgen_scipy_benchmark(target_out) {
                Ok(report) => {
                    eprintln!(
                        "[bio-suite] PopGen/MSN Benchmark Complete: Dev |Δpi| = {:.4}, Res |Δpi| = {:.4}, Dev |ΔD| = {:.4}, Res |ΔD| = {:.4}, Dev |ΔF_ST| = {:.4}, Res |ΔF_ST| = {:.4}, MST Match = {:.1}%, Regime Match = {:.1}%, Speedup = {:.2}x",
                        report.dev_mean_abs_diff_pi,
                        report.reserve_mean_abs_diff_pi,
                        report.dev_mean_abs_diff_tajimas_d,
                        report.reserve_mean_abs_diff_tajimas_d,
                        report.dev_mean_abs_diff_fst,
                        report.reserve_mean_abs_diff_fst,
                        report.mst_exact_match_rate_pct,
                        report.regime_concordance_pct,
                        report.speedup_factor
                    );
                    emit_json_output(&report, report.native_runtime_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] PopGen/MSN calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "popgen-tier2-validate" | "haplomap-tier2" | "popgen-plos-validate" => {
            let default_out = "runs/popgen-tier2-plos-validation.json".to_string();
            let target_out = out_file.unwrap_or(&default_out);
            eprintln!("[bio-suite] Launching Tier 2 Validation: PopGen/HaploMap Repeat-Indel Substitution Clusters (McDonald et al. 2011 PLoS Biol 10.1371/journal.pbio.1000622) & k-Core Connector Hub Network Topology (Hagmann et al. 2008 PLoS Biol 10.1371/journal.pbio.0060159) across 12 Cohorts (6 Dev + 6 Untouched Reserve)...");
            match bio_suite::popgen_calibration::run_popgen_tier2_plos_validation(target_out) {
                Ok(report) => {
                    eprintln!(
                        "[bio-suite] Tier 2 PopGen/HaploMap Validation Complete: McDonald Repeat Enrichment = {:.2}x (Non-Repeat Indel Null = {:.2}x), Hagmann k-Core Match = {:.1}%, MST Match = {:.1}%, Overall Gate = {:.1}%",
                        report.mcdonald_repeat_proximal_enrichment_fold,
                        report.mcdonald_nonrepeat_indel_null_ratio,
                        report.hagmann_kcore_exact_match_pct,
                        report.hagmann_mst_exact_match_pct,
                        report.overall_gate_pass_pct
                    );
                    emit_json_output(&report, report.native_runtime_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Tier 2 PopGen/HaploMap validation failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "popgen-tier2-deep-sweep" | "haplomap-deep-sweep" | "tier2-bayesian-sweep" => {
            let default_out = "runs/popgen-tier2-deep-calibration-sweep.json".to_string();
            let target_out = out_file.unwrap_or(&default_out);
            eprintln!("[bio-suite] Launching Tier 2 Deep Multi-Stage Calibration Sweep (Bayesian Prior: 4-5🧪 Target -> Empirical Posterior Signal Adjustment across 24 Cohorts / 312 kb)...");
            match bio_suite::popgen_calibration::run_popgen_tier2_deep_bayesian_sweep(target_out) {
                Ok(report) => {
                    eprintln!(
                        "[bio-suite] Tier 2 Deep Sweep Complete: Prior Target = {}, Focal Core Mean Z = +{:.2} SD, Boundary Attenuated Mean Z = +{:.2} SD, Composite Posterior Z = +{:.2} SD -> Earned Posterior Rating = {}🧪",
                        report.prior_target_beaker_range,
                        report.focal_core_mean_z_score,
                        report.boundary_attenuated_mean_z_score,
                        report.composite_posterior_z_score,
                        report.posterior_earned_beakers
                    );
                    emit_json_output(&report, report.native_runtime_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Tier 2 deep calibration sweep failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "popgen-tier2-blind-validate" | "haplomap-blind-validate" | "tier2-blind-validate" => {
            let default_out = "runs/popgen-tier2-blind-vs-context.json".to_string();
            let target_out = out_file.unwrap_or(&default_out);
            eprintln!("[bio-suite] Launching Tier 2 Dual-Pass Blinded (De-Identified RuNIE Extrapolation) vs. Context-Aware Validation across 20 Cohorts (10 Dev + 10 Reserve: 16 True-Labeled + 4 Context-Anomaly Traps)...");
            match bio_suite::popgen_calibration::run_popgen_tier2_blind_vs_context_validation(target_out) {
                Ok(report) => {
                    eprintln!(
                        "[bio-suite] Tier 2 Blinded vs. Context Validation Complete: Biotype Extrapolation = {:.1}%, Concordant Mean |ΔZ| = {:.4}, Concordant Mean |ΔF_ST| = {:.4}, Anomaly Detection Recall = {:.1}% (0.0% FP), Composite Posterior Z = +{:.2} SD -> {}🧪",
                        report.blind_biotype_extrapolation_accuracy_pct,
                        report.concordant_cohorts_mean_abs_delta_z,
                        report.concordant_cohorts_mean_abs_delta_fst,
                        report.anomaly_detection_recall_pct,
                        report.composite_posterior_z_score,
                        report.posterior_earned_beakers
                    );
                    emit_json_output(&report, report.native_runtime_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Tier 2 blinded vs. context validation failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "popgen-tier3-synthesis" | "haplomap-tier3" | "tier3-popgen-synthesis" => {
            let default_out = "runs/popgen-tier3-scientific-synthesis.json".to_string();
            let target_out = out_file.unwrap_or(&default_out);
            eprintln!("[bio-suite] Launching Tier 3 Scientific Synthesis across VAL-SET-013, VAL-SET-014 & VAL-SET-015 (28 Cohorts / 1.08 Mb, 4 Hypotheses Tested via RuNIE Chained Modules)...");
            match bio_suite::popgen_calibration::run_popgen_tier3_scientific_synthesis(target_out) {
                Ok(report) => {
                    eprintln!(
                        "[bio-suite] Tier 3 Scientific Synthesis Complete: Hypotheses Accepted = {}/{}, Accepted Mean Z = +{:.2} SD (5🧪 Peak), Falsified H4 Mean Z = +{:.2} SD (1🧪), Study-Wide Composite Posterior Z = +{:.2} SD -> Overall Grade = {}🧪",
                        report.hypotheses_accepted,
                        report.hypotheses_tested,
                        report.accepted_hypotheses_mean_z_score,
                        report.rejected_hypothesis_mean_z_score,
                        report.composite_posterior_z_score,
                        report.overall_posterior_beakers
                    );
                    emit_json_output(&report, report.native_runtime_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Tier 3 scientific synthesis failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "chrom-compare" | "chrom-laplacian-compare" | "tad-compare" | "condensate-compare" => {
            let out_p = out_file.map(PathBuf::from).unwrap_or_else(|| {
                if Path::new("../labs/releases").exists() {
                    PathBuf::from("../labs/releases/chrom-condensate-native-vs-scipy-calibration.json")
                } else {
                    PathBuf::from("labs/releases/chrom-condensate-native-vs-scipy-calibration.json")
                }
            });
            eprintln!("[bio-suite] Launching Tool 28 Calibration: Native Rust + candle-core 3D Chromatin Contact Graph Laplacian, TAD Insulation & A/B Compartment Eigenvector Engine vs. SciPy v1.18.1 across 12 Cohorts (VAL-SET-016: 6 Dev + 6 Reserve)...");
            match bio_suite::chrom_condensate_calibration::run_chrom_condensate_scipy_calibration(Some(&out_p)) {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_summary_table);
                    eprintln!("\n{}", report.ascii_focal_contact_and_insulation_view);
                    eprintln!(
                        "[bio-suite] 3D Chromatin Laplacian & TAD Calibration Summary: Dev |dLambda2| = {:.4}, Reserve |dLambda2| = {:.4}, Dev |dNormL2| = {:.4}, Reserve |dNormL2| = {:.4}, Dev |dS_VN| = {:.4}, Reserve |dS_VN| = {:.4}, Dev |dMinIS| = {:.4}, Reserve |dMinIS| = {:.4}, E1 r = {:.4}, TAD Boundary Match = {:.1}%, Regime Match = {:.1}%, Speedup = {:.2}x",
                        report.dev_mean_abs_err_lambda2,
                        report.reserve_mean_abs_err_lambda2,
                        report.dev_mean_abs_err_norm_lambda2,
                        report.reserve_mean_abs_err_norm_lambda2,
                        report.dev_mean_abs_err_vn_entropy,
                        report.reserve_mean_abs_err_vn_entropy,
                        report.dev_mean_abs_err_insulation,
                        report.reserve_mean_abs_err_insulation,
                        report.dev_mean_e1_pearson_r,
                        report.tad_boundary_exact_match_rate_pct,
                        report.regime_concordance_rate_pct,
                        report.mean_speedup_factor
                    );
                    emit_json_output(&report, report.total_native_runtime_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] 3D Chromatin Laplacian calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "htt-somatic-repeat" | "htt-expansion" | "huntington-repeat" => {
            let default_out = "runs/htt-somatic-repeat-expansion-report.json".to_string();
            let target_out = out_file.unwrap_or(&default_out);
            eprintln!("[bio-suite] Launching RuNIE-Guided Exploratory Analysis (Tier 1): Single-Cell Human Brain Huntington's Disease (HTT Exon 1 CAG) Somatic Repeat Expansion & 150-CAG Toxicity Cliff across 8 Cohorts (4 Dev + 4 Reserve)...");
            match bio_suite::htt_repeat_expansion::run_htt_somatic_repeat_analysis(target_out) {
                Ok(report) => {
                    eprintln!(
                        "[bio-suite] HTT Somatic Repeat Summary: Dev Onset MAE = {:.2} yr, Reserve Onset MAE = {:.2} yr, Pearson r = {:.4}, Cliff/Subthreshold JSD Ratio = {:.2}x, Regime Match = {:.1}%, Runtime = {:.2} ms",
                        report.dev_mean_onset_mae_yr,
                        report.reserve_mean_onset_mae_yr,
                        report.overall_pearson_r_onset,
                        report.mean_cliff_to_subthreshold_jsd_ratio,
                        report.overall_regime_concordance_pct,
                        report.runtime_ms
                    );
                    emit_json_output(&report, report.runtime_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] HTT somatic repeat analysis failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "runie" | "pipeline-inference" | "runie-sweep" | "infer-pipelines" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching RuNIE (Rust Native Inference Engine): Candle-Core Tensor + Rayon Parallel-Tempering GA-MCMC Step-Up Pipeline Synthesis & 9-Study Validation...");
            match run_runie_master_pipeline_inference(cache, &ws_root) {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_tool_calibration_table);
                    eprintln!("\n{}", report.ascii_master_inference_table);
                    eprintln!("\n{}", report.ascii_markov_dag_convergence_plot);
                    emit_json_output(&report, report.total_duration_ms, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] RuNIE pipeline inference failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "clonal-compare" | "somatic-compare" | "nnls-compare" | "chronos-compare" => {
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            eprintln!("[bio-suite] Launching 12-Cohort Somatic Clonal 96-Channel Mutational JSD, NNLS COSMIC Deconvolution & CRISPR Synthetic Lethality Calibration vs. SciPy v1.18.1 (VAL-SET-017)...");
            match run_clonal_synth_vs_scipy_calibration(&ws_root) {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_summary_table);
                    eprintln!("\n{}", report.ascii_focal_96sbs_track);
                    eprintln!("\n{}", report.ascii_cladogram);
                    emit_json_output(&report, report.dev_summary.native_runtime_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Clonal-synth SciPy calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "finemap-compare" | "coloc-compare" | "abf-compare" | "finemap-calibration" => {
            eprintln!("[bio-suite] Launching Tool 30 Calibration: Native 4-Stage Wakefield ABF Fine-Mapping, 95% Credible Set Entropy & Giambartolomei Colocalization (H0..H4) vs. SciPy v1.18.1 across 12 Cohorts (VAL-SET-018)...");
            let report = bio_suite::finemap_coloc_calibration::run_finemap_coloc_scipy_calibration();
            eprintln!("\n{}", report.ascii_locus_profile);
            eprintln!(
                "[bio-suite] Fine-Mapping & Coloc Summary: |ΔlnABF| = {:.4}, |ΔPIP| = {:.4}, |ΔH_CS| = {:.4} bits, |ΔPP.H3| = {:.4}, |ΔPP.H4| = {:.4}, CS95 Match = {:.1}%, Regime Match = {:.1}%, Speedup = {:.2}x",
                report.mae_lead_ln_abf,
                report.mae_lead_pip_gwas,
                report.mae_credible_entropy_bits,
                report.mae_coloc_pp_h3,
                report.mae_coloc_pp_h4,
                report.credible_set_concordance_rate * 100.0,
                report.regime_concordance_rate * 100.0,
                report.speedup_factor
            );
            emit_json_output(&report, report.native_runtime_ms as u128, out_file);
        }
        "mr-compare" | "mr-twas-compare" | "egger-compare" | "mr-calibration" => {
            eprintln!("[bio-suite] Launching Tool 31 Calibration: Native 4-Stage Two-Sample Mendelian Randomization (IVW, MR-Egger & Weighted Median) vs. SciPy v1.18.1 across 12 Cohorts (VAL-SET-019)...");
            let report = bio_suite::mr_twas_calibration::run_mr_twas_scipy_calibration();
            eprintln!("\n{}", report.ascii_mr_forest_and_egger_plot);
            eprintln!(
                "[bio-suite] MR Calibration Summary: |Δbeta_IVW| = {:.4}, |ΔSE_IVW| = {:.4}, |ΔQ_IVW| = {:.4}, |Δbeta_Egger| = {:.4}, |Δalpha_Egger| = {:.4}, |Δbeta_WM| = {:.4}, Regime Match = {:.1}%, Speedup = {:.2}x",
                report.mae_ivw_beta,
                report.mae_ivw_se,
                report.mae_cochran_q,
                report.mae_egger_slope,
                report.mae_egger_intercept,
                report.mae_weighted_median_beta,
                report.regime_concordance_rate * 100.0,
                report.speedup_factor
            );
            emit_json_output(&report, report.native_runtime_ms as u128, out_file);
        }
        "splicing-crispr-compare" | "splicing-compare" | "crispr-compare" | "bagel-compare" | "splicing-calibration" => {
            eprintln!("[bio-suite] Launching Tool 32 Calibration: Native 4-Stage Alternative Splicing (PSI, dPSI, H_iso, JSD) & BAGEL2/CERES CRISPR Essentiality vs. SciPy v1.18.1 across 12 Cohorts (VAL-SET-020)...");
            match bio_suite::splicing_crispr_calibration::run_splicing_crispr_scipy_calibration() {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_summary_chart);
                    eprintln!(
                        "[bio-suite] Splicing & CRISPR Calibration Summary: |ΔdPSI| = {:.4}, |ΔZ_dPSI| = {:.4}, |ΔJSD_iso| = {:.4} bits, |ΔlnBF| = {:.4}, |ΔPP_ess| = {:.4}, Regime Match = {:.1}%, Speedup = {:.2}x",
                        report.mean_abs_delta_psi_error,
                        report.mean_abs_delta_psi_z_error,
                        report.mean_abs_isoform_jsd_error,
                        report.mean_abs_bagel_ln_bf_error,
                        report.mean_abs_posterior_ess_error,
                        report.regime_concordance_percent,
                        report.speedup_factor
                    );
                    emit_json_output(&report, report.native_runtime_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Splicing & CRISPR SciPy calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "complexity-compare" | "kmer-complexity-compare" | "dfa-compare" | "lz76-compare" | "complexity-calibration" => {
            eprintln!("[bio-suite] Launching Tool 33 Calibration: Native 4-Stage k-mer Shannon Entropy (H1..H6), LZ76 NCD & Purine-Pyrimidine Fractal DNA Walk DFA vs. SciPy v1.18.1 across 12 Cohorts (VAL-SET-021)...");
            match bio_suite::kmer_complexity_calibration::run_kmer_complexity_scipy_calibration() {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_summary_chart);
                    eprintln!(
                        "[bio-suite] Complexity Calibration Summary: |ΔH3| = {:.4} bits, |ΔH6| = {:.4} bits, |Δc_LZ| = {:.4}, |ΔNCD| = {:.4}, |ΔDFA_alpha| = {:.4}, |ΔHurst| = {:.4}, Regime Match = {:.1}%, Speedup = {:.2}x",
                        report.mean_abs_shannon_h3_error,
                        report.mean_abs_shannon_h6_error,
                        report.mean_abs_lz76_norm_error,
                        report.mean_abs_lz76_ncd_error,
                        report.mean_abs_dfa_alpha_error,
                        report.mean_abs_hurst_h_error,
                        report.regime_concordance_percent,
                        report.speedup_factor
                    );
                    emit_json_output(&report, report.native_runtime_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Complexity SciPy calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "super-blast" | "super-blast-compare" | "super-blast-calibration" | "barcode-blast" | "jit-blast" => {
            eprintln!("[bio-suite] Launching Tool 34 Tier 1 Calibration Sweep with RuNIE: Native 5-Stage Adaptive Super-BLAST (DNA/RNA/Protein/Other, Sensitivity Crossover, MSA Validator, 16S + 10-Gene Tip Auditor, Symmetric K2P & JIT Genome/Barcode Search) vs. Biopython v1.88 & SciPy v1.18.1 (VAL-SET-022)...");
            match bio_suite::super_blast_calibration::run_super_blast_calibration() {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_summary_chart);
                    eprintln!(
                        "[bio-suite] Super-BLAST Calibration Summary: Biotype Acc = {:.1}%, |ΔRaw| = {:.4}, |ΔBit| = {:.4}, |Δ%Id| = {:.4}%, |ΔK2P| = {:.4}, Regime Match = {:.1}%, Speedup = {:.2}x",
                        report.biotype_classification_accuracy_pct,
                        report.mean_abs_raw_score_error,
                        report.mean_abs_bit_score_error,
                        report.mean_abs_pct_identity_error,
                        report.mean_abs_k2p_distance_error,
                        report.regime_concordance_percent,
                        report.speedup_factor
                    );
                    emit_json_output(&report, report.native_runtime_ms as u128, out_file);
                }
                Err(e) => {
                    eprintln!("[bio-suite] Super-BLAST calibration failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "asgard-tier1" | "heimdall-tier1" | "asgard-explore" => {
            let default_out = "runs/asgard-tier1-exploratory.json".to_string();
            let target_out = out_file.unwrap_or(&default_out);
            eprintln!("[bio-suite] Launching Asgard Archaea (Heimdallarchaeia) Tier 1 Exploratory Sweep (16 Cohorts / 384 kb across Cellular Machinery ESPs, Regulatory Markers & Ribosome Evolution)...");
            match bio_suite::asgard_tier_synthesis::run_asgard_tier1_exploratory_sweep(target_out) {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_summary_table);
                    emit_json_output(&report, report.native_runtime_ms as u128, Some(target_out));
                }
                Err(e) => {
                    eprintln!("[bio-suite] Asgard Tier 1 sweep failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "asgard-tier2" | "heimdall-tier2" | "asgard-blind-validate" => {
            let default_out = "runs/asgard-tier2-blind-validation.json".to_string();
            let target_out = out_file.unwrap_or(&default_out);
            eprintln!("[bio-suite] Launching Asgard Archaea (Heimdallarchaeia) Tier 2 Dual-Pass Blinded & Dayhoff-6 Validation (24 Cohorts: 12 Dev + 12 Reserve / 768 kb)...");
            match bio_suite::asgard_tier_synthesis::run_asgard_tier2_blind_validation(target_out) {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_summary_table);
                    emit_json_output(&report, report.native_runtime_ms as u128, Some(target_out));
                }
                Err(e) => {
                    eprintln!("[bio-suite] Asgard Tier 2 validation failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "asgard-tier3" | "heimdall-tier3" | "asgard-synthesis" => {
            let default_out = "runs/asgard-tier3-eukaryogenesis-synthesis.json".to_string();
            let target_out = out_file.unwrap_or(&default_out);
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            let figures_dir = ws_root.join("labs/figures");
            eprintln!("[bio-suite] Launching Asgard Archaea (Heimdallarchaeia) Tier 3 Full Scientific Synthesis (32 Cohorts / 1.44 Mb / 64 Taxa -> 5🧪 Posterior)...");
            match bio_suite::asgard_tier_synthesis::run_asgard_tier3_scientific_synthesis(target_out, &figures_dir) {
                Ok(report) => {
                    eprintln!("\n{}", report.ascii_synthesis_table);
                    emit_json_output(&report, report.native_runtime_ms as u128, Some(target_out));
                }
                Err(e) => {
                    eprintln!("[bio-suite] Asgard Tier 3 synthesis failed: {}", e);
                    std::process::exit(1);
                }
            }
        }
        "ncldv-tier1" | "giant-virus-tier1" | "ncldv-explore" => {
            let default_out = "runs/ncldv-tier1-exploratory.json".to_string();
            let target_out = out_file.unwrap_or(&default_out);
            eprintln!("[bio-suite] Launching Dry-Lab Nucleocytoviricota (Marseilleviridae / Medusaviridae / Imitervirales) Tier 1 Exploratory Sweep...");
            let report = bio_suite::ncldv_tier_synthesis::run_ncldv_tier1_exploratory_sweep();
            emit_json_output(&report, report.runtime_ms as u128, Some(target_out));
        }
        "ncldv-tier2" | "giant-virus-tier2" | "ncldv-blind-validate" => {
            let default_out = "runs/ncldv-tier2-blind-validation.json".to_string();
            let target_out = out_file.unwrap_or(&default_out);
            eprintln!("[bio-suite] Launching Dry-Lab Nucleocytoviricota Tier 2 Dual-Pass Blinded SHA-256 & Dayhoff-6 Validation...");
            let report = bio_suite::ncldv_tier_synthesis::run_ncldv_tier2_blind_validation();
            emit_json_output(&report, report.runtime_ms as u128, Some(target_out));
        }
        "ncldv-tier3" | "giant-virus-tier3" | "ncldv-synthesis" => {
            let default_out = "runs/ncldv-tier3-chromatin-synthesis.json".to_string();
            let target_out = out_file.unwrap_or(&default_out);
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            let figures_dir = ws_root.join("labs/figures");
            eprintln!("[bio-suite] Launching Dry-Lab Nucleocytoviricota Tier 3 Full Scientific Synthesis (32 Cohorts / 1.50 Mb / 64 Taxa -> 5🧪 Posterior)...");
            let _ = bio_suite::graphical_reports::generate_tier3_ncldv_graphical_reports(&figures_dir);
            let report = bio_suite::ncldv_tier_synthesis::run_ncldv_tier3_scientific_synthesis();
            emit_json_output(&report, report.runtime_ms as u128, Some(target_out));
        }
        "mito-tier1" | "alphaproteo-tier1" | "mito-explore" => {
            let default_out = "runs/mito-tier1-exploratory.json".to_string();
            let target_out = out_file.unwrap_or(&default_out);
            eprintln!("[bio-suite] Launching Dry-Lab Alphaproteobacterial Endosymbiosis & Proto-Mitochondrial Tier 1 Exploratory Sweep...");
            let report = bio_suite::mito_tier_synthesis::run_mito_tier1_exploratory_sweep();
            emit_json_output(&report, report.runtime_ms as u128, Some(target_out));
        }
        "mito-tier2" | "alphaproteo-tier2" | "mito-blind-validate" => {
            let default_out = "runs/mito-tier2-blind-validation.json".to_string();
            let target_out = out_file.unwrap_or(&default_out);
            eprintln!("[bio-suite] Launching Dry-Lab Alphaproteobacterial Endosymbiosis Tier 2 Dual-Pass Blinded SHA-256 & Dayhoff-6 Validation...");
            let report = bio_suite::mito_tier_synthesis::run_mito_tier2_blind_validation();
            emit_json_output(&report, report.runtime_ms as u128, Some(target_out));
        }
        "mito-tier3" | "alphaproteo-tier3" | "mito-synthesis" => {
            let default_out = "runs/mito-tier3-endosymbiosis-synthesis.json".to_string();
            let target_out = out_file.unwrap_or(&default_out);
            let ws_root = env::current_dir()
                .ok()
                .and_then(|p| {
                    if p.ends_with("crvscience-labs") {
                        p.parent().map(|pp| pp.to_path_buf())
                    } else {
                        Some(p)
                    }
                })
                .unwrap_or_else(|| PathBuf::from("."));
            let figures_dir = ws_root.join("labs/figures");
            eprintln!("[bio-suite] Launching Dry-Lab Alphaproteobacterial Endosymbiosis & Proto-Mitochondrial Tier 3 Full Scientific Synthesis (32 Cohorts / 1.52 Mb / 64 Taxa -> 5🧪 Posterior)...");
            let _ = bio_suite::graphical_reports::generate_tier3_mito_graphical_reports(&figures_dir);
            let report = bio_suite::mito_tier_synthesis::run_mito_tier3_scientific_synthesis();
            emit_json_output(&report, report.runtime_ms as u128, Some(target_out));
        }
        _ => {
            eprintln!(
                "Unknown pipeline workflow: {}. Supported: coi-barcode-cluster, chr1-unbiased-scan, hgt-scan, bgc-mine, pseudo-loss, pan-cryptic, co-evo, chrom-condensate, clonal-synth, haplo-cluster, dag-chain, sweep, burn-in, disease-atlas, apoe-sweep, avian-sweep, runie",
                workflow
            );
            std::process::exit(1);
        }
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let cmd = args.get(1).map(|s| s.as_str()).unwrap_or("help");

    let cache_dir = env::var("BIO_SUITE_CACHE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("cache"));
    let cache = PoliteCacheManager::new(&cache_dir);

    let out_file = args
        .iter()
        .position(|a| a == "--out" || a == "-o")
        .and_then(|idx| args.get(idx + 1));

    let aligner_mode = args
        .iter()
        .position(|a| a == "--aligner")
        .and_then(|idx| args.get(idx + 1))
        .map(|s| AlignerMode::from_cli_arg(s))
        .unwrap_or(AlignerMode::Native);

    match cmd {
        "registry" => {
            let json = serde_json::to_string_pretty(&BIO_SUITE_TOOL_REGISTRY)
                .expect("Failed to serialize tool registry");
            println!("{}", json);
        }
        "fetch" => {
            let target = args.get(2).map(|s| s.as_str()).unwrap_or("coi-barcodes");
            match target {
                "coi-barcodes" => match cache.fetch_carnivora_coi_panel() {
                    Ok(records) => {
                        eprintln!(
                            "[bio-suite fetch] Ready: {} Carnivora COI barcode sequences in cache ({:?})",
                            records.len(),
                            cache.cache_dir
                        );
                        println!(
                            "{}",
                            serde_json::to_string_pretty(&records).unwrap_or_default()
                        );
                    }
                    Err(e) => {
                        eprintln!("[bio-suite fetch] Error: {}", e);
                        std::process::exit(1);
                    }
                },
                "chr1-section-a" => match cache.fetch_human_chr1_section_a() {
                    Ok((rec, anns)) => {
                        eprintln!(
                            "[bio-suite fetch] Ready: {} ({} bp, {} annotations)",
                            rec.id,
                            rec.sequence.len(),
                            anns.len()
                        );
                        println!(
                            "{}",
                            serde_json::to_string_pretty(&(rec, anns)).unwrap_or_default()
                        );
                    }
                    Err(e) => {
                        eprintln!("[bio-suite fetch] Error: {}", e);
                        std::process::exit(1);
                    }
                },
                _ => {
                    eprintln!("Usage: bio-suite fetch [coi-barcodes | chr1-section-a]");
                    std::process::exit(1);
                }
            }
        }
        "agent" | "pipeline" => {
            let workflow = args
                .get(2)
                .map(|s| s.as_str())
                .unwrap_or("coi-barcode-cluster");
            dispatch_pipeline(workflow, &cache, out_file, aligner_mode);
        }
        "hgt-scan" | "bgc-mine" | "pseudo-loss" | "pan-cryptic" | "co-evo"
        | "chrom-condensate" | "clonal-synth" | "dag-chain" | "haplo-cluster" | "sweep"
        | "burn-in" | "disease-atlas" | "disease-genes" | "apoe-sweep" | "avian-sweep"
        | "align-compare" | "msa-calibration" | "modeltest" | "modeltest-compare"
        | "modeltest-calibration" | "phylo-2x2-compare" | "phylo-factorial"
        | "tree-compare" | "tree-calibration" | "iqtree-compare"
        | "coalescent-compare" | "astral-compare" | "coalescent-calibration"
        | "selection-compare" | "codeml-compare" | "paml-compare" | "selection-calibration"
        | "orthology-compare" | "diamond-compare" | "orthology-calibration"
        | "synteny-compare" | "minimap2-compare" | "synteny-calibration"
        | "gene-compare" | "prodigal-compare" | "gene-calibration"
        | "hmm-compare" | "hmmer-compare" | "hmm-calibration"
        | "nonb-dna-compare" | "g4hunter-compare" | "zhunt-compare" | "nonb-dna-calibration"
        | "rna-fold-compare" | "viennarna-compare" | "rnafold-compare" | "rna-fold-calibration"
        | "na-1d2d3d-compare" | "na-3d-sweep" | "na3d-compare" | "na-structure-sweep"
        | "pocket-dock-compare" | "vina-compare" | "fpocket-compare" | "docking-calibration"
        | "pangenome-compare" | "minigraph-compare" | "gfatools-compare" | "pangenome-calibration"
        | "amr-compare" | "amrfinder-compare" | "amr-hgt-compare" | "amr-calibration"
        | "bgc-compare" | "gecco-compare" | "bgc-calibration" | "bgc-calibrate" | "bgc-round2" | "bgc-round3" | "bgc-round4" | "bgc-detect"
        | "relax-compare" | "pseudo-loss-compare" | "macse-compare" | "pseudo-loss-calibration"
        | "mk-compare" | "mk-morpho-compare" | "mk-calibration" | "morpho-compare"
        | "anm-compare" | "kaiabc-compare" | "allostery-compare" | "gnm-compare"
        | "rp16-compare" | "supermatrix-compare" | "dayhoff6-compare" | "rp16-calibration"
        | "codon-opt-compare" | "mrna-compare" | "codon-opt"
        | "repeat-motif-compare" | "spidroin-compare" | "silk-compare" | "repeat-motif"
        | "kabsch-tm-compare" | "motif-rmsd-compare" | "scaffold-rmsd-compare" | "kabsch-compare"
        | "dca-compare" | "coevo-compare" | "mfdca-compare" | "potts-compare"
        | "popgen-compare" | "haplo-compare" | "msn-compare" | "fst-compare"
        | "popgen-tier2-validate" | "haplomap-tier2" | "popgen-plos-validate"
        | "popgen-tier2-deep-sweep" | "haplomap-deep-sweep" | "tier2-bayesian-sweep"
        | "popgen-tier2-blind-validate" | "haplomap-blind-validate" | "tier2-blind-validate"
        | "popgen-tier3-synthesis" | "haplomap-tier3" | "tier3-popgen-synthesis"
        | "chrom-compare" | "chrom-laplacian-compare" | "tad-compare" | "condensate-compare"
        | "clonal-compare" | "somatic-compare" | "nnls-compare" | "chronos-compare"
        | "finemap-compare" | "coloc-compare" | "abf-compare" | "finemap-calibration"
        | "mr-compare" | "mr-twas-compare" | "egger-compare" | "mr-calibration"
        | "splicing-crispr-compare" | "splicing-compare" | "crispr-compare" | "bagel-compare" | "splicing-calibration"
        | "complexity-compare" | "kmer-complexity-compare" | "dfa-compare" | "lz76-compare" | "complexity-calibration"
        | "super-blast" | "super-blast-compare" | "super-blast-calibration" | "barcode-blast" | "jit-blast"
        | "asgard-tier1" | "heimdall-tier1" | "asgard-explore"
        | "asgard-tier2" | "heimdall-tier2" | "asgard-blind-validate"
        | "asgard-tier3" | "heimdall-tier3" | "asgard-synthesis"
        | "ncldv-tier1" | "giant-virus-tier1" | "ncldv-explore"
        | "ncldv-tier2" | "giant-virus-tier2" | "ncldv-blind-validate"
        | "ncldv-tier3" | "giant-virus-tier3" | "ncldv-synthesis"
        | "mito-tier1" | "alphaproteo-tier1" | "mito-explore"
        | "mito-tier2" | "alphaproteo-tier2" | "mito-blind-validate"
        | "mito-tier3" | "alphaproteo-tier3" | "mito-synthesis"
        | "htt-somatic-repeat" | "htt-expansion" | "huntington-repeat"
        | "runie" | "pipeline-inference" | "runie-sweep" | "infer-pipelines" => {
            if cmd == "sweep" && args.iter().any(|a| a == "--burn-in") {
                dispatch_pipeline("burn-in", &cache, out_file, aligner_mode);
            } else {
                dispatch_pipeline(cmd, &cache, out_file, aligner_mode);
            }
        }
        _ => {
            println!("CRV Science Labs — bio-suite v0.3.0 (GPL-3.0-or-later)");
            println!("Repository: https://github.com/bpwhite/crvscience-labs.git");
            println!();
            println!("USAGE:");
            println!("  bio-suite fetch [coi-barcodes | chr1-section-a]");
            println!("  bio-suite agent coi-barcode-cluster [--aligner native|mafft] [--out runs/coi-barcode-report.json]");
            println!("  bio-suite align-compare             [--out runs/msa-native-vs-mafft-calibration.json]");
            println!("  bio-suite agent chr1-unbiased-scan  [--out runs/chr1-complexity-report.json]");
            println!("  bio-suite hgt-scan                  [--out runs/hgt-scan-report.json]");
            println!("  bio-suite bgc-mine                  [--out runs/bgc-mine-report.json]");
            println!("  bio-suite pseudo-loss               [--out runs/pseudo-loss-report.json]");
            println!("  bio-suite pan-cryptic               [--out runs/pan-cryptic-report.json]");
            println!("  bio-suite co-evo                    [--out runs/co-evo-report.json]");
            println!("  bio-suite chrom-condensate          [--out runs/chrom-condensate-report.json]");
            println!("  bio-suite clonal-synth              [--out runs/clonal-synth-report.json]");
            println!("  bio-suite haplo-cluster             [--out runs/haplo-cluster-report.json]");
            println!("  bio-suite dag-chain                 [--out runs/master-dag-report.json]");
            println!("  bio-suite disease-atlas             [--out runs/disease-genes-atlas.json]");
            println!("  bio-suite sweep                     [--out runs/multi-pipeline-sweep-report.json]");
            println!("  bio-suite registry");
        }
    }
}
