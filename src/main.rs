//! `bio-suite` — CRV Science Labs High-Throughput Bioinformatics, Molecular Systematics
//! & Information-Theoretic Genome Complexity Engine
//!
//! Released under GNU General Public License v3.0 (GPL-3.0-or-later).
//! Repository: https://github.com/bpwhite/crvscience-labs.git

use bio_suite::agent::{
    run_chr1_complexity_agent, run_coi_barcode_cluster_agent_with_mode,
};
use bio_suite::avian_systematics::run_avian_tiered_systematics_sweep;
use bio_suite::dag_pipelines::{
    run_chromatin_condensate_scaffolder_agent, run_coevo_epistatic_escape_agent,
    run_master_dag_chain_agent, run_somatic_clonal_synthetic_lethal_agent,
};
use bio_suite::disease_panels::{build_disease_genes_atlas, run_apoe_tiered_discovery_sweep};
use bio_suite::fetch::PoliteCacheManager;
use bio_suite::genome_target::BIO_SUITE_TOOL_REGISTRY;
use bio_suite::haplotype::run_haplotype_cluster_pipeline_agent;
use bio_suite::modeltest::{
    run_factorial_msa_x_modeltest_calibration, run_modeltest_native_vs_jmodeltest_calibration,
    run_modeltest_with_mode, ModelTestMode,
};
use bio_suite::msa_calibration::run_msa_native_vs_mafft_calibration;
use bio_suite::phylo::{align_sequences_with_mode, AlignerMode};
use bio_suite::pipelines::{
    run_bgc_amp_miner_agent, run_hgt_resistome_scanner_agent, run_pangenome_cryptic_target_agent,
    run_pseudogene_trait_loss_agent,
};
use bio_suite::sweep::{
    run_multi_pipeline_discovery_sweep, run_procedural_burnin_calibration_sweep,
};
use bio_suite::tree_builder::run_tree_builder_native_vs_iqtree_calibration;
use serde::Serialize;
use std::env;
use std::fs;
use std::path::PathBuf;

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
        _ => {
            eprintln!(
                "Unknown pipeline workflow: {}. Supported: coi-barcode-cluster, chr1-unbiased-scan, hgt-scan, bgc-mine, pseudo-loss, pan-cryptic, co-evo, chrom-condensate, clonal-synth, haplo-cluster, dag-chain, sweep, burn-in, disease-atlas, apoe-sweep, avian-sweep",
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
        | "tree-compare" | "tree-calibration" | "iqtree-compare" => {
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
