//! `bio-suite` — CRV Science Labs High-Throughput Bioinformatics, Molecular Systematics
//! & Information-Theoretic Genome Complexity Engine
//!
//! Released under GNU General Public License v3.0 (GPL-3.0-or-later).
//! Repository: https://github.com/bpwhite/crvscience-labs.git

use bio_suite::agent::{run_chr1_complexity_agent, run_coi_barcode_cluster_agent};
use bio_suite::fetch::PoliteCacheManager;
use bio_suite::genome_target::BIO_SUITE_TOOL_REGISTRY;
use bio_suite::pipelines::{
    run_bgc_amp_miner_agent, run_hgt_resistome_scanner_agent, run_pangenome_cryptic_target_agent,
    run_pseudogene_trait_loss_agent,
};
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

fn dispatch_pipeline(workflow: &str, cache: &PoliteCacheManager, out_file: Option<&String>) {
    match workflow {
        "coi-barcode-cluster" => {
            eprintln!("[bio-suite] Launching Pipeline 1: coi-barcode-cluster (500 Rayon bootstrap replicates)...");
            match run_coi_barcode_cluster_agent(cache, 500) {
                Ok(report) => emit_json_output(&report, report.total_duration_ms, out_file),
                Err(e) => {
                    eprintln!("[bio-suite] Pipeline failed: {}", e);
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
        _ => {
            eprintln!(
                "Unknown pipeline workflow: {}. Supported: coi-barcode-cluster, chr1-unbiased-scan, hgt-scan, bgc-mine, pseudo-loss, pan-cryptic",
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
            dispatch_pipeline(workflow, &cache, out_file);
        }
        "hgt-scan" | "bgc-mine" | "pseudo-loss" | "pan-cryptic" => {
            dispatch_pipeline(cmd, &cache, out_file);
        }
        _ => {
            println!("CRV Science Labs — bio-suite v0.2.0 (GPL-3.0-or-later)");
            println!("Repository: https://github.com/bpwhite/crvscience-labs.git");
            println!();
            println!("USAGE:");
            println!("  bio-suite fetch [coi-barcodes | chr1-section-a]");
            println!("  bio-suite agent coi-barcode-cluster [--out runs/coi-barcode-report.json]");
            println!("  bio-suite agent chr1-unbiased-scan  [--out runs/chr1-complexity-report.json]");
            println!("  bio-suite hgt-scan                  [--out runs/hgt-scan-report.json]");
            println!("  bio-suite bgc-mine                  [--out runs/bgc-mine-report.json]");
            println!("  bio-suite pseudo-loss               [--out runs/pseudo-loss-report.json]");
            println!("  bio-suite pan-cryptic               [--out runs/pan-cryptic-report.json]");
            println!("  bio-suite registry");
        }
    }
}
