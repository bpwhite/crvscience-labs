//! `bio-suite-supervisor` — Standalone CLI & Interactive Multi-Tab GUI Binary
//!
//! Oversees autonomous multi-threaded data fetch, Super-BLAST annotation validation & self-healing
//! loop-back hygiene, full-length sequence alignment, 32-chain RuNIE factorial analysis, locked-format
//! figure generation (`.svg` + `.png` + `.pdf` modeled on Preprints 18 & 21), and arXiv LaTeX (`.tex`)
//! preprint compilation.

use bio_suite::fetch::PoliteCacheManager;
use bio_suite::supervisor::run_supervisor_cli;
use std::env;
use std::path::PathBuf;

fn main() {
    let args: Vec<String> = env::args().collect();
    let cache_dir = env::current_dir()
        .ok()
        .map(|p| {
            if p.join("cache").exists() {
                p.join("cache")
            } else {
                PathBuf::from("cache")
            }
        })
        .unwrap_or_else(|| PathBuf::from("cache"));
    let cache = PoliteCacheManager::new(cache_dir);

    if let Err(e) = run_supervisor_cli(&args, &cache) {
        eprintln!("[bio-suite-supervisor] Error: {}", e);
        std::process::exit(1);
    }
}
