//! `runie-align-view` — Standalone Rust-Native Sequence Alignment Viewer & Paper Snapshot Binary
//!
//! Launches the interactive AliView-modeled desktop GUI (`runie-align-view [optional_fasta]`)
//! or executes headless pipeline/paper snapshots (`runie-align-view --input <fasta> --region 120:210 --snapshot-svg fig2.svg --snapshot-png fig2.png`).

use bio_suite::align_viewer::run_align_viewer_cli;
use bio_suite::fetch::PoliteCacheManager;
use std::env;
use std::path::PathBuf;

fn main() {
    let args: Vec<String> = env::args().collect();
    let cache_dir = env::var("BIO_SUITE_CACHE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|_| {
            if PathBuf::from("cache").exists() {
                PathBuf::from("cache")
            } else {
                PathBuf::from("crvscience-labs/cache")
            }
        });
    let cache = PoliteCacheManager::new(&cache_dir);

    if let Err(e) = run_align_viewer_cli(&args, &cache) {
        eprintln!("[runie-align-view] Error: {}", e);
        std::process::exit(1);
    }
}
