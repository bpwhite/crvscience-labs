//! Dual-Mode CLI & GUI Entrypoint for `bio-suite-seqnet`.
//!
//! Supports:
//! 1. `run` — Live NCBI Entrez + SQLite incremental multi-locus neural network run
//!    in either `--mode headless` (console) or `--mode gui` (console + live `egui` monitor).
//! 2. `demo-sim` — Offline 2-pass Coleoptera simulation demonstrating Pass 1 cold-start
//!    quarantine/chimera convergence and Pass 2 warm-start SQLite instant consensus.
//! 3. `inspect-db` — Audits registered `database_id` partitions and soft-quarantine history.

use anyhow::Result;
use clap::{Parser, Subcommand, ValueEnum};
use std::path::PathBuf;

use bio_suite_seqnet::config::{AnalysisLevel, ConsensusConfig, ExecutionMode, OtuRank, build_locus_panel};
use bio_suite_seqnet::db::SeqNetDatabase;
use bio_suite_seqnet::entrez::generate_coleoptera_benchmark_cohort;
use bio_suite_seqnet::gui::launch_gui_with_worker;
use bio_suite_seqnet::neurons::{EngineControlFlags, run_neural_network_schedule};

#[derive(Debug, Clone, Copy, ValueEnum)]
enum ModeArg {
    Headless,
    Gui,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum RankArg {
    Population,
    Species,
    Genus,
    Family,
}

#[derive(Parser, Debug)]
#[command(
    name = "bio-suite-seqnet",
    version,
    about = "Pure-Rust 3-Stage Hybrid Sequence Analysis Neural Network (SQLite Cache + Chimeric OTU + NEXUS)"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Execute the 3-stage sequence analysis neural network on a target organism group
    Run {
        /// Organism group (e.g., "Coleoptera")
        #[arg(long, default_value = "Coleoptera")]
        taxon: String,

        /// Number of nuclear loci in target panel
        #[arg(long, default_value_t = 18)]
        nuclear_loci: usize,

        /// Number of mitochondrial loci in target panel
        #[arg(long, default_value_t = 7)]
        mito_loci: usize,

        /// Target number of retained OTU terminals
        #[arg(long, default_value_t = 100)]
        target_otus: usize,

        /// Taxonomic rank governing OTU grouping rules
        #[arg(long, value_enum, default_value_t = RankArg::Genus)]
        otu_rank: RankArg,

        /// Minimum separate species terminals required per retained Genus
        #[arg(long, default_value_t = 2)]
        min_species_per_genus: usize,

        /// Minimum per-OTU locus coverage ratio (e.g. 0.90 = >=90% of target loci represented)
        #[arg(long, default_value_t = 0.90)]
        min_coverage: f64,

        /// Maximum feedback epochs before finalizing consensus
        #[arg(long, default_value_t = 8)]
        max_epochs: usize,

        /// Path to SQLite cache database file
        #[arg(long, default_value = "coleoptera_seqnet_cache.sqlite")]
        db_path: PathBuf,

        /// Explicit logical `database_id` stored in SQLite `database_registry`
        #[arg(long, default_value = "CRV-COLEOPTERA-ML25-DB01")]
        db_id: String,

        /// Output path for clean partitioned NEXUS supermatrix
        #[arg(long, default_value = "coleoptera_25loci_optimal.nex")]
        out: PathBuf,

        /// Execution mode: `headless` (terminal only) or `gui` (terminal + desktop monitor)
        #[arg(long, value_enum, default_value_t = ModeArg::Headless)]
        mode: ModeArg,

        /// Disable live NCBI Entrez calls and operate strictly on existing SQLite cache
        #[arg(long, default_value_t = false)]
        offline_cache_only: bool,

        /// Explicit locus panel by canonical name, comma-separated (e.g. `SSU_18S,LSU_28S,COI,rrnL_16S`);
        /// overrides --nuclear-loci / --mito-loci
        #[arg(long, value_delimiter = ',')]
        loci: Vec<String>,

        /// Stratified sampling taxa, comma-separated (e.g. beetle families); one Entrez query per taxon x locus
        #[arg(long, value_delimiter = ',')]
        sampling_taxa: Vec<String>,

        /// Entrez retmax per (taxon, locus) query on a cold start
        #[arg(long, default_value_t = 400)]
        retmax: usize,

        /// Gap-refill species directives per epoch
        #[arg(long, default_value_t = 6)]
        refill_directives: usize,

        /// Minimum locus coverage for a species to receive a gap-refill directive
        #[arg(long, default_value_t = 0.55)]
        refill_min_coverage: f64,
    },

    /// Run a 2-pass offline Coleoptera benchmark (Pass 1 cold start -> Pass 2 cached delta sync)
    DemoSim {
        #[arg(long, default_value = "demo_coleoptera_seqnet.sqlite")]
        db_path: PathBuf,

        #[arg(long, default_value = "CRV-DEMO-COLEOPTERA-DB01")]
        db_id: String,

        #[arg(long, default_value = "demo_coleoptera_25loci.nex")]
        out: PathBuf,

        #[arg(long, value_enum, default_value_t = ModeArg::Headless)]
        mode: ModeArg,
    },

    /// Inspect registered `database_id` partitions and soft-quarantine ledger in a SQLite file
    InspectDb {
        #[arg(long)]
        db_path: PathBuf,

        #[arg(long, default_value = "CRV-COLEOPTERA-ML25-DB01")]
        db_id: String,

        #[arg(long, default_value = "Coleoptera")]
        taxon: String,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Run {
            taxon,
            nuclear_loci,
            mito_loci,
            target_otus,
            otu_rank,
            min_species_per_genus,
            min_coverage,
            max_epochs,
            db_path,
            db_id,
            out,
            mode,
            offline_cache_only,
            loci,
            sampling_taxa,
            retmax,
            refill_directives,
            refill_min_coverage,
        } => {
            let (nuclear_loci, mito_loci) = if loci.is_empty() {
                (nuclear_loci, mito_loci)
            } else {
                let panel = bio_suite_seqnet::config::build_named_locus_panel(&loci, true);
                let nuc = panel
                    .iter()
                    .filter(|l| l.locus_class == bio_suite_seqnet::config::LocusClass::Nuclear)
                    .count();
                (nuc, panel.len() - nuc)
            };
            let cfg = ConsensusConfig {
                locus_names: loci,
                sampling_taxa,
                retmax_per_stratum: retmax,
                max_refill_directives: refill_directives,
                refill_min_coverage,
                taxon_group: taxon,
                analysis_level: AnalysisLevel::MultiLocusPhylogenetic,
                target_loci_count: nuclear_loci + mito_loci,
                target_nuclear_loci: nuclear_loci,
                target_mito_loci: mito_loci,
                target_otu_count: target_otus,
                otu_rank: match otu_rank {
                    RankArg::Population => OtuRank::Population,
                    RankArg::Species => OtuRank::Species,
                    RankArg::Genus => OtuRank::Genus,
                    RankArg::Family => OtuRank::Family,
                },
                min_species_per_genus,
                min_locus_coverage_ratio: min_coverage,
                max_epochs,
                sqlite_path: db_path,
                database_id: db_id,
                nexus_output_path: out,
                allow_live_entrez: !offline_cache_only,
                execution_mode: match mode {
                    ModeArg::Headless => ExecutionMode::Headless,
                    ModeArg::Gui => ExecutionMode::GuiMonitor,
                },
                ..ConsensusConfig::default()
            };

            println!(
                "=== bio-suite-seqnet | Taxon: {} | DB ID: {} | Target: {} loci ({} nuc + {} mito), >= {:.0}% fill ===",
                cfg.taxon_group,
                cfg.database_id,
                cfg.target_loci_count,
                cfg.target_nuclear_loci,
                cfg.target_mito_loci,
                cfg.min_locus_coverage_ratio * 100.0
            );

            match cfg.execution_mode {
                ExecutionMode::Headless => {
                    let res = run_neural_network_schedule(
                        &cfg,
                        None,
                        None,
                        EngineControlFlags::new(),
                    )?;
                    print_run_summary(&res);
                }
                ExecutionMode::GuiMonitor => {
                    launch_gui_with_worker(cfg, None)?;
                }
            }
        }

        Commands::DemoSim {
            db_path,
            db_id,
            out,
            mode,
        } => {
            let panel = build_locus_panel(18, 7, true);
            let wave1 = generate_coleoptera_benchmark_cohort(&panel, false);

            let cfg = ConsensusConfig {
                taxon_group: "Coleoptera".to_string(),
                target_loci_count: 25,
                target_nuclear_loci: 18,
                target_mito_loci: 7,
                target_otu_count: 100,
                otu_rank: OtuRank::Genus,
                min_species_per_genus: 2,
                min_locus_coverage_ratio: 0.90,
                sqlite_path: db_path,
                database_id: db_id,
                nexus_output_path: out,
                allow_live_entrez: false,
                execution_mode: match mode {
                    ModeArg::Headless => ExecutionMode::Headless,
                    ModeArg::Gui => ExecutionMode::GuiMonitor,
                },
                ..ConsensusConfig::default()
            };

            match cfg.execution_mode {
                ExecutionMode::Headless => {
                    println!(">>> [Pass 1: Cold-Start Coleoptera 25-Locus Run]");
                    let res1 = run_neural_network_schedule(
                        &cfg,
                        Some(wave1),
                        None,
                        EngineControlFlags::new(),
                    )?;
                    print_run_summary(&res1);

                    println!("\n>>> [Pass 2: Warm-Start Cached + Incremental Delta Sequence Check]");
                    let wave2_delta = generate_coleoptera_benchmark_cohort(&panel, true);
                    let res2 = run_neural_network_schedule(
                        &cfg,
                        Some(wave2_delta),
                        None,
                        EngineControlFlags::new(),
                    )?;
                    print_run_summary(&res2);
                }
                ExecutionMode::GuiMonitor => {
                    let wave2 = generate_coleoptera_benchmark_cohort(&panel, true);
                    launch_gui_with_worker(cfg, Some(wave2))?;
                }
            }
        }

        Commands::InspectDb {
            db_path,
            db_id,
            taxon,
        } => {
            let db = SeqNetDatabase::open_or_init(&db_path, &db_id, "Inspector")?;
            println!("Registered logical databases in {:?}:", db_path);
            for (id, label, opened) in db.list_registered_databases()? {
                println!("  - {} ({}) [last opened: {}]", id, label, opened);
            }
            let qmap = db.load_quarantine_ledger(&taxon)?;
            println!(
                "\nSoft-Quarantine Ledger for `{}` (database_id={}): {} flagged entries",
                taxon,
                db_id,
                qmap.len()
            );
            for q in qmap.values() {
                println!(
                    "  * {} | {} | {} | quarantined={} | strikes={} | skep_prior={:.3} | reason={}",
                    q.accession_version,
                    q.species_binomial,
                    q.locus_canonical,
                    q.currently_quarantined,
                    q.prior_quarantine_count,
                    q.skepticism_prior,
                    q.last_reason_code
                );
            }
        }
    }

    Ok(())
}

fn print_run_summary(res: &bio_suite_seqnet::neurons::SeqNetRunResult) {
    println!(
        "=== Summary [Run {} | DB {}] ===",
        res.audit.run_id, res.audit.database_id
    );
    println!(
        "  Converged: {} in {} epoch(s) | Consensus: {:.4} | Matrix Fill: {:.2}%",
        res.audit.converged,
        res.audit.epochs_executed,
        res.audit.final_consensus_score,
        res.audit.final_matrix_occupancy * 100.0
    );
    println!(
        "  Retained Species Terminals: {} | Soft-Quarantined (Excluded from NEXUS): {} | Rehabilitated: {}",
        res.otus.len(),
        res.audit.quarantined_count,
        res.audit.rehabilitated_count
    );
    println!("  Clean NEXUS Output: {}", res.audit.nexus_output_path);
}
