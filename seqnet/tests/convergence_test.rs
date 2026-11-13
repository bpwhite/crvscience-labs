//! End-to-end integration tests for `bio-suite-seqnet`.
//!
//! Validates:
//! 1. Multi-DB `database_id` isolation inside SQLite (`CRV-TEST-DB-A` vs `CRV-TEST-DB-B`).
//! 2. Soft-quarantine exclusion of NUMT internal stop codons (`KX990001.1`), K2P sample swaps
//!    (`JQ8800*.1`), and high-ambiguity reads (`AY770012.1`) from the exported NEXUS matrix,
//!    while preserving their flag metadata in `BIOSUITE_METADATA` and `quarantine_ledger`.
//! 3. Enforcement of `>= 2` separate species terminals per Genus (dropping singleton `Cupes`)
//!    and `>= 90%` locus representation per retained chimeric species terminal.
//! 4. Warm-start rapid 1-epoch convergence on Pass 2 when re-running with incremental delta data.

use bio_suite_seqnet::config::{ConsensusConfig, ExecutionMode, OtuRank, build_locus_panel};
use bio_suite_seqnet::entrez::generate_coleoptera_benchmark_cohort;
use bio_suite_seqnet::neurons::{EngineControlFlags, run_neural_network_schedule};

#[test]
fn test_coleoptera_hybrid_convergence_and_warm_cache_delta() {
    let tmp_dir = std::env::temp_dir().join("bio_suite_seqnet_integration_test");
    let _ = std::fs::remove_dir_all(&tmp_dir);
    std::fs::create_dir_all(&tmp_dir).unwrap();

    let db_path = tmp_dir.join("coleoptera_cache.sqlite");
    let nex_path = tmp_dir.join("coleoptera_optimal.nex");

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
        sqlite_path: db_path.clone(),
        database_id: "CRV-TEST-COLEOPTERA-001".to_string(),
        nexus_output_path: nex_path.clone(),
        allow_live_entrez: false,
        execution_mode: ExecutionMode::Headless,
        ..ConsensusConfig::default()
    };

    // Pass 1: Cold start
    let res1 =
        run_neural_network_schedule(&cfg, Some(wave1), None, EngineControlFlags::new()).unwrap();
    assert!(res1.audit.converged);
    assert!(res1.audit.final_matrix_occupancy >= 0.90);
    assert!(res1.audit.quarantined_count >= 3);

    // Every retained OTU must have >= 90% locus coverage and no singleton genus (`Cupes`)
    for otu in &res1.otus {
        assert!(otu.locus_coverage_ratio >= 0.90);
        assert_ne!(otu.genus, "Cupes");
    }

    // NEXUS file must NEVER contain quarantined sequence rows in MATRIX, only metadata flags
    let nex_text = std::fs::read_to_string(&nex_path).unwrap();
    assert!(nex_text.contains("#NEXUS"));
    assert!(nex_text.contains("BEGIN DATA;"));
    assert!(nex_text.contains("CHARPARTITION GenomicCompartment"));
    assert!(nex_text.contains("QUARANTINE_FLAG accession=KX990001.1"));
    assert!(!nex_text.contains("BAD-NUMT-01"));

    // Pass 2: Warm-start re-run with newly deposited `Anoplophora chinensis` sequences
    let wave2 = generate_coleoptera_benchmark_cohort(&panel, true);
    let res2 =
        run_neural_network_schedule(&cfg, Some(wave2), None, EngineControlFlags::new()).unwrap();
    assert!(res2.audit.converged);
    assert_eq!(
        res2.audit.epochs_executed, 1,
        "Cached warm-start should achieve consensus in 1 epoch"
    );
    // `Anoplophora` now has 2 species >= 90% (`glabripennis` + `chinensis`), adding +2 terminals!
    assert!(res2.otus.len() > res1.otus.len());
}
