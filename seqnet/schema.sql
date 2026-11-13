-- ============================================================================
-- bio-suite-seqnet: Multi-Database SQLite Schema (GPL-3.0)
-- Supports isolated/portable SQLite caches keyed by `database_id`, incremental
-- NCBI Entrez `mdat` watermarks, non-destructive soft-quarantine history,
-- persistent hybrid neural weights, and chimeric OTU accession provenance.
-- ============================================================================

PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS database_registry (
    database_id       TEXT PRIMARY KEY,
    label             TEXT NOT NULL,
    schema_version    INTEGER NOT NULL DEFAULT 1,
    created_at        TEXT NOT NULL,
    last_opened_at    TEXT NOT NULL,
    description       TEXT NOT NULL DEFAULT ''
);

CREATE TABLE IF NOT EXISTS entrez_sync_watermarks (
    database_id            TEXT NOT NULL,
    taxon_group            TEXT NOT NULL,
    locus_canonical        TEXT NOT NULL,
    last_sync_mdat         TEXT NOT NULL, -- YYYY/MM/DD format for NCBI Entrez mindate
    cached_accession_count INTEGER NOT NULL DEFAULT 0,
    last_run_id            TEXT NOT NULL,
    updated_at             TEXT NOT NULL,
    PRIMARY KEY (database_id, taxon_group, locus_canonical),
    FOREIGN KEY (database_id) REFERENCES database_registry(database_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS sequences_cache (
    database_id       TEXT NOT NULL,
    accession_version TEXT NOT NULL,
    ncbi_taxid        INTEGER NOT NULL DEFAULT 0,
    organism_name     TEXT NOT NULL,
    genus             TEXT NOT NULL,
    species_binomial  TEXT NOT NULL,
    family            TEXT NOT NULL DEFAULT '',
    voucher_id        TEXT NOT NULL DEFAULT '',
    isolate_code      TEXT NOT NULL DEFAULT '',
    country_locality  TEXT NOT NULL DEFAULT '',
    locus_canonical   TEXT NOT NULL,
    locus_class       TEXT NOT NULL CHECK (locus_class IN ('nuclear', 'mitochondrial')),
    sequence_raw      TEXT NOT NULL,
    seq_length        INTEGER NOT NULL,
    gc_fraction       REAL NOT NULL,
    ambig_fraction    REAL NOT NULL,
    shannon_entropy   REAL NOT NULL,
    reading_frame_ok  INTEGER NOT NULL DEFAULT 1,
    ncbi_mod_date     TEXT NOT NULL,
    fetched_at        TEXT NOT NULL,
    PRIMARY KEY (database_id, accession_version, locus_canonical),
    FOREIGN KEY (database_id) REFERENCES database_registry(database_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_seq_taxon_locus
    ON sequences_cache(database_id, genus, species_binomial, locus_canonical);

CREATE INDEX IF NOT EXISTS idx_seq_locus
    ON sequences_cache(database_id, locus_canonical);

-- Soft-quarantine ledger: sequences are NEVER permanently deleted.
-- Flagging `prior_quarantine_count` applies a Bayesian skepticism prior on future
-- runs while permitting rehabilitation if denser taxon sampling resolves an outlier.
CREATE TABLE IF NOT EXISTS quarantine_ledger (
    database_id             TEXT NOT NULL,
    accession_version       TEXT NOT NULL,
    locus_canonical         TEXT NOT NULL,
    taxon_group             TEXT NOT NULL,
    genus                   TEXT NOT NULL,
    species_binomial        TEXT NOT NULL,
    currently_quarantined   INTEGER NOT NULL CHECK (currently_quarantined IN (0, 1)),
    prior_quarantine_count  INTEGER NOT NULL DEFAULT 0,
    prior_pass_count        INTEGER NOT NULL DEFAULT 0,
    skepticism_prior        REAL NOT NULL DEFAULT 0.0,
    last_reason_code        TEXT NOT NULL,
    within_locus_k2p_dist   REAL NOT NULL DEFAULT 0.0,
    neural_consensus_score  REAL NOT NULL DEFAULT 0.0,
    last_evaluated_run_id   TEXT NOT NULL,
    last_evaluated_at       TEXT NOT NULL,
    audit_notes             TEXT NOT NULL DEFAULT '',
    PRIMARY KEY (database_id, accession_version, locus_canonical),
    FOREIGN KEY (database_id) REFERENCES database_registry(database_id) ON DELETE CASCADE
);

CREATE INDEX IF NOT EXISTS idx_quarantine_group
    ON quarantine_ledger(database_id, taxon_group, currently_quarantined);

-- Stores learned continuous neural consensus weights per (database_id, taxon_group)
-- so subsequent incremental runs warm-start near optimal consensus immediately.
CREATE TABLE IF NOT EXISTS neural_weights_store (
    database_id       TEXT NOT NULL,
    taxon_group       TEXT NOT NULL,
    weights_json      TEXT NOT NULL,
    epochs_trained    INTEGER NOT NULL DEFAULT 0,
    last_loss         REAL NOT NULL DEFAULT 0.0,
    updated_at        TEXT NOT NULL,
    PRIMARY KEY (database_id, taxon_group),
    FOREIGN KEY (database_id) REFERENCES database_registry(database_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS neural_runs (
    run_id                   TEXT PRIMARY KEY,
    database_id              TEXT NOT NULL,
    taxon_group              TEXT NOT NULL,
    analysis_level           TEXT NOT NULL,
    target_loci_count        INTEGER NOT NULL,
    target_nuclear_loci      INTEGER NOT NULL,
    target_mito_loci         INTEGER NOT NULL,
    target_otu_count         INTEGER NOT NULL,
    otu_rank_level           TEXT NOT NULL,
    min_species_per_genus    INTEGER NOT NULL,
    min_locus_coverage_ratio REAL NOT NULL,
    epochs_executed          INTEGER NOT NULL,
    cache_hits               INTEGER NOT NULL DEFAULT 0,
    live_delta_fetched       INTEGER NOT NULL DEFAULT 0,
    quarantined_count        INTEGER NOT NULL DEFAULT 0,
    rehabilitated_count      INTEGER NOT NULL DEFAULT 0,
    converged                INTEGER NOT NULL CHECK (converged IN (0, 1)),
    final_consensus_score    REAL NOT NULL,
    final_matrix_occupancy   REAL NOT NULL,
    nexus_output_path        TEXT NOT NULL,
    started_at               TEXT NOT NULL,
    completed_at             TEXT NOT NULL,
    FOREIGN KEY (database_id) REFERENCES database_registry(database_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS artificial_otus (
    otu_id                   TEXT PRIMARY KEY,
    run_id                   TEXT NOT NULL,
    database_id              TEXT NOT NULL,
    terminal_label           TEXT NOT NULL,
    genus                    TEXT NOT NULL,
    species_binomial         TEXT NOT NULL,
    family                   TEXT NOT NULL DEFAULT '',
    is_chimeric              INTEGER NOT NULL CHECK (is_chimeric IN (0, 1)),
    distinct_specimen_count  INTEGER NOT NULL,
    distinct_accession_count INTEGER NOT NULL,
    loci_filled              INTEGER NOT NULL,
    loci_target              INTEGER NOT NULL,
    locus_coverage_ratio     REAL NOT NULL,
    mean_neural_confidence   REAL NOT NULL,
    created_at               TEXT NOT NULL,
    FOREIGN KEY (run_id) REFERENCES neural_runs(run_id) ON DELETE CASCADE,
    FOREIGN KEY (database_id) REFERENCES database_registry(database_id) ON DELETE CASCADE
);

CREATE TABLE IF NOT EXISTS artificial_otu_accessions (
    otu_id                TEXT NOT NULL,
    locus_canonical       TEXT NOT NULL,
    locus_class           TEXT NOT NULL,
    accession_version     TEXT NOT NULL,
    voucher_id            TEXT NOT NULL DEFAULT '',
    organism_name         TEXT NOT NULL,
    raw_seq_length        INTEGER NOT NULL,
    aligned_seq_length    INTEGER NOT NULL,
    neural_confidence     REAL NOT NULL,
    within_locus_k2p_dist REAL NOT NULL,
    aligned_sequence      TEXT NOT NULL,
    PRIMARY KEY (otu_id, locus_canonical),
    FOREIGN KEY (otu_id) REFERENCES artificial_otus(otu_id) ON DELETE CASCADE
);
