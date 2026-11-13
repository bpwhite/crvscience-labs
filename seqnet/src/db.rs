//! Multi-database SQLite persistence engine (`database_id`-scoped), incremental Entrez
//! modification-date (`mdat`) watermarks, non-destructive soft-quarantine ledger,
//! persistent neural weight store, and chimeric OTU provenance recording.

use anyhow::{bail, Context, Result};
use chrono::Utc;
use rusqlite::{params, Connection, OptionalExtension};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

use crate::config::{ConsensusConfig, LocusClass};

const SCHEMA_SQL: &str = include_str!("../schema.sql");

/// Raw or cached nucleotide sequence record with extracted taxonomic and voucher metadata.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceRecord {
    pub accession_version: String,
    pub ncbi_taxid: u64,
    pub organism_name: String,
    pub genus: String,
    pub species_binomial: String,
    pub family: String,
    pub voucher_id: String,
    pub isolate_code: String,
    pub country_locality: String,
    pub locus_canonical: String,
    pub locus_class: LocusClass,
    pub sequence_raw: String,
    pub seq_length: usize,
    pub gc_fraction: f64,
    pub ambig_fraction: f64,
    pub shannon_entropy: f64,
    pub reading_frame_ok: bool,
    /// Modification date in `YYYY/MM/DD` format for NCBI Entrez `mdat` watermarks.
    pub ncbi_mod_date: String,
}

/// Historical soft-quarantine state stored in `quarantine_ledger` (never permanently deleted).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QuarantineEntry {
    pub accession_version: String,
    pub locus_canonical: String,
    pub taxon_group: String,
    pub genus: String,
    pub species_binomial: String,
    pub currently_quarantined: bool,
    pub prior_quarantine_count: u32,
    pub prior_pass_count: u32,
    pub skepticism_prior: f64,
    pub last_reason_code: String,
    pub within_locus_k2p_dist: f64,
    pub neural_consensus_score: f64,
    pub last_evaluated_run_id: String,
    pub audit_notes: String,
}

/// Summary statistics recorded at the end of a neural network convergence run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunAuditRecord {
    pub run_id: String,
    pub database_id: String,
    pub taxon_group: String,
    pub epochs_executed: usize,
    pub cache_hits: usize,
    pub live_delta_fetched: usize,
    pub quarantined_count: usize,
    pub rehabilitated_count: usize,
    pub converged: bool,
    pub final_consensus_score: f64,
    pub final_matrix_occupancy: f64,
    pub nexus_output_path: String,
    pub started_at: String,
    pub completed_at: String,
}

/// Persistent SQLite repository bound to a verified `database_id`.
pub struct SeqNetDatabase {
    conn: Connection,
    pub database_id: String,
}

impl SeqNetDatabase {
    /// Opens or initializes a SQLite database file and enforces `database_id` identity
    /// so multiple distinct taxonomic caches never cross-contaminate by accident.
    pub fn open_or_init(path: &Path, expected_db_id: &str, label: &str) -> Result<Self> {
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("failed to create directory {:?}", parent))?;
            }
        }
        let conn = Connection::open(path)
            .with_context(|| format!("failed to open SQLite DB at {:?}", path))?;

        conn.execute_batch(
            "PRAGMA journal_mode = WAL;
             PRAGMA synchronous = NORMAL;
             PRAGMA foreign_keys = ON;",
        )?;
        conn.execute_batch(SCHEMA_SQL)
            .context("failed to apply bio-suite-seqnet SQLite schema")?;

        let now = Utc::now().to_rfc3339();

        // Check whether this SQLite file already hosts a registered `database_id`.
        let existing_ids: Vec<String> = {
            let mut stmt = conn.prepare("SELECT database_id FROM database_registry")?;
            let rows = stmt.query_map([], |r| r.get::<_, String>(0))?;
            rows.collect::<std::result::Result<Vec<_>, _>>()?
        };

        if existing_ids.is_empty() {
            conn.execute(
                "INSERT INTO database_registry
                 (database_id, label, schema_version, created_at, last_opened_at, description)
                 VALUES (?1, ?2, 1, ?3, ?3, ?4)",
                params![
                    expected_db_id,
                    label,
                    now,
                    format!("bio-suite-seqnet multi-locus cache [{}]", expected_db_id)
                ],
            )?;
        } else if !existing_ids.iter().any(|id| id == expected_db_id) {
            // Multiple logical database IDs are explicitly supported inside either separate
            // SQLite files or partitioned inside a shared SQLite file. Register the new ID cleanly.
            conn.execute(
                "INSERT INTO database_registry
                 (database_id, label, schema_version, created_at, last_opened_at, description)
                 VALUES (?1, ?2, 1, ?3, ?3, ?4)",
                params![
                    expected_db_id,
                    label,
                    now,
                    format!("Additional logical partition [{}]", expected_db_id)
                ],
            )?;
        } else {
            conn.execute(
                "UPDATE database_registry SET last_opened_at = ?1 WHERE database_id = ?2",
                params![now, expected_db_id],
            )?;
        }

        Ok(Self {
            conn,
            database_id: expected_db_id.to_string(),
        })
    }

    /// Lists all logical `database_id` partitions registered in this SQLite file.
    pub fn list_registered_databases(&self) -> Result<Vec<(String, String, String)>> {
        let mut stmt = self.conn.prepare(
            "SELECT database_id, label, last_opened_at
             FROM database_registry ORDER BY created_at ASC",
        )?;
        let rows = stmt.query_map([], |r| Ok((r.get(0)?, r.get(1)?, r.get(2)?)))?;
        rows.collect::<std::result::Result<Vec<_>, _>>()
            .map_err(Into::into)
    }

    /// Retrieves the last NCBI Entrez modification-date (`mdat`) watermark for incremental queries.
    pub fn get_sync_watermark(&self, taxon_group: &str, locus: &str) -> Result<Option<String>> {
        self.conn
            .query_row(
                "SELECT last_sync_mdat FROM entrez_sync_watermarks
                 WHERE database_id = ?1 AND taxon_group = ?2 AND locus_canonical = ?3",
                params![self.database_id, taxon_group, locus],
                |r| r.get::<_, String>(0),
            )
            .optional()
            .map_err(Into::into)
    }

    /// Updates the incremental Entrez `mdat` watermark after live querying completes.
    pub fn upsert_sync_watermark(
        &self,
        taxon_group: &str,
        locus: &str,
        mdat_yyyy_mm_dd: &str,
        cached_count: usize,
        run_id: &str,
    ) -> Result<()> {
        let now = Utc::now().to_rfc3339();
        self.conn.execute(
            "INSERT INTO entrez_sync_watermarks
             (database_id, taxon_group, locus_canonical, last_sync_mdat, cached_accession_count, last_run_id, updated_at)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7)
             ON CONFLICT(database_id, taxon_group, locus_canonical) DO UPDATE SET
                 last_sync_mdat = excluded.last_sync_mdat,
                 cached_accession_count = excluded.cached_accession_count,
                 last_run_id = excluded.last_run_id,
                 updated_at = excluded.updated_at",
            params![
                self.database_id,
                taxon_group,
                locus,
                mdat_yyyy_mm_dd,
                cached_count as i64,
                run_id,
                now
            ],
        )?;
        Ok(())
    }

    /// Bulk-inserts or updates newly fetched sequence records inside a single WAL transaction.
    pub fn upsert_sequences(&mut self, records: &[SequenceRecord]) -> Result<usize> {
        if records.is_empty() {
            return Ok(0);
        }
        let now = Utc::now().to_rfc3339();
        let tx = self.conn.transaction()?;
        {
            let mut stmt = tx.prepare_cached(
                "INSERT INTO sequences_cache (
                    database_id, accession_version, ncbi_taxid, organism_name, genus,
                    species_binomial, family, voucher_id, isolate_code, country_locality,
                    locus_canonical, locus_class, sequence_raw, seq_length, gc_fraction,
                    ambig_fraction, shannon_entropy, reading_frame_ok, ncbi_mod_date, fetched_at
                 ) VALUES (
                    ?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10,
                    ?11, ?12, ?13, ?14, ?15, ?16, ?17, ?18, ?19, ?20
                 )
                 ON CONFLICT(database_id, accession_version, locus_canonical) DO UPDATE SET
                    ncbi_taxid = excluded.ncbi_taxid,
                    organism_name = excluded.organism_name,
                    genus = excluded.genus,
                    species_binomial = excluded.species_binomial,
                    family = excluded.family,
                    voucher_id = excluded.voucher_id,
                    isolate_code = excluded.isolate_code,
                    country_locality = excluded.country_locality,
                    sequence_raw = excluded.sequence_raw,
                    seq_length = excluded.seq_length,
                    gc_fraction = excluded.gc_fraction,
                    ambig_fraction = excluded.ambig_fraction,
                    shannon_entropy = excluded.shannon_entropy,
                    reading_frame_ok = excluded.reading_frame_ok,
                    ncbi_mod_date = excluded.ncbi_mod_date,
                    fetched_at = excluded.fetched_at",
            )?;

            for rec in records {
                stmt.execute(params![
                    self.database_id,
                    rec.accession_version,
                    rec.ncbi_taxid as i64,
                    rec.organism_name,
                    rec.genus,
                    rec.species_binomial,
                    rec.family,
                    rec.voucher_id,
                    rec.isolate_code,
                    rec.country_locality,
                    rec.locus_canonical,
                    rec.locus_class.as_str(),
                    rec.sequence_raw,
                    rec.seq_length as i64,
                    rec.gc_fraction,
                    rec.ambig_fraction,
                    rec.shannon_entropy,
                    if rec.reading_frame_ok { 1 } else { 0 },
                    rec.ncbi_mod_date,
                    now,
                ])?;
            }
        }
        tx.commit()?;
        Ok(records.len())
    }

    /// Loads all cached sequences matching the target locus panel for warm-start processing.
    pub fn load_cached_sequences(&self, loci: &[String]) -> Result<Vec<SequenceRecord>> {
        if loci.is_empty() {
            return Ok(Vec::new());
        }
        let mut stmt = self.conn.prepare(
            "SELECT accession_version, ncbi_taxid, organism_name, genus, species_binomial,
                    family, voucher_id, isolate_code, country_locality, locus_canonical,
                    locus_class, sequence_raw, seq_length, gc_fraction, ambig_fraction,
                    shannon_entropy, reading_frame_ok, ncbi_mod_date
             FROM sequences_cache
             WHERE database_id = ?1
             ORDER BY genus ASC, species_binomial ASC, locus_canonical ASC",
        )?;
        let locus_set: std::collections::HashSet<&str> = loci.iter().map(|s| s.as_str()).collect();
        let rows = stmt.query_map(params![self.database_id], |r| {
            let lclass_str: String = r.get(10)?;
            let lclass = if lclass_str == "mitochondrial" {
                LocusClass::Mitochondrial
            } else {
                LocusClass::Nuclear
            };
            Ok(SequenceRecord {
                accession_version: r.get(0)?,
                ncbi_taxid: r.get::<_, i64>(1)? as u64,
                organism_name: r.get(2)?,
                genus: r.get(3)?,
                species_binomial: r.get(4)?,
                family: r.get(5)?,
                voucher_id: r.get(6)?,
                isolate_code: r.get(7)?,
                country_locality: r.get(8)?,
                locus_canonical: r.get(9)?,
                locus_class: lclass,
                sequence_raw: r.get(11)?,
                seq_length: r.get::<_, i64>(12)? as usize,
                gc_fraction: r.get(13)?,
                ambig_fraction: r.get(14)?,
                shannon_entropy: r.get(15)?,
                reading_frame_ok: r.get::<_, i64>(16)? != 0,
                ncbi_mod_date: r.get(17)?,
            })
        })?;

        let mut out = Vec::new();
        for row in rows {
            let rec = row?;
            if locus_set.contains(rec.locus_canonical.as_str()) {
                out.push(rec);
            }
        }
        Ok(out)
    }

    /// Loads the historical soft-quarantine ledger keyed by `(accession_version, locus_canonical)`.
    pub fn load_quarantine_ledger(
        &self,
        taxon_group: &str,
    ) -> Result<HashMap<(String, String), QuarantineEntry>> {
        let mut stmt = self.conn.prepare(
            "SELECT accession_version, locus_canonical, taxon_group, genus, species_binomial,
                    currently_quarantined, prior_quarantine_count, prior_pass_count,
                    skepticism_prior, last_reason_code, within_locus_k2p_dist,
                    neural_consensus_score, last_evaluated_run_id, audit_notes
             FROM quarantine_ledger
             WHERE database_id = ?1 AND taxon_group = ?2",
        )?;
        let rows = stmt.query_map(params![self.database_id, taxon_group], |r| {
            Ok(QuarantineEntry {
                accession_version: r.get(0)?,
                locus_canonical: r.get(1)?,
                taxon_group: r.get(2)?,
                genus: r.get(3)?,
                species_binomial: r.get(4)?,
                currently_quarantined: r.get::<_, i64>(5)? != 0,
                prior_quarantine_count: r.get::<_, i64>(6)? as u32,
                prior_pass_count: r.get::<_, i64>(7)? as u32,
                skepticism_prior: r.get(8)?,
                last_reason_code: r.get(9)?,
                within_locus_k2p_dist: r.get(10)?,
                neural_consensus_score: r.get(11)?,
                last_evaluated_run_id: r.get(12)?,
                audit_notes: r.get(13)?,
            })
        })?;

        let mut map = HashMap::new();
        for row in rows {
            let entry = row?;
            map.insert(
                (
                    entry.accession_version.clone(),
                    entry.locus_canonical.clone(),
                ),
                entry,
            );
        }
        Ok(map)
    }

    /// Persists updated soft-quarantine states at the end of a run. Previously quarantined
    /// accessions remain in SQLite with updated `prior_quarantine_count` / `skepticism_prior`
    /// rather than being deleted.
    pub fn save_quarantine_updates(&mut self, entries: &[QuarantineEntry]) -> Result<()> {
        let now = Utc::now().to_rfc3339();
        let tx = self.conn.transaction()?;
        {
            let mut stmt = tx.prepare_cached(
                "INSERT INTO quarantine_ledger (
                    database_id, accession_version, locus_canonical, taxon_group,
                    genus, species_binomial, currently_quarantined, prior_quarantine_count,
                    prior_pass_count, skepticism_prior, last_reason_code,
                    within_locus_k2p_dist, neural_consensus_score, last_evaluated_run_id,
                    last_evaluated_at, audit_notes
                 ) VALUES (
                    ?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15, ?16
                 )
                 ON CONFLICT(database_id, accession_version, locus_canonical) DO UPDATE SET
                    currently_quarantined = excluded.currently_quarantined,
                    prior_quarantine_count = excluded.prior_quarantine_count,
                    prior_pass_count = excluded.prior_pass_count,
                    skepticism_prior = excluded.skepticism_prior,
                    last_reason_code = excluded.last_reason_code,
                    within_locus_k2p_dist = excluded.within_locus_k2p_dist,
                    neural_consensus_score = excluded.neural_consensus_score,
                    last_evaluated_run_id = excluded.last_evaluated_run_id,
                    last_evaluated_at = excluded.last_evaluated_at,
                    audit_notes = excluded.audit_notes",
            )?;

            for e in entries {
                stmt.execute(params![
                    self.database_id,
                    e.accession_version,
                    e.locus_canonical,
                    e.taxon_group,
                    e.genus,
                    e.species_binomial,
                    if e.currently_quarantined { 1 } else { 0 },
                    e.prior_quarantine_count as i64,
                    e.prior_pass_count as i64,
                    e.skepticism_prior,
                    e.last_reason_code,
                    e.within_locus_k2p_dist,
                    e.neural_consensus_score,
                    e.last_evaluated_run_id,
                    now,
                    e.audit_notes,
                ])?;
            }
        }
        tx.commit()?;
        Ok(())
    }

    /// Loads serialized hybrid neural network weights for warm-start consensus scoring.
    pub fn load_neural_weights(&self, taxon_group: &str) -> Result<Option<(String, usize)>> {
        self.conn
            .query_row(
                "SELECT weights_json, epochs_trained FROM neural_weights_store
                 WHERE database_id = ?1 AND taxon_group = ?2",
                params![self.database_id, taxon_group],
                |r| Ok((r.get::<_, String>(0)?, r.get::<_, i64>(1)? as usize)),
            )
            .optional()
            .map_err(Into::into)
    }

    /// Saves updated hybrid neural network weights after feedback epochs complete.
    pub fn save_neural_weights(
        &self,
        taxon_group: &str,
        weights_json: &str,
        epochs_trained: usize,
        last_loss: f64,
    ) -> Result<()> {
        let now = Utc::now().to_rfc3339();
        self.conn.execute(
            "INSERT INTO neural_weights_store
             (database_id, taxon_group, weights_json, epochs_trained, last_loss, updated_at)
             VALUES (?1, ?2, ?3, ?4, ?5, ?6)
             ON CONFLICT(database_id, taxon_group) DO UPDATE SET
                 weights_json = excluded.weights_json,
                 epochs_trained = neural_weights_store.epochs_trained + excluded.epochs_trained,
                 last_loss = excluded.last_loss,
                 updated_at = excluded.updated_at",
            params![
                self.database_id,
                taxon_group,
                weights_json,
                epochs_trained as i64,
                last_loss,
                now
            ],
        )?;
        Ok(())
    }

    /// Records the final run metadata and every artificial/chimeric OTU + constituent NCBI accession.
    pub fn save_run_and_otus(
        &mut self,
        cfg: &ConsensusConfig,
        audit: &RunAuditRecord,
        otus: &[crate::otu::AssembledSpeciesOtu],
    ) -> Result<()> {
        if audit.database_id != self.database_id {
            bail!(
                "database_id mismatch: run has {}, connection bound to {}",
                audit.database_id,
                self.database_id
            );
        }
        let now = Utc::now().to_rfc3339();
        let tx = self.conn.transaction()?;

        tx.execute(
            "INSERT INTO neural_runs (
                run_id, database_id, taxon_group, analysis_level, target_loci_count,
                target_nuclear_loci, target_mito_loci, target_otu_count, otu_rank_level,
                min_species_per_genus, min_locus_coverage_ratio, epochs_executed,
                cache_hits, live_delta_fetched, quarantined_count, rehabilitated_count,
                converged, final_consensus_score, final_matrix_occupancy,
                nexus_output_path, started_at, completed_at
             ) VALUES (
                ?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12,
                ?13, ?14, ?15, ?16, ?17, ?18, ?19, ?20, ?21, ?22
             )",
            params![
                audit.run_id,
                self.database_id,
                cfg.taxon_group,
                cfg.analysis_level.as_str(),
                cfg.target_loci_count as i64,
                cfg.target_nuclear_loci as i64,
                cfg.target_mito_loci as i64,
                cfg.target_otu_count as i64,
                cfg.otu_rank.as_str(),
                cfg.min_species_per_genus as i64,
                cfg.min_locus_coverage_ratio,
                audit.epochs_executed as i64,
                audit.cache_hits as i64,
                audit.live_delta_fetched as i64,
                audit.quarantined_count as i64,
                audit.rehabilitated_count as i64,
                if audit.converged { 1 } else { 0 },
                audit.final_consensus_score,
                audit.final_matrix_occupancy,
                audit.nexus_output_path,
                audit.started_at,
                audit.completed_at,
            ],
        )?;

        {
            let mut otu_stmt = tx.prepare_cached(
                "INSERT INTO artificial_otus (
                    otu_id, run_id, database_id, terminal_label, genus, species_binomial,
                    family, is_chimeric, distinct_specimen_count, distinct_accession_count,
                    loci_filled, loci_target, locus_coverage_ratio, mean_neural_confidence,
                    created_at
                 ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15)",
            )?;

            let mut acc_stmt = tx.prepare_cached(
                "INSERT INTO artificial_otu_accessions (
                    otu_id, locus_canonical, locus_class, accession_version, voucher_id,
                    organism_name, raw_seq_length, aligned_seq_length, neural_confidence,
                    within_locus_k2p_dist, aligned_sequence
                 ) VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11)",
            )?;

            for (idx, otu) in otus.iter().enumerate() {
                let otu_id = format!("{}_OTU_{:04}", audit.run_id, idx + 1);
                otu_stmt.execute(params![
                    otu_id,
                    audit.run_id,
                    self.database_id,
                    otu.terminal_label,
                    otu.genus,
                    otu.species_binomial,
                    otu.family,
                    if otu.is_chimeric { 1 } else { 0 },
                    otu.distinct_specimen_count as i64,
                    otu.distinct_accession_count as i64,
                    otu.loci_filled as i64,
                    otu.loci_target as i64,
                    otu.locus_coverage_ratio,
                    otu.mean_neural_confidence,
                    now,
                ])?;

                for slot in otu.locus_slots.values() {
                    let aligned_seq = slot
                        .aligned_sequence
                        .as_deref()
                        .unwrap_or(&slot.sequence_raw);
                    acc_stmt.execute(params![
                        otu_id,
                        slot.locus_canonical,
                        slot.locus_class.as_str(),
                        slot.accession_version,
                        slot.voucher_id,
                        slot.organism_name,
                        slot.sequence_raw.len() as i64,
                        aligned_seq.len() as i64,
                        slot.neural_confidence,
                        slot.within_locus_k2p_dist,
                        aligned_seq,
                    ])?;
                }
            }
        }

        tx.commit()?;
        Ok(())
    }
}
