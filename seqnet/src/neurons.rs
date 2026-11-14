//! Core 3-Stage Hybrid Sequence Analysis Neural Network (`neurons.rs`).
//!
//! Orchestrates cyclic message-passing and feedback loops across:
//! - **Stage 1 (Input Neurons)**: SQLite warm-start cache loader, NCBI Entrez `mdat` delta
//!   fetcher, and targeted gap-refill fetcher driven by Stage 3 feedback directives.
//! - **Stage 2 (Hybrid Processing Neurons)**: Deterministic bioinformatic quarantine gates
//!   (taxonomy/voucher parser, ambiguity/length filter, CDS reading-frame & NUMT stop-codon
//!   detector, parallel 4-mer cosine + Kimura 2-Parameter within-locus outlier/swap detector,
//!   non-destructive soft-quarantine re-evaluator) coupled with `HybridConsensusWeights`
//!   continuous neural scoring and online contrastive weight adaptation.
//! - **Stage 3 (Output & Valuation Neurons)**: Species-terminal-per-genus chimeric OTU
//!   assembler (`>= min_species_per_genus`, `>= 90%` per-OTU locus occupancy), between-locus
//!   coverage evaluator, parallel Gotoh MSA, and clean partitioned NEXUS supermatrix exporter.

use anyhow::Result;
use chrono::Utc;
use crossbeam_channel::Sender;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

use crate::align_nexus::{align_all_loci_parallel, write_partitioned_nexus, LocusPartitionMeta};
use crate::config::{ConsensusConfig, LocusClass, LocusSpec};

/// Sync-watermark key for the per-stratum mitogenome query (mitogenome mode).
const MITOGENOME_WATERMARK: &str = "__mitogenome__";
use crate::db::{QuarantineEntry, RunAuditRecord, SeqNetDatabase, SequenceRecord};
use crate::entrez::EntrezClient;
use crate::neural_tensor::{HybridConsensusWeights, SeqFeatureVec};
use crate::otu::{assemble_and_value_otus, AssembledSpeciesOtu, GapRefillDirective};

/// Candidate sequence passing Stage 2 Processing Neurons into Stage 3 OTU Synthesis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidatedCandidate {
    pub record: SequenceRecord,
    pub features: SeqFeatureVec,
    pub within_locus_k2p_dist: f64,
    pub kmer_cosine_sim: f64,
    pub neural_confidence: f64,
    pub rehabilitated_from_quarantine: bool,
}

/// Interactive thread-safe control flags shared between the GUI Monitor and the neural engine.
#[derive(Debug, Default)]
pub struct EngineControlFlags {
    pub paused: AtomicBool,
    pub force_finalize: AtomicBool,
}

impl EngineControlFlags {
    pub fn new() -> Arc<Self> {
        Arc::new(Self::default())
    }
}

/// Real-time telemetry snapshot streamed to both the console logger and the `egui` monitor.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NeuronTelemetrySnapshot {
    pub database_id: String,
    pub run_id: String,
    pub taxon_group: String,
    pub epoch: usize,
    pub max_epochs: usize,
    pub active_stage: String,
    pub cache_hits: usize,
    pub live_delta_fetched: usize,
    pub total_sequences_evaluated: usize,
    pub clean_candidates_count: usize,
    pub quarantined_count: usize,
    pub rehabilitated_count: usize,
    pub retained_otus_count: usize,
    pub retained_genera_count: usize,
    pub chimeric_otus_count: usize,
    pub matrix_occupancy: f64,
    pub global_consensus_score: f64,
    pub neural_bce_loss: f64,
    pub converged: bool,
    pub finished: bool,
    pub locus_names: Vec<String>,
    pub retained_otus_preview: Vec<AssembledSpeciesOtu>,
    pub quarantine_preview: Vec<QuarantineEntry>,
    pub feedback_log: Vec<String>,
}

/// Complete final result returned when the neural network finishes its schedule.
#[derive(Debug, Clone)]
pub struct SeqNetRunResult {
    pub audit: RunAuditRecord,
    pub otus: Vec<AssembledSpeciesOtu>,
    pub partitions: Vec<LocusPartitionMeta>,
    pub quarantine_entries: Vec<QuarantineEntry>,
}

/// Computes a 256-bin normalized 4-mer frequency vector for fast alignment-free homology checks.
pub fn compute_4mer_profile(seq: &str) -> [f32; 256] {
    let mut vec = [0.0f32; 256];
    let bytes = seq.as_bytes();
    if bytes.len() < 4 {
        return vec;
    }
    let mut valid = 0usize;
    for win in bytes.windows(4) {
        let mut code = 0usize;
        let mut ok = true;
        for &b in win {
            let bits = match b.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' | b'U' => 3,
                _ => {
                    ok = false;
                    break;
                }
            };
            code = (code << 2) | bits;
        }
        if ok {
            vec[code] += 1.0;
            valid += 1;
        }
    }
    if valid > 0 {
        let mut norm_sq = 0.0f32;
        for v in &vec {
            norm_sq += v * v;
        }
        let inv_norm = 1.0 / norm_sq.sqrt().max(1e-8);
        for v in &mut vec {
            *v *= inv_norm;
        }
    }
    vec
}

#[inline]
fn norm_nuc(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'U' => b'T',
        other => other,
    }
}

#[inline]
fn k2p_from_counts(compared: usize, transitions: usize, transversions: usize) -> f64 {
    if compared < 20 {
        return 0.50;
    }
    let p = transitions as f64 / (compared as f64);
    let q = transversions as f64 / (compared as f64);
    let term1 = (1.0 - 2.0 * p - q).max(0.01);
    let term2 = (1.0 - 2.0 * q).max(0.01);
    (-0.5 * term1.ln() - 0.25 * term2.ln()).clamp(0.0, 1.5)
}

fn unaligned_k2p(ba: &[u8], bb: &[u8]) -> f64 {
    let overlap = ba.len().min(bb.len());
    if overlap < 30 {
        return 0.50;
    }
    let mut transitions = 0usize;
    let mut transversions = 0usize;
    let mut compared = 0usize;
    for i in 0..overlap {
        let a = norm_nuc(ba[i]);
        let b = norm_nuc(bb[i]);
        if !matches!(a, b'A' | b'C' | b'G' | b'T') || !matches!(b, b'A' | b'C' | b'G' | b'T') {
            continue;
        }
        compared += 1;
        if a != b {
            if matches!((a, b), (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')) {
                transitions += 1;
            } else {
                transversions += 1;
            }
        }
    }
    k2p_from_counts(compared, transitions, transversions)
}

/// Computes Kimura 2-Parameter (K2P, Kimura 1980) genetic distance over homologous aligned
/// positions using a prebuilt k-mer `AnchorIndex` (or building one when `None`) plus banded
/// semi-global Needleman-Wunsch alignment across both strands, eliminating false K2P saturation
/// from 5'/3' start offsets, indels, or reverse-complemented submissions.
pub fn compute_k2p_with_index(
    index: Option<&crate::region::AnchorIndex>,
    seq_a: &str,
    seq_b: &str,
) -> f64 {
    let ba = seq_a.as_bytes();
    let bb = seq_b.as_bytes();
    if ba.len().min(bb.len()) < 30 {
        return 0.50;
    }

    let owned_idx;
    let idx = match index {
        Some(i) => i,
        None => {
            owned_idx = crate::region::AnchorIndex::new(ba);
            &owned_idx
        }
    };

    if let Some((is_rev, offset, hit_center)) = idx.best_diagonal(bb) {
        let rc_buf;
        let sb: &[u8] = if is_rev {
            rc_buf = crate::region::revcomp(bb);
            &rc_buf
        } else {
            bb
        };

        let ov_a0 = offset.max(0) as usize;
        let ov_a1 = ((offset + sb.len() as i64).max(0) as usize).min(ba.len());
        let ov_len = ov_a1.saturating_sub(ov_a0);
        if ov_len >= 30 {
            const MAX_K2P_WIN: usize = 2400;
            const BAND: i64 = 36;
            let win = ov_len.min(MAX_K2P_WIN);
            let a0 = hit_center
                .saturating_sub(win / 2)
                .max(ov_a0)
                .min(ov_a1 - win);
            let a1 = a0 + win;
            let b0_raw = a0 as i64 - offset;
            let b0 = (b0_raw - BAND).max(0) as usize;
            let b1 = (b0_raw + win as i64 + BAND).clamp(b0 as i64, sb.len() as i64) as usize;
            let sa = &ba[a0..a1];
            let sbs = &sb[b0..b1];
            let n = sa.len();
            let m = sbs.len();

            if m >= 20 {
                let d0 = b0_raw - b0 as i64;
                let width = (2 * BAND + 1) as usize;
                let mut bt = vec![0u8; (n + 1) * width];
                let mut row_lo = vec![0usize; n + 1];
                let mut row_hi = vec![0usize; n + 1];
                let neg = i32::MIN / 4;
                let mut dp_prev = vec![0i32; m + 1];
                let mut dp_curr = vec![neg; m + 1];

                row_lo[0] = 0;
                row_hi[0] = m;

                let mut best_end_score = neg;
                let mut best_end_i = n;
                let mut best_end_j = 0usize;

                for i in 1..=n {
                    let lo = (i as i64 + d0 - BAND).clamp(0, m as i64) as usize;
                    let hi = (i as i64 + d0 + BAND).clamp(lo as i64, m as i64) as usize;
                    row_lo[i] = lo;
                    row_hi[i] = hi;
                    dp_curr[lo..=hi].fill(neg);
                    let ua = norm_nuc(sa[i - 1]);
                    let valid_a = matches!(ua, b'A' | b'C' | b'G' | b'T');
                    for j in lo..=hi {
                        if j == 0 {
                            dp_curr[0] = 0;
                            bt[i * width] = 1;
                            continue;
                        }
                        let ub = norm_nuc(sbs[j - 1]);
                        let valid_b = matches!(ub, b'A' | b'C' | b'G' | b'T');
                        let sub = if valid_a && valid_b {
                            if ua == ub { 2 } else { -3 }
                        } else {
                            0
                        };
                        let mut best_v = neg;
                        let mut dir = 0u8;
                        if (j - 1) >= row_lo[i - 1]
                            && (j - 1) <= row_hi[i - 1]
                            && dp_prev[j - 1] > neg / 2
                        {
                            best_v = dp_prev[j - 1] + sub;
                            dir = 0;
                        }
                        if j >= row_lo[i - 1] && j <= row_hi[i - 1] && dp_prev[j] > neg / 2 {
                            let v = dp_prev[j] - 4;
                            if v > best_v {
                                best_v = v;
                                dir = 1;
                            }
                        }
                        if j > lo && dp_curr[j - 1] > neg / 2 {
                            let v = dp_curr[j - 1] - 4;
                            if v > best_v {
                                best_v = v;
                                dir = 2;
                            }
                        }
                        dp_curr[j] = best_v;
                        bt[i * width + (j - lo)] = dir;
                        if (i == n || j == m) && best_v > best_end_score {
                            best_end_score = best_v;
                            best_end_i = i;
                            best_end_j = j;
                        }
                    }
                    std::mem::swap(&mut dp_prev, &mut dp_curr);
                }

                if best_end_score > neg / 2 {
                    let mut ci = best_end_i;
                    let mut cj = best_end_j;
                    let mut transitions = 0usize;
                    let mut transversions = 0usize;
                    let mut compared = 0usize;
                    while ci > 0 && cj > 0 {
                        let lo = row_lo[ci];
                        let dir = bt[ci * width + (cj - lo)];
                        match dir {
                            0 => {
                                let a = norm_nuc(sa[ci - 1]);
                                let b = norm_nuc(sbs[cj - 1]);
                                if matches!(a, b'A' | b'C' | b'G' | b'T')
                                    && matches!(b, b'A' | b'C' | b'G' | b'T')
                                {
                                    compared += 1;
                                    if a != b {
                                        if matches!(
                                            (a, b),
                                            (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')
                                        ) {
                                            transitions += 1;
                                        } else {
                                            transversions += 1;
                                        }
                                    }
                                }
                                ci -= 1;
                                cj -= 1;
                            }
                            1 => ci -= 1,
                            _ => cj -= 1,
                        }
                    }
                    if compared >= 20 {
                        return k2p_from_counts(compared, transitions, transversions);
                    }
                }
            }
        }
    }

    unaligned_k2p(ba, bb)
}

/// Computes Kimura 2-Parameter (K2P, Kimura 1980) genetic distance on aligned homologous positions.
pub fn compute_approx_k2p(seq_a: &str, seq_b: &str) -> f64 {
    compute_k2p_with_index(None, seq_a, seq_b)
}

/// Executes the full 3-stage sequence analysis neural network schedule until consensus
/// converges on an optimal dataset or `max_epochs` is reached.
pub fn run_neural_network_schedule(
    cfg: &ConsensusConfig,
    seed_records: Option<Vec<SequenceRecord>>,
    telemetry_tx: Option<Sender<NeuronTelemetrySnapshot>>,
    controls: Arc<EngineControlFlags>,
) -> Result<SeqNetRunResult> {
    let started_at = Utc::now().to_rfc3339();
    // Millisecond timestamp + per-process sequence: back-to-back runs (e.g. a cold pass followed
    // by a warm delta pass) must not collide on the `neural_runs.run_id` primary key.
    static RUN_SEQ: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let run_id = format!(
        "RUN-{}-{}-{:03}",
        cfg.taxon_group.to_ascii_uppercase().replace(' ', "_"),
        Utc::now().format("%Y%m%dT%H%M%S%.3fZ"),
        RUN_SEQ.fetch_add(1, Ordering::Relaxed)
    );

    let panel = cfg.locus_panel();
    if panel.is_empty() {
        anyhow::bail!("empty locus panel (check `locus_names`: {:?})", cfg.locus_names);
    }
    let locus_names: Vec<String> = panel.iter().map(|l| l.canonical_name.clone()).collect();
    let locus_map: HashMap<String, &LocusSpec> = panel
        .iter()
        .map(|l| (l.canonical_name.clone(), l))
        .collect();

    let mut db = SeqNetDatabase::open_or_init(
        &cfg.sqlite_path,
        &cfg.database_id,
        &format!("{} Multi-Locus Cache", cfg.taxon_group),
    )?;

    let mut feedback_log: Vec<String> = Vec::new();
    feedback_log.push(format!(
        "[Init] Bound SQLite DB `{}` ({:?}) for {} ({} nuclear + {} mito loci, target >= {:.0}% per-OTU locus fill)",
        db.database_id,
        cfg.sqlite_path,
        cfg.taxon_group,
        cfg.target_nuclear_loci,
        cfg.target_mito_loci,
        cfg.min_locus_coverage_ratio * 100.0
    ));

    // Load historical soft-quarantine ledger & warm-start neural consensus weights if cached!
    let quarantine_map = db.load_quarantine_ledger(&cfg.taxon_group)?;
    let (mut consensus_net, prior_epochs) =
        if let Some((w_json, ep)) = db.load_neural_weights(&cfg.taxon_group)? {
            let parsed = serde_json::from_str::<HybridConsensusWeights>(&w_json)
                .unwrap_or_default();
            feedback_log.push(format!(
                "[Stage 2 Warm-Start] Loaded pre-trained neural consensus weights ({} prior epochs) from `{}`",
                ep, db.database_id
            ));
            (parsed, ep)
        } else {
            (HybridConsensusWeights::default(), 0)
        };
    let _ = prior_epochs;

    // Optional seed ingestion (for offline benchmarks / local imports)
    if let Some(ref seeds) = seed_records {
        db.upsert_sequences(seeds)?;
    }

    // Stage 1 Initial Input Cluster: Warm-Start Cache + Incremental NCBI Entrez `mdat` Delta Sync
    let cached_initial = db.load_cached_sequences(&locus_names)?;
    let cache_hits = cached_initial.len();
    let mut live_delta_fetched = 0usize;
    let today_mdat = Utc::now().format("%Y/%m/%d").to_string();
    let mut attempted_refills: HashSet<(String, String)> = HashSet::new();
    let min_required_loci =
        ((panel.len() as f64 * cfg.min_locus_coverage_ratio) - 1e-9).ceil().max(1.0) as usize;

    if cfg.allow_live_entrez {
        let client = EntrezClient::new(cfg);
        let mut newly_fetched = Vec::new();
        let strata = cfg.sampling_strata();
        let mut offline = false;
        let mut consecutive_failures = 0usize;
        // Mitogenome mode: mito loci are sliced from mitogenome records (one query per stratum).
        let mito_panel: Vec<LocusSpec> = panel
            .iter()
            .filter(|l| l.locus_class == LocusClass::Mitochondrial)
            .cloned()
            .collect();
        let use_mitogenomes = cfg.mitogenome_loci && !mito_panel.is_empty();
        'strata: for stratum in &strata {
            if use_mitogenomes {
                let watermark = db.get_sync_watermark(stratum, MITOGENOME_WATERMARK)?;
                let retmax = if watermark.is_some() {
                    cfg.mitogenome_retmax_per_stratum.min(150)
                } else {
                    cfg.mitogenome_retmax_per_stratum
                };
                match client.fetch_mitogenome_delta(stratum, &mito_panel, watermark.as_deref(), retmax) {
                    Ok(batch) => {
                        consecutive_failures = 0;
                        live_delta_fetched += batch.len();
                        newly_fetched.extend(batch);
                        let _ = db.upsert_sync_watermark(
                            stratum,
                            MITOGENOME_WATERMARK,
                            &today_mdat,
                            cache_hits + live_delta_fetched,
                            &run_id,
                        );
                    }
                    Err(err) => {
                        feedback_log.push(format!(
                            "[Stage 1 Entrez Notice] `{}` / mitogenomes live query failed after retries ({}); skipped",
                            stratum, err
                        ));
                        consecutive_failures += 1;
                        if consecutive_failures >= 4 {
                            offline = true;
                            break 'strata;
                        }
                    }
                }
            }
            for locus in &panel {
                if use_mitogenomes && locus.locus_class == LocusClass::Mitochondrial {
                    continue;
                }
                // Watermarks are kept per (stratum, locus) so each query only pulls its own delta.
                let watermark = db.get_sync_watermark(stratum, &locus.canonical_name)?;
                let retmax = if watermark.is_some() {
                    cfg.retmax_per_stratum.min(150)
                } else {
                    cfg.retmax_per_stratum
                };
                match client.fetch_locus_delta(stratum, locus, watermark.as_deref(), retmax) {
                    Ok(batch) => {
                        consecutive_failures = 0;
                        live_delta_fetched += batch.len();
                        newly_fetched.extend(batch);
                        let _ = db.upsert_sync_watermark(
                            stratum,
                            &locus.canonical_name,
                            &today_mdat,
                            cache_hits + live_delta_fetched,
                            &run_id,
                        );
                    }
                    Err(err) => {
                        // Skip only this (stratum, locus) after retries; give up on live fetching
                        // altogether only when Entrez keeps failing.
                        feedback_log.push(format!(
                            "[Stage 1 Entrez Notice] `{}` / `{}` live query failed after retries ({}); skipped",
                            stratum, locus.canonical_name, err
                        ));
                        consecutive_failures += 1;
                        if consecutive_failures >= 4 {
                            offline = true;
                            break 'strata;
                        }
                        continue;
                    }
                }
            }
            if strata.len() > 1 {
                println!(
                    "  [Stage 1] {:<24} cumulative fetched: {}",
                    stratum, live_delta_fetched
                );
            }
        }
        // Mitogenome mode: query nuclear loci by name for top-coverage mitogenome species per stratum,
        // considering both warm-start cache and fresh fetches, capped to avoid unbounded queries.
        if use_mitogenomes && !offline {
            let nuclear: Vec<&LocusSpec> = panel
                .iter()
                .filter(|l| l.locus_class == LocusClass::Nuclear)
                .collect();
            let mut have: HashSet<(String, String)> = HashSet::new();
            let mut mito_loci_by_sp: HashMap<String, HashSet<String>> = HashMap::new();
            let mut sp_stratum: HashMap<String, String> = HashMap::new();
            for r in cached_initial.iter().chain(newly_fetched.iter()) {
                have.insert((r.species_binomial.clone(), r.locus_canonical.clone()));
                if r.locus_class == LocusClass::Mitochondrial {
                    mito_loci_by_sp
                        .entry(r.species_binomial.clone())
                        .or_default()
                        .insert(r.locus_canonical.clone());
                }
                if !r.family.is_empty() {
                    let matches_strata = strata.iter().any(|s| s.eq_ignore_ascii_case(&r.family));
                    let entry = sp_stratum.entry(r.species_binomial.clone()).or_default();
                    let cur_matches = !entry.is_empty()
                        && strata.iter().any(|s| s.eq_ignore_ascii_case(entry));
                    if entry.is_empty()
                        || (!cur_matches && (matches_strata || r.locus_class == LocusClass::Mitochondrial))
                        || (matches_strata && r.locus_class == LocusClass::Mitochondrial)
                    {
                        *entry = r.family.clone();
                    }
                }
            }
            for rec in &mut newly_fetched {
                if let Some(parent_stratum) = sp_stratum.get(&rec.species_binomial) {
                    if !parent_stratum.is_empty() {
                        let fam_ok = !rec.family.is_empty()
                            && (strata.is_empty() || strata.iter().any(|s| s.eq_ignore_ascii_case(&rec.family)));
                        if !fam_ok {
                            rec.family = parent_stratum.clone();
                        }
                    }
                }
            }
            let mut by_stratum: std::collections::BTreeMap<String, Vec<(usize, String)>> =
                std::collections::BTreeMap::new();
            for (sp, mito_loci) in &mito_loci_by_sp {
                let mito_count = mito_loci.len();
                if mito_count + nuclear.len() < min_required_loci {
                    continue;
                }
                let st = sp_stratum.get(sp).cloned().unwrap_or_default();
                by_stratum.entry(st).or_default().push((mito_count, sp.clone()));
            }
            let mut mito_species = Vec::new();
            for (_st, mut candidates) in by_stratum {
                candidates.sort_by(|a, b| b.0.cmp(&a.0).then_with(|| a.1.cmp(&b.1)));
                candidates.truncate(24);
                for (_cnt, sp) in candidates {
                    mito_species.push(sp);
                }
            }
            let mut anchored = Vec::new();
            for sp in &mito_species {
                let parent_stratum = sp_stratum.get(sp).map(String::as_str).unwrap_or("");
                for locus in &nuclear {
                    if have.contains(&(sp.clone(), locus.canonical_name.clone())) {
                        continue;
                    }
                    if !attempted_refills.insert((sp.clone(), locus.canonical_name.clone())) {
                        continue;
                    }
                    if let Ok(mut found) = client.fetch_locus_delta(sp, locus, None, 5) {
                        if !parent_stratum.is_empty() {
                            for rec in &mut found {
                                let fam_ok = !rec.family.is_empty()
                                    && (strata.is_empty() || strata.iter().any(|s| s.eq_ignore_ascii_case(&rec.family)));
                                if !fam_ok {
                                    rec.family = parent_stratum.to_string();
                                }
                            }
                        }
                        anchored.extend(found);
                    }
                }
            }
            println!(
                "  [Stage 1] mitogenome-anchored nuclear queries: {} species x {} loci -> {} records",
                mito_species.len(),
                nuclear.len(),
                anchored.len()
            );
            live_delta_fetched += anchored.len();
            newly_fetched.extend(anchored);
        }
        if offline {
            println!("  [Stage 1] Entrez unavailable; continuing on SQLite cache");
        }
        if !newly_fetched.is_empty() {
            db.upsert_sequences(&newly_fetched)?;
            feedback_log.push(format!(
                "[Stage 1 Delta Sync] Ingested {} new/modified NCBI accessions into `{}`",
                newly_fetched.len(),
                db.database_id
            ));
        }
    }

    let mut working_pool = db.load_cached_sequences(&locus_names)?;
    let mut pending_directives: Vec<GapRefillDirective> = Vec::new();

    let mut best_otus: Vec<AssembledSpeciesOtu> = Vec::new();
    let mut latest_quarantine_vec: Vec<QuarantineEntry> = Vec::new();
    let mut final_consensus_score = 0.0;
    let mut final_occupancy = 0.0;
    let mut converged = false;
    let mut epochs_executed = 0usize;
    let mut total_rehabilitated = 0usize;
    let mut final_training_samples: Vec<(SeqFeatureVec, f64, f64)> = Vec::new();

    // ========================================================================
    // 3-STAGE CYCLIC NEURAL PROCESSING & FEEDBACK SCHEDULE
    // ========================================================================
    for epoch in 1..=cfg.max_epochs {
        while controls.paused.load(Ordering::Relaxed) {
            std::thread::sleep(std::time::Duration::from_millis(80));
        }
        epochs_executed = epoch;

        // Stage 1 Feedback Refill: If Epoch > 1 and Stage 3 requested targeted species/locus fills
        if epoch > 1 && cfg.allow_live_entrez && !pending_directives.is_empty() {
            let client = EntrezClient::new(cfg);
            let strata = cfg.sampling_strata();
            let mut sp_stratum_pool: HashMap<String, String> = HashMap::new();
            let mut pool_have: HashSet<(String, String)> = HashSet::new();
            for r in &working_pool {
                pool_have.insert((r.species_binomial.clone(), r.locus_canonical.clone()));
                if !r.family.is_empty() {
                    let matches_strata = strata.iter().any(|s| s.eq_ignore_ascii_case(&r.family));
                    let entry = sp_stratum_pool.entry(r.species_binomial.clone()).or_default();
                    let cur_matches = !entry.is_empty()
                        && strata.iter().any(|s| s.eq_ignore_ascii_case(entry));
                    if entry.is_empty()
                        || (!cur_matches && (matches_strata || r.locus_class == LocusClass::Mitochondrial))
                        || (matches_strata && r.locus_class == LocusClass::Mitochondrial)
                    {
                        *entry = r.family.clone();
                    }
                }
            }
            let mut refill_records = Vec::new();
            for dir in pending_directives.iter().take(cfg.max_refill_directives) {
                // In mitogenome mode a species' mito loci come from its mitogenome; per-locus
                // refills would mix specimens, so only nuclear gaps are refilled.
                let refillable = dir.missing_loci.iter().filter(|l| {
                    let key = (dir.species_binomial.clone(), (*l).clone());
                    !pool_have.contains(&key)
                        && !attempted_refills.contains(&key)
                        && !(cfg.mitogenome_loci
                            && locus_map
                                .get(*l)
                                .map(|s| s.locus_class == LocusClass::Mitochondrial)
                                .unwrap_or(false))
                });
                let parent_stratum = sp_stratum_pool
                    .get(&dir.species_binomial)
                    .map(String::as_str)
                    .unwrap_or("");
                let chosen: Vec<String> = refillable
                    .take(cfg.max_refill_loci_per_directive)
                    .cloned()
                    .collect();
                for missing_loc in chosen {
                    attempted_refills.insert((dir.species_binomial.clone(), missing_loc.clone()));
                    if let Some(spec) = locus_map.get(&missing_loc) {
                        if let Ok(mut found) = client.fetch_locus_delta(
                            &dir.species_binomial,
                            spec,
                            None,
                            25,
                        ) {
                            if !parent_stratum.is_empty() {
                                for rec in &mut found {
                                    let fam_ok = !rec.family.is_empty()
                                        && (strata.is_empty()
                                            || strata.iter().any(|s| s.eq_ignore_ascii_case(&rec.family)));
                                    if !fam_ok {
                                        rec.family = parent_stratum.to_string();
                                    }
                                }
                            }
                            refill_records.extend(found);
                        }
                    }
                }
            }
            if !refill_records.is_empty() {
                live_delta_fetched += refill_records.len();
                db.upsert_sequences(&refill_records)?;
                working_pool = db.load_cached_sequences(&locus_names)?;
                feedback_log.push(format!(
                    "[Epoch {} Stage 1 Refill] Pulled {} targeted gap-fill sequences from NCBI",
                    epoch,
                    refill_records.len()
                ));
            }
        }

        // Precompute cross-locus voucher linkage counts across the pool
        let mut voucher_locus_counts: HashMap<&str, HashSet<&str>> = HashMap::new();
        for rec in &working_pool {
            if !rec.voucher_id.is_empty() {
                voucher_locus_counts
                    .entry(rec.voucher_id.as_str())
                    .or_default()
                    .insert(rec.locus_canonical.as_str());
            }
        }

        // Group pool by locus for parallel Stage 2 Within-Locus Neuron evaluation via `rayon`
        let mut pool_by_locus: HashMap<&str, Vec<&SequenceRecord>> = HashMap::new();
        for rec in &working_pool {
            pool_by_locus
                .entry(rec.locus_canonical.as_str())
                .or_default()
                .push(rec);
        }

        // Compute locus-level 4-mer centroids, median GC, and within-locus K2P cohesion in parallel
        let locus_evals: Vec<(
            String,
            Vec<(SequenceRecord, SeqFeatureVec, f64, f64, Option<&'static str>)>,
        )> = panel
            .par_iter()
            .map(|locus| {
                let Some(recs) = pool_by_locus.get(locus.canonical_name.as_str()) else {
                    return (locus.canonical_name.clone(), Vec::new());
                };
                if recs.is_empty() {
                    return (locus.canonical_name.clone(), Vec::new());
                }

                // Compute 4-mer profiles & locus centroid
                let profiles: Vec<[f32; 256]> = recs
                    .iter()
                    .map(|r| compute_4mer_profile(&r.sequence_raw))
                    .collect();
                let mut centroid = [0.0f32; 256];
                let mut valid_count = 0usize;
                for (r, p) in recs.iter().zip(profiles.iter()) {
                    if r.reading_frame_ok && r.ambig_fraction <= cfg.max_ambig_fraction {
                        for b in 0..256 {
                            centroid[b] += p[b];
                        }
                        valid_count += 1;
                    }
                }
                if valid_count == 0 {
                    for p in &profiles {
                        for b in 0..256 {
                            centroid[b] += p[b];
                        }
                    }
                }
                let c_norm: f32 = centroid.iter().map(|v| v * v).sum::<f32>().sqrt().max(1e-8);
                for v in &mut centroid {
                    *v /= c_norm;
                }

                // Median GC & representative consensus anchor sequence
                let mut gc_vals: Vec<f64> = recs.iter().map(|r| r.gc_fraction).collect();
                gc_vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let gc_med = gc_vals[gc_vals.len() / 2];

                // Pick sequence with highest k-mer similarity to centroid as locus reference anchor
                let sims: Vec<f64> = profiles
                    .iter()
                    .map(|p| {
                        p.iter()
                            .zip(centroid.iter())
                            .map(|(a, b)| (*a as f64) * (*b as f64))
                            .sum::<f64>()
                            .clamp(0.0, 1.0)
                    })
                    .collect();

                let anchor_idx = sims
                    .iter()
                    .enumerate()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
                    .map(|(i, _)| i)
                    .unwrap_or(0);
                let anchor_seq = &recs[anchor_idx].sequence_raw;
                let k2p_index = crate::region::AnchorIndex::new(anchor_seq.as_bytes());

                let k2p_dists: Vec<f64> = recs
                    .iter()
                    .map(|r| compute_k2p_with_index(Some(&k2p_index), anchor_seq, &r.sequence_raw))
                    .collect();

                let mut sorted_k2p = k2p_dists.clone();
                sorted_k2p.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let k2p_med = sorted_k2p[sorted_k2p.len() / 2];
                let mut abs_devs: Vec<f64> =
                    sorted_k2p.iter().map(|d| (d - k2p_med).abs()).collect();
                abs_devs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                // Robust MAD floor (`0.035` prevents zero-variance division on conserved nuclear genes)
                let k2p_mad = (abs_devs[abs_devs.len() / 2] * 1.4826).max(0.035);

                // Region coherence: map every record onto a centroid-representative long anchor and
                // keep only records covering the locus core (stretch shared by most records), so one
                // locus column never mixes non-overlapping fragments (e.g. D1-D2 vs D2-D3 28S).
                let region_cov: Vec<Option<f64>> = if cfg.region_min_core_coverage > 0.0 {
                    let mut clean_lens: Vec<usize> = recs
                        .iter()
                        .filter(|r| r.reading_frame_ok && r.ambig_fraction <= cfg.max_ambig_fraction)
                        .map(|r| r.seq_length)
                        .collect();
                    if clean_lens.is_empty() {
                        clean_lens = recs.iter().map(|r| r.seq_length).collect();
                    }
                    clean_lens.sort_unstable();
                    let p90_len = clean_lens
                        .get((clean_lens.len() * 9 / 10).min(clean_lens.len().saturating_sub(1)))
                        .copied()
                        .unwrap_or(0);
                    let min_anchor_len = (p90_len * 3) / 4;
                    let region_anchor = recs
                        .iter()
                        .enumerate()
                        .filter(|(_, r)| {
                            r.reading_frame_ok
                                && r.ambig_fraction <= cfg.max_ambig_fraction
                                && r.seq_length >= min_anchor_len
                        })
                        .max_by(|(ia, _), (ib, _)| {
                            sims[*ia]
                                .partial_cmp(&sims[*ib])
                                .unwrap_or(std::cmp::Ordering::Equal)
                        })
                        .or_else(|| {
                            recs.iter()
                                .enumerate()
                                .filter(|(_, r)| {
                                    r.reading_frame_ok && r.ambig_fraction <= cfg.max_ambig_fraction
                                })
                                .max_by_key(|(_, r)| r.seq_length)
                        })
                        .or_else(|| recs.iter().enumerate().max_by_key(|(_, r)| r.seq_length))
                        .map(|(_, r)| r.sequence_raw.as_bytes())
                        .unwrap_or(&[]);
                    let index = crate::region::AnchorIndex::new(region_anchor);
                    let spans: Vec<Option<crate::region::AnchorSpan>> =
                        recs.iter().map(|r| index.map(r.sequence_raw.as_bytes())).collect();
                    let mapped: Vec<crate::region::AnchorSpan> = spans.iter().flatten().copied().collect();
                    match crate::region::core_window(&mapped, index.len, cfg.region_core_support) {
                        // On a real (non-repetitive) anchor, a record with no k-mer diagonal at all
                        // is not homologous to the locus (e.g. a predicted mRNA or ITS-only record
                        // matched by a keyword) and gets zero core coverage. Repetitive anchors
                        // (complexity < 0.8) cannot place fragments, so the gate stays off.
                        Some(core) if mapped.len() >= 6 && index.complexity() >= 0.8 => spans
                            .iter()
                            .map(|s| Some(s.map(|s| crate::region::core_coverage(&s, core)).unwrap_or(0.0)))
                            .collect(),
                        _ => vec![None; recs.len()],
                    }
                } else {
                    vec![None; recs.len()]
                };

                let mut per_seq = Vec::with_capacity(recs.len());
                for idx in 0..recs.len() {
                    let rec = recs[idx];
                    let kmer_sim = sims[idx];
                    let k2p = k2p_dists[idx];
                    let k2p_z = ((k2p - k2p_med) / k2p_mad).max(0.0);
                    let gc_z = (rec.gc_fraction - gc_med) / 0.06;

                    // Deterministic Stage 2 Quarantine Gate Checks:
                    let hard_fail: Option<&'static str> = if rec.seq_length < locus.min_expected_bp
                        || rec.seq_length > locus.max_expected_bp * 2
                    {
                        Some("SEQUENCE_LENGTH_OUT_OF_BOUNDS")
                    } else if rec.ambig_fraction > cfg.max_ambig_fraction {
                        Some("EXCESS_AMBIGUOUS_BASES")
                    } else if locus.protein_coding && !rec.reading_frame_ok {
                        Some("INTERNAL_STOP_CODON_NUMT")
                    } else if region_cov[idx].map(|c| c < cfg.region_min_core_coverage).unwrap_or(false) {
                        Some("REGION_OUTSIDE_LOCUS_CORE")
                    } else if k2p_z > cfg.max_within_locus_zscore && kmer_sim < 0.72 {
                        Some("WITHIN_LOCUS_K2P_SAMPLE_SWAP")
                    } else {
                        None
                    };

                    let len_span = (locus.max_expected_bp.saturating_sub(locus.min_expected_bp))
                        .max(100) as f64;
                    let len_adequacy = ((rec.seq_length.saturating_sub(locus.min_expected_bp))
                        as f64
                        / len_span)
                        .clamp(0.25, 1.0);
                    let ambig_clean = (1.0 - rec.ambig_fraction / 0.10).clamp(0.0, 1.0);
                    let rf_score = if rec.reading_frame_ok { 1.0 } else { 0.0 };
                    let gc_homog = (-0.5 * gc_z * gc_z).exp();
                    let k2p_cohesion = (-0.35 * k2p_z).exp();
                    let v_span = voucher_locus_counts
                        .get(rec.voucher_id.as_str())
                        .map(|s| s.len())
                        .unwrap_or(0) as f64
                        / (panel.len().max(1) as f64);
                    let entropy_norm = (rec.shannon_entropy / 2.0).clamp(0.0, 1.0);

                    let feats: SeqFeatureVec = [
                        len_adequacy,
                        ambig_clean,
                        rf_score,
                        gc_homog,
                        kmer_sim,
                        k2p_cohesion,
                        v_span,
                        entropy_norm,
                    ];

                    per_seq.push(((*rec).clone(), feats, k2p, kmer_sim, hard_fail));
                }

                (locus.canonical_name.clone(), per_seq)
            })
            .collect();

        // Stage 2 Hybrid Neural Consensus Scoring + Soft-Quarantine Re-Evaluation + Online Training Batch
        let mut clean_candidates: Vec<ValidatedCandidate> = Vec::new();
        let mut training_samples: Vec<(SeqFeatureVec, f64, f64)> = Vec::new();
        let mut epoch_quarantine_entries: Vec<QuarantineEntry> = Vec::new();
        let mut epoch_rehabilitated = 0usize;

        for (_, eval_rows) in locus_evals {
            for (rec, feats, k2p, kmer_sim, hard_fail) in eval_rows {
                let key = (rec.accession_version.clone(), rec.locus_canonical.clone());
                let prev_q = quarantine_map.get(&key);
                let prior_q_count = prev_q.map(|q| q.prior_quarantine_count).unwrap_or(0);
                let prior_pass_count = prev_q.map(|q| q.prior_pass_count).unwrap_or(0);
                let was_previously_quarantined = prev_q
                    .map(|q| q.currently_quarantined || q.prior_quarantine_count > 0)
                    .unwrap_or(false);

                // Soft-quarantine skepticism prior: decays gently if rehabilitated repeatedly,
                // never permanently locks out a sequence if new cohort evidence vindicates it.
                let net_strikes = (prior_q_count as f64 - 0.5 * prior_pass_count as f64).max(0.0);
                let skep_prior = 1.0 - (-0.38 * net_strikes).exp();

                let (neural_prob, _) = consensus_net.forward(&feats, skep_prior);

                let quarantine_reason: Option<&'static str> = if let Some(r) = hard_fail {
                    Some(r)
                } else if neural_prob < cfg.neural_pass_threshold {
                    Some("LOW_HYBRID_NEURAL_CONSENSUS")
                } else {
                    None
                };

                if let Some(reason) = quarantine_reason {
                    // Add negative contrastive training signal on deterministic gate failures
                    if hard_fail.is_some() {
                        training_samples.push((feats, skep_prior, 0.0));
                    }
                    epoch_quarantine_entries.push(QuarantineEntry {
                        accession_version: rec.accession_version.clone(),
                        locus_canonical: rec.locus_canonical.clone(),
                        taxon_group: cfg.taxon_group.clone(),
                        genus: rec.genus.clone(),
                        species_binomial: rec.species_binomial.clone(),
                        currently_quarantined: true,
                        prior_quarantine_count: prior_q_count + 1,
                        prior_pass_count,
                        skepticism_prior: (1.0 - (-0.38 * (net_strikes + 1.0)).exp()),
                        last_reason_code: reason.to_string(),
                        within_locus_k2p_dist: k2p,
                        neural_consensus_score: neural_prob,
                        last_evaluated_run_id: run_id.clone(),
                        audit_notes: format!(
                            "Epoch {} soft-quarantined (kmer_sim={:.3}, k2p={:.4})",
                            epoch, kmer_sim, k2p
                        ),
                    });
                } else {
                    // Passed both deterministic bioinformatic gates and hybrid neural threshold!
                    if kmer_sim > 0.90 && k2p < 0.15 {
                        training_samples.push((feats, skep_prior, 1.0));
                    }
                    let rehab = was_previously_quarantined;
                    if rehab {
                        epoch_rehabilitated += 1;
                        epoch_quarantine_entries.push(QuarantineEntry {
                            accession_version: rec.accession_version.clone(),
                            locus_canonical: rec.locus_canonical.clone(),
                            taxon_group: cfg.taxon_group.clone(),
                            genus: rec.genus.clone(),
                            species_binomial: rec.species_binomial.clone(),
                            currently_quarantined: false,
                            prior_quarantine_count: prior_q_count,
                            prior_pass_count: prior_pass_count + 1,
                            skepticism_prior: (1.0 - (-0.38 * (net_strikes - 0.5).max(0.0)).exp()),
                            last_reason_code: "REHABILITATED_BY_COHORT_CONSENSUS".to_string(),
                            within_locus_k2p_dist: k2p,
                            neural_consensus_score: neural_prob,
                            last_evaluated_run_id: run_id.clone(),
                            audit_notes: format!(
                                "Rehabilitated on Epoch {} (prior_q={})",
                                epoch, prior_q_count
                            ),
                        });
                    }
                    clean_candidates.push(ValidatedCandidate {
                        record: rec,
                        features: feats,
                        within_locus_k2p_dist: k2p,
                        kmer_cosine_sim: kmer_sim,
                        neural_confidence: neural_prob,
                        rehabilitated_from_quarantine: rehab,
                    });
                }
            }
        }

        // Weights stay fixed within a run so every epoch scores sequences by the same rule
        // (per-epoch self-training drifted the boundary toward rejection on real cohorts).
        // The final epoch's labels train the persisted weights once, after the loop.
        let bce_loss = consensus_net.evaluate_loss(&training_samples);
        final_training_samples = training_samples;

        // Stage 3 Output & Between-Locus Valuation Neurons
        let valuation = assemble_and_value_otus(cfg, &panel, &clean_candidates);
        let mean_otu_conf = if valuation.retained_otus.is_empty() {
            0.0
        } else {
            valuation
                .retained_otus
                .iter()
                .map(|o| o.mean_neural_confidence)
                .sum::<f64>()
                / (valuation.retained_otus.len() as f64)
        };

        // Composite Global Consensus Objective combining matrix occupancy (>=90% goal) & neural confidence
        let epoch_consensus = if valuation.retained_otus.is_empty() {
            0.0
        } else {
            0.60 * valuation.matrix_occupancy + 0.40 * mean_otu_conf
        };

        let chimeric_count = valuation
            .retained_otus
            .iter()
            .filter(|o| o.is_chimeric)
            .count();
        let active_q_count = epoch_quarantine_entries
            .iter()
            .filter(|q| q.currently_quarantined)
            .count();

        feedback_log.push(format!(
            "[Epoch {}/{}] Clean={} | Quarantined={} (Rehab={}) | RetainedOTUs={} ({} genera, {} chimeric) | Occupancy={:.2}% | Consensus={:.4} (loss={:.4})",
            epoch,
            cfg.max_epochs,
            clean_candidates.len(),
            active_q_count,
            epoch_rehabilitated,
            valuation.retained_otus.len(),
            valuation.retained_genera_count,
            chimeric_count,
            valuation.matrix_occupancy * 100.0,
            epoch_consensus,
            bce_loss
        ));

        if cfg.execution_mode == crate::config::ExecutionMode::Headless
            || cfg.execution_mode == crate::config::ExecutionMode::GuiMonitor
        {
            println!(
                "  Epoch {:02}/{} | Clean: {:4} | Quarantined: {:3} | OTUs: {:3} ({:2} genera, {:3} chimeric) | Matrix Fill: {:6.2}% | Consensus: {:.4}",
                epoch,
                cfg.max_epochs,
                clean_candidates.len(),
                active_q_count,
                valuation.retained_otus.len(),
                valuation.retained_genera_count,
                chimeric_count,
                valuation.matrix_occupancy * 100.0,
                epoch_consensus
            );
        }

        let prev_consensus = final_consensus_score;
        let prev_otu_count = best_otus.len();
        let otus_still_growing = epoch >= 2 && valuation.retained_otus.len() > prev_otu_count;
        best_otus = valuation.retained_otus;
        latest_quarantine_vec = epoch_quarantine_entries;
        final_occupancy = valuation.matrix_occupancy;
        final_consensus_score = epoch_consensus;
        total_rehabilitated = epoch_rehabilitated;
        pending_directives = valuation.gap_directives;

        // Check convergence criteria:
        // 1. Retained OTUs non-empty and matrix occupancy >= min_locus_coverage_ratio (e.g. 90%)
        // 2. Either `epoch_consensus >= optimal_consensus_target` OR consensus plateaued (`|delta| < 0.002` after >=2 epochs, or 1st epoch warm-start hit optimal!)
        let warm_start_instant = cache_hits > 0
            && prior_epochs > 0
            && final_occupancy >= cfg.min_locus_coverage_ratio
            && epoch_consensus >= cfg.optimal_consensus_target;

        // A consensus plateau only counts once gap-refill has stopped adding OTUs
        let plateaued = epoch >= 2
            && !otus_still_growing
            && final_occupancy >= cfg.min_locus_coverage_ratio
            && ((epoch_consensus - prev_consensus).abs() < 0.0025
                || epoch_consensus >= cfg.optimal_consensus_target);

        if warm_start_instant || plateaued || controls.force_finalize.load(Ordering::Relaxed) {
            converged = final_occupancy >= cfg.min_locus_coverage_ratio && !best_otus.is_empty();
            feedback_log.push(format!(
                "[Consensus Reached] Optimal dataset converged at Epoch {} (Matrix Fill {:.2}%, Consensus {:.4})",
                epoch,
                final_occupancy * 100.0,
                final_consensus_score
            ));
        }

        if let Some(ref tx) = telemetry_tx {
            let _ = tx.send(NeuronTelemetrySnapshot {
                database_id: db.database_id.clone(),
                run_id: run_id.clone(),
                taxon_group: cfg.taxon_group.clone(),
                epoch,
                max_epochs: cfg.max_epochs,
                active_stage: if converged {
                    "Stage 3: MSA & NEXUS Export".to_string()
                } else {
                    "Stage 2 <-> Stage 3 Feedback Loop".to_string()
                },
                cache_hits,
                live_delta_fetched,
                total_sequences_evaluated: working_pool.len(),
                clean_candidates_count: clean_candidates.len(),
                quarantined_count: active_q_count,
                rehabilitated_count: total_rehabilitated,
                retained_otus_count: best_otus.len(),
                retained_genera_count: valuation.retained_genera_count,
                chimeric_otus_count: chimeric_count,
                matrix_occupancy: final_occupancy,
                global_consensus_score: final_consensus_score,
                neural_bce_loss: bce_loss,
                converged,
                finished: false,
                locus_names: locus_names.clone(),
                retained_otus_preview: best_otus.clone(),
                quarantine_preview: latest_quarantine_vec.clone(),
                feedback_log: feedback_log.clone(),
            });
        }

        // Priors stay as loaded at run start: re-evaluating the same record in a later epoch of
        // this run must not count as a new strike (that compounding drove quarantine counts up
        // every epoch). The ledger update is persisted once, after the final epoch.

        if converged {
            break;
        }
    }

    // ========================================================================
    // FINAL STAGE 3 ALIGNMENT, SQLITE PERSISTENCE & CLEAN NEXUS EXPORT
    // ========================================================================
    let partitions = align_all_loci_parallel(&panel, &mut best_otus);

    db.save_quarantine_updates(&latest_quarantine_vec)?;
    // One learning step per run, on the final epoch's labels; persisted for warm starts.
    let _ = consensus_net.online_feedback_step(&final_training_samples);
    if let Ok(w_json) = serde_json::to_string(&consensus_net) {
        let _ = db.save_neural_weights(&cfg.taxon_group, &w_json, epochs_executed, 0.0);
    }

    let active_quarantined_total = latest_quarantine_vec
        .iter()
        .filter(|q| q.currently_quarantined)
        .count();

    let completed_at = Utc::now().to_rfc3339();
    let audit = RunAuditRecord {
        run_id: run_id.clone(),
        database_id: db.database_id.clone(),
        taxon_group: cfg.taxon_group.clone(),
        epochs_executed,
        cache_hits,
        live_delta_fetched,
        quarantined_count: active_quarantined_total,
        rehabilitated_count: total_rehabilitated,
        converged,
        final_consensus_score,
        final_matrix_occupancy: final_occupancy,
        nexus_output_path: cfg.nexus_output_path.to_string_lossy().to_string(),
        started_at,
        completed_at,
    };

    db.save_run_and_otus(cfg, &audit, &best_otus)?;
    write_partitioned_nexus(cfg, &audit, &best_otus, &partitions, &latest_quarantine_vec)?;

    feedback_log.push(format!(
        "[Output Complete] Exported clean NEXUS supermatrix ({} OTUs x {} loci) -> {:?}",
        best_otus.len(),
        partitions.len(),
        cfg.nexus_output_path
    ));

    if let Some(ref tx) = telemetry_tx {
        let chimeric_count = best_otus.iter().filter(|o| o.is_chimeric).count();
        let genera_set: HashSet<&str> = best_otus.iter().map(|o| o.genus.as_str()).collect();
        let _ = tx.send(NeuronTelemetrySnapshot {
            database_id: db.database_id.clone(),
            run_id,
            taxon_group: cfg.taxon_group.clone(),
            epoch: epochs_executed,
            max_epochs: cfg.max_epochs,
            active_stage: "Completed â€” NEXUS & SQLite Synced".to_string(),
            cache_hits,
            live_delta_fetched,
            total_sequences_evaluated: working_pool.len(),
            clean_candidates_count: working_pool.len().saturating_sub(active_quarantined_total),
            quarantined_count: active_quarantined_total,
            rehabilitated_count: total_rehabilitated,
            retained_otus_count: best_otus.len(),
            retained_genera_count: genera_set.len(),
            chimeric_otus_count: chimeric_count,
            matrix_occupancy: final_occupancy,
            global_consensus_score: final_consensus_score,
            neural_bce_loss: 0.0,
            converged,
            finished: true,
            locus_names,
            retained_otus_preview: best_otus.clone(),
            quarantine_preview: latest_quarantine_vec.clone(),
            feedback_log,
        });
    }

    Ok(SeqNetRunResult {
        audit,
        otus: best_otus,
        partitions,
        quarantine_entries: latest_quarantine_vec,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pseudo_dna(len: usize, mut seed: u64) -> String {
        (0..len)
            .map(|_| {
                seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ["A", "C", "G", "T"][(seed >> 33) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn homology_k2p_handles_shift_indel_and_revcomp_without_false_saturation() {
        let anchor = pseudo_dna(600, 9001);
        // +3 bp shifted homolog with 10% point mutations
        let mut shifted_bytes = anchor.as_bytes()[3..].to_vec();
        for i in (5..shifted_bytes.len()).step_by(10) {
            shifted_bytes[i] = match shifted_bytes[i] {
                b'A' => b'G',
                b'G' => b'A',
                b'C' => b'T',
                _ => b'C',
            };
        }
        let shifted = String::from_utf8(shifted_bytes).unwrap();
        let k2p_shift = compute_approx_k2p(&anchor, &shifted);
        assert!(
            k2p_shift < 0.25,
            "3-bp shifted homolog must not saturate K2P (got {k2p_shift})"
        );

        // Internal 6-bp deletion homolog
        let mut deleted = anchor[..210].to_string();
        deleted.push_str(&anchor[216..]);
        let k2p_indel = compute_approx_k2p(&anchor, &deleted);
        assert!(
            k2p_indel < 0.05,
            "6-bp internal indel homolog must stay aligned (got {k2p_indel})"
        );

        // Reverse-complement homolog
        let rc = String::from_utf8(crate::region::revcomp(shifted.as_bytes())).unwrap();
        let k2p_rc = compute_approx_k2p(&anchor, &rc);
        assert!(
            k2p_rc < 0.25,
            "reverse-complement homolog must align cleanly (got {k2p_rc})"
        );

        // Unrelated sequence must still show high divergence
        let unrelated = pseudo_dna(600, 424242);
        let k2p_unrel = compute_approx_k2p(&anchor, &unrelated);
        assert!(
            k2p_unrel >= 0.80,
            "unrelated decoy must yield high K2P (got {k2p_unrel})"
        );
    }
}
