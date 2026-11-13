//! Stage 3 Chimeric OTU Synthesis & Between-Locus Coverage Valuation (`otu.rs`).
//!
//! Constructs species-level OTU terminals organized under Genus-level constraints
//! (`>= min_species_per_genus` separate species terminals per retained Genus, never collapsing
//! an entire genus into a single hybrid row). Within each species terminal, merges multiple
//! conspecific specimens/vouchers when needed to maximize locus coverage (`>= 90%`), while
//! prioritizing shared `voucher_id` linkage to minimize unnecessary chimerism.

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap, HashSet};

use crate::config::{ConsensusConfig, LocusClass, LocusSpec, OtuRank};
use crate::neurons::ValidatedCandidate;

/// Single locus assignment inside an assembled species-level OTU terminal.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OtuLocusSlot {
    pub locus_canonical: String,
    pub locus_class: LocusClass,
    pub accession_version: String,
    pub voucher_id: String,
    pub organism_name: String,
    pub sequence_raw: String,
    pub aligned_sequence: Option<String>,
    pub neural_confidence: f64,
    pub within_locus_k2p_dist: f64,
}

/// Species-level terminal OTU (potentially chimeric across multiple conspecific specimens).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AssembledSpeciesOtu {
    /// Sanitized NEXUS taxon label, e.g. `Carabus_nemoralis`
    pub terminal_label: String,
    pub genus: String,
    pub species_binomial: String,
    pub family: String,
    /// `true` when >= 2 distinct specimens/vouchers contributed sequences across loci.
    pub is_chimeric: bool,
    pub distinct_specimen_count: usize,
    pub distinct_accession_count: usize,
    pub loci_filled: usize,
    pub loci_target: usize,
    pub locus_coverage_ratio: f64,
    pub mean_neural_confidence: f64,
    /// Map from `locus_canonical` -> chosen clean accession slot.
    pub locus_slots: BTreeMap<String, OtuLocusSlot>,
}

/// Feedback signal emitted by Stage 3 Between-Locus Valuation back to Stage 1 Input Neurons.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GapRefillDirective {
    pub genus: String,
    pub species_binomial: String,
    pub current_coverage_ratio: f64,
    pub missing_loci: Vec<String>,
    pub priority_score: f64,
}

/// Summary report from Stage 3 OTU synthesis & cross-locus valuation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OtuValuationReport {
    pub retained_otus: Vec<AssembledSpeciesOtu>,
    pub dropped_low_coverage_species: Vec<(String, f64)>,
    pub dropped_singleton_genera: Vec<(String, usize)>,
    pub gap_directives: Vec<GapRefillDirective>,
    pub matrix_occupancy: f64,
    pub retained_genera_count: usize,
}

/// Synthesizes species-level chimeric OTU terminals from non-quarantined Stage 2 candidates,
/// enforces `min_locus_coverage_ratio` (e.g. `0.90`) and `min_species_per_genus` (e.g. `2`),
/// and generates targeted gap-refill feedback directives for near-threshold taxa.
pub fn assemble_and_value_otus(
    cfg: &ConsensusConfig,
    panel: &[LocusSpec],
    clean_candidates: &[ValidatedCandidate],
) -> OtuValuationReport {
    let loci_target = panel.len().max(1);
    let panel_names: Vec<String> = panel.iter().map(|l| l.canonical_name.clone()).collect();

    // Group clean candidates by `species_binomial`
    let mut by_species: BTreeMap<String, Vec<&ValidatedCandidate>> = BTreeMap::new();
    for cand in clean_candidates {
        by_species
            .entry(cand.record.species_binomial.clone())
            .or_default()
            .push(cand);
    }

    let mut draft_species_otus: Vec<AssembledSpeciesOtu> = Vec::new();
    let mut dropped_low_coverage_species: Vec<(String, f64)> = Vec::new();
    let mut gap_directives: Vec<GapRefillDirective> = Vec::new();

    for (sp_binomial, cands) in by_species {
        let genus = cands[0].record.genus.clone();
        let family = cands[0].record.family.clone();

        // Count how many distinct loci each non-empty `voucher_id` spans for this species
        // so we prefer assembling around a primary museum voucher before grafting chimeric loci.
        let mut voucher_locus_span: HashMap<&str, HashSet<&str>> = HashMap::new();
        for c in &cands {
            if !c.record.voucher_id.is_empty() {
                voucher_locus_span
                    .entry(c.record.voucher_id.as_str())
                    .or_default()
                    .insert(c.record.locus_canonical.as_str());
            }
        }

        // Group candidates for this species by locus
        let mut by_locus: HashMap<&str, Vec<&ValidatedCandidate>> = HashMap::new();
        for c in &cands {
            by_locus
                .entry(c.record.locus_canonical.as_str())
                .or_default()
                .push(c);
        }

        let mut locus_slots = BTreeMap::new();
        let mut specimens_used: HashSet<String> = HashSet::new();
        let mut accessions_used: HashSet<String> = HashSet::new();
        let mut conf_sum = 0.0;

        for locus_name in &panel_names {
            if let Some(options) = by_locus.get(locus_name.as_str()) {
                // Select best accession maximizing:
                //   neural_confidence + 0.08 * (voucher_span / loci_target) + 0.03 * length_ratio
                let best = options
                    .iter()
                    .copied()
                    .max_by(|a, b| {
                        let span_a = voucher_locus_span
                            .get(a.record.voucher_id.as_str())
                            .map(|s| s.len())
                            .unwrap_or(0) as f64
                            / (loci_target as f64);
                        let span_b = voucher_locus_span
                            .get(b.record.voucher_id.as_str())
                            .map(|s| s.len())
                            .unwrap_or(0) as f64
                            / (loci_target as f64);
                        let score_a = a.neural_confidence + 0.08 * span_a;
                        let score_b = b.neural_confidence + 0.08 * span_b;
                        score_a
                            .partial_cmp(&score_b)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .unwrap();

                let spec_key = if !best.record.voucher_id.is_empty() {
                    best.record.voucher_id.clone()
                } else if !best.record.isolate_code.is_empty() {
                    format!("ISO:{}", best.record.isolate_code)
                } else {
                    format!("ACC:{}", best.record.accession_version)
                };
                specimens_used.insert(spec_key);
                accessions_used.insert(best.record.accession_version.clone());
                conf_sum += best.neural_confidence;

                locus_slots.insert(
                    locus_name.clone(),
                    OtuLocusSlot {
                        locus_canonical: locus_name.clone(),
                        locus_class: best.record.locus_class,
                        accession_version: best.record.accession_version.clone(),
                        voucher_id: best.record.voucher_id.clone(),
                        organism_name: best.record.organism_name.clone(),
                        sequence_raw: best.record.sequence_raw.clone(),
                        aligned_sequence: None,
                        neural_confidence: best.neural_confidence,
                        within_locus_k2p_dist: best.within_locus_k2p_dist,
                    },
                );
            }
        }

        let loci_filled = locus_slots.len();
        let coverage = loci_filled as f64 / (loci_target as f64);

        if coverage + 1e-9 >= cfg.min_locus_coverage_ratio {
            let distinct_specs = specimens_used.len();
            let distinct_accs = accessions_used.len();
            let safe_label = sp_binomial
                .chars()
                .map(|ch| if ch.is_ascii_alphanumeric() { ch } else { '_' })
                .collect::<String>();

            draft_species_otus.push(AssembledSpeciesOtu {
                terminal_label: safe_label,
                genus,
                species_binomial: sp_binomial,
                family,
                is_chimeric: distinct_specs > 1,
                distinct_specimen_count: distinct_specs,
                distinct_accession_count: distinct_accs,
                loci_filled,
                loci_target,
                locus_coverage_ratio: coverage,
                mean_neural_confidence: if loci_filled > 0 {
                    conf_sum / (loci_filled as f64)
                } else {
                    0.0
                },
                locus_slots,
            });
        } else {
            dropped_low_coverage_species.push((sp_binomial.clone(), coverage));
            // If within striking distance (`refill_min_coverage`, default 55%), request targeted Stage 1 refill
            if coverage + 1e-9 >= cfg.refill_min_coverage {
                let missing: Vec<String> = panel_names
                    .iter()
                    .filter(|ln| !locus_slots.contains_key(*ln))
                    .cloned()
                    .collect();
                gap_directives.push(GapRefillDirective {
                    genus,
                    species_binomial: sp_binomial,
                    current_coverage_ratio: coverage,
                    missing_loci: missing,
                    priority_score: coverage,
                });
            }
        }
    }

    // Enforce hierarchical OTU rank rules (e.g., `OtuRank::Genus` requires `>= min_species_per_genus`
    // separate species terminals per Genus).
    let mut by_genus: BTreeMap<String, Vec<AssembledSpeciesOtu>> = BTreeMap::new();
    for otu in draft_species_otus {
        by_genus.entry(otu.genus.clone()).or_default().push(otu);
    }

    let mut retained_otus = Vec::new();
    let mut dropped_singleton_genera = Vec::new();
    let mut retained_genera_count = 0usize;

    let min_sp_req = match cfg.otu_rank {
        OtuRank::Genus => cfg.min_species_per_genus.max(1),
        _ => 1,
    };

    for (genus, mut sp_list) in by_genus {
        if sp_list.len() >= min_sp_req {
            retained_genera_count += 1;
            // Sort by coverage * neural confidence descending
            sp_list.sort_by(|a, b| {
                let sa = a.locus_coverage_ratio * a.mean_neural_confidence;
                let sb = b.locus_coverage_ratio * b.mean_neural_confidence;
                sb.partial_cmp(&sa).unwrap_or(std::cmp::Ordering::Equal)
            });
            retained_otus.extend(sp_list);
        } else {
            dropped_singleton_genera.push((genus, sp_list.len()));
        }
    }

    // Trim to `target_otu_count` while preserving `>= min_species_per_genus` integrity
    if retained_otus.len() > cfg.target_otu_count && cfg.target_otu_count >= min_sp_req {
        retained_otus.truncate(cfg.target_otu_count);
    }

    let matrix_occupancy = if retained_otus.is_empty() {
        0.0
    } else {
        let filled_cells: usize = retained_otus.iter().map(|o| o.loci_filled).sum();
        filled_cells as f64 / ((retained_otus.len() * loci_target) as f64)
    };

    gap_directives.sort_by(|a, b| {
        b.priority_score
            .partial_cmp(&a.priority_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    OtuValuationReport {
        retained_otus,
        dropped_low_coverage_species,
        dropped_singleton_genera,
        gap_directives,
        matrix_occupancy,
        retained_genera_count,
    }
}
