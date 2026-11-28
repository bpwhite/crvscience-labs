//! Coding-ortholog QC, genetic-code translation, RefSeq `db=gene` acquisition,
//! and post-MSA alignment row integrity gates shared across `bio-suite-seqnet`
//! and `bio-suite-supervisor`.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::time::Duration;

fn default_true() -> bool {
    true
}
fn default_genetic_code() -> u8 {
    11
}
fn default_min_bp() -> usize {
    300
}
fn default_max_bp() -> usize {
    4500
}
fn default_min_aa() -> usize {
    80
}
fn default_target_n() -> usize {
    12
}
fn default_min_anchor_sim() -> f64 {
    0.15
}
fn default_min_row_cov() -> f64 {
    0.45
}
fn default_min_cons_id() -> f64 {
    0.18
}
fn default_min_cons_d6() -> f64 {
    0.35
}
fn default_pos_max() -> usize {
    usize::MAX
}

/// Source strategy for functional coding ortholog acquisition.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum CodingOrthologSource {
    #[default]
    RefSeqGene,
    EntrezNuccore,
    PinnedAccessions,
}

/// Required amino-acid motif rule (supports `.` single-residue wildcards and OR alternatives).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteinMotifRule {
    pub name: String,
    #[serde(default)]
    pub pos_min: usize,
    #[serde(default = "default_pos_max")]
    pub pos_max: usize,
    pub patterns: Vec<String>,
}

/// Reusable configuration for frame/motif/homology auditing of a protein-coding ortholog cohort.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CodingOrthologConfig {
    pub cache_key: String,
    #[serde(default = "default_genetic_code")]
    pub genetic_code_table: u8,
    #[serde(default = "default_true")]
    pub preserve_frame0: bool,
    #[serde(default = "default_min_bp")]
    pub min_bp: usize,
    #[serde(default = "default_max_bp")]
    pub max_bp: usize,
    #[serde(default)]
    pub max_domain_bp: Option<usize>,
    #[serde(default = "default_min_aa")]
    pub min_aa_len: usize,
    #[serde(default)]
    pub max_internal_stops: usize,
    #[serde(default = "default_target_n")]
    pub target_n: usize,
    #[serde(default)]
    pub pin_species_substr: Option<String>,
    #[serde(default)]
    pub blocked_organisms: Vec<String>,
    #[serde(default)]
    pub required_motifs: Vec<ProteinMotifRule>,
    #[serde(default = "default_min_anchor_sim")]
    pub min_anchor_aa_kmer_sim: f64,
    #[serde(default = "default_min_row_cov")]
    pub min_row_coverage: f64,
    #[serde(default = "default_min_cons_id")]
    pub min_consensus_aa_identity: f64,
    #[serde(default = "default_min_cons_d6")]
    pub min_consensus_dayhoff6_sim: f64,
}

impl Default for CodingOrthologConfig {
    fn default() -> Self {
        Self {
            cache_key: String::new(),
            genetic_code_table: default_genetic_code(),
            preserve_frame0: true,
            min_bp: default_min_bp(),
            max_bp: default_max_bp(),
            max_domain_bp: None,
            min_aa_len: default_min_aa(),
            max_internal_stops: 0,
            target_n: default_target_n(),
            pin_species_substr: None,
            blocked_organisms: Vec::new(),
            required_motifs: Vec::new(),
            min_anchor_aa_kmer_sim: default_min_anchor_sim(),
            min_row_coverage: default_min_row_cov(),
            min_consensus_aa_identity: default_min_cons_id(),
            min_consensus_dayhoff6_sim: default_min_cons_d6(),
        }
    }
}

/// Candidate coding sequence evaluated by `audit_coding_ortholog_cohort`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CodingOrthologCandidate {
    pub accession: String,
    pub organism: String,
    pub title: String,
    #[serde(default)]
    pub subname: String,
    pub sequence_raw: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct LegacyCachedSummary {
    pub accession: String,
    pub organism: String,
    pub title: String,
    #[serde(default)]
    pub subname: String,
    #[serde(default)]
    pub slen: usize,
}

/// Quarantine event emitted by coding-ortholog or post-MSA alignment row gates.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrthologQuarantineHit {
    pub id: String,
    pub accession: String,
    pub organism: String,
    pub reason: String,
    pub detail: String,
}

/// Aggregate result returned by coding-ortholog QC and RefSeq `db=gene` cohort acquisition.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CodingOrthologAuditReport {
    pub accepted: Vec<CodingOrthologCandidate>,
    pub quarantined: Vec<OrthologQuarantineHit>,
}

pub fn clean_acgt(seq: &str) -> String {
    seq.bytes()
        .filter_map(|b| match b.to_ascii_uppercase() {
            b'A' | b'C' | b'G' | b'T' => Some(b as char),
            b'U' => Some('T'),
            _ => None,
        })
        .collect()
}

pub fn translate_codon(cod: &[u8], genetic_code_table: u8) -> char {
    if cod.len() < 3 {
        return 'X';
    }
    let c = [
        cod[0].to_ascii_uppercase(),
        cod[1].to_ascii_uppercase(),
        cod[2].to_ascii_uppercase(),
    ];
    match (&c, genetic_code_table) {
        (b"TAA" | b"TAG", _) => '*',
        (b"TGA", 2 | 5) => 'W',
        (b"TGA", _) => '*',
        (b"ATA", 2 | 5) => 'M',
        (b"AGA" | b"AGG", 2) => '*',
        (b"AGA" | b"AGG", 5) => 'S',
        (b"TTT" | b"TTC", _) => 'F',
        (b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG", _) => 'L',
        (b"ATT" | b"ATC" | b"ATA", _) => 'I',
        (b"ATG", _) => 'M',
        (b"GTT" | b"GTC" | b"GTA" | b"GTG", _) => 'V',
        (b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC", _) => 'S',
        (b"CCT" | b"CCC" | b"CCA" | b"CCG", _) => 'P',
        (b"ACT" | b"ACC" | b"ACA" | b"ACG", _) => 'T',
        (b"GCT" | b"GCC" | b"GCA" | b"GCG", _) => 'A',
        (b"TAT" | b"TAC", _) => 'Y',
        (b"CAT" | b"CAC", _) => 'H',
        (b"CAA" | b"CAG", _) => 'Q',
        (b"AAT" | b"AAC", _) => 'N',
        (b"AAA" | b"AAG", _) => 'K',
        (b"GAT" | b"GAC", _) => 'D',
        (b"GAA" | b"GAG", _) => 'E',
        (b"TGT" | b"TGC", _) => 'C',
        (b"TGG", _) => 'W',
        (b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG", _) => 'R',
        (b"GGT" | b"GGC" | b"GGA" | b"GGG", _) => 'G',
        _ => 'X',
    }
}

pub fn trim_frame0_cds_to_stop(seq: &str, max_bp: usize, genetic_code_table: u8) -> String {
    let clean = clean_acgt(seq);
    let cap = max_bp - (max_bp % 3);
    let limit = clean.len().min(cap);
    let b = clean.as_bytes();
    let mut end = 0usize;
    while end + 3 <= limit {
        if translate_codon(&b[end..end + 3], genetic_code_table) == '*' {
            break;
        }
        end += 3;
    }
    clean[..end].to_string()
}

pub fn translate_cds_with_table(nt: &str, genetic_code_table: u8, max_codons: usize) -> String {
    let b = nt.as_bytes();
    let limit = (b.len() / 3).min(max_codons);
    let mut aa = String::with_capacity(limit);
    for c in 0..limit {
        aa.push(translate_codon(&b[c * 3..c * 3 + 3], genetic_code_table));
    }
    aa
}

/// Extracts the best sense ORF; when `preserve_frame0` is true (e.g., RefSeq/GenBank annotated
/// bacterial CDS that may initiate at `GTG`/`TTG`/`ATT`), translates strictly from codon 0.
pub fn extract_best_sense_orf_with_table(
    seq: &str,
    max_bp: usize,
    genetic_code_table: u8,
    preserve_frame0: bool,
) -> String {
    if preserve_frame0 {
        return trim_frame0_cds_to_stop(seq, max_bp, genetic_code_table);
    }
    let clean = clean_acgt(seq);
    let b = clean.as_bytes();
    if b.len() < 30 {
        return trim_frame0_cds_to_stop(&clean, max_bp, genetic_code_table);
    }
    let mut best_score = 0usize;
    let mut best_span = (0usize, 0usize);

    for frame in 0..3 {
        let mut run_start = frame;
        let mut first_atg: Option<usize> = None;
        let mut pos = frame;
        while pos + 3 <= b.len() {
            let cod = &b[pos..pos + 3];
            let is_stop = translate_codon(cod, genetic_code_table) == '*';
            if !is_stop && first_atg.is_none() && cod == b"ATG" {
                first_atg = Some(pos);
            }
            if is_stop || pos + 6 > b.len() {
                let run_end = if is_stop { pos } else { pos + 3 };
                let (cand_s, cand_len, bonus) = if let Some(atg) = first_atg {
                    let atg_codons = (run_end.saturating_sub(atg)) / 3;
                    let raw_codons = (run_end.saturating_sub(run_start)) / 3;
                    if atg_codons * 2 >= raw_codons {
                        (atg, run_end.saturating_sub(atg), atg_codons + 25)
                    } else {
                        (run_start, run_end.saturating_sub(run_start), raw_codons)
                    }
                } else {
                    (
                        run_start,
                        run_end.saturating_sub(run_start),
                        (run_end.saturating_sub(run_start)) / 3,
                    )
                };
                if bonus > best_score && cand_len >= 60 {
                    best_score = bonus;
                    best_span = (cand_s, run_end);
                }
                run_start = pos + 3;
                first_atg = None;
            }
            pos += 3;
        }
    }
    let (s, e) = if best_span.1 > best_span.0 {
        best_span
    } else {
        (0, b.len() - (b.len() % 3))
    };
    let cap = max_bp - (max_bp % 3);
    let end = (s + cap).min(e);
    let trimmed_len = (end.saturating_sub(s)) - ((end.saturating_sub(s)) % 3);
    clean[s..s + trimmed_len].to_string()
}

pub fn dayhoff6_class(aa: u8) -> u8 {
    match aa.to_ascii_uppercase() {
        b'A' | b'G' | b'P' | b'S' | b'T' => 0,
        b'D' | b'E' | b'N' | b'Q' => 1,
        b'H' | b'K' | b'R' => 2,
        b'I' | b'L' | b'M' | b'V' => 3,
        b'F' | b'W' | b'Y' => 4,
        b'C' => 5,
        _ => 6,
    }
}

fn dayhoff6_3mer_profile(aa: &str) -> [f32; 216] {
    let mut vec = [0.0f32; 216];
    let classes: Vec<u8> = aa
        .bytes()
        .map(dayhoff6_class)
        .filter(|&c| c < 6)
        .collect();
    for w in classes.windows(3) {
        let idx = (w[0] as usize) * 36 + (w[1] as usize) * 6 + (w[2] as usize);
        vec[idx] += 1.0;
    }
    let norm = vec.iter().map(|v| v * v).sum::<f32>().sqrt();
    if norm > 1e-8 {
        for v in &mut vec {
            *v /= norm;
        }
    }
    vec
}

fn cosine_216(a: &[f32; 216], b: &[f32; 216]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(&x, &y)| (x as f64) * (y as f64))
        .sum::<f64>()
        .clamp(0.0, 1.0)
}

pub fn dayhoff6_3mer_cosine(aa_a: &str, aa_b: &str) -> f64 {
    let pa = dayhoff6_3mer_profile(aa_a);
    let pb = dayhoff6_3mer_profile(aa_b);
    cosine_216(&pa, &pb)
}

fn wildcard_match(window: &[u8], pattern: &[u8]) -> bool {
    if window.len() != pattern.len() {
        return false;
    }
    window.iter().zip(pattern.iter()).all(|(&w, &p)| {
        let pu = p.to_ascii_uppercase();
        pu == b'.' || pu == w.to_ascii_uppercase()
    })
}

pub fn matches_protein_motif_rules(aa: &str, rules: &[ProteinMotifRule]) -> Result<(), String> {
    let b = aa.as_bytes();
    for rule in rules {
        let mut matched = false;
        for pat in &rule.patterns {
            let pb = pat.as_bytes();
            if pb.is_empty() || pb.len() > b.len() || rule.pos_min + pb.len() > b.len() {
                continue;
            }
            let max_start = b.len().saturating_sub(pb.len()).min(rule.pos_max);
            if rule.pos_min <= max_start
                && (rule.pos_min..=max_start).any(|idx| wildcard_match(&b[idx..idx + pb.len()], pb))
            {
                matched = true;
                break;
            }
        }
        if !matched {
            return Err(format!(
                "missing motif '{}' in residues {}..{}",
                rule.name, rule.pos_min, rule.pos_max
            ));
        }
    }
    Ok(())
}

fn norm_sp2(org: &str) -> String {
    org.split(|c: char| !c.is_ascii_alphanumeric())
        .filter(|t| !t.is_empty())
        .take(2)
        .map(|t| t.to_ascii_lowercase())
        .collect::<Vec<_>>()
        .join(" ")
}

fn count_internal_stops(aa: &str) -> usize {
    let bytes = aa.as_bytes();
    if bytes.len() <= 1 {
        return 0;
    }
    bytes[..bytes.len() - 1]
        .iter()
        .filter(|&&ch| ch == b'*')
        .count()
}

pub fn audit_coding_ortholog_cohort(
    candidates: &[CodingOrthologCandidate],
    cfg: &CodingOrthologConfig,
) -> CodingOrthologAuditReport {
    let mut survivors: Vec<(CodingOrthologCandidate, String)> = Vec::new();
    let mut quarantined: Vec<OrthologQuarantineHit> = Vec::new();
    let blocked_low: Vec<String> = cfg
        .blocked_organisms
        .iter()
        .map(|s| s.to_ascii_lowercase())
        .collect();

    for cand in candidates {
        let org_low = cand.organism.to_ascii_lowercase();
        if blocked_low
            .iter()
            .any(|b| !b.is_empty() && org_low.contains(b))
        {
            quarantined.push(OrthologQuarantineHit {
                id: cand.accession.clone(),
                accession: cand.accession.clone(),
                organism: cand.organism.clone(),
                reason: "BLOCKED_ORGANISM".to_string(),
                detail: "matched organism blocklist".to_string(),
            });
            continue;
        }
        let clean = clean_acgt(&cand.sequence_raw);
        if clean.len() < cfg.min_bp || clean.len() > cfg.max_bp {
            quarantined.push(OrthologQuarantineHit {
                id: cand.accession.clone(),
                accession: cand.accession.clone(),
                organism: cand.organism.clone(),
                reason: "SEQUENCE_LENGTH_OUT_OF_BOUNDS".to_string(),
                detail: format!("{} bp not in {}..{}", clean.len(), cfg.min_bp, cfg.max_bp),
            });
            continue;
        }
        let orf_cap = cfg.max_domain_bp.unwrap_or(cfg.max_bp);
        let orf_nt = extract_best_sense_orf_with_table(
            &clean,
            orf_cap,
            cfg.genetic_code_table,
            cfg.preserve_frame0,
        );
        let aa = translate_cds_with_table(&orf_nt, cfg.genetic_code_table, orf_nt.len() / 3);
        if aa.len() < cfg.min_aa_len {
            quarantined.push(OrthologQuarantineHit {
                id: cand.accession.clone(),
                accession: cand.accession.clone(),
                organism: cand.organism.clone(),
                reason: "TRUNCATED_OR_EARLY_STOP".to_string(),
                detail: format!("ORF {} aa < {} aa", aa.len(), cfg.min_aa_len),
            });
            continue;
        }
        let stops = count_internal_stops(&aa);
        if stops > cfg.max_internal_stops {
            quarantined.push(OrthologQuarantineHit {
                id: cand.accession.clone(),
                accession: cand.accession.clone(),
                organism: cand.organism.clone(),
                reason: "INTERNAL_STOP_CODON_NUMT".to_string(),
                detail: format!("{} internal stops > {}", stops, cfg.max_internal_stops),
            });
            continue;
        }
        if let Err(why) = matches_protein_motif_rules(&aa, &cfg.required_motifs) {
            quarantined.push(OrthologQuarantineHit {
                id: cand.accession.clone(),
                accession: cand.accession.clone(),
                organism: cand.organism.clone(),
                reason: "MISSING_ORTHOLOG_MOTIF".to_string(),
                detail: why,
            });
            continue;
        }
        survivors.push((cand.clone(), aa));
    }

    if survivors.len() < 3 || cfg.min_anchor_aa_kmer_sim <= 0.0 {
        return CodingOrthologAuditReport {
            accepted: survivors.into_iter().map(|(c, _)| c).collect(),
            quarantined,
        };
    }

    let profiles: Vec<[f32; 216]> = survivors
        .iter()
        .map(|(_, aa)| dayhoff6_3mer_profile(aa))
        .collect();
    let mut centroid = [0.0f32; 216];
    for p in &profiles {
        for (i, &v) in p.iter().enumerate() {
            centroid[i] += v;
        }
    }
    let cnorm = centroid.iter().map(|v| v * v).sum::<f32>().sqrt();
    if cnorm > 1e-8 {
        for v in &mut centroid {
            *v /= cnorm;
        }
    }

    let anchor_idx = cfg
        .pin_species_substr
        .as_deref()
        .and_then(|pin| {
            let p_low = pin.to_ascii_lowercase();
            survivors
                .iter()
                .position(|(c, _)| c.organism.to_ascii_lowercase().contains(&p_low))
        })
        .unwrap_or_else(|| {
            profiles
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| {
                    cosine_216(a, &centroid)
                        .partial_cmp(&cosine_216(b, &centroid))
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(i, _)| i)
                .unwrap_or(0)
        });

    let mut reference = centroid;
    if cfg.pin_species_substr.is_some() {
        let anc = profiles[anchor_idx];
        for i in 0..216 {
            reference[i] = 0.5 * centroid[i] + 0.5 * anc[i];
        }
        let rnorm = reference.iter().map(|v| v * v).sum::<f32>().sqrt();
        if rnorm > 1e-8 {
            for v in &mut reference {
                *v /= rnorm;
            }
        }
    }

    let mut accepted = Vec::with_capacity(survivors.len());
    for (idx, (cand, _aa)) in survivors.into_iter().enumerate() {
        let sim_ref = cosine_216(&profiles[idx], &reference);
        let sim_cent = cosine_216(&profiles[idx], &centroid);
        if sim_ref < cfg.min_anchor_aa_kmer_sim && sim_cent < cfg.min_anchor_aa_kmer_sim {
            quarantined.push(OrthologQuarantineHit {
                id: cand.accession.clone(),
                accession: cand.accession.clone(),
                organism: cand.organism.clone(),
                reason: "NON_HOMOLOGOUS_PROTEIN_OUTLIER".to_string(),
                detail: format!(
                    "Dayhoff6 3-mer cosine {:.3} < {:.3}",
                    sim_ref.max(sim_cent),
                    cfg.min_anchor_aa_kmer_sim
                ),
            });
        } else {
            accepted.push(cand);
        }
    }

    CodingOrthologAuditReport {
        accepted,
        quarantined,
    }
}

/// Post-MSA codon alignment row integrity check across `(id, aa_row, nt_row)` tuples.
pub fn audit_aligned_codon_rows(
    rows: &[(String, String, String)],
    min_row_coverage: f64,
    min_consensus_aa_identity: f64,
    min_consensus_dayhoff6_sim: f64,
) -> (Vec<usize>, Vec<OrthologQuarantineHit>) {
    if rows.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let cols = rows.iter().map(|(_, aa, _)| aa.len()).max().unwrap_or(0);
    if cols == 0 {
        return (Vec::new(), Vec::new());
    }

    let mut consensus: Vec<u8> = vec![b'-'; cols];
    for c in 0..cols {
        let mut counts = [0u16; 26];
        let mut best_aa = b'-';
        let mut best_cnt = 0u16;
        for (_, aa, _) in rows {
            if let Some(&b) = aa.as_bytes().get(c) {
                let bu = b.to_ascii_uppercase();
                if bu.is_ascii_uppercase() && bu != b'X' {
                    let idx = (bu - b'A') as usize;
                    counts[idx] = counts[idx].saturating_add(1);
                    if counts[idx] > best_cnt {
                        best_cnt = counts[idx];
                        best_aa = bu;
                    }
                }
            }
        }
        consensus[c] = best_aa;
    }

    let mut kept = Vec::with_capacity(rows.len());
    let mut quarantined = Vec::new();
    for (ri, (id, aa, _nt)) in rows.iter().enumerate() {
        let ab = aa.as_bytes();
        let mut non_gap = 0usize;
        let mut compared = 0usize;
        let mut aa_matches = 0usize;
        let mut d6_matches = 0usize;
        for c in 0..cols {
            let r_ch = ab.get(c).copied().unwrap_or(b'-').to_ascii_uppercase();
            if r_ch != b'-' && r_ch != b'?' {
                non_gap += 1;
                let cons_ch = consensus[c];
                if cons_ch != b'-' && r_ch != b'X' {
                    compared += 1;
                    if r_ch == cons_ch {
                        aa_matches += 1;
                    }
                    let rc = dayhoff6_class(r_ch);
                    let cc = dayhoff6_class(cons_ch);
                    if rc < 6 && rc == cc {
                        d6_matches += 1;
                    }
                }
            }
        }
        let cov = non_gap as f64 / cols as f64;
        if cov < min_row_coverage {
            quarantined.push(OrthologQuarantineHit {
                id: id.clone(),
                accession: id.clone(),
                organism: String::new(),
                reason: "LOW_ALIGNMENT_CORE_COVERAGE".to_string(),
                detail: format!("coverage {:.3} < {:.3}", cov, min_row_coverage),
            });
            continue;
        }
        let aa_id = if compared > 0 {
            aa_matches as f64 / compared as f64
        } else {
            0.0
        };
        let d6_sim = if compared > 0 {
            d6_matches as f64 / compared as f64
        } else {
            0.0
        };
        if compared >= 20
            && (aa_id < min_consensus_aa_identity || d6_sim < min_consensus_dayhoff6_sim)
        {
            quarantined.push(OrthologQuarantineHit {
                id: id.clone(),
                accession: id.clone(),
                organism: String::new(),
                reason: "NON_HOMOLOGOUS_ALIGNED_ROW".to_string(),
                detail: format!(
                    "AA id {:.3} (min {:.3}) or Dayhoff6 sim {:.3} (min {:.3})",
                    aa_id, min_consensus_aa_identity, d6_sim, min_consensus_dayhoff6_sim
                ),
            });
            continue;
        }
        // Sliding 20-codon local window homology gate to catch chimeric or locally frame-shifted tracts
        if cols >= 32 {
            let win = 20usize;
            let mut local_fail: Option<(usize, f64, f64)> = None;
            let margin = (cols / 10).min(8);
            let hi = cols.saturating_sub(win + margin).max(margin);
            for s0 in (margin..=hi).step_by(5) {
                if s0 + win > cols {
                    break;
                }
                let mut w_cmp = 0usize;
                let mut w_id = 0usize;
                let mut w_d6 = 0usize;
                for c in s0..(s0 + win) {
                    let r_ch = ab.get(c).copied().unwrap_or(b'-').to_ascii_uppercase();
                    let cons_ch = consensus[c];
                    if r_ch != b'-' && r_ch != b'?' && r_ch != b'X' && cons_ch != b'-' {
                        w_cmp += 1;
                        if r_ch == cons_ch {
                            w_id += 1;
                        }
                        let rc = dayhoff6_class(r_ch);
                        let cc = dayhoff6_class(cons_ch);
                        if rc < 6 && rc == cc {
                            w_d6 += 1;
                        }
                    }
                }
                if w_cmp >= 15 {
                    let loc_id = w_id as f64 / w_cmp as f64;
                    let loc_d6 = w_d6 as f64 / w_cmp as f64;
                    let req_d6 = (min_consensus_dayhoff6_sim * 0.55).clamp(0.22, 0.35);
                    if loc_id < 0.12 && loc_d6 < req_d6 {
                        local_fail = Some((s0 + 1, loc_id, loc_d6));
                        break;
                    }
                }
            }
            if let Some((col1, lid, ld6)) = local_fail {
                quarantined.push(OrthologQuarantineHit {
                    id: id.clone(),
                    accession: id.clone(),
                    organism: String::new(),
                    reason: "LOCAL_WINDOW_MISALIGNMENT_OR_CHIMERA".to_string(),
                    detail: format!(
                        "20-codon window at aa {} failed local score (AA id {:.2}, Dayhoff6 {:.2})",
                        col1, lid, ld6
                    ),
                });
                continue;
            }
        }
        kept.push(ri);
    }

    (kept, quarantined)
}

fn url_encode(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 32);
    for &b in s.as_bytes() {
        if b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_' | b'.' | b'~') {
            out.push(b as char);
        } else if b == b' ' {
            out.push_str("%20");
        } else {
            out.push_str(&format!("%{:02X}", b));
        }
    }
    out
}

fn entrez_get_polite(url: &str) -> Option<String> {
    let api_key = std::env::var("NCBI_API_KEY")
        .ok()
        .map(|v| v.trim().to_string())
        .filter(|v| !v.is_empty());
    let full = if let Some(ref k) = api_key {
        let sep = if url.contains('?') { '&' } else { '?' };
        format!("{}{}api_key={}", url, sep, url_encode(k))
    } else {
        url.to_string()
    };
    let rate_ms = if api_key.is_some() { 140 } else { 360 };
    for attempt in 0..4u64 {
        let wait = if attempt == 0 {
            rate_ms
        } else {
            900u64 << (attempt - 1)
        };
        std::thread::sleep(Duration::from_millis(wait));
        if let Ok(resp) = ureq::get(&full)
            .set(
                "User-Agent",
                "crvscience-bio-suite-seqnet/0.1.0 (+https://crvscience.com/labs; GPL-3.0)",
            )
            .timeout(Duration::from_secs(25))
            .call()
        {
            if let Ok(txt) = resp.into_string() {
                return Some(txt);
            }
        }
    }
    None
}

fn load_cached_candidates(cache_path: &Path) -> Option<Vec<CodingOrthologCandidate>> {
    let raw = fs::read_to_string(cache_path).ok()?;
    if let Ok(cands) = serde_json::from_str::<Vec<CodingOrthologCandidate>>(&raw) {
        if !cands.is_empty() {
            return Some(cands);
        }
    }
    if let Ok(legacy) = serde_json::from_str::<Vec<(LegacyCachedSummary, String)>>(&raw) {
        if !legacy.is_empty() {
            return Some(
                legacy
                    .into_iter()
                    .map(|(s, seq)| CodingOrthologCandidate {
                        accession: s.accession,
                        organism: s.organism,
                        title: s.title,
                        subname: s.subname,
                        sequence_raw: seq,
                    })
                    .collect(),
            );
        }
    }
    None
}

/// Fetches RefSeq `db=gene` orthologs, validates via `audit_coding_ortholog_cohort`,
/// caches verified candidates at `{cache_dir}/{cfg.cache_key}.json`, and returns audit output.
pub fn fetch_refseq_gene_orthologs(
    queries: &[String],
    cfg: &CodingOrthologConfig,
    cache_dir: &Path,
) -> CodingOrthologAuditReport {
    let _ = fs::create_dir_all(cache_dir);
    let cache_path = cache_dir.join(format!("{}.json", cfg.cache_key));
    if let Some(cached) = load_cached_candidates(&cache_path) {
        let mut rep = audit_coding_ortholog_cohort(&cached, cfg);
        rep.accepted.truncate(cfg.target_n.max(1));
        if !rep.accepted.is_empty() {
            return rep;
        }
    }

    let blocked_low: Vec<String> = cfg
        .blocked_organisms
        .iter()
        .map(|s| s.to_ascii_lowercase())
        .collect();
    let mut by_sp: HashMap<String, (u8, String, u64, u64, u8, String, String)> = HashMap::new();

    for q in queries {
        let esearch = format!(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&retmode=json&retmax=120&tool=bio-suite-seqnet&email=labs@crvscience.com&term={}",
            url_encode(q)
        );
        let Some(search_txt) = entrez_get_polite(&esearch) else { continue };
        let Ok(search_json) = serde_json::from_str::<serde_json::Value>(&search_txt) else { continue };
        let ids: Vec<String> = search_json
            .get("esearchresult")
            .and_then(|r| r.get("idlist"))
            .and_then(|v| v.as_array())
            .map(|arr| {
                arr.iter()
                    .filter_map(|v| v.as_str().map(|s| s.to_string()))
                    .collect()
            })
            .unwrap_or_default();
        for chunk in ids.chunks(60) {
            let esum = format!(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&retmode=json&tool=bio-suite-seqnet&email=labs@crvscience.com&id={}",
                chunk.join(",")
            );
            let Some(sum_txt) = entrez_get_polite(&esum) else { continue };
            let Ok(sum_json) = serde_json::from_str::<serde_json::Value>(&sum_txt) else { continue };
            let Some(res_obj) = sum_json.get("result").and_then(|v| v.as_object()) else { continue };
            for id in chunk {
                let Some(it) = res_obj.get(id) else { continue };
                let org = it
                    .get("organism")
                    .and_then(|o| o.get("scientificname"))
                    .and_then(|v| v.as_str())
                    .unwrap_or("");
                if org.is_empty()
                    || blocked_low
                        .iter()
                        .any(|b| !b.is_empty() && org.to_ascii_lowercase().contains(b))
                {
                    continue;
                }
                let sp2 = norm_sp2(org);
                if sp2.is_empty() {
                    continue;
                }
                let Some(gi) = it
                    .get("genomicinfo")
                    .and_then(|v| v.as_array())
                    .and_then(|arr| arr.first())
                else {
                    continue;
                };
                let acc = gi.get("chraccver").and_then(|v| v.as_str()).unwrap_or("");
                let (Some(chrstart), Some(chrstop)) = (
                    gi.get("chrstart").and_then(|v| v.as_u64()),
                    gi.get("chrstop").and_then(|v| v.as_u64()),
                ) else {
                    continue;
                };
                if acc.is_empty() {
                    continue;
                }
                let lo = chrstart.min(chrstop);
                let hi = chrstart.max(chrstop);
                let len = (hi - lo + 1) as usize;
                if len < cfg.min_bp || len > cfg.max_bp {
                    continue;
                }
                let prio: u8 = if acc.starts_with("NC_") {
                    3
                } else if acc.starts_with("NZ_CP") || acc.starts_with("NZ_AP") || acc.starts_with("NZ_LR") {
                    2
                } else if acc.starts_with("NZ_") {
                    1
                } else {
                    0
                };
                let gname = it.get("name").and_then(|v| v.as_str()).unwrap_or("");
                let gdesc = it.get("description").and_then(|v| v.as_str()).unwrap_or("");
                let title = format!("{} {} ({})", org, gdesc, gname);
                let strand: u8 = if chrstart <= chrstop { 1 } else { 2 };
                let replace = by_sp.get(&sp2).map(|ex| prio > ex.0).unwrap_or(true);
                if replace {
                    by_sp.insert(
                        sp2,
                        (prio, acc.to_string(), lo + 1, hi + 1, strand, org.to_string(), title),
                    );
                }
            }
        }
    }

    let mut ranked: Vec<(String, u8, String, u64, u64, u8, String, String)> = by_sp
        .into_iter()
        .map(|(sp2, (p, a, s, e, st, o, t))| (sp2, p, a, s, e, st, o, t))
        .collect();
    ranked.sort_by(|a, b| {
        let pin_a = cfg
            .pin_species_substr
            .as_deref()
            .map(|p| a.6.to_ascii_lowercase().contains(&p.to_ascii_lowercase()))
            .unwrap_or(false);
        let pin_b = cfg
            .pin_species_substr
            .as_deref()
            .map(|p| b.6.to_ascii_lowercase().contains(&p.to_ascii_lowercase()))
            .unwrap_or(false);
        pin_b
            .cmp(&pin_a)
            .then_with(|| b.1.cmp(&a.1))
            .then_with(|| a.0.cmp(&b.0))
    });

    let mut fetched = Vec::with_capacity(cfg.target_n * 2);
    for (_sp2, _prio, acc, s1, e1, st, org, title) in ranked {
        if fetched.len() >= cfg.target_n * 2 {
            break;
        }
        let efetch = format!(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={}&seq_start={}&seq_stop={}&strand={}&rettype=fasta&retmode=text&tool=bio-suite-seqnet&email=labs@crvscience.com",
            url_encode(&acc),
            s1,
            e1,
            st
        );
        let Some(fa) = entrez_get_polite(&efetch) else { continue };
        let mut body = String::with_capacity(fa.len());
        for line in fa.lines() {
            if !line.starts_with('>') {
                body.push_str(line.trim());
            }
        }
        let clean = clean_acgt(&body);
        if !clean.is_empty() {
            fetched.push(CodingOrthologCandidate {
                accession: acc,
                organism: org,
                title,
                subname: String::new(),
                sequence_raw: clean,
            });
        }
    }

    let mut rep = audit_coding_ortholog_cohort(&fetched, cfg);
    rep.accepted.truncate(cfg.target_n.max(1));
    if !rep.accepted.is_empty() {
        let legacy: Vec<(LegacyCachedSummary, String)> = rep
            .accepted
            .iter()
            .map(|c| {
                (
                    LegacyCachedSummary {
                        accession: c.accession.clone(),
                        organism: c.organism.clone(),
                        title: c.title.clone(),
                        subname: c.subname.clone(),
                        slen: c.sequence_raw.len(),
                    },
                    c.sequence_raw.clone(),
                )
            })
            .collect();
        if let Ok(json_s) = serde_json::to_string(&legacy) {
            let _ = fs::write(&cache_path, json_s);
        }
    }
    rep
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aa_to_codons(aa: &str, start_gtg: bool) -> String {
        let mut nt = String::with_capacity(aa.len() * 3);
        for (i, ch) in aa.bytes().enumerate() {
            let c = if i == 0 && start_gtg {
                "GTG"
            } else {
                match ch {
                    b'A' => "GCT",
                    b'C' => "TGT",
                    b'D' => "GAT",
                    b'E' => "GAA",
                    b'F' => "TTT",
                    b'G' => "GGT",
                    b'H' => "CAT",
                    b'I' => "ATT",
                    b'K' => "AAA",
                    b'L' => "CTG",
                    b'M' => "ATG",
                    b'N' => "AAT",
                    b'P' => "CCG",
                    b'Q' => "CAG",
                    b'R' => "CGT",
                    b'S' => "TCT",
                    b'T' => "ACC",
                    b'V' => "GTG",
                    b'W' => "TGG",
                    b'Y' => "TAT",
                    _ => "GCT",
                }
            };
            nt.push_str(c);
        }
        nt
    }

    #[test]
    fn frame0_preserves_gtg_start_without_internal_atg_shift() {
        // Frame 0 begins at GTG followed by an out-of-frame ATG at position +4 (GT-GAT-G -> ATG shifted by 1 bp)
        let mut nt = "GTGGTGATGAAAGCTCTGGAAACCATTCCGCAGCCG".to_string();
        nt.push_str(&"GCTGAACTGAAAGGTATTCCGCTGACC".repeat(12));
        let orf = extract_best_sense_orf_with_table(&nt, 1500, 11, true);
        assert!(orf.starts_with("GTGGTGATG"), "Expected GTG start frame 0, got {}", &orf[..12]);
        let aa = translate_cds_with_table(&orf, 11, orf.len() / 3);
        assert!(aa.starts_with("VVMKALETIPQP"));
        assert!(!aa.contains('*'));
    }

    #[test]
    fn cyp102a_motif_and_dayhoff6_guard_accepts_orthologs_and_rejects_operon() {
        let mut base = vec![b'A'; 455];
        for (i, &b) in b"MTIKEMPQPKTFGELKNLPL".iter().enumerate() {
            base[i] = b;
        }
        // Plant I-helix AGHET at 264 and Cys-loop FGNGQRACIG at 393
        for (i, &b) in b"AGHETTS".iter().enumerate() {
            base[264 + i] = b;
        }
        for (i, &b) in b"FGNGQRACIG".iter().enumerate() {
            base[393 + i] = b;
        }
        // Realistic pattern entropy across body
        let fill = b"DEKLPQSTVYIARGNHWFC";
        for idx in 25..450 {
            if !(264..272).contains(&idx) && !(393..403).contains(&idx) {
                base[idx] = fill[(idx * 7 + 3) % fill.len()];
            }
        }
        let base_s = String::from_utf8(base.clone()).unwrap();
        let mut strep = base.clone();
        strep[395] = b'T';
        strep[397] = b'V';
        let strep_s = String::from_utf8(strep).unwrap();
        let unrelated = "WWWWWWWWWWYYYYYYYYYYMMMMMMMMMMCCCCCCCCCC".repeat(12);

        let cands = vec![
            CodingOrthologCandidate {
                accession: "NC_017138.1".into(),
                organism: "Bacillus megaterium WSH-002".into(),
                title: "P450 BM3".into(),
                subname: String::new(),
                sequence_raw: aa_to_codons(&base_s, false),
            },
            CodingOrthologCandidate {
                accession: "NC_000964.3".into(),
                organism: "Bacillus subtilis 168".into(),
                title: "cypD".into(),
                subname: String::new(),
                sequence_raw: aa_to_codons(&base_s, true),
            },
            CodingOrthologCandidate {
                accession: "NC_003155.5".into(),
                organism: "Streptomyces avermitilis".into(),
                title: "CYP102".into(),
                subname: String::new(),
                sequence_raw: aa_to_codons(&strep_s, false),
            },
            CodingOrthologCandidate {
                accession: "GQ334372.1".into(),
                organism: "Bacillus sp.".into(),
                title: "unrelated xylanase".into(),
                subname: String::new(),
                sequence_raw: aa_to_codons(&unrelated[..455], false),
            },
        ];
        let cfg = CodingOrthologConfig {
            cache_key: "test".into(),
            genetic_code_table: 11,
            preserve_frame0: true,
            min_bp: 1200,
            max_bp: 3600,
            max_domain_bp: Some(1395),
            min_aa_len: 420,
            pin_species_substr: Some("megaterium".into()),
            required_motifs: vec![
                ProteinMotifRule {
                    name: "P450_I_helix".into(),
                    pos_min: 225,
                    pos_max: 315,
                    patterns: vec!["AGHET".into(), "SGHET".into(), "AGFET".into()],
                },
                ProteinMotifRule {
                    name: "P450_Cys_loop".into(),
                    pos_min: 350,
                    pos_max: 435,
                    patterns: vec!["F..G...C.G".into()],
                },
            ],
            ..CodingOrthologConfig::default()
        };
        let rep = audit_coding_ortholog_cohort(&cands, &cfg);
        assert_eq!(rep.accepted.len(), 3);
        assert_eq!(rep.quarantined.len(), 1);
        assert_eq!(rep.quarantined[0].accession, "GQ334372.1");
    }

    #[test]
    fn post_msa_row_guard_rejects_all_gap_and_frameshift_rows() {
        let cons = "MKETGPIPQPKTFGPLGNLPLIDKDKPTLSLIKLAEEQGPIFQIHTPAGTTIVVSGHEL";
        let ortho = "MKETSPVPQPKTFGPLGNLPLLDKDKPTLSLIKLAEEQGPIFQIHTPAGTTIVVSGHEL";
        let all_gap = "-".repeat(cons.len());
        let shifted = "WWWWCCCCYYYYHHHHWWWWCCCCYYYYHHHHWWWWCCCCYYYYHHHHWWWWCCCCYYYYH";
        let rows = vec![
            ("FAM-04_1".into(), cons.into(), aa_to_codons(cons, false)),
            ("FAM-04_2".into(), cons.into(), aa_to_codons(cons, false)),
            ("FAM-04_3".into(), ortho.into(), aa_to_codons(ortho, false)),
            ("FAM-04_4".into(), ortho.into(), aa_to_codons(ortho, false)),
            ("KC839476.1".into(), all_gap.clone(), "---".repeat(cons.len())),
            ("DM173719.1".into(), shifted.into(), aa_to_codons(shifted, false)),
        ];
        let (kept, q) = audit_aligned_codon_rows(&rows, 0.45, 0.18, 0.35);
        assert_eq!(kept, vec![0, 1, 2, 3]);
        assert_eq!(q.len(), 2);
        assert_eq!(q[0].reason, "LOW_ALIGNMENT_CORE_COVERAGE");
        assert_eq!(q[1].reason, "NON_HOMOLOGOUS_ALIGNED_ROW");
    }
}
