//! Native Biosynthetic Gene Cluster (BGC) Detection Engine (`bio-suite bgc-mine` / `bgc-compare`).
//!
//! Calibrated head-to-head against `GECCO v0.11.0` (conditional-random-field BGC detector over
//! the full Pfam profile library) and scored against experimentally curated `MIBiG 4.0` loci.
//!
//! Architectural stages (all native Rust; no external binary is needed for detection):
//! 1. Plan-7 profile construction from real Pfam seed alignments (34 biosynthetic signature
//!    families: PKS, NRPS, terpene, RiPP, siderophore and tailoring domains) via the native
//!    Henikoff-weighted builder in [`crate::profile_hmm`].
//! 2. Native 4-stage ORF graph gene calling ([`crate::gene_caller`]).
//! 3. Parallel MSV-prefiltered Plan-7 domain scan with Pfam curated gathering (GA) cut-offs.
//! 4. Architecture-aware rule-based cluster calling: core-gene definitions that require
//!    multi-domain architectures (so fatty-acid synthase and acyl-CoA synthetase housekeeping
//!    genes are not mistaken for modular PKS/NRPS), plus windowed accessory-gene support.
//! 5. Head-to-head calibration harness versus `GECCO` and `MIBiG` ground truth.

use crate::gene_caller::{run_native_prokaryotic_gene_caller, translate_nt_to_aa};
use crate::profile_hmm::{
    build_native_plan7_from_msa, digitize_protein_seq, p7_bg_null_one, p7_gmsv,
    score_sequence_plan7_native, P7Hmm, P7Profile,
};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

const LN2: f32 = std::f32::consts::LN_2;
/// MSV bit-score prefilter (HMMER's MSV filter passes at roughly this level for true homologs).
const MSV_PREFILTER_BITS: f32 = 7.0;
/// Smoothing window (genes on each side) for accessory-gene support.
const SUPPORT_WINDOW_GENES: usize = 4;
/// Minimum smoothed support score required for a core gene to seed a cluster.
const SEED_MIN_SUPPORT: f32 = 1.8;
/// Gene-to-gene gap (bp) tolerated while extending a cluster.
const MAX_EXTEND_GAP_BP: usize = 5000;
/// Smoothed support needed for a non-biosynthetic neighbour to be absorbed into a cluster.
const EXTEND_MIN_SUPPORT: f32 = 0.9;
/// Flank padding (bp) applied to every called cluster.
const FLANK_PAD_BP: usize = 2000;

/// Biosynthetic signature domain families (Pfam accession, short name, accessory weight).
/// Core-gene roles are decided by domain *architecture* in [`classify_core_gene`], not by weight.
pub const BGC_FAMILIES: &[(&str, &str, f32)] = &[
    ("PF00109", "ketoacyl-synt", 0.8),
    ("PF02801", "Ketoacyl-synt_C", 0.3),
    ("PF00698", "Acyl_transf_1", 0.7),
    ("PF08659", "KR", 0.5),
    ("PF14765", "PS-DH", 0.6),
    ("PF00550", "PP-binding", 0.4),
    ("PF16197", "KAsynt_C_assoc", 0.3),
    ("PF00501", "AMP-binding", 0.6),
    ("PF00668", "Condensation", 1.0),
    ("PF00975", "Thioesterase", 0.4),
    ("PF13193", "AMP-binding_C", 0.2),
    ("PF00195", "Chal_sti_synt_N", 0.8),
    ("PF02797", "Chal_sti_synt_C", 0.8),
    ("PF03936", "Terpene_synth_C", 1.0),
    ("PF01397", "Terpene_synth", 1.0),
    ("PF13243", "SQHop_cyclase_C", 0.6),
    ("PF13249", "SQHop_cyclase_N", 0.6),
    ("PF00494", "SQS_PSY", 0.9),
    ("PF05147", "LANC_like", 1.0),
    ("PF04738", "Lant_dehyd", 1.0),
    ("PF04055", "Radical_SAM", 0.3),
    ("PF13186", "SPASM", 0.4),
    ("PF02624", "YcaO", 1.0),
    ("PF04183", "IucA_IucC", 1.0),
    ("PF06276", "FhuF", 0.4),
    ("PF00067", "p450", 0.35),
    ("PF01494", "FAD_binding_3", 0.35),
    ("PF08241", "Methyltransf_11", 0.25),
    ("PF00535", "Glycos_transf_2", 0.2),
    ("PF00155", "Aminotran_1_2", 0.2),
    ("PF01370", "Epimerase", 0.2),
    ("PF13489", "Methyltransf_23", 0.25),
    ("PF00201", "UDPGT", 0.3),
    ("PF01636", "APH", 0.1),
];

/// Additional biosynthetic marker families used by rule set V3 only. They were shortlisted by
/// enrichment inside curated MIBiG loci on the *development* genomes (never on the reserve set)
/// and then checked for biological sense: NRPS terminal reductase and MbtH partner protein,
/// SARP pathway regulators, fungal PKS starter-unit transacylase, bacterial terpene synthase,
/// siderophore-pathway enzymes, melanin (tyrosinase/MelC1) and deoxysugar enzymes.
pub const BGC_FAMILIES_EXTRA: &[(&str, &str, f32)] = &[
    ("PF19086", "Terpene_syn_C_2", 1.0),
    ("PF07993", "NAD_binding_4", 0.5),
    ("PF03621", "MbtH", 0.8),
    ("PF03704", "BTAD", 0.5),
    ("PF16073", "SAT", 0.8),
    ("PF13434", "Lys_Orn_oxgnase", 0.6),
    ("PF00425", "Chorismate_bind", 0.4),
    ("PF00264", "Tyrosinase", 0.6),
    ("PF06236", "MelC1", 0.8),
    ("PF08541", "ACP_syn_III_C", 0.3),
    ("PF08545", "ACP_syn_III", 0.3),
    ("PF16363", "GDP_Man_Dehyd", 0.4),
    ("PF00108", "Thiolase_N", 0.3),
    ("PF00891", "Methyltransf_2", 0.3),
];

/// Revision-4 families, shortlisted on the development genomes' still-missed loci (RiPP precursors,
/// lasso-peptide and methanobactin-type modification enzymes, bacteriocin processing peptidase,
/// PQQ-D recognition-element protein, carotenoid cyclase, taurine-type dioxygenase).
pub const BGC_FAMILIES_EXTRA2: &[(&str, &str, f32)] = &[
    ("PF19402", "RamS", 1.0),
    ("PF13471", "Transglut_core3", 1.0),
    ("PF09221", "Bacteriocin_IId", 1.0),
    ("PF05114", "MbnB_TglH_ChrH", 1.0),
    ("PF03412", "Peptidase_C39", 0.5),
    ("PF05402", "PqqD", 0.4),
    ("PF05834", "Lycopene_cycl", 0.6),
    ("PF02668", "TauD", 0.3),
];

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DomainFamilyModel {
    pub accession: String,
    pub name: String,
    pub ga_bits: f32,
    pub seed_sequences: usize,
    pub model_length: usize,
    pub accessory_weight: f32,
    /// `seed-built` (native Henikoff builder on the Pfam seed) or `official-pfam` (released HMM).
    #[serde(default)]
    pub model_source: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DomainHit {
    pub accession: String,
    pub bitscore: f32,
    pub i_evalue: f64,
    /// Domain p-value (`i_evalue / n_targets`), comparable to GECCO's `-p` filter.
    #[serde(default)]
    pub pvalue: f64,
    pub ali_from: usize,
    pub ali_to: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BgcGene {
    pub seq_id: String,
    pub gene_index: usize,
    pub beg: usize,
    pub end: usize,
    pub strand: char,
    pub length_aa: usize,
    pub hits: Vec<DomainHit>,
    pub core_class: Option<String>,
    pub support_weight: f32,
    pub smoothed_support: f32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BgcCluster {
    pub cluster_id: String,
    pub seq_id: String,
    pub start: usize,
    pub end: usize,
    pub gene_count: usize,
    pub core_gene_count: usize,
    pub classes: Vec<String>,
    pub domains: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct NativeGenomeScan {
    pub seq_id: String,
    pub genome_bp: usize,
    pub genes_called: usize,
    pub genes: Vec<BgcGene>,
    pub clusters: Vec<BgcCluster>,
    pub clusters_v1: Vec<BgcCluster>,
    pub gene_call_ms: u128,
    pub domain_scan_ms: u128,
    pub cluster_call_ms: u128,
}

// ---------------------------------------------------------------------------------------------
// Stage 1: Pfam seed parsing & Plan-7 profile construction
// ---------------------------------------------------------------------------------------------

/// Parse a Stockholm 1.0 alignment (interleaved blocks supported) into `(name, aligned)` rows,
/// converting `.` gap characters to `-`, plus the `GA` gathering cut-off from `#=GF GA`.
pub fn parse_stockholm(text: &str) -> (Vec<(String, String)>, f32) {
    let mut order: Vec<String> = Vec::new();
    let mut rows: HashMap<String, String> = HashMap::new();
    let mut ga = 20.0f32;
    for line in text.lines() {
        if line.starts_with("#=GF GA") {
            if let Some(v) = line.split_whitespace().nth(2) {
                ga = v.trim_end_matches(';').parse().unwrap_or(20.0);
            }
            continue;
        }
        if line.is_empty() || line.starts_with('#') || line.starts_with("//") {
            continue;
        }
        let mut it = line.split_whitespace();
        let (Some(name), Some(seq)) = (it.next(), it.next()) else {
            continue;
        };
        if !rows.contains_key(name) {
            order.push(name.to_string());
        }
        rows.entry(name.to_string())
            .or_default()
            .push_str(&seq.replace('.', "-"));
    }
    (
        order
            .into_iter()
            .map(|n| {
                let s = rows.remove(&n).unwrap_or_default();
                (n, s)
            })
            .collect(),
        ga,
    )
}

/// Where the Plan-7 profiles come from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModelSource {
    /// Built natively from the Pfam seed alignment (`tools/gecco/pfam_seed/<acc>.sto`).
    SeedBuilt,
    /// Released Pfam HMMs exported to HMMER3/f text (`tools/gecco/pfam_official/<acc>.hmm`);
    /// any family missing from that directory falls back to the seed-built model.
    OfficialPfam,
}

/// Build all 34 native Plan-7 models from `tools/gecco/pfam_seed/<acc>.sto`.
pub fn build_domain_models(
    seed_dir: &Path,
) -> Result<Vec<(DomainFamilyModel, P7Hmm)>, String> {
    build_domain_models_from(seed_dir, None)
}

/// Build the 34 models, optionally preferring released Pfam HMMs from `official_dir`.
pub fn build_domain_models_from(
    seed_dir: &Path,
    official_dir: Option<&Path>,
) -> Result<Vec<(DomainFamilyModel, P7Hmm)>, String> {
    build_domain_models_ext(seed_dir, official_dir, false)
}

/// Build the 34 published families, plus the 14 V3 families when `include_extra` is set.
pub fn build_domain_models_ext(
    seed_dir: &Path,
    official_dir: Option<&Path>,
    include_extra: bool,
) -> Result<Vec<(DomainFamilyModel, P7Hmm)>, String> {
    build_domain_models_level(seed_dir, official_dir, if include_extra { 1 } else { 0 })
}

/// Build models at a family level: 0 = 34 published, 1 = + 14 revision-3 families,
/// 2 = + 8 revision-4 families.
pub fn build_domain_models_level(
    seed_dir: &Path,
    official_dir: Option<&Path>,
    level: u8,
) -> Result<Vec<(DomainFamilyModel, P7Hmm)>, String> {
    let mut out = Vec::new();
    let extra1: &[(&str, &str, f32)] = if level >= 1 { BGC_FAMILIES_EXTRA } else { &[] };
    let extra2: &[(&str, &str, f32)] = if level >= 2 { BGC_FAMILIES_EXTRA2 } else { &[] };
    for (acc, name, w) in BGC_FAMILIES.iter().chain(extra1.iter()).chain(extra2.iter()) {
        let path = seed_dir.join(format!("{}.sto", acc));
        let txt = fs::read_to_string(&path)
            .map_err(|e| format!("Missing Pfam seed {}: {}", path.display(), e))?;
        let (rows, ga) = parse_stockholm(&txt);
        let alen = rows.first().map(|r| r.1.len()).unwrap_or(0);
        let rows: Vec<(String, String)> = rows.into_iter().filter(|r| r.1.len() == alen).collect();
        let msa: Vec<(&str, &str)> = rows.iter().map(|(a, b)| (a.as_str(), b.as_str())).collect();

        let official = official_dir
            .map(|d| d.join(format!("{}.hmm", acc)))
            .filter(|p| p.exists())
            .and_then(|p| fs::read_to_string(p).ok())
            .and_then(|t| P7Hmm::parse_hmmer3_ascii(&t).ok());
        let (hmm, source) = match official {
            Some(h) => (h, "official-pfam"),
            None => (build_native_plan7_from_msa(name, acc, name, &msa, 0.5)?, "seed-built"),
        };
        out.push((
            DomainFamilyModel {
                accession: acc.to_string(),
                name: name.to_string(),
                ga_bits: ga,
                seed_sequences: msa.len(),
                model_length: hmm.m,
                accessory_weight: *w,
                model_source: source.to_string(),
            },
            hmm,
        ));
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------
// Stage 2-3: gene calling and domain scanning
// ---------------------------------------------------------------------------------------------

pub fn read_fasta_records(path: &Path) -> Result<Vec<(String, String)>, String> {
    let txt = fs::read_to_string(path).map_err(|e| format!("{}: {}", path.display(), e))?;
    let mut recs = Vec::new();
    let mut id = String::new();
    let mut seq = String::new();
    for line in txt.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if !id.is_empty() {
                recs.push((id.clone(), seq.to_ascii_uppercase()));
            }
            id = h.split_whitespace().next().unwrap_or("").to_string();
            seq = String::new();
        } else {
            seq.push_str(line.trim());
        }
    }
    if !id.is_empty() {
        recs.push((id, seq.to_ascii_uppercase()));
    }
    Ok(recs)
}

fn scan_gene_domains(
    aa: &str,
    models: &[(DomainFamilyModel, P7Hmm)],
    prefilters: &[P7Profile],
    n_targets: f64,
) -> Vec<DomainHit> {
    let dsq = digitize_protein_seq(aa);
    if dsq.len() < 31 {
        return Vec::new();
    }
    let l = dsq.len() - 1;
    let nullsc = p7_bg_null_one(l);
    let mut hits = Vec::new();
    for (fi, (fam, hmm)) in models.iter().enumerate() {
        let usc = p7_gmsv(&dsq, l, &prefilters[fi]);
        if (usc - nullsc) / LN2 < MSV_PREFILTER_BITS {
            continue;
        }
        let res = score_sequence_plan7_native(hmm, "g", "", aa, n_targets, n_targets);
        if let Some(best) = res
            .domains
            .iter()
            .max_by(|a, b| a.bitscore.partial_cmp(&b.bitscore).unwrap())
        {
            // Keep every candidate down to CANDIDATE_GA_FRACTION of GA (or a GECCO-style p-value
            // pass); the active cut-off policy is applied afterwards by `filter_hits`.
            let pvalue = best.i_evalue / n_targets.max(1.0);
            if best.bitscore >= fam.ga_bits * CANDIDATE_GA_FRACTION || pvalue <= 1e-9 {
                hits.push(DomainHit {
                    accession: fam.accession.clone(),
                    bitscore: best.bitscore,
                    i_evalue: best.i_evalue,
                    pvalue,
                    ali_from: best.iali,
                    ali_to: best.jali,
                });
            }
        }
    }
    hits
}

/// Fraction of the Pfam GA cut-off above which sub-threshold candidate hits are retained.
const CANDIDATE_GA_FRACTION: f32 = 0.6;

/// Domain-hit acceptance policy. `Ga` is the policy used for the published benchmark (#59).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CutoffPolicy {
    /// bit-score >= Pfam curated gathering threshold.
    Ga,
    /// bit-score >= 85% of GA.
    Ga085,
    /// bit-score >= 70% of GA.
    Ga070,
    /// domain p-value <= 1e-9 (GECCO's default `-p` filter).
    P1e9,
    /// GA pass **or** p-value <= 1e-9.
    GaOrP1e9,
}

impl CutoffPolicy {
    pub fn label(&self) -> &'static str {
        match self {
            CutoffPolicy::Ga => "GA",
            CutoffPolicy::Ga085 => "0.85xGA",
            CutoffPolicy::Ga070 => "0.70xGA",
            CutoffPolicy::P1e9 => "p<=1e-9",
            CutoffPolicy::GaOrP1e9 => "GA or p<=1e-9",
        }
    }
    pub fn all() -> [CutoffPolicy; 5] {
        [
            CutoffPolicy::Ga,
            CutoffPolicy::Ga085,
            CutoffPolicy::Ga070,
            CutoffPolicy::P1e9,
            CutoffPolicy::GaOrP1e9,
        ]
    }
    pub fn accepts(&self, h: &DomainHit, ga: f32) -> bool {
        match self {
            CutoffPolicy::Ga => h.bitscore >= ga,
            CutoffPolicy::Ga085 => h.bitscore >= ga * 0.85,
            CutoffPolicy::Ga070 => h.bitscore >= ga * 0.70,
            CutoffPolicy::P1e9 => h.pvalue <= 1e-9,
            CutoffPolicy::GaOrP1e9 => h.bitscore >= ga || h.pvalue <= 1e-9,
        }
    }
}

/// Return a copy of `genes` with only the hits accepted by `policy`.
pub fn filter_hits(
    genes: &[BgcGene],
    ga: &HashMap<String, f32>,
    policy: CutoffPolicy,
) -> Vec<BgcGene> {
    genes
        .iter()
        .map(|g| {
            let mut c = g.clone();
            c.hits
                .retain(|h| policy.accepts(h, ga.get(&h.accession).copied().unwrap_or(f32::MAX)));
            c
        })
        .collect()
}

// ---------------------------------------------------------------------------------------------
// Stage 4: architecture-aware cluster calling
// ---------------------------------------------------------------------------------------------

fn has(g: &BgcGene, acc: &str) -> bool {
    g.hits.iter().any(|h| h.accession == acc)
}

/// Rule-set revision. `V1` is the original pre-diagnosis rule set; `V2` adds the NRPS-like
/// (stand-alone adenylation + carrier protein) rule and relaxes type III PKS to either domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RuleSet {
    V1,
    V2,
    /// V2 plus the 14 additional marker families, terminal-reductase / MbtH NRPS rules,
    /// bacterial terpene synthase, fungal starter-unit PKS and a melanin rule.
    V3,
    /// Revision 3 plus: stand-alone terpene synthases seed alone, RiPP precursor / lasso-peptide /
    /// methanobactin-type markers, stricter NRPS-like adjacency, and a seed floor for the weak
    /// single-marker RiPP calls (YcaO, radical-SAM/SPASM) that produced housekeeping false positives.
    V4,
}

/// Core biosynthetic gene definitions. Multi-domain *architecture* is required so that primary
/// metabolism look-alikes do not seed clusters:
/// - modular type I PKS needs a ketosynthase fused to an acyltransferase/reductase/dehydratase
///   (type II fatty-acid synthase genes `fabB`/`fabF` carry a lone ketosynthase),
/// - NRPS needs a condensation domain, or an adenylation domain fused to a carrier protein
///   (acyl-CoA synthetases `fadD`/`acs` lack a carrier domain).
pub fn classify_core_gene(g: &BgcGene, rules: RuleSet) -> Option<&'static str> {
    let v3 = matches!(rules, RuleSet::V3 | RuleSet::V4);
    let v4 = rules == RuleSet::V4;
    if has(g, "PF00109")
        && (has(g, "PF00698") || has(g, "PF08659") || has(g, "PF14765") || has(g, "PF16197")
            || (v3 && has(g, "PF16073")))
    {
        return Some("PKS-I");
    }
    if has(g, "PF00668")
        || (has(g, "PF00501") && has(g, "PF00550"))
        || (v3 && has(g, "PF00501") && has(g, "PF07993"))
    {
        return Some("NRPS");
    }
    let t3 = if rules != RuleSet::V1 {
        has(g, "PF00195") || has(g, "PF02797")
    } else {
        has(g, "PF00195") && has(g, "PF02797")
    };
    if t3 {
        return Some("PKS-III");
    }
    if has(g, "PF03936")
        || has(g, "PF01397")
        || (v3 && has(g, "PF19086"))
    {
        return Some("Terpene");
    }
    if has(g, "PF13243") && has(g, "PF13249") {
        // Squalene-hopene-type cyclase pair: kept as a multi-gene terpene core before revision 4.
        return Some(if v4 { "Terpene-cyclase" } else { "Terpene" });
    }
    if v4
        && (has(g, "PF19402") || has(g, "PF13471") || has(g, "PF09221") || has(g, "PF05114"))
    {
        // Precursor-peptide / lasso-peptide / methanobactin-type RiPP marker genes are specific.
        return Some("RiPP");
    }
    if has(g, "PF05147") || has(g, "PF04738") {
        return Some("RiPP");
    }
    if has(g, "PF02624") || (has(g, "PF04055") && has(g, "PF13186")) {
        // YcaO and radical-SAM/SPASM also occur in primary metabolism; from revision 4 they are
        // treated as weak RiPP evidence that needs neighbouring support.
        return Some(if v4 { "RiPP-weak" } else { "RiPP" });
    }
    if has(g, "PF04183") {
        return Some("Siderophore");
    }
    if v3 && has(g, "PF06236") {
        return Some("Melanin");
    }
    None
}

/// Tunable cluster-calling parameters. Default reproduces the published (#59) behaviour exactly.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ClusterParams {
    pub support_window: usize,
    pub seed_min: f32,
    pub extend_gap: usize,
    pub extend_min: f32,
    pub pad: usize,
}

impl Default for ClusterParams {
    fn default() -> Self {
        ClusterParams {
            support_window: SUPPORT_WINDOW_GENES,
            seed_min: SEED_MIN_SUPPORT,
            extend_gap: MAX_EXTEND_GAP_BP,
            extend_min: EXTEND_MIN_SUPPORT,
            pad: FLANK_PAD_BP,
        }
    }
}

pub fn call_clusters_for_contig(
    seq_id: &str,
    genes: &mut [BgcGene],
    weights: &HashMap<String, f32>,
    rules: RuleSet,
) -> Vec<BgcCluster> {
    call_clusters_with_params(seq_id, genes, weights, rules, &ClusterParams::default())
}
pub fn call_clusters_with_params(
    seq_id: &str,
    genes: &mut [BgcGene],
    weights: &HashMap<String, f32>,
    rules: RuleSet,
    prm: &ClusterParams,
) -> Vec<BgcCluster> {
    let n = genes.len();
    if n == 0 {
        return Vec::new();
    }
    // Per-gene support weight (capped at 1.0).
    for g in genes.iter_mut() {
        g.core_class = classify_core_gene(g, rules).map(|s| s.to_string());
        let mut w: f32 = g
            .hits
            .iter()
            .map(|h| weights.get(&h.accession).copied().unwrap_or(0.0))
            .sum();
        if g.core_class.is_some() {
            w = w.max(1.0);
        }
        g.support_weight = w.min(1.0);
    }
    // Type II PKS: two ketosynthase-bearing genes within a few genes of each other.
    let ks_lone: Vec<usize> = (0..n)
        .filter(|&i| has(&genes[i], "PF00109") && genes[i].core_class.is_none())
        .collect();
    let mut t2pks: HashSet<usize> = HashSet::new();
    for (a, &i) in ks_lone.iter().enumerate() {
        for &j in ks_lone.iter().skip(a + 1) {
            if j - i <= 3 {
                t2pks.insert(i);
                t2pks.insert(j);
            }
        }
    }
    for &i in &t2pks {
        genes[i].core_class = Some("PKS-II".to_string());
        genes[i].support_weight = 1.0;
    }
    // V2: NRPS-like clusters -- a stand-alone adenylation gene with a carrier-protein gene within
    // three genes (e.g. enterobactin `entE` + `entB`), where neither carries a core architecture.
    if rules != RuleSet::V1 {
        // Revision 4 requires the partner gene to be adjacent (operon-like packing); earlier
        // revisions allowed up to three genes either side.
        let partner_reach: usize = if rules == RuleSet::V4 { 1 } else { 3 };
        let a_only: Vec<usize> = (0..n)
            .filter(|&i| has(&genes[i], "PF00501") && genes[i].core_class.is_none())
            .collect();
        for i in a_only {
            let lo = i.saturating_sub(partner_reach);
            let hi = (i + partner_reach).min(n - 1);
            let pcp = (lo..=hi).find(|&j| {
                j != i
                    && (has(&genes[j], "PF00550")
                        || (matches!(rules, RuleSet::V3 | RuleSet::V4) && has(&genes[j], "PF03621")))
            });
            if let Some(j) = pcp {
                genes[i].core_class = Some("NRPS-like".to_string());
                genes[i].support_weight = 1.0;
                if genes[j].core_class.is_none() {
                    genes[j].support_weight = genes[j].support_weight.max(0.6);
                }
            }
        }
    }
    // Windowed smoothed support.
    let sw: Vec<f32> = genes.iter().map(|g| g.support_weight).collect();
    for i in 0..n {
        let lo = i.saturating_sub(prm.support_window);
        let hi = (i + prm.support_window).min(n - 1);
        genes[i].smoothed_support = sw[lo..=hi].iter().sum();
    }
    // Seed & extend.
    let mut regions: Vec<(usize, usize)> = Vec::new();
    for i in 0..n {
        let seed_min = match genes[i].core_class.as_deref() {
            Some("NRPS-like") => 0.99,
            Some("Melanin") => 1.5f32.min(prm.seed_min),
            // Revision 4: stand-alone terpene synthases and precursor-type RiPP markers are
            // specific enough to seed on their own; weak RiPP markers need real neighbour support.
            Some("Terpene") if rules == RuleSet::V4 => 0.99,
            Some("RiPP") if rules == RuleSet::V4 => 0.99,
            Some("RiPP-weak") => 1.8f32.max(prm.seed_min),
            _ => prm.seed_min,
        };
        if genes[i].core_class.is_none() || genes[i].smoothed_support < seed_min {
            continue;
        }
        let (mut lo, mut hi) = (i, i);
        while lo > 0 {
            let p = &genes[lo - 1];
            let gap = genes[lo].beg.saturating_sub(p.end);
            if gap <= prm.extend_gap
                && (p.support_weight > 0.0 || p.smoothed_support >= prm.extend_min)
            {
                lo -= 1;
            } else {
                break;
            }
        }
        while hi + 1 < n {
            let q = &genes[hi + 1];
            let gap = q.beg.saturating_sub(genes[hi].end);
            if gap <= prm.extend_gap
                && (q.support_weight > 0.0 || q.smoothed_support >= prm.extend_min)
            {
                hi += 1;
            } else {
                break;
            }
        }
        regions.push((lo, hi));
    }
    regions.sort();
    let mut merged: Vec<(usize, usize)> = Vec::new();
    for (lo, hi) in regions {
        if let Some(last) = merged.last_mut() {
            let gap = genes[lo].beg.saturating_sub(genes[last.1].end);
            if lo <= last.1 + 1 || gap <= prm.extend_gap {
                last.1 = last.1.max(hi);
                continue;
            }
        }
        merged.push((lo, hi));
    }
    merged
        .into_iter()
        .enumerate()
        .map(|(k, (lo, hi))| {
            let slice = &genes[lo..=hi];
            let classes: BTreeSet<String> =
                slice.iter().filter_map(|g| g.core_class.clone()).collect();
            let domains: BTreeSet<String> = slice
                .iter()
                .flat_map(|g| g.hits.iter().map(|h| h.accession.clone()))
                .collect();
            BgcCluster {
                cluster_id: format!("{}_native_{}", seq_id, k + 1),
                seq_id: seq_id.to_string(),
                start: genes[lo].beg.saturating_sub(prm.pad).max(1),
                end: genes[hi].end + prm.pad,
                gene_count: slice.len(),
                core_gene_count: slice.iter().filter(|g| g.core_class.is_some()).count(),
                classes: classes.into_iter().collect(),
                domains: domains.into_iter().collect(),
            }
        })
        .collect()
}

/// Run the complete native pipeline on one FASTA genome (multi-record supported).
pub fn run_native_bgc_scan(
    fasta: &Path,
    models: &[(DomainFamilyModel, P7Hmm)],
) -> Result<NativeGenomeScan, String> {
    let recs = read_fasta_records(fasta)?;
    let prefilters: Vec<P7Profile> = models
        .iter()
        .map(|(_, h)| P7Profile::configure(h, 400, true))
        .collect();
    let weights: HashMap<String, f32> = models
        .iter()
        .map(|(f, _)| (f.accession.clone(), f.accessory_weight))
        .collect();
    let ga_map: HashMap<String, f32> = models
        .iter()
        .map(|(f, _)| (f.accession.clone(), f.ga_bits))
        .collect();

    let mut all_genes: Vec<BgcGene> = Vec::new();
    let mut clusters: Vec<BgcCluster> = Vec::new();
    let mut clusters_v1: Vec<BgcCluster> = Vec::new();
    let (mut gene_ms, mut dom_ms, mut cl_ms) = (0u128, 0u128, 0u128);
    let mut genome_bp = 0usize;
    let mut genes_called = 0usize;
    let seq_label = recs.first().map(|r| r.0.clone()).unwrap_or_default();

    for (id, dna) in &recs {
        genome_bp += dna.len();
        let t0 = Instant::now();
        let out = run_native_prokaryotic_gene_caller(id, dna, false);
        gene_ms += t0.elapsed().as_millis();
        genes_called += out.predicted_genes.len();

        let t1 = Instant::now();
        let n_targets = out.predicted_genes.len().max(1) as f64;
        let mut genes: Vec<BgcGene> = out
            .predicted_genes
            .par_iter()
            .enumerate()
            .map(|(i, pg)| {
                let aa = translate_nt_to_aa(dna, pg.beg, pg.end, pg.strand);
                let aa: String = aa.trim_end_matches('*').to_string();
                let hits = scan_gene_domains(&aa, models, &prefilters, n_targets);
                BgcGene {
                    seq_id: id.clone(),
                    gene_index: i,
                    beg: pg.beg,
                    end: pg.end,
                    strand: pg.strand,
                    length_aa: aa.len(),
                    hits,
                    core_class: None,
                    support_weight: 0.0,
                    smoothed_support: 0.0,
                }
            })
            .collect();
        dom_ms += t1.elapsed().as_millis();
        genes.sort_by_key(|g| g.beg);
        // Published-benchmark policy: curated Pfam gathering cut-off.
        let mut genes = filter_hits(&genes, &ga_map, CutoffPolicy::Ga);

        let t2 = Instant::now();
        let mut genes_v1 = genes.clone();
        let mut cl1 = call_clusters_for_contig(id, &mut genes_v1, &weights, RuleSet::V1);
        let mut cl = call_clusters_for_contig(id, &mut genes, &weights, RuleSet::V2);
        cl_ms += t2.elapsed().as_millis();
        clusters_v1.append(&mut cl1);
        clusters.append(&mut cl);
        all_genes.extend(genes);
    }
    Ok(NativeGenomeScan {
        seq_id: seq_label,
        genome_bp,
        genes_called,
        genes: all_genes,
        clusters,
        clusters_v1,
        gene_call_ms: gene_ms,
        domain_scan_ms: dom_ms,
        cluster_call_ms: cl_ms,
    })
}

// ---------------------------------------------------------------------------------------------
// Stage 5: head-to-head calibration vs GECCO and MIBiG
// ---------------------------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Interval {
    pub seq_id: String,
    pub start: usize,
    pub end: usize,
    pub label: String,
}

fn overlap_bp(a_s: usize, a_e: usize, b_s: usize, b_e: usize) -> usize {
    let s = a_s.max(b_s);
    let e = a_e.min(b_e);
    if e >= s {
        e - s + 1
    } else {
        0
    }
}

/// Fraction of `target` base pairs covered by the union of `preds` (same contig).
fn covered_fraction(target: &Interval, preds: &[Interval]) -> f64 {
    let mut segs: Vec<(usize, usize)> = preds
        .iter()
        .filter(|p| p.seq_id == target.seq_id)
        .filter_map(|p| {
            let o = overlap_bp(target.start, target.end, p.start, p.end);
            if o > 0 {
                Some((p.start.max(target.start), p.end.min(target.end)))
            } else {
                None
            }
        })
        .collect();
    segs.sort();
    let mut total = 0usize;
    let mut cur: Option<(usize, usize)> = None;
    for (s, e) in segs {
        match cur {
            Some((cs, ce)) if s <= ce + 1 => cur = Some((cs, ce.max(e))),
            Some((cs, ce)) => {
                total += ce - cs + 1;
                cur = Some((s, e));
            }
            None => cur = Some((s, e)),
        }
    }
    if let Some((cs, ce)) = cur {
        total += ce - cs + 1;
    }
    total as f64 / (target.end - target.start + 1).max(1) as f64
}

fn union_coverage_bp(preds: &[Interval]) -> BTreeMap<String, Vec<(usize, usize)>> {
    let mut by: BTreeMap<String, Vec<(usize, usize)>> = BTreeMap::new();
    for p in preds {
        by.entry(p.seq_id.clone()).or_default().push((p.start, p.end));
    }
    for v in by.values_mut() {
        v.sort();
        let mut m: Vec<(usize, usize)> = Vec::new();
        for (s, e) in v.drain(..) {
            if let Some(last) = m.last_mut() {
                if s <= last.1 + 1 {
                    last.1 = last.1.max(e);
                    continue;
                }
            }
            m.push((s, e));
        }
        *v = m;
    }
    by
}

fn bp_jaccard(a: &[Interval], b: &[Interval]) -> f64 {
    let ua = union_coverage_bp(a);
    let ub = union_coverage_bp(b);
    let (mut inter, mut uni) = (0usize, 0usize);
    let keys: BTreeSet<&String> = ua.keys().chain(ub.keys()).collect();
    for k in keys {
        let va = ua.get(k).cloned().unwrap_or_default();
        let vb = ub.get(k).cloned().unwrap_or_default();
        let la: usize = va.iter().map(|(s, e)| e - s + 1).sum();
        let lb: usize = vb.iter().map(|(s, e)| e - s + 1).sum();
        let mut i_bp = 0usize;
        for (s1, e1) in &va {
            for (s2, e2) in &vb {
                i_bp += overlap_bp(*s1, *e1, *s2, *e2);
            }
        }
        inter += i_bp;
        uni += la + lb - i_bp;
    }
    if uni == 0 {
        1.0
    } else {
        inter as f64 / uni as f64
    }
}

/// Cluster-level matches: pairs whose overlap covers >= 50% of the shorter cluster.
fn cluster_matches(a: &[Interval], b: &[Interval]) -> usize {
    let mut used_b: HashSet<usize> = HashSet::new();
    let mut matched = 0usize;
    for ia in a {
        for (j, ib) in b.iter().enumerate() {
            if used_b.contains(&j) || ia.seq_id != ib.seq_id {
                continue;
            }
            let o = overlap_bp(ia.start, ia.end, ib.start, ib.end);
            let shorter = (ia.end - ia.start + 1).min(ib.end - ib.start + 1);
            if o * 2 >= shorter && o > 0 {
                used_b.insert(j);
                matched += 1;
                break;
            }
        }
    }
    matched
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToolScore {
    pub tool: String,
    pub n_clusters: usize,
    pub total_bp: usize,
    pub mibig_recall_any_overlap: f64,
    pub mibig_recall_half_covered: f64,
    pub mibig_truth_hit: usize,
    pub predicted_overlapping_mibig_frac: f64,
    pub wall_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TruthEntry {
    pub mibig_id: String,
    pub seq_id: String,
    pub start: usize,
    pub end: usize,
    pub classes: String,
    pub quality: String,
    pub native_cover_frac: f64,
    pub gecco_cover_frac: f64,
    pub gecco_sensitive_cover_frac: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DomainConcordance {
    pub native_calls: usize,
    pub gecco_calls: usize,
    pub shared: usize,
    pub jaccard: f64,
    pub family_count_pearson_r: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BgcGenomeComparison {
    pub genome_id: String,
    pub organism: String,
    pub kingdom: String,
    pub regime_note: String,
    pub genome_bp: usize,
    pub genes_called_native: usize,
    pub native: ToolScore,
    pub native_rules_v1: ToolScore,
    pub gecco_default: ToolScore,
    pub gecco_sensitive: ToolScore,
    pub mibig_truth_loci: usize,
    pub native_vs_gecco_bp_jaccard: f64,
    pub native_vs_gecco_matched_clusters: usize,
    pub native_vs_gecco_cluster_jaccard: f64,
    pub domain_concordance: DomainConcordance,
    pub native_gene_call_ms: u128,
    pub native_domain_scan_ms: u128,
    pub native_cluster_call_ms: u128,
    pub truth: Vec<TruthEntry>,
    pub native_clusters: Vec<BgcCluster>,
    pub gecco_default_clusters: Vec<Interval>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BgcCalibrationReport {
    pub pipeline_name: String,
    pub external_tool_version: String,
    pub truth_database: String,
    pub domain_models: Vec<DomainFamilyModel>,
    pub genomes: Vec<BgcGenomeComparison>,
    pub total_genome_bp: usize,
    pub total_mibig_loci: usize,
    pub native_total_recalled_any: usize,
    pub native_v1_total_recalled_any: usize,
    pub gecco_total_recalled_any: usize,
    pub gecco_sensitive_total_recalled_any: usize,
    pub native_total_recalled_half: usize,
    pub gecco_total_recalled_half: usize,
    pub mean_bp_jaccard_native_vs_gecco: f64,
    pub mean_domain_jaccard: f64,
    pub bgc_rich_clusters_native: usize,
    pub bgc_rich_clusters_gecco: usize,
    pub negative_control_clusters_native: usize,
    pub negative_control_clusters_gecco: usize,
    pub signal_rating_beakers: u8,
    pub signal_rating_label: String,
    pub ascii_master_summary_table: String,
    pub total_duration_ms: u128,
}

fn load_mibig_truth(mibig_dir: &Path, accession_base: &str) -> Vec<(String, usize, usize, String, String)> {
    let mut out: Vec<(String, usize, usize, String, String)> = Vec::new();
    let Ok(rd) = fs::read_dir(mibig_dir) else {
        return out;
    };
    for ent in rd.flatten() {
        let p = ent.path();
        if p.extension().and_then(|s| s.to_str()) != Some("json") {
            continue;
        }
        let Ok(txt) = fs::read_to_string(&p) else { continue };
        let Ok(v) = serde_json::from_str::<serde_json::Value>(&txt) else { continue };
        let id = v["accession"].as_str().unwrap_or("").to_string();
        let quality = v["quality"].as_str().unwrap_or("").to_string();
        let classes: Vec<String> = v["biosynthesis"]["classes"]
            .as_array()
            .map(|a| a.iter().filter_map(|c| c["class"].as_str().map(String::from)).collect())
            .unwrap_or_default();
        if let Some(loci) = v["loci"].as_array() {
            for l in loci {
                let acc = l["accession"].as_str().unwrap_or("");
                if acc.split('.').next() != Some(accession_base) {
                    continue;
                }
                let s = l["location"]["from"].as_u64().unwrap_or(0) as usize;
                let e = l["location"]["to"].as_u64().unwrap_or(0) as usize;
                if e > s {
                    out.push((id.clone(), s.max(1), e, classes.join("/"), quality.clone()));
                }
            }
        }
    }
    out.sort_by_key(|t| (t.1, t.2));
    // Collapse records sharing identical coordinates (several MIBiG entries can describe the same locus).
    out.dedup_by(|a, b| a.1 == b.1 && a.2 == b.2);
    out
}

fn parse_gecco_clusters(path: &Path) -> Vec<Interval> {
    let Ok(txt) = fs::read_to_string(path) else { return Vec::new() };
    let mut out = Vec::new();
    for (i, line) in txt.lines().enumerate() {
        if i == 0 || line.is_empty() {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 7 {
            continue;
        }
        out.push(Interval {
            seq_id: c[0].to_string(),
            start: c[2].parse().unwrap_or(0),
            end: c[3].parse().unwrap_or(0),
            label: c[6].to_string(),
        });
    }
    out
}

fn gene_key(seq: &str, strand: char, beg: usize, end: usize) -> String {
    // Stop-codon coordinate is shared between callers even when start codons differ.
    if strand == '+' {
        format!("{}:+:{}", seq, end)
    } else {
        format!("{}:-:{}", seq, beg)
    }
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    if n < 2.0 {
        return 1.0;
    }
    let ma = a.iter().sum::<f64>() / n;
    let mb = b.iter().sum::<f64>() / n;
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for i in 0..a.len() {
        sab += (a[i] - ma) * (b[i] - mb);
        saa += (a[i] - ma).powi(2);
        sbb += (b[i] - mb).powi(2);
    }
    if saa == 0.0 || sbb == 0.0 {
        0.0
    } else {
        sab / (saa * sbb).sqrt()
    }
}

fn domain_concordance(
    native: &NativeGenomeScan,
    gecco_features: &Path,
    families: &HashSet<String>,
) -> DomainConcordance {
    let mut nset: HashSet<(String, String)> = HashSet::new();
    for g in &native.genes {
        for h in &g.hits {
            nset.insert((gene_key(&g.seq_id, g.strand, g.beg, g.end), h.accession.clone()));
        }
    }
    let mut gset: HashSet<(String, String)> = HashSet::new();
    if let Ok(txt) = fs::read_to_string(gecco_features) {
        for (i, line) in txt.lines().enumerate() {
            if i == 0 {
                continue;
            }
            let c: Vec<&str> = line.split('\t').collect();
            if c.len() < 8 || !families.contains(c[5]) {
                continue;
            }
            let strand = c[4].chars().next().unwrap_or('+');
            let beg: usize = c[2].parse().unwrap_or(0);
            let end: usize = c[3].parse().unwrap_or(0);
            gset.insert((gene_key(c[0], strand, beg, end), c[5].to_string()));
        }
    }
    let shared = nset.intersection(&gset).count();
    let union = nset.union(&gset).count();
    let mut fam_n: BTreeMap<&str, f64> = BTreeMap::new();
    let mut fam_g: BTreeMap<&str, f64> = BTreeMap::new();
    for f in families {
        fam_n.insert(f.as_str(), 0.0);
        fam_g.insert(f.as_str(), 0.0);
    }
    for (_, a) in &nset {
        if let Some(v) = fam_n.get_mut(a.as_str()) {
            *v += 1.0;
        }
    }
    for (_, a) in &gset {
        if let Some(v) = fam_g.get_mut(a.as_str()) {
            *v += 1.0;
        }
    }
    let va: Vec<f64> = fam_n.values().copied().collect();
    let vb: Vec<f64> = fam_g.values().copied().collect();
    DomainConcordance {
        native_calls: nset.len(),
        gecco_calls: gset.len(),
        shared,
        jaccard: if union == 0 { 1.0 } else { shared as f64 / union as f64 },
        family_count_pearson_r: pearson(&va, &vb),
    }
}

fn run_gecco(
    gecco_exe: &Path,
    fasta: &Path,
    out_dir: &Path,
    threshold: f64,
) -> Result<u128, String> {
    let t0 = Instant::now();
    let st = Command::new(gecco_exe)
        .env("PYTHONIOENCODING", "utf-8")
        .env("PYTHONUTF8", "1")
        .args([
            "run",
            "-g",
            &fasta.to_string_lossy(),
            "-o",
            &out_dir.to_string_lossy(),
            "-j",
            "0",
            "--force-tsv",
            "-m",
            &format!("{}", threshold),
            "--no-color",
            "--no-markup",
            "-q",
        ])
        .output()
        .map_err(|e| format!("Failed to launch GECCO: {}", e))?;
    if !st.status.success() {
        return Err(format!("GECCO failed: {}", String::from_utf8_lossy(&st.stderr)));
    }
    Ok(t0.elapsed().as_millis())
}

fn score_tool(
    tool: &str,
    preds: &[Interval],
    truth: &[(String, usize, usize, String, String)],
    seq_id: &str,
    wall_ms: u128,
) -> ToolScore {
    let mut any = 0usize;
    let mut half = 0usize;
    for t in truth {
        let ti = Interval { seq_id: seq_id.to_string(), start: t.1, end: t.2, label: t.0.clone() };
        let f = covered_fraction(&ti, preds);
        if f > 0.0 {
            any += 1;
        }
        if f >= 0.5 {
            half += 1;
        }
    }
    let overl = preds
        .iter()
        .filter(|p| truth.iter().any(|t| overlap_bp(p.start, p.end, t.1, t.2) > 0))
        .count();
    ToolScore {
        tool: tool.to_string(),
        n_clusters: preds.len(),
        total_bp: union_coverage_bp(preds).values().flatten().map(|(s, e)| e - s + 1).sum(),
        mibig_recall_any_overlap: if truth.is_empty() { 0.0 } else { any as f64 / truth.len() as f64 },
        mibig_recall_half_covered: if truth.is_empty() { 0.0 } else { half as f64 / truth.len() as f64 },
        mibig_truth_hit: any,
        predicted_overlapping_mibig_frac: if preds.is_empty() { 0.0 } else { overl as f64 / preds.len() as f64 },
        wall_ms,
    }
}

pub fn run_bgc_calibration_benchmark(ws_root: &Path) -> Result<BgcCalibrationReport, String> {
    let t_start = Instant::now();
    let labs = if ws_root.ends_with("crvscience-labs") {
        ws_root.to_path_buf()
    } else {
        ws_root.join("crvscience-labs")
    };
    let seed_dir = labs.join("tools/gecco/pfam_seed");
    let mibig_dir = labs.join("tools/gecco/data/mibig_json_4.0");
    let gecco_exe: PathBuf = labs.join("tools/gecco/venv/Scripts/gecco.exe");
    let work = labs.join("runs/gecco_work");
    let cache = labs.join("cache/bgc");

    let models = build_domain_models(&seed_dir)?;
    let family_set: HashSet<String> = models.iter().map(|(f, _)| f.accession.clone()).collect();

    // (accession, organism, kingdom, regime)
    let specs: Vec<(&str, &str, &str, &str)> = vec![
        ("AL645882.2", "Streptomyces coelicolor A3(2)", "Prokaryota", "BGC-rich actinomycete (15 MIBiG loci)"),
        ("CP000560.1", "Bacillus velezensis FZB42", "Prokaryota", "BGC-rich Firmicute (6 MIBiG loci)"),
        ("AE004091.2", "Pseudomonas aeruginosa PAO1", "Prokaryota", "Proteobacterial pathogen (7 MIBiG loci)"),
        ("U00096.3", "Escherichia coli K-12 MG1655", "Prokaryota", "BGC-poor near-negative control (enterobactin only)"),
        ("BN001302.1", "Aspergillus nidulans FGSC A4 chr II", "Eukaryota", "Out-of-regime: spliced fungal genes, bacterial gene callers"),
    ];

    let mut genomes = Vec::new();
    for (acc, organism, kingdom, regime) in specs {
        let fasta = cache.join(format!("{}.fasta", acc));
        let base = acc.split('.').next().unwrap_or(acc);

        let native = run_native_bgc_scan(&fasta, &models)?;
        let native_wall = native.gene_call_ms + native.domain_scan_ms + native.cluster_call_ms;

        // Persist per-gene native domain calls for provenance and failure-mode auditing.
        {
            let dump_dir = labs.join("runs/bgc_work");
            let _ = fs::create_dir_all(&dump_dir);
            let mut s = String::from("seq_id\tbeg\tend\tstrand\tlength_aa\tcore_class\tsupport\tsmoothed\tdomains(acc:bits)\n");
            for g in &native.genes {
                let d: Vec<String> = g
                    .hits
                    .iter()
                    .map(|h| format!("{}:{:.1}", h.accession, h.bitscore))
                    .collect();
                s.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{}\n",
                    g.seq_id,
                    g.beg,
                    g.end,
                    g.strand,
                    g.length_aa,
                    g.core_class.clone().unwrap_or_default(),
                    g.support_weight,
                    g.smoothed_support,
                    d.join(",")
                ));
            }
            let _ = fs::write(dump_dir.join(format!("{}.native_genes.tsv", acc)), s);
        }

        let dir_def = work.join("m080").join(acc);
        // Reuse a completed GECCO run (and its recorded wall time) when BGC_REUSE_GECCO=1.
        let def_done = dir_def.join(format!("{}.clusters.tsv", acc)).exists();
        let wall_def = if std::env::var("BGC_REUSE_GECCO").is_ok() && def_done {
            fs::read_to_string(work.join("timings_m080.tsv"))
                .ok()
                .and_then(|t| {
                    t.lines()
                        .find(|l| l.starts_with(acc))
                        .and_then(|l| l.split('\t').nth(1).and_then(|v| v.trim().parse::<f64>().ok()))
                })
                .map(|s| (s * 1000.0) as u128)
                .unwrap_or(0)
        } else {
            run_gecco(&gecco_exe, &fasta, &dir_def, 0.8)?
        };
        let dir_sens = work.join(acc);
        let wall_sens = if dir_sens.join(format!("{}.clusters.tsv", acc)).exists() {
            0
        } else {
            run_gecco(&gecco_exe, &fasta, &dir_sens, 0.3)?
        };

        let gecco_def = parse_gecco_clusters(&dir_def.join(format!("{}.clusters.tsv", acc)));
        let gecco_sens = parse_gecco_clusters(&dir_sens.join(format!("{}.clusters.tsv", acc)));
        let native_iv: Vec<Interval> = native
            .clusters
            .iter()
            .map(|c| Interval {
                seq_id: c.seq_id.clone(),
                start: c.start,
                end: c.end,
                label: c.classes.join(";"),
            })
            .collect();

        let native_v1_iv: Vec<Interval> = native
            .clusters_v1
            .iter()
            .map(|c| Interval {
                seq_id: c.seq_id.clone(),
                start: c.start,
                end: c.end,
                label: c.classes.join(";"),
            })
            .collect();

        let truth = load_mibig_truth(&mibig_dir, base);
        let seq_id = native.seq_id.clone();
        let truth_entries: Vec<TruthEntry> = truth
            .iter()
            .map(|t| {
                let ti = Interval { seq_id: seq_id.clone(), start: t.1, end: t.2, label: t.0.clone() };
                TruthEntry {
                    mibig_id: t.0.clone(),
                    seq_id: seq_id.clone(),
                    start: t.1,
                    end: t.2,
                    classes: t.3.clone(),
                    quality: t.4.clone(),
                    native_cover_frac: covered_fraction(&ti, &native_iv),
                    gecco_cover_frac: covered_fraction(&ti, &gecco_def),
                    gecco_sensitive_cover_frac: covered_fraction(&ti, &gecco_sens),
                }
            })
            .collect();

        let matched = cluster_matches(&native_iv, &gecco_def);
        let cl_union = native_iv.len() + gecco_def.len() - matched;
        let dc = domain_concordance(
            &native,
            &dir_sens.join(format!("{}.features.tsv", acc)),
            &family_set,
        );

        genomes.push(BgcGenomeComparison {
            genome_id: acc.to_string(),
            organism: organism.to_string(),
            kingdom: kingdom.to_string(),
            regime_note: regime.to_string(),
            genome_bp: native.genome_bp,
            genes_called_native: native.genes_called,
            native: score_tool("native (rules v2)", &native_iv, &truth, &seq_id, native_wall),
            native_rules_v1: score_tool("native (rules v1)", &native_v1_iv, &truth, &seq_id, native_wall),
            gecco_default: score_tool("GECCO v0.11.0 (-m 0.8)", &gecco_def, &truth, &seq_id, wall_def),
            gecco_sensitive: score_tool("GECCO v0.11.0 (-m 0.3)", &gecco_sens, &truth, &seq_id, wall_sens),
            mibig_truth_loci: truth.len(),
            native_vs_gecco_bp_jaccard: bp_jaccard(&native_iv, &gecco_def),
            native_vs_gecco_matched_clusters: matched,
            native_vs_gecco_cluster_jaccard: if cl_union == 0 { 1.0 } else { matched as f64 / cl_union as f64 },
            domain_concordance: dc,
            native_gene_call_ms: native.gene_call_ms,
            native_domain_scan_ms: native.domain_scan_ms,
            native_cluster_call_ms: native.cluster_call_ms,
            truth: truth_entries,
            native_clusters: native.clusters.clone(),
            gecco_default_clusters: gecco_def.clone(),
        });
    }

    let sum = |f: &dyn Fn(&BgcGenomeComparison) -> usize| genomes.iter().map(f).sum::<usize>();
    let total_truth = sum(&|g| g.mibig_truth_loci);
    let n_any = sum(&|g| g.native.mibig_truth_hit);
    let n1_any = sum(&|g| g.native_rules_v1.mibig_truth_hit);
    let g_any = sum(&|g| g.gecco_default.mibig_truth_hit);
    let s_any = sum(&|g| g.gecco_sensitive.mibig_truth_hit);
    let half = |sel: &dyn Fn(&BgcGenomeComparison) -> &ToolScore| -> usize {
        genomes
            .iter()
            .map(|g| (sel(g).mibig_recall_half_covered * g.mibig_truth_loci as f64).round() as usize)
            .sum()
    };
    let rich: Vec<&BgcGenomeComparison> = genomes.iter().filter(|g| g.genome_id != "U00096.3" && g.kingdom == "Prokaryota").collect();
    let neg = genomes.iter().find(|g| g.genome_id == "U00096.3");
    let mean_bp = genomes.iter().map(|g| g.native_vs_gecco_bp_jaccard).sum::<f64>() / genomes.len() as f64;
    let mean_dom = genomes.iter().map(|g| g.domain_concordance.jaccard).sum::<f64>() / genomes.len() as f64;

    let mut table = String::from(
        "==================================================================================================================================================\n\
         TOOL COMPARISON #15: NATIVE ARCHITECTURE-AWARE BGC ENGINE vs. GECCO v0.11.0 (CRF) -- MIBiG 4.0 GROUND TRUTH\n\
         ==================================================================================================================================================\n\
         Genome       | Organism                          | Mb    | MIBiG | Nat v2/v1 hit | GECCO.8 hit/half | GECCO.3 hit | Nat/G.8/G.3 clusters | bp-Jacc | Dom-Jacc\n\
         -------------+-----------------------------------+-------+-------+---------------+------------------+-------------+----------------------+---------+---------\n",
    );
    for g in &genomes {
        table.push_str(&format!(
            "{:<12} | {:<33} | {:>5.2} | {:>5} | {:>5}/{:<7} | {:>8}/{:<8} | {:>11} | {:>6}/{:>5}/{:>6}  | {:>7.3} | {:>7.3}\n",
            g.genome_id,
            g.organism.chars().take(33).collect::<String>(),
            g.genome_bp as f64 / 1e6,
            g.mibig_truth_loci,
            g.native.mibig_truth_hit,
            g.native_rules_v1.mibig_truth_hit,
            g.gecco_default.mibig_truth_hit,
            (g.gecco_default.mibig_recall_half_covered * g.mibig_truth_loci as f64).round() as usize,
            g.gecco_sensitive.mibig_truth_hit,
            g.native.n_clusters,
            g.gecco_default.n_clusters,
            g.gecco_sensitive.n_clusters,
            g.native_vs_gecco_bp_jaccard,
            g.domain_concordance.jaccard,
        ));
    }
    table.push_str(&format!(
        "-------------+-----------------------------------+-------+-------+---------------+------------------+-------------+----------------------+---------+---------\n\
         TOTAL        | 5 genomes (4 bacterial, 1 fungal) | {:>5.2} | {:>5} | {:>5}/{:<7} | {:>8}/{:<8} | {:>11} |                      | {:>7.3} | {:>7.3}\n\
         ==================================================================================================================================================",
        genomes.iter().map(|g| g.genome_bp).sum::<usize>() as f64 / 1e6,
        total_truth,
        n_any,
        n1_any,
        g_any,
        half(&|g| &g.gecco_default),
        s_any,
        mean_bp,
        mean_dom,
    ));

    Ok(BgcCalibrationReport {
        pipeline_name: "bgc-compare (Native Architecture-Aware Biosynthetic Gene Cluster Engine vs. GECCO v0.11.0)".to_string(),
        external_tool_version: "GECCO v0.11.0 (pyrodigal 3.7.1, pyhmmer 0.12.3, full Pfam CRF)".to_string(),
        truth_database: "MIBiG 4.0 (mibig_json_4.0)".to_string(),
        domain_models: models.iter().map(|(f, _)| f.clone()).collect(),
        total_genome_bp: genomes.iter().map(|g| g.genome_bp).sum(),
        total_mibig_loci: total_truth,
        native_total_recalled_any: n_any,
        native_v1_total_recalled_any: n1_any,
        gecco_total_recalled_any: g_any,
        gecco_sensitive_total_recalled_any: s_any,
        native_total_recalled_half: half(&|g| &g.native),
        gecco_total_recalled_half: half(&|g| &g.gecco_default),
        mean_bp_jaccard_native_vs_gecco: mean_bp,
        mean_domain_jaccard: mean_dom,
        bgc_rich_clusters_native: rich.iter().map(|g| g.native.n_clusters).sum(),
        bgc_rich_clusters_gecco: rich.iter().map(|g| g.gecco_default.n_clusters).sum(),
        negative_control_clusters_native: neg.map(|g| g.native.n_clusters).unwrap_or(0),
        negative_control_clusters_gecco: neg.map(|g| g.gecco_default.n_clusters).unwrap_or(0),
        signal_rating_beakers: 3,
        signal_rating_label: "3 Comparative Synthesis (provisional; finalised after review)".to_string(),
        genomes,
        ascii_master_summary_table: table,
        total_duration_ms: t_start.elapsed().as_millis(),
    })
}

// ---------------------------------------------------------------------------------------------
// Calibration round 2: model-source x cut-off-policy ablation, BGC-free negative controls and a
// held-out genome set scored once with a configuration frozen on the development set.
// ---------------------------------------------------------------------------------------------

/// Genes and translated proteins for one genome, computed once and reused for every model source.
pub struct GenomeProteins {
    pub seq_id: String,
    pub genome_bp: usize,
    pub genes_called: usize,
    pub gene_call_ms: u128,
    /// `(contig, [(beg, end, strand, protein)])`
    pub contigs: Vec<(String, Vec<(usize, usize, char, String)>)>,
}

pub fn call_genome_proteins(fasta: &Path) -> Result<GenomeProteins, String> {
    let recs = read_fasta_records(fasta)?;
    let mut contigs = Vec::new();
    let (mut bp, mut n, mut ms) = (0usize, 0usize, 0u128);
    let label = recs.first().map(|r| r.0.clone()).unwrap_or_default();
    for (id, dna) in &recs {
        bp += dna.len();
        let t0 = Instant::now();
        let out = run_native_prokaryotic_gene_caller(id, dna, false);
        ms += t0.elapsed().as_millis();
        n += out.predicted_genes.len();
        let genes: Vec<(usize, usize, char, String)> = out
            .predicted_genes
            .iter()
            .map(|pg| {
                let aa = translate_nt_to_aa(dna, pg.beg, pg.end, pg.strand);
                (pg.beg, pg.end, pg.strand, aa.trim_end_matches('*').to_string())
            })
            .collect();
        contigs.push((id.clone(), genes));
    }
    Ok(GenomeProteins { seq_id: label, genome_bp: bp, genes_called: n, gene_call_ms: ms, contigs })
}

/// Domain-scan cached proteins; returns candidate hits (before any cut-off policy) and elapsed ms.
pub fn scan_proteins(
    gp: &GenomeProteins,
    models: &[(DomainFamilyModel, P7Hmm)],
) -> (Vec<BgcGene>, u128) {
    let prefilters: Vec<P7Profile> =
        models.iter().map(|(_, h)| P7Profile::configure(h, 400, true)).collect();
    let t0 = Instant::now();
    let mut all = Vec::new();
    for (id, genes) in &gp.contigs {
        let n_targets = genes.len().max(1) as f64;
        let mut g: Vec<BgcGene> = genes
            .par_iter()
            .enumerate()
            .map(|(i, (beg, end, strand, aa))| BgcGene {
                seq_id: id.clone(),
                gene_index: i,
                beg: *beg,
                end: *end,
                strand: *strand,
                length_aa: aa.len(),
                hits: scan_gene_domains(aa, models, &prefilters, n_targets),
                core_class: None,
                support_weight: 0.0,
                smoothed_support: 0.0,
            })
            .collect();
        g.sort_by_key(|x| x.beg);
        all.extend(g);
    }
    (all, t0.elapsed().as_millis())
}

/// Call clusters for genes spanning several contigs (genes must be sorted by start within contig).
pub fn clusters_from_genes(
    genes: &[BgcGene],
    weights: &HashMap<String, f32>,
    rules: RuleSet,
) -> Vec<BgcCluster> {
    let mut by: BTreeMap<String, Vec<BgcGene>> = BTreeMap::new();
    for g in genes {
        by.entry(g.seq_id.clone()).or_default().push(g.clone());
    }
    let mut out = Vec::new();
    for (id, mut v) in by {
        out.extend(call_clusters_for_contig(&id, &mut v, weights, rules));
    }
    out
}

fn domain_concordance_genes(
    genes: &[BgcGene],
    gecco_features: &Path,
    families: &HashSet<String>,
) -> DomainConcordance {
    let mut nset: HashSet<(String, String)> = HashSet::new();
    for g in genes {
        for h in &g.hits {
            nset.insert((gene_key(&g.seq_id, g.strand, g.beg, g.end), h.accession.clone()));
        }
    }
    let mut gset: HashSet<(String, String)> = HashSet::new();
    if let Ok(txt) = fs::read_to_string(gecco_features) {
        for (i, line) in txt.lines().enumerate() {
            if i == 0 {
                continue;
            }
            let c: Vec<&str> = line.split('\t').collect();
            if c.len() < 8 || !families.contains(c[5]) {
                continue;
            }
            let strand = c[4].chars().next().unwrap_or('+');
            let beg: usize = c[2].parse().unwrap_or(0);
            let end: usize = c[3].parse().unwrap_or(0);
            gset.insert((gene_key(c[0], strand, beg, end), c[5].to_string()));
        }
    }
    let shared = nset.intersection(&gset).count();
    let union = nset.union(&gset).count();
    DomainConcordance {
        native_calls: nset.len(),
        gecco_calls: gset.len(),
        shared,
        jaccard: if union == 0 { 1.0 } else { shared as f64 / union as f64 },
        family_count_pearson_r: 0.0,
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConfigRow {
    pub model_source: String,
    pub policy: String,
    pub clusters: usize,
    pub cluster_bp: usize,
    pub loci_any: usize,
    pub loci_half: usize,
    pub clusters_overlapping_mibig: usize,
    pub domain_native_calls: usize,
    pub domain_gecco_calls: usize,
    pub domain_shared: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Round2Genome {
    pub genome_id: String,
    pub organism: String,
    pub kingdom: String,
    pub set: String,
    pub role: String,
    pub genome_bp: usize,
    pub mibig_loci: usize,
    pub gecco: ToolScore,
    pub gecco_domain_calls_total: usize,
    pub gene_call_ms: u128,
    pub domain_scan_ms_seed_built: u128,
    pub domain_scan_ms_official: u128,
    pub configs: Vec<ConfigRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConfigAggregate {
    pub model_source: String,
    pub policy: String,
    pub genomes: usize,
    pub loci: usize,
    pub loci_any: usize,
    pub loci_half: usize,
    pub clusters: usize,
    pub clusters_overlapping_mibig: usize,
    pub domain_native_calls: usize,
    pub domain_gecco_calls: usize,
    pub domain_shared: usize,
    pub domain_precision_vs_gecco: f64,
    pub domain_recall_vs_gecco: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Round2Report {
    pub pipeline_name: String,
    pub rules: String,
    pub selection_rule: String,
    pub negative_cluster_budget: usize,
    pub chosen_model_source: String,
    pub chosen_policy: String,
    pub chosen_is_baseline: bool,
    pub development_set: Vec<String>,
    pub heldout_set: Vec<String>,
    pub genomes: Vec<Round2Genome>,
    pub development_aggregates: Vec<ConfigAggregate>,
    pub development_selection_scope: Vec<ConfigAggregate>,
    pub negative_control_aggregates: Vec<ConfigAggregate>,
    pub heldout_aggregates: Vec<ConfigAggregate>,
    pub heldout_bacterial_aggregates: Vec<ConfigAggregate>,
    pub gecco_development: ToolAggregate,
    pub gecco_heldout: ToolAggregate,
    pub gecco_heldout_bacterial: ToolAggregate,
    pub gecco_negative_controls: ToolAggregate,
    pub total_duration_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToolAggregate {
    pub genomes: usize,
    pub loci: usize,
    pub loci_any: usize,
    pub loci_half: usize,
    pub clusters: usize,
    pub clusters_overlapping_mibig: usize,
}

fn aggregate_config(
    rows: &[&Round2Genome],
    model_source: &str,
    policy: &str,
) -> ConfigAggregate {
    let mut a = ConfigAggregate {
        model_source: model_source.to_string(),
        policy: policy.to_string(),
        genomes: rows.len(),
        loci: 0,
        loci_any: 0,
        loci_half: 0,
        clusters: 0,
        clusters_overlapping_mibig: 0,
        domain_native_calls: 0,
        domain_gecco_calls: 0,
        domain_shared: 0,
        domain_precision_vs_gecco: 0.0,
        domain_recall_vs_gecco: 0.0,
    };
    for g in rows {
        a.loci += g.mibig_loci;
        if let Some(c) = g.configs.iter().find(|c| c.model_source == model_source && c.policy == policy) {
            a.loci_any += c.loci_any;
            a.loci_half += c.loci_half;
            a.clusters += c.clusters;
            a.clusters_overlapping_mibig += c.clusters_overlapping_mibig;
            a.domain_native_calls += c.domain_native_calls;
            a.domain_gecco_calls += c.domain_gecco_calls;
            a.domain_shared += c.domain_shared;
        }
    }
    a.domain_precision_vs_gecco =
        if a.domain_native_calls == 0 { 0.0 } else { a.domain_shared as f64 / a.domain_native_calls as f64 };
    a.domain_recall_vs_gecco =
        if a.domain_gecco_calls == 0 { 0.0 } else { a.domain_shared as f64 / a.domain_gecco_calls as f64 };
    a
}

fn aggregate_gecco(rows: &[&Round2Genome]) -> ToolAggregate {
    let mut t = ToolAggregate {
        genomes: rows.len(),
        loci: 0,
        loci_any: 0,
        loci_half: 0,
        clusters: 0,
        clusters_overlapping_mibig: 0,
    };
    for g in rows {
        t.loci += g.mibig_loci;
        t.loci_any += g.gecco.mibig_truth_hit;
        t.loci_half += (g.gecco.mibig_recall_half_covered * g.mibig_loci as f64).round() as usize;
        t.clusters += g.gecco.n_clusters;
        t.clusters_overlapping_mibig +=
            (g.gecco.predicted_overlapping_mibig_frac * g.gecco.n_clusters as f64).round() as usize;
    }
    t
}

/// Calibration round 2. Rules are frozen at `RuleSet::V2`; only the model source and the
/// domain cut-off policy vary, and the configuration is chosen **on the development set only**.
pub fn run_bgc_calibration_round2(ws_root: &Path) -> Result<Round2Report, String> {
    let t_start = Instant::now();
    let labs = if ws_root.ends_with("crvscience-labs") {
        ws_root.to_path_buf()
    } else {
        ws_root.join("crvscience-labs")
    };
    let seed_dir = labs.join("tools/gecco/pfam_seed");
    let official_dir = labs.join("tools/gecco/pfam_official");
    let mibig_dir = labs.join("tools/gecco/data/mibig_json_4.0");
    let gecco_exe = labs.join("tools/gecco/venv/Scripts/gecco.exe");
    let work = labs.join("runs/gecco_work");
    let cache = labs.join("cache/bgc");

    let seed_models = build_domain_models_from(&seed_dir, None)?;
    let official_models = build_domain_models_from(&seed_dir, Some(&official_dir))?;
    let n_official = official_models.iter().filter(|(f, _)| f.model_source == "official-pfam").count();
    let weights: HashMap<String, f32> =
        seed_models.iter().map(|(f, _)| (f.accession.clone(), f.accessory_weight)).collect();
    let ga_map: HashMap<String, f32> =
        seed_models.iter().map(|(f, _)| (f.accession.clone(), f.ga_bits)).collect();
    let family_set: HashSet<String> = seed_models.iter().map(|(f, _)| f.accession.clone()).collect();
    eprintln!("[bgc-calibrate] official Pfam HMMs loaded for {}/34 families", n_official);

    // (accession, organism, kingdom, set, role)
    let specs: Vec<(&str, &str, &str, &str, &str)> = vec![
        ("AL645882.2", "Streptomyces coelicolor A3(2)", "Prokaryota", "development", "bgc-rich"),
        ("CP000560.1", "Bacillus velezensis FZB42", "Prokaryota", "development", "bgc-rich"),
        ("AE004091.2", "Pseudomonas aeruginosa PAO1", "Prokaryota", "development", "bgc-rich"),
        ("U00096.3", "Escherichia coli K-12 MG1655", "Prokaryota", "development", "bgc-poor"),
        ("BN001302.1", "Aspergillus nidulans FGSC A4 chr II", "Eukaryota", "development", "out-of-regime"),
        ("NC_000908.2", "Mycoplasmoides genitalium G37", "Prokaryota", "development", "negative-control"),
        ("NC_000915.1", "Helicobacter pylori 26695", "Prokaryota", "development", "negative-control"),
        ("NC_000907.1", "Haemophilus influenzae Rd KW20", "Prokaryota", "development", "negative-control"),
        ("NC_001136.10", "Saccharomyces cerevisiae S288C chr IV", "Eukaryota", "development", "negative-control"),
        ("AP009493.1", "Streptomyces griseus NBRC 13350", "Prokaryota", "held-out", "bgc-rich"),
        ("CP000667.1", "Salinispora tropica CNB-440", "Prokaryota", "held-out", "bgc-rich"),
        ("BX470251.1", "Photorhabdus laumondii TTO1", "Prokaryota", "held-out", "bgc-rich"),
        ("CP000113.1", "Myxococcus xanthus DK 1622", "Prokaryota", "held-out", "bgc-rich"),
        ("CP000085.1", "Burkholderia thailandensis E264 chr II", "Prokaryota", "held-out", "bgc-rich"),
        ("CM000574.1", "Fusarium graminearum PH-1 chr 1", "Eukaryota", "held-out", "out-of-regime"),
    ];

    let mut genomes: Vec<Round2Genome> = Vec::new();
    for (acc, organism, kingdom, set, role) in &specs {
        let fasta = cache.join(format!("{}.fasta", acc));
        let base = acc.split('.').next().unwrap_or(acc);
        let dir = work.join("m080").join(acc);
        if !dir.join(format!("{}.clusters.tsv", acc)).exists() {
            run_gecco(&gecco_exe, &fasta, &dir, 0.8)?;
        }
        let gecco_iv = parse_gecco_clusters(&dir.join(format!("{}.clusters.tsv", acc)));
        let features = dir.join(format!("{}.features.tsv", acc));
        let truth = load_mibig_truth(&mibig_dir, base);

        eprintln!("[bgc-calibrate] {} ({}): gene calling...", acc, role);
        let gp = call_genome_proteins(&fasta)?;
        let seq_id = gp.seq_id.clone();
        let gecco_score = score_tool("GECCO v0.11.0 (-m 0.8)", &gecco_iv, &truth, &seq_id, 0);
        let gecco_dom_total = domain_concordance_genes(&[], &features, &family_set).gecco_calls;

        let mut configs: Vec<ConfigRow> = Vec::new();
        let mut scan_ms = [0u128; 2];
        for (si, (src_name, models)) in
            [("seed-built", &seed_models), ("official-pfam", &official_models)].iter().enumerate()
        {
            eprintln!("[bgc-calibrate] {}: domain scan ({})...", acc, src_name);
            let (raw, ms) = scan_proteins(&gp, models);
            scan_ms[si] = ms;
            for policy in CutoffPolicy::all() {
                let filtered = filter_hits(&raw, &ga_map, policy);
                let cl = clusters_from_genes(&filtered, &weights, RuleSet::V2);
                let iv: Vec<Interval> = cl
                    .iter()
                    .map(|c| Interval { seq_id: c.seq_id.clone(), start: c.start, end: c.end, label: c.classes.join(";") })
                    .collect();
                let sc = score_tool(src_name, &iv, &truth, &seq_id, 0);
                let dc = domain_concordance_genes(&filtered, &features, &family_set);
                configs.push(ConfigRow {
                    model_source: src_name.to_string(),
                    policy: policy.label().to_string(),
                    clusters: sc.n_clusters,
                    cluster_bp: sc.total_bp,
                    loci_any: sc.mibig_truth_hit,
                    loci_half: (sc.mibig_recall_half_covered * truth.len() as f64).round() as usize,
                    clusters_overlapping_mibig: (sc.predicted_overlapping_mibig_frac * sc.n_clusters as f64).round() as usize,
                    domain_native_calls: dc.native_calls,
                    domain_gecco_calls: dc.gecco_calls,
                    domain_shared: dc.shared,
                });
            }
        }
        genomes.push(Round2Genome {
            genome_id: acc.to_string(),
            organism: organism.to_string(),
            kingdom: kingdom.to_string(),
            set: set.to_string(),
            role: role.to_string(),
            genome_bp: gp.genome_bp,
            mibig_loci: truth.len(),
            gecco: gecco_score,
            gecco_domain_calls_total: gecco_dom_total,
            gene_call_ms: gp.gene_call_ms,
            domain_scan_ms_seed_built: scan_ms[0],
            domain_scan_ms_official: scan_ms[1],
            configs,
        });
    }

    // ---- aggregates ----
    let dev: Vec<&Round2Genome> = genomes.iter().filter(|g| g.set == "development").collect();
    let held: Vec<&Round2Genome> = genomes.iter().filter(|g| g.set == "held-out").collect();
    let held_bact: Vec<&Round2Genome> = held.iter().copied().filter(|g| g.kingdom == "Prokaryota").collect();
    let neg: Vec<&Round2Genome> = dev.iter().copied().filter(|g| g.role == "negative-control").collect();
    // Selection scope: bacterial development genomes (rich or poor); negatives gate the choice.
    let dev_rich: Vec<&Round2Genome> = dev
        .iter()
        .copied()
        .filter(|g| g.role == "bgc-rich" || g.role == "bgc-poor")
        .collect();

    let mut dev_agg = Vec::new();
    let mut scope_agg = Vec::new();
    let mut neg_agg = Vec::new();
    let mut held_agg = Vec::new();
    let mut held_bact_agg = Vec::new();
    for src in ["seed-built", "official-pfam"] {
        for p in CutoffPolicy::all() {
            dev_agg.push(aggregate_config(&dev, src, p.label()));
            scope_agg.push(aggregate_config(&dev_rich, src, p.label()));
            neg_agg.push(aggregate_config(&neg, src, p.label()));
            held_agg.push(aggregate_config(&held, src, p.label()));
            held_bact_agg.push(aggregate_config(&held_bact, src, p.label()));
        }
    }

    // Pre-registered selection: among configurations whose total cluster count on the BGC-free
    // development negative controls is within the budget, maximise curated loci overlapped on the
    // bacterial development genomes; break ties by fewer clusters, then by baseline preference.
    const NEG_BUDGET: usize = 2;
    let mut best: Option<(usize, usize, usize, String, String)> = None; // (loci, -clusters via inverse, baseline flag)
    for (i, s) in scope_agg.iter().enumerate() {
        let n = &neg_agg[i];
        if n.clusters > NEG_BUDGET {
            continue;
        }
        let baseline = (s.model_source == "seed-built" && s.policy == "GA") as usize;
        let key = (s.loci_any, usize::MAX - s.clusters, baseline, s.model_source.clone(), s.policy.clone());
        if best.as_ref().map(|b| (b.0, b.1, b.2) < (key.0, key.1, key.2)).unwrap_or(true) {
            best = Some(key);
        }
    }
    let (chosen_src, chosen_pol) = best
        .map(|b| (b.3, b.4))
        .unwrap_or(("seed-built".to_string(), "GA".to_string()));

    Ok(Round2Report {
        pipeline_name: "bgc-calibrate (round 2: model source x cut-off policy, negative controls, held-out genomes)".to_string(),
        rules: "RuleSet::V2 (frozen; not modified in this round)".to_string(),
        selection_rule: format!(
            "Among configurations with <= {} total clusters on the BGC-free development negative controls, \
             maximise curated MIBiG loci overlapped on the bacterial development genomes; ties -> fewer clusters -> baseline. \
             Held-out genomes are scored only after the choice is made.",
            NEG_BUDGET
        ),
        negative_cluster_budget: NEG_BUDGET,
        chosen_is_baseline: chosen_src == "seed-built" && chosen_pol == "GA",
        chosen_model_source: chosen_src,
        chosen_policy: chosen_pol,
        development_set: dev.iter().map(|g| g.genome_id.clone()).collect(),
        heldout_set: held.iter().map(|g| g.genome_id.clone()).collect(),
        development_aggregates: dev_agg,
        development_selection_scope: scope_agg,
        negative_control_aggregates: neg_agg,
        heldout_aggregates: held_agg,
        heldout_bacterial_aggregates: held_bact_agg,
        gecco_development: aggregate_gecco(&dev),
        gecco_heldout: aggregate_gecco(&held),
        gecco_heldout_bacterial: aggregate_gecco(&held_bact),
        gecco_negative_controls: aggregate_gecco(&neg),
        genomes,
        total_duration_ms: t_start.elapsed().as_millis(),
    })
}

// ---------------------------------------------------------------------------------------------
// Calibration round 3: rule set V3 + parameter search on the development genomes, then a single
// scoring pass on a reserve set that is not touched until the configuration is frozen on disk.
// ---------------------------------------------------------------------------------------------

pub fn clusters_from_genes_with(
    genes: &[BgcGene],
    weights: &HashMap<String, f32>,
    rules: RuleSet,
    prm: &ClusterParams,
) -> Vec<BgcCluster> {
    let mut by: BTreeMap<String, Vec<BgcGene>> = BTreeMap::new();
    for g in genes {
        by.entry(g.seq_id.clone()).or_default().push(g.clone());
    }
    let mut out = Vec::new();
    for (id, mut v) in by {
        out.extend(call_clusters_with_params(&id, &mut v, weights, rules, prm));
    }
    out
}

/// Keep only hits from the 34 published families (used to evaluate V2 on an extended scan).
pub fn restrict_to_published(genes: &[BgcGene]) -> Vec<BgcGene> {
    let pub_set: HashSet<&str> = BGC_FAMILIES.iter().map(|f| f.0).collect();
    genes
        .iter()
        .map(|g| {
            let mut c = g.clone();
            c.hits.retain(|h| pub_set.contains(h.accession.as_str()));
            c
        })
        .collect()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LocusRow {
    pub genome_id: String,
    pub mibig_id: String,
    pub classes: String,
    pub start: usize,
    pub end: usize,
    pub coverage_fraction: f64,
    pub best_iou: f64,
}

/// `(loci overlapped, loci >= half covered, mean best-IoU over overlapped loci, per-locus rows)`.
fn locus_metrics(
    genome_id: &str,
    seq_id: &str,
    truth: &[(String, usize, usize, String, String)],
    preds: &[Interval],
) -> (usize, usize, f64, Vec<LocusRow>) {
    let (mut any, mut half) = (0usize, 0usize);
    let (mut iou_sum, mut iou_n) = (0.0f64, 0usize);
    let mut rows = Vec::new();
    for t in truth {
        let ti = Interval { seq_id: seq_id.to_string(), start: t.1, end: t.2, label: t.0.clone() };
        let cov = covered_fraction(&ti, preds);
        let mut best = 0.0f64;
        for p in preds.iter().filter(|p| p.seq_id == seq_id) {
            let o = overlap_bp(t.1, t.2, p.start, p.end);
            if o > 0 {
                let uni = (t.2 - t.1 + 1) + (p.end - p.start + 1) - o;
                best = best.max(o as f64 / uni as f64);
            }
        }
        if cov > 0.0 {
            any += 1;
            iou_sum += best;
            iou_n += 1;
        }
        if cov >= 0.5 {
            half += 1;
        }
        rows.push(LocusRow {
            genome_id: genome_id.to_string(),
            mibig_id: t.0.clone(),
            classes: t.3.clone(),
            start: t.1,
            end: t.2,
            coverage_fraction: cov,
            best_iou: best,
        });
    }
    (any, half, if iou_n == 0 { 0.0 } else { iou_sum / iou_n as f64 }, rows)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct R3Config {
    pub label: String,
    pub rules: String,
    pub params: ClusterParams,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct R3Score {
    pub label: String,
    pub loci: usize,
    pub loci_any: usize,
    pub loci_half: usize,
    pub mean_iou: f64,
    pub clusters: usize,
    pub clusters_overlapping: usize,
    pub negative_clusters: usize,
    pub overlap_share: f64,
    pub recall: f64,
    pub f1: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct R3GenomeResult {
    pub genome_id: String,
    pub organism: String,
    pub kingdom: String,
    pub role: String,
    pub genome_bp: usize,
    pub mibig_loci: usize,
    /// `(label, loci_any, loci_half, clusters, clusters_overlapping, mean_iou)`
    pub configs: Vec<(String, usize, usize, usize, usize, f64)>,
    pub gecco: (usize, usize, usize, usize, f64),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Round3Report {
    pub pipeline_name: String,
    pub protocol: String,
    pub grid_size: usize,
    pub objective: String,
    pub negative_budget: usize,
    pub frozen: R3Config,
    pub dev_top_configs: Vec<(R3Config, R3Score)>,
    pub dev_baselines: Vec<R3Score>,
    pub dev_gecco: R3Score,
    pub reserve_scores: Vec<R3Score>,
    pub reserve_gecco: R3Score,
    pub reserve_genomes: Vec<R3GenomeResult>,
    pub dev_genomes: Vec<R3GenomeResult>,
    pub reserve_loci_frozen: Vec<LocusRow>,
    pub reserve_loci_gecco: Vec<LocusRow>,
    pub dev_loci_frozen: Vec<LocusRow>,
    pub families_added: Vec<String>,
    pub total_duration_ms: u128,
}

struct GenomeCtx {
    acc: String,
    organism: String,
    kingdom: String,
    role: String,
    seq_id: String,
    genome_bp: usize,
    truth: Vec<(String, usize, usize, String, String)>,
    gecco_iv: Vec<Interval>,
    /// Official-Pfam extended-scan candidate hits, GA-filtered.
    genes_ext: Vec<BgcGene>,
}

fn score_cfg(
    label: &str,
    ctxs: &[&GenomeCtx],
    clusters_for: &dyn Fn(&GenomeCtx) -> Vec<BgcCluster>,
) -> R3Score {
    let (mut loci, mut any, mut half, mut ncl, mut ovl, mut neg) = (0usize, 0usize, 0usize, 0usize, 0usize, 0usize);
    let (mut iou_sum, mut iou_n) = (0.0f64, 0usize);
    for c in ctxs {
        let cl = clusters_for(c);
        let iv: Vec<Interval> = cl
            .iter()
            .map(|x| Interval { seq_id: x.seq_id.clone(), start: x.start, end: x.end, label: x.classes.join(";") })
            .collect();
        if c.role == "negative-control" {
            neg += iv.len();
            continue;
        }
        let (a, h, mi, _) = locus_metrics(&c.acc, &c.seq_id, &c.truth, &iv);
        loci += c.truth.len();
        any += a;
        half += h;
        iou_sum += mi * a as f64;
        iou_n += a;
        ncl += iv.len();
        ovl += iv.iter().filter(|p| c.truth.iter().any(|t| overlap_bp(p.start, p.end, t.1, t.2) > 0)).count();
    }
    let recall = if loci == 0 { 0.0 } else { any as f64 / loci as f64 };
    let share = if ncl == 0 { 0.0 } else { ovl as f64 / ncl as f64 };
    R3Score {
        label: label.to_string(),
        loci,
        loci_any: any,
        loci_half: half,
        mean_iou: if iou_n == 0 { 0.0 } else { iou_sum / iou_n as f64 },
        clusters: ncl,
        clusters_overlapping: ovl,
        negative_clusters: neg,
        overlap_share: share,
        recall,
        f1: if recall + share == 0.0 { 0.0 } else { 2.0 * recall * share / (recall + share) },
    }
}

fn score_gecco(label: &str, ctxs: &[&GenomeCtx]) -> R3Score {
    let (mut loci, mut any, mut half, mut ncl, mut ovl, mut neg) = (0usize, 0usize, 0usize, 0usize, 0usize, 0usize);
    let (mut iou_sum, mut iou_n) = (0.0f64, 0usize);
    for c in ctxs {
        if c.role == "negative-control" {
            neg += c.gecco_iv.len();
            continue;
        }
        let (a, h, mi, _) = locus_metrics(&c.acc, &c.seq_id, &c.truth, &c.gecco_iv);
        loci += c.truth.len();
        any += a;
        half += h;
        iou_sum += mi * a as f64;
        iou_n += a;
        ncl += c.gecco_iv.len();
        ovl += c.gecco_iv.iter().filter(|p| c.truth.iter().any(|t| overlap_bp(p.start, p.end, t.1, t.2) > 0)).count();
    }
    let recall = if loci == 0 { 0.0 } else { any as f64 / loci as f64 };
    let share = if ncl == 0 { 0.0 } else { ovl as f64 / ncl as f64 };
    R3Score {
        label: label.to_string(),
        loci,
        loci_any: any,
        loci_half: half,
        mean_iou: if iou_n == 0 { 0.0 } else { iou_sum / iou_n as f64 },
        clusters: ncl,
        clusters_overlapping: ovl,
        negative_clusters: neg,
        overlap_share: share,
        recall,
        f1: if recall + share == 0.0 { 0.0 } else { 2.0 * recall * share / (recall + share) },
    }
}

#[allow(clippy::too_many_arguments)]
fn build_ctx(
    spec: &(&str, &str, &str, &str),
    labs: &Path,
    mibig_dir: &Path,
    gecco_exe: &Path,
    models_ext: &[(DomainFamilyModel, P7Hmm)],
    ga_map: &HashMap<String, f32>,
    cache_dir: Option<&Path>,
    keep_gp: bool,
) -> Result<(GenomeCtx, Option<GenomeProteins>), String> {
    let (acc, organism, kingdom, role) = *spec;
    let work = labs.join("runs/gecco_work");
    let fasta = labs.join("cache/bgc").join(format!("{}.fasta", acc));
    let dir = work.join("m080").join(acc);
    if !dir.join(format!("{}.clusters.tsv", acc)).exists() {
        run_gecco(gecco_exe, &fasta, &dir, 0.8)?;
    }
    let gecco_iv = parse_gecco_clusters(&dir.join(format!("{}.clusters.tsv", acc)));
    let truth = load_mibig_truth(mibig_dir, acc.split('.').next().unwrap_or(acc));

    let cache_file = cache_dir.map(|d| d.join(format!("{}.raw.json", acc)));
    let mut gp_out: Option<GenomeProteins> = None;
    let (seq_id, genome_bp, raw): (String, usize, Vec<BgcGene>) = match cache_file
        .as_ref()
        .filter(|p| p.exists())
        .and_then(|p| fs::read_to_string(p).ok())
        .and_then(|t| serde_json::from_str::<(String, usize, Vec<BgcGene>)>(&t).ok())
    {
        Some(x) => x,
        None => {
            eprintln!("[bgc-round3] {} ({}): gene calling + extended domain scan...", acc, role);
            let gp = call_genome_proteins(&fasta)?;
            let (raw, _) = scan_proteins(&gp, models_ext);
            let tuple = (gp.seq_id.clone(), gp.genome_bp, raw);
            if let Some(p) = cache_file.as_ref() {
                let _ = fs::create_dir_all(p.parent().unwrap());
                let _ = fs::write(p, serde_json::to_string(&tuple).unwrap_or_default());
            }
            if keep_gp {
                gp_out = Some(gp);
            }
            tuple
        }
    };
    let genes_ext = filter_hits(&raw, ga_map, CutoffPolicy::Ga);
    Ok((
        GenomeCtx {
            acc: acc.to_string(),
            organism: organism.to_string(),
            kingdom: kingdom.to_string(),
            role: role.to_string(),
            seq_id,
            genome_bp,
            truth,
            gecco_iv,
            genes_ext,
        },
        gp_out,
    ))
}

pub fn run_bgc_calibration_round3(ws_root: &Path) -> Result<Round3Report, String> {
    let t_start = Instant::now();
    let labs = if ws_root.ends_with("crvscience-labs") {
        ws_root.to_path_buf()
    } else {
        ws_root.join("crvscience-labs")
    };
    let seed_dir = labs.join("tools/gecco/pfam_seed");
    let official_dir = labs.join("tools/gecco/pfam_official");
    let mibig_dir = labs.join("tools/gecco/data/mibig_json_4.0");
    let gecco_exe = labs.join("tools/gecco/venv/Scripts/gecco.exe");
    let cache_dir = labs.join("runs/bgc_work/cache_v3");

    let models_ext = build_domain_models_ext(&seed_dir, Some(&official_dir), true)?;
    let models_seed34 = build_domain_models_ext(&seed_dir, None, false)?;
    let weights_ext: HashMap<String, f32> =
        models_ext.iter().map(|(f, _)| (f.accession.clone(), f.accessory_weight)).collect();
    let weights_34: HashMap<String, f32> =
        models_seed34.iter().map(|(f, _)| (f.accession.clone(), f.accessory_weight)).collect();
    let ga_ext: HashMap<String, f32> =
        models_ext.iter().map(|(f, _)| (f.accession.clone(), f.ga_bits)).collect();
    let ga_34: HashMap<String, f32> =
        models_seed34.iter().map(|(f, _)| (f.accession.clone(), f.ga_bits)).collect();

    // ------------------------------- development genomes -------------------------------
    let dev_specs: Vec<(&str, &str, &str, &str)> = vec![
        ("AL645882.2", "Streptomyces coelicolor A3(2)", "Prokaryota", "bgc-rich"),
        ("CP000560.1", "Bacillus velezensis FZB42", "Prokaryota", "bgc-rich"),
        ("AE004091.2", "Pseudomonas aeruginosa PAO1", "Prokaryota", "bgc-rich"),
        ("U00096.3", "Escherichia coli K-12 MG1655", "Prokaryota", "bgc-poor"),
        ("AP009493.1", "Streptomyces griseus NBRC 13350", "Prokaryota", "bgc-rich"),
        ("CP000667.1", "Salinispora tropica CNB-440", "Prokaryota", "bgc-rich"),
        ("BX470251.1", "Photorhabdus laumondii TTO1", "Prokaryota", "bgc-rich"),
        ("CP000113.1", "Myxococcus xanthus DK 1622", "Prokaryota", "bgc-rich"),
        ("CP000085.1", "Burkholderia thailandensis E264 chr II", "Prokaryota", "bgc-rich"),
        ("BN001302.1", "Aspergillus nidulans FGSC A4 chr II", "Eukaryota", "out-of-regime"),
        ("CM000574.1", "Fusarium graminearum PH-1 chr 1", "Eukaryota", "out-of-regime"),
        ("NC_000908.2", "Mycoplasmoides genitalium G37", "Prokaryota", "negative-control"),
        ("NC_000915.1", "Helicobacter pylori 26695", "Prokaryota", "negative-control"),
        ("NC_000907.1", "Haemophilus influenzae Rd KW20", "Prokaryota", "negative-control"),
        ("NC_001136.10", "Saccharomyces cerevisiae S288C chr IV", "Eukaryota", "negative-control"),
    ];
    let mut dev: Vec<GenomeCtx> = Vec::new();
    for s in &dev_specs {
        let (c, _) = build_ctx(s, &labs, &mibig_dir, &gecco_exe, &models_ext, &ga_ext, Some(&cache_dir), false)?;
        dev.push(c);
    }
    // Selection scope: bacterial genomes with curated loci, plus every BGC-free negative control.
    let dev_scope: Vec<&GenomeCtx> = dev
        .iter()
        .filter(|c| (c.kingdom == "Prokaryota" && c.role != "out-of-regime") || c.role == "negative-control")
        .collect();
    let dev_all: Vec<&GenomeCtx> = dev.iter().collect();

    // Pre-compute published-family-restricted genes for V2 candidates.
    let dev_restricted: HashMap<String, Vec<BgcGene>> =
        dev.iter().map(|c| (c.acc.clone(), restrict_to_published(&c.genes_ext))).collect();

    // ------------------------------- parameter grid -------------------------------
    let mut grid: Vec<(RuleSet, ClusterParams)> = Vec::new();
    let default = ClusterParams::default();
    grid.push((RuleSet::V2, default));
    grid.push((RuleSet::V3, default));
    for rules in [RuleSet::V2, RuleSet::V3] {
        for sw in [3usize, 4, 6] {
            for sm in [1.2f32, 1.5, 1.8, 2.2] {
                for eg in [3000usize, 5000, 8000] {
                    for pad in [0usize, 1000, 2000, 4000] {
                        let prm = ClusterParams { support_window: sw, seed_min: sm, extend_gap: eg, extend_min: 0.9, pad };
                        if prm != default {
                            grid.push((rules, prm));
                        }
                    }
                }
            }
        }
    }
    let grid_size = grid.len();
    eprintln!("[bgc-round3] development grid: {} configurations", grid_size);

    let run_cfg = |rules: RuleSet, prm: &ClusterParams, c: &GenomeCtx| -> Vec<BgcCluster> {
        if rules == RuleSet::V3 {
            clusters_from_genes_with(&c.genes_ext, &weights_ext, rules, prm)
        } else {
            clusters_from_genes_with(&dev_restricted[&c.acc], &weights_34, rules, prm)
        }
    };

    const NEG_BUDGET: usize = 1;
    let mut scored: Vec<(R3Config, R3Score)> = Vec::new();
    for (rules, prm) in &grid {
        let label = format!(
            "{:?} w={} seed>={} gap={} pad={}",
            rules, prm.support_window, prm.seed_min, prm.extend_gap, prm.pad
        );
        let sc = score_cfg(&label, &dev_scope, &|c| run_cfg(*rules, prm, c));
        scored.push((R3Config { label, rules: format!("{:?}", rules), params: *prm }, sc));
    }
    // Pre-registered objective: negative-control budget, then F1 of (locus recall, MIBiG-overlap
    // share), then mean best-IoU of overlapped loci, then fewer clusters; ties keep grid order.
    let mut order: Vec<usize> = (0..scored.len()).filter(|&i| scored[i].1.negative_clusters <= NEG_BUDGET).collect();
    if order.is_empty() {
        order = (0..scored.len()).collect();
    }
    order.sort_by(|&a, &b| {
        let (sa, sb) = (&scored[a].1, &scored[b].1);
        sb.f1
            .partial_cmp(&sa.f1)
            .unwrap()
            .then(sb.mean_iou.partial_cmp(&sa.mean_iou).unwrap())
            .then(sa.clusters.cmp(&sb.clusters))
            .then(a.cmp(&b))
    });
    let frozen = scored[order[0]].0.clone();
    let dev_top: Vec<(R3Config, R3Score)> = order.iter().take(8).map(|&i| scored[i].clone()).collect();
    let frozen_rules = if frozen.rules == "V3" { RuleSet::V3 } else { RuleSet::V2 };

    // Persist the frozen configuration BEFORE any reserve genome is touched.
    let frozen_path = labs.join("runs/bgc-round3-frozen.json");
    fs::write(&frozen_path, serde_json::to_string_pretty(&frozen).unwrap_or_default())
        .map_err(|e| e.to_string())?;
    eprintln!("[bgc-round3] FROZEN on development set: {} (written to {})", frozen.label, frozen_path.display());

    let dev_baselines = vec![
        score_cfg("published #59: V2 default (34 families)", &dev_scope, &|c| {
            run_cfg(RuleSet::V2, &default, c)
        }),
        score_cfg("V3 default parameters", &dev_scope, &|c| run_cfg(RuleSet::V3, &default, c)),
        score_cfg("frozen", &dev_scope, &|c| run_cfg(frozen_rules, &frozen.params, c)),
    ];
    let dev_gecco = score_gecco("GECCO v0.11.0 (-m 0.8)", &dev_scope);

    // ------------------------------- reserve genomes -------------------------------
    let reserve_specs: Vec<(&str, &str, &str, &str)> = vec![
        ("FN667742.1", "Xenorhabdus nematophila ATCC 19061", "Prokaryota", "bgc-rich"),
        ("FN554889.1", "Streptomyces scabiei 87.22", "Prokaryota", "bgc-rich"),
        ("BA000030.4", "Streptomyces avermitilis MA-4680", "Prokaryota", "bgc-rich"),
        ("AM420293.1", "Saccharopolyspora erythraea NRRL 2338", "Prokaryota", "bgc-rich"),
        ("CP000850.1", "Salinispora arenicola CNS-205", "Prokaryota", "bgc-rich"),
        ("NC_002528.1", "Buchnera aphidicola APS", "Prokaryota", "negative-control"),
        ("NC_000117.1", "Chlamydia trachomatis D/UW-3", "Prokaryota", "negative-control"),
        ("NC_000919.1", "Treponema pallidum Nichols", "Prokaryota", "negative-control"),
        ("NC_003421.2", "Schizosaccharomyces pombe chr III", "Eukaryota", "negative-control"),
    ];
    let mut reserve: Vec<GenomeCtx> = Vec::new();
    let mut reserve_seed34: HashMap<String, Vec<BgcGene>> = HashMap::new();
    for s in &reserve_specs {
        // Reserve genomes are scanned only now. One gene-calling pass is shared by both model sets.
        let (c, gp_opt) = build_ctx(s, &labs, &mibig_dir, &gecco_exe, &models_ext, &ga_ext, Some(&cache_dir), true)?;
        let gp = match gp_opt {
            Some(gp) => gp,
            None => call_genome_proteins(&labs.join("cache/bgc").join(format!("{}.fasta", s.0)))?,
        };
        let (raw34, _) = scan_proteins(&gp, &models_seed34);
        reserve_seed34.insert(s.0.to_string(), filter_hits(&raw34, &ga_34, CutoffPolicy::Ga));
        reserve.push(c);
    }
    let reserve_scope: Vec<&GenomeCtx> = reserve.iter().collect();
    let reserve_restricted: HashMap<String, Vec<BgcGene>> =
        reserve.iter().map(|c| (c.acc.clone(), restrict_to_published(&c.genes_ext))).collect();

    let res_published = |c: &GenomeCtx| -> Vec<BgcCluster> {
        clusters_from_genes_with(&reserve_seed34[&c.acc], &weights_34, RuleSet::V2, &default)
    };
    let res_official34 = |c: &GenomeCtx| -> Vec<BgcCluster> {
        clusters_from_genes_with(&reserve_restricted[&c.acc], &weights_34, RuleSet::V2, &default)
    };
    let res_v3_default = |c: &GenomeCtx| -> Vec<BgcCluster> {
        clusters_from_genes_with(&c.genes_ext, &weights_ext, RuleSet::V3, &default)
    };
    let res_frozen = |c: &GenomeCtx| -> Vec<BgcCluster> {
        if frozen_rules == RuleSet::V3 {
            clusters_from_genes_with(&c.genes_ext, &weights_ext, RuleSet::V3, &frozen.params)
        } else {
            clusters_from_genes_with(&reserve_restricted[&c.acc], &weights_34, RuleSet::V2, &frozen.params)
        }
    };
    let reserve_scores = vec![
        score_cfg("published #59 (seed-built, V2, default)", &reserve_scope, &res_published),
        score_cfg("official Pfam models, V2, default", &reserve_scope, &res_official34),
        score_cfg("V3 rules, default parameters", &reserve_scope, &res_v3_default),
        score_cfg("FROZEN configuration", &reserve_scope, &res_frozen),
    ];
    let reserve_gecco = score_gecco("GECCO v0.11.0 (-m 0.8)", &reserve_scope);

    // ------------------------------- per-genome & per-locus detail -------------------------------
    let genome_results = |ctxs: &[GenomeCtx], extra: &[(&str, &dyn Fn(&GenomeCtx) -> Vec<BgcCluster>)]| -> Vec<R3GenomeResult> {
        ctxs.iter()
            .map(|c| {
                let mut configs = Vec::new();
                for (label, f) in extra {
                    let cl = f(c);
                    let iv: Vec<Interval> = cl
                        .iter()
                        .map(|x| Interval { seq_id: x.seq_id.clone(), start: x.start, end: x.end, label: x.classes.join(";") })
                        .collect();
                    let (a, h, mi, _) = locus_metrics(&c.acc, &c.seq_id, &c.truth, &iv);
                    let ov = iv.iter().filter(|p| c.truth.iter().any(|t| overlap_bp(p.start, p.end, t.1, t.2) > 0)).count();
                    configs.push((label.to_string(), a, h, iv.len(), ov, mi));
                }
                let (ga, gh, gmi, _) = locus_metrics(&c.acc, &c.seq_id, &c.truth, &c.gecco_iv);
                let gov = c.gecco_iv.iter().filter(|p| c.truth.iter().any(|t| overlap_bp(p.start, p.end, t.1, t.2) > 0)).count();
                R3GenomeResult {
                    genome_id: c.acc.clone(),
                    organism: c.organism.clone(),
                    kingdom: c.kingdom.clone(),
                    role: c.role.clone(),
                    genome_bp: c.genome_bp,
                    mibig_loci: c.truth.len(),
                    configs,
                    gecco: (ga, gh, c.gecco_iv.len(), gov, gmi),
                }
            })
            .collect()
    };
    let dev_pub = |c: &GenomeCtx| run_cfg(RuleSet::V2, &default, c);
    let dev_frz = |c: &GenomeCtx| run_cfg(frozen_rules, &frozen.params, c);
    let dev_genomes = genome_results(&dev, &[("published V2 default (official models)", &dev_pub), ("frozen", &dev_frz)]);
    let reserve_genomes = genome_results(
        &reserve,
        &[
            ("published #59", &res_published),
            ("official V2", &res_official34),
            ("V3 default", &res_v3_default),
            ("frozen", &res_frozen),
        ],
    );
    let mut reserve_loci_frozen = Vec::new();
    let mut reserve_loci_gecco = Vec::new();
    for c in &reserve {
        let cl = res_frozen(c);
        let iv: Vec<Interval> = cl
            .iter()
            .map(|x| Interval { seq_id: x.seq_id.clone(), start: x.start, end: x.end, label: x.classes.join(";") })
            .collect();
        reserve_loci_frozen.extend(locus_metrics(&c.acc, &c.seq_id, &c.truth, &iv).3);
        reserve_loci_gecco.extend(locus_metrics(&c.acc, &c.seq_id, &c.truth, &c.gecco_iv).3);
    }
    let mut dev_loci_frozen = Vec::new();
    for c in dev.iter().filter(|c| c.role != "negative-control") {
        let cl = dev_frz(c);
        let iv: Vec<Interval> = cl
            .iter()
            .map(|x| Interval { seq_id: x.seq_id.clone(), start: x.start, end: x.end, label: x.classes.join(";") })
            .collect();
        dev_loci_frozen.extend(locus_metrics(&c.acc, &c.seq_id, &c.truth, &iv).3);
    }
    let _ = dev_all;

    Ok(Round3Report {
        pipeline_name: "bgc-round3 (rule set V3 + parameter search; reserve-set evaluation)".to_string(),
        protocol: "Development = all 15 genomes examined in rounds 1-2. A grid of rule sets and cluster parameters is scored on the \
                   development scope only; the winner is written to disk; the 9 reserve genomes are gene-called and scored afterwards."
            .to_string(),
        grid_size,
        objective: format!(
            "<= {} cluster(s) on development negative controls; then maximise F1 of (locus recall, MIBiG-overlap share); \
             then mean best-IoU; then fewer clusters; ties keep grid order (defaults first).",
            NEG_BUDGET
        ),
        negative_budget: NEG_BUDGET,
        frozen,
        dev_top_configs: dev_top,
        dev_baselines,
        dev_gecco,
        reserve_scores,
        reserve_gecco,
        reserve_genomes,
        dev_genomes,
        reserve_loci_frozen,
        reserve_loci_gecco,
        dev_loci_frozen,
        families_added: BGC_FAMILIES_EXTRA.iter().map(|f| format!("{} {}", f.0, f.1)).collect(),
        total_duration_ms: t_start.elapsed().as_millis(),
    })
}

// ---------------------------------------------------------------------------------------------
// Production entry point: frozen V3 configuration on a user-supplied genome.
// ---------------------------------------------------------------------------------------------

/// Cluster parameters frozen by calibration round 3 (development-set grid search).
pub const FROZEN_V3_PARAMS: ClusterParams = ClusterParams {
    support_window: 6,
    seed_min: 1.2,
    extend_gap: 5000,
    extend_min: 0.9,
    pad: 1000,
};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DetectReport {
    pub genome: String,
    pub genome_bp: usize,
    pub genes_called: usize,
    pub rules: String,
    pub model_source: String,
    pub families: usize,
    pub params: ClusterParams,
    pub clusters: Vec<BgcCluster>,
    pub gene_call_ms: u128,
    pub domain_scan_ms: u128,
}

/// Detect biosynthetic gene clusters in a FASTA genome with the frozen V4 configuration
/// (released Pfam profiles for 55 of 56 families, Pfam GA cut-offs, V4 architecture rules).
pub fn run_bgc_detect(ws_root: &Path, fasta: &Path) -> Result<DetectReport, String> {
    let labs = if ws_root.ends_with("crvscience-labs") {
        ws_root.to_path_buf()
    } else {
        ws_root.join("crvscience-labs")
    };
    let seed_dir = labs.join("tools/gecco/pfam_seed");
    let official_dir = labs.join("tools/gecco/pfam_official");
    let models = build_domain_models_level(&seed_dir, Some(&official_dir), 2)?;
    let weights: HashMap<String, f32> =
        models.iter().map(|(f, _)| (f.accession.clone(), f.accessory_weight)).collect();
    let ga: HashMap<String, f32> = models.iter().map(|(f, _)| (f.accession.clone(), f.ga_bits)).collect();
    let n_official = models.iter().filter(|(f, _)| f.model_source == "official-pfam").count();

    let gp = call_genome_proteins(fasta)?;
    let (raw, scan_ms) = scan_proteins(&gp, &models);
    let genes = filter_hits(&raw, &ga, CutoffPolicy::Ga);
    let clusters = clusters_from_genes_with(&genes, &weights, RuleSet::V4, &FROZEN_V3_PARAMS);
    Ok(DetectReport {
        genome: gp.seq_id.clone(),
        genome_bp: gp.genome_bp,
        genes_called: gp.genes_called,
        rules: "V4".to_string(),
        model_source: format!("official-pfam ({} of {} families), remainder seed-built", n_official, models.len()),
        families: models.len(),
        params: FROZEN_V3_PARAMS,
        clusters,
        gene_call_ms: gp.gene_call_ms,
        domain_scan_ms: scan_ms,
    })
}

// ---------------------------------------------------------------------------------------------
// Calibration round 4: revision 4 rules/families, wider grid, second (independent) reserve set.
// ---------------------------------------------------------------------------------------------

fn restrict_to_families(genes: &[BgcGene], keep: &HashSet<String>) -> Vec<BgcGene> {
    genes
        .iter()
        .map(|g| {
            let mut c = g.clone();
            c.hits.retain(|h| keep.contains(&h.accession));
            c
        })
        .collect()
}

pub fn run_bgc_calibration_round4(ws_root: &Path) -> Result<Round3Report, String> {
    let t_start = Instant::now();
    let labs = if ws_root.ends_with("crvscience-labs") {
        ws_root.to_path_buf()
    } else {
        ws_root.join("crvscience-labs")
    };
    let seed_dir = labs.join("tools/gecco/pfam_seed");
    let official_dir = labs.join("tools/gecco/pfam_official");
    let mibig_dir = labs.join("tools/gecco/data/mibig_json_4.0");
    let gecco_exe = labs.join("tools/gecco/venv/Scripts/gecco.exe");
    let cache_dir = labs.join("runs/bgc_work/cache_v4");

    let models56 = build_domain_models_level(&seed_dir, Some(&official_dir), 2)?;
    let models34s = build_domain_models_level(&seed_dir, None, 0)?;
    let w56: HashMap<String, f32> = models56.iter().map(|(f, _)| (f.accession.clone(), f.accessory_weight)).collect();
    let ga56: HashMap<String, f32> = models56.iter().map(|(f, _)| (f.accession.clone(), f.ga_bits)).collect();
    let w34: HashMap<String, f32> = models34s.iter().map(|(f, _)| (f.accession.clone(), f.accessory_weight)).collect();
    let ga34: HashMap<String, f32> = models34s.iter().map(|(f, _)| (f.accession.clone(), f.ga_bits)).collect();
    let fam48: HashSet<String> = BGC_FAMILIES
        .iter()
        .chain(BGC_FAMILIES_EXTRA.iter())
        .map(|f| f.0.to_string())
        .collect();
    let w48: HashMap<String, f32> = w56.iter().filter(|(k, _)| fam48.contains(*k)).map(|(k, v)| (k.clone(), *v)).collect();
    let n_off = models56.iter().filter(|(f, _)| f.model_source == "official-pfam").count();
    eprintln!("[bgc-round4] {} families, released Pfam profiles for {}", models56.len(), n_off);

    // ---------------- development genomes (everything examined so far) ----------------
    let dev_specs: Vec<(&str, &str, &str, &str)> = vec![
        ("AL645882.2", "Streptomyces coelicolor A3(2)", "Prokaryota", "bgc-rich"),
        ("CP000560.1", "Bacillus velezensis FZB42", "Prokaryota", "bgc-rich"),
        ("AE004091.2", "Pseudomonas aeruginosa PAO1", "Prokaryota", "bgc-rich"),
        ("U00096.3", "Escherichia coli K-12 MG1655", "Prokaryota", "bgc-poor"),
        ("AP009493.1", "Streptomyces griseus NBRC 13350", "Prokaryota", "bgc-rich"),
        ("CP000667.1", "Salinispora tropica CNB-440", "Prokaryota", "bgc-rich"),
        ("BX470251.1", "Photorhabdus laumondii TTO1", "Prokaryota", "bgc-rich"),
        ("CP000113.1", "Myxococcus xanthus DK 1622", "Prokaryota", "bgc-rich"),
        ("CP000085.1", "Burkholderia thailandensis E264 chr II", "Prokaryota", "bgc-rich"),
        ("FN667742.1", "Xenorhabdus nematophila ATCC 19061", "Prokaryota", "bgc-rich"),
        ("FN554889.1", "Streptomyces scabiei 87.22", "Prokaryota", "bgc-rich"),
        ("BA000030.4", "Streptomyces avermitilis MA-4680", "Prokaryota", "bgc-rich"),
        ("AM420293.1", "Saccharopolyspora erythraea NRRL 2338", "Prokaryota", "bgc-rich"),
        ("CP000850.1", "Salinispora arenicola CNS-205", "Prokaryota", "bgc-rich"),
        ("BN001302.1", "Aspergillus nidulans FGSC A4 chr II", "Eukaryota", "out-of-regime"),
        ("CM000574.1", "Fusarium graminearum PH-1 chr 1", "Eukaryota", "out-of-regime"),
        ("NC_000908.2", "Mycoplasmoides genitalium G37", "Prokaryota", "negative-control"),
        ("NC_000915.1", "Helicobacter pylori 26695", "Prokaryota", "negative-control"),
        ("NC_000907.1", "Haemophilus influenzae Rd KW20", "Prokaryota", "negative-control"),
        ("NC_001136.10", "Saccharomyces cerevisiae S288C chr IV", "Eukaryota", "negative-control"),
        ("NC_002528.1", "Buchnera aphidicola APS", "Prokaryota", "negative-control"),
        ("NC_000117.1", "Chlamydia trachomatis D/UW-3", "Prokaryota", "negative-control"),
        ("NC_000919.1", "Treponema pallidum Nichols", "Prokaryota", "negative-control"),
        ("NC_003421.2", "Schizosaccharomyces pombe chr III", "Eukaryota", "negative-control"),
    ];
    let mut dev: Vec<GenomeCtx> = Vec::new();
    for s in &dev_specs {
        let (c, _) = build_ctx(s, &labs, &mibig_dir, &gecco_exe, &models56, &ga56, Some(&cache_dir), false)?;
        dev.push(c);
    }
    let dev_scope: Vec<&GenomeCtx> = dev
        .iter()
        .filter(|c| (c.kingdom == "Prokaryota" && c.role != "out-of-regime") || c.role == "negative-control")
        .collect();

    // ---------------- grid ----------------
    let default = ClusterParams::default();
    let mut grid: Vec<(RuleSet, ClusterParams)> = vec![(RuleSet::V3, FROZEN_V3_PARAMS), (RuleSet::V4, FROZEN_V3_PARAMS)];
    for rules in [RuleSet::V3, RuleSet::V4] {
        for sw in [4usize, 6, 8, 10] {
            for sm in [0.8f32, 1.0, 1.2, 1.5, 1.8] {
                for eg in [3000usize, 5000] {
                    for pad in [0usize, 1000, 2000] {
                        let prm = ClusterParams { support_window: sw, seed_min: sm, extend_gap: eg, extend_min: 0.9, pad };
                        if prm != FROZEN_V3_PARAMS {
                            grid.push((rules, prm));
                        }
                    }
                }
            }
        }
    }
    let grid_size = grid.len();
    eprintln!("[bgc-round4] development grid: {} configurations", grid_size);

    let run_cfg = |rules: RuleSet, prm: &ClusterParams, c: &GenomeCtx| -> Vec<BgcCluster> {
        match rules {
            RuleSet::V4 => clusters_from_genes_with(&c.genes_ext, &w56, rules, prm),
            RuleSet::V3 => clusters_from_genes_with(&restrict_to_families(&c.genes_ext, &fam48), &w48, rules, prm),
            _ => clusters_from_genes_with(&restrict_to_published(&c.genes_ext), &w34, rules, prm),
        }
    };

    const NEG_BUDGET: usize = 2;
    let mut scored: Vec<(R3Config, R3Score)> = Vec::new();
    for (rules, prm) in &grid {
        let label = format!("{:?} w={} seed>={} gap={} pad={}", rules, prm.support_window, prm.seed_min, prm.extend_gap, prm.pad);
        let sc = score_cfg(&label, &dev_scope, &|c| run_cfg(*rules, prm, c));
        scored.push((R3Config { label, rules: format!("{:?}", rules), params: *prm }, sc));
    }
    let mut order: Vec<usize> = (0..scored.len()).filter(|&i| scored[i].1.negative_clusters <= NEG_BUDGET).collect();
    if order.is_empty() {
        order = (0..scored.len()).collect();
    }
    order.sort_by(|&a, &b| {
        let (sa, sb) = (&scored[a].1, &scored[b].1);
        sb.f1
            .partial_cmp(&sa.f1)
            .unwrap()
            .then(sb.mean_iou.partial_cmp(&sa.mean_iou).unwrap())
            .then(sa.clusters.cmp(&sb.clusters))
            .then(a.cmp(&b))
    });
    // The frozen candidate must be a revision-4 configuration; revision 3 is kept as the comparison.
    let best_v4 = order.iter().copied().find(|&i| scored[i].0.rules == "V4").unwrap_or(order[0]);
    let frozen = scored[best_v4].0.clone();
    let dev_top: Vec<(R3Config, R3Score)> = order.iter().take(8).map(|&i| scored[i].clone()).collect();
    let frozen_rules = RuleSet::V4;

    let frozen_path = labs.join("runs/bgc-round4-frozen.json");
    fs::write(&frozen_path, serde_json::to_string_pretty(&frozen).unwrap_or_default()).map_err(|e| e.to_string())?;
    eprintln!("[bgc-round4] FROZEN on development set: {} (written to {})", frozen.label, frozen_path.display());

    let dev_baselines = vec![
        score_cfg("published #59: revision 2, 34 families", &dev_scope, &|c| run_cfg(RuleSet::V2, &default, c)),
        score_cfg("revision 3, parameters frozen in #60", &dev_scope, &|c| run_cfg(RuleSet::V3, &FROZEN_V3_PARAMS, c)),
        score_cfg("revision 4, parameters frozen in #60", &dev_scope, &|c| run_cfg(RuleSet::V4, &FROZEN_V3_PARAMS, c)),
        score_cfg("revision 4, newly frozen", &dev_scope, &|c| run_cfg(frozen_rules, &frozen.params, c)),
    ];
    let dev_gecco = score_gecco("GECCO v0.11.0 (-m 0.8)", &dev_scope);

    // ---------------- second reserve set (scanned only now) ----------------
    let reserve_specs: Vec<(&str, &str, &str, &str)> = vec![
        ("CP000076.1", "Pseudomonas protegens Pf-5", "Prokaryota", "bgc-rich"),
        ("CP000075.1", "Pseudomonas syringae B728a", "Prokaryota", "bgc-rich"),
        ("CP001037.1", "Nostoc punctiforme PCC 73102", "Prokaryota", "bgc-rich"),
        ("CP005094.1", "Azotobacter vinelandii CA", "Prokaryota", "bgc-rich"),
        ("AL123456.3", "Mycobacterium tuberculosis H37Rv", "Prokaryota", "bgc-rich"),
        ("FR845719.1", "Streptomyces venezuelae ATCC 10712", "Prokaryota", "bgc-rich"),
        ("CP007155.1", "Kutzneria albida DSM 43870", "Prokaryota", "bgc-rich"),
        ("CM000913.1", "Streptomyces clavuligerus ATCC 27064", "Prokaryota", "bgc-rich"),
        ("CM000172.1", "Aspergillus fumigatus Af293 chr 4", "Eukaryota", "out-of-regime"),
        ("NC_000963.1", "Rickettsia prowazekii Madrid E", "Prokaryota", "negative-control"),
        ("NC_001318.1", "Borreliella burgdorferi B31", "Prokaryota", "negative-control"),
        ("NC_000909.1", "Methanocaldococcus jannaschii (archaeon)", "Prokaryota", "negative-control"),
    ];
    let mut reserve: Vec<GenomeCtx> = Vec::new();
    let mut res_seed34: HashMap<String, Vec<BgcGene>> = HashMap::new();
    for s in &reserve_specs {
        let (c, gp_opt) = build_ctx(s, &labs, &mibig_dir, &gecco_exe, &models56, &ga56, Some(&cache_dir), true)?;
        let gp = match gp_opt {
            Some(gp) => gp,
            None => call_genome_proteins(&labs.join("cache/bgc").join(format!("{}.fasta", s.0)))?,
        };
        let (raw34, _) = scan_proteins(&gp, &models34s);
        res_seed34.insert(s.0.to_string(), filter_hits(&raw34, &ga34, CutoffPolicy::Ga));
        reserve.push(c);
    }
    let reserve_scope: Vec<&GenomeCtx> = reserve.iter().filter(|c| c.role != "out-of-regime").collect();
    let reserve_all: Vec<&GenomeCtx> = reserve.iter().collect();

    let res_pub = |c: &GenomeCtx| clusters_from_genes_with(&res_seed34[&c.acc], &w34, RuleSet::V2, &default);
    let res_r3 = |c: &GenomeCtx| clusters_from_genes_with(&restrict_to_families(&c.genes_ext, &fam48), &w48, RuleSet::V3, &FROZEN_V3_PARAMS);
    let res_r4_old = |c: &GenomeCtx| clusters_from_genes_with(&c.genes_ext, &w56, RuleSet::V4, &FROZEN_V3_PARAMS);
    let res_r4 = |c: &GenomeCtx| clusters_from_genes_with(&c.genes_ext, &w56, frozen_rules, &frozen.params);
    let reserve_scores = vec![
        score_cfg("published #59 (revision 2, seed profiles)", &reserve_scope, &res_pub),
        score_cfg("revision 3 as frozen in #60", &reserve_scope, &res_r3),
        score_cfg("revision 4, #60 parameters", &reserve_scope, &res_r4_old),
        score_cfg("revision 4, newly frozen", &reserve_scope, &res_r4),
    ];
    let reserve_gecco = score_gecco("GECCO v0.11.0 (-m 0.8)", &reserve_scope);

    let genome_results = |ctxs: &[GenomeCtx], extra: &[(&str, &dyn Fn(&GenomeCtx) -> Vec<BgcCluster>)]| -> Vec<R3GenomeResult> {
        ctxs.iter()
            .map(|c| {
                let mut configs = Vec::new();
                for (label, f) in extra {
                    let cl = f(c);
                    let iv: Vec<Interval> = cl
                        .iter()
                        .map(|x| Interval { seq_id: x.seq_id.clone(), start: x.start, end: x.end, label: x.classes.join(";") })
                        .collect();
                    let (a, h, mi, _) = locus_metrics(&c.acc, &c.seq_id, &c.truth, &iv);
                    let ov = iv.iter().filter(|p| c.truth.iter().any(|t| overlap_bp(p.start, p.end, t.1, t.2) > 0)).count();
                    configs.push((label.to_string(), a, h, iv.len(), ov, mi));
                }
                let (ga, gh, gmi, _) = locus_metrics(&c.acc, &c.seq_id, &c.truth, &c.gecco_iv);
                let gov = c.gecco_iv.iter().filter(|p| c.truth.iter().any(|t| overlap_bp(p.start, p.end, t.1, t.2) > 0)).count();
                R3GenomeResult {
                    genome_id: c.acc.clone(),
                    organism: c.organism.clone(),
                    kingdom: c.kingdom.clone(),
                    role: c.role.clone(),
                    genome_bp: c.genome_bp,
                    mibig_loci: c.truth.len(),
                    configs,
                    gecco: (ga, gh, c.gecco_iv.len(), gov, gmi),
                }
            })
            .collect()
    };
    let d_r3 = |c: &GenomeCtx| run_cfg(RuleSet::V3, &FROZEN_V3_PARAMS, c);
    let d_r4 = |c: &GenomeCtx| run_cfg(frozen_rules, &frozen.params, c);
    let dev_genomes = genome_results(&dev, &[("revision 3 (#60)", &d_r3), ("revision 4", &d_r4)]);
    let reserve_genomes = genome_results(
        &reserve,
        &[("published #59", &res_pub), ("revision 3 (#60)", &res_r3), ("revision 4 old params", &res_r4_old), ("revision 4", &res_r4)],
    );
    let (mut rl_f, mut rl_g, mut dl_f) = (Vec::new(), Vec::new(), Vec::new());
    for c in &reserve_all {
        let iv: Vec<Interval> = res_r4(c)
            .iter()
            .map(|x| Interval { seq_id: x.seq_id.clone(), start: x.start, end: x.end, label: x.classes.join(";") })
            .collect();
        rl_f.extend(locus_metrics(&c.acc, &c.seq_id, &c.truth, &iv).3);
        rl_g.extend(locus_metrics(&c.acc, &c.seq_id, &c.truth, &c.gecco_iv).3);
    }
    for c in dev.iter().filter(|c| c.role != "negative-control") {
        let iv: Vec<Interval> = d_r4(c)
            .iter()
            .map(|x| Interval { seq_id: x.seq_id.clone(), start: x.start, end: x.end, label: x.classes.join(";") })
            .collect();
        dl_f.extend(locus_metrics(&c.acc, &c.seq_id, &c.truth, &iv).3);
    }

    Ok(Round3Report {
        pipeline_name: "bgc-round4 (revision 4 rules/families, wider grid, second reserve set)".to_string(),
        protocol: "Development = the 24 genomes examined in rounds 1-3. A grid of rule revisions and cluster parameters is scored on the \
                   development scope only; the best revision-4 configuration is written to disk; the 12 second-reserve genomes are \
                   gene-called and scored afterwards."
            .to_string(),
        grid_size,
        objective: format!(
            "<= {} cluster(s) across the development negative controls; then maximise the balanced score of (locus recall, MIBiG-overlap share); \
             then mean best-IoU; then fewer clusters. Revision 3 with the #60 parameters is scored alongside as the comparison.",
            NEG_BUDGET
        ),
        negative_budget: NEG_BUDGET,
        frozen,
        dev_top_configs: dev_top,
        dev_baselines,
        dev_gecco,
        reserve_scores,
        reserve_gecco,
        reserve_genomes,
        dev_genomes,
        reserve_loci_frozen: rl_f,
        reserve_loci_gecco: rl_g,
        dev_loci_frozen: dl_f,
        families_added: BGC_FAMILIES_EXTRA2.iter().map(|f| format!("{} {}", f.0, f.1)).collect(),
        total_duration_ms: t_start.elapsed().as_millis(),
    })
}
