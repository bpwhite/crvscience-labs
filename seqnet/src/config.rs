//! Configuration, locus panel definitions, and runtime targets for `bio-suite-seqnet`.

use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// Execution mode for the 3-stage sequence analysis neural network.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExecutionMode {
    /// Pure terminal execution with structured epoch progress tables.
    Headless,
    /// Dual console telemetry + live `eframe`/`egui` desktop monitor window.
    GuiMonitor,
}

/// High-level phylogenetic or population-genetic study design.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum AnalysisLevel {
    SingleLocusBarcode,
    MultiLocusPhylogenetic,
    PhylogenomicSupermatrix,
    PopulationCoalescent,
}

impl AnalysisLevel {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::SingleLocusBarcode => "single-locus-barcode",
            Self::MultiLocusPhylogenetic => "multi-locus-phylogenetic",
            Self::PhylogenomicSupermatrix => "phylogenomic-supermatrix",
            Self::PopulationCoalescent => "population-coalescent",
        }
    }
}

/// Taxonomic rank at which OTU sampling rules are organized.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum OtuRank {
    Population,
    Species,
    Genus,
    Family,
}

impl OtuRank {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Population => "population",
            Self::Species => "species",
            Self::Genus => "genus",
            Self::Family => "family",
        }
    }
}

/// Genomic compartment of a target locus, governing genetic code & evolutionary rate priors.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum LocusClass {
    Nuclear,
    Mitochondrial,
}

impl LocusClass {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Nuclear => "nuclear",
            Self::Mitochondrial => "mitochondrial",
        }
    }
}

/// Canonical marker specification used by Input, Processing, and Output neuron clusters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LocusSpec {
    pub canonical_name: String,
    pub locus_class: LocusClass,
    pub protein_coding: bool,
    /// NCBI genetic code table (1 = Standard Nuclear, 5 = Invertebrate Mitochondrial).
    pub genetic_code_table: u8,
    pub min_expected_bp: usize,
    pub max_expected_bp: usize,
    /// Synonyms recognized in NCBI Entrez queries and GenBank feature annotations.
    pub entrez_synonyms: Vec<String>,
}

/// Convergence and quality-control hyperparameters for the 3-stage neuron schedule.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConsensusConfig {
    /// Target organism clade (e.g., `"Coleoptera"`).
    pub taxon_group: String,
    /// Analysis design (e.g., `MultiLocusPhylogenetic`).
    pub analysis_level: AnalysisLevel,
    /// Total target loci (`target_nuclear_loci + target_mito_loci`).
    pub target_loci_count: usize,
    pub target_nuclear_loci: usize,
    pub target_mito_loci: usize,
    /// Target count of final retained OTU terminals (e.g., `100`).
    pub target_otu_count: usize,
    /// Organizational OTU rank (e.g., `Genus` with separate species terminals per genus).
    pub otu_rank: OtuRank,
    /// Minimum distinct species terminals required to retain a Genus (e.g., `2`).
    pub min_species_per_genus: usize,
    /// Minimum fraction of target loci every retained OTU terminal must possess (default `0.90`).
    pub min_locus_coverage_ratio: f64,
    /// Minimum neural consensus score (`[0.0, 1.0]`) to pass Stage 2 Processing without quarantine.
    pub neural_pass_threshold: f64,
    /// Global consensus score required to declare an "Optimal Data Set" and halt early.
    pub optimal_consensus_target: f64,
    /// Maximum feedback schedule epochs before exporting best available consensus set.
    pub max_epochs: usize,
    /// Maximum ambiguous nucleotide fraction (`N`, `?`, IUPAC degenerate bases) tolerated.
    pub max_ambig_fraction: f64,
    /// Maximum robust z-score on within-locus K2P distance before flagging paralog/swap.
    pub max_within_locus_zscore: f64,
    /// SQLite database file path.
    pub sqlite_path: PathBuf,
    /// Explicit logical `database_id` stored inside `database_registry`.
    pub database_id: String,
    /// Destination path for clean partitioned NEXUS supermatrix.
    pub nexus_output_path: PathBuf,
    /// Optional NCBI Entrez API key (`10 req/s` when set, `3 req/s` default).
    pub ncbi_api_key: Option<String>,
    /// Contact email passed to NCBI E-utilities per usage guidelines.
    pub ncbi_email: String,
    /// Whether to perform live NCBI Entrez delta queries (`false` = cache/offline-only mode).
    pub allow_live_entrez: bool,
    /// Runtime UI mode (`Headless` vs `GuiMonitor`).
    pub execution_mode: ExecutionMode,
    /// Explicit locus panel by canonical name (e.g. `["SSU_18S", "LSU_28S", "COI", "rrnL_16S"]`).
    /// When empty, the first `target_nuclear_loci` + `target_mito_loci` panel entries are used.
    pub locus_names: Vec<String>,
    /// Stratified Stage 1 sampling: one Entrez query per `(stratum, locus)` (e.g. beetle families),
    /// so sampled species overlap across loci. When empty, `taxon_group` is queried as one stratum.
    pub sampling_taxa: Vec<String>,
    /// Entrez `retmax` per `(stratum, locus)` query on a cold start (no watermark yet).
    pub retmax_per_stratum: usize,
    /// Stage 3 -> Stage 1 feedback limits: species directives and missing loci refilled per epoch.
    pub max_refill_directives: usize,
    pub max_refill_loci_per_directive: usize,
    /// Minimum locus coverage for a species to receive a gap-refill directive (default `0.55`).
    pub refill_min_coverage: f64,
    /// Region coherence gate: a record must cover this fraction of the locus core (default `0.6`;
    /// `0` disables the gate). See `region.rs`.
    pub region_min_core_coverage: f64,
    /// Fraction of mapped records that must cover an anchor position for it to be in the core.
    pub region_core_support: f64,
    /// Mitogenome mode: mitochondrial loci are sliced from complete/near-complete mitogenome
    /// records (CDS/rRNA features) instead of per-locus title queries, so every mito locus of a
    /// species comes from one specimen. Nuclear loci are still queried per locus.
    pub mitogenome_loci: bool,
    /// Entrez `retmax` for the per-stratum mitogenome query (mitogenome mode only).
    pub mitogenome_retmax_per_stratum: usize,
}

impl ConsensusConfig {
    /// The locus panel this configuration runs on.
    pub fn locus_panel(&self) -> Vec<LocusSpec> {
        if self.locus_names.is_empty() {
            build_locus_panel(self.target_nuclear_loci, self.target_mito_loci, true)
        } else {
            build_named_locus_panel(&self.locus_names, true)
        }
    }

    /// Entrez query strata (`sampling_taxa`, or the whole `taxon_group`).
    pub fn sampling_strata(&self) -> Vec<String> {
        if self.sampling_taxa.is_empty() {
            vec![self.taxon_group.clone()]
        } else {
            self.sampling_taxa.clone()
        }
    }
}

/// Builds a panel from canonical locus names (see `build_locus_panel` for the catalogue).
/// Unknown names are skipped; panel order follows `names`.
pub fn build_named_locus_panel(names: &[String], invert_mito: bool) -> Vec<LocusSpec> {
    let catalogue = build_locus_panel(25, 15, invert_mito);
    names
        .iter()
        .filter_map(|n| catalogue.iter().find(|l| l.canonical_name.eq_ignore_ascii_case(n)).cloned())
        .collect()
}

impl Default for ConsensusConfig {
    fn default() -> Self {
        Self {
            taxon_group: "Coleoptera".to_string(),
            analysis_level: AnalysisLevel::MultiLocusPhylogenetic,
            target_loci_count: 25,
            target_nuclear_loci: 18,
            target_mito_loci: 7,
            target_otu_count: 100,
            otu_rank: OtuRank::Genus,
            min_species_per_genus: 2,
            min_locus_coverage_ratio: 0.90,
            neural_pass_threshold: 0.62,
            optimal_consensus_target: 0.88,
            max_epochs: 8,
            max_ambig_fraction: 0.05,
            max_within_locus_zscore: 3.2,
            sqlite_path: PathBuf::from("coleoptera_seqnet_cache.sqlite"),
            database_id: "CRV-COLEOPTERA-ML25-DB01".to_string(),
            nexus_output_path: PathBuf::from("coleoptera_25loci_optimal.nex"),
            ncbi_api_key: std::env::var("NCBI_API_KEY").ok(),
            ncbi_email: std::env::var("NCBI_EMAIL")
                .unwrap_or_else(|_| "labs@crvscience.com".to_string()),
            allow_live_entrez: true,
            execution_mode: ExecutionMode::Headless,
            locus_names: Vec::new(),
            sampling_taxa: Vec::new(),
            retmax_per_stratum: 400,
            max_refill_directives: 6,
            max_refill_loci_per_directive: 4,
            refill_min_coverage: 0.55,
            region_min_core_coverage: 0.6,
            region_core_support: 0.5,
            mitogenome_loci: false,
            mitogenome_retmax_per_stratum: 60,
        }
    }
}

/// Builds a curated panel of `n_nuclear` nuclear + `n_mito` mitochondrial phylogenetic markers
/// widely sequenced across Coleoptera & Metazoa (McKenna et al. / Che et al. / standard insect panels),
/// scaling cleanly up to 100+ single-copy nuclear orthologs (USCO/BUSCO) when requested.
pub fn build_locus_panel(n_nuclear: usize, n_mito: usize, invert_mito: bool) -> Vec<LocusSpec> {
    let mito_code = if invert_mito { 5 } else { 2 };

    let core_mito: &[(&str, bool, usize, usize, &[&str])] = &[
        ("COI", true, 500, 1600, &["COI", "COX1", "CO1", "cytochrome c oxidase subunit I", "cytochrome oxidase subunit I"]),
        ("COII", true, 450, 800, &["COII", "COX2", "CO2", "cytochrome c oxidase subunit II"]),
        ("COIII", true, 500, 850, &["COIII", "COX3", "CO3", "cytochrome c oxidase subunit III"]),
        ("CYTB", true, 600, 1200, &["CYTB", "COB", "cytochrome b"]),
        ("ND1", true, 600, 1050, &["ND1", "NAD1", "NADH dehydrogenase subunit 1"]),
        ("ND2", true, 600, 1100, &["ND2", "NAD2", "NADH dehydrogenase subunit 2"]),
        ("ND4", true, 750, 1450, &["ND4", "NAD4", "NADH dehydrogenase subunit 4"]),
        ("ND5", true, 900, 1800, &["ND5", "NAD5", "NADH dehydrogenase subunit 5"]),
        ("ATP6", true, 450, 750, &["ATP6", "ATP synthase F0 subunit 6"]),
        ("rrnL_16S", false, 400, 1450, &["16S", "rrnL", "l-rRNA", "16S ribosomal RNA"]),
        ("rrnS_12S", false, 300, 950, &["12S", "rrnS", "s-rRNA", "12S ribosomal RNA"]),
        ("ND3", true, 280, 420, &["ND3", "NAD3", "NADH dehydrogenase subunit 3"]),
        ("ND4L", true, 220, 350, &["ND4L", "NAD4L"]),
        ("ND6", true, 350, 600, &["ND6", "NAD6"]),
        ("ATP8", true, 120, 220, &["ATP8", "ATP synthase F0 subunit 8"]),
    ];

    let core_nuclear: &[(&str, bool, usize, usize, &[&str])] = &[
        // "SSU"/"LSU" alone also name mitochondrial 12S/16S rRNA, so they are not used as synonyms
        ("SSU_18S", false, 600, 2100, &["18S", "18S ribosomal RNA", "18S rRNA"]),
        ("LSU_28S", false, 600, 3800, &["28S", "28S ribosomal RNA", "28S rRNA"]),
        ("H3", true, 280, 400, &["H3", "histone H3"]),
        ("EF1a", true, 500, 1300, &["EF-1a", "EF1a", "elongation factor 1-alpha"]),
        ("CAD", true, 600, 2600, &["CAD", "carbamoyl-phosphate synthetase"]),
        ("Wg", true, 350, 600, &["wingless", "Wg", "Wnt1"]),
        ("ArgK", true, 500, 1100, &["ArgK", "arginine kinase"]),
        ("RNApolII", true, 600, 2200, &["RPB1", "RNA polymerase II largest subunit"]),
        ("RPB2", true, 600, 2100, &["RPB2", "RNA polymerase II second largest subunit"]),
        ("Topo", true, 550, 1500, &["topoisomerase I", "Topo1"]),
        ("MPSR1", true, 500, 1400, &["MPSR1", "muscle protein 20"]),
        ("NaK", true, 600, 1600, &["NaK", "sodium potassium ATPase alpha subunit"]),
        ("PEPCK", true, 500, 1500, &["PEPCK", "phosphoenolpyruvate carboxykinase"]),
        ("PGD", true, 450, 1300, &["PGD", "6-phosphogluconate dehydrogenase"]),
        ("TPI", true, 400, 900, &["TPI", "triosephosphate isomerase"]),
        ("GAPDH", true, 500, 1100, &["GAPDH", "glyceraldehyde-3-phosphate dehydrogenase"]),
        ("IDH", true, 500, 1200, &["IDH", "isocitrate dehydrogenase"]),
        ("AAT", true, 500, 1300, &["AAT", "GOT1", "aspartate aminotransferase"]),
        ("DDC", true, 500, 1400, &["DDC", "dopa decarboxylase"]),
        ("Enolase", true, 500, 1300, &["Eno", "enolase"]),
        ("ITS2", false, 250, 800, &["ITS2", "internal transcribed spacer 2"]),
        ("H4", true, 220, 350, &["H4", "histone H4"]),
        ("Actin5C", true, 500, 1200, &["Act5C", "actin"]),
        ("AlphaTub", true, 500, 1400, &["alpha-tubulin", "TUBA"]),
        ("Hsp70", true, 600, 1900, &["Hsp70", "heat shock protein 70"]),
    ];

    let mut panel = Vec::with_capacity(n_nuclear + n_mito);

    for idx in 0..n_nuclear {
        if idx < core_nuclear.len() {
            let (name, coding, min_bp, max_bp, syns) = core_nuclear[idx];
            panel.push(LocusSpec {
                canonical_name: name.to_string(),
                locus_class: LocusClass::Nuclear,
                protein_coding: coding,
                genetic_code_table: 1,
                min_expected_bp: min_bp,
                max_expected_bp: max_bp,
                entrez_synonyms: syns.iter().map(|s| (*s).to_string()).collect(),
            });
        } else {
            let usco_id = format!("USCO_NUC_{:03}", idx + 1);
            panel.push(LocusSpec {
                canonical_name: usco_id.clone(),
                locus_class: LocusClass::Nuclear,
                protein_coding: true,
                genetic_code_table: 1,
                min_expected_bp: 450,
                max_expected_bp: 1800,
                entrez_synonyms: vec![usco_id],
            });
        }
    }

    for idx in 0..n_mito {
        if idx < core_mito.len() {
            let (name, coding, min_bp, max_bp, syns) = core_mito[idx];
            panel.push(LocusSpec {
                canonical_name: name.to_string(),
                locus_class: LocusClass::Mitochondrial,
                protein_coding: coding,
                genetic_code_table: mito_code,
                min_expected_bp: min_bp,
                max_expected_bp: max_bp,
                entrez_synonyms: syns.iter().map(|s| (*s).to_string()).collect(),
            });
        } else {
            let mito_id = format!("MITO_LOC_{:02}", idx + 1);
            panel.push(LocusSpec {
                canonical_name: mito_id.clone(),
                locus_class: LocusClass::Mitochondrial,
                protein_coding: true,
                genetic_code_table: mito_code,
                min_expected_bp: 300,
                max_expected_bp: 1500,
                entrez_synonyms: vec![mito_id],
            });
        }
    }

    panel
}
