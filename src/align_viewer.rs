//! `align_viewer` — Rust-Native Interactive Sequence Alignment Viewer & Publication Snapshot Engine
//!
//! Modeled after AliView (Larsson 2014, https://ormbunkar.se/aliview/) and purpose-built for the
//! CRV Science Labs `bio-suite` & RuNIE (Rust Native Inference Engine) phylogenomic ecosystem.
//!
//! Capabilities:
//! 1. **Interactive Native Desktop GUI (`bio-suite align-view` / `runie-align-view`)**:
//!    - Native top Menu Bar (`File`, `Edit`, `View`, `Translation`, `Features & Annotations`,
//!      `RuNIE Pipeline`, `Snapshots & Export`, `Help`)
//!    - Interactive Ribbon Toolbar, Global Bird's-Eye Minimap with draggable viewport lens,
//!      Feature/Annotation Sidebar (click any gene, CDS, catalytic motif, or complexity region
//!      to zoom/select), Conservation + Shannon Entropy track, Consensus row, and live hover telemetry.
//!    - Smooth zoom from 1 px/column (whole-mitogenome / multi-locus overview) up to 24 px/column
//!      residue-level inspection.
//! 2. **Multi-Alphabet & Genetic-Code Translation Engine**:
//!    - DNA / RNA Nucleotide view
//!    - Amino Acid Translated Overlay (3-base codon span)
//!    - Amino Acid Compact view (1 column per codon)
//!    - Dual DNA + Amino Acid split-row view
//!    - RuNIE `RY3` Purine/Pyrimidine recoding (`R`/`Y`/`-`)
//!    - RuNIE `Dayhoff-6` amino-acid functional group recoding (`1`..`6`)
//!    - Supports Reading Frames +1, +2, +3, Reverse Complement, and NCBI Genetic Codes
//!      1 (Standard), 2 (Vertebrate Mitochondrial), 5 (Invertebrate Mitochondrial), and 11 (Bacterial/Plastid).
//! 3. **Automated Pipeline & Manuscript Snapshot Engine (`--snapshot-svg`, `--snapshot-png`, `--snapshot-json`)**:
//!    - Zooms to any coordinate range (`--region start:end`), named feature (`--feature <name>`),
//!      or automatically selects the highest-information focal motif window (`--auto-focal`)
//!      and exports pixel-locked publication SVGs (CRV Science Labs Figure 2 standard) and
//!      high-resolution raster PNG snapshots.

use crate::fetch::{FastaRecord, GenomicAnnotation, PoliteCacheManager};
use crate::phylo::{align_sequences_with_mode, AlignerMode};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

/// Sequence alphabet / display mode in the alignment viewer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ViewAlphabetMode {
    /// Standard 4-base DNA/RNA nucleotide view (`A`, `C`, `G`, `T`/`U`, `-`, `N`)
    Nucleotide,
    /// Translated Amino Acid 1-letter code spanned across 3-nucleotide codon columns
    AminoAcidOverlay,
    /// Compact Amino Acid view (1 column = 1 translated codon)
    AminoAcidCompact,
    /// Dual row per taxon: DNA nucleotide on top + translated Amino Acid below
    DualDnaAndAminoAcid,
    /// RuNIE RY3 Purine (`R`) / Pyrimidine (`Y`) / Gap (`-`) recoding
    Ry3Recoded,
    /// RuNIE Dayhoff-6 6-state amino-acid functional group recoding (`1`..`6`)
    Dayhoff6Recoded,
}

impl ViewAlphabetMode {
    pub fn from_cli_arg(s: &str) -> Self {
        match s.trim().to_ascii_lowercase().as_str() {
            "aa" | "protein" | "amino" | "amino-overlay" => ViewAlphabetMode::AminoAcidOverlay,
            "aa-compact" | "compact-aa" | "codon" | "protein-compact" => {
                ViewAlphabetMode::AminoAcidCompact
            }
            "dual" | "both" | "dna+aa" => ViewAlphabetMode::DualDnaAndAminoAcid,
            "ry3" | "ry" | "purine-pyrimidine" => ViewAlphabetMode::Ry3Recoded,
            "dayhoff6" | "dayhoff" | "d6" => ViewAlphabetMode::Dayhoff6Recoded,
            _ => ViewAlphabetMode::Nucleotide,
        }
    }

    pub fn label(&self) -> &'static str {
        match self {
            ViewAlphabetMode::Nucleotide => "DNA / Nucleotide",
            ViewAlphabetMode::AminoAcidOverlay => "AA Translation (3-bp Codon)",
            ViewAlphabetMode::AminoAcidCompact => "AA Compact (1 Col/Codon)",
            ViewAlphabetMode::DualDnaAndAminoAcid => "Dual DNA + AA Track",
            ViewAlphabetMode::Ry3Recoded => "RuNIE RY3 (R/Y/-)",
            ViewAlphabetMode::Dayhoff6Recoded => "RuNIE Dayhoff-6",
        }
    }
}

/// Color palette for nucleotide and amino-acid residues.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AlignmentColorScheme {
    /// AliView default high-contrast nucleotide & Taylor/ClustalX amino-acid palette
    AliViewDefault,
    /// ClustalX chemical group coloring
    ClustalX,
    /// Kyte-Doolittle hydrophobicity gradient (warm red hydrophobic -> cool blue hydrophilic)
    Hydrophobicity,
    /// RuNIE Conservation & Parsimony-Informative site highlighter
    RunieConservation,
    /// Okabe-Ito color-blind accessible palette
    OkabeItoSafe,
}

impl AlignmentColorScheme {
    pub fn from_cli_arg(s: &str) -> Self {
        match s.trim().to_ascii_lowercase().as_str() {
            "clustal" | "clustalx" => AlignmentColorScheme::ClustalX,
            "hydro" | "hydrophobicity" | "kd" => AlignmentColorScheme::Hydrophobicity,
            "runie" | "conservation" | "pi" => AlignmentColorScheme::RunieConservation,
            "okabe" | "okabe-ito" | "cb" | "colorblind" => AlignmentColorScheme::OkabeItoSafe,
            _ => AlignmentColorScheme::AliViewDefault,
        }
    }

    pub fn label(&self) -> &'static str {
        match self {
            AlignmentColorScheme::AliViewDefault => "AliView Default",
            AlignmentColorScheme::ClustalX => "ClustalX Chemistry",
            AlignmentColorScheme::Hydrophobicity => "Kyte-Doolittle Hydrophobicity",
            AlignmentColorScheme::RunieConservation => "RuNIE Conservation / PI",
            AlignmentColorScheme::OkabeItoSafe => "Okabe-Ito Color-Blind Safe",
        }
    }
}

/// Character display mode for residues matching the consensus or first reference taxon.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConsensusMatchMode {
    /// Show every residue letter explicitly
    ShowAllCharacters,
    /// Show `.` for residues identical to the column majority consensus
    DotsForConsensus,
    /// Show `.` for residues identical to Taxon #1 (Reference sequence)
    DotsForReferenceTaxon,
    /// Dim invariant columns and highlight only variable / parsimony-informative sites
    HighlightVariableOnly,
}

impl ConsensusMatchMode {
    pub fn label(&self) -> &'static str {
        match self {
            ConsensusMatchMode::ShowAllCharacters => "All Residues",
            ConsensusMatchMode::DotsForConsensus => "Dots (.) vs Consensus",
            ConsensusMatchMode::DotsForReferenceTaxon => "Dots (.) vs Ref Taxon #1",
            ConsensusMatchMode::HighlightVariableOnly => "Highlight Variable / PI Only",
        }
    }
}

/// Supported NCBI Genetic Code translation tables.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GeneticCodeTable {
    /// NCBI Table 1: Standard Nuclear Code
    Standard1,
    /// NCBI Table 2: Vertebrate Mitochondrial Code (AGA/AGG=*, ATA=M, TGA=W)
    VertebrateMito2,
    /// NCBI Table 5: Invertebrate Mitochondrial Code (AGA/AGG=S, ATA=M, TGA=W)
    InvertebrateMito5,
    /// NCBI Table 11: Bacterial, Archaeal & Plant Plastid Code
    BacterialPlastid11,
}

impl GeneticCodeTable {
    pub fn from_id(id: u8) -> Self {
        match id {
            2 => GeneticCodeTable::VertebrateMito2,
            5 => GeneticCodeTable::InvertebrateMito5,
            11 => GeneticCodeTable::BacterialPlastid11,
            _ => GeneticCodeTable::Standard1,
        }
    }

    pub fn table_id(&self) -> u8 {
        match self {
            GeneticCodeTable::Standard1 => 1,
            GeneticCodeTable::VertebrateMito2 => 2,
            GeneticCodeTable::InvertebrateMito5 => 5,
            GeneticCodeTable::BacterialPlastid11 => 11,
        }
    }

    pub fn label(&self) -> &'static str {
        match self {
            GeneticCodeTable::Standard1 => "Table 1: Standard Nuclear",
            GeneticCodeTable::VertebrateMito2 => "Table 2: Vertebrate Mitochondrial",
            GeneticCodeTable::InvertebrateMito5 => "Table 5: Invertebrate Mitochondrial",
            GeneticCodeTable::BacterialPlastid11 => "Table 11: Bacterial / Plastid",
        }
    }

    /// Translate a 3-nucleotide codon (`[u8; 3]`, uppercase ASCII) into a 1-letter amino acid.
    pub fn translate_codon(&self, c1: u8, c2: u8, c3: u8) -> u8 {
        if c1 == b'-' && c2 == b'-' && c3 == b'-' {
            return b'-';
        }
        if c1 == b'-' || c2 == b'-' || c3 == b'-' {
            return b'X';
        }
        let u1 = if c1 == b'U' { b'T' } else { c1 };
        let u2 = if c2 == b'U' { b'T' } else { c2 };
        let u3 = if c3 == b'U' { b'T' } else { c3 };
        let codon = [u1, u2, u3];

        // Table-specific overrides over Standard Table 1
        match self {
            GeneticCodeTable::VertebrateMito2 => match &codon {
                b"AGA" | b"AGG" => return b'*',
                b"ATA" => return b'M',
                b"TGA" => return b'W',
                _ => {}
            },
            GeneticCodeTable::InvertebrateMito5 => match &codon {
                b"AGA" | b"AGG" => return b'S',
                b"ATA" => return b'M',
                b"TGA" => return b'W',
                _ => {}
            },
            GeneticCodeTable::Standard1 | GeneticCodeTable::BacterialPlastid11 => {}
        }

        match &codon {
            b"TTT" | b"TTC" => b'F',
            b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => b'L',
            b"ATT" | b"ATC" | b"ATA" => b'I',
            b"ATG" => b'M',
            b"GTT" | b"GTC" | b"GTA" | b"GTG" => b'V',
            b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => b'S',
            b"CCT" | b"CCC" | b"CCA" | b"CCG" => b'P',
            b"ACT" | b"ACC" | b"ACA" | b"ACG" => b'T',
            b"GCT" | b"GCC" | b"GCA" | b"GCG" => b'A',
            b"TAT" | b"TAC" => b'Y',
            b"TAA" | b"TAG" | b"TGA" => b'*',
            b"CAT" | b"CAC" => b'H',
            b"CAA" | b"CAG" => b'Q',
            b"AAT" | b"AAC" => b'N',
            b"AAA" | b"AAG" => b'K',
            b"GAT" | b"GAC" => b'D',
            b"GAA" | b"GAG" => b'E',
            b"TGT" | b"TGC" => b'C',
            b"TGG" => b'W',
            b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => b'R',
            b"GGT" | b"GGC" | b"GGA" | b"GGG" => b'G',
            _ => b'X',
        }
    }
}

/// Convert a nucleotide to RuNIE `RY3` Purine/Pyrimidine code (`R`, `Y`, `-`, `N`).
pub fn recode_ry3(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' | b'G' => b'R',
        b'C' | b'T' | b'U' => b'Y',
        b'-' | b'.' => b'-',
        _ => b'N',
    }
}

/// Convert an amino acid 1-letter residue to RuNIE `Dayhoff-6` functional state (`1`..`6`)
/// and its group label:
/// `1`: AGPST (Small/neutral), `2`: DENQ (Acidic/amide), `3`: HKR (Basic),
/// `4`: ILMV (Hydrophobic aliphatic), `5`: FYW (Aromatic), `6`: C (Cysteine).
pub fn recode_dayhoff6(aa: u8) -> u8 {
    match aa.to_ascii_uppercase() {
        b'A' | b'G' | b'P' | b'S' | b'T' => b'1',
        b'D' | b'E' | b'N' | b'Q' => b'2',
        b'H' | b'K' | b'R' => b'3',
        b'I' | b'L' | b'M' | b'V' => b'4',
        b'F' | b'Y' | b'W' => b'5',
        b'C' => b'6',
        b'-' => b'-',
        b'*' => b'*',
        _ => b'X',
    }
}

/// Amino acid chemical property group name for status-bar inspection.
pub fn amino_acid_full_name_and_chemistry(aa: u8) -> (&'static str, &'static str, f64) {
    match aa.to_ascii_uppercase() {
        b'A' => ("Alanine (Ala)", "Hydrophobic Aliphatic (Dayhoff-1)", 1.8),
        b'R' => ("Arginine (Arg)", "Positively Charged Basic (Dayhoff-3)", -4.5),
        b'N' => ("Asparagine (Asn)", "Polar Uncharged Amide (Dayhoff-2)", -3.5),
        b'D' => ("Aspartate (Asp)", "Negatively Charged Acidic (Dayhoff-2)", -3.5),
        b'C' => ("Cysteine (Cys)", "Disulfide Thiol (Dayhoff-6)", 2.5),
        b'Q' => ("Glutamine (Gln)", "Polar Uncharged Amide (Dayhoff-2)", -3.5),
        b'E' => ("Glutamate (Glu)", "Negatively Charged Acidic (Dayhoff-2)", -3.5),
        b'G' => ("Glycine (Gly)", "Conformational Special (Dayhoff-1)", -0.4),
        b'H' => ("Histidine (His)", "Aromatic Imidazole Basic (Dayhoff-3)", -3.2),
        b'I' => ("Isoleucine (Ile)", "Hydrophobic Branched (Dayhoff-4)", 4.5),
        b'L' => ("Leucine (Leu)", "Hydrophobic Aliphatic (Dayhoff-4)", 3.8),
        b'K' => ("Lysine (Lys)", "Positively Charged Basic (Dayhoff-3)", -3.9),
        b'M' => ("Methionine (Met)", "Hydrophobic Thioether / Start (Dayhoff-4)", 1.9),
        b'F' => ("Phenylalanine (Phe)", "Hydrophobic Aromatic (Dayhoff-5)", 2.8),
        b'P' => ("Proline (Pro)", "Cyclic Imino Turn Inducer (Dayhoff-1)", -1.6),
        b'S' => ("Serine (Ser)", "Polar Hydroxyl (Dayhoff-1)", -0.8),
        b'T' => ("Threonine (Thr)", "Polar Hydroxyl (Dayhoff-1)", -0.7),
        b'W' => ("Tryptophan (Trp)", "Bulky Indole Aromatic (Dayhoff-5)", -0.9),
        b'Y' => ("Tyrosine (Tyr)", "Phenolic Aromatic (Dayhoff-5)", -1.3),
        b'V' => ("Valine (Val)", "Hydrophobic Branched (Dayhoff-4)", 4.2),
        b'*' => ("Stop Codon (Ter)", "Translation Termination", 0.0),
        b'-' => ("Alignment Gap (-)", "Indel / Missing Column", 0.0),
        _ => ("Ambiguous / Unknown (X)", "Uncertain Residue", 0.0),
    }
}

/// Single aligned taxon sequence entry in the viewer workspace.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ViewerSequenceEntry {
    pub index: usize,
    pub id: String,
    pub taxon: String,
    pub common_name: String,
    pub clade_or_family: String,
    pub aligned_dna: Vec<u8>,
    pub gc_percent: f64,
    pub gap_percent: f64,
    pub rcfv_taxon: f64,
    pub adna_5prime_ct_count: usize,
    pub adna_3prime_ga_count: usize,
}

/// Per-column conservation, entropy, consensus, and parsimony telemetry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ColumnStats {
    pub col_index_1based: usize,
    pub consensus_dna: u8,
    pub consensus_aa: u8,
    pub conservation_fraction: f64,
    pub shannon_entropy_bits: f64,
    pub gap_fraction: f64,
    pub is_variable: bool,
    pub is_parsimony_informative: bool,
    pub codon_omega_estimate: f64,
}

/// Annotated genomic feature, domain, motif, or sliding-window signal on the alignment coordinate axis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentFeatureAnnotation {
    pub feature_id: String,
    pub feature_type: String, // "CDS", "Catalytic_Motif", "Conserved_Block", "Positive_Selection", "Hypervariable_Loop", "aDNA_Deam_Zone", "G4_CpG"
    pub label: String,
    pub start_col_1based: usize,
    pub end_col_1based: usize,
    pub strand: char,
    pub score: f64,
    pub color_rgb: (u8, u8, u8),
    pub description: String,
}

/// Complete Alignment Workspace loaded into the viewer or snapshot generator.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentWorkspace {
    pub title: String,
    pub source_path: String,
    pub alignment_method: String,
    pub sequences: Vec<ViewerSequenceEntry>,
    pub alignment_length: usize,
    pub column_stats: Vec<ColumnStats>,
    pub consensus_dna_string: String,
    pub consensus_aa_string: String,
    pub features: Vec<AlignmentFeatureAnnotation>,
    pub conserved_sites: usize,
    pub variable_sites: usize,
    pub singleton_sites: usize,
    pub parsimony_informative_sites: usize,
    pub mean_gc_percent: f64,
    pub mean_gap_percent: f64,
    pub total_rcfv: f64,
    pub ti_tv_ratio: f64,
    pub detected_protein_input: bool,
}
impl AlignmentWorkspace {
    /// Load an alignment workspace from a FASTA/CLUSTAL/JSON file on disk, optionally accompanied by a GFF3/JSON annotation file.
    pub fn from_file(
        path: &Path,
        annotation_path: Option<&Path>,
        genetic_code: GeneticCodeTable,
        reading_frame: usize,
        auto_align_if_unequal: bool,
    ) -> Result<Self, String> {
        let raw = fs::read_to_string(path)
            .map_err(|e| format!("Failed to read alignment file {:?}: {}", path, e))?;
        let file_stem = path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("Alignment")
            .to_string();

        let records = parse_any_sequence_format(&raw, &file_stem)?;
        if records.is_empty() {
            return Err(format!("No sequences found in {:?}", path));
        }

        let mut ws = Self::from_fasta_records(
            &file_stem,
            &path.display().to_string(),
            &records,
            genetic_code,
            reading_frame,
            auto_align_if_unequal,
        );

        if let Some(ann_p) = annotation_path {
            if let Ok(ext_anns) = parse_external_annotations(ann_p, ws.alignment_length) {
                let mut combined = ext_anns;
                combined.extend(ws.features);
                ws.features = combined;
            }
        }
        Ok(ws)
    }

    /// Load one of the built-in empirical `bio-suite` cohorts from `crvscience-labs/cache/`.
    pub fn from_builtin_cohort(
        cohort_name: &str,
        cache: &PoliteCacheManager,
        genetic_code: GeneticCodeTable,
    ) -> Result<Self, String> {
        match cohort_name.trim().to_ascii_lowercase().as_str() {
            "chr1" | "chr1-section-a" | "human-chr1" => {
                let (rec, anns) = cache.fetch_human_chr1_section_a()?;
                // Slice multiple overlapping functional windows or load cached FASTAs if single sequence
                let mut records = vec![rec.clone()];
                // Also include any homologous mammal/primate loci in cache if present, or sub-windows for multi-row inspection
                if rec.sequence.len() >= 3600 {
                    for (i, label) in [
                        (0, "SAMD11 Promoter & CpG Island Window"),
                        (600, "NOC2L Antisense Exonic Window"),
                        (1200, "KLHL17 Catalytic BTB/Kelch Window"),
                        (1800, "PLEKHN1 PH-Domain Coding Window"),
                        (2400, "PERM1 Regulatory Enhancer Window"),
                    ] {
                        let end = (i + 1200).min(rec.sequence.len());
                        records.push(FastaRecord {
                            id: format!("NC_000001.11_w{}", i),
                            accession: format!("NC_000001.11:{}-{}", i + 1, end),
                            taxon: format!("Homo sapiens ({})", label),
                            common_name: label.to_string(),
                            family: "Hominidae".to_string(),
                            gene: "1p36.33".to_string(),
                            description: label.to_string(),
                            sequence: rec.sequence[i..end].to_string(),
                            cached_from_disk: true,
                        });
                    }
                    records[0].sequence.truncate(1200);
                }
                let mut ws = Self::from_fasta_records(
                    "Human GRCh38 Chr1p36.33 Annotated Functional Panel",
                    "cache/NC_000001.11_chr1_section_a.fasta",
                    &records,
                    genetic_code,
                    1,
                    true,
                );
                for a in anns {
                    let s = a.start.max(1).min(ws.alignment_length);
                    let e = a.end.max(s).min(ws.alignment_length);
                    ws.features.insert(
                        0,
                        AlignmentFeatureAnnotation {
                            feature_id: a.feature_id,
                            feature_type: a.feature_type.clone(),
                            label: format!("{} ({})", a.symbol, a.feature_type),
                            start_col_1based: s,
                            end_col_1based: e,
                            strand: a.strand,
                            score: 100.0,
                            color_rgb: (56, 189, 248),
                            description: a.description,
                        },
                    );
                }
                Ok(ws)
            }
            "hymenoptera" | "bees-ants" => {
                let mut records = Vec::new();
                if let Ok(entries) = fs::read_dir(&cache.cache_dir) {
                    for ent in entries.flatten() {
                        let fname = ent.file_name().to_string_lossy().to_string();
                        if fname.starts_with("HYMEN_") && fname.ends_with(".fasta") {
                            if let Ok(raw) = fs::read_to_string(ent.path()) {
                                if let Ok(mut parsed) = parse_any_sequence_format(&raw, &fname) {
                                    for r in &mut parsed {
                                        if r.sequence.len() > 1800 {
                                            r.sequence.truncate(1800);
                                        }
                                    }
                                    records.extend(parsed);
                                }
                            }
                        }
                    }
                }
                if records.is_empty() {
                    let coi = cache.fetch_carnivora_coi_panel()?;
                    return Ok(Self::from_fasta_records(
                        "Carnivora COX1 Barcode Reference Alignment (16 Taxa, 4 Families)",
                        "cache/carnivora_coi_panel",
                        &coi,
                        GeneticCodeTable::VertebrateMito2,
                        1,
                        true,
                    ));
                }
                Ok(Self::from_fasta_records(
                    "Hymenoptera Eusociality & Venom Clade Alignment",
                    "cache/HYMEN_*.fasta",
                    &records,
                    GeneticCodeTable::InvertebrateMito5,
                    1,
                    true,
                ))
            }
            _ => {
                let records = cache.fetch_carnivora_coi_panel()?;
                Ok(Self::from_fasta_records(
                    "Carnivora COX1 Mitochondrial Barcode Alignment (16 Taxa, 1,545 bp)",
                    "cache/carnivora_coi_panel",
                    &records,
                    GeneticCodeTable::VertebrateMito2,
                    1,
                    true,
                ))
            }
        }
    }

    /// Construct an `AlignmentWorkspace` from a slice of `FastaRecord`s.
    pub fn from_fasta_records(
        title: &str,
        source_path: &str,
        records: &[FastaRecord],
        genetic_code: GeneticCodeTable,
        reading_frame: usize,
        auto_align_if_unequal: bool,
    ) -> Self {
        let min_len = records.iter().map(|r| r.sequence.len()).min().unwrap_or(0);
        let max_len = records.iter().map(|r| r.sequence.len()).max().unwrap_or(0);

        // Detect whether input is already aligned or needs Gotoh/progressive alignment
        let (alignment_method, aligned_pairs): (String, Vec<(FastaRecord, String)>) =
            if min_len == max_len || !auto_align_if_unequal || records.len() <= 1 {
                // Pad to max_len if slight trailing length differences exist
                let pairs = records
                    .iter()
                    .map(|r| {
                        let mut s = r.sequence.trim().to_ascii_uppercase();
                        if s.len() < max_len {
                            s.push_str(&"-".repeat(max_len - s.len()));
                        }
                        (r.clone(), s)
                    })
                    .collect();
                ("Pre-Aligned FASTA Matrix".to_string(), pairs)
            } else if max_len <= 6000 && records.len() <= 64 {
                let (method, aligned, _) = align_sequences_with_mode(records, AlignerMode::Native);
                let pairs = records
                    .iter()
                    .zip(aligned.into_iter())
                    .map(|(rec, aln)| (rec.clone(), aln.aligned_seq.to_ascii_uppercase()))
                    .collect();
                (method, pairs)
            } else {
                // Fast codon-preserving/center-anchored alignment for very long genome tracks
                let pairs = records
                    .iter()
                    .map(|r| {
                        let mut s = r.sequence.trim().to_ascii_uppercase();
                        if s.len() < max_len {
                            s.push_str(&"-".repeat(max_len - s.len()));
                        }
                        (r.clone(), s)
                    })
                    .collect();
                ("Direct Coordinate-Locked Alignment".to_string(), pairs)
            };

        let alignment_length = aligned_pairs.first().map(|(_, s)| s.len()).unwrap_or(0);

        // Detect if sequences are protein (contain amino acid letters beyond A,C,G,T,U,N,-)
        let mut non_dna_chars = 0usize;
        let mut total_alpha_chars = 0usize;
        for (_, seq) in &aligned_pairs {
            for &b in seq.as_bytes() {
                if b.is_ascii_alphabetic() {
                    total_alpha_chars += 1;
                    if !matches!(
                        b,
                        b'A' | b'C' | b'G' | b'T' | b'U' | b'N' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M'
                    ) {
                        non_dna_chars += 1;
                    }
                }
            }
        }
        let detected_protein_input =
            total_alpha_chars > 0 && (non_dna_chars as f64 / total_alpha_chars as f64) > 0.12;

        // Compute global base frequencies for RCFV
        let mut global_counts = [0usize; 4];
        let mut global_total = 0usize;
        for (_, seq) in &aligned_pairs {
            for &b in seq.as_bytes() {
                match b {
                    b'A' => {
                        global_counts[0] += 1;
                        global_total += 1;
                    }
                    b'C' => {
                        global_counts[1] += 1;
                        global_total += 1;
                    }
                    b'G' => {
                        global_counts[2] += 1;
                        global_total += 1;
                    }
                    b'T' | b'U' => {
                        global_counts[3] += 1;
                        global_total += 1;
                    }
                    _ => {}
                }
            }
        }
        let mean_freqs: [f64; 4] = if global_total > 0 {
            [
                global_counts[0] as f64 / global_total as f64,
                global_counts[1] as f64 / global_total as f64,
                global_counts[2] as f64 / global_total as f64,
                global_counts[3] as f64 / global_total as f64,
            ]
        } else {
            [0.25; 4]
        };

        let mut sequences = Vec::with_capacity(aligned_pairs.len());
        let mut total_rcfv = 0.0;

        for (idx, (rec, aln_str)) in aligned_pairs.into_iter().enumerate() {
            let bytes = aln_str.into_bytes();
            let mut counts = [0usize; 4];
            let mut gaps = 0usize;
            for &b in &bytes {
                match b {
                    b'A' => counts[0] += 1,
                    b'C' => counts[1] += 1,
                    b'G' => counts[2] += 1,
                    b'T' | b'U' => counts[3] += 1,
                    b'-' | b'.' => gaps += 1,
                    _ => {}
                }
            }
            let ungapped = (counts[0] + counts[1] + counts[2] + counts[3]).max(1);
            let gc_percent = ((counts[1] + counts[2]) as f64 * 100.0) / ungapped as f64;
            let gap_percent = (gaps as f64 * 100.0) / bytes.len().max(1) as f64;

            let mut rcfv_taxon = 0.0;
            for k in 0..4 {
                let f = counts[k] as f64 / ungapped as f64;
                rcfv_taxon += (f - mean_freqs[k]).abs();
            }
            rcfv_taxon /= 4.0;
            total_rcfv += rcfv_taxon;

            // Count terminal 25-bp ancient-DNA deamination candidates (5' T at C-consensus, 3' A at G-consensus)
            let term_win = bytes.len().min(25);
            let adna_5prime_ct_count = bytes[..term_win].iter().filter(|&&b| b == b'T').count();
            let adna_3prime_ga_count = bytes[bytes.len().saturating_sub(term_win)..]
                .iter()
                .filter(|&&b| b == b'A')
                .count();

            sequences.push(ViewerSequenceEntry {
                index: idx,
                id: rec.id,
                taxon: if rec.taxon.is_empty() {
                    rec.accession.clone()
                } else {
                    rec.taxon
                },
                common_name: rec.common_name,
                clade_or_family: if rec.family.is_empty() {
                    "Cohort".to_string()
                } else {
                    rec.family
                },
                aligned_dna: bytes,
                gc_percent,
                gap_percent,
                rcfv_taxon,
                adna_5prime_ct_count,
                adna_3prime_ga_count,
            });
        }

        if !sequences.is_empty() {
            total_rcfv /= sequences.len() as f64;
        }

        let mut ws = Self {
            title: title.to_string(),
            source_path: source_path.to_string(),
            alignment_method,
            sequences,
            alignment_length,
            column_stats: Vec::new(),
            consensus_dna_string: String::new(),
            consensus_aa_string: String::new(),
            features: Vec::new(),
            conserved_sites: 0,
            variable_sites: 0,
            singleton_sites: 0,
            parsimony_informative_sites: 0,
            mean_gc_percent: 0.0,
            mean_gap_percent: 0.0,
            total_rcfv,
            ti_tv_ratio: 2.0,
            detected_protein_input,
        };
        ws.recompute_telemetry_and_auto_features(genetic_code, reading_frame);
        ws
    }

    /// Recompute per-column conservation, Shannon entropy, consensus sequence, parsimony-informative sites,
    /// and automatically scan for biological features/motifs across the alignment.
    pub fn recompute_telemetry_and_auto_features(
        &mut self,
        genetic_code: GeneticCodeTable,
        reading_frame: usize,
    ) {
        let n_seq = self.sequences.len().max(1);
        let len = self.alignment_length;
        let frame_offset = reading_frame.saturating_sub(1) % 3;

        let mut col_stats = Vec::with_capacity(len);
        let mut consensus_dna = Vec::with_capacity(len);
        let mut conserved = 0usize;
        let mut variable = 0usize;
        let mut singleton = 0usize;
        let mut parsimony_inf = 0usize;

        for col in 0..len {
            let mut freq: BTreeMap<u8, usize> = BTreeMap::new();
            let mut gaps = 0usize;
            for seq in &self.sequences {
                let b = seq.aligned_dna.get(col).copied().unwrap_or(b'-');
                if b == b'-' || b == b'.' {
                    gaps += 1;
                } else {
                    *freq.entry(b).or_insert(0) += 1;
                }
            }
            let non_gap = n_seq.saturating_sub(gaps);
            let (top_base, top_count) = freq
                .iter()
                .max_by_key(|(_, &c)| c)
                .map(|(&b, &c)| (b, c))
                .unwrap_or((b'-', 0));

            consensus_dna.push(top_base);
            let conservation_fraction = if non_gap > 0 {
                top_count as f64 / non_gap as f64
            } else {
                0.0
            };

            let mut entropy = 0.0;
            if non_gap > 0 {
                for &c in freq.values() {
                    if c > 0 {
                        let p = c as f64 / non_gap as f64;
                        entropy -= p * p.log2();
                    }
                }
            }

            let distinct_states = freq.len();
            let states_ge_2 = freq.values().filter(|&&c| c >= 2).count();
            let is_var = distinct_states > 1;
            let is_pi = states_ge_2 >= 2;

            if distinct_states <= 1 && non_gap > 0 {
                conserved += 1;
            } else if is_var {
                variable += 1;
                if is_pi {
                    parsimony_inf += 1;
                } else {
                    singleton += 1;
                }
            }

            col_stats.push(ColumnStats {
                col_index_1based: col + 1,
                consensus_dna: top_base,
                consensus_aa: b' ',
                conservation_fraction,
                shannon_entropy_bits: entropy,
                gap_fraction: gaps as f64 / n_seq as f64,
                is_variable: is_var,
                is_parsimony_informative: is_pi,
                codon_omega_estimate: 0.25,
            });
        }

        // Compute codon-level consensus AA and approximate dN/dS (omega) across triplets
        let mut consensus_aa_bytes = Vec::with_capacity(len / 3 + 1);
        let mut pos = frame_offset;
        while pos + 2 < len {
            let c1 = consensus_dna[pos];
            let c2 = consensus_dna[pos + 1];
            let c3 = consensus_dna[pos + 2];
            let aa = if self.detected_protein_input {
                c1
            } else {
                genetic_code.translate_codon(c1, c2, c3)
            };
            consensus_aa_bytes.push(aa);

            // Estimate nonsynonymous vs synonymous polymorphism at this codon across taxa
            let mut syn_diffs = 0usize;
            let mut nonsyn_diffs = 0usize;
            if !self.detected_protein_input {
                for seq in &self.sequences {
                    let s1 = seq.aligned_dna[pos];
                    let s2 = seq.aligned_dna[pos + 1];
                    let s3 = seq.aligned_dna[pos + 2];
                    if s1 == b'-' || s2 == b'-' || s3 == b'-' {
                        continue;
                    }
                    if (s1, s2, s3) != (c1, c2, c3) {
                        let taxon_aa = genetic_code.translate_codon(s1, s2, s3);
                        if taxon_aa == aa {
                            syn_diffs += 1;
                        } else {
                            nonsyn_diffs += 1;
                        }
                    }
                }
            }
            let omega = if syn_diffs == 0 && nonsyn_diffs == 0 {
                0.15
            } else {
                ((nonsyn_diffs as f64 + 0.25) / ((syn_diffs as f64 + 0.5) * 2.8)).min(3.5)
            };
            for k in 0..3 {
                col_stats[pos + k].consensus_aa = aa;
                col_stats[pos + k].codon_omega_estimate = omega;
            }
            pos += 3;
        }

        // Compute pairwise Transition / Transversion (Ti/Tv) ratio vs consensus
        let mut transitions = 0usize;
        let mut transversions = 0usize;
        for col in 0..len {
            let cb = consensus_dna[col];
            if !matches!(cb, b'A' | b'C' | b'G' | b'T') {
                continue;
            }
            for seq in &self.sequences {
                let sb = seq.aligned_dna[col];
                if sb == cb || !matches!(sb, b'A' | b'C' | b'G' | b'T') {
                    continue;
                }
                let pair = (cb.min(sb), cb.max(sb));
                if pair == (b'A', b'G') || pair == (b'C', b'T') {
                    transitions += 1;
                } else {
                    transversions += 1;
                }
            }
        }
        self.ti_tv_ratio = if transversions > 0 {
            transitions as f64 / transversions as f64
        } else {
            2.0
        };

        self.mean_gc_percent = if self.sequences.is_empty() {
            0.0
        } else {
            self.sequences.iter().map(|s| s.gc_percent).sum::<f64>() / self.sequences.len() as f64
        };
        self.mean_gap_percent = if self.sequences.is_empty() {
            0.0
        } else {
            self.sequences.iter().map(|s| s.gap_percent).sum::<f64>() / self.sequences.len() as f64
        };

        self.conserved_sites = conserved;
        self.variable_sites = variable;
        self.singleton_sites = singleton;
        self.parsimony_informative_sites = parsimony_inf;
        self.consensus_dna_string = String::from_utf8_lossy(&consensus_dna).to_string();
        self.consensus_aa_string = String::from_utf8_lossy(&consensus_aa_bytes).to_string();
        self.column_stats = col_stats;

        if self.features.is_empty() {
            self.features = self.auto_detect_genomic_features(genetic_code, reading_frame);
        }
    }

    /// Automatically scan the alignment for conserved catalytic domains, ORFs, positive selection
    /// windows (`omega > 1.0`), hypervariable indel loops, G-quadruplexes, and terminal aDNA zones.
    pub fn auto_detect_genomic_features(
        &self,
        genetic_code: GeneticCodeTable,
        reading_frame: usize,
    ) -> Vec<AlignmentFeatureAnnotation> {
        let mut features = Vec::new();
        let len = self.alignment_length;
        if len == 0 {
            return features;
        }
        let frame_offset = reading_frame.saturating_sub(1) % 3;

        // 1. Primary Full-Alignment Coding / Locus Feature
        features.push(AlignmentFeatureAnnotation {
            feature_id: "FEAT-FULL-01".to_string(),
            feature_type: "Locus_Overview".to_string(),
            label: format!("Full Alignment Span (1..{} bp)", len),
            start_col_1based: 1,
            end_col_1based: len,
            strand: '+',
            score: 100.0,
            color_rgb: (56, 189, 248),
            description: format!(
                "Complete {} alignment ({} taxa, {} PI sites, {:.1}% GC)",
                self.title,
                self.sequences.len(),
                self.parsimony_informative_sites,
                self.mean_gc_percent
            ),
        });

        // 2. Scan 45-bp sliding windows for Ultraconserved Catalytic Anchors (>= 90% mean conservation)
        let win = 45usize.min(len);
        if win >= 15 {
            let mut best_cons_windows: Vec<(usize, usize, f64)> = Vec::new();
            let mut step = 0usize;
            while step + win <= len {
                let slice = &self.column_stats[step..step + win];
                let mean_cons =
                    slice.iter().map(|c| c.conservation_fraction).sum::<f64>() / win as f64;
                let mean_gaps = slice.iter().map(|c| c.gap_fraction).sum::<f64>() / win as f64;
                if mean_cons >= 0.88 && mean_gaps <= 0.10 {
                    best_cons_windows.push((step + 1, step + win, mean_cons));
                }
                step += 30;
            }
            best_cons_windows.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
            for (idx, (s, e, score)) in best_cons_windows.into_iter().take(4).enumerate() {
                let aa_motif = self.extract_consensus_aa_window(s, e, genetic_code, frame_offset);
                features.push(AlignmentFeatureAnnotation {
                    feature_id: format!("ANCHOR-{:02}", idx + 1),
                    feature_type: "Conserved_Domain".to_string(),
                    label: format!("Conserved Core #{} ({}..{} bp, {:.0}% id)", idx + 1, s, e, score * 100.0),
                    start_col_1based: s,
                    end_col_1based: e,
                    strand: '+',
                    score: score * 100.0,
                    color_rgb: (16, 185, 129),
                    description: format!(
                        "Ultraconserved anchor window ({:.1}% identity; translated consensus motif: {})",
                        score * 100.0,
                        aa_motif
                    ),
                });
            }
        }

        // 3. Scan for High-Diagnostic / Parsimony-Informative Barcode & dN/dS Windows (45 bp)
        if win >= 15 {
            let mut diag_windows: Vec<(usize, usize, usize, f64)> = Vec::new();
            let mut step = 0usize;
            while step + win <= len {
                let slice = &self.column_stats[step..step + win];
                let pi_count = slice.iter().filter(|c| c.is_parsimony_informative).count();
                let mean_omega =
                    slice.iter().map(|c| c.codon_omega_estimate).sum::<f64>() / win as f64;
                if pi_count >= 4 {
                    diag_windows.push((step + 1, step + win, pi_count, mean_omega));
                }
                step += 30;
            }
            diag_windows.sort_by(|a, b| b.2.cmp(&a.2));
            for (idx, (s, e, pi_cnt, omega)) in diag_windows.into_iter().take(4).enumerate() {
                let aa_motif = self.extract_consensus_aa_window(s, e, genetic_code, frame_offset);
                features.push(AlignmentFeatureAnnotation {
                    feature_id: format!("DIAG-PI-{:02}", idx + 1),
                    feature_type: "Diagnostic_Motif".to_string(),
                    label: format!("Diagnostic Clade Window #{} ({}..{} bp, {} PI)", idx + 1, s, e, pi_cnt),
                    start_col_1based: s,
                    end_col_1based: e,
                    strand: '+',
                    score: pi_cnt as f64,
                    color_rgb: (168, 85, 247),
                    description: format!(
                        "High-resolution phylogenetic window ({} parsimony-informative sites, mean omega={:.2}; AA: {})",
                        pi_cnt, omega, aa_motif
                    ),
                });
            }
        }

        // 4. Scan translated consensus for canonical functional protein motifs (PROSITE-style patterns)
        let aa_bytes = self.consensus_aa_string.as_bytes();
        for (codon_idx, window) in aa_bytes.windows(6).enumerate() {
            let start_bp = frame_offset + codon_idx * 3 + 1;
            let end_bp = (start_bp + 17).min(len);
            // Check H-X-X-H / H-E-X-X-H (Heme / Metalloprotease catalytic site)
            if window[0] == b'H' && (window[3] == b'H' || window[4] == b'H') {
                features.push(AlignmentFeatureAnnotation {
                    feature_id: format!("MOTIF-HEME-{}", start_bp),
                    feature_type: "Catalytic_Motif".to_string(),
                    label: format!(
                        "His-Ligand Catalytic Motif ({}) [{}..{}]",
                        String::from_utf8_lossy(window),
                        start_bp,
                        end_bp
                    ),
                    start_col_1based: start_bp,
                    end_col_1based: end_bp,
                    strand: '+',
                    score: 95.0,
                    color_rgb: (245, 158, 11),
                    description: "Conserved Histidine coordination / catalytic motif".to_string(),
                });
                if features.len() > 14 {
                    break;
                }
            } else if window[0] == b'G' && window[2] == b'G' && window[5] == b'G' {
                features.push(AlignmentFeatureAnnotation {
                    feature_id: format!("MOTIF-GXGXXG-{}", start_bp),
                    feature_type: "Catalytic_Motif".to_string(),
                    label: format!(
                        "G-x-G-x-x-G Nucleotide-Binding Loop ({}) [{}..{}]",
                        String::from_utf8_lossy(window),
                        start_bp,
                        end_bp
                    ),
                    start_col_1based: start_bp,
                    end_col_1based: end_bp,
                    strand: '+',
                    score: 92.0,
                    color_rgb: (236, 72, 153),
                    description: "Conserved glycine-rich Rossmann/P-loop motif".to_string(),
                });
                if features.len() > 14 {
                    break;
                }
            }
        }

        // 5. Terminal Ancient-DNA (aDNA) 5' and 3' Deamination Audit Windows
        if len >= 50 {
            features.push(AlignmentFeatureAnnotation {
                feature_id: "ADNA-5PRIME".to_string(),
                feature_type: "aDNA_Deam_Zone".to_string(),
                label: "5' Terminal aDNA C->T Audit Zone (1..25 bp)".to_string(),
                start_col_1based: 1,
                end_col_1based: 25,
                strand: '+',
                score: 80.0,
                color_rgb: (239, 68, 68),
                description: "RuNIE 5' terminal cytosine deamination (C->T) inspection & masking window".to_string(),
            });
            features.push(AlignmentFeatureAnnotation {
                feature_id: "ADNA-3PRIME".to_string(),
                feature_type: "aDNA_Deam_Zone".to_string(),
                label: format!("3' Terminal aDNA G->A Audit Zone ({}..{} bp)", len - 24, len),
                start_col_1based: len - 24,
                end_col_1based: len,
                strand: '-',
                score: 80.0,
                color_rgb: (239, 68, 68),
                description: "RuNIE 3' terminal guanine complement deamination (G->A) inspection window".to_string(),
            });
        }

        features.sort_by_key(|f| f.start_col_1based);
        features
    }

    pub fn extract_consensus_aa_window(
        &self,
        start_1based: usize,
        end_1based: usize,
        genetic_code: GeneticCodeTable,
        frame_offset: usize,
    ) -> String {
        let mut out = String::new();
        let s0 = start_1based.saturating_sub(1);
        let e0 = end_1based.min(self.alignment_length);
        let mut p = s0;
        while (p.saturating_sub(frame_offset)) % 3 != 0 && p < e0 {
            p += 1;
        }
        let bytes = self.consensus_dna_string.as_bytes();
        while p + 2 < e0 {
            let aa = if self.detected_protein_input {
                bytes[p]
            } else {
                genetic_code.translate_codon(bytes[p], bytes[p + 1], bytes[p + 2])
            };
            out.push(aa as char);
            p += 3;
        }
        out
    }

    /// Translate a specific taxon row (`seq_idx`) at 0-based nucleotide column `col` under the given
    /// `genetic_code` and `reading_frame` (1, 2, or 3). Returns `(codon_position_0_1_2, translated_aa_byte, [c1, c2, c3])`.
    pub fn codon_and_aa_at(
        &self,
        seq_idx: usize,
        col: usize,
        genetic_code: GeneticCodeTable,
        reading_frame: usize,
    ) -> (usize, u8, [u8; 3]) {
        let seq = match self.sequences.get(seq_idx) {
            Some(s) => &s.aligned_dna,
            None => return (0, b'-', [b'-'; 3]),
        };
        if self.detected_protein_input {
            let aa = seq.get(col).copied().unwrap_or(b'-');
            return (1, aa, [aa, aa, aa]);
        }
        let frame_offset = reading_frame.saturating_sub(1) % 3;
        if col < frame_offset {
            let b = seq.get(col).copied().unwrap_or(b'-');
            return (0, b' ', [b, b'-', b'-']);
        }
        let rel = col - frame_offset;
        let codon_pos = rel % 3;
        let codon_start = col - codon_pos;
        if codon_start + 2 >= seq.len() {
            let b = seq.get(col).copied().unwrap_or(b'-');
            return (codon_pos, b' ', [b, b'-', b'-']);
        }
        let c1 = seq[codon_start];
        let c2 = seq[codon_start + 1];
        let c3 = seq[codon_start + 2];
        let aa = genetic_code.translate_codon(c1, c2, c3);
        (codon_pos, aa, [c1, c2, c3])
    }

    /// Search for a nucleotide or amino-acid motif across all sequences or the consensus.
    /// Returns matching `(seq_idx, start_col_1based, end_col_1based, matched_string)` tuples.
    pub fn search_motif(
        &self,
        query: &str,
        genetic_code: GeneticCodeTable,
        reading_frame: usize,
    ) -> Vec<(usize, usize, usize, String)> {
        let q = query.trim().to_ascii_uppercase();
        if q.is_empty() {
            return Vec::new();
        }
        let mut hits = Vec::new();
        let frame_offset = reading_frame.saturating_sub(1) % 3;

        for (s_idx, seq) in self.sequences.iter().enumerate() {
            let dna_str = String::from_utf8_lossy(&seq.aligned_dna);
            for (pos, _) in dna_str.match_indices(&q) {
                hits.push((s_idx, pos + 1, pos + q.len(), q.clone()));
                if hits.len() >= 100 {
                    return hits;
                }
            }
            // Also search translated amino-acid sequence in the active reading frame
            if !self.detected_protein_input && seq.aligned_dna.len() >= frame_offset + 3 {
                let mut aa_seq = Vec::new();
                let mut p = frame_offset;
                while p + 2 < seq.aligned_dna.len() {
                    aa_seq.push(genetic_code.translate_codon(
                        seq.aligned_dna[p],
                        seq.aligned_dna[p + 1],
                        seq.aligned_dna[p + 2],
                    ));
                    p += 3;
                }
                let aa_str = String::from_utf8_lossy(&aa_seq);
                for (aa_pos, _) in aa_str.match_indices(&q) {
                    let start_bp = frame_offset + aa_pos * 3 + 1;
                    let end_bp = start_bp + q.len() * 3 - 1;
                    hits.push((s_idx, start_bp, end_bp, q.clone()));
                    if hits.len() >= 100 {
                        return hits;
                    }
                }
            }
        }
        hits
    }

    /// Select the highest-information focal window of width `window_bp` (default 60 bp = 20 codons)
    /// for automatic publication Figure 2 snapshots.
    pub fn select_best_focal_window(&self, window_bp: usize) -> (usize, usize) {
        let len = self.alignment_length;
        if len <= window_bp {
            return (1, len.max(1));
        }
        let w = (window_bp / 3 * 3).max(15).min(len);
        let mut best_start = 1usize;
        let mut best_score = -1.0e9f64;

        let mut s = 0usize;
        while s + w <= len {
            let slice = &self.column_stats[s..s + w];
            let mean_gaps = slice.iter().map(|c| c.gap_fraction).sum::<f64>() / w as f64;
            let pi_count = slice.iter().filter(|c| c.is_parsimony_informative).count();
            let cons_count = slice
                .iter()
                .filter(|c| c.conservation_fraction >= 0.85)
                .count();
            // Favor windows with low gaps, clear conserved anchor flanking, and high parsimony-informative density
            let score = (pi_count as f64) * 3.5 + (cons_count as f64) * 0.8 - mean_gaps * 120.0;
            if score > best_score {
                best_score = score;
                best_start = s + 1;
            }
            s += 3;
        }
        (best_start, (best_start + w - 1).min(len))
    }

    /// RuNIE Pipeline Operator: Mask ancient-DNA terminal 5' C->T and 3' G->A deamination sites (within 25 bp of ends).
    pub fn apply_adna_deamination_mask(&mut self, genetic_code: GeneticCodeTable, reading_frame: usize) -> usize {
        let len = self.alignment_length;
        let term = len.min(25);
        let cons = self.consensus_dna_string.as_bytes().to_vec();
        let mut masked_count = 0usize;

        for seq in &mut self.sequences {
            for col in 0..term {
                if cons.get(col) == Some(&b'C') && seq.aligned_dna.get(col) == Some(&b'T') {
                    seq.aligned_dna[col] = b'N';
                    masked_count += 1;
                }
            }
            for col in len.saturating_sub(term)..len {
                if cons.get(col) == Some(&b'G') && seq.aligned_dna.get(col) == Some(&b'A') {
                    seq.aligned_dna[col] = b'N';
                    masked_count += 1;
                }
            }
        }
        self.recompute_telemetry_and_auto_features(genetic_code, reading_frame);
        masked_count
    }

    /// RuNIE Pipeline Operator: Strip alignment columns exceeding `max_gap_fraction` (e.g. > 0.50).
    pub fn strip_gap_heavy_columns(
        &mut self,
        max_gap_fraction: f64,
        genetic_code: GeneticCodeTable,
        reading_frame: usize,
    ) -> usize {
        let keep_cols: Vec<usize> = self
            .column_stats
            .iter()
            .enumerate()
            .filter(|(_, c)| c.gap_fraction <= max_gap_fraction)
            .map(|(idx, _)| idx)
            .collect();
        let removed = self.alignment_length.saturating_sub(keep_cols.len());
        if removed == 0 || keep_cols.is_empty() {
            return 0;
        }
        for seq in &mut self.sequences {
            let mut filtered = Vec::with_capacity(keep_cols.len());
            for &c in &keep_cols {
                if let Some(&b) = seq.aligned_dna.get(c) {
                    filtered.push(b);
                }
            }
            seq.aligned_dna = filtered;
        }
        self.alignment_length = keep_cols.len();
        self.features.clear();
        self.recompute_telemetry_and_auto_features(genetic_code, reading_frame);
        removed
    }

    /// RuNIE Pipeline Operator: Reverse-complement all sequences in the alignment workspace.
    pub fn reverse_complement_all(&mut self, genetic_code: GeneticCodeTable, reading_frame: usize) {
        if self.detected_protein_input {
            return;
        }
        for seq in &mut self.sequences {
            seq.aligned_dna.reverse();
            for b in &mut seq.aligned_dna {
                *b = match *b {
                    b'A' => b'T',
                    b'T' | b'U' => b'A',
                    b'C' => b'G',
                    b'G' => b'C',
                    b'R' => b'Y',
                    b'Y' => b'R',
                    other => other,
                };
            }
        }
        self.features.clear();
        self.recompute_telemetry_and_auto_features(genetic_code, reading_frame);
    }

    /// Export current alignment (or a 1-based column range `start_1based..=end_1based`) as FASTA string.
    pub fn export_fasta(&self, range_1based: Option<(usize, usize)>, translated_aa: bool, genetic_code: GeneticCodeTable, reading_frame: usize) -> String {
        let (s1, e1) = range_1based.unwrap_or((1, self.alignment_length));
        let s0 = s1.saturating_sub(1).min(self.alignment_length);
        let e0 = e1.min(self.alignment_length);
        let frame_offset = reading_frame.saturating_sub(1) % 3;

        let mut out = String::new();
        for seq in &self.sequences {
            out.push_str(&format!(
                ">{} | {} | {} | cols {}..{}\n",
                seq.id, seq.taxon, seq.clade_or_family, s0 + 1, e0
            ));
            if translated_aa && !self.detected_protein_input {
                let mut p = s0 + frame_offset;
                let mut aa_line = String::new();
                while p + 2 < e0 {
                    let aa = genetic_code.translate_codon(
                        seq.aligned_dna[p],
                        seq.aligned_dna[p + 1],
                        seq.aligned_dna[p + 2],
                    );
                    aa_line.push(aa as char);
                    p += 3;
                }
                out.push_str(&aa_line);
                out.push('\n');
            } else {
                out.push_str(&String::from_utf8_lossy(&seq.aligned_dna[s0..e0]));
                out.push('\n');
            }
        }
        out
    }
}

/// Parse multi-FASTA, CLUSTAL, PHYLIP, or `bio-suite` JSON sequence records.
pub fn parse_any_sequence_format(raw: &str, default_title: &str) -> Result<Vec<FastaRecord>, String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err("Empty sequence file".to_string());
    }

    // 1. Check JSON array of FastaRecord
    if trimmed.starts_with('[') || trimmed.starts_with('{') {
        if let Ok(records) = serde_json::from_str::<Vec<FastaRecord>>(trimmed) {
            if !records.is_empty() {
                return Ok(records);
            }
        }
    }

    // 2. Standard FASTA (`>header`)
    if trimmed.contains('>') {
        let mut records = Vec::new();
        let mut cur_header = String::new();
        let mut cur_seq = String::new();

        for line in trimmed.lines() {
            let l = line.trim();
            if l.is_empty() || l.starts_with(';') {
                continue;
            }
            if let Some(hdr) = l.strip_prefix('>') {
                if !cur_header.is_empty() && !cur_seq.is_empty() {
                    records.push(build_fasta_record_from_header(&cur_header, &cur_seq));
                }
                cur_header = hdr.trim().to_string();
                cur_seq.clear();
            } else {
                for ch in l.chars() {
                    if !ch.is_whitespace() && !ch.is_ascii_digit() {
                        cur_seq.push(ch.to_ascii_uppercase());
                    }
                }
            }
        }
        if !cur_header.is_empty() && !cur_seq.is_empty() {
            records.push(build_fasta_record_from_header(&cur_header, &cur_seq));
        }
        return Ok(records);
    }

    // 3. CLUSTAL / PHYLIP interleaved fallback
    let mut seq_map: Vec<(String, String)> = Vec::new();
    for line in trimmed.lines() {
        let l = line.trim();
        if l.is_empty()
            || l.starts_with("CLUSTAL")
            || l.starts_with("MUSCLE")
            || l.starts_with('#')
            || l.chars().all(|c| c == '*' || c == ':' || c == '.' || c == ' ')
        {
            continue;
        }
        let parts: Vec<&str> = l.split_whitespace().collect();
        if parts.len() >= 2 {
            // Skip PHYLIP numeric header `N L`
            if parts.len() == 2 && parts[0].parse::<usize>().is_ok() && parts[1].parse::<usize>().is_ok() {
                continue;
            }
            let id = parts[0].to_string();
            let seq_part: String = parts[1]
                .chars()
                .filter(|c| !c.is_whitespace() && !c.is_ascii_digit())
                .map(|c| c.to_ascii_uppercase())
                .collect();
            if !seq_part.is_empty() {
                if let Some(existing) = seq_map.iter_mut().find(|(k, _)| k == &id) {
                    existing.1.push_str(&seq_part);
                } else {
                    seq_map.push((id, seq_part));
                }
            }
        }
    }
    if !seq_map.is_empty() {
        return Ok(seq_map
            .into_iter()
            .map(|(id, seq)| build_fasta_record_from_header(&id, &seq))
            .collect());
    }

    // 4. Raw single sequence fallback
    let clean: String = trimmed
        .chars()
        .filter(|c| c.is_ascii_alphabetic() || *c == '-')
        .map(|c| c.to_ascii_uppercase())
        .collect();
    if !clean.is_empty() {
        return Ok(vec![build_fasta_record_from_header(default_title, &clean)]);
    }

    Err("Unrecognized alignment format".to_string())
}

fn build_fasta_record_from_header(header: &str, seq: &str) -> FastaRecord {
    let parts: Vec<&str> = header.split('|').map(|s| s.trim()).collect();
    let first_tok = header.split_whitespace().next().unwrap_or("SEQ").to_string();
    let (id, taxon, family) = if parts.len() >= 3 {
        (
            parts[0].to_string(),
            parts[1].to_string(),
            parts[2].to_string(),
        )
    } else {
        // Extract binomial taxon name if present in header
        let words: Vec<&str> = header.split_whitespace().collect();
        let taxon_guess = if words.len() >= 3 {
            format!("{} {}", words[1], words[2])
        } else {
            header.chars().take(36).collect()
        };
        (first_tok.clone(), taxon_guess, "Alignment".to_string())
    };
    FastaRecord {
        id: id.clone(),
        accession: first_tok,
        taxon,
        common_name: String::new(),
        family,
        gene: String::new(),
        description: header.to_string(),
        sequence: seq.to_string(),
        cached_from_disk: true,
    }
}

/// Parse external GFF3, BED, or JSON (`Vec<GenomicAnnotation>`) feature annotations.
pub fn parse_external_annotations(
    path: &Path,
    max_len: usize,
) -> Result<Vec<AlignmentFeatureAnnotation>, String> {
    let raw = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let trimmed = raw.trim();
    let mut out = Vec::new();

    if trimmed.starts_with('[') {
        if let Ok(anns) = serde_json::from_str::<Vec<GenomicAnnotation>>(trimmed) {
            for a in anns {
                let s = a.start.max(1).min(max_len.max(1));
                let e = a.end.max(s).min(max_len.max(1));
                out.push(AlignmentFeatureAnnotation {
                    feature_id: a.feature_id,
                    feature_type: a.feature_type.clone(),
                    label: format!("{} ({})", a.symbol, a.feature_type),
                    start_col_1based: s,
                    end_col_1based: e,
                    strand: a.strand,
                    score: 100.0,
                    color_rgb: (56, 189, 248),
                    description: a.description,
                });
            }
            return Ok(out);
        }
    }

    // Parse GFF3 / GTF / BED lines
    for (idx, line) in trimmed.lines().enumerate() {
        let l = line.trim();
        if l.is_empty() || l.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = l.split('\t').collect();
        if cols.len() >= 8 {
            let ftype = cols[2].to_string();
            let start = cols[3].parse::<usize>().unwrap_or(1).max(1).min(max_len.max(1));
            let end = cols[4].parse::<usize>().unwrap_or(start).max(start).min(max_len.max(1));
            let strand = cols[6].chars().next().unwrap_or('+');
            let attrs = cols[8];
            let label = attrs
                .split(';')
                .find_map(|kv| {
                    let kv = kv.trim();
                    kv.strip_prefix("Name=")
                        .or_else(|| kv.strip_prefix("gene="))
                        .or_else(|| kv.strip_prefix("ID="))
                        .or_else(|| kv.strip_prefix("product="))
                })
                .unwrap_or(&ftype)
                .to_string();
            out.push(AlignmentFeatureAnnotation {
                feature_id: format!("GFF-{:03}", idx + 1),
                feature_type: ftype.clone(),
                label: format!("{} [{}..{}]", label, start, end),
                start_col_1based: start,
                end_col_1based: end,
                strand,
                score: 100.0,
                color_rgb: (52, 211, 153),
                description: attrs.to_string(),
            });
        } else if cols.len() >= 3 {
            // BED: chrom start end [name]
            if let (Ok(s0), Ok(e0)) = (cols[1].parse::<usize>(), cols[2].parse::<usize>()) {
                let s = (s0 + 1).min(max_len.max(1));
                let e = e0.max(s).min(max_len.max(1));
                let name = cols.get(3).copied().unwrap_or("BED_Feature");
                out.push(AlignmentFeatureAnnotation {
                    feature_id: format!("BED-{:03}", idx + 1),
                    feature_type: "Region".to_string(),
                    label: format!("{} [{}..{}]", name, s, e),
                    start_col_1based: s,
                    end_col_1based: e,
                    strand: '+',
                    score: 100.0,
                    color_rgb: (251, 191, 36),
                    description: name.to_string(),
                });
            }
        }
    }
    Ok(out)
}

/// Compute `(background_rgb, foreground_text_rgb)` for a residue under the active
/// `ViewAlphabetMode`, `AlignmentColorScheme`, and column conservation telemetry.
pub fn residue_cell_colors(
    residue: u8,
    is_amino_acid: bool,
    mode: ViewAlphabetMode,
    scheme: AlignmentColorScheme,
    col_stat: Option<&ColumnStats>,
) -> ((u8, u8, u8), (u8, u8, u8)) {
    let ch = residue.to_ascii_uppercase();
    if ch == b'-' || ch == b'.' || ch == b' ' {
        return ((15, 23, 42), (100, 116, 139));
    }

    // 1. RY3 Mode
    if mode == ViewAlphabetMode::Ry3Recoded && !is_amino_acid {
        return match ch {
            b'R' => ((14, 116, 144), (240, 249, 255)), // Purine (Cyan-Blue)
            b'Y' => ((180, 83, 9), (255, 251, 235)),   // Pyrimidine (Amber-Orange)
            _ => ((30, 41, 59), (148, 163, 184)),
        };
    }

    // 2. Dayhoff-6 Mode
    if mode == ViewAlphabetMode::Dayhoff6Recoded {
        return match ch {
            b'1' => ((16, 185, 129), (6, 46, 34)),   // AGPST
            b'2' => ((168, 85, 247), (250, 245, 255)), // DENQ
            b'3' => ((239, 68, 68), (254, 242, 242)),  // HKR
            b'4' => ((59, 130, 246), (239, 246, 255)), // ILMV
            b'5' => ((6, 182, 212), (8, 51, 68)),      // FYW
            b'6' => ((234, 179, 8), (28, 25, 23)),     // C
            _ => ((30, 41, 59), (148, 163, 184)),
        };
    }

    // 3. RuNIE Conservation & Parsimony-Informative Highlighter
    if scheme == AlignmentColorScheme::RunieConservation {
        if let Some(st) = col_stat {
            if st.is_parsimony_informative {
                if ch != st.consensus_dna && ch != st.consensus_aa {
                    return ((220, 38, 38), (255, 255, 255)); // Variant allele at PI site
                } else {
                    return ((37, 99, 235), (255, 255, 255)); // Consensus allele at PI site
                }
            } else if st.conservation_fraction >= 0.98 {
                return ((30, 41, 59), (148, 163, 184)); // Muted invariant background
            }
        }
    }

    // 4. Nucleotide Coloring (AliView / Okabe-Ito)
    if !is_amino_acid {
        return match scheme {
            AlignmentColorScheme::OkabeItoSafe => match ch {
                b'A' => ((0, 158, 115), (255, 255, 255)),   // Bluish Green
                b'C' => ((0, 114, 178), (255, 255, 255)),   // Blue
                b'G' => ((230, 159, 0), (15, 23, 42)),      // Orange
                b'T' | b'U' => ((213, 94, 0), (255, 255, 255)), // Vermilion
                _ => ((51, 65, 85), (226, 232, 240)),
            },
            _ => match ch {
                // AliView signature nucleotide colors
                b'A' => ((34, 197, 94), (5, 46, 22)),       // A = Bright Green
                b'C' => ((59, 130, 246), (239, 246, 255)),  // C = Vivid Blue
                b'G' => ((245, 158, 11), (28, 25, 23)),     // G = Amber/Gold
                b'T' | b'U' => ((239, 68, 68), (254, 242, 242)), // T = Crimson Red
                _ => ((51, 65, 85), (226, 232, 240)),
            },
        };
    }

    // 5. Amino Acid Coloring (AliView / ClustalX / Hydrophobicity)
    if scheme == AlignmentColorScheme::Hydrophobicity {
        let (_, _, kd) = amino_acid_full_name_and_chemistry(ch);
        // kd ranges from -4.5 (Arg, hydrophilic blue) to +4.5 (Ile, hydrophobic red)
        let norm = ((kd + 4.5) / 9.0).clamp(0.0, 1.0);
        let r = (30.0 + norm * 210.0) as u8;
        let g = (95.0 - (norm - 0.5).abs() * 80.0) as u8;
        let b = (235.0 - norm * 200.0) as u8;
        return ((r, g, b), (255, 255, 255));
    }

    // AliView / ClustalX Amino Acid chemical groups
    match ch {
        // Hydrophobic aliphatic / aromatic: Blue
        b'A' | b'I' | b'L' | b'M' | b'F' | b'W' | b'V' => ((37, 99, 235), (255, 255, 255)),
        // Positively charged basic (Lys, Arg): Red
        b'K' | b'R' => ((220, 38, 38), (255, 255, 255)),
        // Negatively charged acidic (Asp, Glu): Magenta / Purple
        b'D' | b'E' => ((147, 51, 234), (255, 255, 255)),
        // Polar uncharged hydroxyl / amide (Ser, Thr, Asn, Gln): Emerald Green
        b'N' | b'Q' | b'S' | b'T' => ((16, 185, 129), (6, 46, 34)),
        // Cysteine (disulfide bridge): Yellow
        b'C' => ((234, 179, 8), (28, 25, 23)),
        // Glycine: Orange
        b'G' => ((249, 115, 22), (255, 255, 255)),
        // Proline: Gold-Olive
        b'P' => ((202, 138, 4), (255, 255, 255)),
        // Aromatic polar (His, Tyr): Cyan-Teal
        b'H' | b'Y' => ((8, 145, 178), (255, 255, 255)),
        // Stop codon (*)
        b'*' => ((15, 23, 42), (248, 113, 113)),
        _ => ((51, 65, 85), (226, 232, 240)),
    }
}
/// Configuration for taking an alignment snapshot (SVG, PNG, and/or JSON report) for inclusion in papers.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnapshotConfig {
    pub start_col_1based: usize,
    pub end_col_1based: usize,
    pub mode: ViewAlphabetMode,
    pub color_scheme: AlignmentColorScheme,
    pub match_mode: ConsensusMatchMode,
    pub genetic_code: GeneticCodeTable,
    pub reading_frame: usize,
    pub max_taxa_rows: usize,
    pub subtitle_note: String,
}

impl Default for SnapshotConfig {
    fn default() -> Self {
        Self {
            start_col_1based: 1,
            end_col_1based: 60,
            mode: ViewAlphabetMode::DualDnaAndAminoAcid,
            color_scheme: AlignmentColorScheme::AliViewDefault,
            match_mode: ConsensusMatchMode::ShowAllCharacters,
            genetic_code: GeneticCodeTable::VertebrateMito2,
            reading_frame: 1,
            max_taxa_rows: 24,
            subtitle_note: String::new(),
        }
    }
}

fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

/// Generate a publication-grade, pixel-locked SVG multiple sequence alignment snapshot
/// (modeled on CRV Science Labs Figure 2 standards, with zero duplicate top Figure title banner
/// and explicit `<text x="...">` coordinates on every residue cell).
pub fn render_alignment_snapshot_svg(ws: &AlignmentWorkspace, cfg: &SnapshotConfig) -> String {
    let s1 = cfg.start_col_1based.max(1).min(ws.alignment_length.max(1));
    let e1 = cfg
        .end_col_1based
        .max(s1)
        .min(ws.alignment_length.max(1));
    let s0 = s1 - 1;
    let e0 = e1;
    let num_cols = (e0 - s0).max(1);
    let n_rows = ws.sequences.len().min(cfg.max_taxa_rows.max(1));

    let col_w: f64 = if cfg.mode == ViewAlphabetMode::AminoAcidCompact {
        22.0
    } else if num_cols > 90 {
        11.0
    } else if num_cols > 65 {
        14.0
    } else {
        17.0
    };
    let row_h: f64 = match cfg.mode {
        ViewAlphabetMode::DualDnaAndAminoAcid => 30.0,
        _ => 21.0,
    };

    let left_label_w = 270.0;
    let right_pad = 95.0;
    let display_cols = if cfg.mode == ViewAlphabetMode::AminoAcidCompact {
        (num_cols / 3).max(1)
    } else {
        num_cols
    };
    let matrix_w = display_cols as f64 * col_w;
    let svg_w = (left_label_w + matrix_w + right_pad).max(880.0);

    let top_bar_h = 54.0;
    let feat_track_h = 34.0;
    let cons_track_h = 44.0;
    let ruler_h = 24.0;
    let matrix_y0 = top_bar_h + feat_track_h + cons_track_h + ruler_h;
    let matrix_h = n_rows as f64 * row_h;
    let consensus_row_h = 42.0;
    let legend_h = 44.0;
    let svg_h = matrix_y0 + matrix_h + consensus_row_h + legend_h + 16.0;

    let mut svg = String::with_capacity(64_000);
    svg.push_str(&format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {:.0} {:.0}" width="100%" height="100%" font-family="'JetBrains Mono', 'Fira Code', Consolas, monospace">"#,
        svg_w, svg_h
    ));
    svg.push_str("\n  <defs>\n    <linearGradient id=\"hdrGrad\" x1=\"0\" y1=\"0\" x2=\"1\" y2=\"0\">\n      <stop offset=\"0%\" stop-color=\"#0f172a\"/>\n      <stop offset=\"100%\" stop-color=\"#1e293b\"/>\n    </linearGradient>\n  </defs>\n");

    // Background card
    svg.push_str(&format!(
        "  <rect x=\"0\" y=\"0\" width=\"{:.0}\" height=\"{:.0}\" rx=\"10\" fill=\"#f8fafc\" stroke=\"#cbd5e1\" stroke-width=\"1.5\"/>\n",
        svg_w, svg_h
    ));

    // Top metadata header bar (no "Figure N" duplicate banner per Guardrail 2.00.5)
    svg.push_str(&format!(
        "  <rect x=\"10\" y=\"8\" width=\"{:.0}\" height=\"40\" rx=\"7\" fill=\"url(#hdrGrad)\"/>\n",
        svg_w - 20.0
    ));
    let win_pi = ws.column_stats[s0..e0]
        .iter()
        .filter(|c| c.is_parsimony_informative)
        .count();
    let win_cons = ws.column_stats[s0..e0]
        .iter()
        .map(|c| c.conservation_fraction)
        .sum::<f64>()
        / num_cols as f64
        * 100.0;
    svg.push_str(&format!(
        "  <text x=\"22\" y=\"25\" fill=\"#38bdf8\" font-size=\"11.5\" font-weight=\"700\">Panel A | {} — Focal Alignment Window ({}..{} bp of {} aligned bp)</text>\n",
        xml_escape(&ws.title),
        s1,
        e1,
        ws.alignment_length
    ));
    svg.push_str(&format!(
        "  <text x=\"22\" y=\"41\" fill=\"#cbd5e1\" font-size=\"9.8\">View: {} | Palette: {} | Genetic Code: {} (Frame +{}) | Window Conservation: {:.1}% | Window PI Sites: {} | Global RCFV: {:.4}</text>\n",
        xml_escape(cfg.mode.label()),
        xml_escape(cfg.color_scheme.label()),
        xml_escape(cfg.genetic_code.label()),
        cfg.reading_frame,
        win_cons,
        win_pi,
        ws.total_rcfv
    ));

    // Feature & Annotation Track
    let feat_y = top_bar_h;
    svg.push_str(&format!(
        "  <text x=\"18\" y=\"{:.1}\" fill=\"#334155\" font-size=\"10\" font-weight=\"700\">Annotated Features</text>\n",
        feat_y + 18.0
    ));
    svg.push_str(&format!(
        "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"24\" rx=\"4\" fill=\"#e2e8f0\" stroke=\"#cbd5e1\" stroke-width=\"0.8\"/>\n",
        left_label_w,
        feat_y + 3.0,
        matrix_w
    ));

    let overlapping_feats: Vec<&AlignmentFeatureAnnotation> = ws
        .features
        .iter()
        .filter(|f| {
            f.feature_type != "Locus_Overview"
                && f.start_col_1based <= e1
                && f.end_col_1based >= s1
        })
        .collect();

    if overlapping_feats.is_empty() {
        svg.push_str(&format!(
            "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"18\" rx=\"3\" fill=\"#0284c7\" fill-opacity=\"0.18\" stroke=\"#0284c7\" stroke-width=\"1\"/>\n",
            left_label_w + 2.0,
            feat_y + 6.0,
            (matrix_w - 4.0).max(10.0)
        ));
        svg.push_str(&format!(
            "  <text x=\"{:.1}\" y=\"{:.1}\" fill=\"#0369a1\" font-size=\"9.5\" font-weight=\"700\">Coding / Homology Span ({}..{} bp, Frame +{})</text>\n",
            left_label_w + 10.0,
            feat_y + 18.5,
            s1,
            e1,
            cfg.reading_frame
        ));
    } else {
        for f in overlapping_feats.iter().take(4) {
            let fs = f.start_col_1based.max(s1);
            let fe = f.end_col_1based.min(e1);
            let fx = left_label_w
                + ((fs - s1) as f64 / num_cols as f64) * matrix_w;
            let fw = (((fe - fs + 1) as f64 / num_cols as f64) * matrix_w).max(24.0);
            let (r, g, b) = f.color_rgb;
            svg.push_str(&format!(
                "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"18\" rx=\"3\" fill=\"rgb({},{},{})\" fill-opacity=\"0.22\" stroke=\"rgb({},{},{})\" stroke-width=\"1.2\"/>\n",
                fx,
                feat_y + 6.0,
                fw,
                r, g, b,
                r, g, b
            ));
            svg.push_str(&format!(
                "  <text x=\"{:.1}\" y=\"{:.1}\" fill=\"#0f172a\" font-size=\"9\" font-weight=\"700\">{}</text>\n",
                fx + 5.0,
                feat_y + 18.0,
                xml_escape(&f.label)
            ));
        }
    }

    // Conservation & Parsimony-Informative Histogram Track
    let cons_y = top_bar_h + feat_track_h;
    svg.push_str(&format!(
        "  <text x=\"18\" y=\"{:.1}\" fill=\"#334155\" font-size=\"10\" font-weight=\"700\">Conservation &amp; PI</text>\n",
        cons_y + 18.0
    ));
    svg.push_str(&format!(
        "  <text x=\"18\" y=\"{:.1}\" fill=\"#64748b\" font-size=\"8.5\">Bars = Id% | ◆ = PI site</text>\n",
        cons_y + 31.0
    ));
    svg.push_str(&format!(
        "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"36\" fill=\"#f1f5f9\" stroke=\"#cbd5e1\" stroke-width=\"0.8\"/>\n",
        left_label_w, cons_y + 2.0, matrix_w
    ));

    for d_idx in 0..display_cols {
        let col_idx = if cfg.mode == ViewAlphabetMode::AminoAcidCompact {
            s0 + d_idx * 3
        } else {
            s0 + d_idx
        };
        if let Some(st) = ws.column_stats.get(col_idx) {
            let x = left_label_w + d_idx as f64 * col_w;
            let bar_h = (st.conservation_fraction * 28.0).clamp(2.0, 28.0);
            let bar_y = cons_y + 34.0 - bar_h;
            let bar_color = if st.is_parsimony_informative {
                "#7c3aed"
            } else if st.conservation_fraction >= 0.95 {
                "#0ea5e9"
            } else {
                "#94a3b8"
            };
            svg.push_str(&format!(
                "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"{}\"/>\n",
                x + 1.2,
                bar_y,
                (col_w - 2.4).max(2.0),
                bar_h,
                bar_color
            ));
            if st.is_parsimony_informative {
                svg.push_str(&format!(
                    "  <circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"2.2\" fill=\"#dc2626\"/>\n",
                    x + col_w * 0.5,
                    cons_y + 6.5
                ));
            }
        }
    }

    // Coordinate Ruler
    let ruler_y = top_bar_h + feat_track_h + cons_track_h;
    svg.push_str(&format!(
        "  <text x=\"18\" y=\"{:.1}\" fill=\"#475569\" font-size=\"9.5\" font-weight=\"700\">Taxon / Specimen (GC%)</text>\n",
        ruler_y + 15.0
    ));
    for d_idx in 0..display_cols {
        let col_1based = if cfg.mode == ViewAlphabetMode::AminoAcidCompact {
            s1 + d_idx * 3
        } else {
            s1 + d_idx
        };
        let x_center = left_label_w + d_idx as f64 * col_w + col_w * 0.5;
        let is_major = if cfg.mode == ViewAlphabetMode::AminoAcidCompact {
            d_idx % 5 == 0
        } else {
            col_1based % 10 == 0 || d_idx == 0
        };
        if is_major {
            svg.push_str(&format!(
                "  <line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#475569\" stroke-width=\"1.2\"/>\n",
                x_center,
                ruler_y + 14.0,
                x_center,
                ruler_y + 22.0
            ));
            svg.push_str(&format!(
                "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" fill=\"#1e293b\" font-size=\"8.8\" font-weight=\"700\">{}</text>\n",
                x_center,
                ruler_y + 11.5,
                col_1based
            ));
        } else if col_1based % 5 == 0 {
            svg.push_str(&format!(
                "  <line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#94a3b8\" stroke-width=\"0.9\"/>\n",
                x_center,
                ruler_y + 17.0,
                x_center,
                ruler_y + 22.0
            ));
        }
    }

    // Reference sequence (Taxon #0) for DotsForReferenceTaxon mode
    let ref_dna = ws
        .sequences
        .first()
        .map(|s| s.aligned_dna.as_slice())
        .unwrap_or(&[]);

    // Sequence Matrix Rows (Pixel-Locked `<text x="...">`)
    for (r_idx, seq) in ws.sequences.iter().take(n_rows).enumerate() {
        let y = matrix_y0 + r_idx as f64 * row_h;
        let row_bg = if r_idx % 2 == 0 { "#ffffff" } else { "#f8fafc" };
        svg.push_str(&format!(
            "  <rect x=\"10\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"{}\" stroke=\"#e2e8f0\" stroke-width=\"0.5\"/>\n",
            y,
            svg_w - 20.0,
            row_h,
            row_bg
        ));

        // Left Taxon Label + Clade pill
        let short_taxon: String = seq.taxon.chars().take(24).collect();
        let short_id: String = seq.id.chars().take(12).collect();
        svg.push_str(&format!(
            "  <text x=\"16\" y=\"{:.1}\" fill=\"#0f172a\" font-size=\"9.8\" font-weight=\"700\"><tspan fill=\"#0284c7\">{}</tspan> <tspan font-style=\"italic\">{}</tspan></text>\n",
            y + row_h * 0.64,
            xml_escape(&short_id),
            xml_escape(&short_taxon)
        ));
        svg.push_str(&format!(
            "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" fill=\"#64748b\" font-size=\"8.5\">{:.1}%</text>\n",
            left_label_w - 6.0,
            y + row_h * 0.64,
            seq.gc_percent
        ));

        // Residue Cells
        for d_idx in 0..display_cols {
            let col_idx = if cfg.mode == ViewAlphabetMode::AminoAcidCompact {
                s0 + d_idx * 3
            } else {
                s0 + d_idx
            };
            if col_idx >= ws.alignment_length {
                break;
            }
            let x = left_label_w + d_idx as f64 * col_w;
            let x_center = x + col_w * 0.5;
            let col_stat = ws.column_stats.get(col_idx);
            let raw_base = seq.aligned_dna.get(col_idx).copied().unwrap_or(b'-');
            let (codon_pos, aa_byte, _) =
                ws.codon_and_aa_at(r_idx, col_idx, cfg.genetic_code, cfg.reading_frame);

            match cfg.mode {
                ViewAlphabetMode::DualDnaAndAminoAcid => {
                    // Top half: Nucleotide
                    let disp_dna = match cfg.match_mode {
                        ConsensusMatchMode::DotsForConsensus
                            if col_stat.map(|c| c.consensus_dna) == Some(raw_base)
                                && raw_base != b'-' =>
                        {
                            b'.'
                        }
                        ConsensusMatchMode::DotsForReferenceTaxon
                            if r_idx > 0
                                && ref_dna.get(col_idx) == Some(&raw_base)
                                && raw_base != b'-' =>
                        {
                            b'.'
                        }
                        _ => raw_base,
                    };
                    let (bg_dna, fg_dna) = residue_cell_colors(
                        raw_base,
                        false,
                        ViewAlphabetMode::Nucleotide,
                        cfg.color_scheme,
                        col_stat,
                    );
                    svg.push_str(&format!(
                        "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"rgb({},{},{})\"/>\n",
                        x,
                        y + 0.5,
                        col_w - 0.5,
                        row_h * 0.52,
                        bg_dna.0,
                        bg_dna.1,
                        bg_dna.2
                    ));
                    svg.push_str(&format!(
                        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" fill=\"rgb({},{},{})\" font-size=\"8.8\" font-weight=\"700\">{}</text>\n",
                        x_center,
                        y + row_h * 0.38,
                        fg_dna.0,
                        fg_dna.1,
                        fg_dna.2,
                        disp_dna as char
                    ));

                    // Bottom half: Translated Codon AA spanned across the 3 bases of the codon
                    if codon_pos == 0 && col_idx + 2 < e0 {
                        let (bg_aa, fg_aa) = residue_cell_colors(
                            aa_byte,
                            true,
                            ViewAlphabetMode::AminoAcidOverlay,
                            cfg.color_scheme,
                            col_stat,
                        );
                        svg.push_str(&format!(
                            "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"2\" fill=\"rgb({},{},{})\" stroke=\"#0f172a\" stroke-width=\"0.4\"/>\n",
                            x,
                            y + row_h * 0.54,
                            col_w * 3.0 - 0.8,
                            row_h * 0.44,
                            bg_aa.0,
                            bg_aa.1,
                            bg_aa.2
                        ));
                        svg.push_str(&format!(
                            "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" fill=\"rgb({},{},{})\" font-size=\"8.6\" font-weight=\"700\">{}</text>\n",
                            x + col_w * 1.5,
                            y + row_h * 0.86,
                            fg_aa.0,
                            fg_aa.1,
                            fg_aa.2,
                            aa_byte as char
                        ));
                    }
                }
                ViewAlphabetMode::AminoAcidOverlay => {
                    let (bg_aa, fg_aa) = residue_cell_colors(
                        aa_byte,
                        true,
                        cfg.mode,
                        cfg.color_scheme,
                        col_stat,
                    );
                    svg.push_str(&format!(
                        "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"rgb({},{},{})\"/>\n",
                        x,
                        y + 0.5,
                        col_w - 0.5,
                        row_h - 1.0,
                        bg_aa.0,
                        bg_aa.1,
                        bg_aa.2
                    ));
                    if codon_pos == 1 {
                        svg.push_str(&format!(
                            "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" fill=\"rgb({},{},{})\" font-size=\"10\" font-weight=\"700\">{}</text>\n",
                            x_center,
                            y + row_h * 0.68,
                            fg_aa.0,
                            fg_aa.1,
                            fg_aa.2,
                            aa_byte as char
                        ));
                    }
                }
                ViewAlphabetMode::AminoAcidCompact => {
                    let (bg_aa, fg_aa) = residue_cell_colors(
                        aa_byte,
                        true,
                        cfg.mode,
                        cfg.color_scheme,
                        col_stat,
                    );
                    svg.push_str(&format!(
                        "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"rgb({},{},{})\"/>\n",
                        x,
                        y + 0.5,
                        col_w - 0.6,
                        row_h - 1.0,
                        bg_aa.0,
                        bg_aa.1,
                        bg_aa.2
                    ));
                    svg.push_str(&format!(
                        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" fill=\"rgb({},{},{})\" font-size=\"10.2\" font-weight=\"700\">{}</text>\n",
                        x_center,
                        y + row_h * 0.68,
                        fg_aa.0,
                        fg_aa.1,
                        fg_aa.2,
                        aa_byte as char
                    ));
                }
                ViewAlphabetMode::Ry3Recoded => {
                    let ry = recode_ry3(raw_base);
                    let (bg, fg) = residue_cell_colors(ry, false, cfg.mode, cfg.color_scheme, col_stat);
                    svg.push_str(&format!(
                        "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"rgb({},{},{})\"/>\n",
                        x, y + 0.5, col_w - 0.5, row_h - 1.0, bg.0, bg.1, bg.2
                    ));
                    svg.push_str(&format!(
                        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" fill=\"rgb({},{},{})\" font-size=\"9.5\" font-weight=\"700\">{}</text>\n",
                        x_center, y + row_h * 0.68, fg.0, fg.1, fg.2, ry as char
                    ));
                }
                ViewAlphabetMode::Dayhoff6Recoded => {
                    let d6 = recode_dayhoff6(aa_byte);
                    let (bg, fg) = residue_cell_colors(d6, true, cfg.mode, cfg.color_scheme, col_stat);
                    svg.push_str(&format!(
                        "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"rgb({},{},{})\"/>\n",
                        x, y + 0.5, col_w - 0.5, row_h - 1.0, bg.0, bg.1, bg.2
                    ));
                    if codon_pos == 1 || ws.detected_protein_input {
                        svg.push_str(&format!(
                            "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" fill=\"rgb({},{},{})\" font-size=\"9.5\" font-weight=\"700\">{}</text>\n",
                            x_center, y + row_h * 0.68, fg.0, fg.1, fg.2, d6 as char
                        ));
                    }
                }
                ViewAlphabetMode::Nucleotide => {
                    let disp = match cfg.match_mode {
                        ConsensusMatchMode::DotsForConsensus
                            if col_stat.map(|c| c.consensus_dna) == Some(raw_base)
                                && raw_base != b'-' =>
                        {
                            b'.'
                        }
                        ConsensusMatchMode::DotsForReferenceTaxon
                            if r_idx > 0
                                && ref_dna.get(col_idx) == Some(&raw_base)
                                && raw_base != b'-' =>
                        {
                            b'.'
                        }
                        _ => raw_base,
                    };
                    let (bg, fg) =
                        residue_cell_colors(raw_base, false, cfg.mode, cfg.color_scheme, col_stat);
                    svg.push_str(&format!(
                        "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"rgb({},{},{})\"/>\n",
                        x, y + 0.5, col_w - 0.5, row_h - 1.0, bg.0, bg.1, bg.2
                    ));
                    svg.push_str(&format!(
                        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" fill=\"rgb({},{},{})\" font-size=\"9.5\" font-weight=\"700\">{}</text>\n",
                        x_center, y + row_h * 0.68, fg.0, fg.1, fg.2, disp as char
                    ));
                }
            }
        }

        // Right-side clade badge
        let right_x = left_label_w + matrix_w + 8.0;
        let short_clade: String = seq.clade_or_family.chars().take(12).collect();
        svg.push_str(&format!(
            "  <text x=\"{:.1}\" y=\"{:.1}\" fill=\"#475569\" font-size=\"8.8\" font-weight=\"600\">{}</text>\n",
            right_x,
            y + row_h * 0.64,
            xml_escape(&short_clade)
        ));
    }

    // Majority-Rule Consensus Bar at Bottom of Matrix
    let cons_row_y = matrix_y0 + matrix_h + 4.0;
    svg.push_str(&format!(
        "  <rect x=\"10\" y=\"{:.1}\" width=\"{:.1}\" height=\"36\" rx=\"5\" fill=\"#0f172a\"/>\n",
        cons_row_y,
        svg_w - 20.0
    ));
    svg.push_str(&format!(
        "  <text x=\"18\" y=\"{:.1}\" fill=\"#38bdf8\" font-size=\"9.8\" font-weight=\"700\">Majority Consensus (DNA / AA)</text>\n",
        cons_row_y + 22.0
    ));
    for d_idx in 0..display_cols {
        let col_idx = if cfg.mode == ViewAlphabetMode::AminoAcidCompact {
            s0 + d_idx * 3
        } else {
            s0 + d_idx
        };
        if let Some(st) = ws.column_stats.get(col_idx) {
            let x_center = left_label_w + d_idx as f64 * col_w + col_w * 0.5;
            if cfg.mode == ViewAlphabetMode::AminoAcidCompact {
                svg.push_str(&format!(
                    "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" fill=\"#f8fafc\" font-size=\"10\" font-weight=\"700\">{}</text>\n",
                    x_center,
                    cons_row_y + 22.0,
                    st.consensus_aa as char
                ));
            } else {
                svg.push_str(&format!(
                    "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" fill=\"#38bdf8\" font-size=\"9\" font-weight=\"700\">{}</text>\n",
                    x_center,
                    cons_row_y + 15.0,
                    st.consensus_dna as char
                ));
                if (col_idx.saturating_sub(cfg.reading_frame.saturating_sub(1) % 3)) % 3 == 1 {
                    svg.push_str(&format!(
                        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" fill=\"#fbbf24\" font-size=\"8.8\" font-weight=\"700\">{}</text>\n",
                        x_center,
                        cons_row_y + 29.0,
                        st.consensus_aa as char
                    ));
                }
            }
        }
    }

    // Bottom Color Key Legend
    let leg_y = cons_row_y + 44.0;
    svg.push_str(&format!(
        "  <text x=\"18\" y=\"{:.1}\" fill=\"#334155\" font-size=\"9.2\" font-weight=\"700\">AliView / RuNIE Residue Color Key:</text>\n",
        leg_y + 16.0
    ));
    let legend_items: &[(&str, &str)] = &[
        ("A (Adenine)", "#22c55e"),
        ("C (Cytosine)", "#3b82f6"),
        ("G (Guanine)", "#f59e0b"),
        ("T/U (Thymine)", "#ef4444"),
        ("AA Hydrophobic (AILMFWV)", "#2563eb"),
        ("AA Basic (KR)", "#dc2626"),
        ("AA Acidic (DE)", "#9333ea"),
        ("AA Polar (NQST)", "#10b981"),
        ("Cys (C)", "#eab308"),
    ];
    let mut lx = 235.0;
    for (lbl, col_hex) in legend_items {
        svg.push_str(&format!(
            "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"11\" height=\"11\" rx=\"2\" fill=\"{}\"/>\n",
            lx,
            leg_y + 6.5,
            col_hex
        ));
        svg.push_str(&format!(
            "  <text x=\"{:.1}\" y=\"{:.1}\" fill=\"#334155\" font-size=\"8.6\" font-weight=\"600\">{}</text>\n",
            lx + 15.0,
            leg_y + 15.5,
            lbl
        ));
        lx += (lbl.len() as f64) * 5.6 + 26.0;
    }

    svg.push_str("</svg>\n");
    svg
}
/// Pure-Rust 32-bit BGRA software framebuffer used for both the live Win32 interactive GUI
/// (`StretchDIBits`) and headless raster PNG snapshot generation (`--snapshot-png`).
pub struct SoftwareFrameBuffer {
    pub width: usize,
    pub height: usize,
    /// Pixels stored as `0x00RRGGBB` (`u32`, compatible with Win32 `BI_RGB` top-down 32-bit DIB)
    pub pixels: Vec<u32>,
}

impl SoftwareFrameBuffer {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width: width.max(1),
            height: height.max(1),
            pixels: vec![0x000B0F19; width.max(1) * height.max(1)],
        }
    }

    pub fn resize(&mut self, width: usize, height: usize) {
        let w = width.max(1);
        let h = height.max(1);
        if self.width != w || self.height != h {
            self.width = w;
            self.height = h;
            self.pixels.resize(w * h, 0x000B0F19);
        }
    }

    #[inline]
    pub fn pack_rgb(r: u8, g: u8, b: u8) -> u32 {
        ((r as u32) << 16) | ((g as u32) << 8) | (b as u32)
    }

    pub fn clear(&mut self, rgb: (u8, u8, u8)) {
        let c = Self::pack_rgb(rgb.0, rgb.1, rgb.2);
        self.pixels.fill(c);
    }

    pub fn fill_rect(&mut self, x: i32, y: i32, w: i32, h: i32, rgb: (u8, u8, u8)) {
        if w <= 0 || h <= 0 {
            return;
        }
        let x0 = x.max(0) as usize;
        let y0 = y.max(0) as usize;
        let x1 = (x + w).max(0).min(self.width as i32) as usize;
        let y1 = (y + h).max(0).min(self.height as i32) as usize;
        if x0 >= x1 || y0 >= y1 {
            return;
        }
        let c = Self::pack_rgb(rgb.0, rgb.1, rgb.2);
        for row in y0..y1 {
            let base = row * self.width;
            self.pixels[base + x0..base + x1].fill(c);
        }
    }

    pub fn stroke_rect(&mut self, x: i32, y: i32, w: i32, h: i32, rgb: (u8, u8, u8)) {
        self.fill_rect(x, y, w, 1, rgb);
        self.fill_rect(x, y + h - 1, w, 1, rgb);
        self.fill_rect(x, y, 1, h, rgb);
        self.fill_rect(x + w - 1, y, 1, h, rgb);
    }

    /// Draw a single ASCII character at `(x, y)` using our crisp built-in 5x7 bitmap font scaled by `scale` (1 or 2).
    pub fn draw_char(&mut self, x: i32, y: i32, ch: char, rgb: (u8, u8, u8), scale: i32) {
        let glyph = glyph_5x7(ch);
        let c = Self::pack_rgb(rgb.0, rgb.1, rgb.2);
        let s = scale.max(1);
        for (row_idx, &row_bits) in glyph.iter().enumerate() {
            for col_idx in 0..5 {
                if (row_bits >> (4 - col_idx)) & 1 == 1 {
                    for dy in 0..s {
                        let py = y + (row_idx as i32) * s + dy;
                        if py < 0 || py >= self.height as i32 {
                            continue;
                        }
                        for dx in 0..s {
                            let px = x + (col_idx as i32) * s + dx;
                            if px >= 0 && px < self.width as i32 {
                                self.pixels[(py as usize) * self.width + (px as usize)] = c;
                            }
                        }
                    }
                }
            }
        }
    }

    /// Draw a text string at `(x, y)` with `scale` (`1` = 6x8 px per char, `2` = 12x16 px per char),
    /// optionally truncated to `max_px_width`.
    pub fn draw_text(
        &mut self,
        x: i32,
        y: i32,
        text: &str,
        rgb: (u8, u8, u8),
        scale: i32,
        max_px_width: Option<i32>,
    ) {
        let s = scale.max(1);
        let step = 6 * s;
        let mut cx = x;
        let limit_x = max_px_width.map(|w| x + w).unwrap_or(i32::MAX);
        for ch in text.chars() {
            if cx + 5 * s > limit_x {
                break;
            }
            self.draw_char(cx, y, ch, rgb, s);
            cx += step;
        }
    }

    /// Encode this framebuffer as a valid, standards-compliant 24-bit RGB PNG byte vector
    /// using pure Rust (uncompressed zlib blocks + Adler-32 + IEEE CRC-32).
    pub fn encode_png(&self) -> Vec<u8> {
        let w = self.width as u32;
        let h = self.height as u32;

        // Build raw PNG scanlines: each row has 1 filter byte (0) + 3 * width RGB bytes
        let row_len = 1 + self.width * 3;
        let mut raw_data = Vec::with_capacity(row_len * self.height);
        for y in 0..self.height {
            raw_data.push(0u8); // Filter type 0 (None)
            let base = y * self.width;
            for x in 0..self.width {
                let p = self.pixels[base + x];
                raw_data.push(((p >> 16) & 0xFF) as u8);
                raw_data.push(((p >> 8) & 0xFF) as u8);
                raw_data.push((p & 0xFF) as u8);
            }
        }

        // Wrap raw_data in a valid Zlib stream using uncompressed DEFLATE stored blocks (BTYPE=00)
        let mut zlib_stream = Vec::with_capacity(raw_data.len() + (raw_data.len() / 65535 + 2) * 6 + 6);
        zlib_stream.push(0x78); // CMF: CM=8 (Deflate), CINFO=7 (32K window)
        zlib_stream.push(0x01); // FLG: FCHECK so 0x7801 % 31 == 0

        let chunks: Vec<&[u8]> = raw_data.chunks(65_535).collect();
        for (i, chunk) in chunks.iter().enumerate() {
            let is_last = i + 1 == chunks.len();
            zlib_stream.push(if is_last { 0x01 } else { 0x00 });
            let len = chunk.len() as u16;
            let nlen = !len;
            zlib_stream.extend_from_slice(&len.to_le_bytes());
            zlib_stream.extend_from_slice(&nlen.to_le_bytes());
            zlib_stream.extend_from_slice(chunk);
        }

        // Adler-32 checksum of raw_data
        let mut s1: u32 = 1;
        let mut s2: u32 = 0;
        for &b in &raw_data {
            s1 = (s1 + b as u32) % 65_521;
            s2 = (s2 + s1) % 65_521;
        }
        let adler = (s2 << 16) | s1;
        zlib_stream.extend_from_slice(&adler.to_be_bytes());

        let mut png = Vec::with_capacity(zlib_stream.len() + 128);
        // PNG Signature
        png.extend_from_slice(&[137, 80, 78, 71, 13, 10, 26, 10]);

        // IHDR Chunk (13 bytes)
        let mut ihdr = Vec::with_capacity(13);
        ihdr.extend_from_slice(&w.to_be_bytes());
        ihdr.extend_from_slice(&h.to_be_bytes());
        ihdr.push(8); // Bit depth = 8
        ihdr.push(2); // Color type = 2 (Truecolor RGB)
        ihdr.push(0); // Compression = 0
        ihdr.push(0); // Filter = 0
        ihdr.push(0); // Interlace = 0
        write_png_chunk(&mut png, b"IHDR", &ihdr);

        // IDAT Chunk
        write_png_chunk(&mut png, b"IDAT", &zlib_stream);

        // IEND Chunk
        write_png_chunk(&mut png, b"IEND", &[]);
        png
    }
}

fn write_png_chunk(out: &mut Vec<u8>, chunk_type: &[u8; 4], data: &[u8]) {
    let len = data.len() as u32;
    out.extend_from_slice(&len.to_be_bytes());
    out.extend_from_slice(chunk_type);
    out.extend_from_slice(data);
    let mut crc: u32 = 0xFFFF_FFFF;
    for &b in chunk_type.iter().chain(data.iter()) {
        crc ^= b as u32;
        for _ in 0..8 {
            if crc & 1 != 0 {
                crc = (crc >> 1) ^ 0xEDB8_8320;
            } else {
                crc >>= 1;
            }
        }
    }
    out.extend_from_slice(&(!crc).to_be_bytes());
}

/// Crisp 5x7 pixel bitmap font for all printable ASCII characters.
fn glyph_5x7(ch: char) -> [u8; 7] {
    match ch {
        'A' | 'a' => [0b01110, 0b10001, 0b10001, 0b11111, 0b10001, 0b10001, 0b10001],
        'B' | 'b' => [0b11110, 0b10001, 0b10001, 0b11110, 0b10001, 0b10001, 0b11110],
        'C' | 'c' => [0b01110, 0b10001, 0b10000, 0b10000, 0b10000, 0b10001, 0b01110],
        'D' | 'd' => [0b11100, 0b10010, 0b10001, 0b10001, 0b10001, 0b10010, 0b11100],
        'E' | 'e' => [0b11111, 0b10000, 0b10000, 0b11110, 0b10000, 0b10000, 0b11111],
        'F' | 'f' => [0b11111, 0b10000, 0b10000, 0b11110, 0b10000, 0b10000, 0b10000],
        'G' | 'g' => [0b01110, 0b10001, 0b10000, 0b10111, 0b10001, 0b10001, 0b01110],
        'H' | 'h' => [0b10001, 0b10001, 0b10001, 0b11111, 0b10001, 0b10001, 0b10001],
        'I' | 'i' => [0b01110, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110],
        'J' | 'j' => [0b00111, 0b00010, 0b00010, 0b00010, 0b10010, 0b10010, 0b01100],
        'K' | 'k' => [0b10001, 0b10010, 0b10100, 0b11000, 0b10100, 0b10010, 0b10001],
        'L' | 'l' => [0b10000, 0b10000, 0b10000, 0b10000, 0b10000, 0b10000, 0b11111],
        'M' | 'm' => [0b10001, 0b11011, 0b10101, 0b10101, 0b10001, 0b10001, 0b10001],
        'N' | 'n' => [0b10001, 0b11001, 0b10101, 0b10011, 0b10001, 0b10001, 0b10001],
        'O' | 'o' => [0b01110, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01110],
        'P' | 'p' => [0b11110, 0b10001, 0b10001, 0b11110, 0b10000, 0b10000, 0b10000],
        'Q' | 'q' => [0b01110, 0b10001, 0b10001, 0b10001, 0b10101, 0b10010, 0b01101],
        'R' | 'r' => [0b11110, 0b10001, 0b10001, 0b11110, 0b10100, 0b10010, 0b10001],
        'S' | 's' => [0b01111, 0b10000, 0b10000, 0b01110, 0b00001, 0b00001, 0b11110],
        'T' | 't' => [0b11111, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100],
        'U' | 'u' => [0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01110],
        'V' | 'v' => [0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01010, 0b00100],
        'W' | 'w' => [0b10001, 0b10001, 0b10001, 0b10101, 0b10101, 0b10101, 0b01010],
        'X' | 'x' => [0b10001, 0b10001, 0b01010, 0b00100, 0b01010, 0b10001, 0b10001],
        'Y' | 'y' => [0b10001, 0b10001, 0b01010, 0b00100, 0b00100, 0b00100, 0b00100],
        'Z' | 'z' => [0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b10000, 0b11111],
        '0' => [0b01110, 0b10001, 0b10011, 0b10101, 0b11001, 0b10001, 0b01110],
        '1' => [0b00100, 0b01100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110],
        '2' => [0b01110, 0b10001, 0b00001, 0b00010, 0b00100, 0b01000, 0b11111],
        '3' => [0b11110, 0b00001, 0b00001, 0b01110, 0b00001, 0b00001, 0b11110],
        '4' => [0b00010, 0b00110, 0b01010, 0b10010, 0b11111, 0b00010, 0b00010],
        '5' => [0b11111, 0b10000, 0b10000, 0b11110, 0b00001, 0b00001, 0b11110],
        '6' => [0b01110, 0b10000, 0b10000, 0b11110, 0b10001, 0b10001, 0b01110],
        '7' => [0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b01000, 0b01000],
        '8' => [0b01110, 0b10001, 0b10001, 0b01110, 0b10001, 0b10001, 0b01110],
        '9' => [0b01110, 0b10001, 0b10001, 0b01111, 0b00001, 0b00001, 0b01110],
        '-' => [0b00000, 0b00000, 0b00000, 0b11111, 0b00000, 0b00000, 0b00000],
        '+' => [0b00000, 0b00100, 0b00100, 0b11111, 0b00100, 0b00100, 0b00000],
        '.' => [0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b01100, 0b01100],
        ':' => [0b00000, 0b01100, 0b01100, 0b00000, 0b01100, 0b01100, 0b00000],
        '/' => [0b00001, 0b00010, 0b00010, 0b00100, 0b01000, 0b01000, 0b10000],
        '|' => [0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100],
        '(' | '[' => [0b00110, 0b01000, 0b01000, 0b01000, 0b01000, 0b01000, 0b00110],
        ')' | ']' => [0b01100, 0b00010, 0b00010, 0b00010, 0b00010, 0b00010, 0b01100],
        '%' => [0b11001, 0b11010, 0b00010, 0b00100, 0b01000, 0b01011, 0b10011],
        '*' => [0b00000, 0b10101, 0b01110, 0b11111, 0b01110, 0b10101, 0b00000],
        '_' => [0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b11111],
        '=' => [0b00000, 0b00000, 0b11111, 0b00000, 0b11111, 0b00000, 0b00000],
        '>' => [0b10000, 0b01000, 0b00100, 0b00010, 0b00100, 0b01000, 0b10000],
        '<' => [0b00001, 0b00010, 0b00100, 0b01000, 0b00100, 0b00010, 0b00001],
        '#' => [0b01010, 0b11111, 0b01010, 0b01010, 0b11111, 0b01010, 0b00000],
        _ => [0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b00000],
    }
}
/// Interactive state for the AliView-modeled desktop GUI and raster PNG renderer.
#[derive(Debug, Clone)]
pub struct AlignmentViewerState {
    pub workspace: AlignmentWorkspace,
    pub mode: ViewAlphabetMode,
    pub color_scheme: AlignmentColorScheme,
    pub match_mode: ConsensusMatchMode,
    pub genetic_code: GeneticCodeTable,
    pub reading_frame: usize,
    pub scroll_col: usize,
    pub scroll_row: usize,
    pub col_width_px: i32, // Zoom level: 1..28 px per column
    pub row_height_px: i32,
    pub show_feature_sidebar: bool,
    pub selected_feature_idx: Option<usize>,
    pub selected_col_range_1based: Option<(usize, usize)>,
    pub hover_col_1based: Option<usize>,
    pub hover_row_idx: Option<usize>,
    pub search_query: String,
    pub search_hits: Vec<(usize, usize, usize, String)>,
    pub active_search_hit_idx: usize,
    pub status_message: String,
    pub dragging_minimap: bool,
    pub dragging_selection_start_col: Option<usize>,
}

#[derive(Debug, Clone)]
pub struct UiHitRect {
    pub x: i32,
    pub y: i32,
    pub w: i32,
    pub h: i32,
    pub action_id: u16,
}

pub const ACT_MODE_DNA: u16 = 101;
pub const ACT_MODE_AA_OVERLAY: u16 = 102;
pub const ACT_MODE_AA_COMPACT: u16 = 103;
pub const ACT_MODE_DUAL: u16 = 104;
pub const ACT_MODE_RY3: u16 = 105;
pub const ACT_MODE_DAYHOFF6: u16 = 106;

pub const ACT_COLOR_ALIVIEW: u16 = 201;
pub const ACT_COLOR_CLUSTAL: u16 = 202;
pub const ACT_COLOR_HYDRO: u16 = 203;
pub const ACT_COLOR_RUNIE_PI: u16 = 204;
pub const ACT_COLOR_OKABE: u16 = 205;

pub const ACT_MATCH_ALL: u16 = 301;
pub const ACT_MATCH_DOTS_CONS: u16 = 302;
pub const ACT_MATCH_DOTS_REF: u16 = 303;

pub const ACT_FRAME_1: u16 = 401;
pub const ACT_FRAME_2: u16 = 402;
pub const ACT_FRAME_3: u16 = 403;
pub const ACT_REV_COMP: u16 = 404;

pub const ACT_GCODE_1: u16 = 411;
pub const ACT_GCODE_2: u16 = 412;
pub const ACT_GCODE_5: u16 = 413;
pub const ACT_GCODE_11: u16 = 414;

pub const ACT_ZOOM_IN: u16 = 501;
pub const ACT_ZOOM_OUT: u16 = 502;
pub const ACT_ZOOM_OVERVIEW: u16 = 503;
pub const ACT_ZOOM_STANDARD: u16 = 504;
pub const ACT_ZOOM_FOCAL_BEST: u16 = 505;

pub const ACT_TOGGLE_SIDEBAR: u16 = 601;
pub const ACT_RESCAN_FEATURES: u16 = 602;
pub const ACT_MASK_ADNA: u16 = 603;
pub const ACT_STRIP_GAPS: u16 = 604;

pub const ACT_SNAPSHOT_SVG_PNG: u16 = 701;
pub const ACT_OPEN_FILE: u16 = 801;
pub const ACT_SAVE_FASTA: u16 = 802;
pub const ACT_LOAD_COI: u16 = 811;
pub const ACT_LOAD_CHR1: u16 = 812;
pub const ACT_LOAD_HYMEN: u16 = 813;
pub const ACT_EXIT: u16 = 999;

pub const ACT_FEAT_BASE: u16 = 2000;

impl AlignmentViewerState {
    pub fn new(workspace: AlignmentWorkspace) -> Self {
        let len = workspace.alignment_length.max(1);
        Self {
            workspace,
            mode: ViewAlphabetMode::Nucleotide,
            color_scheme: AlignmentColorScheme::AliViewDefault,
            match_mode: ConsensusMatchMode::ShowAllCharacters,
            genetic_code: GeneticCodeTable::VertebrateMito2,
            reading_frame: 1,
            scroll_col: 0,
            scroll_row: 0,
            col_width_px: 14,
            row_height_px: 20,
            show_feature_sidebar: true,
            selected_feature_idx: Some(0),
            selected_col_range_1based: None,
            hover_col_1based: None,
            hover_row_idx: None,
            search_query: String::new(),
            search_hits: Vec::new(),
            active_search_hit_idx: 0,
            status_message: format!(
                "Loaded alignment ({} bp). Use Menu/Toolbar or Keys: [1-6] Mode | [+/-] Zoom | [0] Bird's-Eye | [F] Best Focal | [F12] Snapshot SVG+PNG",
                len
            ),
            dragging_minimap: false,
            dragging_selection_start_col: None,
        }
    }

    pub fn visible_cols_count(&self, canvas_w: i32) -> usize {
        let sidebar_w = if self.show_feature_sidebar { 245 } else { 0 };
        let label_w = 215;
        let matrix_w = (canvas_w - sidebar_w - label_w - 16).max(64);
        let step_px = self.col_width_px.max(1);
        let cols = (matrix_w / step_px).max(1) as usize;
        if self.mode == ViewAlphabetMode::AminoAcidCompact {
            cols * 3
        } else {
            cols
        }
    }

    pub fn clamp_scroll(&mut self, canvas_w: i32, canvas_h: i32) {
        let vis_cols = self.visible_cols_count(canvas_w);
        let max_scroll_col = self.workspace.alignment_length.saturating_sub(vis_cols / 2);
        if self.scroll_col > max_scroll_col {
            self.scroll_col = max_scroll_col;
        }
        let vis_rows = ((canvas_h - 210).max(40) / self.row_height_px.max(12)) as usize;
        let max_scroll_row = self.workspace.sequences.len().saturating_sub(vis_rows.max(1));
        if self.scroll_row > max_scroll_row {
            self.scroll_row = max_scroll_row;
        }
    }

    pub fn select_and_zoom_to_feature(&mut self, feat_idx: usize, canvas_w: i32) {
        if let Some(f) = self.workspace.features.get(feat_idx) {
            self.selected_feature_idx = Some(feat_idx);
            let s = f.start_col_1based;
            let e = f.end_col_1based;
            self.selected_col_range_1based = Some((s, e));
            let span = (e.saturating_sub(s) + 1).max(1);
            let sidebar_w = if self.show_feature_sidebar { 245 } else { 0 };
            let avail_w = (canvas_w - sidebar_w - 235).max(200);
            let target_col_px = (avail_w / (span as i32 + 6)).clamp(1, 22);
            self.col_width_px = target_col_px;
            self.scroll_col = s.saturating_sub(4);
            self.status_message = format!(
                "Selected Feature [{}]: {} ({}..{} bp, {} columns) — Press F12 to Snapshot SVG+PNG",
                f.feature_id, f.label, s, e, span
            );
        }
    }

    pub fn execute_action(&mut self, action_id: u16, canvas_w: i32) {
        if action_id >= ACT_FEAT_BASE {
            let idx = (action_id - ACT_FEAT_BASE) as usize;
            self.select_and_zoom_to_feature(idx, canvas_w);
            return;
        }
        match action_id {
            ACT_MODE_DNA => {
                self.mode = ViewAlphabetMode::Nucleotide;
                self.row_height_px = 20;
                self.status_message = "Switched to DNA / Nucleotide View (AliView 4-Base Color Mode)".to_string();
            }
            ACT_MODE_AA_OVERLAY => {
                self.mode = ViewAlphabetMode::AminoAcidOverlay;
                self.row_height_px = 20;
                self.status_message = format!(
                    "Switched to Translated Amino Acid Overlay (Frame +{}, {})",
                    self.reading_frame,
                    self.genetic_code.label()
                );
            }
            ACT_MODE_AA_COMPACT => {
                self.mode = ViewAlphabetMode::AminoAcidCompact;
                self.row_height_px = 20;
                self.status_message = format!(
                    "Switched to Compact Amino Acid View (1 Column per Codon, Frame +{})",
                    self.reading_frame
                );
            }
            ACT_MODE_DUAL => {
                self.mode = ViewAlphabetMode::DualDnaAndAminoAcid;
                self.row_height_px = 30;
                self.status_message = "Switched to Dual DNA + Translated Amino Acid Split Track".to_string();
            }
            ACT_MODE_RY3 => {
                self.mode = ViewAlphabetMode::Ry3Recoded;
                self.row_height_px = 20;
                self.status_message = "Switched to RuNIE RY3 Purine (R) / Pyrimidine (Y) Recoding View".to_string();
            }
            ACT_MODE_DAYHOFF6 => {
                self.mode = ViewAlphabetMode::Dayhoff6Recoded;
                self.row_height_px = 20;
                self.status_message = "Switched to RuNIE Dayhoff-6 Functional State Recoding View (1..6)".to_string();
            }
            ACT_COLOR_ALIVIEW => {
                self.color_scheme = AlignmentColorScheme::AliViewDefault;
                self.status_message = "Color Palette: AliView Default".to_string();
            }
            ACT_COLOR_CLUSTAL => {
                self.color_scheme = AlignmentColorScheme::ClustalX;
                self.status_message = "Color Palette: ClustalX Chemistry".to_string();
            }
            ACT_COLOR_HYDRO => {
                self.color_scheme = AlignmentColorScheme::Hydrophobicity;
                self.status_message = "Color Palette: Kyte-Doolittle Hydrophobicity (Red=Hydrophobic, Blue=Hydrophilic)".to_string();
            }
            ACT_COLOR_RUNIE_PI => {
                self.color_scheme = AlignmentColorScheme::RunieConservation;
                self.status_message = "Color Palette: RuNIE Parsimony-Informative & Conservation Highlighter".to_string();
            }
            ACT_COLOR_OKABE => {
                self.color_scheme = AlignmentColorScheme::OkabeItoSafe;
                self.status_message = "Color Palette: Okabe-Ito Color-Blind Safe".to_string();
            }
            ACT_MATCH_ALL => {
                self.match_mode = ConsensusMatchMode::ShowAllCharacters;
                self.status_message = "Display Mode: Show All Characters".to_string();
            }
            ACT_MATCH_DOTS_CONS => {
                self.match_mode = if self.match_mode == ConsensusMatchMode::DotsForConsensus {
                    ConsensusMatchMode::ShowAllCharacters
                } else {
                    ConsensusMatchMode::DotsForConsensus
                };
                self.status_message = format!("Display Mode: {}", self.match_mode.label());
            }
            ACT_MATCH_DOTS_REF => {
                self.match_mode = if self.match_mode == ConsensusMatchMode::DotsForReferenceTaxon {
                    ConsensusMatchMode::ShowAllCharacters
                } else {
                    ConsensusMatchMode::DotsForReferenceTaxon
                };
                self.status_message = format!("Display Mode: {}", self.match_mode.label());
            }
            ACT_FRAME_1 | ACT_FRAME_2 | ACT_FRAME_3 => {
                self.reading_frame = (action_id - ACT_FRAME_1 + 1) as usize;
                self.workspace
                    .recompute_telemetry_and_auto_features(self.genetic_code, self.reading_frame);
                self.status_message = format!("Reading Frame set to +{}", self.reading_frame);
            }
            ACT_REV_COMP => {
                self.workspace
                    .reverse_complement_all(self.genetic_code, self.reading_frame);
                self.status_message = "Reverse-complemented all sequences in workspace".to_string();
            }
            ACT_GCODE_1 => {
                self.genetic_code = GeneticCodeTable::Standard1;
                self.workspace
                    .recompute_telemetry_and_auto_features(self.genetic_code, self.reading_frame);
                self.status_message = "Genetic Code set to Table 1 (Standard Nuclear)".to_string();
            }
            ACT_GCODE_2 => {
                self.genetic_code = GeneticCodeTable::VertebrateMito2;
                self.workspace
                    .recompute_telemetry_and_auto_features(self.genetic_code, self.reading_frame);
                self.status_message = "Genetic Code set to Table 2 (Vertebrate Mitochondrial)".to_string();
            }
            ACT_GCODE_5 => {
                self.genetic_code = GeneticCodeTable::InvertebrateMito5;
                self.workspace
                    .recompute_telemetry_and_auto_features(self.genetic_code, self.reading_frame);
                self.status_message = "Genetic Code set to Table 5 (Invertebrate Mitochondrial)".to_string();
            }
            ACT_GCODE_11 => {
                self.genetic_code = GeneticCodeTable::BacterialPlastid11;
                self.workspace
                    .recompute_telemetry_and_auto_features(self.genetic_code, self.reading_frame);
                self.status_message = "Genetic Code set to Table 11 (Bacterial / Plastid)".to_string();
            }
            ACT_ZOOM_IN => {
                self.col_width_px = (self.col_width_px + 2).min(26);
                self.status_message = format!("Zoomed In: {} px/column", self.col_width_px);
            }
            ACT_ZOOM_OUT => {
                self.col_width_px = (self.col_width_px - 2).max(1);
                self.status_message = format!("Zoomed Out: {} px/column", self.col_width_px);
            }
            ACT_ZOOM_OVERVIEW => {
                self.col_width_px = 2;
                self.scroll_col = 0;
                self.status_message = "Bird's-Eye Whole-Alignment Overview (2 px/col — scan entire locus)".to_string();
            }
            ACT_ZOOM_STANDARD => {
                self.col_width_px = 14;
                self.status_message = "Standard Residue Reading View (14 px/col)".to_string();
            }
            ACT_ZOOM_FOCAL_BEST => {
                let (s, e) = self.workspace.select_best_focal_window(60);
                self.selected_col_range_1based = Some((s, e));
                self.scroll_col = s.saturating_sub(1);
                self.col_width_px = 16;
                self.status_message = format!(
                    "Zoomed to Highest-Information Focal Motif Window: {}..{} bp (ready for paper snapshot F12)",
                    s, e
                );
            }
            ACT_TOGGLE_SIDEBAR => {
                self.show_feature_sidebar = !self.show_feature_sidebar;
            }
            ACT_RESCAN_FEATURES => {
                self.workspace.features = self
                    .workspace
                    .auto_detect_genomic_features(self.genetic_code, self.reading_frame);
                self.status_message = format!(
                    "Rescanned alignment: detected {} functional features & diagnostic windows",
                    self.workspace.features.len()
                );
            }
            ACT_MASK_ADNA => {
                let masked = self
                    .workspace
                    .apply_adna_deamination_mask(self.genetic_code, self.reading_frame);
                self.status_message = format!(
                    "RuNIE aDNA Deamination Filter: masked {} terminal 5' C->T / 3' G->A sites to N",
                    masked
                );
            }
            ACT_STRIP_GAPS => {
                let stripped = self.workspace.strip_gap_heavy_columns(
                    0.50,
                    self.genetic_code,
                    self.reading_frame,
                );
                self.status_message = format!(
                    "Stripped {} gap-heavy columns (>50% gaps); alignment length now {} bp",
                    stripped, self.workspace.alignment_length
                );
            }
            ACT_SNAPSHOT_SVG_PNG => {
                match self.capture_current_viewport_snapshots(canvas_w, None, None) {
                    Ok((svg_p, png_p)) => {
                        self.status_message = format!(
                            "Saved Paper Snapshots -> {} & {}",
                            svg_p.display(),
                            png_p.display()
                        );
                    }
                    Err(e) => {
                        self.status_message = format!("Snapshot error: {}", e);
                    }
                }
            }
            _ => {}
        }
    }

    /// Save both publication SVG and high-DPI PNG snapshots of the currently selected region
    /// (or current visible viewport) to `runs/snapshots/` (or custom output paths).
    pub fn capture_current_viewport_snapshots(
        &self,
        canvas_w: i32,
        custom_svg_path: Option<&Path>,
        custom_png_path: Option<&Path>,
    ) -> Result<(PathBuf, PathBuf), String> {
        let (s1, e1) = if let Some((sel_s, sel_e)) = self.selected_col_range_1based {
            (sel_s.min(sel_e), sel_s.max(sel_e))
        } else {
            let vis = self.visible_cols_count(canvas_w).clamp(15, 90);
            let s = self.scroll_col + 1;
            let e = (s + vis - 1).min(self.workspace.alignment_length.max(1));
            (s, e)
        };

        let snap_cfg = SnapshotConfig {
            start_col_1based: s1,
            end_col_1based: e1,
            mode: self.mode,
            color_scheme: self.color_scheme,
            match_mode: self.match_mode,
            genetic_code: self.genetic_code,
            reading_frame: self.reading_frame,
            max_taxa_rows: 32,
            subtitle_note: String::new(),
        };

        let out_dir = PathBuf::from("runs/snapshots");
        let _ = fs::create_dir_all(&out_dir);

        let svg_path = custom_svg_path
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| out_dir.join(format!("alignment_snapshot_{}_{}.svg", s1, e1)));
        let png_path = custom_png_path
            .map(|p| p.to_path_buf())
            .unwrap_or_else(|| out_dir.join(format!("alignment_snapshot_{}_{}.png", s1, e1)));

        if let Some(p) = svg_path.parent() {
            let _ = fs::create_dir_all(p);
        }
        if let Some(p) = png_path.parent() {
            let _ = fs::create_dir_all(p);
        }

        let svg_str = render_alignment_snapshot_svg(&self.workspace, &snap_cfg);
        fs::write(&svg_path, svg_str).map_err(|e| e.to_string())?;

        // Render high-res raster PNG of the snapshot window
        let mut snap_state = self.clone();
        snap_state.scroll_col = s1.saturating_sub(1);
        snap_state.selected_col_range_1based = Some((s1, e1));
        snap_state.col_width_px = 16;
        let png_w = (260 + ((e1 - s1 + 1) as usize) * 16 + 260).clamp(1180, 1920);
        let png_h = (220 + self.workspace.sequences.len().min(32) * (snap_state.row_height_px as usize)).clamp(560, 1200);
        let mut fb = SoftwareFrameBuffer::new(png_w, png_h);
        let _ = snap_state.render_interactive_frame(&mut fb);
        let png_bytes = fb.encode_png();
        fs::write(&png_path, png_bytes).map_err(|e| e.to_string())?;

        Ok((svg_path, png_path))
    }

    /// Render the complete AliView-style interactive interface into `fb` and return clickable `UiHitRect`s.
    pub fn render_interactive_frame(&self, fb: &mut SoftwareFrameBuffer) -> Vec<UiHitRect> {
        let mut hits = Vec::with_capacity(64);
        let w = fb.width as i32;
        let h = fb.height as i32;

        // 1. Base background
        fb.clear((11, 15, 25));

        // 2. Top Ribbon Toolbar (y = 0..42)
        fb.fill_rect(0, 0, w, 42, (15, 23, 42));
        fb.fill_rect(0, 41, w, 1, (51, 65, 85));

        let toolbar_buttons: &[(&str, u16, bool, (u8, u8, u8))] = &[
            (
                "DNA",
                ACT_MODE_DNA,
                self.mode == ViewAlphabetMode::Nucleotide,
                (34, 197, 94),
            ),
            (
                "AA-3bp",
                ACT_MODE_AA_OVERLAY,
                self.mode == ViewAlphabetMode::AminoAcidOverlay,
                (59, 130, 246),
            ),
            (
                "AA-1col",
                ACT_MODE_AA_COMPACT,
                self.mode == ViewAlphabetMode::AminoAcidCompact,
                (99, 102, 241),
            ),
            (
                "Dual DNA+AA",
                ACT_MODE_DUAL,
                self.mode == ViewAlphabetMode::DualDnaAndAminoAcid,
                (16, 185, 129),
            ),
            (
                "RY3",
                ACT_MODE_RY3,
                self.mode == ViewAlphabetMode::Ry3Recoded,
                (6, 182, 212),
            ),
            (
                "Dayhoff6",
                ACT_MODE_DAYHOFF6,
                self.mode == ViewAlphabetMode::Dayhoff6Recoded,
                (245, 158, 11),
            ),
            (
                "Fr+1",
                ACT_FRAME_1,
                self.reading_frame == 1,
                (56, 189, 248),
            ),
            (
                "Fr+2",
                ACT_FRAME_2,
                self.reading_frame == 2,
                (56, 189, 248),
            ),
            (
                "Fr+3",
                ACT_FRAME_3,
                self.reading_frame == 3,
                (56, 189, 248),
            ),
            (
                "Dots(.)",
                ACT_MATCH_DOTS_CONS,
                self.match_mode == ConsensusMatchMode::DotsForConsensus,
                (168, 85, 247),
            ),
            ("Zoom-", ACT_ZOOM_OUT, false, (148, 163, 184)),
            ("Zoom+", ACT_ZOOM_IN, false, (148, 163, 184)),
            (
                "Overview",
                ACT_ZOOM_OVERVIEW,
                self.col_width_px <= 2,
                (56, 189, 248),
            ),
            ("BestFocal", ACT_ZOOM_FOCAL_BEST, false, (236, 72, 153)),
            ("Mask-aDNA", ACT_MASK_ADNA, false, (239, 68, 68)),
            ("Snap SVG+PNG", ACT_SNAPSHOT_SVG_PNG, false, (16, 185, 129)),
        ];

        let mut bx = 8;
        for &(label, act, is_active, accent) in toolbar_buttons {
            let bw = (label.len() as i32) * 6 + 14;
            if bx + bw > w - 6 {
                break;
            }
            let bg = if is_active {
                accent
            } else {
                (30, 41, 59)
            };
            let fg = if is_active {
                (11, 15, 25)
            } else {
                (226, 232, 240)
            };
            fb.fill_rect(bx, 7, bw, 26, bg);
            fb.stroke_rect(bx, 7, bw, 26, accent);
            fb.draw_text(bx + 7, 16, label, fg, 1, Some(bw - 8));
            hits.push(UiHitRect {
                x: bx,
                y: 7,
                w: bw,
                h: 26,
                action_id: act,
            });
            bx += bw + 5;
        }

        // 3. Global Alignment Bird's-Eye Minimap & Viewport Lens (y = 42..78)
        let minimap_y = 42;
        let minimap_h = 36;
        fb.fill_rect(0, minimap_y, w, minimap_h, (9, 13, 22));
        fb.fill_rect(0, minimap_y + minimap_h - 1, w, 1, (30, 41, 59));
        fb.draw_text(
            10,
            minimap_y + 7,
            "GLOBAL LOCUS MAP:",
            (56, 189, 248),
            1,
            Some(120),
        );
        let len = self.workspace.alignment_length.max(1);
        fb.draw_text(
            10,
            minimap_y + 20,
            &format!("L={} bp | {} Taxa", len, self.workspace.sequences.len()),
            (148, 163, 184),
            1,
            Some(135),
        );

        let map_x0 = 150;
        let map_w = (w - map_x0 - 14).max(100);
        fb.fill_rect(map_x0, minimap_y + 5, map_w, 24, (15, 23, 42));
        fb.stroke_rect(map_x0, minimap_y + 5, map_w, 24, (51, 65, 85));

        // Draw conservation sparkline & feature blocks across the full alignment length
        for px in 0..map_w {
            let col = ((px as usize) * len) / (map_w as usize);
            if let Some(st) = self.workspace.column_stats.get(col) {
                let bh = (st.conservation_fraction * 18.0) as i32;
                let col_rgb = if st.is_parsimony_informative {
                    (168, 85, 247)
                } else if st.gap_fraction > 0.35 {
                    (239, 68, 68)
                } else {
                    (14, 165, 233)
                };
                fb.fill_rect(map_x0 + px, minimap_y + 27 - bh, 1, bh, col_rgb);
            }
        }
        // Highlight current visible viewport lens on minimap
        let vis_cols = self.visible_cols_count(w);
        let lens_x = map_x0 + ((self.scroll_col as i64 * map_w as i64) / len as i64) as i32;
        let lens_w = (((vis_cols as i64 * map_w as i64) / len as i64) as i32)
            .max(8)
            .min(map_w - (lens_x - map_x0));
        fb.stroke_rect(lens_x, minimap_y + 4, lens_w, 26, (52, 211, 153));
        fb.stroke_rect(lens_x + 1, minimap_y + 5, (lens_w - 2).max(2), 24, (52, 211, 153));

        // 4. Left Feature / Annotation Selector Sidebar (if visible)
        let content_y0 = 78;
        let status_h = 28;
        let consensus_bar_h = 34;
        let content_h = (h - content_y0 - status_h).max(80);
        let sidebar_w = if self.show_feature_sidebar { 245 } else { 0 };

        if self.show_feature_sidebar {
            fb.fill_rect(0, content_y0, sidebar_w, content_h, (13, 19, 33));
            fb.fill_rect(sidebar_w - 1, content_y0, 1, content_h, (51, 65, 85));
            fb.draw_text(
                8,
                content_y0 + 8,
                "FEATURES & ANNOTATIONS (CLICK)",
                (56, 189, 248),
                1,
                Some(sidebar_w - 12),
            );

            let mut fy = content_y0 + 24;
            for (f_idx, feat) in self.workspace.features.iter().enumerate() {
                if fy + 28 > content_y0 + content_h - 4 {
                    break;
                }
                let is_sel = self.selected_feature_idx == Some(f_idx);
                let card_bg = if is_sel {
                    (30, 58, 95)
                } else {
                    (19, 28, 46)
                };
                fb.fill_rect(6, fy, sidebar_w - 12, 26, card_bg);
                fb.fill_rect(6, fy, 4, 26, feat.color_rgb);
                if is_sel {
                    fb.stroke_rect(6, fy, sidebar_w - 12, 26, feat.color_rgb);
                }
                fb.draw_text(
                    14,
                    fy + 5,
                    &feat.label,
                    (241, 245, 249),
                    1,
                    Some(sidebar_w - 22),
                );
                fb.draw_text(
                    14,
                    fy + 15,
                    &format!(
                        "{} | {}..{} bp ({})",
                        feat.feature_type,
                        feat.start_col_1based,
                        feat.end_col_1based,
                        feat.strand
                    ),
                    (148, 163, 184),
                    1,
                    Some(sidebar_w - 22),
                );
                hits.push(UiHitRect {
                    x: 6,
                    y: fy,
                    w: sidebar_w - 12,
                    h: 26,
                    action_id: ACT_FEAT_BASE + (f_idx as u16),
                });
                fy += 29;
            }
        }

        // 5. Main Alignment Matrix Area (Right of Sidebar)
        let main_x0 = sidebar_w;
        let label_col_w = 215;
        let matrix_x0 = main_x0 + label_col_w;
        let matrix_w = (w - matrix_x0 - 6).max(60);

        // 5a. Feature Banner Track + Conservation Track + Coordinate Ruler (y = content_y0 .. content_y0 + 62)
        let track_h = 62;
        fb.fill_rect(main_x0, content_y0, w - main_x0, track_h, (15, 23, 42));
        fb.fill_rect(main_x0, content_y0 + track_h - 1, w - main_x0, 1, (51, 65, 85));

        fb.draw_text(
            main_x0 + 8,
            content_y0 + 8,
            "CONSERVATION / ENTROPY",
            (52, 211, 153),
            1,
            Some(label_col_w - 12),
        );
        fb.draw_text(
            main_x0 + 8,
            content_y0 + 22,
            &format!(
                "PI={} | GC={:.1}% | Ti/Tv={:.2}",
                self.workspace.parsimony_informative_sites,
                self.workspace.mean_gc_percent,
                self.workspace.ti_tv_ratio
            ),
            (148, 163, 184),
            1,
            Some(label_col_w - 12),
        );
        fb.draw_text(
            main_x0 + 8,
            content_y0 + 44,
            "TAXON / SPECIMEN",
            (56, 189, 248),
            1,
            Some(label_col_w - 12),
        );

        let col_px = self.col_width_px.max(1);
        let num_disp_cols = (matrix_w / col_px).max(1) as usize;

        for d_idx in 0..num_disp_cols {
            let col_idx = if self.mode == ViewAlphabetMode::AminoAcidCompact {
                self.scroll_col + d_idx * 3
            } else {
                self.scroll_col + d_idx
            };
            if col_idx >= self.workspace.alignment_length {
                break;
            }
            let cx = matrix_x0 + (d_idx as i32) * col_px;
            let col_1based = col_idx + 1;

            // Highlight selected column range if active
            if let Some((sel_s, sel_e)) = self.selected_col_range_1based {
                let lo = sel_s.min(sel_e);
                let hi = sel_s.max(sel_e);
                if col_1based >= lo && col_1based <= hi {
                    fb.fill_rect(cx, content_y0, col_px, track_h, (30, 58, 95));
                }
            }

            // Conservation bar
            if let Some(st) = self.workspace.column_stats.get(col_idx) {
                let bh = (st.conservation_fraction * 26.0) as i32;
                let bar_c = if st.is_parsimony_informative {
                    (168, 85, 247)
                } else if st.conservation_fraction >= 0.95 {
                    (56, 189, 248)
                } else {
                    (100, 116, 139)
                };
                fb.fill_rect(
                    cx,
                    content_y0 + 30 - bh,
                    (col_px - 1).max(1),
                    bh.max(1),
                    bar_c,
                );
            }
        }

        // Second pass for ruler ticks & coordinate numbers so background highlights never clip multi-digit labels
        for d_idx in 0..num_disp_cols {
            let col_idx = if self.mode == ViewAlphabetMode::AminoAcidCompact {
                self.scroll_col + d_idx * 3
            } else {
                self.scroll_col + d_idx
            };
            if col_idx >= self.workspace.alignment_length {
                break;
            }
            let cx = matrix_x0 + (d_idx as i32) * col_px;
            let col_1based = col_idx + 1;

            let tick_interval = if self.mode == ViewAlphabetMode::AminoAcidCompact {
                15
            } else if col_px >= 12 {
                10
            } else if col_px >= 4 {
                25
            } else {
                100
            };
            let is_major_tick = if self.mode == ViewAlphabetMode::AminoAcidCompact {
                d_idx == 0 || (d_idx % 5 == 0)
            } else {
                col_1based % tick_interval == 0 || d_idx == 0
            };
            if is_major_tick {
                fb.fill_rect(cx, content_y0 + 52, 1, 9, (226, 232, 240));
                fb.draw_text(
                    cx + 2,
                    content_y0 + 42,
                    &col_1based.to_string(),
                    (226, 232, 240),
                    1,
                    Some(54),
                );
            } else if col_px >= 10 && col_1based % 5 == 0 {
                fb.fill_rect(cx, content_y0 + 56, 1, 5, (100, 116, 139));
            }
        }

        // 5b. Sequence Rows
        let rows_y0 = content_y0 + track_h;
        let rows_avail_h = (content_h - track_h - consensus_bar_h).max(24);
        let row_h = self.row_height_px.max(12);
        let max_vis_rows = (rows_avail_h / row_h).max(1) as usize;
        let ref_dna = self
            .workspace
            .sequences
            .first()
            .map(|s| s.aligned_dna.as_slice())
            .unwrap_or(&[]);

        for r_vis in 0..max_vis_rows {
            let seq_idx = self.scroll_row + r_vis;
            let seq = match self.workspace.sequences.get(seq_idx) {
                Some(s) => s,
                None => break,
            };
            let ry = rows_y0 + (r_vis as i32) * row_h;
            let is_hover_row = self.hover_row_idx == Some(seq_idx);
            let label_bg = if is_hover_row {
                (30, 41, 59)
            } else if seq_idx % 2 == 0 {
                (15, 23, 42)
            } else {
                (11, 17, 30)
            };
            fb.fill_rect(main_x0, ry, label_col_w, row_h, label_bg);
            fb.fill_rect(main_x0, ry + row_h - 1, label_col_w, 1, (30, 41, 59));

            let label_str = format!("{:02} {}", seq_idx + 1, seq.taxon);
            fb.draw_text(
                main_x0 + 6,
                ry + (row_h - 8) / 2,
                &label_str,
                if is_hover_row {
                    (56, 189, 248)
                } else {
                    (226, 232, 240)
                },
                1,
                Some(label_col_w - 12),
            );

            // Draw cells across visible columns
            for d_idx in 0..num_disp_cols {
                let col_idx = if self.mode == ViewAlphabetMode::AminoAcidCompact {
                    self.scroll_col + d_idx * 3
                } else {
                    self.scroll_col + d_idx
                };
                if col_idx >= self.workspace.alignment_length {
                    break;
                }
                let cx = matrix_x0 + (d_idx as i32) * col_px;
                let col_stat = self.workspace.column_stats.get(col_idx);
                let raw_base = seq.aligned_dna.get(col_idx).copied().unwrap_or(b'-');
                let (codon_pos, aa_byte, _) = self.workspace.codon_and_aa_at(
                    seq_idx,
                    col_idx,
                    self.genetic_code,
                    self.reading_frame,
                );

                match self.mode {
                    ViewAlphabetMode::DualDnaAndAminoAcid => {
                        let top_h = row_h / 2;
                        let bot_h = row_h - top_h;
                        let (bg_dna, fg_dna) = residue_cell_colors(
                            raw_base,
                            false,
                            ViewAlphabetMode::Nucleotide,
                            self.color_scheme,
                            col_stat,
                        );
                        fb.fill_rect(cx, ry, (col_px - 1).max(1), top_h - 1, bg_dna);
                        if col_px >= 8 {
                            fb.draw_char(
                                cx + (col_px - 5) / 2,
                                ry + (top_h - 7) / 2,
                                raw_base as char,
                                fg_dna,
                                1,
                            );
                        }
                        let (bg_aa, fg_aa) = residue_cell_colors(
                            aa_byte,
                            true,
                            ViewAlphabetMode::AminoAcidOverlay,
                            self.color_scheme,
                            col_stat,
                        );
                        fb.fill_rect(cx, ry + top_h, (col_px - 1).max(1), bot_h - 1, bg_aa);
                        if col_px >= 8 && codon_pos == 1 {
                            fb.draw_char(
                                cx + (col_px - 5) / 2,
                                ry + top_h + (bot_h - 7) / 2,
                                aa_byte as char,
                                fg_aa,
                                1,
                            );
                        }
                    }
                    ViewAlphabetMode::AminoAcidOverlay => {
                        let (bg, fg) = residue_cell_colors(
                            aa_byte,
                            true,
                            self.mode,
                            self.color_scheme,
                            col_stat,
                        );
                        fb.fill_rect(cx, ry, (col_px - 1).max(1), row_h - 1, bg);
                        if col_px >= 8 && (codon_pos == 1 || self.workspace.detected_protein_input) {
                            fb.draw_char(
                                cx + (col_px - 5) / 2,
                                ry + (row_h - 7) / 2,
                                aa_byte as char,
                                fg,
                                1,
                            );
                        }
                    }
                    ViewAlphabetMode::AminoAcidCompact => {
                        let (bg, fg) = residue_cell_colors(
                            aa_byte,
                            true,
                            self.mode,
                            self.color_scheme,
                            col_stat,
                        );
                        fb.fill_rect(cx, ry, (col_px - 1).max(1), row_h - 1, bg);
                        if col_px >= 8 {
                            fb.draw_char(
                                cx + (col_px - 5) / 2,
                                ry + (row_h - 7) / 2,
                                aa_byte as char,
                                fg,
                                1,
                            );
                        }
                    }
                    ViewAlphabetMode::Ry3Recoded => {
                        let ry_b = recode_ry3(raw_base);
                        let (bg, fg) =
                            residue_cell_colors(ry_b, false, self.mode, self.color_scheme, col_stat);
                        fb.fill_rect(cx, ry, (col_px - 1).max(1), row_h - 1, bg);
                        if col_px >= 8 {
                            fb.draw_char(
                                cx + (col_px - 5) / 2,
                                ry + (row_h - 7) / 2,
                                ry_b as char,
                                fg,
                                1,
                            );
                        }
                    }
                    ViewAlphabetMode::Dayhoff6Recoded => {
                        let d6 = recode_dayhoff6(aa_byte);
                        let (bg, fg) =
                            residue_cell_colors(d6, true, self.mode, self.color_scheme, col_stat);
                        fb.fill_rect(cx, ry, (col_px - 1).max(1), row_h - 1, bg);
                        if col_px >= 8 && (codon_pos == 1 || self.workspace.detected_protein_input) {
                            fb.draw_char(
                                cx + (col_px - 5) / 2,
                                ry + (row_h - 7) / 2,
                                d6 as char,
                                fg,
                                1,
                            );
                        }
                    }
                    ViewAlphabetMode::Nucleotide => {
                        let disp = match self.match_mode {
                            ConsensusMatchMode::DotsForConsensus
                                if col_stat.map(|c| c.consensus_dna) == Some(raw_base)
                                    && raw_base != b'-' =>
                            {
                                b'.'
                            }
                            ConsensusMatchMode::DotsForReferenceTaxon
                                if seq_idx > 0
                                    && ref_dna.get(col_idx) == Some(&raw_base)
                                    && raw_base != b'-' =>
                            {
                                b'.'
                            }
                            _ => raw_base,
                        };
                        let (bg, fg) = residue_cell_colors(
                            raw_base,
                            false,
                            self.mode,
                            self.color_scheme,
                            col_stat,
                        );
                        fb.fill_rect(cx, ry, (col_px - 1).max(1), row_h - 1, bg);
                        if col_px >= 8 {
                            fb.draw_char(
                                cx + (col_px - 5) / 2,
                                ry + (row_h - 7) / 2,
                                disp as char,
                                fg,
                                1,
                            );
                        }
                    }
                }
            }
        }

        // 5c. Consensus Row at Bottom of Matrix
        let cons_y = h - status_h - consensus_bar_h;
        fb.fill_rect(main_x0, cons_y, w - main_x0, consensus_bar_h, (15, 23, 42));
        fb.fill_rect(main_x0, cons_y, w - main_x0, 1, (56, 189, 248));
        fb.draw_text(
            main_x0 + 8,
            cons_y + 7,
            "MAJORITY CONSENSUS (DNA)",
            (56, 189, 248),
            1,
            Some(label_col_w - 12),
        );
        fb.draw_text(
            main_x0 + 8,
            cons_y + 19,
            "TRANSLATED CONSENSUS (AA)",
            (251, 191, 36),
            1,
            Some(label_col_w - 12),
        );

        for d_idx in 0..num_disp_cols {
            let col_idx = if self.mode == ViewAlphabetMode::AminoAcidCompact {
                self.scroll_col + d_idx * 3
            } else {
                self.scroll_col + d_idx
            };
            if let Some(st) = self.workspace.column_stats.get(col_idx) {
                let cx = matrix_x0 + (d_idx as i32) * col_px;
                if col_px >= 8 {
                    fb.draw_char(
                        cx + (col_px - 5) / 2,
                        cons_y + 7,
                        st.consensus_dna as char,
                        (56, 189, 248),
                        1,
                    );
                    if self.mode == ViewAlphabetMode::AminoAcidCompact
                        || (col_idx.saturating_sub(self.reading_frame.saturating_sub(1) % 3)) % 3 == 1
                    {
                        fb.draw_char(
                            cx + (col_px - 5) / 2,
                            cons_y + 19,
                            st.consensus_aa as char,
                            (251, 191, 36),
                            1,
                        );
                    }
                }
            }
        }

        // 6. Bottom Status & Live Residue Hover Telemetry Bar
        let stat_y = h - status_h;
        fb.fill_rect(0, stat_y, w, status_h, (6, 10, 19));
        fb.fill_rect(0, stat_y, w, 1, (51, 65, 85));

        let status_line = if let (Some(h_col1), Some(h_row)) =
            (self.hover_col_1based, self.hover_row_idx)
        {
            let col0 = h_col1.saturating_sub(1);
            if let (Some(seq), Some(st)) = (
                self.workspace.sequences.get(h_row),
                self.workspace.column_stats.get(col0),
            ) {
                let base = seq.aligned_dna.get(col0).copied().unwrap_or(b'-');
                let (_, aa, codon) = self.workspace.codon_and_aa_at(
                    h_row,
                    col0,
                    self.genetic_code,
                    self.reading_frame,
                );
                let (aa_name, aa_chem, _) = amino_acid_full_name_and_chemistry(aa);
                format!(
                    "Col {}/{} | Taxon: {} | Base: {} | Codon: {}{}{} -> {} ({}, {}) | Id: {:.0}% | H={:.2}b | dN/dS={:.2}",
                    h_col1,
                    self.workspace.alignment_length,
                    seq.taxon,
                    base as char,
                    codon[0] as char,
                    codon[1] as char,
                    codon[2] as char,
                    aa as char,
                    aa_name,
                    aa_chem,
                    st.conservation_fraction * 100.0,
                    st.shannon_entropy_bits,
                    st.codon_omega_estimate
                )
            } else {
                self.status_message.clone()
            }
        } else {
            self.status_message.clone()
        };

        fb.draw_text(10, stat_y + 10, &status_line, (226, 232, 240), 1, Some(w - 20));
        hits
    }
}
/// Structured JSON report emitted when `bio-suite align-view` is run in pipeline/snapshot mode.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentViewerSnapshotReport {
    pub tool_name: String,
    pub alignment_title: String,
    pub source_path: String,
    pub alignment_method: String,
    pub taxa_count: usize,
    pub alignment_length_bp: usize,
    pub conserved_sites: usize,
    pub variable_sites: usize,
    pub singleton_sites: usize,
    pub parsimony_informative_sites: usize,
    pub mean_gc_percent: f64,
    pub mean_gap_percent: f64,
    pub total_rcfv: f64,
    pub ti_tv_ratio: f64,
    pub snapshot_start_1based: usize,
    pub snapshot_end_1based: usize,
    pub snapshot_mode: String,
    pub color_scheme: String,
    pub genetic_code: String,
    pub reading_frame: usize,
    pub consensus_dna_window: String,
    pub consensus_aa_window: String,
    pub detected_features: Vec<AlignmentFeatureAnnotation>,
    pub snapshot_svg_path: Option<String>,
    pub snapshot_png_path: Option<String>,
    pub duration_ms: u128,
}

/// CLI & Pipeline Entrypoint for `bio-suite align-view` / `runie-align-view`.
pub fn run_align_viewer_cli(args: &[String], cache: &PoliteCacheManager) -> Result<(), String> {
    let t0 = Instant::now();

    let get_arg = |flag: &str| -> Option<String> {
        args.iter()
            .position(|a| a == flag)
            .and_then(|idx| args.get(idx + 1))
            .cloned()
    };

    let explicit_cohort = get_arg("--cohort");
    let input_path = get_arg("--input")
        .or_else(|| get_arg("-i"))
        .or_else(|| get_arg("--fasta"))
        .or_else(|| {
            if explicit_cohort.is_some() {
                return None;
            }
            // Allow positional FASTA file path: `runie-align-view path/to/file.fasta`
            // (Skip arguments whose preceding token is a flag like `--snapshot-svg` or `--gff`)
            args.iter()
                .enumerate()
                .skip(1)
                .find(|(idx, a)| {
                    let prev_is_value_flag = *idx > 1
                        && matches!(
                            args[*idx - 1].as_str(),
                            "--cohort"
                                | "--mode"
                                | "--color"
                                | "--palette"
                                | "--genetic-code"
                                | "--table"
                                | "--frame"
                                | "--region"
                                | "-r"
                                | "--feature"
                                | "--window"
                                | "--strip-gaps"
                                | "--gff"
                                | "--annotations"
                                | "--snapshot-svg"
                                | "--svg"
                                | "--snapshot-png"
                                | "--png"
                                | "--out"
                                | "-o"
                                | "--snapshot-json"
                        );
                    !prev_is_value_flag
                        && !a.starts_with('-')
                        && *a != "align-view"
                        && *a != "viewer"
                        && *a != "aliview"
                        && *a != "align-snapshot"
                        && Path::new(a).exists()
                })
                .map(|(_, a)| a.clone())
        });

    let gff_path = get_arg("--gff")
        .or_else(|| get_arg("--annotations"))
        .map(PathBuf::from);
    let cohort_name = explicit_cohort.unwrap_or_else(|| "coi-barcodes".to_string());
    let mode = get_arg("--mode")
        .map(|s| ViewAlphabetMode::from_cli_arg(&s))
        .unwrap_or(ViewAlphabetMode::DualDnaAndAminoAcid);
    let color_scheme = get_arg("--color")
        .or_else(|| get_arg("--palette"))
        .map(|s| AlignmentColorScheme::from_cli_arg(&s))
        .unwrap_or(AlignmentColorScheme::AliViewDefault);
    let genetic_code = get_arg("--genetic-code")
        .or_else(|| get_arg("--table"))
        .and_then(|s| s.parse::<u8>().ok())
        .map(GeneticCodeTable::from_id)
        .unwrap_or(GeneticCodeTable::VertebrateMito2);
    let reading_frame = get_arg("--frame")
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or(1)
        .clamp(1, 3);

    let mut match_mode = ConsensusMatchMode::ShowAllCharacters;
    if args.iter().any(|a| a == "--dots" || a == "--dots-consensus") {
        match_mode = ConsensusMatchMode::DotsForConsensus;
    } else if args.iter().any(|a| a == "--dots-ref") {
        match_mode = ConsensusMatchMode::DotsForReferenceTaxon;
    }

    let mut workspace = if let Some(ref p_str) = input_path {
        AlignmentWorkspace::from_file(
            Path::new(p_str),
            gff_path.as_deref(),
            genetic_code,
            reading_frame,
            !args.iter().any(|a| a == "--no-auto-align"),
        )?
    } else {
        AlignmentWorkspace::from_builtin_cohort(&cohort_name, cache, genetic_code)?
    };

    if args.iter().any(|a| a == "--mask-adna") {
        workspace.apply_adna_deamination_mask(genetic_code, reading_frame);
    }
    if let Some(gap_thr) = get_arg("--strip-gaps").and_then(|s| s.parse::<f64>().ok()) {
        workspace.strip_gap_heavy_columns(gap_thr, genetic_code, reading_frame);
    }

    // Determine selected region from `--region start:end`, `--feature <name>`, or `--auto-focal`
    let mut selected_range: Option<(usize, usize)> = None;
    if let Some(reg_str) = get_arg("--region").or_else(|| get_arg("-r")) {
        let parts: Vec<&str> = reg_str.split(|c| c == ':' || c == '-' || c == '.').filter(|s| !s.is_empty()).collect();
        if parts.len() >= 2 {
            if let (Ok(s), Ok(e)) = (parts[0].parse::<usize>(), parts[1].parse::<usize>()) {
                let s1 = s.max(1).min(workspace.alignment_length.max(1));
                let e1 = e.max(s1).min(workspace.alignment_length.max(1));
                selected_range = Some((s1, e1));
            }
        }
    } else if let Some(feat_query) = get_arg("--feature") {
        let fq = feat_query.to_ascii_lowercase();
        if let Some(f) = workspace.features.iter().find(|f| {
            f.feature_id.to_ascii_lowercase().contains(&fq)
                || f.label.to_ascii_lowercase().contains(&fq)
                || f.feature_type.to_ascii_lowercase().contains(&fq)
        }) {
            selected_range = Some((f.start_col_1based, f.end_col_1based));
        }
    } else if args.iter().any(|a| a == "--auto-focal") {
        let win_w = get_arg("--window")
            .and_then(|s| s.parse::<usize>().ok())
            .unwrap_or(60);
        selected_range = Some(workspace.select_best_focal_window(win_w));
    }

    let snap_svg = get_arg("--snapshot-svg").or_else(|| get_arg("--svg"));
    let snap_png = get_arg("--snapshot-png").or_else(|| get_arg("--png"));
    let snap_json = get_arg("--out").or_else(|| get_arg("-o")).or_else(|| get_arg("--snapshot-json"));
    let is_snapshot_cmd = args.iter().any(|a| {
        a == "align-snapshot"
            || a == "--headless"
            || a == "--snapshot"
            || a == "--snapshot-svg"
            || a == "--snapshot-png"
            || a == "--svg"
            || a == "--png"
    });

    if is_snapshot_cmd || snap_svg.is_some() || snap_png.is_some() || snap_json.is_some() {
        let (s1, e1) = selected_range.unwrap_or_else(|| workspace.select_best_focal_window(60));
        let mut state = AlignmentViewerState::new(workspace.clone());
        state.mode = mode;
        state.color_scheme = color_scheme;
        state.match_mode = match_mode;
        state.genetic_code = genetic_code;
        state.reading_frame = reading_frame;
        state.selected_col_range_1based = Some((s1, e1));
        state.scroll_col = s1.saturating_sub(1);

        let (svg_out, png_out) = state.capture_current_viewport_snapshots(
            1440,
            snap_svg.as_deref().map(Path::new),
            snap_png.as_deref().map(Path::new),
        )?;

        let s0 = s1.saturating_sub(1);
        let e0 = e1.min(state.workspace.alignment_length);
        let cons_dna_win = state
            .workspace
            .consensus_dna_string
            .get(s0..e0)
            .unwrap_or("")
            .to_string();
        let cons_aa_win = state.workspace.extract_consensus_aa_window(
            s1,
            e1,
            genetic_code,
            reading_frame.saturating_sub(1) % 3,
        );

        let report = AlignmentViewerSnapshotReport {
            tool_name: "bio-suite align-view (Rust-Native AliView-Modeled Alignment Viewer & Snapshot Engine)".to_string(),
            alignment_title: state.workspace.title.clone(),
            source_path: state.workspace.source_path.clone(),
            alignment_method: state.workspace.alignment_method.clone(),
            taxa_count: state.workspace.sequences.len(),
            alignment_length_bp: state.workspace.alignment_length,
            conserved_sites: state.workspace.conserved_sites,
            variable_sites: state.workspace.variable_sites,
            singleton_sites: state.workspace.singleton_sites,
            parsimony_informative_sites: state.workspace.parsimony_informative_sites,
            mean_gc_percent: (state.workspace.mean_gc_percent * 100.0).round() / 100.0,
            mean_gap_percent: (state.workspace.mean_gap_percent * 100.0).round() / 100.0,
            total_rcfv: (state.workspace.total_rcfv * 10000.0).round() / 10000.0,
            ti_tv_ratio: (state.workspace.ti_tv_ratio * 1000.0).round() / 1000.0,
            snapshot_start_1based: s1,
            snapshot_end_1based: e1,
            snapshot_mode: mode.label().to_string(),
            color_scheme: color_scheme.label().to_string(),
            genetic_code: genetic_code.label().to_string(),
            reading_frame,
            consensus_dna_window: cons_dna_win,
            consensus_aa_window: cons_aa_win,
            detected_features: state.workspace.features.clone(),
            snapshot_svg_path: Some(svg_out.display().to_string()),
            snapshot_png_path: Some(png_out.display().to_string()),
            duration_ms: t0.elapsed().as_millis(),
        };

        let json_str = serde_json::to_string_pretty(&report).unwrap_or_default();
        if let Some(ref jp) = snap_json {
            if let Some(parent) = Path::new(jp).parent() {
                let _ = fs::create_dir_all(parent);
            }
            fs::write(jp, &json_str).map_err(|e| e.to_string())?;
        }
        eprintln!(
            "[bio-suite align-view] Snapshot captured ({}..{} bp of {} bp) in {} ms -> SVG: {} | PNG: {}",
            s1,
            e1,
            report.alignment_length_bp,
            report.duration_ms,
            svg_out.display(),
            png_out.display()
        );
        if snap_json.is_none() {
            println!("{}", json_str);
        }
        if !args.iter().any(|a| a == "--gui") {
            return Ok(());
        }
    }

    // Launch Interactive Native Desktop GUI
    let mut state = AlignmentViewerState::new(workspace);
    state.mode = mode;
    state.color_scheme = color_scheme;
    state.match_mode = match_mode;
    state.genetic_code = genetic_code;
    state.reading_frame = reading_frame;
    if let Some((s1, e1)) = selected_range {
        state.selected_col_range_1based = Some((s1, e1));
        state.scroll_col = s1.saturating_sub(1);
    }

    #[cfg(target_os = "windows")]
    {
        win32_gui::run_interactive_win32_viewer(state, cache.cache_dir.clone())
    }
    #[cfg(not(target_os = "windows"))]
    {
        let (svg_out, png_out) = state.capture_current_viewport_snapshots(1440, None, None)?;
        eprintln!(
            "Non-Windows environment detected; saved snapshots to {:?} and {:?}",
            svg_out, png_out
        );
        Ok(())
    }
}

#[cfg(target_os = "windows")]
mod win32_gui {
    use super::*;
    use std::ffi::OsStr;
    use std::os::windows::ffi::OsStrExt;
    use std::ptr::{null, null_mut};

    type HWND = *mut core::ffi::c_void;
    type HMENU = *mut core::ffi::c_void;
    type HINSTANCE = *mut core::ffi::c_void;
    type HDC = *mut core::ffi::c_void;
    type HICON = *mut core::ffi::c_void;
    type HCURSOR = *mut core::ffi::c_void;
    type HBRUSH = *mut core::ffi::c_void;
    type LRESULT = isize;
    type WPARAM = usize;
    type LPARAM = isize;
    type UINT = u32;
    type DWORD = u32;
    type WORD = u16;
    type BOOL = i32;

    const WS_OVERLAPPEDWINDOW: DWORD = 0x00CF_0000;
    const WS_VISIBLE: DWORD = 0x1000_0000;
    const CS_HREDRAW: UINT = 0x0002;
    const CS_VREDRAW: UINT = 0x0001;
    const IDC_ARROW: *const u16 = 32512 as *const u16;
    const CW_USEDEFAULT: i32 = 0x8000_0000u32 as i32;

    const WM_DESTROY: UINT = 0x0002;
    const WM_SIZE: UINT = 0x0005;
    const WM_PAINT: UINT = 0x000F;
    const WM_CLOSE: UINT = 0x0010;
    const WM_ERASEBKGND: UINT = 0x0014;
    const WM_KEYDOWN: UINT = 0x0100;
    const WM_COMMAND: UINT = 0x0111;
    const WM_MOUSEMOVE: UINT = 0x0200;
    const WM_LBUTTONDOWN: UINT = 0x0201;
    const WM_LBUTTONUP: UINT = 0x0202;
    const WM_MOUSEWHEEL: UINT = 0x020A;

    const VK_ESCAPE: usize = 0x1B;
    const VK_LEFT: usize = 0x25;
    const VK_UP: usize = 0x26;
    const VK_RIGHT: usize = 0x27;
    const VK_DOWN: usize = 0x28;
    const VK_PRIOR: usize = 0x21; // PageUp
    const VK_NEXT: usize = 0x22;  // PageDown
    const VK_HOME: usize = 0x24;
    const VK_END: usize = 0x23;
    const VK_F1: usize = 0x70;
    const VK_F2: usize = 0x71;
    const VK_F3: usize = 0x72;
    const VK_F4: usize = 0x73;
    const VK_F5: usize = 0x74;
    const VK_F12: usize = 0x7B;
    const VK_OEM_PLUS: usize = 0xBB;
    const VK_OEM_MINUS: usize = 0xBD;
    const VK_ADD: usize = 0x6B;
    const VK_SUBTRACT: usize = 0x6D;

    const MF_STRING: UINT = 0x0000;
    const MF_POPUP: UINT = 0x0010;
    const MF_SEPARATOR: UINT = 0x0800;
    const BI_RGB: DWORD = 0;
    const DIB_RGB_COLORS: UINT = 0;
    const SRCCOPY: DWORD = 0x00CC_0020;

    #[repr(C)]
    struct WNDCLASSEXW {
        cb_size: UINT,
        style: UINT,
        lpfn_wnd_proc: Option<unsafe extern "system" fn(HWND, UINT, WPARAM, LPARAM) -> LRESULT>,
        cb_cls_extra: i32,
        cb_wnd_extra: i32,
        h_instance: HINSTANCE,
        h_icon: HICON,
        h_cursor: HCURSOR,
        hbr_background: HBRUSH,
        lpsz_menu_name: *const u16,
        lpsz_class_name: *const u16,
        h_icon_sm: HICON,
    }

    #[repr(C)]
    struct POINT {
        x: i32,
        y: i32,
    }

    #[repr(C)]
    struct MSG {
        hwnd: HWND,
        message: UINT,
        w_param: WPARAM,
        l_param: LPARAM,
        time: DWORD,
        pt: POINT,
    }

    #[repr(C)]
    struct RECT {
        left: i32,
        top: i32,
        right: i32,
        bottom: i32,
    }

    #[repr(C)]
    struct PAINTSTRUCT {
        hdc: HDC,
        f_erase: BOOL,
        rc_paint: RECT,
        f_restore: BOOL,
        f_inc_update: BOOL,
        rgb_reserved: [u8; 32],
    }

    #[repr(C)]
    struct BITMAPINFOHEADER {
        bi_size: DWORD,
        bi_width: i32,
        bi_height: i32,
        bi_planes: WORD,
        bi_bit_count: WORD,
        bi_compression: DWORD,
        bi_size_image: DWORD,
        bi_x_pels_per_meter: i32,
        bi_y_pels_per_meter: i32,
        bi_clr_used: DWORD,
        bi_clr_important: DWORD,
    }

    #[repr(C)]
    struct BITMAPINFO {
        bmi_header: BITMAPINFOHEADER,
        bmi_colors: [u32; 1],
    }

    #[repr(C)]
    struct OPENFILENAMEW {
        l_struct_size: DWORD,
        hwnd_owner: HWND,
        h_instance: HINSTANCE,
        lpstr_filter: *const u16,
        lpstr_custom_filter: *mut u16,
        n_max_cust_filter: DWORD,
        n_filter_index: DWORD,
        lpstr_file: *mut u16,
        n_max_file: DWORD,
        lpstr_file_title: *mut u16,
        n_max_file_title: DWORD,
        lpstr_initial_dir: *const u16,
        lpstr_title: *const u16,
        flags: DWORD,
        n_file_offset: WORD,
        n_file_extension: WORD,
        lpstr_def_ext: *const u16,
        l_cust_data: LPARAM,
        lpfn_hook: *const core::ffi::c_void,
        lp_template_name: *const u16,
        pv_reserved: *mut core::ffi::c_void,
        dw_reserved: DWORD,
        flags_ex: DWORD,
    }

    #[link(name = "kernel32")]
    extern "system" {
        fn GetModuleHandleW(lp_module_name: *const u16) -> HINSTANCE;
    }

    #[link(name = "user32")]
    extern "system" {
        fn RegisterClassExW(lpwcx: *const WNDCLASSEXW) -> WORD;
        fn CreateWindowExW(
            dw_ex_style: DWORD,
            lp_class_name: *const u16,
            lp_window_name: *const u16,
            dw_style: DWORD,
            x: i32,
            y: i32,
            n_width: i32,
            n_height: i32,
            h_wnd_parent: HWND,
            h_menu: HMENU,
            h_instance: HINSTANCE,
            lp_param: *mut core::ffi::c_void,
        ) -> HWND;
        fn DefWindowProcW(hwnd: HWND, msg: UINT, wparam: WPARAM, lparam: LPARAM) -> LRESULT;
        fn GetMessageW(lp_msg: *mut MSG, hwnd: HWND, min: UINT, max: UINT) -> BOOL;
        fn TranslateMessage(lp_msg: *const MSG) -> BOOL;
        fn DispatchMessageW(lp_msg: *const MSG) -> LRESULT;
        fn PostQuitMessage(n_exit_code: i32);
        fn DestroyWindow(hwnd: HWND) -> BOOL;
        fn LoadCursorW(h_instance: HINSTANCE, lp_cursor_name: *const u16) -> HCURSOR;
        fn BeginPaint(hwnd: HWND, lp_paint: *mut PAINTSTRUCT) -> HDC;
        fn EndPaint(hwnd: HWND, lp_paint: *const PAINTSTRUCT) -> BOOL;
        fn InvalidateRect(hwnd: HWND, lp_rect: *const RECT, b_erase: BOOL) -> BOOL;
        fn GetClientRect(hwnd: HWND, lp_rect: *mut RECT) -> BOOL;
        fn CreateMenu() -> HMENU;
        fn CreatePopupMenu() -> HMENU;
        fn AppendMenuW(h_menu: HMENU, u_flags: UINT, u_id_new_item: usize, lp_new_item: *const u16) -> BOOL;
        fn SetMenu(hwnd: HWND, h_menu: HMENU) -> BOOL;
        fn SetWindowTextW(hwnd: HWND, lp_string: *const u16) -> BOOL;
    }

    #[link(name = "gdi32")]
    extern "system" {
        fn StretchDIBits(
            hdc: HDC,
            x_dest: i32,
            y_dest: i32,
            dest_width: i32,
            dest_height: i32,
            x_src: i32,
            y_src: i32,
            src_width: i32,
            src_height: i32,
            lp_bits: *const core::ffi::c_void,
            lpbmi: *const BITMAPINFO,
            i_usage: UINT,
            rop: DWORD,
        ) -> i32;
    }

    #[link(name = "comdlg32")]
    extern "system" {
        fn GetOpenFileNameW(lpofn: *mut OPENFILENAMEW) -> BOOL;
    }

    fn to_wide(s: &str) -> Vec<u16> {
        OsStr::new(s).encode_wide().chain(std::iter::once(0)).collect()
    }

    struct Win32ViewerContext {
        state: AlignmentViewerState,
        fb: SoftwareFrameBuffer,
        hits: Vec<UiHitRect>,
        cache_dir: PathBuf,
    }

    static mut GLOBAL_VIEWER_CTX: *mut Win32ViewerContext = null_mut();

    unsafe fn build_top_menubar(features: &[AlignmentFeatureAnnotation]) -> HMENU {
        let menubar = CreateMenu();

        // 1. File Menu
        let file_menu = CreatePopupMenu();
        AppendMenuW(file_menu, MF_STRING, ACT_OPEN_FILE as usize, to_wide("&Open FASTA / Alignment...\tCtrl+O").as_ptr());
        AppendMenuW(file_menu, MF_SEPARATOR, 0, null());
        AppendMenuW(file_menu, MF_STRING, ACT_LOAD_COI as usize, to_wide("Load Cohort: &Carnivora COX1 Barcodes (16 Taxa)").as_ptr());
        AppendMenuW(file_menu, MF_STRING, ACT_LOAD_CHR1 as usize, to_wide("Load Cohort: &Human GRCh38 Chr1p36.33 Panel").as_ptr());
        AppendMenuW(file_menu, MF_STRING, ACT_LOAD_HYMEN as usize, to_wide("Load Cohort: H&ymenoptera Eusociality Clade").as_ptr());
        AppendMenuW(file_menu, MF_SEPARATOR, 0, null());
        AppendMenuW(file_menu, MF_STRING, ACT_SAVE_FASTA as usize, to_wide("&Save Current Alignment (FASTA)...\tCtrl+S").as_ptr());
        AppendMenuW(file_menu, MF_STRING, ACT_SNAPSHOT_SVG_PNG as usize, to_wide("Capture Paper &Snapshot (SVG + PNG)\tF12").as_ptr());
        AppendMenuW(file_menu, MF_SEPARATOR, 0, null());
        AppendMenuW(file_menu, MF_STRING, ACT_EXIT as usize, to_wide("E&xit\tAlt+F4").as_ptr());
        AppendMenuW(menubar, MF_POPUP, file_menu as usize, to_wide("&File").as_ptr());

        // 2. View Menu
        let view_menu = CreatePopupMenu();
        AppendMenuW(view_menu, MF_STRING, ACT_MODE_DNA as usize, to_wide("&1. DNA / Nucleotide View (4-Base)\tKey 1").as_ptr());
        AppendMenuW(view_menu, MF_STRING, ACT_MODE_AA_OVERLAY as usize, to_wide("&2. Amino Acid Translation (3-bp Codon Span)\tKey 2").as_ptr());
        AppendMenuW(view_menu, MF_STRING, ACT_MODE_AA_COMPACT as usize, to_wide("&3. Amino Acid Compact (1 Column / Codon)\tKey 3").as_ptr());
        AppendMenuW(view_menu, MF_STRING, ACT_MODE_DUAL as usize, to_wide("&4. Dual DNA + Amino Acid Split Track\tKey 4").as_ptr());
        AppendMenuW(view_menu, MF_STRING, ACT_MODE_RY3 as usize, to_wide("&5. RuNIE RY3 Purine/Pyrimidine Recoding\tKey 5").as_ptr());
        AppendMenuW(view_menu, MF_STRING, ACT_MODE_DAYHOFF6 as usize, to_wide("&6. RuNIE Dayhoff-6 Functional Group Recoding\tKey 6").as_ptr());
        AppendMenuW(view_menu, MF_SEPARATOR, 0, null());
        AppendMenuW(view_menu, MF_STRING, ACT_MATCH_ALL as usize, to_wide("Show &All Residue Characters\tKey A").as_ptr());
        AppendMenuW(view_menu, MF_STRING, ACT_MATCH_DOTS_CONS as usize, to_wide("Show &Dots (.) for Consensus Matches\tKey D").as_ptr());
        AppendMenuW(view_menu, MF_STRING, ACT_MATCH_DOTS_REF as usize, to_wide("Show Dots (.) for &Reference Taxon #1 Matches\tKey R").as_ptr());
        AppendMenuW(view_menu, MF_SEPARATOR, 0, null());
        AppendMenuW(view_menu, MF_STRING, ACT_ZOOM_IN as usize, to_wide("Zoom &In (+)\t+").as_ptr());
        AppendMenuW(view_menu, MF_STRING, ACT_ZOOM_OUT as usize, to_wide("Zoom &Out (-)\t-").as_ptr());
        AppendMenuW(view_menu, MF_STRING, ACT_ZOOM_OVERVIEW as usize, to_wide("&Bird's-Eye Whole-Alignment Overview (2 px/col)\tKey 0").as_ptr());
        AppendMenuW(view_menu, MF_STRING, ACT_ZOOM_STANDARD as usize, to_wide("&Standard Residue Reading View (14 px/col)").as_ptr());
        AppendMenuW(view_menu, MF_STRING, ACT_ZOOM_FOCAL_BEST as usize, to_wide("Jump to &Best Focal Diagnostic Window (60 bp)\tKey F").as_ptr());
        AppendMenuW(menubar, MF_POPUP, view_menu as usize, to_wide("&View").as_ptr());

        // 3. Color Palette Menu
        let color_menu = CreatePopupMenu();
        AppendMenuW(color_menu, MF_STRING, ACT_COLOR_ALIVIEW as usize, to_wide("&AliView Default Palette (DNA + Taylor AA)").as_ptr());
        AppendMenuW(color_menu, MF_STRING, ACT_COLOR_CLUSTAL as usize, to_wide("&ClustalX Amino Acid Chemistry").as_ptr());
        AppendMenuW(color_menu, MF_STRING, ACT_COLOR_HYDRO as usize, to_wide("&Kyte-Doolittle Hydrophobicity (Warm-to-Cool)").as_ptr());
        AppendMenuW(color_menu, MF_STRING, ACT_COLOR_RUNIE_PI as usize, to_wide("&RuNIE Parsimony-Informative & Conservation Highlight").as_ptr());
        AppendMenuW(color_menu, MF_STRING, ACT_COLOR_OKABE as usize, to_wide("&Okabe-Ito Color-Blind Safe Palette").as_ptr());
        AppendMenuW(menubar, MF_POPUP, color_menu as usize, to_wide("&Color Palette").as_ptr());

        // 4. Translation & Genetic Code Menu
        let trans_menu = CreatePopupMenu();
        AppendMenuW(trans_menu, MF_STRING, ACT_FRAME_1 as usize, to_wide("Reading Frame &+1\tF1").as_ptr());
        AppendMenuW(trans_menu, MF_STRING, ACT_FRAME_2 as usize, to_wide("Reading Frame +&2\tF2").as_ptr());
        AppendMenuW(trans_menu, MF_STRING, ACT_FRAME_3 as usize, to_wide("Reading Frame +&3\tF3").as_ptr());
        AppendMenuW(trans_menu, MF_STRING, ACT_REV_COMP as usize, to_wide("&Reverse Complement Alignment").as_ptr());
        AppendMenuW(trans_menu, MF_SEPARATOR, 0, null());
        AppendMenuW(trans_menu, MF_STRING, ACT_GCODE_1 as usize, to_wide("NCBI Table &1: Standard Nuclear Code").as_ptr());
        AppendMenuW(trans_menu, MF_STRING, ACT_GCODE_2 as usize, to_wide("NCBI Table &2: Vertebrate Mitochondrial Code").as_ptr());
        AppendMenuW(trans_menu, MF_STRING, ACT_GCODE_5 as usize, to_wide("NCBI Table &5: Invertebrate Mitochondrial Code").as_ptr());
        AppendMenuW(trans_menu, MF_STRING, ACT_GCODE_11 as usize, to_wide("NCBI Table 1&1: Bacterial / Archaeal / Plastid Code").as_ptr());
        AppendMenuW(menubar, MF_POPUP, trans_menu as usize, to_wide("&Translation").as_ptr());

        // 5. Features & Annotations Menu
        let feat_menu = CreatePopupMenu();
        AppendMenuW(feat_menu, MF_STRING, ACT_TOGGLE_SIDEBAR as usize, to_wide("Toggle Feature &Annotation Sidebar\tF4").as_ptr());
        AppendMenuW(feat_menu, MF_STRING, ACT_RESCAN_FEATURES as usize, to_wide("Re-Scan Alignment for &Motifs & Domains\tF5").as_ptr());
        AppendMenuW(feat_menu, MF_SEPARATOR, 0, null());
        for (idx, f) in features.iter().take(15).enumerate() {
            let item_txt = format!("Zoom to [{}] {}", f.feature_id, f.label);
            AppendMenuW(feat_menu, MF_STRING, (ACT_FEAT_BASE + idx as u16) as usize, to_wide(&item_txt).as_ptr());
        }
        AppendMenuW(menubar, MF_POPUP, feat_menu as usize, to_wide("F&eatures").as_ptr());

        // 6. RuNIE Pipeline Menu
        let runie_menu = CreatePopupMenu();
        AppendMenuW(runie_menu, MF_STRING, ACT_MASK_ADNA as usize, to_wide("Mask Terminal 5' C->T / 3' G->A &aDNA Deamination").as_ptr());
        AppendMenuW(runie_menu, MF_STRING, ACT_STRIP_GAPS as usize, to_wide("Strip &Gap-Heavy Alignment Columns (>50% Gaps)").as_ptr());
        AppendMenuW(runie_menu, MF_SEPARATOR, 0, null());
        AppendMenuW(runie_menu, MF_STRING, ACT_SNAPSHOT_SVG_PNG as usize, to_wide("Export Publication &Figure 2 Snapshot (SVG + PNG)\tF12").as_ptr());
        AppendMenuW(menubar, MF_POPUP, runie_menu as usize, to_wide("&RuNIE Pipeline").as_ptr());

        menubar
    }

    unsafe fn prompt_open_fasta_dialog(hwnd: HWND) -> Option<PathBuf> {
        let mut file_buf = vec![0u16; 520];
        let filter = "Alignment & Sequence Files (*.fasta;*.fa;*.afa;*.aln;*.json)\0*.fasta;*.fa;*.afa;*.aln;*.json\0All Files (*.*)\0*.*\0\0"
            .encode_utf16()
            .collect::<Vec<u16>>();
        let title = to_wide("Open Sequence Alignment (FASTA / CLUSTAL / bio-suite JSON)");
        let mut ofn = OPENFILENAMEW {
            l_struct_size: core::mem::size_of::<OPENFILENAMEW>() as DWORD,
            hwnd_owner: hwnd,
            h_instance: null_mut(),
            lpstr_filter: filter.as_ptr(),
            lpstr_custom_filter: null_mut(),
            n_max_cust_filter: 0,
            n_filter_index: 1,
            lpstr_file: file_buf.as_mut_ptr(),
            n_max_file: file_buf.len() as DWORD,
            lpstr_file_title: null_mut(),
            n_max_file_title: 0,
            lpstr_initial_dir: null(),
            lpstr_title: title.as_ptr(),
            flags: 0x0008_0000 | 0x0000_1000 | 0x0000_0008, // OFN_EXPLORER | OFN_FILEMUSTEXIST | OFN_NOCHANGEDIR
            n_file_offset: 0,
            n_file_extension: 0,
            lpstr_def_ext: null(),
            l_cust_data: 0,
            lpfn_hook: null(),
            lp_template_name: null(),
            pv_reserved: null_mut(),
            dw_reserved: 0,
            flags_ex: 0,
        };
        if GetOpenFileNameW(&mut ofn) != 0 {
            let len = file_buf.iter().position(|&c| c == 0).unwrap_or(0);
            if len > 0 {
                let s = String::from_utf16_lossy(&file_buf[..len]);
                return Some(PathBuf::from(s));
            }
        }
        None
    }

    unsafe extern "system" fn wnd_proc(
        hwnd: HWND,
        msg: UINT,
        wparam: WPARAM,
        lparam: LPARAM,
    ) -> LRESULT {
        let ctx_ptr = GLOBAL_VIEWER_CTX;
        if ctx_ptr.is_null() {
            return DefWindowProcW(hwnd, msg, wparam, lparam);
        }
        let ctx = &mut *ctx_ptr;

        match msg {
            WM_ERASEBKGND => 1,
            WM_SIZE => {
                let mut rc = RECT { left: 0, top: 0, right: 1280, bottom: 800 };
                GetClientRect(hwnd, &mut rc);
                let w = (rc.right - rc.left).max(320) as usize;
                let h = (rc.bottom - rc.top).max(240) as usize;
                ctx.fb.resize(w, h);
                ctx.state.clamp_scroll(w as i32, h as i32);
                ctx.hits = ctx.state.render_interactive_frame(&mut ctx.fb);
                InvalidateRect(hwnd, null(), 0);
                0
            }
            WM_COMMAND => {
                let cmd_id = (wparam & 0xFFFF) as u16;
                match cmd_id {
                    ACT_EXIT => {
                        DestroyWindow(hwnd);
                    }
                    ACT_OPEN_FILE => {
                        if let Some(path) = prompt_open_fasta_dialog(hwnd) {
                            match AlignmentWorkspace::from_file(
                                &path,
                                None,
                                ctx.state.genetic_code,
                                ctx.state.reading_frame,
                                true,
                            ) {
                                Ok(ws) => {
                                    ctx.state = AlignmentViewerState::new(ws);
                                    let win_title = to_wide(&format!(
                                        "CRV Science Labs — bio-suite AlignView | {} ({} Taxa, {} bp)",
                                        ctx.state.workspace.title,
                                        ctx.state.workspace.sequences.len(),
                                        ctx.state.workspace.alignment_length
                                    ));
                                    SetWindowTextW(hwnd, win_title.as_ptr());
                                }
                                Err(e) => {
                                    ctx.state.status_message = format!("Open error: {}", e);
                                }
                            }
                        }
                    }
                    ACT_LOAD_COI | ACT_LOAD_CHR1 | ACT_LOAD_HYMEN => {
                        let cohort = match cmd_id {
                            ACT_LOAD_CHR1 => "chr1-section-a",
                            ACT_LOAD_HYMEN => "hymenoptera",
                            _ => "coi-barcodes",
                        };
                        let cache = PoliteCacheManager::new(&ctx.cache_dir);
                        if let Ok(ws) = AlignmentWorkspace::from_builtin_cohort(
                            cohort,
                            &cache,
                            ctx.state.genetic_code,
                        ) {
                            ctx.state = AlignmentViewerState::new(ws);
                        }
                    }
                    ACT_SAVE_FASTA => {
                        let out_dir = PathBuf::from("runs/snapshots");
                        let _ = fs::create_dir_all(&out_dir);
                        let out_p = out_dir.join("exported_alignment.fasta");
                        let fasta_txt = ctx.state.workspace.export_fasta(
                            ctx.state.selected_col_range_1based,
                            false,
                            ctx.state.genetic_code,
                            ctx.state.reading_frame,
                        );
                        let _ = fs::write(&out_p, fasta_txt);
                        ctx.state.status_message = format!("Saved FASTA alignment to {}", out_p.display());
                    }
                    _ => {
                        let w = ctx.fb.width as i32;
                        ctx.state.execute_action(cmd_id, w);
                    }
                }
                ctx.hits = ctx.state.render_interactive_frame(&mut ctx.fb);
                InvalidateRect(hwnd, null(), 0);
                0
            }
            WM_LBUTTONDOWN => {
                let mx = (lparam & 0xFFFF) as i16 as i32;
                let my = ((lparam >> 16) & 0xFFFF) as i16 as i32;
                let w = ctx.fb.width as i32;

                // 1. Check toolbar & feature sidebar buttons
                let mut matched_act = None;
                for h in &ctx.hits {
                    if mx >= h.x && mx <= h.x + h.w && my >= h.y && my <= h.y + h.h {
                        matched_act = Some(h.action_id);
                        break;
                    }
                }
                if let Some(act) = matched_act {
                    ctx.state.execute_action(act, w);
                } else if my >= 42 && my <= 78 && mx >= 150 {
                    // 2. Clicked Global Minimap -> jump to coordinate
                    ctx.state.dragging_minimap = true;
                    let map_w = (w - 164).max(100);
                    let frac = ((mx - 150) as f64 / map_w as f64).clamp(0.0, 1.0);
                    ctx.state.scroll_col = ((frac * ctx.state.workspace.alignment_length as f64) as usize)
                        .saturating_sub(ctx.state.visible_cols_count(w) / 2);
                    ctx.state.clamp_scroll(w, ctx.fb.height as i32);
                } else {
                    // 3. Clicked inside alignment matrix -> start column selection
                    let sidebar_w = if ctx.state.show_feature_sidebar { 245 } else { 0 };
                    let matrix_x0 = sidebar_w + 215;
                    if mx >= matrix_x0 && my >= 78 {
                        let d_col = ((mx - matrix_x0) / ctx.state.col_width_px.max(1)).max(0) as usize;
                        let col_1based = if ctx.state.mode == ViewAlphabetMode::AminoAcidCompact {
                            ctx.state.scroll_col + d_col * 3 + 1
                        } else {
                            ctx.state.scroll_col + d_col + 1
                        }
                        .min(ctx.state.workspace.alignment_length.max(1));
                        ctx.state.dragging_selection_start_col = Some(col_1based);
                        ctx.state.selected_col_range_1based = Some((col_1based, col_1based));
                    }
                }

                ctx.hits = ctx.state.render_interactive_frame(&mut ctx.fb);
                InvalidateRect(hwnd, null(), 0);
                0
            }
            WM_LBUTTONUP => {
                ctx.state.dragging_minimap = false;
                ctx.state.dragging_selection_start_col = None;
                0
            }
            WM_MOUSEMOVE => {
                let mx = (lparam & 0xFFFF) as i16 as i32;
                let my = ((lparam >> 16) & 0xFFFF) as i16 as i32;
                let w = ctx.fb.width as i32;
                let h = ctx.fb.height as i32;

                if ctx.state.dragging_minimap {
                    let map_w = (w - 164).max(100);
                    let frac = ((mx - 150) as f64 / map_w as f64).clamp(0.0, 1.0);
                    ctx.state.scroll_col = ((frac * ctx.state.workspace.alignment_length as f64) as usize)
                        .saturating_sub(ctx.state.visible_cols_count(w) / 2);
                    ctx.state.clamp_scroll(w, h);
                }

                let sidebar_w = if ctx.state.show_feature_sidebar { 245 } else { 0 };
                let matrix_x0 = sidebar_w + 215;
                let rows_y0 = 78 + 62;
                if mx >= matrix_x0 && my >= rows_y0 && my < h - 62 {
                    let d_col = ((mx - matrix_x0) / ctx.state.col_width_px.max(1)).max(0) as usize;
                    let col_1based = if ctx.state.mode == ViewAlphabetMode::AminoAcidCompact {
                        ctx.state.scroll_col + d_col * 3 + 1
                    } else {
                        ctx.state.scroll_col + d_col + 1
                    };
                    let row_idx =
                        ctx.state.scroll_row + ((my - rows_y0) / ctx.state.row_height_px.max(12)).max(0) as usize;
                    if col_1based <= ctx.state.workspace.alignment_length
                        && row_idx < ctx.state.workspace.sequences.len()
                    {
                        ctx.state.hover_col_1based = Some(col_1based);
                        ctx.state.hover_row_idx = Some(row_idx);
                        if let Some(start_c) = ctx.state.dragging_selection_start_col {
                            ctx.state.selected_col_range_1based = Some((start_c, col_1based));
                        }
                    }
                } else {
                    ctx.state.hover_col_1based = None;
                    ctx.state.hover_row_idx = None;
                }

                ctx.hits = ctx.state.render_interactive_frame(&mut ctx.fb);
                InvalidateRect(hwnd, null(), 0);
                0
            }
            WM_MOUSEWHEEL => {
                let delta = ((wparam >> 16) & 0xFFFF) as i16;
                let key_flags = (wparam & 0xFFFF) as usize;
                let w = ctx.fb.width as i32;
                let h = ctx.fb.height as i32;
                // If Ctrl (0x0008) is held, zoom in/out; otherwise scroll horizontally
                if (key_flags & 0x0008) != 0 {
                    if delta > 0 {
                        ctx.state.execute_action(ACT_ZOOM_IN, w);
                    } else {
                        ctx.state.execute_action(ACT_ZOOM_OUT, w);
                    }
                } else {
                    let step = (ctx.state.visible_cols_count(w) / 8).max(3);
                    if delta > 0 {
                        ctx.state.scroll_col = ctx.state.scroll_col.saturating_sub(step);
                    } else {
                        ctx.state.scroll_col += step;
                    }
                    ctx.state.clamp_scroll(w, h);
                }
                ctx.hits = ctx.state.render_interactive_frame(&mut ctx.fb);
                InvalidateRect(hwnd, null(), 0);
                0
            }
            WM_KEYDOWN => {
                let w = ctx.fb.width as i32;
                let h = ctx.fb.height as i32;
                let step = (ctx.state.visible_cols_count(w) / 6).max(3);
                match wparam {
                    VK_LEFT => ctx.state.scroll_col = ctx.state.scroll_col.saturating_sub(step),
                    VK_RIGHT => ctx.state.scroll_col += step,
                    VK_UP => ctx.state.scroll_row = ctx.state.scroll_row.saturating_sub(1),
                    VK_DOWN => ctx.state.scroll_row += 1,
                    VK_PRIOR => ctx.state.scroll_col = ctx.state.scroll_col.saturating_sub(step * 5),
                    VK_NEXT => ctx.state.scroll_col += step * 5,
                    VK_HOME => ctx.state.scroll_col = 0,
                    VK_END => ctx.state.scroll_col = ctx.state.workspace.alignment_length,
                    VK_ESCAPE => {
                        ctx.state.selected_col_range_1based = None;
                        ctx.state.status_message = "Cleared region selection".to_string();
                    }
                    0x31 => ctx.state.execute_action(ACT_MODE_DNA, w),
                    0x32 => ctx.state.execute_action(ACT_MODE_AA_OVERLAY, w),
                    0x33 => ctx.state.execute_action(ACT_MODE_AA_COMPACT, w),
                    0x34 => ctx.state.execute_action(ACT_MODE_DUAL, w),
                    0x35 => ctx.state.execute_action(ACT_MODE_RY3, w),
                    0x36 => ctx.state.execute_action(ACT_MODE_DAYHOFF6, w),
                    0x30 => ctx.state.execute_action(ACT_ZOOM_OVERVIEW, w),
                    0x46 => ctx.state.execute_action(ACT_ZOOM_FOCAL_BEST, w),
                    0x44 => ctx.state.execute_action(ACT_MATCH_DOTS_CONS, w),
                    0x52 => ctx.state.execute_action(ACT_MATCH_DOTS_REF, w),
                    0x41 => ctx.state.execute_action(ACT_MATCH_ALL, w),
                    VK_OEM_PLUS | VK_ADD => ctx.state.execute_action(ACT_ZOOM_IN, w),
                    VK_OEM_MINUS | VK_SUBTRACT => ctx.state.execute_action(ACT_ZOOM_OUT, w),
                    VK_F1 => ctx.state.execute_action(ACT_FRAME_1, w),
                    VK_F2 => ctx.state.execute_action(ACT_FRAME_2, w),
                    VK_F3 => ctx.state.execute_action(ACT_FRAME_3, w),
                    VK_F4 => ctx.state.execute_action(ACT_TOGGLE_SIDEBAR, w),
                    VK_F5 => ctx.state.execute_action(ACT_RESCAN_FEATURES, w),
                    VK_F12 => ctx.state.execute_action(ACT_SNAPSHOT_SVG_PNG, w),
                    _ => {}
                }
                ctx.state.clamp_scroll(w, h);
                ctx.hits = ctx.state.render_interactive_frame(&mut ctx.fb);
                InvalidateRect(hwnd, null(), 0);
                0
            }
            WM_PAINT => {
                let mut ps: PAINTSTRUCT = core::mem::zeroed();
                let hdc = BeginPaint(hwnd, &mut ps);
                let w = ctx.fb.width as i32;
                let h = ctx.fb.height as i32;
                let bmi = BITMAPINFO {
                    bmi_header: BITMAPINFOHEADER {
                        bi_size: core::mem::size_of::<BITMAPINFOHEADER>() as DWORD,
                        bi_width: w,
                        bi_height: -h, // Top-down DIB
                        bi_planes: 1,
                        bi_bit_count: 32,
                        bi_compression: BI_RGB,
                        bi_size_image: 0,
                        bi_x_pels_per_meter: 0,
                        bi_y_pels_per_meter: 0,
                        bi_clr_used: 0,
                        bi_clr_important: 0,
                    },
                    bmi_colors: [0],
                };
                StretchDIBits(
                    hdc,
                    0,
                    0,
                    w,
                    h,
                    0,
                    0,
                    w,
                    h,
                    ctx.fb.pixels.as_ptr() as *const core::ffi::c_void,
                    &bmi,
                    DIB_RGB_COLORS,
                    SRCCOPY,
                );
                EndPaint(hwnd, &ps);
                0
            }
            WM_CLOSE => {
                DestroyWindow(hwnd);
                0
            }
            WM_DESTROY => {
                PostQuitMessage(0);
                0
            }
            _ => DefWindowProcW(hwnd, msg, wparam, lparam),
        }
    }

    pub fn run_interactive_win32_viewer(
        state: AlignmentViewerState,
        cache_dir: PathBuf,
    ) -> Result<(), String> {
        unsafe {
            let h_instance = GetModuleHandleW(null());
            let class_name = to_wide("CrvBioSuiteAlignViewWndClass");

            let wc = WNDCLASSEXW {
                cb_size: core::mem::size_of::<WNDCLASSEXW>() as UINT,
                style: CS_HREDRAW | CS_VREDRAW,
                lpfn_wnd_proc: Some(wnd_proc),
                cb_cls_extra: 0,
                cb_wnd_extra: 0,
                h_instance,
                h_icon: null_mut(),
                h_cursor: LoadCursorW(null_mut(), IDC_ARROW),
                hbr_background: null_mut(),
                lpsz_menu_name: null(),
                lpsz_class_name: class_name.as_ptr(),
                h_icon_sm: null_mut(),
            };
            RegisterClassExW(&wc);

            let init_w = 1420usize;
            let init_h = 840usize;
            let mut fb = SoftwareFrameBuffer::new(init_w, init_h);
            let hits = state.render_interactive_frame(&mut fb);
            let menubar = build_top_menubar(&state.workspace.features);

            let win_title = to_wide(&format!(
                "CRV Science Labs — bio-suite AlignView (AliView-Modeled Rust Alignment Inspector) | {} ({} Taxa, {} bp)",
                state.workspace.title,
                state.workspace.sequences.len(),
                state.workspace.alignment_length
            ));

            let mut ctx = Box::new(Win32ViewerContext {
                state,
                fb,
                hits,
                cache_dir,
            });
            GLOBAL_VIEWER_CTX = &mut *ctx as *mut Win32ViewerContext;

            let hwnd = CreateWindowExW(
                0,
                class_name.as_ptr(),
                win_title.as_ptr(),
                WS_OVERLAPPEDWINDOW | WS_VISIBLE,
                CW_USEDEFAULT,
                CW_USEDEFAULT,
                init_w as i32,
                init_h as i32,
                null_mut(),
                menubar,
                h_instance,
                null_mut(),
            );
            if hwnd.is_null() {
                GLOBAL_VIEWER_CTX = null_mut();
                return Err("Failed to create Win32 window for AlignView".to_string());
            }

            SetMenu(hwnd, menubar);

            let mut msg: MSG = core::mem::zeroed();
            while GetMessageW(&mut msg, null_mut(), 0, 0) > 0 {
                TranslateMessage(&msg);
                DispatchMessageW(&msg);
            }

            GLOBAL_VIEWER_CTX = null_mut();
            Ok(())
        }
    }
}
