//! Tier 2 (Engines 2 & 3): Structural/Functional Genomics (`bio-suite genome`)
//! & Genetic Target Identification / Drug Discovery (`bio-suite target`)
//! plus the GPL-3.0 Verified Tool Registry.

use crate::complexity::{round4, SlidingWindowRecord};
use crate::fetch::GenomicAnnotation;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotatedComplexityRegion {
    pub region_id: String,
    pub window_id: String,
    pub start: usize,
    pub end: usize,
    pub structural_complexity_z: f64,
    pub gc_percent: f64,
    pub cpg_obs_exp: f64,
    pub g4hunter_score: f64,
    pub zdna_score: f64,
    pub lz76_normalized: f64,
    pub dfa_alpha: f64,
    pub overlapping_cds_count: usize,
    pub overlapping_repeat_count: usize,
    pub overlapping_features: Vec<String>,
    pub narrative_summary: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FineMappedVariantCredibleMember {
    pub rsid: String,
    pub position: usize,
    pub beta: f64,
    pub se: f64,
    pub z_score: f64,
    pub wakefield_abf: f64,
    pub posterior_inclusion_prob: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ToolAdapterEntry {
    pub engine: &'static str,
    pub category: &'static str,
    pub tool_name: &'static str,
    pub binary_cmd: &'static str,
    pub spdx_license: &'static str,
    pub gpl3_compatible: bool,
    pub native_rust_fallback: &'static str,
    pub description: &'static str,
}

/// Complete GPL-3.0-audited tool registry across all 4 `bio-suite` engines.
pub const BIO_SUITE_TOOL_REGISTRY: &[ToolAdapterEntry] = &[
    // Engine 1: Molecular Systematics & Phylogenomics (`bio-suite phylo`)
    ToolAdapterEntry {
        engine: "phylo",
        category: "Orthology & Homology Inference",
        tool_name: "OrthoFinder / DIAMOND / BUSCO",
        binary_cmd: "orthofinder",
        spdx_license: "GPL-3.0-or-later / MIT",
        gpl3_compatible: true,
        native_rust_fallback: "Reciprocal Best-Hit (RBH) k-mer homology & NCBI Ortholog Fetcher",
        description: "Identifies single-copy orthologs, paralogs, and conserved benchmark loci across proteomes and coding sequences.",
    },
    ToolAdapterEntry {
        engine: "phylo",
        category: "Multiple Sequence Alignment & Trimming",
        tool_name: "MAFFT / MUSCLE5 / trimAl",
        binary_cmd: "mafft",
        spdx_license: "BSD-3-Clause / GPL-3.0",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite parallel codon-anchored progressive MSA & entropy block trimmer",
        description: "Aligns coding or non-coding loci while preserving triplet reading frames and trimming saturated sites.",
    },
    ToolAdapterEntry {
        engine: "phylo",
        category: "ML & Bayesian Tree Inference",
        tool_name: "IQ-TREE 2 / RAxML-ng / RevBayes / MrBayes",
        binary_cmd: "iqtree2",
        spdx_license: "GPL-2.0-or-later / GPL-3.0",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite K2P / TN93 Neighbor-Joining & UPGMA with Rayon parallel bootstrap",
        description: "Infers phylogenetic topologies, substitution model parameters, and clade bootstrap/posterior support.",
    },
    ToolAdapterEntry {
        engine: "phylo",
        category: "Multispecies Coalescent & Introgression",
        tool_name: "ASTRAL-IV / Dsuite / SNaQ",
        binary_cmd: "astral4",
        spdx_license: "Apache-2.0 / GPL-3.0 / MIT",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite ABBA-BABA Patterson's D & ABGD Barcode Gap MOTU Delimiter",
        description: "Resolves incomplete lineage sorting (ILS), hybridization, and DNA barcoding species boundaries.",
    },
    ToolAdapterEntry {
        engine: "phylo",
        category: "Molecular Clock & Selection (dN/dS)",
        tool_name: "HyPhy (MEME/BUSTED/aBSREL) / TreePL / BEAST 2",
        binary_cmd: "hyphy",
        spdx_license: "MIT / GPL-3.0 / LGPL-2.1",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite Nei-Gojobori dN/dS & synapomorphic codon state reconstructor",
        description: "Detects episodic diversifying selection and reconstructs ancestral nucleotide/amino-acid states.",
    },
    // Engine 2: Structural, Population & Functional Genomics (`bio-suite genome`)
    ToolAdapterEntry {
        engine: "genome",
        category: "Long-Read Assembly & Phasing",
        tool_name: "hifiasm / Verkko / Flye",
        binary_cmd: "hifiasm",
        spdx_license: "MIT / BSD-3-Clause",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite overlap-layout contig statistics & N50/L50 assembly evaluator",
        description: "Haplotype-resolved telomere-to-telomere (T2T) assembly from PacBio HiFi, ONT Ultra-Long, and Hi-C reads.",
    },
    ToolAdapterEntry {
        engine: "genome",
        category: "Pangenome Graphs & Structural Variants",
        tool_name: "Minigraph-Cactus / PGGB / vg / Sniffles2",
        binary_cmd: "vg",
        spdx_license: "MIT",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite indel/inversion breakpoint & variant spectrum analyzer",
        description: "Constructs reference-bias-free variation graphs and calls SNPs, indels, and structural variants.",
    },
    ToolAdapterEntry {
        engine: "genome",
        category: "Ab Initio Annotation & Repeat Masking (GPL-3 Drop-In)",
        tool_name: "Helixer + AUGUSTUS + EDTA + Liftoff",
        binary_cmd: "helixer",
        spdx_license: "GPL-3.0-or-later",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite ORF/CDS extractor & tandem-repeat / GFF3 coordinate intersector",
        description: "100% GPL-3.0 replacement for restricted BRAKER3/RepeatMasker stacks for TE annotation and gene calling.",
    },
    ToolAdapterEntry {
        engine: "genome",
        category: "Synteny & Whole-Genome Alignment (GPL-3 Drop-In)",
        tool_name: "Minimap2 / MUMmer4 / SyRI / GENESPACE",
        binary_cmd: "minimap2",
        spdx_license: "MIT / GPL-3.0",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite k-mer synteny dotplot & collinear block detector",
        description: "Detects chromosomal collinearity, inversions, and segmental duplications (replaces academic-only MCScanX).",
    },
    ToolAdapterEntry {
        engine: "genome",
        category: "Non-B DNA & Chromatin Architecture",
        tool_name: "G4Hunter / ZHunt / Cooler / MACS3",
        binary_cmd: "macs3",
        spdx_license: "GPL-3.0 / BSD-3-Clause",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite native Rayon G4Hunter G-Quadruplex, Z-DNA & CpG island scanner",
        description: "Scans for G-quadruplexes, Z-DNA alternating YR tracts, i-motifs, and chromatin accessibility peaks.",
    },
    // Engine 3: Genetic Target Identification & Drug Discovery (`bio-suite target`)
    ToolAdapterEntry {
        engine: "target",
        category: "Statistical Fine-Mapping & Colocalization",
        tool_name: "SuSiE / coloc / PolyFun",
        binary_cmd: "Rscript",
        spdx_license: "MIT / GPL-3.0",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite Wakefield ABF Bayesian fine-mapping PIP & coloc PP.H4 engine",
        description: "Narrows GWAS loci to 95% causal credible sets and colocalizes eQTL/pQTL signals (replaces restricted FINEMAP).",
    },
    ToolAdapterEntry {
        engine: "target",
        category: "Mendelian Randomization & TWAS",
        tool_name: "TwoSampleMR / S-PrediXcan / FUSION",
        binary_cmd: "Rscript",
        spdx_license: "MIT / GPL-3.0",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite Inverse-Variance Weighted (IVW) & MR-Egger causal estimator",
        description: "Estimates causal effects of genetically predicted gene/protein levels on disease phenotypes.",
    },
    ToolAdapterEntry {
        engine: "target",
        category: "CRISPR Screens & Isoform Splicing (GPL-3 Drop-In)",
        tool_name: "MAGeCK / BAGEL2 / rMATS-turbo / Whippet",
        binary_cmd: "rmats.py",
        spdx_license: "BSD-3 / GPL-3.0 / MIT",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite essentiality ranker & cryptic splice-site entropy scorer",
        description: "Identifies synthetic lethal CRISPR targets and disease-specific alternative splicing (replaces academic MAJIQ).",
    },
    ToolAdapterEntry {
        engine: "target",
        category: "Biomolecular Structure & Virtual Screening (GPL-3 Drop-In)",
        tool_name: "Boltz-1 / OpenFold / ESMFold / fpocket / AutoDock Vina / GNINA / OpenMM",
        binary_cmd: "vina",
        spdx_license: "MIT / Apache-2.0 / GPL-2.0+",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite pocket hydropathy/volume profiler & ligand physicochemical scorer",
        description: "100% open-source MIT/Apache-2.0 replacement for non-commercial AlphaFold 3 plus cryptic pocket & docking engines.",
    },
    // Engine 4: Complexity Theory & Information-Theoretic Genomics (`bio-suite complexity`)
    ToolAdapterEntry {
        engine: "complexity",
        category: "Shannon Entropy & Mutual Information",
        tool_name: "bio-suite native / Jellyfish / KMC3",
        binary_cmd: "bio-suite",
        spdx_license: "GPL-3.0-or-later",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite Rayon sliding-window k-mer Shannon entropy (H1..H6), GC3 & JSD engine",
        description: "Quantifies positional codon information content, k-mer block entropy, and Jensen-Shannon divergence.",
    },
    ToolAdapterEntry {
        engine: "complexity",
        category: "Algorithmic (Kolmogorov) Complexity",
        tool_name: "bio-suite native LZ76/LZ77 / GeCo3 / Jarvis",
        binary_cmd: "bio-suite",
        spdx_license: "GPL-3.0-or-later",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite native Lempel-Ziv 1976 (LZ76) exhaustive dictionary parser & NCD",
        description: "Approximates Kolmogorov sequence complexity and normalized compression distance across genomic regimes.",
    },
    ToolAdapterEntry {
        engine: "complexity",
        category: "Long-Range Correlations & Fractal DNA Walks",
        tool_name: "bio-suite native CGR / Hurst / DFA",
        binary_cmd: "bio-suite",
        spdx_license: "GPL-3.0-or-later",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite Purine-Pyrimidine DNA walk Detrended Fluctuation Analysis (DFA alpha) & R/S Hurst",
        description: "Measures scale-free long-range power-law correlations across coding, intronic, and intergenic regions.",
    },
    ToolAdapterEntry {
        engine: "complexity",
        category: "Formal Grammars & Topological Network Complexity",
        tool_name: "HMMER / Infernal / igraph / graph-tool",
        binary_cmd: "hmmsearch",
        spdx_license: "BSD-3-Clause / GPL-2.0+ / LGPL-3.0",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite spectral graph Laplacian & Von Neumann entropy calculator",
        description: "Parses Chomsky hierarchy RNA/protein grammars and quantifies 3D chromatin / phylogenetic network entropy.",
    },
    // Pipeline 3 (`hgt-scan`) & Pipeline 4 (`bgc-mine`) & Pipeline 5 (`pseudo-loss`) & Pipeline 6 (`pan-cryptic`) Adapters
    ToolAdapterEntry {
        engine: "genome",
        category: "Horizontal Gene Transfer & Mobile Resistome (Pipeline 3: hgt-scan)",
        tool_name: "tRNAscan-SE 2.0 / VirSorter2 / PhiSpy / ISEScan / MOB-suite / AMRFinderPlus / RGI",
        binary_cmd: "amrfinder",
        spdx_license: "GPL-3.0 / Public-Domain / Apache-2.0",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 256-dim tetranucleotide KL-divergence, 64-codon RSCU & attB/resistome scanner",
        description: "Detects horizontally transferred genomic islands via k-mer KL-divergence, tRNA/integrase flanks, and CARD/AMRFinderPlus cargo.",
    },
    ToolAdapterEntry {
        engine: "target",
        category: "Cryptic BGC, Self-Resistance & AMP Mining (Pipeline 4: bgc-mine)",
        tool_name: "Prodigal-gv / antiSMASH 7 / GECCO / BiG-SCAPE / ARTS / Macrel / ESMFold",
        binary_cmd: "antismash",
        spdx_license: "GPL-3.0 / AGPL-3.0 / BSD-3-Clause / MIT",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite MIBiG domain-Jaccard dereplicator, ARTS paralog scanner & Eisenberg mu_H AMP profiler",
        description: "Mines cryptic NRPS/PKS/RiPP clusters, dereplicates against MIBiG, and pairs scaffolds with duplicated self-resistance targets.",
    },
    ToolAdapterEntry {
        engine: "phylo",
        category: "Syntenic Pseudogenization & Selection Relaxation (Pipeline 5: pseudo-loss)",
        tool_name: "TOGA / Miniprot / MACSE / HyPhy RELAX / bcftools / gnomAD pLoF",
        binary_cmd: "hyphy",
        spdx_license: "MIT / GPL-3.0 / CeCILL",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite codon frameshift/stop-gain caller, Nei-Gojobori dN/dS & HyPhy RELAX k estimator",
        description: "Distinguishes true evolutionary gene loss (omega -> 1.0) from sequencing artifacts and identifies safe human knockout drug targets.",
    },
    ToolAdapterEntry {
        engine: "complexity",
        category: "Pangenome Graph Entropy & Cryptic Allosteric Pockets (Pipeline 6: pan-cryptic)",
        tool_name: "odgi / panacus / tskit / Modkit / IsoQuant / Boltz-1 / BioEmu / PocketMiner",
        binary_cmd: "odgi",
        spdx_license: "MIT / GPL-3.0 / Apache-2.0",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite GFA walk Shannon entropy (H_walk), Ollivier-Ricci curvature, H_iso & RIN betweenness engine",
        description: "Links pangenome variation graph bubble entropy and 5mC/scIso-Seq phase-switching to isoform-specific cryptic allosteric pockets.",
    },
];

/// Intersect high-complexity sliding windows with genomic annotations (genes, CDS, repeats)
/// to produce human- and machine-readable structural complexity discoveries (e.g., Region X32).
pub fn intersect_high_complexity_windows(
    windows: &[SlidingWindowRecord],
    annotations: &[GenomicAnnotation],
    z_threshold: f64,
) -> Vec<AnnotatedComplexityRegion> {
    let mut ranked: Vec<&SlidingWindowRecord> = windows
        .iter()
        .filter(|w| w.structural_complexity_z >= z_threshold)
        .collect();
    ranked.sort_by(|a, b| {
        b.structural_complexity_z
            .partial_cmp(&a.structural_complexity_z)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    ranked
        .into_iter()
        .enumerate()
        .map(|(idx, win)| {
            let overlaps: Vec<&GenomicAnnotation> = annotations
                .iter()
                .filter(|ann| ann.start < win.end && ann.end > win.start)
                .collect();

            let cds_count = overlaps
                .iter()
                .filter(|a| a.feature_type == "CDS" || a.feature_type == "gene")
                .count();
            let repeat_count = overlaps
                .iter()
                .filter(|a| a.feature_type == "repeat_region")
                .count();
            let feature_labels: Vec<String> = overlaps
                .iter()
                .map(|a| format!("{} [{} {}:{}..{}]", a.symbol, a.feature_type, a.strand, a.start, a.end))
                .collect();

            let region_id = format!("Region-X{:02}", idx + 1);
            let narrative = format!(
                "{} ({}..{} bp, Z = {:+.2}) overlaps {} protein-coding locus/loci and {} repeat element(s) ({}), exhibiting elevated GC ({:.1}%), G4Hunter propensity ({:.2}), and normalized Lempel-Ziv complexity ({:.3}).",
                region_id,
                win.start,
                win.end,
                win.structural_complexity_z,
                cds_count,
                repeat_count,
                if feature_labels.is_empty() {
                    "unannotated intergenic space".to_string()
                } else {
                    feature_labels.join("; ")
                },
                win.gc_percent,
                win.g4hunter_score,
                win.lz76_normalized
            );

            AnnotatedComplexityRegion {
                region_id,
                window_id: win.window_id.clone(),
                start: win.start,
                end: win.end,
                structural_complexity_z: win.structural_complexity_z,
                gc_percent: win.gc_percent,
                cpg_obs_exp: win.cpg_obs_exp,
                g4hunter_score: win.g4hunter_score,
                zdna_score: win.zdna_score,
                lz76_normalized: win.lz76_normalized,
                dfa_alpha: win.dfa_alpha,
                overlapping_cds_count: cds_count,
                overlapping_repeat_count: repeat_count,
                overlapping_features: feature_labels,
                narrative_summary: narrative,
            }
        })
        .collect()
}

/// Wakefield Approximate Bayes Factor (ABF) Bayesian fine-mapping credible set calculator (`bio-suite target`).
pub fn compute_wakefield_fine_mapping(
    variants: &[(String, usize, f64, f64)], // (rsid, pos, beta, se)
    prior_variance_w: f64,
) -> Vec<FineMappedVariantCredibleMember> {
    let mut raw: Vec<(String, usize, f64, f64, f64, f64)> = Vec::with_capacity(variants.len());
    for (rsid, pos, beta, se) in variants {
        let v = se * se;
        let z = beta / se.max(1e-9);
        let r = prior_variance_w / (v + prior_variance_w);
        let log_abf = 0.5 * (1.0 - r).ln() + 0.5 * r * z * z;
        let abf = log_abf.exp();
        raw.push((rsid.clone(), *pos, *beta, *se, z, abf));
    }
    let sum_abf: f64 = raw.iter().map(|x| x.5).sum::<f64>().max(1e-12);
    let mut out: Vec<FineMappedVariantCredibleMember> = raw
        .into_iter()
        .map(|(rsid, pos, beta, se, z, abf)| FineMappedVariantCredibleMember {
            rsid,
            position: pos,
            beta: round4(beta),
            se: round4(se),
            z_score: round4(z),
            wakefield_abf: round4(abf),
            posterior_inclusion_prob: round4(abf / sum_abf),
        })
        .collect();
    out.sort_by(|a, b| {
        b.posterior_inclusion_prob
            .partial_cmp(&a.posterior_inclusion_prob)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    out
}
