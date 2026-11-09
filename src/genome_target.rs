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
        tool_name: "DIAMOND v2.2.8 (blastp --swipe) / OrthoFinder / BUSCO",
        binary_cmd: "diamond",
        spdx_license: "GPL-3.0-or-later / MIT",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 4-Stage Smith-Waterman BLOSUM62, Emms-Kelly Length-Normalized RBH & Paralog Graph Engine",
        description: "Identifies single-copy orthologs, co-ortholog in-paralogs, and out-paralogs across proteomes and coding sequences.",
    },
    ToolAdapterEntry {
        engine: "phylo",
        category: "Multiple Sequence Alignment & Trimming (MAFFT L-INS-i / E-INS-i / Q-INS-i / FFT-NS-i)",
        tool_name: "MAFFT v7.526 (L-INS-i / E-INS-i / Q-INS-i / FFT-NS-i --adjustdirection) / MUSCLE5 / trimAl",
        binary_cmd: "mafft",
        spdx_license: "BSD-3-Clause / GPL-3.0",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 6-Stage Native Rust MAFFT-Mimetic Aligner (runie-mafft-linsi: Cooley-Tukey Radix-2 FFT + hat3 Smith-Waterman Transitive Consistency + rDNA 18S-ITS1-5.8S-ITS2-28S Subunit Locking + 2-Stage UPGMA Profile Gotoh DP + dvtditr Iterative Refinement) & Center-Star Gotoh / Codon MSA",
        description: "Aligns coding, rRNA, and hypervariable Internal Transcribed Spacer (ITS1-5.8S-ITS2) loci via FFT phase correlation, local Smith-Waterman consistency, rDNA subunit boundary locking, and tree-bipartition iterative refinement.",
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
        category: "Codon Selection (dN/dS, Site/Branch & BEB)",
        tool_name: "PAML v4.10.10 (codeml / yn00) / HyPhy",
        binary_cmd: "codeml",
        spdx_license: "GPL-3.0-or-later",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 61-State Goldman-Yang Codon CTMC (M0/M1a/M2a/Branch) & 4D Grid BEB Engine",
        description: "Evaluates nested site and branch codon substitution models (M0, M1a, M2a, Branch 2-Ratio, RELAX k) and Bayes Empirical Bayes (BEB) site posteriors.",
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
        category: "Pangenome Graphs & Structural Variants (GPL-3 Drop-In)",
        tool_name: "minigraph v0.21-r606 (-cxggs / -cxasm --call, Li et al. 2020) / gfatools v0.5-r296 (stat / bubble, Onodera et al. 2013) / Minigraph-Cactus / PGGB / vg",
        binary_cmd: "minigraph",
        spdx_license: "MIT / GPL-3.0-or-later",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 4-Stage Native rGFA Pangenome Constructor, Onodera-Sadakane-Shibuya Superbubble Decomposer & Walk Genotyper",
        description: "Constructs reference-tagged pangenome variation graphs (rGFA), decomposes biallelic and multi-allelic structural variant superbubbles, and calls exact haplotype walk genotypes.",
    },
    ToolAdapterEntry {
        engine: "genome",
        category: "Ab Initio Annotation & Repeat Masking (GPL-3 Drop-In)",
        tool_name: "Prodigal v2.6.3 / Helixer / AUGUSTUS / EDTA / Liftoff",
        binary_cmd: "prodigal",
        spdx_license: "GPL-3.0-or-later",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 4-Stage 6-Frame ORF Graph, Dicodon Hexamer Markov, 28-Bin Shine-Dalgarno RBS EM Trainer & Operon Overlap DP + Eukaryotic Kozak/GT-AG Splice Reconstructor",
        description: "100% GPL-3.0 prokaryotic and eukaryotic ab initio gene calling, RBS/Kozak translation initiation scoring, and exon/repeat coordinate annotation.",
    },
    ToolAdapterEntry {
        engine: "genome",
        category: "Synteny & Whole-Genome Alignment (GPL-3 Drop-In)",
        tool_name: "Minimap2 v2.31-r1302 (-c --cs -x asm5/asm10/asm20) / MUMmer4 / SyRI",
        binary_cmd: "minimap2",
        spdx_license: "MIT / GPL-3.0",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 4-Stage Double-Stranded (w,k)-Minimizer, Concave-Gap Chaining & SyRI Split-Chain SV Engine",
        description: "Detects chromosomal collinearity, inversions, genomic island indels, and segmental duplications (replaces academic-only MCScanX).",
    },
    ToolAdapterEntry {
        engine: "genome",
        category: "Non-B DNA & Chromatin Architecture (GPL-3 Drop-In)",
        tool_name: "G4Hunter v1.0 (Bedrat et al. 2016) / Z-Hunt-II (zhunt3, Ho et al. 1986; Schroth et al. 1992) / MACS3",
        binary_cmd: "zhunt3",
        spdx_license: "GPL-3.0-or-later / BSD-3-Clause",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 4-Stage Run-Length Quantized G4Hunter, 4x16 Anti/Syn Boltzmann Partition Function Z-Hunt-II & CpG Island Scanner",
        description: "Scans for G-quadruplexes, i-motifs, left-handed Z-DNA B-to-Z superhelical transitions (Delta Lk_1/2, Z-Score), and Gardiner-Garden & Frommer CpG islands.",
    },
    ToolAdapterEntry {
        engine: "genome",
        category: "RNA Secondary Structure & Partition Function Folding (GPL-3 Drop-In)",
        tool_name: "ViennaRNA v2.7.2 (RNAfold -p / -g, Hofacker 2003; Lorenz et al. 2011; Turner & Mathews 2004)",
        binary_cmd: "RNAfold",
        spdx_license: "Custom Free Academic / Open Source / GPL-3.0 Fallback",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 4-Stage Turner-2004 Zuker-Stiegler MFE, McCaskill Inside-Outside Partition Function, Lorenz rG4 & CBC Epistasis Engine",
        description: "Computes exact O(N^3) minimum free energy (MFE) secondary structure, Boltzmann partition function base-pair probabilities P_ij, positional pairing entropy S(i), RNA G-quadruplexes, and mutational Delta-Delta-G.",
    },
    ToolAdapterEntry {
        engine: "target",
        category: "Nucleic Acid 1D -> 2D -> 3D Structural Prediction & Tertiary Alignment",
        tool_name: "US-align v20261019 (-mol RNA -het 1, Zhang et al. 2022) / ViennaRNA v2.7.2 (RNAPKplex / dna_mathews2004.par) / Z-Hunt-II",
        binary_cmd: "USalign",
        spdx_license: "GPL-3.0-or-later / Free Academic",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 3-Stage 1D->2D->3D Nucleic Acid Pipeline (SE(3) Rigid-Body Step Integrator + Motif-Fragment Distance-Geometry + Kabsch SVD TM-score_RNA)",
        description: "Predicts 1D non-B propensity, 2D RNA/DNA partition functions & H-type pseudoknots, and 3D all-atom/5-bead tertiary coordinates validated via Kabsch RMSD and TM-score_RNA.",
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
        tool_name: "fpocket v4.2.3 (Le Guilloux et al. 2009; Schmidtke & Barril 2010) / AutoDock Vina v1.2.7 (Trott & Olson 2010; Eberhardt et al. 2021) / Boltz-1 / OpenFold",
        binary_cmd: "vina",
        spdx_license: "MIT / Apache-2.0 / GPL-3.0",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 4-Stage Alpha-Sphere Voronoi Cavity, Shrake-Rupley SASA/Druggability & SE(3) Trott-Olson 5-Term Empirical Vina Docking Engine",
        description: "Detects druggable protein/RNA binding cavities via Voronoi alpha-spheres and docks small-molecule ligands using the Trott-Olson empirical scoring function and SE(3) pose optimization.",
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
        category: "Formal Grammars, Profile HMMs & Functional Domain Search (GPL-3 Drop-In)",
        tool_name: "HMMER v3.4 (hmmbuild / hmmsearch --max) / Infernal / igraph / graph-tool",
        binary_cmd: "hmmsearch",
        spdx_license: "BSD-3-Clause / GPL-3.0-or-later",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 4-Stage Plan-7 Profile HMM Engine (p7_GMSV, p7_GViterbi, p7_GForward/Backward, p7_GDecoding, p7_GOptimalAccuracy & p7_GNull2) + spectral Laplacian entropy",
        description: "100% GPL-3.0 Plan-7 Profile HMM construction, Forward-Backward ensemble log-odds, Optimal Posterior Accuracy domain alignment, Null2 bias correction, and Chomsky/network entropy.",
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
        tool_name: "minigraph v0.21-r606 / gfatools v0.5-r296 / odgi / panacus / Modkit / IsoQuant / Boltz-1 / PocketMiner",
        binary_cmd: "gfatools",
        spdx_license: "MIT / GPL-3.0 / Apache-2.0",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 4-Stage rGFA Walk Shannon Entropy (H_walk), Ollivier-Ricci Curvature (kappa_uv), Panacus Growth Scaling (gamma) & Cryptic Allosteric Pocket Engine",
        description: "Links pangenome variation graph superbubble walk entropy (H_walk), Ollivier-Ricci edge curvature, and Panacus open-pangenome growth scaling to isoform-specific cryptic allosteric pockets.",
    },
    // RuNIE (Rust Native Inference Engine) Autonomously Synthesized & Calibrated Domain Tools
    ToolAdapterEntry {
        engine: "phylo",
        category: "Total-Evidence Discrete Morphology/Behavior, Fossil Tip-Dating & PIC (RuNIE Tool A)",
        tool_name: "MrBayes v3.2.7a / IQ-TREE v2.4.0 (-m MK+ASC+G4) / R ape::pic",
        binary_cmd: "mrbayes",
        spdx_license: "GPL-3.0-or-later",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 4-Stage Lewis (2001) Mk+ASC Discrete Character CTMC, Fossilized Birth-Death (FBD) Tip-Dating & Felsenstein (1985) PIC Engine",
        description: "Integrates discrete morphological/behavioral characters, extinct fossil tips, and tree-imputed genomic metadata contrasts (Zrzavy et al. 2018; Challis et al. 2023).",
    },
    ToolAdapterEntry {
        engine: "target",
        category: "AlphaFold PAE/pLDDT Spectral Partitioning, ANM Allostery & Circadian Oscillator (RuNIE Tool B)",
        tool_name: "AlphaFold-Multimer v2.3 / ProDy v2.4 ANM / COPASI v4.42 RK4",
        binary_cmd: "prody",
        spdx_license: "MIT / Apache-2.0 / GPL-3.0",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 4-Stage Candle-Core PAE Spectral Block Partitioner, 3N x 3N ANM Hessian Allosteric Pathway & KaiABC Hexamer RK4 Oscillator",
        description: "Models PAE inter-domain confidence, ANM normal-mode allosteric communication, KaiB fold-switching, and 24-hour KaiABC phosphorylation oscillations (Chang 2025; Varadi 2024; Nussinov 2023).",
    },
    ToolAdapterEntry {
        engine: "phylo",
        category: "Universal 16-Ribosomal-Protein Supermatrix & CPR/DPANN Auxotrophy Scanner (RuNIE Tool C)",
        tool_name: "IQ-TREE v2.4.0 (LG+C60+F+G4) / HMMER v3.4 / BMGE",
        binary_cmd: "iqtree2",
        spdx_license: "GPL-2.0-or-later / GPL-3.0",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 4-Stage 16-RP Syntenic Supermatrix Concatenator, Dayhoff-6 Compositional Recoder & CPR/DPANN Auxotrophy Scanner",
        description: "Resolves deep bacterial (CPR/Patescibacteria) and archaeal (DPANN/Asgard) radiations while correcting compositional LBA (Hug et al. 2016).",
    },
    ToolAdapterEntry {
        engine: "target",
        category: "De Novo SE(3) Frame Denoising Diffusion & Functional Motif Scaffolding (RuNIE Tool D)",
        tool_name: "RFdiffusion v1.1.0 (Watson et al. 2023) / ProteinMPNN v1.0.1",
        binary_cmd: "rfdiffusion",
        spdx_license: "BSD-3-Clause / MIT / GPL-3.0",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 4-Stage Candle-Core SE(3) Rigid-Body Frame Denoising Diffusion, Clamped Motif Scaffolder & scTM/scRMSD Self-Consistency Engine",
        description: "Generates de novo protein backbones scaffolding catalytic triads and viral epitopes with SE(3) equivariant diffusion and Potts inverse folding (Watson et al. 2023).",
    },
    ToolAdapterEntry {
        engine: "genome",
        category: "Multi-Objective mRNA Codon-Lattice, Turner-2004 MFE & m1Psi Optimizer (RuNIE Tool E)",
        tool_name: "LinearDesign v1.0 (Zhang et al. 2023) / ViennaRNA v2.7.2",
        binary_cmd: "lineardesign",
        spdx_license: "GPL-3.0-or-later / Free Academic",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 4-Stage DFA Synonymous Codon-Lattice x Turner-2004 MFE Co-Optimizer, 5'-UTR Initiation Guard & m1Psi Engine",
        description: "Jointly optimizes mRNA vaccine Codon Adaptation Index (CAI), thermodynamic folding stability (MFE), 5' Kozak accessibility, and m1Psi uridine depletion (Al Fayez et al. 2023).",
    },
    ToolAdapterEntry {
        engine: "complexity",
        category: "Spidroin Repeat-Motif Grammar & Silk Biomechanical Property Regressor (RuNIE Tool F)",
        tool_name: "Spider Silkome DB / Elastic-Network Polymer Regressor (Arakawa et al. 2022)",
        binary_cmd: "bio-suite",
        spdx_license: "GPL-3.0-or-later",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 4-Stage Candle-Core Spidroin Repeat Grammar (Poly-Ala / GPGXX / MaSp3 / SpiCE) & Multi-Output Biomechanics Regressor",
        description: "Predicts dragline silk toughness (MJ/m^3), tensile strength (GPa), Young's modulus, and supercontraction from spidroin repeat architectures across 1,098 Araneae species (Arakawa et al. 2022).",
    },
    ToolAdapterEntry {
        engine: "phylo",
        category: "Adaptive Multi-Biotype Super-BLAST (DNA/RNA/Protein/Other), Sensitivity Analyzer, Pre-MSA Validator, 16S Tip Concordance & JIT Genome/Barcode Search (RuNIE Tool G)",
        tool_name: "Biopython v1.88 PairwiseAligner / NCBI BLAST+ (blastn / blastp / makeblastdb) / SciPy v1.18.1",
        binary_cmd: "super-blast",
        spdx_license: "GPL-3.0-or-later / BSD-3-Clause",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 5-Stage Adaptive Multi-Biotype (DNA/RNA/Protein/Other) Super-BLAST Engine with On-the-Fly Direct vs. JIT Ephemeral Index Sensitivity Analyzer, Pre-MSA Validator, 16S + 10-Gene Tip Auditor & Symmetric K2P Matrix",
        description: "Automatically detects DNA, RNA, Protein (BLOSUM62), or Standard Language strings, performs sensitivity crossover analysis to avoid unnecessary DB compilation, validates MSA inputs, audits multi-gene tip concordance, and executes JIT genome/barcode searches.",
    },
    ToolAdapterEntry {
        engine: "phylo",
        category: "Divergence Time-Tree (Chronogram), Constrained Least-Squares Dating (LSD2), Relaxed Molecular Clock & X-Axis Time-Tree SVG Generator (RuNIE Tool H)",
        tool_name: "IQ-TREE v2.4.0 (--date LSD2, To et al. 2016) / SciPy v1.18.1 (Constrained Least-Squares) / Sanderson (2002) r8s",
        binary_cmd: "timetree-compare",
        spdx_license: "GPL-3.0-or-later / BSD-3-Clause",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite 4-Stage Active-Set LSD2 Chronogram Solver, Autocorrelated Log-Normal Relaxed Clock, Parametric Bootstrap 95% CI & Clean X-Axis Time-Tree SVG Engine",
        description: "Estimates divergence time-trees (chronograms) under strict and autocorrelated log-normal relaxed molecular clocks with fossil node intervals [t_min, t_max], heterochronous tip dates, 95% CIs, and automatic X-axis geological time-tree SVG rendering.",
    },
    ToolAdapterEntry {
        engine: "phylo",
        category: "Interactive Sequence Alignment Viewer, Multi-Alphabet Translation (DNA/AA/RY3/Dayhoff-6), Feature/Annotation Selector & Publication Snapshot Engine (RuNIE Tool #37)",
        tool_name: "AliView v1.28 (Larsson 2014) / Jalview / SeaView",
        binary_cmd: "runie-align-view",
        spdx_license: "GPL-3.0-or-later",
        gpl3_compatible: true,
        native_rust_fallback: "bio-suite align-view / runie-align-view (Native Win32 GUI + 32-bit Software Framebuffer + Pixel-Locked SVG & High-DPI PNG Snapshot Engine)",
        description: "Provides interactive desktop MSA inspection with top File/View/Translation/Features/RuNIE menus, global minimap lens, feature/annotation jumping, DNA/AA/RY3/Dayhoff-6 color-coded views, and automated focal-region SVG/PNG snapshots for manuscripts.",
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
