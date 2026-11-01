# CRV Science Labs (`bio-suite` v0.3.0)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language: Rust](https://img.shields.io/badge/Language-Rust%202021-orange.svg)](https://www.rust-lang.org/)
[![Portal: crvscience.com/labs](https://img.shields.io/badge/Portal-crvscience.com%2Flabs-00bcd4.svg)](https://crvscience.com/labs)

**`bio-suite`** is the open-source Rust command-line toolkit powering [CRV Science Labs](https://crvscience.com/labs). It provides fast, self-contained pipelines for molecular systematics, comparative genomics, pangenome variation graphs, nucleic acid 1D/2D/3D structural prediction, biomolecular pocket detection and small-molecule docking, sequence complexity analysis, and human disease gene mapping.

Every pipeline runs in native multi-threaded Rust out of the box—requiring zero external software—while offering optional command-line flags and benchmark harnesses to compare results directly against standard public bioinformatics tools (`MAFFT`, `jModelTest`, `IQ-TREE 2`, `ASTRAL`, `PAML codeml`, `DIAMOND`, `Minimap2`, `Prodigal`, `HMMER3`, `G4Hunter`, `Z-Hunt-II`, `ViennaRNA`, `RNAPKplex`, `US-align`, `fpocket`, `AutoDock Vina`, `minigraph`, and `gfatools`).

---

## Quick Start

```bash
git clone https://github.com/bpwhite/crvscience-labs.git
cd crvscience-labs
cargo test
cargo build --release
```

The standalone executable is built at `target/release/bio-suite` (`target/release/bio-suite.exe` on Windows).

---

## Native Tools, Public Counterparts & Validation Roadmap

### 1. Head-to-Head Calibrated Native Engines vs. Public Counterparts (`✅ Validated`)

All 13 core tool comparisons below have been calibrated head-to-head against compiled external reference binaries across real NCBI RefSeq and experimental PDB datasets:

| # | Native Tool (`bio-suite` Command & Module) | Public Counterpart(s) | Validation Status & Concordance Summary | Function of the Tool |
| :---: | :--- | :--- | :--- | :--- |
| **1** | **`align-compare`** (`src/phylo.rs`, `src/msa_calibration.rs`) | `MAFFT v7.526` (`--aligner mafft`), `MUSCLE5`, `trimAl` | ✅ **Validated** (`SP = 1.0000`, `TC = 1.0000` on coding panels; `r = 0.999999` on `K2P` distances; `110×` speedup) | Aligns protein-coding and ribosomal RNA sequences using codon-anchored triplet alignment and 9-mer anchor-segmented Gotoh 3-state affine-gap dynamic programming. |
| **2** | **`modeltest-compare`** & **`phylo-2x2-compare`** (`src/modeltest.rs`) | `jModelTest v2.1.9` / `PhyML v3.0` (`--modeltest jmodeltest`) | ✅ **Validated** (`|\Delta \ln L| = 0.0086` log-likelihood units across 24 models; `100%` top `BIC`/`AICc` family agreement) | Evaluates 24 nucleotide substitution models (`JC`, `F81`, `K80`, `HKY`, `SYM`, `GTR` × `[Uniform, +I, +G4, +I+G4]`) via column-compressed Felsenstein pruning, Padé `[6/6]` matrix exponentials, and `AICc`/`BIC` ranking. |
| **3** | **`tree-compare`** (`src/tree_builder.rs`) | `IQ-TREE v2.4.0` (`--tree-builder iqtree`), `RAxML-ng` | ✅ **Validated** (54-taxon vertebrate mitogenome ML tree: `100%` crown-clade recovery, `RF = 0.039`, `12.4×` speedup) | Infers maximum-likelihood phylogenetic trees using a 4-stage Multi-Distance BIONJ starter pool, exact Yang discrete $\Gamma_4$ site rates, inside-outside pulley branch optimization, and quartet-polished NNI search. |
| **4** | **`coalescent-compare`** (`src/coalescent.rs`) | `ASTRAL v5.7.8` (`--coalescent astral`), `Dsuite` | ✅ **Validated** (1,064 gene trees across 3 cohorts: `100%` topology match `RF = 0.0000`, quartet score `r = 1.000000`) | Infers species trees under incomplete lineage sorting (ILS) from multi-locus gene trees via exact quartet tripartition dynamic programming, Bayesian local posterior probabilities (`LPP`), coalescent branch lengths, and `ABBA-BABA` introgression tests (`Patterson's D`). |
| **5** | **`selection-compare`** (`src/codon_selection.rs`) | `PAML v4.10.10` (`codeml`, `--selection paml`), `HyPhy` | ✅ **Validated** (Goldman-Yang `M0`/`M1a`/`M2a`/Branch across 4 cohorts: `ln L` `r = 0.999999`, `100%` BEB site match) | Evaluates 61-state standard and 60-state vertebrate mitochondrial Goldman-Yang codon CTMC models (`M0`, `M1a`, `M2a`, Branch 2-Ratio, `RELAX k`) and 10,000-cell 4D grid Bayes Empirical Bayes (`BEB`) positive selection site posteriors ($d_N/d_S = \omega$). |
| **6** | **`orthology-compare`** (`src/orthology.rs`) | `DIAMOND v2.2.8` (`blastp --swipe`), `OrthoFinder`, `BUSCO` | ✅ **Validated** (2,646 pairs across 4 cohorts: bit-score `r = 0.99989`, `100%` 1:1 RBH & in/out-paralog classification) | Identifies single-copy orthologs, co-ortholog in-paralogs, and out-paralogs across proteomes using exact Smith-Waterman `BLOSUM62` local alignment, Karlin-Altschul $E$-values, and Emms-Kelly length-normalized reciprocal best hits (`RBH`). |
| **7** | **`synteny-compare`** (`src/synteny.rs`) | `Minimap2 v2.31-r1302` (`-x asm5/asm10/asm20`), `MUMmer4`, `SyRI` | ✅ **Validated** (376.9 kb across 4 cohorts: `100%` SV call concordance, BLAST ID MAE `0.06%`, `de` MAE `0.0004`) | Maps whole-chromosome collinearity, inversions, translocations, and genomic island insertion breakpoints using double-stranded $(w,k)$-minimizer indexing, concave-gap dynamic programming chaining, and split-chain structural variant calling. |
| **8** | **`gene-compare`** (`src/gene_caller.rs`) | `Prodigal v2.6.3`, `AUGUSTUS`, `Helixer` | ✅ **Validated** (396 kb across 6 bacterial genomes: `99.71%` 3′-stop match, `cscore` `r = 0.99995`, `sscore` `r = 0.9988`) | Performs *ab initio* prokaryotic and eukaryotic gene prediction via a 6-frame open reading frame (ORF) graph, dicodon hexamer Markov coding potentials, 28-bin Shine-Dalgarno RBS expectation-maximization training, and operon overlap dynamic programming. |
| **9** | **`hmm-compare`** (`src/profile_hmm.rs`) | `HMMER v3.4` (`hmmbuild` / `hmmsearch --max`), `Infernal` | ✅ **Validated** (5 Pfam/AMR families, 30 targets: `100%` sensitivity/specificity, bitscore `r = 0.99978`, IoU `0.9936`) | Builds Plan-7 profile hidden Markov models and searches sequence databases for conserved protein/RNA domains using multi-segment ungapped `MSV`, Viterbi, Forward-Backward ensemble log-odds, posterior optimal accuracy alignment, and `Null2` composition bias correction. |
| **10** | **`nonb-dna-compare`** (`src/nonb_dna.rs`, `src/complexity.rs`) | `G4Hunter v1.0`, `Z-Hunt-II` (`zhunt3`) | ✅ **Validated** (500 kb across 6 cohorts: G4Hunter `r = 1.000000`, Z-Hunt-II $\Delta\text{Lk}_{1/2}$ `MAE = 0.022`, Z-score `r = 0.9997`) | Scans chromosomes for non-canonical DNA conformations, computing run-length quantized G-quadruplex/i-motif propensity, exact $4\times 16$ *anti*/*syn* Boltzmann partition-function B-to-Z DNA superhelical transitions ($\Delta\text{Lk}_{1/2}$, Z-score), and CpG islands. |
| **11** | **`rna-fold-compare`** & **`na-1d2d3d-compare`** (`src/rna_fold.rs`, `src/na_structure_1d2d3d.rs`) | `ViennaRNA v2.7.2` (`RNAfold -p -g`, `RNAPKplex`, `dna_mathews2004.par`), `US-align v20261019` | ✅ **Validated** (12 RNA families: MFE `r = 0.9992`, `F1 = 0.9854`; 7 PDB 3D benchmarks: mean `RMSD = 1.48 Å`, `TM-score = 0.724`) | Predicts 1D non-B propensity, 2D Turner-2004 MFE/McCaskill partition functions ($P_{ij}$), RNA G-quadruplexes (`rG4`), and H-type pseudoknots, and folds 5-atom coarse-grained 3D `.pdb` coordinates for DNA and RNA validated via Kabsch SVD `RMSD` and `TM-score`. |
| **12** | **`pocket-dock-compare`** (`src/pocket_dock.rs`) | `fpocket v4.2.3`, `AutoDock Vina v1.2.7` | ✅ **Validated** (10 PDB complexes + negative controls: druggability `r = 0.9718`, Vina $\Delta G_{\text{dock}}$ `r = 0.9906`, `MAE = 0.259 kcal/mol`, `100%` sub-`2.0 Å` redocking) | Detects protein binding pockets via 3D Voronoi $\alpha$-sphere clustering and Schmidtke–Barril druggability scoring ($D_{\text{drug}}$, polar/apolar SASA), and performs $\text{SE}(3)$ rigid-body + torsional small-molecule docking under the AutoDock Vina empirical free-energy potential ($\Delta G_{\text{dock}}$). |
| **13** | **`pangenome-compare`** (`src/pangenome.rs`) | `minigraph v0.21-r606` (`-cxggs`, `-cxasm --call`), `gfatools v0.5-r296` (`stat`, `bubble`) | ✅ **Validated** (10 loci, 60 haplotypes, 1.81 Mb: `100%` segment/link match, length `r = 1.000000`, `100%` bubble `F1` across 50 SVs, `100%` walk call match, `16.7×` speedup) | Constructs reference-tagged pangenome variation graphs (`rGFA`), decomposes biallelic and multi-allelic structural variant superbubbles via the Onodera–Sadakane–Shibuya algorithm, calls exact haplotype walk paths, and computes walk Shannon entropy ($H_{\text{walk}}$), Ollivier–Ricci curvature ($\kappa_{uv}$), and Panacus growth scaling ($\gamma$). |

---

### 2. Future Tools & Active External Validation Queue (`🔬 In Validation`)

The following native `bio-suite` engines and `BIO_SUITE_TOOL_REGISTRY` adapters are implemented in Rust and actively undergoing (or queued for) dedicated head-to-head external binary calibration benchmarks:

| Planned Benchmark Command & Module | Public Counterpart(s) Under Validation | Validation Status | Function & Planned Calibration Scope |
| :--- | :--- | :---: | :--- |
| **`hgt-scan` / `amr-compare`** (`src/pipelines.rs`, `src/genome_target.rs`) | `AMRFinderPlus v3.12`, `RGI v6` (`CARD`), `MOB-suite`, `tRNAscan-SE 2.0`, `ISEScan` | 🔬 **In Validation** | Detects horizontally transferred genomic islands and mobile antimicrobial resistance (`AMR`) plasmids via 256-dim tetranucleotide Kullback-Leibler divergence ($D_{KL}$), 64-codon `RSCU` shift, `tRNA`/`attB` integration flanks, and `Inc` replicon typing across ESKAPE pathogens. |
| **`bgc-mine` / `bgc-compare`** (`src/pipelines.rs`, `src/genome_target.rs`) | `antiSMASH 7`, `GECCO`, `BiG-SCAPE`, `ARTS`, `Macrel` | 🔬 **In Validation** | Mines microbial genomes for biosynthetic gene clusters (`NRPS`, `PKS`, `RiPP`), dereplicates against `MIBiG` reference clusters via domain-Jaccard similarity, detects duplicated self-resistance target paralogs, and profiles amphipathic antimicrobial peptides (`Eisenberg` $\mu_H$). |
| **`pseudo-loss` / `relax-compare`** (`src/pipelines.rs`, `src/codon_selection.rs`) | `HyPhy RELAX`, `TOGA`, `MACSE v2`, `bcftools`, `gnomAD v4 pLoF` | 🔬 **In Validation** | Identifies shared frameshifts, premature stop codons, splice-site disruptions, and relaxed purifying selection ($\omega \to 1.0$, `RELAX` selection intensity parameter $k$) across syntenic vertebrate coding alignments paired with human homozygous loss-of-function (`pLoF`) biobank tolerances. |
| **`co-evo` / `dca-compare`** (`src/dag_pipelines.rs`) | `plmDCA`, `CCMpred`, `EVcouplings`, `Infernal v1.1.5` | 🔬 **In Validation** | Computes pairwise mutual information (`MI`), Average Product Correction (`APC`), and $L_2$-regularized pseudolikelihood Potts direct-coupling parameters ($\|J_{ij}\|$) to detect 3D contact pairs and compensatory epistatic resistance networks in proteins and `23S` ribosomal RNA. |
| **`haplo-cluster` / `popgen-compare`** (`src/haplotype.rs`) | `PopART`, `DnaSP v6`, `Network`, `vcftools` | 🔬 **In Validation** | Collapses population sequences into unique haplotypes, computes haplotype/nucleotide diversity ($H_d$, $\pi$), Tajima's $D$, and Fu's $F_S$, builds polymorphic site dot-identity matrices, and constructs Kruskal/Bandelt Minimum Spanning Networks (`MSN`). |
| **`coi-barcode-cluster` / `abgd-compare`** (`src/phylo.rs`, `src/agent.rs`) | `ABGD` (Puillandre et al.), `ASAP`, `MEGA11`, `BOLD` | 🔬 **In Validation** | Computes `K2P`/`TN93` evolutionary distance matrices, Neighbor-Joining/UPGMA trees with parallel bootstrap support, and Automatic Barcode Gap Discovery (`ABGD`) prior intraspecific divergence partitions across `COX1`, `MT-RNR2`, and `ITS1–5.8S–ITS2` panels. |
| **`finemap-compare` (`target` Engine)** (`src/genome_target.rs`, `src/disease_panels.rs`) | `SuSiE`, `coloc`, `PolyFun`, `FINEMAP` | 🔬 **In Validation** | Computes Wakefield Approximate Bayes Factors (`ABF`), 95% causal credible sets, posterior inclusion probabilities (`PIP`), and Bayesian eQTL/pQTL colocalization posterior probabilities (`PP.H0`–`PP.H4`) across human GWAS disease loci. |
| **`mr-twas-compare` (`target` Engine)** (`src/genome_target.rs`) | `TwoSampleMR`, `S-PrediXcan`, `FUSION` | 🔬 **In Validation** | Estimates causal effects of genetically predicted gene expression and circulating protein levels (`pQTL`/`eQTL`) on human clinical phenotypes via Inverse-Variance Weighted (`IVW`), weighted median, and `MR-Egger` pleiotropy-robust regression. |
| **`splicing-crispr-compare` (`target` Engine)** (`src/genome_target.rs`, `src/pipelines.rs`) | `rMATS-turbo`, `Whippet`, `IsoQuant`, `MAGeCK`, `BAGEL2` | 🔬 **In Validation** | Quantifies alternative transcript isoform splicing entropy ($H_{\text{iso}}$, $\Delta\Psi$), MaxEntScan `5′`/`3′` splice-site strength, and synthetic-lethal paralog co-dependencies across CRISPR knockout screens. |
| **`assembly-compare` (`genome` Engine)** (`src/genome_target.rs`) | `hifiasm`, `Verkko`, `Flye`, `QUAST`, `Merqury` | 🔬 **In Validation** | Evaluates overlap-layout contig phasing, `N50`/`NG50`/`L50`/`auN` continuity statistics, $k$-mer completeness (`QV`), and telomere-to-telomere (`T2T`) chromosomal scaffold integrity. |
| **`kmer-complexity-compare` (`complexity` Engine)** (`src/complexity.rs`) | `Jellyfish 2`, `KMC3`, `GeCo3`, `Jarvis` | 🔬 **In Validation** | Benchmarks multi-threaded sliding-window $k$-mer Shannon block entropy ($H_1\dots H_6$), exhaustive Lempel-Ziv (`LZ76`) dictionary parsing, Normalized Compression Distance (`NCD`), and Purine-Pyrimidine DNA walk Detrended Fluctuation Analysis (`DFA` $\alpha$). |
| **`chrom-condensate` & `clonal-synth`** (`src/dag_pipelines.rs`) | `Cooler` / `HiCExplorer`, `SigProfilerExtractor`, `PyClone-VI` | 🔬 **In Validation** | Models 3D chromatin contact graph spectral Laplacian connectivity (`Fiedler` $\lambda_2$, Von Neumann entropy $S_{\text{VN}}$) and reconstructs somatic clonal lineage trees with 96-channel trinucleotide mutational signature Jensen-Shannon divergence (`JSD`). |
| **`disease-atlas` & `sweep`** (`src/disease_panels.rs`, `src/sweep.rs`) | `OMIM`, `ClinVar`, `GWAS Catalog`, `Nextflow` | 🔬 **Active Continuous Calibration** | Maintains the 32-panel (630-gene) human disease and pharmacogenomic pleiotropy catalog and executes the continuous cross-kingdom (`Prokaryota` + `Eukaryota`) `0🧪`–`5🧪` ranking engine calibration sweep. |

---

## How It Works

1. **Polite NCBI Sequence Cache (`src/fetch.rs`)**:
   - Downloads verified RefSeq and GenBank sequences from NCBI Entrez under a strict rate limit (`<= 2.5 requests/sec`) and saves them in a local SHA-256 disk cache (`./cache/`).
   - Once cached, subsequent runs execute offline in milliseconds with zero network traffic.
2. **Native-First Execution with Optional External Tool Benchmarks**:
   - By default, all alignments, substitution model tests, maximum-likelihood trees, species-tree coalescent analyses, codon selection models, orthology searches, whole-genome synteny chains, gene predictions, profile HMM searches, non-B DNA scans, 1D/2D/3D nucleic acid structure predictions, protein pocket docking runs, and rGFA pangenome variation graph decompositions run in native Rust.
   - Pass `--aligner mafft`, `--modeltest jmodeltest`, `--tree-builder iqtree`, `--coalescent astral`, or `--selection paml` (or run any `*-compare` subcommand) to benchmark native outputs head-to-head against external reference binaries.
3. **Structured JSON Telemetry & Calibrated Signal Scoring**:
   - Every command outputs a reproducible `.json` report containing sequence accession metadata, focal alignments, phylogenetic trees, 3D structural superposition metrics, and statistical effect sizes used in the [CRV Science Lab Discoveries Journal](https://crvscience.com/labs#discoveries).

---

## Repository Structure

```text
crvscience-labs/
├── Cargo.toml                   # Rust crate manifest & dependencies
├── src/
│   ├── main.rs                  # CLI entry point & command dispatcher
│   ├── lib.rs                   # Library module exports
│   ├── fetch.rs                 # Polite NCBI Entrez fetcher & SHA-256 local cache
│   ├── phylo.rs                 # Sequence alignment, K2P/TN93 distances & NJ trees
│   ├── msa_calibration.rs       # Native aligner vs. MAFFT v7.526 benchmark harness
│   ├── modeltest.rs             # 24-model Felsenstein pruning vs. jModelTest v2.1.9
│   ├── tree_builder.rs          # 4-stage BIONJ+NNI ML tree builder vs. IQ-TREE v2.4.0
│   ├── coalescent.rs            # Quartet multispecies coalescent & ABBA-BABA vs. ASTRAL v5.7.8
│   ├── codon_selection.rs       # 61/60-state Goldman-Yang codon CTMC & 4D BEB vs. PAML codeml
│   ├── orthology.rs             # Smith-Waterman BLOSUM62 & Emms-Kelly RBH vs. DIAMOND v2.2.8
│   ├── synteny.rs               # Double-stranded (w,k)-minimizer & SyRI SV vs. Minimap2 v2.31
│   ├── gene_caller.rs           # 6-frame ORF graph, dicodon Markov & RBS EM vs. Prodigal v2.6.3
│   ├── profile_hmm.rs           # Plan-7 Profile HMM (MSV/Viterbi/Forward-Backward) vs. HMMER v3.4
│   ├── nonb_dna.rs              # Non-B DNA scanner vs. G4Hunter v1.0 & Z-Hunt-II (zhunt3)
│   ├── rna_fold.rs              # Turner-2004 MFE, McCaskill PF & rG4 vs. ViennaRNA v2.7.2
│   ├── na_structure_1d2d3d.rs   # Unified 1D->2D->3D DNA/RNA folding vs. ViennaRNA, RNAPKplex & US-align
│   ├── pocket_dock.rs           # Alpha-sphere pocket druggability & SE(3) docking vs. fpocket & AutoDock Vina
│   ├── pangenome.rs             # rGFA pangenome graph, superbubble SV & walk entropy vs. minigraph & gfatools
│   ├── haplotype.rs             # Haplotype networks & polymorphic site matrices
│   ├── complexity.rs            # Shannon entropy, KL divergence, RSCU & DNA walk scaling
│   ├── genome_target.rs         # Non-B DNA scanners, fine-mapping & GPL-3.0 tool registry
│   ├── pipelines.rs             # HGT resistome, BGC mining, pseudogene & pangenome modules
│   ├── dag_pipelines.rs         # Co-evolution, chromatin & somatic evolution modules
│   ├── avian_systematics.rs     # Tiered avian mitogenome & rRNA supermatrix systematics
│   ├── disease_panels.rs        # 32-panel (630-gene) human disease catalog
│   └── sweep.rs                 # Cross-kingdom calibration sweep & 0★–5★ scoring gates
└── tests/
    └── core_fidelity.rs         # Automated unit & integration test suite
```

*(Note: Local working directories `cache/`, `runs/`, and `tools/` are git-ignored so the repository remains lightweight.)*

---

## License

Released under the **GNU General Public License v3.0 (`GPL-3.0-or-later`)**. See [LICENSE](LICENSE) for full details.
