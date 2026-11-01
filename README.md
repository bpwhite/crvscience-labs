# CRV Science Labs (`bio-suite` v0.3.0)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language: Rust](https://img.shields.io/badge/Language-Rust%202021-orange.svg)](https://www.rust-lang.org/)
[![Portal: crvscience.com/labs](https://img.shields.io/badge/Portal-crvscience.com%2Flabs-00bcd4.svg)](https://crvscience.com/labs)

**`bio-suite`** is the open-source Rust command-line toolkit powering [CRV Science Labs](https://crvscience.com/labs). It provides fast, self-contained pipelines for molecular systematics, comparative genomics, nucleic acid 1D/2D/3D structural prediction, sequence complexity analysis, and human disease gene mapping.

Every pipeline runs in native multi-threaded Rust out of the box—requiring zero external software—while offering optional command-line flags and benchmark harnesses to compare results directly against standard public bioinformatics tools (`MAFFT`, `jModelTest`, `IQ-TREE 2`, `ASTRAL`, `PAML codeml`, `DIAMOND`, `Minimap2`, `Prodigal`, `HMMER3`, `G4Hunter`, `Z-Hunt-II`, `ViennaRNA`, `RNAPKplex`, and `US-align`).

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

## Native Tools, Public Counterparts & Functions

### 1. Head-to-Head Calibrated Native Engines vs. Public Counterparts

| Native Tool (`bio-suite` Command & Module) | Public Counterpart(s) | Function of the Tool |
| :--- | :--- | :--- |
| **`align-compare`** (`src/phylo.rs`, `src/msa_calibration.rs`) | `MAFFT v7.526` (`--aligner mafft`), `MUSCLE5`, `trimAl` | Aligns protein-coding and ribosomal RNA sequences using codon-anchored triplet alignment and 9-mer anchor-segmented Gotoh 3-state affine-gap dynamic programming. |
| **`modeltest-compare`** (`src/modeltest.rs`) | `jModelTest v2.1.9` / `PhyML` (`--modeltest jmodeltest`) | Evaluates 24 nucleotide substitution models (`JC`, `F81`, `K80`, `HKY`, `SYM`, `GTR` × `[Uniform, +I, +G4, +I+G4]`) via column-compressed Felsenstein pruning, Padé `[6/6]` matrix exponentials, and `AICc`/`BIC` ranking. |
| **`tree-compare`** (`src/tree_builder.rs`) | `IQ-TREE v2.4.0` (`--tree-builder iqtree`), `RAxML-ng` | Infers maximum-likelihood phylogenetic trees using a 4-stage Multi-Distance BIONJ starter pool, exact Yang discrete $\Gamma_4$ site rates, inside-outside pulley branch optimization, and quartet-polished NNI search. |
| **`coalescent-compare`** (`src/coalescent.rs`) | `ASTRAL v5.7.8` (`--coalescent astral`), `Dsuite` | Infers species trees under incomplete lineage sorting (ILS) from multi-locus gene trees via exact quartet tripartition dynamic programming, Bayesian local posterior probabilities (`LPP`), coalescent branch lengths, and `ABBA-BABA` introgression tests (`Patterson's D`). |
| **`selection-compare`** (`src/codon_selection.rs`) | `PAML v4.10.10` (`codeml`, `--selection paml`), `HyPhy` | Evaluates 61-state standard and 60-state vertebrate mitochondrial Goldman-Yang codon CTMC models (`M0`, `M1a`, `M2a`, Branch 2-Ratio, `RELAX k`) and 10,000-cell 4D grid Bayes Empirical Bayes (`BEB`) positive selection site posteriors ($d_N/d_S = \omega$). |
| **`orthology-compare`** (`src/orthology.rs`) | `DIAMOND v2.2.8` (`blastp --swipe`), `OrthoFinder`, `BUSCO` | Identifies single-copy orthologs, co-ortholog in-paralogs, and out-paralogs across proteomes using exact Smith-Waterman `BLOSUM62` local alignment, Karlin-Altschul $E$-values, and Emms-Kelly length-normalized reciprocal best hits (`RBH`). |
| **`synteny-compare`** (`src/synteny.rs`) | `Minimap2 v2.31` (`-x asm5/asm10/asm20`), `MUMmer4`, `SyRI` | Maps whole-chromosome collinearity, inversions, translocations, and genomic island insertion breakpoints using double-stranded $(w,k)$-minimizer indexing, concave-gap dynamic programming chaining, and split-chain structural variant calling. |
| **`gene-compare`** (`src/gene_caller.rs`) | `Prodigal v2.6.3`, `AUGUSTUS`, `Helixer` | Performs *ab initio* prokaryotic and eukaryotic gene prediction via a 6-frame open reading frame (ORF) graph, dicodon hexamer Markov coding potentials, 28-bin Shine-Dalgarno RBS expectation-maximization training, and operon overlap dynamic programming. |
| **`hmm-compare`** (`src/profile_hmm.rs`) | `HMMER v3.4` (`hmmbuild` / `hmmsearch --max`), `Infernal` | Builds Plan-7 profile hidden Markov models and searches sequence databases for conserved protein/RNA domains using multi-segment ungapped `MSV`, Viterbi, Forward-Backward ensemble log-odds, posterior optimal accuracy alignment, and `Null2` composition bias correction. |
| **`nonb-dna-compare`** (`src/nonb_dna.rs`, `src/complexity.rs`) | `G4Hunter v1.0`, `Z-Hunt-II` (`zhunt3`) | Scans chromosomes for non-canonical DNA conformations, computing run-length quantized G-quadruplex/i-motif propensity, exact $4\times 16$ *anti*/*syn* Boltzmann partition-function B-to-Z DNA superhelical transitions ($\Delta\text{Lk}_{1/2}$, Z-score), and CpG islands. |
| **`rna-fold-compare`** (`src/rna_fold.rs`) | `ViennaRNA v2.7.2` (`RNAfold -p -g`) | Predicts RNA secondary structure via $O(N^3)$ Turner 2004 nearest-neighbor Zuker-Stiegler minimum free energy (`MFE`) dynamic programming, McCaskill inside-outside partition functions ($P_{ij}$), positional pairing entropy $S(i)$, RNA G-quadruplexes (`rG4`), and Compensatory Base Changes (`CBC`). |
| **`na-1d2d3d-compare`** (`src/na_structure_1d2d3d.rs`) | `US-align v20261019`, `ViennaRNA v2.7.2` (`RNAPKplex`, `dna_mathews2004.par`) | Unifies `1D → 2D → 3D` nucleic acid structural prediction for both DNA (`B-DNA`, `Z-DNA`, parallel/antiparallel `G4`) and RNA (`A-RNA`, coaxial `tRNA` L-hinges, H-type pseudoknots), generating 5-atom coarse-grained `.pdb` coordinates validated via Kabsch SVD `RMSD`, `TM-score`, and `INF`. |

---

### 2. Core Systematics, Comparative Genomics & Disease Atlas Pipelines

| Native Tool (`bio-suite` Command & Module) | Public Counterpart(s) | Function of the Tool |
| :--- | :--- | :--- |
| **`coi-barcode-cluster`** (`src/phylo.rs`, `src/agent.rs`) | `MEGA11`, `BOLD`, `ABGD`, `iqtree2` | Computes `K2P`/`TN93` evolutionary distance matrices, Neighbor-Joining/UPGMA trees with parallel bootstrap support, and intraspecific vs. interspecific DNA barcoding gaps across `COX1` and rRNA panels. |
| **`haplo-cluster`** (`src/haplotype.rs`) | `PopART`, `DnaSP v6`, `Network` | Collapses population sequences into unique haplotypes, computes haplotype/nucleotide diversity ($H_d$, $\pi$), builds polymorphic site dot-identity matrices, and constructs Kruskal Minimum Spanning Networks (`MSN`). |
| **`chr1-unbiased-scan`** (`src/complexity.rs`, `src/genome_target.rs`) | `Jellyfish`, `KMC3`, `GeCo3`, `MACS3` | Scans chromosomal windows for $k$-mer Shannon entropy ($H_1\dots H_6$), Lempel-Ziv (`LZ76`) algorithmic complexity, and Purine-Pyrimidine DNA walk Detrended Fluctuation Analysis (`DFA` $\alpha$) intersected with RefSeq annotations. |
| **`hgt-scan`** (`src/pipelines.rs`) | `AMRFinderPlus`, `RGI` (`CARD`), `MOB-suite`, `tRNAscan-SE 2.0`, `ISEScan` | Detects horizontally transferred genomic islands and mobile antimicrobial resistance (`AMR`) plasmids via 256-dim tetranucleotide Kullback-Leibler divergence ($D_{KL}$), 64-codon `RSCU` shift, and flanking `tRNA`/integrase attachment sites. |
| **`bgc-mine`** (`src/pipelines.rs`) | `antiSMASH 7`, `GECCO`, `BiG-SCAPE`, `ARTS`, `Macrel` | Mines microbial genomes for biosynthetic gene clusters (`NRPS`, `PKS`, `RiPP`), dereplicates against `MIBiG` reference clusters via domain-Jaccard similarity, detects duplicated self-resistance target paralogs, and profiles amphipathic antimicrobial peptides (`Eisenberg` $\mu_H$). |
| **`pseudo-loss`** (`src/pipelines.rs`) | `TOGA`, `MACSE`, `HyPhy RELAX`, `gnomAD pLoF` | Identifies shared frameshifts, premature stop codons, splice-site disruptions, and relaxed purifying selection ($\omega \to 1.0$, `RELAX` $k$) across syntenic multi-species coding alignments. |
| **`pan-cryptic`** (`src/pipelines.rs`) | `odgi`, `vg`, `panacus`, `IsoQuant`, `fpocket`, `PocketMiner` | Quantifies pangenome variation graph walk Shannon entropy ($H_{\text{walk}}$), Ollivier-Ricci graph curvature, and alternative transcript isoform entropy ($H_{\text{iso}}$) linked to structural pocket geometry. |
| **`co-evo`** (`src/dag_pipelines.rs`) | `CCMpred`, `plmDCA`, `EVcouplings`, `Infernal` | Computes pairwise mutual information (`MI`), Average Product Correction (`APC`), and pseudolikelihood Potts direct-coupling parameters ($\|J_{ij}\|$) to detect compensatory epistatic networks in proteins and ribosomal RNA (`23S PTC`). |
| **`chrom-condensate`** (`src/dag_pipelines.rs`) | `igraph`, `Cooler`, `HiCExplorer`, `ViennaRNA` | Models 3D chromatin contact graph spectral Laplacian entropy and subtelomeric/ribosomal DNA G-quadruplex phase-separation condensate scaffolding. |
| **`clonal-synth`** (`src/dag_pipelines.rs`) | `SigProfiler`, `PyClone-VI`, `MAGeCK`, `BAGEL2` | Reconstructs somatic clonal lineage trees, evaluates trinucleotide mutational signature Jensen-Shannon divergence (`JSD`), and screens synthetic-lethal paralog co-dependencies. |
| **`disease-atlas`** (`src/disease_panels.rs`) | `OMIM`, `ClinVar`, `GWAS Catalog`, `SuSiE`, `coloc` | Builds the 32-panel (630-gene) human disease and pharmacogenomic pathway catalog with Wakefield Approximate Bayes Factor (`ABF`) fine-mapping credible sets and pleiotropic hub rankings. |
| **`dag-chain` & `sweep`** (`src/dag_pipelines.rs`, `src/sweep.rs`) | `Nextflow`, `Snakemake` | Chains multi-stage genomic, phylogenetic, and structural modules and executes the symmetric cross-kingdom (`Prokaryota` + `Eukaryota`) calibration sweep across `0🧪`–`5🧪` signal gates. |

---

## How It Works

1. **Polite NCBI Sequence Cache (`src/fetch.rs`)**:
   - Downloads verified RefSeq and GenBank sequences from NCBI Entrez under a strict rate limit (`<= 2.5 requests/sec`) and saves them in a local SHA-256 disk cache (`./cache/`).
   - Once cached, subsequent runs execute offline in milliseconds with zero network traffic.
2. **Native-First Execution with Optional External Tool Benchmarks**:
   - By default, all alignments, substitution model tests, maximum-likelihood trees, species-tree coalescent analyses, codon selection models, profile HMM searches, and 1D/2D/3D nucleic acid structure predictions run in native Rust.
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
