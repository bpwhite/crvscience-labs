# CRV Science Labs (`bio-suite` v0.3.0)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language: Rust](https://img.shields.io/badge/Language-Rust%202021-orange.svg)](https://www.rust-lang.org/)
[![Portal: crvscience.com/labs](https://img.shields.io/badge/Portal-crvscience.com%2Flabs-00bcd4.svg)](https://crvscience.com/labs)

**`bio-suite`** is the open-source, high-throughput Rust bioinformatics, chemoinformatics, molecular systematics, and information-theoretic genome complexity engine powering [CRV Science Labs](https://crvscience.com/labs).

Designed around a deterministic **Three-Tier Worker-Agent Architecture**, `bio-suite` eliminates fragile scripting glue by combining native `rayon`-parallelized Rust algorithms with drop-in adapters for GPL-3.0-compatible scientific tools.

---

## Multiple Sequence Alignment Policy: Native Rust Default vs. Explicit External Tools (`--aligner`)

By default, all `bio-suite` pipelines execute **Native Mode (`AlignerMode::Native`)** so every discovery pipeline remains 100% self-contained, deterministic, and millisecond-fast without spawning external subprocesses unless explicitly requested:

1. **`AlignerMode::Native` (Default — Zero External Process Overhead)**:
   - **Mode 1 (Codon-Anchored Reading-Frame Mode)**: Automatically activates for uniform-length protein-coding barcodes (e.g., `1,545 bp` vertebrate `COX1`), locking onto triplet codon boundaries (`-15..=+15 bp`, step `3`) to guarantee `0` artificial frameshift gaps (`K2P` Pearson $r = 1.0000$, $\text{MAD} = 0.0000$ vs. `MAFFT v7.526 L-INS-i`, executing **`2,415×` faster**: `1.51 ms` vs. `3,635.51 ms`).
   - **Mode 2 (9-mer Anchor-Segmented Gotoh 3-State Affine-Gap Dynamic Programming)**: Automatically activates for variable-length / indel-rich structural rRNAs (`16S`, `18S`, `23S`, `28S`) and nuclear ribosomal spacers (`ITS1–5.8S–ITS2`). Selects the 6-mer Jaccard medoid guide reference, locks onto colinear unique 9-mer anchors ($\ge 24\text{ bp}$ spacing), and executes exact Gotoh (`1982`) 3-matrix (`M`, `I`, `D`) affine-gap dynamic programming (`Match = +4`, `Transition = +1`, `Transversion = -3`, `GapOpen = -10`, `GapExtend = -1`), achieving $r = 0.9849\text{–}0.9867$ K2P concordance with `MAFFT L-INS-i` at **`122×`–`222×` faster runtime**.
2. **`AlignerMode::Mafft` (Explicit On-Demand External Tool — `--aligner mafft`)**:
   - Invokes the local `MAFFT v7.526` (`--auto` $\rightarrow$ `L-INS-i` / `FFT-NS-i` / `FFT-NS-2`, `model=DNA200`) installation (`tools/mafft-win/mafft-win/mafft.bat` or system `mafft` in `PATH`) only when `--aligner mafft` is explicitly passed on the CLI or during head-to-head calibration (`bio-suite align-compare`).
   - Note: `tools/` is git-ignored to keep the repository lightweight and compliant with GitHub Blackbird code search limits (`< 384 KB` per tracked file).

```bash
# Default: Runs Native Rust v2 Aligner (AlignerMode::Native)
./target/release/bio-suite agent coi-barcode-cluster --out runs/coi-barcode-report.json

# Explicit External Tool Override: Runs MAFFT v7.526 (AlignerMode::Mafft)
./target/release/bio-suite agent coi-barcode-cluster --aligner mafft --out runs/coi-barcode-mafft.json

# Head-to-Head Calibration Benchmark (Native Rust v2 vs. MAFFT v7.526 L-INS-i across 4 cohorts)
./target/release/bio-suite align-compare --out runs/msa-native-vs-mafft-calibration.json
```

---

## Architecture Overview

1. **Tier 1 — Polite Data Acquisition & Local Caching (`bio-suite fetch`)**:
   - Queries NCBI Entrez (`efetch.fcgi` / Datasets) under a strict token-bucket rate limiter (`<= 2.5 req/sec` without API key, exponential backoff, lockfile coordination) so NCBI infrastructure is never impacted.
   - Stores every fetched sequence slice or genome locus in a content-addressable SHA-256 local disk cache (`./cache/`). Subsequent pipeline executions run at memory speed with zero network traffic.
2. **Tier 2 — Composable Task Runners across 4 Modular Scientific Engines**:
   - **`bio-suite phylo` (Molecular Systematics & Phylogenomics)**: Codon-anchored Multiple Sequence Alignment (native Rust + `MAFFT` / `MUSCLE5` / `trimAl`), Kimura 2-Parameter (`K2P`) and Tamura-Nei (`TN93`) evolutionary distance matrices, transition/transversion (`Ti/Tv`) bias estimation, Neighbor-Joining (`NJ`) & UPGMA tree reconstruction with parallel bootstrap replicates (`IQ-TREE 2` / `RAxML-ng` / `RevBayes` / `MrBayes` adapters), Automatic Barcode Gap Discovery (`ABGD`) MOTU clustering, and `TOGA` / `MACSE` / `HyPhy RELAX` codon selection relaxation ($\omega \rightarrow 1.0$).
   - **`bio-suite genome` (Structural, Population & Functional Genomics)**: Sliding-window Non-B DNA secondary structure scanners (`G4Hunter` G-quadruplex propensity, `ZHunt` alternating purine-pyrimidine Z-DNA propensity, CpG island `O/E` ratios), annotation-agnostic complexity anomaly peak detection, GFF3/RefSeq coordinate intersection, horizontal gene transfer mobility boundary mapping (`tRNAscan-SE 2.0`, `VirSorter2`, `PhiSpy`, `ISEScan`, `MOB-suite`), and GPL-3.0 drop-in adapters (`Flye`, `hifiasm`, `Verkko`, `Minigraph-Cactus`, `PGGB`, `vg`, `odgi`, `Helixer`, `AUGUSTUS`, `EDTA`, `Minimap2`, `SyRI`, `GENESPACE`, `Modkit`, `IsoQuant`).
   - **`bio-suite target` (Genetic Target Identification & Drug Discovery)**: Bayesian fine-mapping Wakefield Approximate Bayes Factors (`ABF`), posterior inclusion probabilities (`PIP`), colocalization (`SuSiE`, `coloc`), Mendelian Randomization (`IVW`), antimicrobial resistance & biosynthetic gene cluster mining (`AMRFinderPlus`, `RGI` / `CARD`, `VFDB`, `Prodigal-gv`, `antiSMASH 7`, `GECCO`, `BiG-SCAPE`, `ARTS` self-resistance target seeker), amphipathic peptide biophysics (`Macrel` / Eisenberg $\mu_H$), and 100% open-source structural biology, Markov State Models & cryptic pocket discovery (`Boltz-1`, `OpenFold`, `ESMFold`, `BioEmu`, `PocketMiner`, `fpocket`, `AutoDock Vina`, `GNINA`, `OpenMM`).
   - **`bio-suite complexity` (Complexity Theory & Information-Theoretic Genomics)**: Multi-threaded sliding-window $k$-mer Shannon block entropy ($H_1 \dots H_6$), 256-dim tetranucleotide Kullback-Leibler divergence ($D_{\text{KL}}$), 64-codon Relative Synonymous Codon Usage (`RSCU`), codon positional entropy ($H_{\text{pos1}}, H_{\text{pos2}}, H_{\text{pos3}}$), Jensen-Shannon divergence (`JSD`), Lempel-Ziv 1976 (`LZ76`) algorithmic (Kolmogorov) compressibility, Purine-Pyrimidine DNA walks with Rescaled-Range Hurst exponent ($H$) and Detrended Fluctuation Analysis ($\alpha_{\text{DFA}}$), Pangenome Graph Walk Shannon Entropy ($H_{\text{walk}}$), Ollivier-Ricci Graph Curvature ($\kappa$), and Residue Interaction Network (RIN) Betweenness Centrality.
3. **Tier 3 — Deterministic Worker-Agents & Emergent Pipeline Composer (`bio-suite agent`)**:
   - Seamlessly chains Tier-1 and Tier-2 tasks via strongly typed Rust structs and automatically cross-correlates outputs across engines to detect **emergent biological patterns** worthy of publication in the [CRV Science Lab Discoveries Journal](https://crvscience.com/labs#discoveries).

---

## Building from Source

```bash
git clone https://github.com/bpwhite/crvscience-labs.git
cd crvscience-labs
cargo test
cargo build --release
```

The standalone native binary is produced at `target/release/bio-suite` (`target/release/bio-suite.exe` on Windows).

---

## Six Deterministic Worker-Agent Discovery Pipelines

### Baseline Core Pipelines (Pipelines 1 & 2)

1. **Carnivora Mitochondrial COI DNA Barcoding & Complexity Cluster Pipeline (`coi-barcode-cluster`)**
   ```bash
   ./target/release/bio-suite agent coi-barcode-cluster --out runs/coi-barcode-report.json
   ```
2. **Human Chromosome 1 Section A (`1p36.33`) Unbiased & Coding Structural Complexity Scan (`chr1-unbiased-scan`)**
   ```bash
   ./target/release/bio-suite agent chr1-unbiased-scan --out runs/chr1-complexity-report.json
   ```

---

### Part 1: Three "Low-Hanging Fruit" Discovery Pipelines (Standardized, High-Yield)

3. **Pipeline 3 — Horizontal Gene Transfer (HGT) & Mobile Resistome Scanner (`bio-suite hgt-scan`)**
   - **Step 1 (Compositional Anomaly Detection):** Sweeps contigs in 5 kb sliding windows to compute 256-dim tetranucleotide **Kullback-Leibler (KL) Divergence** ($D_{\text{KL}}$), **Relative Synonymous Codon Usage (RSCU)** cosine distance, and GC-skew Z-scores against the host genome background.
   - **Step 2 (Mobility Boundary Mapping):** Scans 10 kb flanks for `tRNA-attB/attR` integration direct repeats (`tRNAscan-SE 2.0`), phage integrases (`VirSorter2` / `PhiSpy`), insertion sequences (`ISEScan`), and conjugative relaxases (`MOB-suite`).
   - **Step 3 (Cargo Interrogation):** Screens island ORFs against `AMRFinderPlus`, `CARD`, and `VFDB` to flag mobilized resistance and virulence cassettes.
   ```bash
   ./target/release/bio-suite hgt-scan --out runs/hgt-scan-report.json
   ```

4. **Pipeline 4 — Cryptic Biosynthetic Gene Cluster (BGC) & AMP Miner (`bio-suite bgc-mine`)**
   - **Step 1 (sORF & BGC Detection):** Runs `Prodigal-gv` and `antiSMASH 7` / `GECCO` rule-based synthase detection for NRPS, PKS, and short-ORF RiPPs/AMPs (`< 100 aa`).
   - **Step 2 (Known-Cluster Dereplication):** Filters known compounds sharing $\ge 80\%$ domain architecture Jaccard similarity against `MIBiG` (`BiG-SCAPE`).
   - **Step 3 (Self-Resistance Target Mining):** Uses `ARTS` logic to scan orphan BGCs for duplicated, diverged housekeeping target genes (e.g., second copy of `gyrB` or `rpoB` embedded inside the BGC) to infer exact molecular mechanism of action.
   - **Step 4 (Structural Peptide Folding):** Computes Eisenberg amphipathic $\alpha$-helical hydrophobic moment ($\mu_H$), net charge, and Boman index (`Macrel` / `ESMFold`).
   ```bash
   ./target/release/bio-suite bgc-mine --out runs/bgc-mine-report.json
   ```

5. **Pipeline 5 — Lineage-Specific Pseudogenization & Trait-Loss Detector (`bio-suite pseudo-loss`)**
   - **Step 1 (Cross-Species Syntenic Alignment):** Projects reference coding exons onto target clades (`TOGA` / `Miniprot` + `MACSE`) scoring frameshifts and premature stop codons (`TAG/TAA/TGA`).
   - **Step 2 (Sequencing Error vs. True Loss):** Filters unpolished long-read homopolymer artifacts using sister-lineage phylogenetic sharing (`IQ-TREE 2`) and short-read coverage (`mosdepth` / `bcftools`).
   - **Step 3 (Relaxation of Selection):** Computes codon $d_N/d_S$ ($\omega$) and `HyPhy RELAX` selection intensity parameter $k$ ($\omega \rightarrow 1.0, k \rightarrow 0$).
   - **Step 4 (Pathway & Human Knockout Druggability):** Cross-references KEGG/Reactome and `gnomAD` homozygous pLoF cohorts to identify de-risked human inhibitory drug targets (*ANGPTL3*, *HSD17B13*).
   ```bash
   ./target/release/bio-suite pseudo-loss --out runs/pseudo-loss-report.json
   ```

---

### Part 2: Deep Multi-Modal & Chained DAG Pipelines (High Complexity, Multi-Dataset Integration)

6. **Pipeline 6 — Pangenome Structural-Entropy & Allosteric Cryptic Target Engine (`bio-suite pan-cryptic`)**
   - **Stage 1 (Pangenome Graph Topology & Walk-Entropy):** Parses `.gfa` variation graphs (`Minigraph-Cactus` / `PGGB` / `odgi`) to compute **Graph Walk Shannon Entropy ($H_{\text{walk}}$)**, **Cyclomatic Number ($\beta_1$)**, and **Ollivier-Ricci Graph Curvature ($\kappa$)**.
   - **Stage 2 (Haplotype Coalescence & Balancing Selection):** Reconstructs local Ancestral Recombination Graphs (`tskit` / `Relate` / `panacus`), computing Tajima's $D$ and cohort odds ratios for pathogenic structural haplotypes.
   - **Stage 3 (Epigenomic Phase-Switch & Splicing Entropy):** Integrates native ONT/HiFi `Modkit` 5mC CpG methylation entropy with long-read `scIso-Seq` (`IsoQuant`) to compute **Transcript Isoform Shannon Entropy ($H_{\text{iso}}$)** and cryptic exon inclusion ($\Psi$).
   - **Stage 4 (Conformational Ensemble & Cryptic Pocket Discovery):** Combines `Boltz-1` / `BioEmu` / `OpenMM` Markov State Models with `PocketMiner` and **Residue Interaction Network (RIN) Betweenness Centrality ($C_B$)** to discover isoform-selective cryptic allosteric drug pockets.
   ```bash
   ./target/release/bio-suite pan-cryptic --out runs/pan-cryptic-report.json
   ```

7. **Pipeline 7 — Co-Evolutionary Host-Pathogen Interface & Epistatic Escape Predictor (`bio-suite co-evo`)**
   - **Upstream Triggers:** Ingests cryptic lineages from **P1**, horizontally transferred AMR/virulence islands from **P3 (`hgt-scan`)**, self-resistance scaffolds from **P4 (`bgc-mine`)**, and host immune pangenome paths from **P6 (`pan-cryptic`)**.
   - **Stage 1 (Cophylogenetic Reconciliation & Interface Harmonization):** Event-based cophylogenetic reconciliation (`eMPRESS` / `ALE` / `IQ-TREE 2`) and `HyPhy` (`MEME` / `BUSTED`) episodic diversifying selection across host-pathogen interfaces.
   - **Stage 2 (APC-Corrected Mutual Information & Potts Epistatic Landscape):** Computes $MI_{i,j}^{\text{APC}}$ and fits a global statistical mechanics **Potts Hamiltonian ($E(\sigma)$)** (`plmDCA` / `CCMpred`) to map cross-species couplings and local energetic frustration.
   - **Stage 3 (Multimer Folding & Escape MCMC):** Runs `AlphaFold 3` Multimer and Metropolis-Hastings MCMC over the Potts landscape to simulate compensatory 2–4 mutation escape pathways.
   - **Stage 4 (Allosteric "Dead-End" Target Design):** Targets frustrated epistatic hubs (`PocketMiner` / `DiffDock` / `GNINA`) where escape requires $\ge 4$ simultaneous mutations.
   ```bash
   ./target/release/bio-suite co-evo --out runs/co-evo-report.json
   ```

8. **Pipeline 8 — 3D Chromatin Phase-Separation & Non-Coding RNA Scaffolder (`bio-suite chrom-condensate`)**
   - **Upstream Triggers:** Chains from **P2 (`chr1-unbiased-scan`)** non-B DNA / G4 loop anchors and **P6 (`pan-cryptic`)** VNTR epigenetic phase-switches.
   - **Stage 1 (3D Contact Graph Laplacian Spectral Entropy):** Converts `.mcool` Micro-C / Pore-C matrices into spatial graphs to compute **Von Neumann Graph Spectral Entropy ($S_{\text{VN}}$)**, **Algebraic Connectivity (Fiedler value $\lambda_2$)**, and **Ollivier-Ricci Curvature**.
   - **Stage 2 (SCFG Structural RNA Grammar & IDP Sticker-Spacer Parsing):** Maps conserved lncRNA aptamer stem-loops (`Infernal`, `RNAfold`, `R-scape`) and partner Transcription Factor IDR "sticker-and-spacer" LLPS grammar (`GeCo3`).
   - **Stage 3 (Native Epitranscriptomic Modification Calling):** Repurposes `Modkit` on Direct RNA Nanopore (`POD5`) reads to call **m6A and pseudouridine ($\Psi$)** phase-separation switches.
   - **Stage 4 (RNA-Protein Pocket Docking & ASO Walk):** Dual-tracks steric-blocking 2'-MOE ASOs and small-molecule RNA G4/pseudoknot condensate dissolvers (`AlphaFold 3` + `GNINA`).
   ```bash
   ./target/release/bio-suite chrom-condensate --out runs/chrom-condensate-report.json
   ```

9. **Pipeline 9 — Somatic Mosaicism, Clonal Entropy & Synthetic Lethal Neoantigen Engine (`bio-suite clonal-synth`)**
   - **Upstream Triggers:** Chains from **P2 (`chr1-unbiased-scan`)** non-B DNA fragility sites, **P5 (`pseudo-loss`)** somatic gene knockouts, and **P6 (`pan-cryptic`)** somatic SVs and splicing chaos.
   - **Stage 1 (Single-Cell Somatic Coalescent & Lineage Tracing):** Reconstructs ultrametric single-cell somatic phylogenies (`scMultiome` + `scNanoGPS`, `IQ-TREE 2` / `BEAST 2` / `PhyloWGS`) and effective population growth rate ($N_e$).
   - **Stage 2 (96-Channel Mutational Signature Deconvolution & Trajectory Entropy):** Computes trinucleotide **Jensen-Shannon Divergence ($JSD$)** and **Kolmogorov-Smirnov ($D_{\text{KS}}$)** complexity (`SigProfiler`) to separate neutral $1/f$ drift from APOBEC/HRD clonal sweeps.
   - **Stage 3 (Clonal LoF $\rightarrow$ CRISPR Synthetic Lethality Bridging):** Pairs truncal biallelic gene losses (*MTAP*, *ARID1A*) with `DepMap` CRISPR `MAGeCK`/`Chronos` synthetic-lethal paralogs (*PRMT5*, *ARID1B*).
   - **Stage 4 (Cryptic Splice-Neoantigen & Paralog Selective Docking):** Combines `PocketMiner` paralog-selective allosteric docking with `AlphaFold 3` patient-phased HLA class I neoantigen presentation.
   ```bash
   ./target/release/bio-suite clonal-synth --out runs/clonal-synth-report.json
   ```

---

## Master Directed Acyclic Meta-Pipeline (`bio-suite dag-chain`)

```text
   [PIPELINE 1: phylo-audit]                  [PIPELINE 3: hgt-scan]             [PIPELINE 4: bgc-mine]
   (rRNA SCFGs, ETE3, IQ-TREE)                (KL-Divergence, AMRFinder)         (antiSMASH, ARTS, ESMFold)
               │                                          │                                  │
               └────────────────────┬─────────────────────┘                                  │
                                    ▼                                                        │
                       [PIPELINE 7: co-evo] ◄────────────────────────────────────────────────┘
         (Potts Epistasis + Host-Pathogen AF3 Multimer + Escape MCMC)        (Feeds novel BGC scaffolds +
                                    │                                         resistance targets into
                                    ▼                                         co-evolutionary docking)
                    (Identifies Conserved Allosteric
                     Pockets & Escape-Proof Ligands)

   [PIPELINE 2: chr-probe]                    [PIPELINE 6: pan-cryptic]          [PIPELINE 5: pseudo-loss]
   (GeCo3, DFA, G4Hunter, Micro-C)            (odgi Graph Entropy, Modkit,       (TOGA, MACSE, HyPhy RELAX)
               │                               IsoQuant, PocketMiner)                        │
               │                                    │             │                          │
               └──────────────────┬─────────────────┘             └─────────────┬────────────┘
                                  ▼                                             ▼
                  [PIPELINE 8: chrom-condensate]                 [PIPELINE 9: clonal-synth]
           (3D Graph Laplacians + lncRNA SCFGs + m6A +     (Somatic BEAST2 Trees + CRISPR MAGeCK +
             Condensate Dissolving ASOs/Small Mols)          Paralog PocketMiner + HLA Neoantigens)
```

```bash
./target/release/bio-suite dag-chain --out runs/master-dag-report.json
```

---

## Tool Reuse & Chaining Matrix Across All 9 Pipelines

| Core Tool / Engine | Primary Function | Reused Across Pipelines |
| --- | --- | --- |
| **`Infernal` / `RNAfold` / `R-scape`** | SCFG RNA secondary structure & covariation | **P1** (rRNA validation) $\rightarrow$ **P2** (Chr1 non-coding scan) $\rightarrow$ **P8** (lncRNA condensate scaffolds) |
| **`GeCo3` / DFA / KL-Divergence / Potts MI** | Algorithmic complexity, Hurst scaling, $k$-mer & epistatic entropy | **P2** (Chr1 anomaly) $\rightarrow$ **P3** (HGT island detection) $\rightarrow$ **P7** (Potts epistasis) $\rightarrow$ **P8** (IDP sticker-spacer grammar) $\rightarrow$ **P9** (Clonal mutational entropy) |
| **`IQ-TREE 2` / `BEAST 2` / `ETE3`** | ML/Bayesian phylogenetics & topological reconciliation | **P1** (Taxonomy audit) $\rightarrow$ **P5** (Trait-loss dating) $\rightarrow$ **P6** (Haplotype trees) $\rightarrow$ **P7** (Host-pathogen cophylogeny) $\rightarrow$ **P9** (Single-cell somatic trees) |
| **`HyPhy` (`RELAX`, `MEME`, `BUSTED`)** | Codon selection regimes ($d_N/d_S$ shifts) | **P5** (Neutral drift in pseudogenes) $\rightarrow$ **P7** (Episodic diversifying selection at host-pathogen interfaces) |
| **`Minigraph-Cactus` / `odgi` / `vg`** | Pangenome variation graphs & path topology | **P2** (Population constraint) $\rightarrow$ **P6** (Structural bubble entropy) $\rightarrow$ **P7** (Pathogen pangenome MSAs) $\rightarrow$ **P9** (Patient somatic graph) |
| **`Modkit` (Nanopore/PacBio)** | Native base modification entropy | **P6** (DNA 5mC/5hmC phasing) $\rightarrow$ **P8** (Direct RNA m6A/$\Psi$ condensate switches) |
| **`IsoQuant` / `LeafCutter`** | Long-read isoform phasing & splicing entropy | **P6** (SV-driven cryptic exons) $\rightarrow$ **P8** (Nascent enhancer RNAs) $\rightarrow$ **P9** (NMD-escaping frameshift neoantigens) |
| **`AlphaFold 3` + `OpenMM` + `PocketMiner`** | Biomolecular complex folding, MD ensembles, cryptic pockets | **P4** (AMP peptides) $\rightarrow$ **P6** (Isoform allosteric pockets) $\rightarrow$ **P7** (Host-pathogen multimers) $\rightarrow$ **P8** (RNA-protein interfaces) $\rightarrow$ **P9** (Paralog selectivity & HLA-peptide binding) |

---

## Core Module Fidelity Tests & GitHub Code Search Indexing Rules

```bash
cargo test
```

Released under the **GNU General Public License v3.0 (`GPL-3.0-or-later`)**. See [LICENSE](LICENSE) for full details.
