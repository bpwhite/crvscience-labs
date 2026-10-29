# CRV Science Labs (`bio-suite` v0.2.0)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language: Rust](https://img.shields.io/badge/Language-Rust%202021-orange.svg)](https://www.rust-lang.org/)
[![Portal: crvscience.com/labs](https://img.shields.io/badge/Portal-crvscience.com%2Flabs-00bcd4.svg)](https://crvscience.com/labs)

**`bio-suite`** is the open-source, high-throughput Rust bioinformatics, chemoinformatics, molecular systematics, and information-theoretic genome complexity engine powering [CRV Science Labs](https://crvscience.com/labs).

Designed around a deterministic **Three-Tier Worker-Agent Architecture**, `bio-suite` eliminates fragile scripting glue by combining native `rayon`-parallelized Rust algorithms with drop-in adapters for GPL-3.0-compatible scientific tools.

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

### Part 2: Deep Multi-Modal Pipeline (High Complexity, Multi-Dataset Integration)

6. **Pipeline 6 — Pangenome Structural-Entropy & Allosteric Cryptic Target Engine (`bio-suite pan-cryptic`)**
   - **Stage 1 (Pangenome Graph Topology & Walk-Entropy):** Parses `.gfa` variation graphs (`Minigraph-Cactus` / `PGGB` / `odgi`) to compute **Graph Walk Shannon Entropy ($H_{\text{walk}}$)**, **Cyclomatic Number ($\beta_1$)**, and **Ollivier-Ricci Graph Curvature ($\kappa$)**.
   - **Stage 2 (Haplotype Coalescence & Balancing Selection):** Reconstructs local Ancestral Recombination Graphs (`tskit` / `Relate` / `panacus`), computing Tajima's $D$ and cohort odds ratios for pathogenic structural haplotypes.
   - **Stage 3 (Epigenomic Phase-Switch & Splicing Entropy):** Integrates native ONT/HiFi `Modkit` 5mC CpG methylation entropy with long-read `scIso-Seq` (`IsoQuant`) to compute **Transcript Isoform Shannon Entropy ($H_{\text{iso}}$)** and cryptic exon inclusion ($\Psi$).
   - **Stage 4 (Conformational Ensemble & Cryptic Pocket Discovery):** Combines `Boltz-1` / `BioEmu` / `OpenMM` Markov State Models with `PocketMiner` and **Residue Interaction Network (RIN) Betweenness Centrality ($C_B$)** to discover isoform-selective cryptic allosteric drug pockets.
   ```bash
   ./target/release/bio-suite pan-cryptic --out runs/pan-cryptic-report.json
   ```

---

## Summary Comparison of Pipelines 3–6

| Pipeline | Complexity Tier | Primary Data Inputs | Core Analytical Hook | Discovery Deliverable |
| --- | --- | --- | --- | --- |
| **3. `hgt-scan`** | **Low (Easy Win)** | Draft assemblies (`.fasta`) | KL-divergence of $k$-mer spectra + flanking integrase/tRNA calls | Mobilized AMR cassettes and horizontally transferred metabolic islands |
| **4. `bgc-mine`** | **Low (Easy Win)** | Bacterial/Fungal genomes, MAGs | `antiSMASH` + `MIBiG` dereplication + `ARTS` self-resistance duplication | Novel natural product clusters paired with their exact molecular target |
| **5. `pseudo-loss`** | **Low (Easy Win)** | Multi-species alignments, Biobank VCFs | `TOGA`/`MACSE` frameshift detection + `HyPhy RELAX` ($\omega \rightarrow 1.0$) | Convergent trait-loss markers and safe, dispensable human knockout drug targets |
| **6. `pan-cryptic`** | **High (Deep Data)** | Pangenome `.gfa`, ONT 5mC, scIso-Seq, MD trajectories | Graph walk entropy ($H_{\text{walk}}$) $\rightarrow$ Isoform entropy ($H_{\text{iso}}$) $\rightarrow$ `PocketMiner` | Isoform-specific cryptic allosteric drug pockets driven by complex structural variants |

---

## Core Module Fidelity Tests & GitHub Code Search Indexing Rules

```bash
cargo test
```

Released under the **GNU General Public License v3.0 (`GPL-3.0-or-later`)**. See [LICENSE](LICENSE) for full details.
