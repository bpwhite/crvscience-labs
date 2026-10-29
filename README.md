# CRV Science Labs (`bio-suite`)

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
2. **Tier 2 — Composable Task Runners (`bio-suite task`) across 4 Modular Engines**:
   - **`bio-suite phylo` (Molecular Systematics & Phylogenomics)**: Codon-anchored Multiple Sequence Alignment (native Rust + `MAFFT` / `MUSCLE5` / `trimAl`), Kimura 2-Parameter (`K2P`) and Tamura-Nei (`TN93`) evolutionary distance matrices, transition/transversion (`Ti/Tv`) bias estimation, Neighbor-Joining (`NJ`) & UPGMA tree reconstruction with parallel bootstrap replicates (`IQ-TREE 2` / `RAxML-ng` / `RevBayes` / `MrBayes` adapters), and Automatic Barcode Gap Discovery (`ABGD`) MOTU species clustering.
   - **`bio-suite genome` (Structural, Population & Functional Genomics)**: Sliding-window Non-B DNA secondary structure scanners (`G4Hunter` G-quadruplex propensity, `ZHunt` alternating purine-pyrimidine Z-DNA propensity, CpG island `O/E` ratios), annotation-agnostic complexity anomaly peak detection, GFF3/RefSeq coordinate intersection (`intersect_annotations`), and GPL-3.0 drop-in adapters (`Flye`, `hifiasm`, `Verkko`, `Minigraph-Cactus`, `vg`, `Helixer`, `AUGUSTUS`, `EDTA`, `Minimap2`, `SyRI`, `GENESPACE`, `Cooler`, `MACS3`).
   - **`bio-suite target` (Genetic Target Identification & Drug Discovery)**: Bayesian fine-mapping Wakefield Approximate Bayes Factors (`ABF`), posterior inclusion probabilities (`PIP`), colocalization (`SuSiE`, `coloc`), Mendelian Randomization (`IVW`), splicing & CRISPR perturbation (`MAGeCK`, `BAGEL2`, `rMATS-turbo`, `Whippet`), and 100% open-source structural biology & docking (`Boltz-1`, `OpenFold`, `ESMFold`, `fpocket`, `P2Rank`, `AutoDock Vina`, `GNINA`, `OpenMM`).
   - **`bio-suite complexity` (Complexity Theory & Information-Theoretic Genomics)**: Multi-threaded sliding-window $k$-mer Shannon block entropy ($H_1 \dots H_6$), codon positional entropy ($H_{\text{pos1}}, H_{\text{pos2}}, H_{\text{pos3}}$), Jensen-Shannon divergence (`JSD`), Lempel-Ziv 1976 (`LZ76`) algorithmic (Kolmogorov) compressibility, Chaos Game Representation (`CGR`), and Purine-Pyrimidine DNA walks with Rescaled-Range Hurst exponent ($H$) and Detrended Fluctuation Analysis ($\alpha_{\text{DFA}}$).
3. **Tier 3 — Deterministic Worker-Agents & Emergent Pipeline Composer (`bio-suite agent`)**:
   - Seamlessly chains Tier-1 and Tier-2 tasks via strongly typed Rust structs and automatically cross-correlates outputs across engines to detect **emergent biological patterns** worthy of publication in the [CRV Science Lab Discoveries Journal](https://crvscience.com/labs#discoveries).

---

## Building from Source

```bash
git clone https://github.com/bpwhite/crvscience-labs.git
cd crvscience-labs
cargo build --release
```

The standalone native binary is produced at `target/release/bio-suite` (`target/release/bio-suite.exe` on Windows).

---

## Flagship Worker-Agent Pipelines

### 1. Carnivora Mitochondrial COI DNA Barcoding & Complexity Cluster Pipeline
Downloads (or loads from local cache) 14 NCBI mitochondrial Cytochrome c Oxidase Subunit I (`COX1` / `COI`) sequences across 4 Carnivoran families (`Ursidae`, `Felidae`, `Canidae`, `Mustelidae`), performs codon-preserving alignment, computes the $14 \times 14$ K2P/TN93 distance matrix, infers the Neighbor-Joining phylogeny with 500 parallel bootstrap replicates, runs ABGD MOTU barcoding gap clustering, and overlays Tier-4 codon Shannon entropy, CpG depletion, and Lempel-Ziv `LZ76` complexity:

```bash
./target/release/bio-suite agent coi-barcode-cluster --out runs/coi-barcode-report.json
```

### 2. Human Chromosome 1 Section A (`1p36.33`) Unbiased & Coding Structural Complexity Scan
Downloads/caches Human GRCh38 `NC_000001.11:1000001-1060000` (`60 kb` subtelomeric locus), profiles coding sequences (`HES4`, `ISG15`, `AGRN`), runs an annotation-agnostic sliding-window structural complexity scan, and intersects anomalous complexity peaks (`Region-X01`, `Region-X02`, ...) with gene and repeat annotations:

```bash
./target/release/bio-suite agent chr1-unbiased-scan --out runs/chr1-complexity-report.json
```

### 3. Inspect GPL-3.0 Verified Tool & Engine Registry
```bash
./target/release/bio-suite registry
```

---

## License

Released under the **GNU General Public License v3.0 (`GPL-3.0-or-later`)**. See [LICENSE](LICENSE) for full details.
