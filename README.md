# CRV Science Labs (`bio-suite` v0.3.0)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language: Rust](https://img.shields.io/badge/Language-Rust%202021-orange.svg)](https://www.rust-lang.org/)
[![Portal: crvscience.com/labs](https://img.shields.io/badge/Portal-crvscience.com%2Flabs-00bcd4.svg)](https://crvscience.com/labs)

**`bio-suite`** is the open-source Rust command-line toolkit powering [CRV Science Labs](https://crvscience.com/labs). It provides fast, self-contained pipelines for molecular systematics, comparative genomics, sequence complexity analysis, and human disease gene mapping.

Every pipeline runs in native multi-threaded Rust out of the box—requiring zero external software—while offering optional command-line flags to benchmark results directly against standard external tools (`MAFFT`, `jModelTest`, `IQ-TREE 2`, and `ASTRAL`).

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

## How It Works

1. **Polite NCBI Sequence Cache (`src/fetch.rs`)**:
   - Downloads verified RefSeq and GenBank sequences from NCBI Entrez under a strict rate limit (`<= 2.5 requests/sec`) and saves them in a local SHA-256 disk cache (`./cache/`).
   - Once cached, subsequent runs execute offline in milliseconds with zero network traffic.
2. **Native-First Execution with Optional External Tool Benchmarks**:
   - By default, all alignments, substitution model tests, maximum-likelihood trees, and species-tree coalescent analyses run in native Rust.
   - Pass `--aligner mafft`, `--modeltest jmodeltest`, `--tree-builder iqtree`, or `--coalescent astral` to compare native outputs side by side against external benchmark binaries.
3. **Structured JSON Telemetry & Calibrated Signal Scoring**:
   - Every command outputs a reproducible `.json` report containing sequence accession metadata, focal alignments, phylogenetic trees, and statistical effect sizes used in the [CRV Science Lab Discoveries Journal](https://crvscience.com/labs#discoveries).

---

## Current Pipelines & Commands

### 1. Molecular Systematics, Barcoding & Haplotype Networks

| Command | What It Does | Example |
| --- | --- | --- |
| `coi-barcode-cluster` | Aligns mitochondrial `COX1` or ribosomal RNA sequences, computes evolutionary distances (`K2P` / `TN93`), builds Neighbor-Joining trees, and measures species barcoding gaps. | `./target/release/bio-suite agent coi-barcode-cluster --out runs/coi-barcode-report.json` |
| `haplo-cluster` | Builds within-species haplotype networks, polymorphic site matrices, and minimum spanning trees from population sequence panels. | `./target/release/bio-suite haplo-cluster --out runs/haplo-cluster-report.json` |

---

### 2. Native Rust vs. External Tool Benchmarks

| Command | Native Engine (Default) | External Benchmark Tool | What It Compares |
| --- | --- | --- | --- |
| `align-compare` | Codon-anchored + Gotoh affine-gap aligner (`--aligner native`) | `MAFFT v7.526` (`--aligner mafft`) | Pairwise distance correlation, gap placement, and runtime across coding and rRNA cohorts. |
| `modeltest-compare` | 24-model Felsenstein pruning engine (`--modeltest native`) | `jModelTest v2.1.9` / `PhyML` (`--modeltest jmodeltest`) | Log-likelihood scores, AICc/BIC model selection, and downstream tree topology. |
| `tree-compare` | 4-stage BIONJ + pulley NNI ML tree builder (`--tree-builder native`) | `IQ-TREE v2.4.0` (`--tree-builder iqtree`) | Maximum-likelihood scores, branch lengths, and clade recovery across 54 vertebrate species. |
| `coalescent-compare` | Quartet multispecies coalescent & `ABBA-BABA` (`--coalescent native`) | `ASTRAL v5.7.8` (`--coalescent astral`) | Species-tree topology, quartet support (`q1, q2, q3`), local posterior probabilities, and introgression across 1,064 gene trees. |

---

### 3. Structural, Microbial & Comparative Genomics

| Command | What It Does | Example |
| --- | --- | --- |
| `chr1-unbiased-scan` | Scans chromosomal windows for non-canonical DNA structures (G-quadruplexes, Z-DNA, CpG islands) and sequence entropy peaks. | `./target/release/bio-suite agent chr1-unbiased-scan --out runs/chr1-complexity-report.json` |
| `hgt-scan` | Detects horizontally transferred genomic islands and mobile antimicrobial resistance cassettes using oligonucleotide divergence and flanking tRNA/integrase markers. | `./target/release/bio-suite hgt-scan --out runs/hgt-scan-report.json` |
| `bgc-mine` | Scans microbial genomes for biosynthetic gene clusters, filters known reference clusters, and flags duplicated self-resistance target genes. | `./target/release/bio-suite bgc-mine --out runs/bgc-mine-report.json` |
| `pseudo-loss` | Identifies shared frameshifts, premature stop codons, and relaxed purifying selection (`dN/dS`) across multi-species coding alignments. | `./target/release/bio-suite pseudo-loss --out runs/pseudo-loss-report.json` |
| `pan-cryptic` | Measures human pangenome variation graph walk entropy, transcript isoform diversity, and structural pocket geometry across phased haplotypes. | `./target/release/bio-suite pan-cryptic --out runs/pan-cryptic-report.json` |
| `co-evo` | Computes mutual information and Potts direct-coupling networks to identify co-evolving residue pairs in proteins and ribosomal RNA. | `./target/release/bio-suite co-evo --out runs/co-evo-report.json` |
| `chrom-condensate` | Analyzes 3D chromatin contact graph connectivity alongside conserved non-coding RNA secondary structures and modification sites. | `./target/release/bio-suite chrom-condensate --out runs/chrom-condensate-report.json` |
| `clonal-synth` | Reconstructs somatic lineage trees, compares mutational signature divergence, and screens paralog co-dependency candidates. | `./target/release/bio-suite clonal-synth --out runs/clonal-synth-report.json` |

---

### 4. Human Disease Gene Panels & Calibration Sweeps

| Command | What It Does | Example |
| --- | --- | --- |
| `disease-atlas` | Builds the 32-panel (630-gene) human disease and pharmacogenomic pathway catalog, ranking cross-panel pleiotropic hub genes. | `./target/release/bio-suite disease-atlas --out runs/disease-genes-atlas.json` |
| `dag-chain` | Chains outputs across the core genomic, structural, and comparative modules into a combined multi-stage summary report. | `./target/release/bio-suite dag-chain --out runs/master-dag-report.json` |
| `sweep` | Executes the cross-kingdom calibration sweep (including negative controls and multi-locus benchmarks) across prokaryotic and eukaryotic datasets. | `./target/release/bio-suite sweep --out runs/multi-pipeline-sweep-report.json` |

---

## Repository Structure

```text
crvscience-labs/
├── Cargo.toml                 # Rust crate manifest & dependencies
├── src/
│   ├── main.rs                # CLI entry point & command dispatcher
│   ├── fetch.rs               # Polite NCBI Entrez fetcher & SHA-256 local cache
│   ├── phylo.rs               # Sequence alignment, K2P/TN93 distances & NJ trees
│   ├── msa_calibration.rs     # Native aligner vs. MAFFT benchmark harness
│   ├── modeltest.rs           # 24-model Felsenstein pruning vs. jModelTest
│   ├── tree_builder.rs        # 4-stage BIONJ+NNI ML tree builder vs. IQ-TREE 2
│   ├── coalescent.rs          # Quartet multispecies coalescent & ABBA-BABA vs. ASTRAL
│   ├── haplotype.rs           # Haplotype networks & polymorphic site matrices
│   ├── complexity.rs          # Shannon entropy, KL divergence, RSCU & DNA walk scaling
│   ├── genome_target.rs       # Non-B DNA scanners, fine-mapping & tool registry
│   ├── pipelines.rs           # HGT resistome, BGC mining, pseudogene & pangenome modules
│   ├── dag_pipelines.rs       # Co-evolution, chromatin & somatic evolution modules
│   ├── disease_panels.rs      # 32-panel (630-gene) human disease catalog
│   └── sweep.rs               # Cross-kingdom calibration sweep & 0★–5★ scoring gates
└── tests/
    └── core_fidelity.rs       # Automated unit & integration test suite
```

*(Note: Local working directories `cache/`, `runs/`, and `tools/` are git-ignored so the repository remains lightweight.)*

---

## License

Released under the **GNU General Public License v3.0 (`GPL-3.0-or-later`)**. See [LICENSE](LICENSE) for full details.
