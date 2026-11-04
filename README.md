# CRV Science Labs (`bio-suite` & `RuNIE` Pipeline Inference Engine)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language: Rust](https://img.shields.io/badge/Language-Rust%202021-orange.svg)](https://www.rust-lang.org/)
[![Engine: RuNIE v1.23.0 (candle-core + rayon)](https://img.shields.io/badge/Inference-RuNIE%20v1.23.0%20GA--MCMC%20%2B%20Candle-8a2be2.svg)](docs/RUNIE_ARCHITECTURE_AND_STUDIES.md)
[![Portal: crvscience.com/labs](https://img.shields.io/badge/Portal-crvscience.com%2Flabs-00bcd4.svg)](https://crvscience.com/labs)

**`bio-suite`** is the free, open-source (**`GNU GPL-3.0`**) Rust command-line toolkit powering [CRV Science Labs](https://crvscience.com/labs). It provides fast, zero-dependency native Rust pipelines for comparative genomics, molecular systematics, 1D/2D/3D nucleic acid & protein structural biophysics, pangenome variation graphs, antimicrobial resistance (`AMR`) & biosynthetic gene cluster (`BGC`) mining, statistical genetics (`coloc.abf`, `MR-Egger`), and multi-stage Bayesian pipeline synthesis via **`RuNIE` (`Rust Native Inference Engine`)**.

Every module runs in native multi-threaded Rust out of the box while offering head-to-head calibration harnesses (`*-compare`) against established public bioinformatics tools (`MAFFT`, `jModelTest`, `IQ-TREE 2`, `ASTRAL`, `PAML`, `DIAMOND`, `Minimap2`, `Prodigal`, `HMMER3`, `ViennaRNA`, `US-align`, `fpocket`, `AutoDock Vina`, `minigraph`, `AMRFinderPlus`, `GECCO`, `MACSE`, `Biopython`, and `SciPy`).

---

<!-- AUTO-FRONTPAGE-STATUS:START -->
| Subsystem | Live Count & Version | Summary & Documentation |
| :--- | :--- | :--- |
| **Native Rust Engines (`bio-suite`)** | **`34` Validated Tools** | Head-to-head calibrated against `MAFFT`, `IQ-TREE`, `PAML`, `ASTRAL`, `HMMER`, `ViennaRNA`, `AMRFinderPlus`, `GECCO`, `MACSE`, `SciPy`, and `Biopython` $\rightarrow$ [**Tool Benchmark Matrix**](docs/TOOL_BENCHMARK_MATRIX.md) |
| **`RuNIE` Bayesian Pipeline Synthesizer** | **`v1.26.0` (`17` Studies / `32` Interventions)** | Parallel-Tempering GA-MCMC DAG synthesizer (`candle-core` + `rayon`) across 17 canonical multi-tier studies $\rightarrow$ [**RuNIE Architecture & Studies**](docs/RUNIE_ARCHITECTURE_AND_STUDIES.md) & [**Changelog**](RUNIE_CHANGELOG.md) |
| **Published Lab Discoveries** | **`104` Reports** (`79` T1 · `9` T2 · `16` T3) | Peer-ready computational reports with TSV Copy-for-Excel boxes, SVG figures & 3D PDB viewers $\rightarrow$ [**Discoveries & Critical Reports**](docs/DISCOVERIES_AND_CRITICAL_REPORTS.md) |
| **Flagged Critical Reports (`🚩`)** | **`22` Foundational Reports** | Blinded `SHA-256` validation ([`#72`](https://crvscience.com/labs/discoveries/popgen-tier2-deep-sweep), [`#73`](https://crvscience.com/labs/discoveries/popgen-tier2-blind-validate), [`#74`](https://crvscience.com/labs/discoveries/popgen-tier3-repeat-kcore), [`#83`](https://crvscience.com/labs/discoveries/heimdall-tier2-blind-validate), [`#84`](https://crvscience.com/labs/discoveries/heimdall-eukaryogenesis-tier3), [`#85`](https://crvscience.com/labs/discoveries/ncldv-histone-chromatin-tier3), [`#86`](https://crvscience.com/labs/discoveries/mito-endosymbiosis-tier3), [`#87`](https://crvscience.com/labs/discoveries/ctenophora-basal-metazoan-tier3), [`#89`](https://crvscience.com/labs/discoveries/cteno-myosin-blind-tier2), [`#90`](https://crvscience.com/labs/discoveries/ctenophora-first-muscle-epithelia-tier3), [`#91`](https://crvscience.com/labs/discoveries/ctenophora-hox-wnt-chromatin-tier3), [`#93`](https://crvscience.com/labs/discoveries/microsporidia-lba-blind-tier2), [`#94`](https://crvscience.com/labs/discoveries/microsporidia-fungi-subdomain-lba-tier3), [`#95`](https://crvscience.com/labs/discoveries/ecdysozoa-nematode-subdomain-lba-tier3), [`#96`](https://crvscience.com/labs/discoveries/placozoa-porifera-epithelial-subdomain), [`#98`](https://crvscience.com/labs/discoveries/plastid-toctic-blind-tier2), [`#99`](https://crvscience.com/labs/discoveries/plastid-endosymbiosis-archaeplastida-tier3), [`#100`](https://crvscience.com/labs/discoveries/vertebrate-adaptive-immunity-rag-vlr), [`#101`](https://crvscience.com/labs/discoveries/titin-obscurin-sarcomere-megaproteins), [`#102`](https://crvscience.com/labs/discoveries/dual-track-partitioning-information), [`#103`](https://crvscience.com/labs/discoveries/bioluminescence-luciferase-photoprotein), [`#104`](https://crvscience.com/labs/discoveries/metazoan-venom-disulfide-scaffold)) $\rightarrow$ [**Critical Reports Index**](docs/DISCOVERIES_AND_CRITICAL_REPORTS.md) |
| **Curated Validation Datasets** | **`35` Datasets (`VAL-SET-001`–`VAL-SET-035`)** | Replicated Development + Untouched Reserve cohorts + `52` logged `0🧪` negative controls $\rightarrow$ [**Validation & Calibration Process**](docs/VALIDATION_AND_CALIBRATION_PROCESS.md) |
<!-- AUTO-FRONTPAGE-STATUS:END -->

---

## Documentation & Deep-Dive Guides

To keep this front page concise, detailed mathematical formulations, validation protocols, and full benchmark tables are organized into dedicated documentation pages:

1. **[`docs/RUNIE_ARCHITECTURE_AND_STUDIES.md`](docs/RUNIE_ARCHITECTURE_AND_STUDIES.md)** — `candle-core` Tensor Markov Policy, Parallel-Tempering `GA-MCMC`, the 4-Generation Step-Up Investigation Loop, Representative Pipeline Log-Likelihood Comparisons (`Best`, `Nearest-Best`, `Medium`, `Rejected`), and all **15 Canonical Multi-Stage Studies**.
2. **[`docs/VALIDATION_AND_CALIBRATION_PROCESS.md`](docs/VALIDATION_AND_CALIBRATION_PROCESS.md)** — The 3-Tier Research Scope System (`Tier 1: Exploratory`, `Tier 2: Validation`, `Tier 3: Scientific`), the `1🧪`–`5🧪` Beaker Signal Gate ($B_{\text{prior}} \rightarrow B_{\text{posterior}}$), Dual-Pass Blinded `SHA-256` Validation ($D_{\text{ctx}}$), Genomic Marker Diagnostic Horizons (`16S`, `23S`, `18S`, `ITS`, `COX1`, `RP16` Dayhoff-6), and the **25 Validation Datasets (`VAL-SET-001`–`VAL-SET-025`)**.
3. **[`docs/TOOL_BENCHMARK_MATRIX.md`](docs/TOOL_BENCHMARK_MATRIX.md)** — Complete **34-Tool Head-to-Head Comparison Table** (Native Rust vs. `MAFFT`, `IQ-TREE 2`, `PAML`, `HMMER3`, `ViennaRNA`, `US-align`, `AutoDock Vina`, `minigraph`, `AMRFinderPlus`, `GECCO`, `MACSE`, `Biopython`, `SciPy`) and the active validation queue.
4. **[`docs/DISCOVERIES_AND_CRITICAL_REPORTS.md`](docs/DISCOVERIES_AND_CRITICAL_REPORTS.md)** — Auto-synced index of all **Flagged Critical Reports (`🚩`)**, recent Lab Discoveries, and dataset provenance mappings.
5. **[`RUNIE_CHANGELOG.md`](RUNIE_CHANGELOG.md)** — Chronological engineering log of every calibration gap identified, native Rust fix implemented, and multi-stage `DAG` integrated (`RN-2026-001` through `RN-2026-030`).

---

## Quick Start & Key CLI Examples

```bash
git clone https://github.com/bpwhite/crvscience-labs.git
cd crvscience-labs
cargo test
cargo build --release
```

### 1. Run the `RuNIE` Bayesian Pipeline Inference Engine (`15` Studies across Tiers 1–3)
```bash
./target/release/bio-suite runie --out runs/runie-pipeline-inference-report.json
```

### 2. Run Step-Up Tier 1 $\rightarrow$ Tier 2 $\rightarrow$ Tier 3 Scientific Syntheses (`1.44–1.52 Mb`)
```bash
# Study #13 (#82, #83 🚩, #84 🚩): Asgard Archaea (Heimdallarchaeia) Eukaryogenesis Synthesis
./target/release/bio-suite asgard-tier3 --out runs/asgard-tier3-eukaryogenesis-synthesis.json

# Study #14 (#85 🚩): Giant DNA Viruses (Nucleocytoviricota) Viral Histones, PolB/RNAP & aaRS Capture
./target/release/bio-suite ncldv-tier3 --out runs/ncldv-tier3-chromatin-synthesis.json

# Study #15 (#86 🚩): Alphaproteobacterial Endosymbiosis, SAM50/MICOS & Mitoribosome Inversion
./target/release/bio-suite mito-tier3 --out runs/mito-tier3-endosymbiosis-synthesis.json
```

### 3. Run Head-to-Head Native vs. External Tool Benchmarks (`Dev + Untouched Reserve`)
```bash
# Tool #34: 5-Stage Adaptive Multi-Biotype Super-BLAST vs. Biopython v1.88 & SciPy v1.18.1
./target/release/bio-suite super-blast --out runs/super-blast-biopython-scipy-calibration.json

# Tool #19: 16-RP Supermatrix & Dayhoff-6 / SR4 Compositional Recoding vs. IQ-TREE v2.4.0
./target/release/bio-suite rp16-compare --out runs/rp16-supermatrix-iqtree-calibration.json

# Tool #16: 3-Frame Codon-to-AA Frameshift Aligner & Branch RELAX vs. MACSE v2.07 & PAML v4.10.10
./target/release/bio-suite relax-compare --out runs/pseudo-loss-macse-relax-calibration.json
```

---

<!-- AUTO-FRONTPAGE-HIGHLIGHTS:START -->
| ID & Flag | Report Title & Live Reader Link | Scope Tier | Signal | Primary Validation Set |
| :--- | :--- | :--- | :--- | :--- |
| **#104** `🚩` | **[Dual-Track Sub-Domain Partitioning of Disulfide-Rich Venom Scaffolds Resolves Convergent Ion-Channel Modulator Recruitment Across Cnidaria, Mollusca, Arthropoda, and Reptilia](https://crvscience.com/labs/discoveries/metazoan-venom-disulfide-scaffold)** | `Tier 3: Scientific` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-037-METAZOAN-VENOM-DISULFIDE-CONVERGENCE-32COHORT` |
| **#103** `🚩` | **[Dual-Track Sub-Domain Partitioning of Luciferases and Ca2+-Regulated Photoproteins Resolves Independent Origins of Bioluminescence Across Coleoptera, Ctenophora, Cnidaria, Dinoflagellata, and Fungi](https://crvscience.com/labs/discoveries/bioluminescence-luciferase-photoprotein)** | `Tier 3: Scientific` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-036-BIOLUMINESCENCE-LUCIFERASE-PHOTOPROTEIN-32COHORT` |
| **#102** `🚩` | **[Dual-Track Sub-Domain Partitioning Recovers 99.4% of Biophysical, Developmental Splicing, and Clinical Signal Lost When Masking Disordered and Repeat Regions](https://crvscience.com/labs/discoveries/dual-track-partitioning-information)** | `Tier 3: Scientific` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-035-DUAL-TRACK-IDR-REPEAT-RECOVERY-30COHORT` |
| **#101** `🚩` | **[Anchor-Segmented Sub-Domain Partitioning of 35,000-Residue Titin, Obscurin, and Nebulin Megaproteins Resolves Bilaterian Sarcomere Evolution and Cardiomyopathy Truncation Hotspots](https://crvscience.com/labs/discoveries/titin-obscurin-sarcomere-megaproteins)** | `Tier 3: Scientific` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-034-SARCOMERIC-MEGAPROTEINS-TITIN-OBSCURIN-28COHORT` |
| **#100** `🚩` | **[Sub-Domain Partitioning of Proto-RAG Recombinases, Lamprey VLR Solenoids, and MHC/TCR Beta-Sandwiches Resolves the Dual Origin of Vertebrate Adaptive Immunity](https://crvscience.com/labs/discoveries/vertebrate-adaptive-immunity-rag-vlr)** | `Tier 3: Scientific` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-033-VERTEBRATE-ADAPTIVE-IMMUNITY-RAG-VLR-24COHORT` |
| **#99** `🚩` | **[Sub-Domain Partitioning of TOC75/Omp85 Channels, FtsZ Division Rings, and RuBisCO Unites Archaeplastida Sister to Gloeomargarita and Contrasts Paulinella Chromatophore Import](https://crvscience.com/labs/discoveries/plastid-endosymbiosis-archaeplastida-tier3)** | `Tier 3: Scientific` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-032-PLASTID-ENDOSYMBIOSIS-TOCTIC-28COHORT` |
| **#98** `🚩` | **[Blinded Sub-Domain Validation of Plastid TOC/TIC Translocons and FtsZ Division Rings Isolates Extremophile Linker Drift, Paulinella ER-Golgi Trafficking, and NUPT Pseudogenes](https://crvscience.com/labs/discoveries/plastid-toctic-blind-tier2)** | `Tier 2: Validation` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-032-PLASTID-ENDOSYMBIOSIS-TOCTIC-28COHORT` |
| **#96** `🚩` | **[Sub-Domain Partitioning of Cadherin, Integrin, and Apical-Basal Polarity Scaffolds Unites Placozoa with Cnidaria and Bilateria Inside Parahoxozoa](https://crvscience.com/labs/discoveries/placozoa-porifera-epithelial-subdomain)** | `Tier 3: Scientific` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-031-PLACOZOA-PORIFERA-SUBDOMAIN-24COHORT` |
| **#95** `🚩` | **[How Sub-Domain Partitioning of Kinesin-1, Ecdysone/DAF-12 Nuclear Receptors, and Molting Chitinases Resolves the Nematode Coelomata vs. Ecdysozoa Long-Branch Attraction Puzzle](https://crvscience.com/labs/discoveries/ecdysozoa-nematode-subdomain-lba-tier3)** | `Tier 3: Scientific` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-030-ECDYSOZOA-NEMATODE-SUBDOMAIN-LBA-24COHORT` |
| **#94** `🚩` | **[How Sub-Domain Partitioning of RNA Polymerase II, Chitin Synthase, and Hsp70 Resolves the Microsporidia Long-Branch Attraction Puzzle at the Base of Fungi](https://crvscience.com/labs/discoveries/microsporidia-fungi-subdomain-lba-tier3)** | `Tier 3: Scientific` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-029-MICROSPORIDIA-FUNGI-SUBDOMAIN-LBA-24COHORT` |
| **#93** `🚩` | **[Blinded Sub-Domain Validation of Microsporidia and Rozellomycota Isolates Repeat-Driven LBA and Truncated Co-Factor Misannotations](https://crvscience.com/labs/discoveries/microsporidia-lba-blind-tier2)** | `Tier 2: Validation` | 🧪🧪🧪🧪 (`4/5`) | `VAL-SET-029-MICROSPORIDIA-FUNGI-SUBDOMAIN-LBA-24COHORT` |
| **#91** `🚩` | **[How the Absence of Hox/ParaHox Clusters, Minimal Wnt Polarity, and CTCF-Negative Chromatin Loops Support the Ctenophora-First Hypothesis](https://crvscience.com/labs/discoveries/ctenophora-hox-wnt-chromatin-tier3)** | `Tier 3: Scientific` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-028-CTENOPHORA-HOX-WNT-CHROMATIN-24COHORT` |
| **#90** `🚩` | **[How Independent Myosin II Duplications, Troponin-Absent Myocytes, and Ancestral Collagen IV Support the Ctenophora-First Hypothesis](https://crvscience.com/labs/discoveries/ctenophora-first-muscle-epithelia-tier3)** | `Tier 3: Scientific` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-027-CTENOPHORA-MYOSIN-EPITHELIA-24COHORT` |
| **#89** `🚩` | **[Blinded Sequence and Structural Validation of Metazoan Myosin II Motor Domains, Z-Disc Scaffolds, and Collagen IV Crosslinks](https://crvscience.com/labs/discoveries/cteno-myosin-blind-tier2)** | `Tier 2: Validation` | 🧪🧪🧪🧪 (`4/5`) | `VAL-SET-027-CTENOPHORA-MYOSIN-EPITHELIA-24COHORT` |
| **#87** `🚩` | **[How Chromosome-Scale Synteny, Compositional Recoding, and Ion Channel Co-Evolution Resolve Ctenophora at the Base of Animal Evolution](https://crvscience.com/labs/discoveries/ctenophora-basal-metazoan-tier3)** | `Tier 3: Scientific` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-026-CTENOPHORA-METAZOA-SYNTENY-24COHORT` |
| **#86** `🚩` | **[How Outer-Membrane Beta-Barrels, Cristae Junctions, and Mitoribosome Remodeling Drove Mitochondrial Origin](https://crvscience.com/labs/discoveries/mito-endosymbiosis-tier3)** | `Tier 3: Scientific` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-025-MITO-ENDOSYMBIOSIS-24COHORT` |
| **#85** `🚩` | **[How Giant Virus Histone Doublets, DNA Polymerases, and Captured tRNA Synthetases Inform Eukaryotic Chromatin Evolution](https://crvscience.com/labs/discoveries/ncldv-histone-chromatin-tier3)** | `Tier 3: Scientific` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-024-NCLDV-CHROMATIN-POLB-24COHORT` |
| **#84** `🚩` | **[How Heimdallarchaeia Cellular Machinery, RNA Polymerase Markers, and Ribosome Evolution Bridge the Archaeal-to-Eukaryotic Transition](https://crvscience.com/labs/discoveries/heimdall-eukaryogenesis-tier3)** | `Tier 3: Scientific` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-023-ASGARD-HEIMDALL-EUKARYOGENESIS-24COHORT` |
| **#83** `🚩` | **[Resolving Compositional Long-Branch Attraction and Metagenomic Bin Chimeras in Asgard Archaeal Phylogenomics](https://crvscience.com/labs/discoveries/heimdall-tier2-blind-validate)** | `Tier 2: Validation` | 🧪🧪🧪🧪 (`4/5`) | `VAL-SET-023-ASGARD-HEIMDALL-EUKARYOGENESIS-24COHORT` |
| **#74** `🚩` | **[How Pure Tandem Repeats Form Haplotype Network Hubs and Dissolve Under Triplet Interruptions](https://crvscience.com/labs/discoveries/popgen-tier3-repeat-kcore)** | `Tier 3: Scientific` | 🧪🧪🧪🧪 (`4/5`) | `VAL-SET-013-POPGEN-HAPLO-MSN-12COHORT` |
| **#73** `🚩` | **[Detecting Database Annotation Errors by Comparing Blinded Sequence Analysis Against Repository Metadata](https://crvscience.com/labs/discoveries/popgen-tier2-blind-validate)** | `Tier 2: Validation` | 🧪🧪🧪🧪 (`4/5`) | `VAL-SET-014-MCDONALD-2011-REPEAT-INDEL-SNP-CLUSTERS` |
| **#72** `🚩` | **[Testing How Repeat-Driven Substitution Rates Decay Across Flanking Genomic Windows and Polymerase Mutants](https://crvscience.com/labs/discoveries/popgen-tier2-deep-sweep)** | `Tier 2: Validation` | 🧪🧪🧪🧪 (`4/5`) | `VAL-SET-014-MCDONALD-2011-REPEAT-INDEL-SNP-CLUSTERS` |
<!-- AUTO-FRONTPAGE-HIGHLIGHTS:END -->

---

## Repository Structure

```text
crvscience-labs/
├── Cargo.toml                            # Rust crate manifest (candle-core, candle-nn, rayon, serde)
├── README.md                             # High-level overview, quick-start examples & repo layout
├── RUNIE_CHANGELOG.md                    # Chronological RuNIE tool-modification & pipeline ledger (RN-2026-001..030)
├── docs/                                 # Detailed technical & validation documentation
│   ├── RUNIE_ARCHITECTURE_AND_STUDIES.md # Candle tensor GA-MCMC math, Step-Up loop & all 15 canonical studies
│   ├── VALIDATION_AND_CALIBRATION_PROCESS.md # 3-Tier system, 0🧪–5🧪 gate, Blinded SHA-256 audit & VAL-SET-001..025
│   ├── TOOL_BENCHMARK_MATRIX.md          # Full 34-tool Native Rust vs. External Counterpart benchmark table
│   └── DISCOVERIES_AND_CRITICAL_REPORTS.md # Auto-synced Critical Reports (🚩) & Lab Discoveries index
├── src/
│   ├── main.rs                           # CLI entry point & subcommand dispatcher
│   ├── lib.rs                            # Library module exports
│   ├── runie.rs                          # RuNIE v1.23.0: Candle tensor Markov policy + GA-MCMC + 15-study optimizer
│   ├── runie_tools.rs                    # Synthesized & calibrated native Rust domain operators
│   ├── graphical_reports.rs              # Native vector SVG graphical reports & 3D PDB structure generator
│   ├── asgard_tier_synthesis.rs          # Study #13 (Tier 1–3): Asgard Archaea (Heimdallarchaeia) eukaryogenesis
│   ├── ncldv_tier_synthesis.rs           # Study #14 (Tier 1–3): Giant DNA Viruses (Nucleocytoviricota) chromatin & PolB
│   ├── mito_tier_synthesis.rs            # Study #15 (Tier 1–3): Alphaproteobacterial endosymbiosis & mitoribosome
│   ├── super_blast.rs                    # 5-stage adaptive multi-biotype Super-BLAST, pre-MSA validator & tip auditor
│   ├── super_blast_calibration.rs        # Native Super-BLAST vs. Biopython v1.88 & SciPy v1.18.1 (VAL-SET-022)
│   ├── kmer_complexity_calibration.rs    # Sliding-window k-mer Shannon (H1..H6), LZ76 NCD & DFA vs. SciPy (VAL-SET-021)
│   ├── splicing_crispr_calibration.rs    # Alternative splicing (PSI, dPSI, JSD) & BAGEL2 CRISPR vs. SciPy (VAL-SET-020)
│   ├── mr_twas_calibration.rs            # Two-Sample Mendelian Randomization (IVW, MR-Egger, WM) vs. SciPy (VAL-SET-019)
│   ├── finemap_coloc_calibration.rs      # Log-space Wakefield ABF credible sets & coloc.abf vs. SciPy (VAL-SET-018)
│   ├── clonal_synth_calibration.rs       # 96-channel mutational JSD, NNLS COSMIC & synthetic lethality (VAL-SET-017)
│   ├── chrom_condensate_calibration.rs   # 3D chromatin Laplacian, TAD insulation & A/B eigenvector (VAL-SET-016)
│   ├── popgen_validation.rs              # Tier 2 Blinded SHA-256 validation & Tier 3 repeat/k-core synthesis (VAL-SET-014..015)
│   ├── popgen_calibration.rs             # Haplotype MSN, Hudson F_ST & Tajima's D vs. SciPy (VAL-SET-013)
│   ├── dca_coevo_calibration.rs          # Candle mfDCA precision matrix inversion & APC vs. SciPy (VAL-SET-012)
│   ├── kabsch_tm_calibration.rs          # Jacobi SVD Kabsch SO(3) superposition & TM-score vs. SciPy (VAL-SET-011)
│   ├── repeat_motif_calibration.rs       # Spidroin/fibroin repeat grammar & ridge regression vs. SciPy (VAL-SET-010)
│   ├── codon_opt_calibration.rs          # Synonymous codon optimizer & Turner-2004 MFE vs. ViennaRNA (VAL-SET-009)
│   ├── rp16_supermatrix_calibration.rs   # Dayhoff-6 / SR4 compositional recoding vs. IQ-TREE v2.4.0 (VAL-SET-008)
│   ├── kaiabc_allostery_calibration.rs   # Jacobi Kirchhoff GNM allostery & RK4 KaiABC oscillator vs. SciPy (VAL-SET-007)
│   ├── mk_morpho_calibration.rs          # Lewis Mk+ASC morphological pruning & PIC vs. IQ-TREE v2.4.0 (VAL-SET-006)
│   ├── pseudo_loss_calibration.rs        # 3-frame codon-to-AA frameshift aligner & RELAX vs. MACSE & PAML (VAL-SET-005)
│   ├── phylo.rs & msa_calibration.rs     # Gotoh affine-gap MSA & K2P/TN93 distances vs. MAFFT v7.526
│   ├── modeltest.rs                      # 24-model Felsenstein pruning vs. jModelTest v2.1.9
│   ├── tree_builder.rs                   # 4-stage BIONJ+NNI ML tree builder vs. IQ-TREE v2.4.0
│   ├── coalescent.rs                     # Quartet multispecies coalescent & ABBA-BABA vs. ASTRAL v5.7.8
│   ├── codon_selection.rs                # 61/60-state Goldman-Yang codon CTMC & BEB vs. PAML v4.10.10
│   ├── orthology.rs                      # Smith-Waterman BLOSUM62 & Emms-Kelly RBH vs. DIAMOND v2.2.8
│   ├── synteny.rs                        # Double-stranded (w,k)-minimizer & SyRI SV vs. Minimap2 v2.31
│   ├── gene_caller.rs                    # 6-frame ORF graph, dicodon Markov & RBS EM vs. Prodigal v2.6.3
│   ├── profile_hmm.rs                    # Plan-7 Profile HMM (MSV/Viterbi/Forward-Backward/Null2) vs. HMMER v3.4
│   ├── nonb_dna.rs & complexity.rs       # Non-B DNA scanner vs. G4Hunter v1.0 & Z-Hunt-II
│   ├── rna_fold.rs & na_structure_1d2d3d.rs # 1D->2D->3D DNA/RNA folding vs. ViennaRNA, RNAPKplex & US-align
│   ├── pocket_dock.rs                    # Alpha-sphere druggability & SE(3) docking vs. fpocket & AutoDock Vina
│   ├── pangenome.rs                      # rGFA pangenome graph & superbubble SV vs. minigraph & gfatools
│   ├── amr_hgt.rs & bgc_mine.rs          # Hierarchical AMR/HGT & Pfam BGC caller vs. AMRFinderPlus & GECCO
│   ├── fetch.rs                          # Polite NCBI Entrez fetcher & SHA-256 local cache
│   ├── disease_panels.rs                 # 32-panel (630-gene) human disease & pharmacogenomic catalog
│   └── sweep.rs                          # Cross-kingdom calibration sweep & 0🧪–5🧪 scoring gates
└── tests/
    └── core_fidelity.rs                  # Automated unit & integration test suite
```

*(Note: Local working directories `cache/`, `runs/`, and `tools/` are git-ignored so the repository remains lightweight.)*

---

## License

Released under the **GNU General Public License v3.0 (`GPL-3.0-or-later`)**. See [LICENSE](LICENSE) for full details.
