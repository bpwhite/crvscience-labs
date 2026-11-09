# CRV Science Labs (`bio-suite` & `RuNIE` Pipeline Inference Engine)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language: Rust](https://img.shields.io/badge/Language-Rust%202021-orange.svg)](https://www.rust-lang.org/)
[![Engine: RuNIE v2.0.1 (candle-core + rayon)](https://img.shields.io/badge/Inference-RuNIE%20v2.0.1%20GA--MCMC%20%2B%20ModelTest-8a2be2.svg)](docs/RUNIE_ARCHITECTURE_AND_STUDIES.md)
[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23969275.svg)](https://doi.org/10.5281/zenodo.23969275)
[![Figshare DOI](https://img.shields.io/badge/Figshare%20DOI-10.6084%2Fm9.figshare.34919469-f472b6.svg)](https://doi.org/10.6084/m9.figshare.34919469)
[![Portal: crvscience.com/labs](https://img.shields.io/badge/Portal-crvscience.com%2Flabs-00bcd4.svg)](https://crvscience.com/labs)

**`bio-suite`** is the free, open-source (**`GNU GPL-3.0`**) Rust command-line toolkit powering [CRV Science Labs](https://crvscience.com/labs). It provides fast, zero-dependency native Rust pipelines for comparative genomics, molecular systematics, 1D/2D/3D nucleic acid & protein structural biophysics, pangenome variation graphs, antimicrobial resistance (`AMR`) & biosynthetic gene cluster (`BGC`) mining, statistical genetics (`coloc.abf`, `MR-Egger`), and multi-stage Bayesian pipeline synthesis via **`RuNIE` (`Rust Native Inference Engine`)**.

Every module runs in native multi-threaded Rust out of the box while offering head-to-head calibration harnesses (`*-compare`) against established public bioinformatics tools (`MAFFT`, `jModelTest`, `IQ-TREE 2`, `ASTRAL`, `PAML`, `BEAST`, `DIAMOND`, `Minimap2`, `Prodigal`, `HMMER3`, `ViennaRNA`, `US-align`, `fpocket`, `AutoDock Vina`, `minigraph`, `AMRFinderPlus`, `GECCO`, `MACSE`, `Biopython`, and `SciPy`).

---

<!-- AUTO-FRONTPAGE-STATUS:START -->
| Subsystem | Live Count & Version | Summary & Documentation |
| :--- | :--- | :--- |
| **Native Rust Engines (`bio-suite`)** | **`36` Validated Tools (`23.19×` Speedup)** | ✅ **Validated**: Head-to-head calibrated against `MAFFT`, `jModelTest`, `IQ-TREE 2`, `PAML`, `ASTRAL`, `BEAST`, `HMMER3`, `ViennaRNA`, `AMRFinderPlus`, `GECCO`, `MACSE`, `SciPy`, and `Biopython` $\rightarrow$ [**Tool Benchmark Matrix**](docs/TOOL_BENCHMARK_MATRIX.md) |
| **`RuNIE` Bayesian Pipeline Synthesizer** | **`v2.0.1` (Stage 4 `ModelTest` + 3-Way Engine)** | ✅ **Validated**: 4-Chain Parallel-Tempering GA-MCMC DAG synthesizer (`candle-core` + `rayon`), Stage 4 `ModelTest` (`BIC`/`AICc`/`RCFV <= 0.050`) gate, 3-Way Native/External/Hybrid execution & `RuNIE Pipeline Flowchart` $\rightarrow$ [**RuNIE Architecture**](docs/RUNIE_ARCHITECTURE_AND_STUDIES.md) & [**Changelog**](RUNIE_CHANGELOG.md) |
| **Empirically Verified Preprints** | **`23` Preprints (`PREPRINT-01`–`PREPRINT-23`, `v2.0.1`)** | ✅ **Validated**: Cached NCBI RefSeq/GenBank sequence syntheses with single-column arXiv PDFs, LaTeX/BibTeX bundles, and Master Concept Zenodo DOIs + Figshare mirrors $\rightarrow$ [**Verified Preprints Portal**](https://crvscience.com/labs#preprints) |
| **Lab Discoveries Directory** | **`134` Total** (`23` Preprints · `35` Tool Cal. · `76` Pipeline Benchmarks) | Full index of empirical preprints, head-to-head tool calibrations, and multi-stage `RuNIE` pipeline benchmarks $\rightarrow$ [**Discoveries & Validation Scope**](docs/DISCOVERIES_AND_CRITICAL_REPORTS.md) |
| **Curated Validation Datasets** | **`62` Datasets (`VAL-SET-001`–`VAL-SET-062`)** | Replicated Development + Untouched Reserve cohorts + `54` logged `0🧪` negative controls $\rightarrow$ [**Validation & Calibration Process**](docs/VALIDATION_AND_CALIBRATION_PROCESS.md) |
<!-- AUTO-FRONTPAGE-STATUS:END -->

---

## Documentation & Deep-Dive Guides

To keep this front page concise, detailed mathematical formulations, validation protocols, and full benchmark tables are organized into dedicated documentation pages:

1. **[`docs/RUNIE_ARCHITECTURE_AND_STUDIES.md`](docs/RUNIE_ARCHITECTURE_AND_STUDIES.md)** — `candle-core` Tensor Markov Policy, Parallel-Tempering `GA-MCMC`, the 4-Generation Step-Up Investigation Loop, Representative Pipeline Log-Likelihood Comparisons (`Best`, `Nearest-Best`, `Medium`, `Rejected`), the **20 Empirically Verified Preprints (`PREPRINT-01`–`PREPRINT-20`)**, and all **15 Canonical Multi-Stage Studies**.
2. **[`docs/VALIDATION_AND_CALIBRATION_PROCESS.md`](docs/VALIDATION_AND_CALIBRATION_PROCESS.md)** — The 3-Tier Research Scope System (`Tier 1: Exploratory`, `Tier 2: Validation`, `Tier 3: Scientific`), the `1🧪`–`5🧪` Beaker Signal Gate ($B_{\text{prior}} \rightarrow B_{\text{posterior}}$), Dual-Pass Blinded `SHA-256` Validation ($D_{\text{ctx}}$), Genomic Marker Diagnostic Horizons (`16S`, `23S`, `18S`, `ITS`, `COX1`, `RP16` Dayhoff-6), and the **59 Validation Datasets (`VAL-SET-001`–`VAL-SET-059`)**.
3. **[`docs/TOOL_BENCHMARK_MATRIX.md`](docs/TOOL_BENCHMARK_MATRIX.md)** — Complete **36-Tool Head-to-Head Comparison Table** (Native Rust vs. `MAFFT`, `IQ-TREE 2`, `PAML`, `BEAST`, `HMMER3`, `ViennaRNA`, `US-align`, `AutoDock Vina`, `minigraph`, `AMRFinderPlus`, `GECCO`, `MACSE`, `Biopython`, `SciPy`) and the active validation queue.
4. **[`docs/DISCOVERIES_AND_CRITICAL_REPORTS.md`](docs/DISCOVERIES_AND_CRITICAL_REPORTS.md)** — Auto-synced index of all **20 Empirically Verified Preprints**, recent Lab Discoveries, and dataset provenance mappings.
5. **[`RUNIE_CHANGELOG.md`](RUNIE_CHANGELOG.md)** — Chronological engineering log of every calibration gap identified, native Rust fix implemented, and multi-stage `DAG` integrated (`RN-2026-001` through `RN-2026-036`).

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
# Tool #1 (v4 Upgrade): 6-Stage Native Rust MAFFT-Mimetic Aligner (L-INS-i / E-INS-i / Q-INS-i / FFT-NS-i) vs. Center-Star Gotoh & MAFFT v7.526
./target/release/bio-suite align-compare --out runs/msa-native-vs-mafft-calibration.json
./target/release/bio-suite mafft-native --strategy auto --out runs/dikarya-its-aligned.fasta

# Tool #36: Coalescent Bayesian Skyline MCMC Sampler vs. BEAST v1.10.4
./target/release/bio-suite mcmc-skyline --out runs/mcmc-skyline-beast-calibration.json

# Tool #35: Relaxed-Clock Constrained Least-Squares Time-Tree (LSD2) vs. IQ-TREE v2.4.0 & SciPy
./target/release/bio-suite timetree-compare --out runs/timetree-lsd2-scipy-calibration.json

# Tool #34: 5-Stage Adaptive Multi-Biotype Super-BLAST vs. Biopython v1.88 & SciPy v1.18.1
./target/release/bio-suite super-blast --out runs/super-blast-biopython-scipy-calibration.json

# Tool #19: 16-RP Supermatrix & Dayhoff-6 / SR4 Compositional Recoding vs. IQ-TREE v2.4.0
./target/release/bio-suite rp16-compare --out runs/rp16-supermatrix-iqtree-calibration.json
```

---

## Empirically Verified Preprints (`PREPRINT-01`–`PREPRINT-20`)

<!-- AUTO-FRONTPAGE-HIGHLIGHTS:START -->
| ID & Status | Verified Preprint Title & Live Reader Link | Master Concept Zenodo DOI | Signal | Primary Validation Set |
| :--- | :--- | :--- | :--- | :--- |
| **#134** `✅ Verified Preprint` | **[Global Phylogenomic Radiation of Bemisia Whiteflies and Allied Aleyrodidae (Hemiptera: Sternorrhyncha), Native Rust Quartet Coalescent Species-Tree Inference, Symbiotic and Cytochrome P450 Selection, and Cross-Regional vs. Multi-Locus Demographic Skylines](https://crvscience.com/labs/discoveries/bemisia-whitefly-global-radiation)** | [`10.5281/zenodo.24093136`](https://doi.org/10.5281/zenodo.24093136) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-062-BEMISIA-GLOBAL-RADIATION-544RECS` |
| **#133** `✅ Verified Preprint` | **[Global Phylogenomic Radiation of Anopheles Mosquitoes (Culicidae: Anophelinae), Native Rust Quartet Coalescent Species-Tree Inference, Chemosensory and Insecticide-Resistance Selection, and Cross-Regional vs. Multi-Locus Demographic Skylines](https://crvscience.com/labs/discoveries/anopheles-global-radiation-native-biosuite)** | [`10.5281/zenodo.24078829`](https://doi.org/10.5281/zenodo.24078829) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-061-ANOPHELES-GLOBAL-RADIATION-544RECS` |
| **#132** `✅ Verified Preprint` | **[Phylogenomic Evolution of Asian Horseshoe Bats (Rhinolophidae), Quartet Coalescent Species-Tree Inference (*ASTRAL vs. Native Rust), Acoustic and Sarbecovirus Receptor Co-Evolution, and Single-Species Karst Population Dynamics](https://crvscience.com/labs/discoveries/rhinolophidae-asia-radiation-astral-skyline)** | `Verified` | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-060-RHINOLOPHIDAE-ASIA-RADIATION-212RECS` |
| **#131** `✅ Verified Preprint` | **[Phylogenomic Evolution of the North American Rodent Radiation, Zoonotic Reservoir Competence, and Coalescent Population Dynamics Across Sciuromorpha, Castorimorpha, and Myomorpha Under FreeRate and Bayesian Skyline Models](https://crvscience.com/labs/discoveries/rodent-north-america-radiation-skyline)** | [`10.5281/zenodo.24066498`](https://doi.org/10.5281/zenodo.24066498) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-059-RODENT-NORTH-AMERICA-RADIATION-212RECS` |
| **#130** `✅ Verified Preprint` | **[Phylogenomic Evolution of Trans-Atlantic Dispersal, Ecosystem Function, and Coalescent Population Dynamics Across African and Neotropical Primates (Anthropoidea: Platyrrhini) Under FreeRate and Bayesian Skyline Models](https://crvscience.com/labs/discoveries/primate-africa-americas-radiation-skyline)** | [`10.5281/zenodo.24061361`](https://doi.org/10.5281/zenodo.24061361) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-058-PRIMATE-AFRICA-AMERICAS-RADIATION-212RECS` |
| **#129** `✅ Verified Preprint` | **[Phylogenomic Evolution of Ecosystem Function and Coalescent Population Dynamics Across Six Families of Beetles (Coleoptera) Under FreeRate and Bayesian Skyline Models](https://crvscience.com/labs/discoveries/coleoptera-ecosystem-function-skyline)** | [`10.5281/zenodo.24058170`](https://doi.org/10.5281/zenodo.24058170) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-057-COLEOPTERA-FUNCTION-SKYLINE-212RECS` |
| **#128** `✅ Verified Preprint` | **[RuNIE Long-Form Bayesian MCMC Skyline & the 5-Stage Neanderthal Population Crash (650 to 39.2 ka BP): Coupled 95% Credible Chronogram, Pristine Phylogram, and Effective Population Size Trajectory Benchmarked Against *BEAST across an 87,227-bp Alignment](https://crvscience.com/labs/discoveries/neanderthal-mcmc-skyline-crash)** | [`10.5281/zenodo.24055713`](https://doi.org/10.5281/zenodo.24055713) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-056-NEANDERTHAL-MCMC-SKYLINE-64TAXA-1.08MB` |
| **#127** `✅ Verified Preprint` | **[Multifactorial N×L Phylogenomic Sensitivity & Clock-Time Curve: 3-Way RuNIE Comparison (100% Native Rust vs. 100% External Binaries vs. Unconstrained Hybrid) Mapping Minimum Alignment-Length Horizons, the 90% MAFFT Alignment Bottleneck, and Saturation Plateaus Across 100 Species](https://crvscience.com/labs/discoveries/multifactorial-nxl-clock-curve-3way)** | [`10.5281/zenodo.24049859`](https://doi.org/10.5281/zenodo.24049859) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-055-MULTIFACTORIAL-NXL-SENSITIVITY-3WAY-100TAXA-1.53MB` |
| **#126** `✅ Verified Preprint` | **[Pancrustacea-to-Hexapoda Tree of Life (538.4 to 0.0 Ma): A 3-Way RuNIE Pipeline Comparison (100% Native Rust vs. 100% External Binaries vs. Unconstrained Hybrid Optimal ln L) Resolving Remipedia Sisterhood, Epipodite-to-Wing Homology, and Terrestrialization Across 100 Species](https://crvscience.com/labs/discoveries/pancrustacea-hexapoda-3way-runie)** | [`10.5281/zenodo.24045575`](https://doi.org/10.5281/zenodo.24045575) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-054-PANCRUSTACEA-HEXAPODA-ORIGIN-3WAY-122RECS-1.53MB` |
| **#125** `✅ Verified Preprint` | **[Mosquito Tree of Life: Triassic–Cenozoic Culicidae Phylogenomics (242.0 to 0.0 Ma; Crown Culicidae 145.0 Ma), Single Origin of Hematophagy, Salivary D7 & Apyrase Recruitment, and Secondary Loss in Toxorhynchites](https://crvscience.com/labs/discoveries/culicidae-phylogenomics-hematophagy-d7)** | [`10.5281/zenodo.24039858`](https://doi.org/10.5281/zenodo.24039858) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-053-CULICIDAE-PHYLOGENOMICS-HEMATOPHAGY-D7-120RECS-1.64MB` |
| **#124** `✅ Verified Preprint` | **[Beetle Tree of Life: Permian–Cenozoic Coleoptera Phylogenomics (308.0 to 0.0 Ma; Crown Coleoptera 292.0 Ma), Suborder Divergence, and Microbial PCWDE Horizontal Gene Transfer (GH28, GH45, GH5)](https://crvscience.com/labs/discoveries/coleoptera-phylogenomics-pcwde-hgt)** | [`10.5281/zenodo.24031323`](https://doi.org/10.5281/zenodo.24031323) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-052-COLEOPTERA-PHYLOGENOMICS-PCWDE-HGT-101RECS-1.38MB` |
| **#123** `✅ Verified Preprint` | **[Lepidoptera Phylogenomics: Devonian–Cenozoic Time-Tree (363.76 Ma), Angiosperm-Coupled Ditrysian Radiation, and COX1 Barcode Saturation Horizon](https://crvscience.com/labs/discoveries/lepidoptera-ditrysia-timetree-cox1-gap)** | [`10.5281/zenodo.23994485`](https://doi.org/10.5281/zenodo.23994485) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-051-LEPIDOPTERA-DITRYSIA-TIMETREE-COX1-72TAXA-1.92MB` |
| **#122** `✅ Verified Preprint` | **[Hymenoptera Phylogenomics: Permian–Cenozoic Time-Tree (281.23 Ma), Single Parasitoid Origin, Eusocial Caste Selection, and COX1 Horizon](https://crvscience.com/labs/discoveries/hymenoptera-radiation-eusociality-timetree)** | [`10.5281/zenodo.23992890`](https://doi.org/10.5281/zenodo.23992890) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-050-HYMENOPTERA-RADIATION-EUSOCIALITY-TIMETREE-105TAXA-1.53MB` |
| **#121** `✅ Verified Preprint` | **[Synapsida and Mammaliaformes Phylogenomics: Carboniferous–Cenozoic Time-Tree, Dual Middle-Ear Detachment, Lactation Selection, and COX1 Horizon](https://crvscience.com/labs/discoveries/synapsida-mammaliaformes-timetree)** | [`10.5281/zenodo.23991261`](https://doi.org/10.5281/zenodo.23991261) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-049-SYNAPSIDA-MAMMALIAFORMES-CROWN-MAMMALIA-TIMETREE-66TAXA-1.82MB` |
| **#120** `✅ Verified Preprint` | **[Sarcopterygii and Tetrapoda Phylogenomics: Devonian–Carboniferous Time-Tree, Dipnoi Sisterhood, Autopod Selection, and COX1 Horizon](https://crvscience.com/labs/discoveries/tetrapod-fin-limb-timetree)** | [`10.5281/zenodo.23988727`](https://doi.org/10.5281/zenodo.23988727) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-048-SARCOPTERYGII-TETRAPODA-FIN-TO-LIMB-TIMETREE-64TAXA-1.20MB` |
| **#119** `✅ Verified Preprint` | **[Paleoproterozoic Primary Plastid Endosymbiosis, Basal Gloeomargarita Sisterhood, and the 2536.2–0 Ma Archaeplastida vs. Paulinella Time-Tree Chronogram](https://crvscience.com/labs/discoveries/plastid-cyanobacteria-timetree)** | [`10.5281/zenodo.23979910`](https://doi.org/10.5281/zenodo.23979910) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-047-ARCHAEPLASTIDA-PLASTID-CYANOBACTERIA-TIMETREE-64TAXA-1.76MB` |
| **#118** `✅ Verified Preprint` | **[The Metazoan Cambrian Explosion Tree: Full-Length Phylogenomics, Hox/Pax6/Biomineralization/Opsin Selection, and Cryogenian–Cambrian Time-Tree Dating Across 56 Metazoan and Choanoflagellata Mitogenomes and 50 Nuclear Loci](https://crvscience.com/labs/discoveries/cambrian-explosion-bilaterian-tree)** | [`10.5281/zenodo.23977970`](https://doi.org/10.5281/zenodo.23977970) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-046-METAZOA-CAMBRIAN-EXPLOSION-TOTAL-EVIDENCE-68TAXA-1.84MB` |
| **#117** `✅ Verified Preprint` | **[Phylogenomic Divergence and Post-K–Pg Adaptive Radiation Across Mammalia Under FreeRate, Quartet Coalescent, and Relaxed-Clock Models](https://crvscience.com/labs/discoveries/mammalia-kpg-supertree-chronogram)** | [`10.5281/zenodo.23976160`](https://doi.org/10.5281/zenodo.23976160) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-045-MAMMALIA-PLACENTALIA-KPG-SUPERTREE-66TAXA-1.72MB` |
| **#114** `✅ Verified Preprint` | **[The Basal Metazoan Time-Tree: Resolving Ctenophora Sister-Group Topology, Poriferan Paraphyly vs. Monophyly, and Cnidarian–Bilaterian Divergence Across 64 Holozoan and Metazoan Genomes (1.64 Mb)](https://crvscience.com/labs/discoveries/cteno-porifera-cnidaria-timetree)** | [`10.5281/zenodo.23974572`](https://doi.org/10.5281/zenodo.23974572) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-044-BASAL-METAZOA-CTENO-PORIFERA-CNIDARIA-TIMETREE-64TAXA-1.64MB` |
| **#113** `✅ Verified Preprint` | **[Out of Africa Through the L3-to-M/N Mitogenome Bottleneck: Time-Tree Chronograms, Haplotype k-Core Hubs, and OXPHOS Purifying Selection Across 58 Hominin Genomes](https://crvscience.com/labs/discoveries/human-ooa-l3-mn-mitogenome)** | [`10.5281/zenodo.23970546`](https://doi.org/10.5281/zenodo.23970546) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-043-HUMAN-OOA-L3-MN-MITOGENOME-58TAXA-1.26MB` |
| **#110** `✅ Verified Preprint` | **[RuNIE (Rust Native Inference Engine): An Open-Source Architecture for Autonomous Bioinformatics Pipeline Synthesis, Blinded Validation, and Bayesian Scope-Depth Auditing](https://crvscience.com/labs/discoveries/runie-automated-bioinformatics-preprint)** | [`10.5281/zenodo.23969275`](https://doi.org/10.5281/zenodo.23969275) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-001 through VAL-SET-043 (43 Curated Dev/Reserve Datasets)` |
| **#109** `✅ Verified Preprint` | **[Reconstructing the Archosauria-to-Crown-Aves Evolutionary Time-Tree (246.5 to 0.0 Ma) Across 62 Extant, Holocene Ancient-DNA, and Mesozoic Fossil Taxa: Paravian Body-Mass Miniaturization, the K–Pg Neoaves Radiation, and Convergent Ratite Flight Loss](https://crvscience.com/labs/discoveries/theropoda-crown-aves-timetree-kpg)** | [`10.5281/zenodo.23971260`](https://doi.org/10.5281/zenodo.23971260) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-042-THEROPODA-AVES-TIMETREE-62TAXA-1.54MB` |
| **#108** `✅ Verified Preprint` | **[Reconstructing the Primate Evolutionary Time-Tree (73.4 to 0.0 Ma) Across 56 Extant and Fossil Taxa: K–Pg Radiation, the Hominoid Rate Slowdown, and Pleistocene Ancient-DNA Tip-Dating](https://crvscience.com/labs/discoveries/primate-evolution-timetree-chronogram)** | [`10.5281/zenodo.23971921`](https://doi.org/10.5281/zenodo.23971921) | 🧪🧪🧪🧪🧪 (`5/5`) | `VAL-SET-041-PRIMATE-TIMETREE-56TAXA-1.42MB` |
<!-- AUTO-FRONTPAGE-HIGHLIGHTS:END -->

---

## Repository Structure

```text
crvscience-labs/
├── Cargo.toml                            # Rust crate manifest (candle-core, candle-nn, rayon, serde)
├── README.md                             # High-level overview, quick-start examples & repo layout
├── RUNIE_CHANGELOG.md                    # Chronological RuNIE tool-modification & pipeline ledger (RN-2026-001..036)
├── docs/                                 # Detailed technical & validation documentation
│   ├── RUNIE_ARCHITECTURE_AND_STUDIES.md # Candle tensor GA-MCMC math, Step-Up loop & 20 Verified Preprints
│   ├── VALIDATION_AND_CALIBRATION_PROCESS.md # 3-Tier system, 0🧪–5🧪 gate, Blinded SHA-256 audit & VAL-SET-001..059
│   ├── TOOL_BENCHMARK_MATRIX.md          # Full 36-tool Native Rust vs. External Counterpart benchmark table
│   └── DISCOVERIES_AND_CRITICAL_REPORTS.md # Auto-synced Verified Preprints & Lab Discoveries index
├── src/
│   ├── main.rs                           # CLI entry point & subcommand dispatcher
│   ├── lib.rs                            # Library module exports
│   ├── runie.rs                          # RuNIE v2.0.1: Candle tensor Markov policy + GA-MCMC + Stage 4 ModelTest
│   ├── runie_tools.rs                    # Synthesized & calibrated native Rust domain operators
│   ├── mcmc_skyline.rs                   # Tool #36: De novo Bayesian coalescent skyline MCMC & relaxed clock vs. BEAST
│   ├── timetree_calibration.rs           # Tool #35: Constrained least-squares molecular clock (LSD2) vs. IQ-TREE & SciPy
│   ├── align_viewer.rs                   # Native multi-track sequence alignment & structural motif viewer
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
│   ├── popgen_calibration.rs             # Haplotype MSN, Hudson F_ST & Tajima's D vs. SciPy (VAL-SET-013)
│   ├── dca_coevo_calibration.rs          # Candle mfDCA precision matrix inversion & APC vs. SciPy (VAL-SET-012)
│   ├── kabsch_tm_calibration.rs          # Jacobi SVD Kabsch SO(3) superposition & TM-score vs. SciPy (VAL-SET-011)
│   ├── repeat_motif_calibration.rs       # Spidroin/fibroin repeat grammar & ridge regression vs. SciPy (VAL-SET-010)
│   ├── codon_opt_calibration.rs          # Synonymous codon optimizer & Turner-2004 MFE vs. ViennaRNA (VAL-SET-009)
│   ├── rp16_supermatrix_calibration.rs   # Dayhoff-6 / SR4 compositional recoding vs. IQ-TREE v2.4.0 (VAL-SET-008)
│   ├── kaiabc_allostery_calibration.rs   # Jacobi Kirchhoff GNM allostery & RK4 KaiABC oscillator vs. SciPy (VAL-SET-007)
│   ├── mk_morpho_calibration.rs          # Lewis Mk+ASC morphological pruning & PIC vs. IQ-TREE v2.4.0 (VAL-SET-006)
│   ├── pseudo_loss_calibration.rs        # 3-frame codon-to-AA frameshift aligner & RELAX vs. MACSE & PAML (VAL-SET-005)
│   ├── supervisor.rs                     # Autonomous Multi-Threaded Pipeline Supervisor & Self-Healing Loop-Back Engine
│   ├── figure_generator.rs               # Locked Preprint 18 & 21 Figure (.svg + .png + .pdf) & LaTeX (.tex) Generator
│   ├── bio_suite_gui.rs                  # 5-Tab Native Win32 Desktop GUI (Pipeline DAG, Super-BLAST Audit, AlignView, Trees, 3D PDB)
│   ├── mafft_native.rs                   # 6-stage Native Rust MAFFT-Mimetic aligner (FFT-NS + hat3 SW consistency + 5.8S ITS lock + 2-Tree UPGMA DP + dvtditr)
│   ├── phylo.rs & msa_calibration.rs     # 3-way MSA calibration (MAFFT-Mimetic v4 vs. Center-Star Gotoh v3 vs. MAFFT v7.526) & K2P/TN93 distances
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
│   ├── sweep.rs                          # Cross-kingdom calibration sweep & 0🧪–5🧪 scoring gates
│   └── bin/                              # Dedicated 3-way RuNIE & clade phylogenomic audit binaries
└── tests/
    └── core_fidelity.rs                  # Automated unit & integration test suite
```

*(Note: Local working directories `cache/`, `runs/`, and `tools/` are git-ignored so the repository remains lightweight.)*

---

## Autonomous Pipeline Supervisor (`bio-suite-supervisor`) & 5-Tab Desktop GUI

`bio-suite-supervisor` oversees end-to-end study execution—from NCBI Entrez / FASTA ingestion and `Super-BLAST` self-healing annotation hygiene to multi-threaded full-length alignment, 32-chain `RuNIE` factorial sweeps, 3D `.pdb` homology superposition, locked Preprint 18 & 21 figure generation (`.svg` + `.png` + `.pdf`), and LaTeX (`main.tex`) preprint compilation:

```powershell
# Execute an autonomous 3-Way Multifactorial run and export all 5 GUI tab snapshots
cargo run --release --bin bio-suite-supervisor -- run `
  --preset anopheles-3way `
  --design 3way-multifactorial `
  --gui-snapshot-png runs/runie-3way-multifactorial-anopheles-mosquitos/figures/gui-overview.png

# Launch the interactive 5-Tab Win32 Desktop GUI
cargo run --release --bin bio-suite-supervisor -- gui --run-id runie-3way-multifactorial-anopheles-mosquitos

# List or inspect managed runs in crvscience-labs/runs/
cargo run --release --bin bio-suite-supervisor -- list
cargo run --release --bin bio-suite-supervisor -- inspect --run-id runie-3way-multifactorial-anopheles-mosquitos
```

See [**`docs/SUPERVISOR_AND_GUI_GUIDE.md`**](docs/SUPERVISOR_AND_GUI_GUIDE.md) for the complete command-line and GUI manual.

---

## How to Cite

If you use **`crvscience-labs` (`bio-suite`)**, **`RuNIE` (`Rust Native Inference Engine`)**, or any of the curated validation datasets (`VAL-SET-001`–`VAL-SET-059`) in your research, software, or teaching, please cite our architecture preprint:

> **White, B. P.** (2026). *RuNIE (Rust Native Inference Engine): An Open-Source Architecture for Autonomous Bioinformatics Pipeline Synthesis, Blinded Validation, and Bayesian Scope-Depth Auditing*. CRV Science Labs Preprint Monograph (`PREPRINT-01`, `v2.0.1`, 26 pp.).  
> - **Master Concept Zenodo DOI**: [`10.5281/zenodo.23969275`](https://doi.org/10.5281/zenodo.23969275)  
> - **Figshare Archive**: [`https://figshare.com/articles/preprint/_/34919469`](https://figshare.com/articles/preprint/_/34919469)  
> - **Canonical Web Edition**: [`https://crvscience.com/labs/preprints/runie-automated-bioinformatics-preprint`](https://crvscience.com/labs/preprints/runie-automated-bioinformatics-preprint)  
> - **Single-Column PDF**: [`https://crvscience.com/labs/pdfs/runie-automated-bioinformatics-preprint.pdf`](https://crvscience.com/labs/pdfs/runie-automated-bioinformatics-preprint.pdf)

### BibTeX

```bibtex
@article{white2026runie,
  author       = {White, Bryan P.},
  title        = {{RuNIE (Rust Native Inference Engine): An Open-Source Architecture for Autonomous Bioinformatics Pipeline Synthesis, Blinded Validation, and Bayesian Scope-Depth Auditing}},
  year         = {2026},
  month        = nov,
  publisher    = {CRV Science Labs / Zenodo},
  version      = {2.0.1},
  doi          = {10.5281/zenodo.23969275},
  url          = {https://doi.org/10.5281/zenodo.23969275},
  note         = {Preprint Monograph PREPRINT-01; Master Concept Zenodo DOI: 10.5281/zenodo.23969275; Web Edition: https://crvscience.com/labs/preprints/runie-automated-bioinformatics-preprint}
}
```

---

## License

Released under the **GNU General Public License v3.0 (`GPL-3.0-or-later`)**. See [LICENSE](LICENSE) for full details.
