# `bio-suite-seqnet` — Pure-Rust 3-Stage Hybrid Sequence Analysis Neural Network

`bio-suite-seqnet` is a standalone Pure-Rust crate (GPL-3.0, compatible with `crvscience-labs` / `RuNIE`) implementing a live 3-stage cluster of sequence analysis neurons with multi-database SQLite caching, non-destructive soft-quarantine feedback loops, species-terminal-per-genus chimeric OTU construction, and clean partitioned NEXUS supermatrix export.

## 1. Three-Stage Hybrid Neuron Architecture

| Cluster Stage | Neuron Units | Bioinformatic & Neural Responsibilities |
| :--- | :--- | :--- |
| **Stage 1: Input Cluster** | `CacheWarmStartNeuron`, `EntrezDeltaFetchNeuron`, `FeedbackGapRefillNeuron` | Binds to a logical `database_id` inside SQLite ([`rusqlite` bundled](https://docs.rs/rusqlite/)). Loads cached sequences and learned neural weights immediately, then queries [NCBI Entrez E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/) (`esearch.fcgi` + `efetch.fcgi` via [`ureq`](https://docs.rs/crate/ureq/2.12.1)) using `datetype=mdat` (`mindate`/`maxdate`) watermarks so repeated runs only fetch newly added or updated accessions. Executes targeted species/locus refill queries when triggered by Stage 3 feedback. |
| **Stage 2: Hybrid Processing Cluster** | `MetadataTaxonomyNeuron`, `SequenceHygieneNeuron`, `ReadingFrameStopCodonNeuron`, `WithinLocusHomologyNeuron`, `SoftQuarantineReevalNeuron`, `HybridNeuralConsensusNeuron` | Runs parallel [`rayon`](https://github.com/huggingface/candle/blob/main/Cargo.toml) within-locus evaluations: rejects unresolved `sp.`/`cf.` names, checks nucleotide length & ambiguity (`N%`), translates CDS frames under NCBI Genetic Code 1 (nuclear) and Code 5 (invertebrate mitochondrial) to catch NUMT pseudogenes, computes 256-bin 4-mer cosine similarity + Kimura 2-Parameter (K2P) robust MAD z-scores to catch sample swaps, and scores 8 normalized features through a 2-layer differentiable MLP (`8 -> 16 -> 1`, with optional [`candle-core`](https://docs.rs/candle-core/) tensor execution). Soft-quarantined sequences stay in `quarantine_ledger` with a decaying `skepticism_prior` rather than permanent deletion. |
| **Stage 3: Output & Valuation Cluster** | `ChimericSpeciesOtuNeuron`, `BetweenLocusValuationNeuron`, `GotohMsaNeuron`, `NexusExportNeuron` | Assembles separate **species-level terminals per Genus**, enforcing `>= min_species_per_genus` (default `2`) and `>= 90%` locus coverage per retained OTU (`min_locus_coverage_ratio = 0.90`). Pools multiple conspecific specimens/vouchers into chimeric species terminals when a single specimen lacks full coverage (recording every NCBI accession and voucher in `artificial_otu_accessions`), emits `GapRefillDirective` feedback for near-threshold species, runs parallel affine-gap Gotoh MSA, and writes a clean partitioned [NEXUS supermatrix](https://www.phylo.org/tools/nexus_document.pdf). |

## 2. Multi-Database SQLite Caching & Soft-Quarantine Rules (`schema.sql`)

* **Explicit `database_id` Isolation (`database_registry`)**: Every run binds to an explicit `database_id` (for example, `CRV-COLEOPTERA-ML25-DB01`). You can maintain separate `.sqlite` files per taxonomic project or host multiple logical `database_id` partitions inside a single SQLite file without cross-contaminating sync watermarks or neural weights.
* **Rapid Incremental Re-Convergence (`entrez_sync_watermarks` + `neural_weights_store`)**: After Pass 1 caches sequences and stores trained `HybridConsensusWeights` in SQLite, re-running the same taxonomic group checks NCBI `nuccore` strictly for records modified since `last_sync_mdat` and warm-starts the neural layer, converging in **1 epoch** when no conflicts arise.
* **Non-Destructive Soft Quarantine (`quarantine_ledger`)**:
  * Flagged records (`INTERNAL_STOP_CODON_NUMT`, `WITHIN_LOCUS_K2P_SAMPLE_SWAP`, `EXCESS_AMBIGUOUS_BASES`, `LOW_HYBRID_NEURAL_CONSENSUS`) are **strictly excluded** from `BEGIN TAXA;` and `BEGIN DATA;` in the exported NEXUS file.
  * Every flagged accession is logged in SQLite (`prior_quarantine_count`, `prior_pass_count`, `skepticism_prior`) and emitted inside the `BEGIN BIOSUITE_METADATA;` audit block of the NEXUS file (`QUARANTINE_FLAG accession=...`) so you see exactly what was withheld for that group.
  * On subsequent runs, `SoftQuarantineReevalNeuron` re-checks previously quarantined accessions against expanded taxon sampling; if denser congeneric coverage confirms an accession is a valid divergent lineage rather than a contaminant, it can rehabilitate the sequence (`currently_quarantined = 0`).

## 3. Dual Execution Modes: Headless Console vs. Console + GUI Monitor

Build the crate with standard Cargo:

```bash
cargo build --release
# Optional: enable RuNIE candle-core tensor bridge
cargo build --release --features candle-backend
# Optional: build strictly headless without desktop window dependencies
cargo build --release --no-default-features
```

### Mode A: Headless Console (`--mode headless`)

Runs multi-threaded console telemetry and writes the SQLite cache and clean partitioned NEXUS file:

```bash
export NCBI_API_KEY="your_ncbi_key"
export NCBI_EMAIL="labs@crvscience.com"

./target/release/bio-suite-seqnet run \
  --taxon "Coleoptera" \
  --nuclear-loci 18 \
  --mito-loci 7 \
  --target-otus 100 \
  --otu-rank genus \
  --min-species-per-genus 2 \
  --min-coverage 0.90 \
  --db-path runs/coleoptera_ml25.sqlite \
  --db-id "CRV-COLEOPTERA-ML25-DB01" \
  --out runs/coleoptera_25loci_optimal.nex \
  --mode headless
```

### Mode B: Console + Live Desktop GUI Monitor (`--mode gui`)

Launches the [`eframe`](https://docs.rs/eframe/latest/eframe/) / [`egui`](https://docs.rs/egui/latest/egui/) real-time monitor window alongside terminal logging, featuring Pause/Resume & Force-Finalize controls, 3-stage neuron activation gauges, a live **OTU x Locus 90% Occupancy Heatmap**, a **Chimeric OTU Provenance Inspector**, and a **Soft-Quarantine Ledger** viewer:

```bash
./target/release/bio-suite-seqnet run \
  --taxon "Coleoptera" \
  --nuclear-loci 18 \
  --mito-loci 7 \
  --target-otus 100 \
  --otu-rank genus \
  --min-species-per-genus 2 \
  --min-coverage 0.90 \
  --db-path runs/coleoptera_ml25.sqlite \
  --db-id "CRV-COLEOPTERA-ML25-DB01" \
  --out runs/coleoptera_25loci_optimal.nex \
  --mode gui
```

### Offline 2-Pass Benchmark & Database Audit Subcommands

```bash
# Run offline 2-pass Coleoptera benchmark (Pass 1 cold start -> Pass 2 cached instant consensus)
./target/release/bio-suite-seqnet demo-sim --mode headless
./target/release/bio-suite-seqnet demo-sim --mode gui

# Audit registered database IDs and soft-quarantined accessions in SQLite
./target/release/bio-suite-seqnet inspect-db \
  --db-path demo_coleoptera_seqnet.sqlite \
  --db-id "CRV-DEMO-COLEOPTERA-DB01" \
  --taxon "Coleoptera"

# Execute integration tests
cargo test
```
