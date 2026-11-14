# seqnet change log

Changes to `bio-suite-seqnet` since it was ingested from the `bio-suite-seqnet.tar` archive. Newest first. Each entry records what changed, why, and how it was verified.

---

## 2026-11-13: Mitogenome mode (commit `b689370`)

**Why:** To reproduce the PREPRINT-18 dataset dimensions (64 species x 16 loci, mitochondrial loci from complete mitogenomes). Per-locus title queries cannot do this, because mitogenome records are titled "... mitochondrion, complete genome" and never name individual genes.

### Added
- **`ConsensusConfig::mitogenome_loci`** (default `false`) and **`mitogenome_retmax_per_stratum`** (default `60`).
  - When enabled, stage 1 sends one Entrez query per stratum for complete or partial mitochondrial genomes. The query uses GenBank only, because RefSeq `NC_` records duplicate the GenBank submission.
  - Each panel mitochondrial locus is sliced from the record's CDS/rRNA feature table. All mito loci for a species therefore come from one specimen.
  - Sync watermarks are stored under the key `__mitogenome__` for each stratum.
- **`entrez::build_mitogenome_query`**, **`parse_mitogenome_insdseq`** and **`extract_location`**.
  - Supported locations: `a..b`, `complement(...)`, `join(...)`/`order(...)` (nested), and partial `<`/`>` ends.
  - Remote locations (`ACC:a..b`) are rejected.
  - Each locus is taken from the first matching feature in a record and is kept only if its length is within the locus spec's bounds.
- **`entrez::mito_feature_canonical`** maps `gene`/`product` qualifier spellings (`cox1`, `COI`, `cytochrome c oxidase subunit I`, `nad4L`, `l-rRNA`, `12S ribosomal RNA`, ...) to canonical panel loci.
- **Mitogenome-anchored nuclear queries.** After the stratum loop, nuclear panel loci are queried by species name for every species that has a mitogenome. Stratum queries for nuclear loci rarely sample the same species as the mitogenomes: in the first run, 28S landed on 0 of the 64 kept species.

### Changed
- **Gap-refill skips mito loci in mitogenome mode.** Refilling one gene from a different specimen would break the one-specimen-per-species guarantee.
- **Entrez client refactor.** The client now has shared `with_retries`, `esearch_accessions` and `efetch_insdseq` helpers. Mitogenome efetch uses 40 records per POST; locus efetch stays at 100.
- **`parse_insdseq_xml` shares a new `build_record`** for record metadata and species-name screening.

### Fixed
- **XML entities in tag values are now decoded** (`&lt;`, `&gt;`, `&amp;`, `&quot;`, `&apos;`). Partial feature locations arrive as `&lt;1487..3019`, which made 206 of 367 mitogenomes lose their COI (most COI annotations are 5'-partial).
- **Species-name screen now also rejects `gen.` and `nr.`.** Before this, records such as "Nemopterinae gen. 1 sp. 1" became the binomial "Nemopterinae gen.".

### Verified
- **New unit tests:**
  - `feature_locations_cover_complement_join_and_partials`
  - `mito_feature_names_map_to_panel_loci`
  - `mitogenome_record_is_sliced_into_panel_loci`, whose fixture uses an escaped `&lt;` location.
- All seqnet tests and `convergence_test` pass.
- **End-to-end run:** `bio-suite seqnet-runie --benchmark benchmarks/coleoptera_preprint18_design.json`.
  - Input: 64 species x 16 loci, 15,931 bp supermatrix, 89.7% matrix fill.
  - The supermatrix ML tree recovered 5/5 scored published clades.
  - The quartet species tree recovered 4/5.

---

## 2026-11-13: Ingestion as RuNIE stage 1 (commit `32174d5`)

**Why:** seqnet became stage 1 of RuNIE pipelines and supersedes `fetch::PoliteCacheManager` for new pipelines. As delivered, the archive did not compile, and it returned no data from live NCBI queries. Its tests passed only on synthetic data.

### Packaging
- Added as the workspace crate `seqnet/`. `bio-suite` depends on it with `default-features = false`, so no GUI is built.

### Fixed (archive as delivered)
- **Did not compile.** In `align_nexus.rs`, the center-star reference sequence selection had a type error.
- **Run-ID collisions.** Back-to-back runs in one process collided on the `neural_runs.run_id` primary key. The ID now includes milliseconds and a per-process counter.
- **Live fetches parsed to nothing.**
  - `efetch rettype=gb` returns GBSeq XML, which the INSDSeq parser cannot read. It now uses `rettype=gbc`.
  - Stage 1 now fails loudly if a non-empty response contains no `<INSDSeq>`.
- **Wrong taxon ID.** The first `db_xref` was used as the taxon ID even when it was a `BOLD:` reference. The parser now finds the `taxon:` entry.
- **Quarantine strikes compounded across epochs.** Strikes are now read once per run, so the filter no longer tightens until it rejects nearly everything.
- **Neural consensus self-training collapse.**
  - Training batches are now class-balanced.
  - Weights stay fixed within a run, with one feedback step after the loop.
- **Premature convergence.** A consensus plateau no longer counts while gap-refill is still adding OTUs.
- **One failed Entrez pair aborted all strata.** Requests now retry with backoff. A failing (stratum, locus) pair is skipped, and stage 1 goes offline only after 4 consecutive failures.

### Added
- **Explicit locus panel** (`locus_names`, `build_named_locus_panel`) and **stratified stage-1 sampling** (`sampling_taxa`, with watermarks per stratum and locus).
- **Configurable refill limits:** `max_refill_directives`, `max_refill_loci_per_directive` and `refill_min_coverage`.
- **Region-coherence gate** (`region.rs`, `region_min_core_coverage`, `region_core_support`). It quarantines records that fall outside the locus's shared core region, such as non-overlapping COI fragments, with the reason `REGION_OUTSIDE_LOCUS_CORE`. It applies only when at least 6 records map and the anchor sequence has a complexity of at least 0.8.
- **Query hygiene.** Locus synonyms are matched in `[Title]` only. Queries also:
  - enforce the compartment with `mitochondrion[Filter]` / `NOT mitochondrion[Filter]`;
  - exclude mRNA (`NOT biomol_mrna`) for rRNA loci, and RefSeq predicted models (`NOT srcdb_refseq_model`) and WGS records for all loci;
  - drop the bare `SSU`/`LSU` synonyms, which also name mitochondrial rRNA.
- **CLI flags** on `run`: `--loci`, `--sampling-taxa`, `--retmax`, `--refill-directives` and `--refill-min-coverage`.

### Not verified
- The default-feature GUI build (`eframe`) has not been compiled or tested since ingestion.
