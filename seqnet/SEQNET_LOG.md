# seqnet change log

Changes to `bio-suite-seqnet` since it was ingested from the `bio-suite-seqnet.tar` archive. Newest first. Each entry records what changed, why, and how it was verified.

---

## 2026-11-14: Compact mitogenome bounds & per-locus fallback for Nematoda (PREPRINT-27)

**Why:** Preparing `bio-suite-seqnet` Stage 1 for phylum-wide Nematoda (`crvscience-labs/benchmarks/nematoda_preprint27_design.json`), where Chromadorean mitogenomes lack `atp8` (`~6.5–15 kb`) and several agricultural pest families (`Hoplolaimidae`, `Tylenchulidae`, `Anguinidae`, `Criconematidae`, `Trichodoridae`) have sparse complete mitogenomes on NCBI but dense nuclear `SSU_18S` / `LSU_28S` and single-locus `COI` / mitochondrial rRNA submissions.

### Changed
- **Widened `build_locus_panel` bounds & rRNA title synonyms (`seqnet/src/config.rs`):** Adjusted minimum length bounds for partial/divergent metazoan markers (`COI` `350..1650`, `SSU_18S` `520..2100`, `LSU_28S` `500..3800`, etc.) and added `"small subunit ribosomal RNA"` / `"large subunit ribosomal RNA"` title synonyms.
- **Compact/multipartite mitogenome length window (`seqnet/src/entrez.rs`):** Expanded `build_mitogenome_query` from `10000:30000[SLEN]` to `6500:30000[SLEN]` to capture compact or multipartite nematode mitochondrial circles.
- **Per-locus mitochondrial fallback when mitogenomes are scarce (`seqnet/src/neurons.rs`):** When `use_mitogenomes` returns fewer than `12` distinct species for a stratum, automatically supplements with per-locus `COI`, `COII`, `CYTB`, `rrnL_16S`, and `rrnS_12S` queries alongside nuclear loci while leaving mitogenome-rich strata unchanged.

### Verified
- Added `nematoda_family_strata_and_compact_mitogenome_bounds` unit test in `seqnet/src/entrez.rs` and verified via `cargo test --workspace --release`.

---

## 2026-11-14: Isoptera (PREPRINT-26) validation sweep & homology/stratum fixes

**Why:** Replicated the PREPRINT-25 phylogenomic pipeline on the global radiation of Isoptera (64 species x 16 loci across 14 Blattodea/Termitidae taxonomic strata, pre-registered in `crvscience-labs/benchmarks/isoptera_preprint26_design.json`) using `bio-suite-seqnet` as Stage 1 ingest, and fixed six bugs surfaced during cold fetch, Stage 2 homology QC, and Stage 2–4 phylogenomic validation.

### Fixed
- **Epifamily `-oidae` false stratum match and `Termitidae` subfamily resolution (`src/entrez.rs`).**
  - `extract_family_or_stratum` previously stopped at the first taxon suffix matching `"idae"` or `"oidae"`, causing Blattodea epifamily `Termitoidae` to shadow downstream Termitidae subfamilies (`Macrotermitinae`, `Apicotermitinae`, `Cubitermitinae`, `Termitinae`, `Nasutitermitinae`, `Syntermitinae`) and NCBI's informal `"Cubitermes group"` clade to fail lookup.
  - Updated `extract_family_or_stratum` to ignore `Termitoidae` (`!ends_with("oidae")`), normalize `"Cubitermes group"` -> `"Cubitermitinae"` prior to `-inae` checks, and scan lineage tokens in reverse order prioritizing exact `sampling_taxa` matches via `parse_mitogenome_insdseq_with_strata` and `parse_insdseq_xml_with_strata`.
- **Extended mitochondrial feature synonyms (`src/entrez.rs`).**
  - Extended `mito_feature_canonical` to recognize `nadh1`..`nadh6`, `nadh4l`, `rrn16`/`16srrn`/`16sribosomalrnarrnl`, and `rrn12`/`12srrn`/`12sribosomalrnarrns`.
- **Mitogenome-anchored nuclear warm-cache, quota ranking & epoch refill dedup (`src/neurons.rs`).**
  - Stage 1 anchored nuclear queries now scan both `cached_initial` and `newly_fetched`, skip candidates unable to reach `min_required_loci`, rank descending by `mito_count` (capped at `24` per stratum), inherit resolved stratum onto nuclear records lacking family metadata, and persist `attempted_refills: HashSet<(String, String)>` across `Epoch > 1` gap refills.
- **Unaligned index-0 K2P false saturation (`src/neurons.rs`).**
  - `compute_approx_k2p` previously compared raw unaligned byte positions `0..min_len`, causing genuine homologs with a 3-bp start offset (*Blattella germanica* vs. *Mastotermes darwiniensis* `COI`, `83.0%` identical at shift `+3`), internal indels (`ND4L`, `rrnL_16S`, `rrnS_12S`), or minus-strand submissions (`GU205094.1`) to saturate at `k2p = 1.5000` despite `kmer_sim = 0.91–0.97`, triggering `488` false-positive quarantines (primarily `LOW_HYBRID_NEURAL_CONSENSUS` on `ND4L`, `COI`, `ATP8`, `ND3`, `COIII`, `ND6`).
  - Added `compute_k2p_with_index` (preserving `compute_approx_k2p` signature) using `AnchorIndex::best_diagonal` (`>= 4` hits across both strands) + banded semi-global Needleman-Wunsch alignment (`match = +2, mismatch = -4, gap = -6`, band `±36 bp`, capped at `2,400 bp`) so only homologous aligned columns enter Kimura 2-Parameter counting.
- **Outgroup region-anchor hijacking & overly strict `K=10` core mapping (`src/neurons.rs`, `src/region.rs`).**
  - Selecting `region_anchor` solely by `max_by_key(seq_length)` picked the minus-strand outgroup cockroach *Blattella germanica* `LSU_28S` (`GU205094.1`, `2,598 bp`, lowest centroid 4-mer similarity `0.8094` vs. `0.9499` for *Pseudacanthotermes spiniger* `EU253824.1`), leaving 238 clean termite `28S` records unmapped and isolating *Sigmella normalis* cockroaches as a spurious core (`303/326` `LSU_28S` quarantined). Simultaneously, exact `K = 10` k-mer seeds rejected divergent short mitochondrial loci (`ATP8` `159 bp`, `ND4L` `288 bp`, `ND6` `492 bp`).
  - Switched `region_anchor` selection to maximize centroid 4-mer similarity `sims[idx]` among clean records `>= p90_clean_len * 3 / 4`, and updated `region::AnchorIndex` to adaptive `k = 7` (`<= 550 bp`) / `k = 8` (`> 550 bp`) across top-3 diagonal bins validated by banded Smith-Waterman local homology (`norm_sw >= 0.20, raw_sw >= 36`) plus exported `pub fn revcomp` and `AnchorIndex::best_diagonal`.

### Added
- **Parallelized `seqnet-runie` Stage 2+3 & `RY3` protein-coding supermatrix (`src/seqnet_runie.rs`, `src/main.rs`).**
  - Parallelized per-locus native MAFFT, BIC ModelTest, and ML gene-tree inference via `rayon::par_iter()`, and added `RY3` third-codon-position recoded supermatrix inference (`supermatrix_ry3.fasta`, `supermatrix_ry3.nwk`, `supermatrix_ry3_ln_l`, `supermatrix_rcfv_4state`, `supermatrix_rcfv_ry3`, `validation_supermatrix_ry3`).
- **`bio-suite-supervisor` Stage 1 `seqnet` integration (`src/supervisor_isoptera.rs`).**
  - Wired `bio-suite-seqnet` Stage 1 (`seqnet_isoptera_preprint26.sqlite`, `CRV-ISOPTERA-PREPRINT26-DB01`) directly into `execute_isoptera_supervisor_study`.

### Verified
- **Unit & integration tests:** All `13` `bio-suite-seqnet` lib unit tests (including 6 new tests across `entrez::tests`, `region::tests`, and `neurons::tests`) + `convergence_test` (`1/1`) + `bio-suite` (`7/7` lib + `20/20` `core_fidelity`) pass in `--release`.
- **End-to-end validation sweep (`bio-suite seqnet-runie --benchmark crvscience-labs/benchmarks/isoptera_preprint26_design.json`):**
  - Pre-fix baseline vs. post-fix (7,807 cached NCBI sequences, `CRV-ISOPTERA-PREPRINT26-DB01`):
    - Clean records: `7,002 -> 7,490` | Quarantined: `805 -> 317` (eliminated `488` false quarantines while keeping genuine `100-N` `LSU_28S` scaffolds, non-overlapping `D1-D2` `28S`, and `INTERNAL_STOP_CODON_NUMT` rejects out).
    - `2`-epoch convergence (`297` OTUs assembled, `64` retained across `14/14` families): matrix occupancy improved `88.88% -> 93.48%`, mean consensus `0.9015 -> 0.9233`.
    - `LSU_28S` recovered from `2` taxa (skipped) to `36 taxa x 2,466 bp` (`lnL = -7780.69`); `ATP8` `63 taxa x 176 bp`; all `16/16` loci retained (`supermatrix 64 taxa x 16,675 bp`; RCFV reduced `90.6%` from `0.0396` 4-state to `0.0037` `RY3`; `msc_normalized_quartet_score = 0.9296`).
    - Cladistic check against the frozen benchmark recovered `5/6` scored clades in `supermatrix_ml`, `supermatrix_ry3_ml`, and `quartet_msc` (`Xylophagodea`, `Isoptera`, `Euisoptera`, `Neoisoptera`, `Termitidae` all `100%`/`95%` bootstrap and `0.9905–1.0000` LPP, plus monophyletic `Cryptocercidae`, `Kalotermitidae`, `Apicotermitinae`, `Cubitermitinae`, `Nasutitermitinae`, `Syntermitinae`). On `non_macrotermitine_termitidae` (`best_jaccard = 0.9565`), the newly recovered 36-taxon nuclear `LSU_28S` signal (`(Sphaerotermes_sphaerothorax, Foraminitermes_valens)` at `100%`) resolved the non-fungus-growing *Sphaerotermes sphaerothorax* (subfamily Sphaerotermitinae, cataloged under `Macrotermitinae` in legacy NCBI taxonomy) inside the non-Macrotermitinae Termitidae grade sister to the `1.00`-supported 9-species fungus-growing `Macrotermitinae sensu stricto`.

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
