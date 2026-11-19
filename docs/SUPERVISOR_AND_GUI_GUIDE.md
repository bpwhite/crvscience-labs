# `bio-suite-supervisor` & Multi-Tab `bio-suite` Desktop GUI Guide

`bio-suite-supervisor` (`src/supervisor.rs`, `src/figure_generator.rs`, `src/bio_suite_gui.rs`) is an autonomous, multi-threaded Rust orchestrator for **CRV Science Labs (`bio-suite` & `RuNIE`)**. It eliminates ad-hoc JavaScript/Python scripts by managing the entire bioinformatics study lifecycle inside compiled Rust:

1. **Flexible Cohort Ingestion**: Accepts NCBI Entrez accession numbers (`--accessions`), local FASTA files (`--fasta`), or built-in gold-standard study presets (`--preset anopheles-3way`, `coleoptera-function-18`, `rhinolophus-asia-21`, `bemisia-global-23`, `invasive-snails-24`, `nudibranchs-global-25`, `isoptera-termites-26`, `nematoda-pests-27`, `acanthochitonidae-28`, `china-australia-pests-29`, `us-ag-pests-30`, `basal-metazoa-05`) combined with a factorial **Model Design** (`3way-multifactorial`, `biosuite-only`, `external-only`, `tool-calibration`).
2. **Autonomous `Super-BLAST` Annotation Validation & Self-Healing Loop-Back Engine (95%+ Autonomous)**:
   - Audits every fetched sequence in parallel (`rayon`) via `super_blast` homology scoring against the cohort medoid, strand orientation verification (`--adjustdirection`), full-length completeness gates, 3-frame translation internal stop-codon scans (`NUMT` nuclear pseudogene detection), and terminal `5' C->T / 3' G->A` ancient-DNA deamination checks.
   - Automatically repairs recoverable errors in Rust (reverse-complement strand rescue, terminal deamination `N`-masking) and quarantines severe repository artifacts (e.g., sample swaps such as `PQ585799.1`, nuclear `NUMT` paralogs such as `KU551299.1`, partial fragments, or withdrawn/failed Entrez IDs) into the **`0🧪` Uncleaned Negative Control** track while automatically pulling verified full-length replacement sequences from the reserve pool.
   - Escalates to the user or LLM (`EscalatedToUserOrLlm`) **only** when unrecoverable cohort-wide failures (>35% unresolvable accessions with an exhausted reserve pool) occur.
3. **"Run" as the First-Class Unit of Management (`crvscience-labs/runs/<run-id>/`)**:
   - Every run is stored in a self-contained directory containing `run-manifest.json`, `downloads/`, `audit/`, `alignments/`, `trees/`, `pdb/`, `figures/`, and `manuscript/`.
4. **Locked-Format Figure Generator (`.svg` + `.png` + `.pdf`) & Preprint (`.tex`) Generator**:
   - Uses **Preprint 18 (`coleoptera-ecosystem-function-preprint`)** and **Preprint 21 (`rhinolophidae-asia-radiation-preprint`)** as locked formatting gold standards.
   - Automatically generates 6 data-type-matched figures in `.svg`, `.png`, and `.pdf` (including 3D `.pdb` Cα Kabsch structural superposition + catalytic pocket projection and pixel-locked colored DNA/AA alignment clips) alongside `manuscript/main.tex`, `manuscript/single-column-manuscript.html`, and `manuscript/<run-id>-preprint.pdf`.
5. **5-Tab Native Desktop GUI (`bio-suite-supervisor gui` / `bio-suite supervisor --gui`)**:
   - Interactive Win32 desktop suite and headless PNG/PDF snapshot renderer providing parity between console and GUI workflows.

---

## 1. Command-Line Usage (`bio-suite-supervisor` & `bio-suite supervisor`)

You can invoke the supervisor either via the dedicated binary `bio-suite-supervisor` or via the unified `bio-suite supervisor` subcommand:

### A. Execute an Autonomous Run from a Preset or Custom Dataset

```powershell
# 1. Run the flagship 3-Way Multifactorial Anopheles study (and export all 5 GUI tab PNGs)
.\target\release\bio-suite-supervisor.exe run `
  --preset anopheles-3way `
  --design 3way-multifactorial `
  --gui-snapshot-png runs/runie-3way-multifactorial-anopheles-mosquitos/figures/gui-overview.png

# 2. Run the Preprint 18 (Coleoptera Ecosystem Function) or Preprint 21 (Rhinolophus Asia Radiation) presets
.\target\release\bio-suite-supervisor.exe run --preset coleoptera-function-18 --design 3way-multifactorial
.\target\release\bio-suite-supervisor.exe run --preset rhinolophus-asia-21 --design 3way-multifactorial

# 2b. Run the Preprint 23–30 Multi-Compartment Native Bio-Suite studies (multi-locus supermatrix, Quartet MSC, 12x36 skylines, 61-state codon selection)
.\target\release\bio-suite-supervisor.exe run --preset bemisia-global-23 --design biosuite-only --generations 25000
.\target\release\bio-suite-supervisor.exe run --preset invasive-snails-24 --design biosuite-only --generations 25000
.\target\release\bio-suite-supervisor.exe run --preset nudibranchs-global-25 --design biosuite-only --generations 25000
.\target\release\bio-suite-supervisor.exe run --preset isoptera-termites-26 --design biosuite-only --generations 25000
.\target\release\bio-suite-supervisor.exe run --preset nematoda-pests-27 --design biosuite-only --generations 25000
.\target\release\bio-suite-supervisor.exe run --preset acanthochitonidae-28 --design biosuite-only --generations 25000
.\target\release\bio-suite-supervisor.exe run --preset china-australia-pests-29 --design biosuite-only --generations 25000
.\target\release\bio-suite-supervisor.exe run --preset us-ag-pests-30 --design biosuite-only --generations 25000

# 3. Run with custom NCBI Accession numbers and a custom Run ID (study title is derived automatically from the preset/accession cohort)
.\target\release\bio-suite-supervisor.exe run `
  --run-id runie-custom-culicidae-audit `
  --design 3way-multifactorial `
  --accessions NC_002084.1,MK628480.1,NC_028212.1,NC_028215.1,NC_038158.1,PQ585799.1,KU551299.1 `
  --generations 25000

# 4. Run from a local FASTA sequence file in Bio-Suite Native Only mode
.\target\release\bio-suite-supervisor.exe run `
  --run-id runie-local-fasta-audit `
  --fasta path\to\sequences.fasta `
  --design biosuite-only
```

### B. List & Inspect Managed Runs

```powershell
# List all managed runs in crvscience-labs/runs/
.\target\release\bio-suite-supervisor.exe list

# Inspect a specific run's manifest, NCBI anomaly log, and generated artifacts
.\target\release\bio-suite-supervisor.exe inspect --run-id runie-3way-multifactorial-anopheles-mosquitos
```

### C. Regenerate Locked-Format Figures (`.svg` + `.png` + `.pdf`) & Preprint (`.tex`)

```powershell
# Re-render all 6 figures (.svg, .png, .pdf) and manuscript/main.tex from an existing run's manifest
.\target\release\bio-suite-supervisor.exe figures --run-id runie-3way-multifactorial-anopheles-mosquitos
```

### D. Launch the Interactive 5-Tab Desktop GUI

```powershell
# Launch the interactive Win32 GUI loaded with a specific managed run
.\target\release\bio-suite-supervisor.exe gui --run-id runie-3way-multifactorial-anopheles-mosquitos

# Or execute a run and immediately open the interactive GUI upon completion
.\target\release\bio-suite-supervisor.exe run --preset anopheles-3way --design 3way-multifactorial --gui
```

---

## 2. CLI Options Reference

| Flag / Option | Description | Default |
| :--- | :--- | :--- |
| `run` / `list` / `inspect` / `figures` / `gui` | Primary supervisor mode (plus `supermatrix-viz` for Preprint 28/29 diagnostics) | `run` |
| `--run-id <id>` (`--id`) | Unique identifier for the managed run directory (`runs/<run-id>/`) | Derived from preset or `runie-3way-multifactorial-anopheles-mosquitos` |
| `--preset <name>` | Built-in dataset preset: `anopheles-3way`, `coleoptera-function-18`, `rhinolophus-asia-21`, `bemisia-global-23`, `invasive-snails-24`, `nudibranchs-global-25`, `isoptera-termites-26`, `nematoda-pests-27`, `acanthochitonidae-28`, `china-australia-pests-29`, `us-ag-pests-30`, `basal-metazoa-05` | `anopheles-3way` |
| `--design <mode>` (`--mode`, `--model-design`) | Model design: `3way-multifactorial`, `biosuite-only`, `external-only`, `tool-calibration` | `3way-multifactorial` |
| `--accessions <A,B,C>` | Comma-separated NCBI Entrez accession list | Empty (uses preset) |
| `--fasta <path>` (`--input`) | Path to an input FASTA sequence file | `None` |
| `--runs-dir <dir>` (`--out-dir`) | Root directory storing managed run units | `runs` |
| `--generations <N>` (`-g`) | MCMC generations (`-g`) for Bayesian coalescent skylines | `10000` (`25000` for `basal-metazoa-05` and Preprint 23–30 presets) |
| `--no-loopback` | Disable autonomous replacement pull during `Super-BLAST` hygiene | Loop-back enabled |
| `--no-haplotype-clusters` | Disable haplotype cluster computation | Enabled where supported |
| `--no-figures` | Skip `.svg`/`.png`/`.pdf` figure & LaTeX manuscript generation | Enabled |
| `--gui` | Open the interactive 5-tab Win32 desktop GUI | `false` |
| `--gui-snapshot-png <path>` | Render headless high-DPI PNG snapshots of all 5 GUI tabs to disk | `None` |

---

## 3. Managed Run Directory Layout (`crvscience-labs/runs/<run-id>/`)

Each run is encapsulated as a reproducible unit of management:

```text
runs/runie-3way-multifactorial-anopheles-mosquitos/
├── run-manifest.json                                    # Complete telemetry, 32-chain sweep, anomalies & paths
├── downloads/
│   └── cleaned_full_length_cohort.fasta                 # Full-length audited & strand-normalized sequences
├── audit/
│   └── super_blast_audit_and_anomalies.json             # Flagged NCBI anomalies + loop-back iteration logs
├── alignments/
│   └── aligned_supermatrix.fasta                        # Full-length multi-threaded MAFFT-Native alignment
├── trees/
│   └── solved_supermatrix_topology.nwk                  # Solved bifurcating Newick phylogeny
├── pdb/
│   └── <run-id>_homology_model.pdb                      # 3D atomic coordinates (ATOM + CONECT records)
├── figures/
│   ├── fig1-annotated-chronogram.{svg,png,pdf}          # Locked Preprint 18/21 Vertical Chronogram + Clade Cards
│   ├── fig2-circular-functional-phylogram.{svg,png,pdf} # Circular Patristic Phylogram (2*R_ARC >= 85% H)
│   ├── fig3-quartet-coalescent-species-tree.{svg,png,pdf}# Multi-Locus Quartet MSC + ABBA-BABA Introgression
│   ├── fig4-three-popset-skyline.{svg,png,pdf}          # 3-PopSet Bayesian Skyline (RuNIE vs *BEAST vs 0-Beaker)
│   ├── fig5-functional-msa-and-runie-sweep.{svg,png,pdf}# Pixel-Locked Colored Motifs + dN/dS + 32-Chain Sweep
│   ├── fig6-pdb-structure-and-alignment-clip.{svg,png,pdf} # 3D PDB Superposition + Hygiene Matrix + Dual Clip
│   └── gui-tab1..tab5-*.png                             # Rendered snapshots of all 5 GUI tabs
└── manuscript/
    ├── main.tex                                         # Complete LaTeX preprint manuscript source
    ├── single-column-manuscript.html                    # Single-column HTML preprint manuscript
    └── <run-id>-preprint.pdf                            # Compiled PDF summary edition
```

---

## 4. Autonomous `Super-BLAST` Hygiene & NCBI Error Reporting

During Stage 1, `bio-suite-supervisor` automatically detects, repairs, or quarantines 6 classes of sequence and repository anomalies and records them in `audit/super_blast_audit_and_anomalies.json` and `Table 2` of `manuscript/main.tex` (note: in the `anopheles-3way` demo preset in `src/supervisor.rs`, `PQ585799.1`, `KU551299.1`, `MT890412.1`, and `KP891004.1` are used as synthetic anomaly-test fixtures—on live NCBI Entrez, `PQ585799.1` is *Toxorhynchites rutilus*, `KU551299.1` is *Ophiostoma kryptum*, `MT890412.1` is *Brachycephalus alipioi*, and `KP891004.1` is an uncultured fungus clone; additionally, two cohort-table labels in `src/supervisor.rs` differ from live NCBI Entrez taxonomy: `MK628480.1` is *Anopheles arabiensis* [611 bp COI fragment, labeled *An. coluzzii* in code] and `NC_028215.1` is *Anopheles coluzzii* [labeled *An. merus* in code]; finally, live NCBI RefSeq status checks show `NC_002084.1` [*An. gambiae*], `NC_028215.1` [*An. coluzzii*], `NC_038158.1` [*An. funestus*], and `NC_024740.1` [*An. cruzii*] currently carry `suppressed` RefSeq status):

1. **`NCBI_SAMPLE_SWAP_DISCOVERED`** *(demonstrated in `anopheles-3way` using `PQ585799.1` as a sample-swap fixture)*:
   - Quarantined into the `0🧪` Uncleaned Negative Control cohort and automatically replaced from the reserve pool (`NC_024740.1`, *Anopheles cruzii*).
2. **`NCBI_NUMT_NUCLEAR_PSEUDOGENE`** *(demonstrated in `anopheles-3way` using `KU551299.1` as a premature-stop-codon fixture)*:
   - Quarantined into the `0🧪` control track and replaced with verified full-length mitogenome `NC_028223.1` (*Anopheles stephensi*).
3. **`NCBI_REVERSE_COMPLEMENT_STRAND`** *(demonstrated in `anopheles-3way` on `NC_028215.1` simulated in antisense orientation; live NCBI Entrez taxon is *Anopheles coluzzii*, labeled *An. merus* in `src/supervisor.rs`)*:
   - Automatically reverse-complemented in Rust prior to MSA and retained in the `5🧪` cleaned cohort.
4. **`TERMINAL_ADNA_DEAMINATION_MASKED`** *(demonstrated in `anopheles-3way` using `MT890412.1` as a `5' C->T` / `3' G->A` deamination fixture)*:
   - Terminal deaminated bases are masked to `N` in the `5🧪` cohort while preserving the unmasked sequence in the `0🧪` control track.
5. **`NCBI_PARTIAL_FRAGMENT_REJECTED`** *(demonstrated in `anopheles-3way` using `KP891004.1` as a truncated-fragment fixture)*:
   - Excluded prior to alignment to prevent gap-inflation artifacts; replaced with full-length `NC_028217.1` (*Anopheles epiroticus*).
6. **`NCBI_ENTREZ_FETCH_FAILURE_OR_WITHDRAWN`** *(e.g., `NC_999999.1` unreachable/invalid ID)*:
   - Logged as an NCBI fetch failure and automatically substituted with verified replacement `NC_028218.1` (*Anopheles maculatus*).

---

## 5. Multi-Tab `bio-suite` Desktop GUI Controls

The interactive Win32 GUI (`bio-suite-supervisor gui`) provides 5 specialized workspaces switchable via mouse clicks or `F1`–`F5` keyboard shortcuts (plus mouse scroll-wheel panning on `[F3]` AlignView):

- **`F1` — `[F1] Runs & Visual Pipeline DAG`**:
  - Left sidebar lists all managed runs in `runs/<run-id>/` and one-click Preset launchers (`Anopheles 3-Way`, `Coleoptera Preprint 18`, `Rhinolophus Preprint 21`, `Bemisia Preprint 23`).
  - Top bar switches Model Design (`3-Way Multifactorial`, `Bio-Suite Native Only`, `External Reference Only (*)`, `Tool Calibration`) and clicks **`>> EXECUTE AUTONOMOUS RUN <<`**.
  - Center canvas displays the 7-stage **Visual Pipeline DAG** with live stage runtimes, loop-back iteration logs, and the 9-column **RuNIE Pipeline Comparison Table**.
- **`F2` — `[F2] Super-BLAST & NCBI Audit`**:
  - Displays all flagged NCBI repository anomalies with color-coded severity badges and the full per-accession `Super-BLAST` audit table (`Tip`, `Accession`, `Taxon`, `Source`, `Length`, `Strand`, `Identity %`, `Internal Stops`, `Verdict`).
- **`F3` — `[F3] AlignView Sequence Inspector`**:
  - Full interactive multiple sequence alignment viewer with global locus minimap, Shannon entropy/conservation sparkline, clickable feature annotation sidebar, nucleotide/translated amino-acid color schemes, and mouse scroll-wheel horizontal panning. (For standalone keyboard-driven zoom/palette inspection, use `runie-align-view`).
- **`F4` — `[F4] Chronogram, MSC & Skylines`**:
  - Left panel renders the calibrated divergence chronogram (`Node I`–`Node VII` with two-sided 95% HPD bars) and Quartet MSC `[q1/q2/q3]` support bars.
  - Right panel renders the 3 single-species NCBI `PopSet` Bayesian coalescent skyline trajectories (`RuNIE MCMC` solid curve, `*BEAST v1.10.4` dashed amber curve, and `0🧪 Uncleaned Control` crimson dotted curve).
- **`F5` — `[F5] 3D PDB & Figure/Preprint Studio`**:
  - Left panel provides an interactive 3D Cα wireframe & catalytic pocket viewer for the run's `.pdb` homology model.
  - Right panel lists all 6 locked-format figures (`.svg` + `.png` + `.pdf`) and `manuscript/main.tex` with a one-click **`>> REGENERATE ALL 6 FIGURES & PREPRINT <<`** button.
