# `bio-suite-supervisor` & Multi-Tab `bio-suite` Desktop GUI Guide

`bio-suite-supervisor` (`src/supervisor.rs`, `src/figure_generator.rs`, `src/bio_suite_gui.rs`) is an autonomous, multi-threaded Rust orchestrator for **CRV Science Labs (`bio-suite` & `RuNIE`)**. It eliminates ad-hoc JavaScript/Python scripts by managing the entire bioinformatics study lifecycle inside compiled Rust:

1. **Flexible Cohort Ingestion**: Accepts NCBI Entrez accession numbers (`--accessions`), local FASTA files (`--fasta`), or built-in gold-standard study presets (`--preset anopheles-3way`, `coleoptera-function-18`, `rhinolophus-asia-21`) combined with a factorial **Model Design** (`3way-multifactorial`, `biosuite-only`, `external-only`, `tool-calibration`).
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

# 3. Run with custom NCBI Accession numbers and a custom Run ID + Title
.\target\release\bio-suite-supervisor.exe run `
  --run-id runie-custom-culicidae-audit `
  --title "Autonomous 3-Way Multifactorial Audit of Custom Culicidae Mitogenomes" `
  --design 3way-multifactorial `
  --accessions NC_002084.1,MK628480.1,NC_028212.1,NC_028215.1,NC_038158.1,PQ585799.1,KU551299.1 `
  --threads 16 `
  --bootstrap 1000

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
| `run` / `list` / `inspect` / `figures` / `gui` | Primary supervisor mode | `run` |
| `--run-id <id>` | Unique identifier for the managed run directory (`runs/<run-id>/`) | Derived from preset or `runie-3way-multifactorial-anopheles-mosquitos` |
| `--title <title>` | Study title rendered in GUI headers, figures, and `manuscript/main.tex` | Preset title |
| `--preset <name>` | Built-in dataset preset: `anopheles-3way`, `coleoptera-function-18`, `rhinolophus-asia-21` | `anopheles-3way` |
| `--design <mode>` | Model design: `3way-multifactorial`, `biosuite-only`, `external-only`, `tool-calibration` | `3way-multifactorial` |
| `--accessions <A,B,C>` | Comma-separated NCBI Entrez accession list | Empty (uses preset) |
| `--fasta <path>` | Path to an input FASTA sequence file | `None` |
| `--runs-dir <dir>` | Root directory storing managed run units | `runs` |
| `--cache-dir <dir>` | Local SHA-256 verified NCBI sequence cache directory | `cache` |
| `--threads <N>` | Number of Rayon worker threads (`0` = all logical CPU cores) | `0` (all cores) |
| `--bootstrap <B>` | Ultrafast bootstrap replicates for ML tree inference | `1000` |
| `--no-loopback` | Disable autonomous replacement pull during `Super-BLAST` hygiene | Loop-back enabled |
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

During Stage 1, `bio-suite-supervisor` automatically detects, repairs, or quarantines 6 classes of sequence and repository anomalies and records them in `audit/super_blast_audit_and_anomalies.json` and `Table 2` of `manuscript/main.tex`:

1. **`NCBI_SAMPLE_SWAP_DISCOVERED`** *(e.g., `PQ585799.1` labeled as Anopheles cruzii but matching Culex quinquefasciatus)*:
   - Quarantined into the `0🧪` Uncleaned Negative Control cohort and automatically replaced from the reserve pool (`NC_024740.1`).
2. **`NCBI_NUMT_NUCLEAR_PSEUDOGENE`** *(e.g., `KU551299.1` carrying internal frame-disrupting premature stop codons)*:
   - Quarantined into the `0🧪` control track and replaced with verified full-length mitogenome `NC_028223.1`.
3. **`NCBI_REVERSE_COMPLEMENT_STRAND`** *(e.g., `NC_028215.1` deposited in antisense orientation)*:
   - Automatically reverse-complemented in Rust prior to MSA and retained in the `5🧪` cleaned cohort.
4. **`TERMINAL_ADNA_DEAMINATION_MASKED`** *(e.g., `MT890412.1` carrying `5' C->T` and `3' G->A` post-mortem deamination overhangs)*:
   - Terminal deaminated bases are masked to `N` in the `5🧪` cohort while preserving the unmasked sequence in the `0🧪` control track.
5. **`NCBI_PARTIAL_FRAGMENT_REJECTED`** *(e.g., `KP891004.1` `412-bp` COI barcoding fragment mixed into a full-organelle cohort)*:
   - Excluded prior to alignment to prevent gap-inflation artifacts; replaced with full-length `NC_028217.1`.
6. **`NCBI_ENTREZ_FETCH_FAILURE_OR_WITHDRAWN`** *(e.g., `NC_999999.1` withdrawn/unreachable record)*:
   - Logged as an NCBI fetch failure and automatically substituted with verified replacement `NC_028218.1`.

---

## 5. Multi-Tab `bio-suite` Desktop GUI Controls

The interactive Win32 GUI (`bio-suite-supervisor gui`) provides 5 specialized workspaces controllable via mouse clicks or keyboard shortcuts:

- **`F1` — `[F1] Runs & Visual Pipeline DAG`**:
  - Left sidebar lists all managed runs in `runs/<run-id>/` and one-click Preset launchers (`Anopheles 3-Way`, `Coleoptera Preprint 18`, `Rhinolophus Preprint 21`).
  - Top bar switches Model Design (`3-Way Multifactorial`, `Bio-Suite Native Only`, `External Reference Only (*)`, `Tool Calibration`) and triggers **`>> EXECUTE AUTONOMOUS RUN <<`** (`R` key).
  - Center canvas displays the 7-stage **Visual Pipeline DAG** with live stage runtimes, loop-back iteration logs, and the 9-column **RuNIE Pipeline Comparison Table**.
- **`F2` — `[F2] Super-BLAST & NCBI Audit`**:
  - Displays all flagged NCBI repository anomalies with color-coded severity badges and the full per-accession `Super-BLAST` audit table (`Tip`, `Accession`, `Taxon`, `Source`, `Length`, `Strand`, `Identity %`, `Internal Stops`, `Verdict`). Use `Up`/`Down` arrows to scroll.
- **`F3` — `[F3] AlignView Sequence Inspector`**:
  - Full interactive multiple sequence alignment viewer with global locus minimap, Shannon entropy/conservation sparkline, clickable feature annotation sidebar, nucleotide/translated amino-acid color schemes (`1`–`6`), zoom (`+`/`-`), bird's-eye view (`0`), and focal window jump (`F`).
- **`F4` — `[F4] Chronogram, MSC & Skylines`**:
  - Left panel renders the calibrated divergence chronogram (`Node I`–`Node VII` with two-sided 95% HPD bars) and Quartet MSC `[q1/q2/q3]` support bars.
  - Right panel renders the 3 single-species NCBI `PopSet` Bayesian coalescent skyline trajectories (`RuNIE MCMC` solid curve, `*BEAST v1.10.4` dashed amber curve with `-6.5 px` visual offset, and `0🧪 Uncleaned Control` crimson dotted curve).
- **`F5` — `[F5] 3D PDB & Figure/Preprint Studio`**:
  - Left panel provides an interactive 3D Cα wireframe & catalytic pocket viewer for the run's `.pdb` homology model (`J`/`L` or `Left`/`Right` to yaw, `I`/`K` to pitch, `0` to reset view).
  - Right panel lists all 6 locked-format figures (`.svg` + `.png` + `.pdf`) and `manuscript/main.tex` with a one-click **`>> REGENERATE ALL 6 FIGURES & PREPRINT <<`** button (`G` or `F12`).
