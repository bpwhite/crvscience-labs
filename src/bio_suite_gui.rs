//! `bio_suite_gui` — Unified Multi-Tab Bioinformatics Suite Desktop GUI (`bio-suite gui` / `bio-suite supervisor --gui`)
//!
//! Upgrades the CRV Science Labs native Win32 framebuffer GUI into a complete 5-workspace
//! Bioinformatics Suite controllable identically via GUI or console:
//! - **Tab 1 (`[F1] Runs & Visual Pipeline Builder`)**:
//!   Manage "Runs" as first-class units (`runs/<run_id>/`), visually configure datasets, presets,
//!   model designs (`3-Way Multifactorial`, `Bio-Suite Only`, `External Only`, `Tool Calibration`),
//!   and pipeline stages on an interactive DAG canvas, and launch multi-threaded autonomous runs.
//! - **Tab 2 (`[F2] Super-BLAST & Data Hygiene Audit`)**:
//!   Inspect all downloaded sequences, completeness checks, strand orientations, flagged NCBI
//!   sample swaps (`PQ585799.1`), `NUMT` pseudogenes, and autonomous loop-back replacements.
//! - **Tab 3 (`[F3] AlignView Sequence Inspector`)**:
//!   Full interactive AliView-modeled alignment viewer (`DNA`, `AA-3bp`, `AA-Compact`, `Dual DNA+AA`,
//!   `RY3`, `Dayhoff-6`, Minimap, Feature Sidebar) bound directly to the active Run's alignment.
//! - **Tab 4 (`[F4] Trees, Chronogram & Skylines`)**:
//!   Interactive visual inspector for the solved calibrated chronogram (`Node I`–`Node VII`),
//!   Quartet MSC `[q1/q2/q3]` bars, and coupled 3-PopSet Bayesian coalescent skyline curves.
//! - **Tab 5 (`[F5] 3D PDB & Locked Figure / Preprint Studio`)**:
//!   Interactive 3D `.pdb` C-alpha & catalytic pocket viewer (with live yaw/pitch rotation) and
//!   one-click generator/inspector for all 6 locked-format figures (`.svg` + `.png` + `.pdf`)
//!   and the arXiv-ready LaTeX preprint (`manuscript/main.tex`).

use crate::align_viewer::{
    AlignmentViewerState, AlignmentWorkspace, GeneticCodeTable, SoftwareFrameBuffer, UiHitRect,
};
use crate::fetch::PoliteCacheManager;
use crate::figure_generator::generate_all_locked_figures_and_preprint;
use crate::supervisor::{
    list_existing_supervisor_runs, run_autonomous_supervisor, ModelDesignMode,
    SupervisorRunManifest, SupervisorRunSpec,
};
use std::fs;
use std::path::{Path, PathBuf};

/// Top-level workspace tab in the `bio-suite` GUI.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BioSuiteGuiTab {
    RunsAndPipelineBuilder,
    SuperBlastDataAudit,
    AlignViewInspector,
    TreesAndSkylines,
    PdbAndFigureStudio,
}

impl BioSuiteGuiTab {
    pub fn label(&self) -> &'static str {
        match self {
            BioSuiteGuiTab::RunsAndPipelineBuilder => "[F1] Runs & Visual Pipeline DAG",
            BioSuiteGuiTab::SuperBlastDataAudit => "[F2] Super-BLAST & NCBI Audit",
            BioSuiteGuiTab::AlignViewInspector => "[F3] AlignView Sequence Inspector",
            BioSuiteGuiTab::TreesAndSkylines => "[F4] Chronogram, MSC & Skylines",
            BioSuiteGuiTab::PdbAndFigureStudio => "[F5] 3D PDB & Figure/Preprint Studio",
        }
    }
}

// GUI Action IDs for top tabs, pipeline builder controls, and studio actions
pub const GUI_ACT_TAB_1: u16 = 3001;
pub const GUI_ACT_TAB_2: u16 = 3002;
pub const GUI_ACT_TAB_3: u16 = 3003;
pub const GUI_ACT_TAB_4: u16 = 3004;
pub const GUI_ACT_TAB_5: u16 = 3005;

pub const GUI_ACT_DESIGN_3WAY: u16 = 3101;
pub const GUI_ACT_DESIGN_BIOSUITE: u16 = 3102;
pub const GUI_ACT_DESIGN_EXTERNAL: u16 = 3103;
pub const GUI_ACT_DESIGN_CALIB: u16 = 3104;

pub const GUI_ACT_PRESET_ANOPHELES: u16 = 3111;
pub const GUI_ACT_PRESET_COLEOPTERA18: u16 = 3112;
pub const GUI_ACT_PRESET_RHINOLOPHUS21: u16 = 3113;
pub const GUI_ACT_PRESET_BEMISIA23: u16 = 3114;

pub const GUI_ACT_TOGGLE_LOOPBACK: u16 = 3121;
pub const GUI_ACT_TOGGLE_RY3: u16 = 3122;
pub const GUI_ACT_EXECUTE_RUN: u16 = 3130;
pub const GUI_ACT_REGEN_FIGURES_TEX: u16 = 3131;

pub const GUI_ACT_PDB_YAW_LEFT: u16 = 3201;
pub const GUI_ACT_PDB_YAW_RIGHT: u16 = 3202;
pub const GUI_ACT_PDB_PITCH_UP: u16 = 3203;
pub const GUI_ACT_PDB_PITCH_DOWN: u16 = 3204;
pub const GUI_ACT_PDB_RESET_VIEW: u16 = 3205;

pub const GUI_ACT_SELECT_RUN_BASE: u16 = 3500;

/// Complete interactive state for the Multi-Tab `bio-suite` GUI.
pub struct BioSuiteGuiState {
    pub active_tab: BioSuiteGuiTab,
    pub spec: SupervisorRunSpec,
    pub active_manifest: SupervisorRunManifest,
    pub available_runs: Vec<SupervisorRunManifest>,
    pub align_viewer: AlignmentViewerState,
    pub runs_root_dir: PathBuf,
    pub cache_dir: PathBuf,
    pub pdb_yaw_rad: f64,
    pub pdb_pitch_rad: f64,
    pub selected_figure_idx: usize,
    pub status_banner: String,
}

impl BioSuiteGuiState {
    pub fn new(
        manifest: SupervisorRunManifest,
        runs_root_dir: PathBuf,
        cache_dir: PathBuf,
    ) -> Self {
        let aln_path = runs_root_dir
            .join(&manifest.run_id)
            .join(&manifest.aligned_fasta_rel_path);
        let ws = AlignmentWorkspace::from_file(
            &aln_path,
            None,
            GeneticCodeTable::InvertebrateMito5,
            1,
            true,
        )
        .unwrap_or_else(|_| {
            let cache = PoliteCacheManager::new(&cache_dir);
            AlignmentWorkspace::from_builtin_cohort(
                "coi-barcodes",
                &cache,
                GeneticCodeTable::InvertebrateMito5,
            )
            .expect("Fallback cohort load")
        });

        let align_viewer = AlignmentViewerState::new(ws);
        let available_runs = list_existing_supervisor_runs(&runs_root_dir);
        let mut spec = SupervisorRunSpec::default();
        spec.run_id = manifest.run_id.clone();
        spec.model_design = manifest.model_design;

        let status_banner = format!(
            "Active Run: {} | Design: {} | {} Cleaned Taxa (L = {} bp) | {} NCBI Anomalies Auto-Resolved | {} Figures + .tex Ready",
            manifest.run_id,
            manifest.model_design.code(),
            manifest.final_cleaned_cohort_count,
            manifest.native_mafft_aligned_cols_l,
            manifest.ncbi_flagged_anomalies.len(),
            manifest.figures.len()
        );

        Self {
            active_tab: BioSuiteGuiTab::RunsAndPipelineBuilder,
            spec,
            active_manifest: manifest,
            available_runs,
            align_viewer,
            runs_root_dir,
            cache_dir,
            pdb_yaw_rad: 0.55,
            pdb_pitch_rad: 0.35,
            selected_figure_idx: 0,
            status_banner,
        }
    }

    pub fn reload_align_viewer_from_active_run(&mut self) {
        let aln_path = self
            .runs_root_dir
            .join(&self.active_manifest.run_id)
            .join(&self.active_manifest.aligned_fasta_rel_path);
        if let Ok(ws) = AlignmentWorkspace::from_file(
            &aln_path,
            None,
            GeneticCodeTable::InvertebrateMito5,
            1,
            true,
        ) {
            self.align_viewer = AlignmentViewerState::new(ws);
        }
    }

    pub fn handle_gui_action(&mut self, action_id: u16, canvas_w: i32) {
        if action_id >= GUI_ACT_SELECT_RUN_BASE && action_id < GUI_ACT_SELECT_RUN_BASE + 100 {
            let idx = (action_id - GUI_ACT_SELECT_RUN_BASE) as usize;
            if let Some(r) = self.available_runs.get(idx).cloned() {
                self.spec.run_id = r.run_id.clone();
                self.spec.model_design = r.model_design;
                self.active_manifest = r;
                self.reload_align_viewer_from_active_run();
                self.status_banner = format!(
                    "Loaded Managed Run '{}' (L = {} aligned columns, {} figures).",
                    self.active_manifest.run_id,
                    self.active_manifest.native_mafft_aligned_cols_l,
                    self.active_manifest.figures.len()
                );
            }
            return;
        }

        match action_id {
            GUI_ACT_TAB_1 => {
                self.active_tab = BioSuiteGuiTab::RunsAndPipelineBuilder;
                self.status_banner =
                    "Switched to [F1] Runs & Visual Pipeline DAG Builder".to_string();
            }
            GUI_ACT_TAB_2 => {
                self.active_tab = BioSuiteGuiTab::SuperBlastDataAudit;
                self.status_banner = format!(
                    "Switched to [F2] Super-BLAST Data Hygiene Audit ({} sequences, {} NCBI anomalies flagged)",
                    self.active_manifest.sequences.len(),
                    self.active_manifest.ncbi_flagged_anomalies.len()
                );
            }
            GUI_ACT_TAB_3 => {
                self.active_tab = BioSuiteGuiTab::AlignViewInspector;
                self.status_banner = format!(
                    "Switched to [F3] AlignView Sequence Inspector (L = {} bp; Keys 1-6 change alphabet)",
                    self.align_viewer.workspace.alignment_length
                );
            }
            GUI_ACT_TAB_4 => {
                self.active_tab = BioSuiteGuiTab::TreesAndSkylines;
                self.status_banner =
                    "Switched to [F4] Calibrated Chronogram, Quartet Coalescent & Bayesian Skylines"
                        .to_string();
            }
            GUI_ACT_TAB_5 => {
                self.active_tab = BioSuiteGuiTab::PdbAndFigureStudio;
                self.status_banner =
                    "Switched to [F5] 3D PDB Structure & Locked Figure / Preprint (.tex) Studio"
                        .to_string();
            }
            GUI_ACT_DESIGN_3WAY => {
                self.spec.model_design = ModelDesignMode::ThreeWayMultifactorial;
                self.status_banner =
                    "Selected Model Design: 3-Way Multifactorial (Native vs. External * vs. Hybrid + 0-Beaker)"
                        .to_string();
            }
            GUI_ACT_DESIGN_BIOSUITE => {
                self.spec.model_design = ModelDesignMode::BioSuiteNativeOnly;
                self.status_banner =
                    "Selected Model Design: 100% Native Rust (bio-suite Only)".to_string();
            }
            GUI_ACT_DESIGN_EXTERNAL => {
                self.spec.model_design = ModelDesignMode::ExternalReferenceOnly;
                self.status_banner =
                    "Selected Model Design: 100% External Reference Binaries (*MAFFT + *IQ-TREE + *BEAST)"
                        .to_string();
            }
            GUI_ACT_DESIGN_CALIB => {
                self.spec.model_design = ModelDesignMode::ToolCalibration;
                self.status_banner =
                    "Selected Model Design: Head-to-Head Native vs. External Tool Calibration"
                        .to_string();
            }
            GUI_ACT_PRESET_ANOPHELES => {
                self.spec.preset_name = Some("anopheles-3way".to_string());
                self.spec.run_id = "runie-3way-multifactorial-anopheles-mosquitos".to_string();
                self.status_banner =
                    "Selected Preset: Runie Run 3-Way Multifactorial Anopheles Mosquitos — Click [▶ EXECUTE AUTONOMOUS RUN]"
                        .to_string();
            }
            GUI_ACT_PRESET_COLEOPTERA18 => {
                self.spec.preset_name = Some("coleoptera-function-18".to_string());
                self.spec.run_id = "runie-3way-coleoptera-function-18".to_string();
                self.status_banner =
                    "Selected Preset: Preprint 18 Gold Standard (Coleoptera 6-Family Functional Phylogenomics & Skyline)"
                        .to_string();
            }
            GUI_ACT_PRESET_RHINOLOPHUS21 => {
                self.spec.preset_name = Some("rhinolophus-asia-21".to_string());
                self.spec.run_id = "runie-3way-rhinolophus-asia-21".to_string();
                self.status_banner =
                    "Selected Preset: Preprint 21 Gold Standard (Asian Rhinolophidae Quartet MSC & Prestin/ACE2)"
                        .to_string();
            }
            GUI_ACT_PRESET_BEMISIA23 => {
                self.spec.preset_name = Some("bemisia-whiteflies-23".to_string());
                self.spec.run_id = "runie-native-bemisia-whiteflies-23".to_string();
                self.spec.model_design = ModelDesignMode::BioSuiteNativeOnly;
                self.spec.mcmc_generations = 25_000;
                self.status_banner =
                    "Selected Preset: Preprint 23 Multi-Compartment (Bemisia / Aleyrodidae 544 Verified Records, Native Bio-Suite)"
                        .to_string();
            }
            GUI_ACT_TOGGLE_LOOPBACK => {
                self.spec.enable_autonomous_loopback = !self.spec.enable_autonomous_loopback;
                self.status_banner = format!(
                    "Autonomous Super-BLAST Loop-Back Self-Healing: {}",
                    if self.spec.enable_autonomous_loopback {
                        "ENABLED (Auto-Repairs & Replaces Failed Accessions)"
                    } else {
                        "DISABLED"
                    }
                );
            }
            GUI_ACT_TOGGLE_RY3 => {
                self.spec.enable_ry3_and_dayhoff6 = !self.spec.enable_ry3_and_dayhoff6;
                self.status_banner = format!(
                    "RuNIE RY3 & Dayhoff-6 Compositional Recoding: {}",
                    if self.spec.enable_ry3_and_dayhoff6 {
                        "ENABLED"
                    } else {
                        "DISABLED"
                    }
                );
            }
            GUI_ACT_EXECUTE_RUN => {
                let cache = PoliteCacheManager::new(&self.cache_dir);
                if let Ok(new_manifest) =
                    run_autonomous_supervisor(&self.spec, &cache, &self.runs_root_dir)
                {
                    self.active_manifest = new_manifest;
                    self.available_runs = list_existing_supervisor_runs(&self.runs_root_dir);
                    self.reload_align_viewer_from_active_run();
                    self.status_banner = format!(
                        "Completed Autonomous Run '{}' in {} ms ({} cleaned taxa, {} NCBI anomalies resolved, 6 figures + main.tex generated)!",
                        self.active_manifest.run_id,
                        self.active_manifest.total_wall_clock_ms,
                        self.active_manifest.final_cleaned_cohort_count,
                        self.active_manifest.ncbi_flagged_anomalies.len()
                    );
                }
            }
            GUI_ACT_REGEN_FIGURES_TEX => {
                let run_dir = self.runs_root_dir.join(&self.active_manifest.run_id);
                let rep_path = run_dir.join("multi-compartment-report.json");
                let regen_ok = if rep_path.is_file() {
                    if let Ok(rep_txt) = fs::read_to_string(&rep_path) {
                        if let Ok(rep) = serde_json::from_str::<
                            crate::supervisor_study::MultiCompartmentStudyReport,
                        >(&rep_txt)
                        {
                            let ws_root = if Path::new("../labs").is_dir() {
                                PathBuf::from("..")
                            } else {
                                PathBuf::from(".")
                            };
                            let labs_fig_dir = ws_root.join("labs").join("figures");
                            crate::supervisor_figures::generate_bemisia_locked_figures(
                                &mut self.active_manifest,
                                &rep,
                                &run_dir,
                                Some(labs_fig_dir.as_path()),
                            )
                            .and_then(|_| {
                                crate::supervisor_manuscript::export_bemisia_manuscript_and_staged_bundle(
                                    &rep,
                                    &self.active_manifest,
                                    &ws_root,
                                    &run_dir,
                                )
                                .map(|_| ())
                            })
                            .is_ok()
                        } else {
                            false
                        }
                    } else {
                        false
                    }
                } else {
                    generate_all_locked_figures_and_preprint(&mut self.active_manifest, &run_dir)
                        .is_ok()
                };
                if regen_ok {
                    self.status_banner = format!(
                        "Regenerated all 6 locked-format figures (.svg + .png + .pdf) and manuscript/main.tex in runs/{}/",
                        self.active_manifest.run_id
                    );
                }
            }
            GUI_ACT_PDB_YAW_LEFT => {
                self.pdb_yaw_rad -= 0.22;
                self.status_banner = format!(
                    "Rotated 3D PDB Structure: Yaw = {:.2} rad, Pitch = {:.2} rad",
                    self.pdb_yaw_rad, self.pdb_pitch_rad
                );
            }
            GUI_ACT_PDB_YAW_RIGHT => {
                self.pdb_yaw_rad += 0.22;
                self.status_banner = format!(
                    "Rotated 3D PDB Structure: Yaw = {:.2} rad, Pitch = {:.2} rad",
                    self.pdb_yaw_rad, self.pdb_pitch_rad
                );
            }
            GUI_ACT_PDB_PITCH_UP => {
                self.pdb_pitch_rad = (self.pdb_pitch_rad + 0.18).clamp(-1.4, 1.4);
                self.status_banner = format!(
                    "Rotated 3D PDB Structure: Yaw = {:.2} rad, Pitch = {:.2} rad",
                    self.pdb_yaw_rad, self.pdb_pitch_rad
                );
            }
            GUI_ACT_PDB_PITCH_DOWN => {
                self.pdb_pitch_rad = (self.pdb_pitch_rad - 0.18).clamp(-1.4, 1.4);
                self.status_banner = format!(
                    "Rotated 3D PDB Structure: Yaw = {:.2} rad, Pitch = {:.2} rad",
                    self.pdb_yaw_rad, self.pdb_pitch_rad
                );
            }
            GUI_ACT_PDB_RESET_VIEW => {
                self.pdb_yaw_rad = 0.55;
                self.pdb_pitch_rad = 0.35;
                self.status_banner = "Reset 3D PDB Structure View Orientation".to_string();
            }
            _ => {
                if self.active_tab == BioSuiteGuiTab::AlignViewInspector {
                    self.align_viewer.execute_action(action_id, canvas_w);
                    self.status_banner = self.align_viewer.status_message.clone();
                }
            }
        }
    }
}
impl BioSuiteGuiState {
    fn draw_pill_button(
        fb: &mut SoftwareFrameBuffer,
        hits: &mut Vec<UiHitRect>,
        x: i32,
        y: i32,
        w: i32,
        h: i32,
        label: &str,
        action_id: u16,
        active: bool,
        accent_rgb: (u8, u8, u8),
    ) {
        let bg = if active {
            accent_rgb
        } else {
            (30, 41, 59)
        };
        let fg = if active {
            (255, 255, 255)
        } else {
            (203, 213, 225)
        };
        fb.fill_rect(x, y, w, h, bg);
        fb.stroke_rect(x, y, w, h, accent_rgb);
        fb.draw_text(x + 8, y + (h - 8) / 2, label, fg, 1, Some(w - 12));
        hits.push(UiHitRect {
            x,
            y,
            w,
            h,
            action_id,
        });
    }

    /// Renders the complete multi-tab `bio-suite` GUI into `fb` and returns all clickable `UiHitRect`s.
    pub fn render_bio_suite_gui_frame(&self, fb: &mut SoftwareFrameBuffer) -> Vec<UiHitRect> {
        let mut hits = Vec::with_capacity(96);
        let w = fb.width as i32;
        let h = fb.height as i32;

        // If on Tab 3 (AlignView Inspector), let AlignView render first, then overlay our top Suite Tab Bar!
        if self.active_tab == BioSuiteGuiTab::AlignViewInspector {
            let mut av_hits = self.align_viewer.render_interactive_frame(fb);
            hits.append(&mut av_hits);
            self.render_top_navigation_bar(fb, &mut hits);
            return hits;
        }

        fb.clear((11, 15, 25));
        self.render_top_navigation_bar(fb, &mut hits);

        match self.active_tab {
            BioSuiteGuiTab::RunsAndPipelineBuilder => {
                self.render_tab1_runs_and_pipeline_dag(fb, &mut hits);
            }
            BioSuiteGuiTab::SuperBlastDataAudit => {
                self.render_tab2_super_blast_audit(fb, &mut hits);
            }
            BioSuiteGuiTab::AlignViewInspector => {}
            BioSuiteGuiTab::TreesAndSkylines => {
                self.render_tab4_trees_and_skylines(fb, &mut hits);
            }
            BioSuiteGuiTab::PdbAndFigureStudio => {
                self.render_tab5_pdb_and_figure_studio(fb, &mut hits);
            }
        }

        // Bottom Status Bar
        fb.fill_rect(0, h - 28, w, 28, (15, 23, 42));
        fb.fill_rect(0, h - 28, w, 1, (51, 65, 85));
        fb.draw_text(
            14,
            h - 18,
            &self.status_banner,
            (56, 189, 248),
            1,
            Some(w - 28),
        );

        hits
    }

    fn render_top_navigation_bar(&self, fb: &mut SoftwareFrameBuffer, hits: &mut Vec<UiHitRect>) {
        let w = fb.width as i32;
        fb.fill_rect(0, 0, w, 38, (15, 23, 42));
        fb.fill_rect(0, 37, w, 1, (56, 189, 248));

        fb.draw_text(
            12,
            14,
            "CRV BIO-SUITE SUPERVISOR",
            (56, 189, 248),
            1,
            Some(165),
        );

        let tabs: &[(BioSuiteGuiTab, u16, (u8, u8, u8))] = &[
            (
                BioSuiteGuiTab::RunsAndPipelineBuilder,
                GUI_ACT_TAB_1,
                (2, 132, 199),
            ),
            (
                BioSuiteGuiTab::SuperBlastDataAudit,
                GUI_ACT_TAB_2,
                (16, 185, 129),
            ),
            (
                BioSuiteGuiTab::AlignViewInspector,
                GUI_ACT_TAB_3,
                (99, 102, 241),
            ),
            (
                BioSuiteGuiTab::TreesAndSkylines,
                GUI_ACT_TAB_4,
                (217, 119, 6),
            ),
            (
                BioSuiteGuiTab::PdbAndFigureStudio,
                GUI_ACT_TAB_5,
                (225, 29, 72),
            ),
        ];

        let mut tx = 180;
        for &(tab, act, col) in tabs {
            let bw = 236;
            Self::draw_pill_button(
                fb,
                hits,
                tx,
                6,
                bw,
                26,
                tab.label(),
                act,
                self.active_tab == tab,
                col,
            );
            tx += bw + 8;
        }
    }

    fn render_tab1_runs_and_pipeline_dag(
        &self,
        fb: &mut SoftwareFrameBuffer,
        hits: &mut Vec<UiHitRect>,
    ) {
        let w = fb.width as i32;
        let h = fb.height as i32;

        // Left Sidebar: Managed Runs Unit Explorer
        let side_w = 320;
        fb.fill_rect(0, 38, side_w, h - 66, (15, 23, 42));
        fb.fill_rect(side_w - 1, 38, 1, h - 66, (51, 65, 85));
        fb.draw_text(
            14,
            52,
            "MANAGED RUN UNITS (runs/<run-id>/)",
            (56, 189, 248),
            1,
            Some(side_w - 24),
        );

        let mut ry = 76;
        for (idx, r) in self.available_runs.iter().take(6).enumerate() {
            let is_cur = r.run_id == self.active_manifest.run_id;
            let bg = if is_cur { (30, 58, 138) } else { (23, 33, 53) };
            let border = if is_cur {
                (56, 189, 248)
            } else {
                (51, 65, 85)
            };
            fb.fill_rect(10, ry, side_w - 20, 54, bg);
            fb.stroke_rect(10, ry, side_w - 20, 54, border);
            fb.draw_text(18, ry + 8, &r.run_id, (255, 255, 255), 1, Some(side_w - 36));
            fb.draw_text(
                18,
                ry + 24,
                &format!(
                    "Design: {} | L={} bp",
                    r.model_design.code(),
                    r.native_mafft_aligned_cols_l
                ),
                (148, 163, 184),
                1,
                Some(side_w - 36),
            );
            fb.draw_text(
                18,
                ry + 38,
                &format!(
                    "Cleaned: {} | Anomalies: {} | Figs: {}",
                    r.final_cleaned_cohort_count,
                    r.ncbi_flagged_anomalies.len(),
                    r.figures.len()
                ),
                (52, 211, 153),
                1,
                Some(side_w - 36),
            );
            hits.push(UiHitRect {
                x: 10,
                y: ry,
                w: side_w - 20,
                h: 54,
                action_id: GUI_ACT_SELECT_RUN_BASE + idx as u16,
            });
            ry += 62;
        }

        // Preset Quick-Loaders in Left Sidebar
        ry += 10;
        fb.draw_text(
            14,
            ry,
            "BUILD NEW RUN FROM PRESET:",
            (251, 191, 36),
            1,
            Some(side_w - 24),
        );
        ry += 18;
        Self::draw_pill_button(
            fb,
            hits,
            12,
            ry,
            side_w - 24,
            26,
            "Preset: Anopheles 3-Way Multifactorial",
            GUI_ACT_PRESET_ANOPHELES,
            self.spec.run_id.contains("anopheles"),
            (2, 132, 199),
        );
        ry += 32;
        Self::draw_pill_button(
            fb,
            hits,
            12,
            ry,
            side_w - 24,
            26,
            "Preset: Coleoptera Preprint 18 Standard",
            GUI_ACT_PRESET_COLEOPTERA18,
            self.spec.run_id.contains("coleoptera"),
            (13, 148, 136),
        );
        ry += 32;
        Self::draw_pill_button(
            fb,
            hits,
            12,
            ry,
            side_w - 24,
            26,
            "Preset: Rhinolophus Preprint 21 Standard",
            GUI_ACT_PRESET_RHINOLOPHUS21,
            self.spec.run_id.contains("rhinolophus"),
            (124, 58, 237),
        );
        ry += 32;
        Self::draw_pill_button(
            fb,
            hits,
            12,
            ry,
            side_w - 24,
            26,
            "Preset: Bemisia Preprint 23 Native (544)",
            GUI_ACT_PRESET_BEMISIA23,
            self.spec.run_id.contains("bemisia"),
            (16, 185, 129),
        );

        // Right Workspace: Model Design Selector + Visual Pipeline DAG Canvas + 9-Col RuNIE Table
        let mx0 = side_w + 16;
        let mw = w - mx0 - 16;

        fb.fill_rect(mx0, 48, mw, 78, (15, 23, 42));
        fb.stroke_rect(mx0, 48, mw, 78, (51, 65, 85));
        fb.draw_text(
            mx0 + 12,
            58,
            &format!(
                "ACTIVE RUN UNIT: {}  |  Worker Threads: {} (Rayon Parallel)",
                self.active_manifest.run_id, self.active_manifest.rayon_worker_threads
            ),
            (255, 255, 255),
            1,
            Some(mw - 24),
        );

        // Model Design Buttons
        let designs: &[(&str, ModelDesignMode, u16, (u8, u8, u8))] = &[
            (
                "3-Way Multifactorial",
                ModelDesignMode::ThreeWayMultifactorial,
                GUI_ACT_DESIGN_3WAY,
                (2, 132, 199),
            ),
            (
                "Bio-Suite Native Only",
                ModelDesignMode::BioSuiteNativeOnly,
                GUI_ACT_DESIGN_BIOSUITE,
                (16, 185, 129),
            ),
            (
                "External Reference Only (*)",
                ModelDesignMode::ExternalReferenceOnly,
                GUI_ACT_DESIGN_EXTERNAL,
                (217, 119, 6),
            ),
            (
                "Tool Calibration",
                ModelDesignMode::ToolCalibration,
                GUI_ACT_DESIGN_CALIB,
                (139, 92, 246),
            ),
        ];
        let mut dx = mx0 + 12;
        for &(lbl, mode, act, col) in designs {
            Self::draw_pill_button(
                fb,
                hits,
                dx,
                84,
                182,
                28,
                lbl,
                act,
                self.spec.model_design == mode,
                col,
            );
            dx += 190;
        }

        // Launch Autonomous Run Button
        Self::draw_pill_button(
            fb,
            hits,
            dx + 8,
            82,
            245,
            32,
            ">> EXECUTE AUTONOMOUS RUN <<",
            GUI_ACT_EXECUTE_RUN,
            true,
            (22, 163, 74),
        );

        // Visual DAG Pipeline Graph Box
        let dag_y = 138;
        let dag_h = 210;
        fb.fill_rect(mx0, dag_y, mw, dag_h, (15, 23, 42));
        fb.stroke_rect(mx0, dag_y, mw, dag_h, (56, 189, 248));
        fb.draw_text(
            mx0 + 14,
            dag_y + 12,
            "VISUAL PIPELINE DAG & AUTONOMOUS SELF-HEALING TELEMETRY (95%+ Autonomous Loop-Back Engine)",
            (56, 189, 248),
            1,
            Some(mw - 28),
        );

        let dag_nodes: &[(&str, String, (u8, u8, u8))] = &[
            (
                "1. Fetch & Input",
                format!(
                    "{} Init Recs\nCache + Entrez",
                    self.active_manifest.initial_requested_count
                ),
                (2, 132, 199),
            ),
            (
                "2. Super-BLAST Gate",
                format!(
                    "{} Anomalies\nLoop-Back: OK",
                    self.active_manifest.ncbi_flagged_anomalies.len()
                ),
                (16, 185, 129),
            ),
            (
                "3. Full-Length MSA",
                format!(
                    "L = {} bp\n{} ms",
                    self.active_manifest.native_mafft_aligned_cols_l,
                    self.active_manifest.stage_msa_ms
                ),
                (99, 102, 241),
            ),
            (
                "4. 32-Chain RuNIE",
                format!(
                    "8 Mod x 4 Reg\nRCFV={:.4}",
                    self.active_manifest.rcfv_ry3_cleaned
                ),
                (217, 119, 6),
            ),
            (
                "5. Clock + MSC + dN/dS",
                format!(
                    "7 Nodes + 4 Fam\n{} ms",
                    self.active_manifest.stage_runie_phylo_ms
                ),
                (14, 165, 233),
            ),
            (
                "6. Skyline + 3D PDB",
                format!(
                    "3 PopSets + PDB\n{} ms",
                    self.active_manifest.stage_skyline_pdb_ms
                ),
                (168, 85, 247),
            ),
            (
                "7. Figures + .tex",
                format!(
                    "{} SVG/PNG/PDF\n+ main.tex",
                    self.active_manifest.figures.len()
                ),
                (225, 29, 72),
            ),
        ];

        let node_w = ((mw - 40) / 7).max(110) - 10;
        for (i, (title, sub, col)) in dag_nodes.iter().enumerate() {
            let nx = mx0 + 16 + (i as i32) * (node_w + 12);
            let ny = dag_y + 42;
            fb.fill_rect(nx, ny, node_w, 86, (23, 33, 53));
            fb.stroke_rect(nx, ny, node_w, 86, *col);
            fb.fill_rect(nx, ny, node_w, 22, *col);
            fb.draw_text(nx + 6, ny + 7, title, (255, 255, 255), 1, Some(node_w - 10));
            for (li, line) in sub.lines().enumerate() {
                fb.draw_text(
                    nx + 6,
                    ny + 34 + (li as i32) * 18,
                    line,
                    (226, 232, 240),
                    1,
                    Some(node_w - 10),
                );
            }
            if i + 1 < dag_nodes.len() {
                fb.draw_line(
                    nx + node_w,
                    ny + 43,
                    nx + node_w + 12,
                    ny + 43,
                    2,
                    (56, 189, 248),
                );
            }
        }

        // Loop-Back Summary Banner inside DAG box
        for (li, lb) in self.active_manifest.loop_back_iterations.iter().enumerate() {
            fb.draw_text(
                mx0 + 16,
                dag_y + 144 + (li as i32) * 20,
                &format!("Loop-Back {}: {}", lb.iteration, lb.summary_note),
                (52, 211, 153),
                1,
                Some(mw - 32),
            );
        }

        // 9-Column RuNIE Pipeline Comparison Table Box
        let tbl_y = dag_y + dag_h + 14;
        let tbl_h = h - tbl_y - 38;
        fb.fill_rect(mx0, tbl_y, mw, tbl_h, (15, 23, 42));
        fb.stroke_rect(mx0, tbl_y, mw, tbl_h, (51, 65, 85));
        fb.draw_text(
            mx0 + 14,
            tbl_y + 12,
            "TABLE 3A. RuNIE PIPELINE COMPARISON TABLE (9-COLUMN LOCKED FORMAT + ARTIFACT BUNDLE)",
            (251, 191, 36),
            1,
            Some(mw - 28),
        );

        for (idx, row) in self
            .active_manifest
            .runie_9col_comparison_table
            .iter()
            .enumerate()
        {
            let ry = tbl_y + 38 + (idx as i32) * 62;
            let c_border = match idx {
                0 => (2, 132, 199),
                1 => (13, 148, 136),
                2 => (217, 119, 6),
                _ => (220, 38, 38),
            };
            fb.fill_rect(mx0 + 12, ry, mw - 24, 54, (23, 33, 53));
            fb.stroke_rect(mx0 + 12, ry, mw - 24, 54, c_border);
            fb.draw_text(
                mx0 + 20,
                ry + 8,
                &format!(
                    "{} | {} | Marginal ln L = {:.2} (Delta = {:.2}) | RCFV = {:.4} | Bootstrap = {} | Z = {:.2} | CPU = {}",
                    row.rank_label,
                    row.model_summary,
                    row.marginal_ln_l,
                    row.delta_ln_l_vs_best,
                    row.rcfv,
                    row.runie_bootstrap,
                    row.posterior_z,
                    row.cpu_time_summary
                ),
                (255, 255, 255),
                1,
                Some(mw - 40),
            );
            fb.draw_text(
                mx0 + 20,
                ry + 28,
                &format!("Chain: {}", row.tool_chain),
                (148, 163, 184),
                1,
                Some(mw - 40),
            );
        }
    }

    fn render_tab2_super_blast_audit(
        &self,
        fb: &mut SoftwareFrameBuffer,
        _hits: &mut Vec<UiHitRect>,
    ) {
        let w = fb.width as i32;
        let h = fb.height as i32;

        // Top Panel: Flagged NCBI Repository & Sequence Anomalies
        fb.fill_rect(16, 48, w - 32, 235, (15, 23, 42));
        fb.stroke_rect(16, 48, w - 32, 235, (220, 38, 38));
        fb.draw_text(
            28,
            60,
            &format!(
                "SUPER-BLAST ANNOTATION VALIDATION & FLAGGED NCBI ANOMALIES ({} Caught & Resolved Autonomously)",
                self.active_manifest.ncbi_flagged_anomalies.len()
            ),
            (248, 113, 113),
            1,
            Some(w - 56),
        );

        for (idx, anom) in self
            .active_manifest
            .ncbi_flagged_anomalies
            .iter()
            .take(6)
            .enumerate()
        {
            let y = 84 + (idx as i32) * 32;
            let col = if anom.severity.contains("CRITICAL") {
                (248, 113, 113)
            } else if anom.severity.contains("WARNING") {
                (251, 191, 36)
            } else {
                (56, 189, 248)
            };
            fb.fill_rect(26, y, w - 52, 28, (23, 33, 53));
            fb.stroke_rect(26, y, w - 52, 28, col);
            fb.draw_text(
                34,
                y + 8,
                &format!(
                    "[{}] {} ({}) -> {}",
                    anom.anomaly_code,
                    anom.accession,
                    anom.declared_organism,
                    anom.autonomous_resolution
                ),
                col,
                1,
                Some(w - 72),
            );
        }

        // Bottom Panel: Complete Audited Sequence Table (Downloads + Modifications + Replacements)
        let tbl_y = 296;
        let tbl_h = h - tbl_y - 38;
        fb.fill_rect(16, tbl_y, w - 32, tbl_h, (15, 23, 42));
        fb.stroke_rect(16, tbl_y, w - 32, tbl_h, (51, 65, 85));
        fb.draw_text(
            28,
            tbl_y + 12,
            "COMPLETE DOWNLOAD & SUPER-BLAST SEQUENCE AUDIT TABLE (TIP | ACCESSION | TAXON | SOURCE | LEN | STRAND | ID% | STOPS | VERDICT)",
            (56, 189, 248),
            1,
            Some(w - 56),
        );

        for (idx, s) in self.active_manifest.sequences.iter().take(22).enumerate() {
            let ry = tbl_y + 36 + (idx as i32) * 20;
            if ry + 18 > tbl_y + tbl_h {
                break;
            }
            let row_col = if s.included_in_cleaned_cohort {
                if s.audit_verdict == "RESCUED_REVCOMP" {
                    (251, 191, 36)
                } else if s.source_origin == "ReserveReplacement" {
                    (52, 211, 153)
                } else {
                    (226, 232, 240)
                }
            } else {
                (248, 113, 113)
            };
            fb.draw_text(
                28,
                ry,
                &format!(
                    "{:<6} {:<4} {:<13} {:<36} {:<18} {:>5} bp  strand={}  Id={:>5.1}%  stops={}  [{}]",
                    s.tip_code,
                    s.sector_code,
                    s.accession,
                    &s.organism[..s.organism.len().min(35)],
                    s.source_origin,
                    s.raw_length_bp,
                    s.strand_orientation,
                    s.pct_identity_vs_centroid,
                    s.internal_stop_codons,
                    s.audit_verdict
                ),
                row_col,
                1,
                Some(w - 56),
            );
        }
    }

    fn render_tab4_trees_and_skylines(
        &self,
        fb: &mut SoftwareFrameBuffer,
        _hits: &mut Vec<UiHitRect>,
    ) {
        let w = fb.width as i32;
        let h = fb.height as i32;
        let half_w = (w - 44) / 2;

        // Left Panel: Calibrated Chronogram & Quartet MSC Support Bars
        fb.fill_rect(16, 48, half_w, h - 86, (15, 23, 42));
        fb.stroke_rect(16, 48, half_w, h - 86, (2, 132, 199));
        fb.draw_text(
            28,
            60,
            "CALIBRATED CHRONOGRAM (NODE I-VII, 95% HPD) & QUARTET MSC [q1/q2/q3]",
            (56, 189, 248),
            1,
            Some(half_w - 24),
        );

        for (idx, n) in self.active_manifest.calibrated_nodes.iter().enumerate() {
            let ny = 104 + (idx as i32) * 92;
            fb.fill_rect(28, ny, half_w - 24, 78, (23, 33, 53));
            fb.stroke_rect(28, ny, half_w - 24, 78, (51, 65, 85));
            fb.draw_text(
                38,
                ny + 8,
                &format!(
                    "{} {} — {} | {:.2} Ma (95% HPD: {:.1}-{:.1} Ma)",
                    n.sector_code, n.node_id, n.clade_name, n.age_ma, n.ci_min_ma, n.ci_max_ma
                ),
                (255, 255, 255),
                1,
                Some(half_w - 44),
            );
            fb.draw_text(
                38,
                ny + 26,
                &format!(
                    "Guild: {} | UFBoot: {:.1}% | LPP: {:.2} ({:.2} CU)",
                    n.trophic_guild, n.ufboot_pct, n.quartet_lpp, n.coalescent_units_cu
                ),
                (148, 163, 184),
                1,
                Some(half_w - 44),
            );

            // Quartet [q1 | q2 | q3] visual bar
            let bw = 220.0_f64;
            let w1 = (bw * n.quartet_q1) as i32;
            let w2 = (bw * n.quartet_q2) as i32;
            let w3 = (bw * n.quartet_q3) as i32;
            fb.fill_rect(38, ny + 48, w1, 16, (2, 132, 199));
            fb.fill_rect(38 + w1, ny + 48, w2, 16, (245, 158, 11));
            fb.fill_rect(38 + w1 + w2, ny + 48, w3, 16, (220, 38, 38));
            fb.draw_text(
                270,
                ny + 52,
                &format!(
                    "Quartet q1={:.2} q2={:.2} q3={:.2}",
                    n.quartet_q1, n.quartet_q2, n.quartet_q3
                ),
                (56, 189, 248),
                1,
                Some(half_w - 280),
            );
        }

        // Right Panel: Coupled 3-PopSet Bayesian Coalescent Skylines
        let rx0 = 28 + half_w;
        fb.fill_rect(rx0, 48, half_w, h - 86, (15, 23, 42));
        fb.stroke_rect(rx0, 48, half_w, h - 86, (16, 185, 129));
        fb.draw_text(
            rx0 + 14,
            60,
            "COUPLED 3-POPSET BAYESIAN COALESCENT SKYLINES (RuNIE vs *BEAST vs 0-Beaker)",
            (52, 211, 153),
            1,
            Some(half_w - 28),
        );

        let colors = [(2, 132, 199), (13, 148, 136), (168, 85, 247)];
        for (p_idx, c) in self.active_manifest.skyline_cohorts.iter().enumerate() {
            let py0 = 92 + (p_idx as i32) * 224;
            let ph = 208;
            let col = colors[p_idx % colors.len()];
            fb.fill_rect(rx0 + 12, py0, half_w - 24, ph, (23, 33, 53));
            fb.stroke_rect(rx0 + 12, py0, half_w - 24, ph, col);
            fb.draw_text(
                rx0 + 22,
                py0 + 8,
                &format!(
                    "{} ({}; N={}) — {}",
                    c.species_name, c.ncbi_popset_id, c.sample_size_n, c.geographic_region
                ),
                (255, 255, 255),
                1,
                Some(half_w - 44),
            );
            fb.draw_text(
                rx0 + 22,
                py0 + 24,
                &format!(
                    "Locality: {} | Tajima D={:.2} | ESS={:.0}",
                    c.country_lat_lon, c.tajima_d, c.ess
                ),
                (148, 163, 184),
                1,
                Some(half_w - 44),
            );

            let gx0 = rx0 + 40;
            let gx1 = rx0 + half_w - 40;
            let gy0 = py0 + 48;
            let gy1 = py0 + ph - 18;
            let max_ne = 5_200_000.0_f64;

            for i in 0..c.epochs_ka.len().saturating_sub(1) {
                let x0 = gx0 + ((c.epochs_ka[i] / 22.0) * (gx1 - gx0) as f64) as i32;
                let x1 = gx0 + ((c.epochs_ka[i + 1] / 22.0) * (gx1 - gx0) as f64) as i32;

                let yr0 =
                    gy1 - ((c.runie_mcmc_ne_median[i] / max_ne) * (gy1 - gy0) as f64) as i32;
                let yr1 =
                    gy1 - ((c.runie_mcmc_ne_median[i + 1] / max_ne) * (gy1 - gy0) as f64) as i32;

                let yb0 =
                    gy1 - ((c.beast_mcc_ne_median[i] / max_ne) * (gy1 - gy0) as f64) as i32 - 5;
                let yb1 =
                    gy1 - ((c.beast_mcc_ne_median[i + 1] / max_ne) * (gy1 - gy0) as f64) as i32 - 5;

                let yz0 =
                    gy1 - ((c.zero_beaker_control_ne[i] / max_ne) * (gy1 - gy0) as f64) as i32;
                let yz1 =
                    gy1 - ((c.zero_beaker_control_ne[i + 1] / max_ne) * (gy1 - gy0) as f64) as i32;

                fb.draw_line(x0, yr0, x1, yr1, 3, col);
                fb.draw_line(x0, yb0, x1, yb1, 2, (245, 158, 11));
                fb.draw_line(x0, yz0, x1, yz1, 2, (239, 68, 68));
                fb.fill_circle(x0, yr0, 3, col);
            }
        }
    }

    fn render_tab5_pdb_and_figure_studio(
        &self,
        fb: &mut SoftwareFrameBuffer,
        hits: &mut Vec<UiHitRect>,
    ) {
        let w = fb.width as i32;
        let h = fb.height as i32;
        let half_w = (w - 44) / 2;

        // Left Panel: Interactive 3D PDB Structure & Active-Site Pocket Viewer
        fb.fill_rect(16, 48, half_w, h - 86, (15, 23, 42));
        fb.stroke_rect(16, 48, half_w, h - 86, (225, 29, 72));
        fb.draw_text(
            28,
            60,
            "INTERACTIVE 3D PDB STRUCTURE & CATALYTIC POCKET VIEWER (.pdb)",
            (251, 113, 133),
            1,
            Some(half_w - 28),
        );

        // 3D Rotation Controls
        Self::draw_pill_button(
            fb,
            hits,
            28,
            82,
            106,
            26,
            "Yaw Left <<",
            GUI_ACT_PDB_YAW_LEFT,
            false,
            (2, 132, 199),
        );
        Self::draw_pill_button(
            fb,
            hits,
            142,
            82,
            106,
            26,
            "Yaw Right >>",
            GUI_ACT_PDB_YAW_RIGHT,
            false,
            (2, 132, 199),
        );
        Self::draw_pill_button(
            fb,
            hits,
            256,
            82,
            106,
            26,
            "Pitch Up ^",
            GUI_ACT_PDB_PITCH_UP,
            false,
            (16, 185, 129),
        );
        Self::draw_pill_button(
            fb,
            hits,
            370,
            82,
            116,
            26,
            "Pitch Down v",
            GUI_ACT_PDB_PITCH_DOWN,
            false,
            (16, 185, 129),
        );
        Self::draw_pill_button(
            fb,
            hits,
            494,
            82,
            116,
            26,
            "Reset 3D View",
            GUI_ACT_PDB_RESET_VIEW,
            false,
            (217, 119, 6),
        );

        if let Some(pdb) = self.active_manifest.pdb_structures.first() {
            fb.draw_text(
                28,
                118,
                &format!(
                    "Target: {} | Template: {} | File: {}",
                    pdb.target_protein, pdb.template_pdb_code, pdb.pdb_file_rel_path
                ),
                (226, 232, 240),
                1,
                Some(half_w - 32),
            );
            fb.draw_text(
                28,
                136,
                &format!(
                    "Kabsch RMSD = {:.3} A | TM-score = {:.4} | Pocket Volume = {:.1} A^3 | Affinity = {:.2} kcal/mol",
                    pdb.kabsch_rmsd_angstrom,
                    pdb.tm_score,
                    pdb.pocket_volume_a3,
                    pdb.docking_affinity_kcal_mol
                ),
                (56, 189, 248),
                1,
                Some(half_w - 32),
            );

            let cx_3d = 16 + half_w / 2;
            let cy_3d = 420;
            let yaw = self.pdb_yaw_rad;
            let pitch = self.pdb_pitch_rad;
            let proj = |x: f64, y: f64, z: f64| -> (i32, i32) {
                let rx = x * yaw.cos() + z * yaw.sin();
                let rz = -x * yaw.sin() + z * yaw.cos();
                let ry = y * pitch.cos() - rz * pitch.sin();
                (cx_3d + (rx * 18.0) as i32, cy_3d - (ry * 14.0) as i32)
            };

            let (pcx, pcy) = proj(
                pdb.pocket_center_xyz.0,
                pdb.pocket_center_xyz.1,
                pdb.pocket_center_xyz.2,
            );
            fb.stroke_circle(pcx, pcy, 76, (245, 158, 11));

            for i in 0..pdb.atoms.len() {
                let a = &pdb.atoms[i];
                let (sx, sy) = proj(a.x, a.y, a.z);
                let (tx, ty) = proj(a.template_x, a.template_y, a.template_z);
                if i + 1 < pdb.atoms.len() {
                    let b = &pdb.atoms[i + 1];
                    let (nsx, nsy) = proj(b.x, b.y, b.z);
                    let (ntx, nty) = proj(b.template_x, b.template_y, b.template_z);
                    fb.draw_line(tx, ty, ntx, nty, 1, (100, 116, 139));
                    let col = match a.secondary_structure.as_str() {
                        "POCKET_ACTIVE" => (244, 63, 94),
                        "HELIX" => (56, 189, 248),
                        "STRAND" => (251, 191, 36),
                        _ => (52, 211, 153),
                    };
                    fb.draw_line(sx, sy, nsx, nsy, 3, col);
                }
                if a.secondary_structure == "POCKET_ACTIVE" {
                    fb.fill_circle(sx, sy, 5, (244, 63, 94));
                    fb.draw_text(
                        sx + 7,
                        sy - 4,
                        &format!("{}{}", a.residue_name_3, a.residue_seq),
                        (253, 164, 175),
                        1,
                        Some(60),
                    );
                } else {
                    fb.fill_circle(sx, sy, 2, (226, 232, 240));
                }
            }

            // Catalytic residue list at bottom of 3D panel
            fb.draw_text(
                28,
                h - 88,
                &format!(
                    "Catalytic Pocket Triad/Residues: {}",
                    pdb.catalytic_residues.join(" | ")
                ),
                (251, 191, 36),
                1,
                Some(half_w - 32),
            );
        }

        // Right Panel: Locked-Format Figure & Preprint (.tex + .pdf) Studio
        let rx0 = 28 + half_w;
        fb.fill_rect(rx0, 48, half_w, h - 86, (15, 23, 42));
        fb.stroke_rect(rx0, 48, half_w, h - 86, (56, 189, 248));
        fb.draw_text(
            rx0 + 14,
            60,
            "LOCKED-FORMAT FIGURE (.SVG + .PNG + .PDF) & PREPRINT (.TEX) STUDIO",
            (56, 189, 248),
            1,
            Some(half_w - 28),
        );

        Self::draw_pill_button(
            fb,
            hits,
            rx0 + 14,
            82,
            half_w - 28,
            30,
            ">> REGENERATE ALL 6 FIGURES (.SVG + .PNG + .PDF) & PREPRINT (main.tex) <<",
            GUI_ACT_REGEN_FIGURES_TEX,
            true,
            (22, 163, 74),
        );

        for (idx, fig) in self.active_manifest.figures.iter().enumerate() {
            let fy = 124 + (idx as i32) * 84;
            fb.fill_rect(rx0 + 14, fy, half_w - 28, 76, (23, 33, 53));
            fb.stroke_rect(rx0 + 14, fy, half_w - 28, 76, (51, 65, 85));
            fb.draw_text(
                rx0 + 22,
                fy + 8,
                &format!(
                    "[Fig {} | {}] {}",
                    fig.figure_number, fig.data_mapping_kind, fig.title
                ),
                (255, 255, 255),
                1,
                Some(half_w - 44),
            );
            fb.draw_text(
                rx0 + 22,
                fy + 26,
                &format!(
                    "Exports: {}  |  {}  |  {} ({}x{} px)",
                    fig.svg_rel_path, fig.png_rel_path, fig.pdf_rel_path, fig.width_px, fig.height_px
                ),
                (52, 211, 153),
                1,
                Some(half_w - 44),
            );
            fb.draw_text(
                rx0 + 22,
                fy + 46,
                &fig.caption_summary,
                (148, 163, 184),
                1,
                Some(half_w - 44),
            );
        }

        // Preprint bundle card at bottom right
        let py = 124 + 6 * 84;
        if py + 64 < h - 40 {
            fb.fill_rect(rx0 + 14, py, half_w - 28, 66, (30, 58, 138));
            fb.stroke_rect(rx0 + 14, py, half_w - 28, 66, (251, 191, 36));
            fb.draw_text(
                rx0 + 22,
                py + 10,
                &format!(
                    "PREPRINT LaTeX & PDF BUNDLE: {}",
                    self.active_manifest.latex_main_tex_rel_path
                ),
                (251, 191, 36),
                1,
                Some(half_w - 44),
            );
            fb.draw_text(
                rx0 + 22,
                py + 30,
                &format!(
                    "HTML: {}  |  PDF: {}",
                    self.active_manifest.html_manuscript_rel_path,
                    self.active_manifest.pdf_preprint_rel_path
                ),
                (255, 255, 255),
                1,
                Some(half_w - 44),
            );
        }
    }
}

/// Renders headless PNG snapshots of the multi-tab `bio-suite` GUI so both CLI verification
/// and documentation can capture exact pixel-locked screenshots of all 5 tabs.
pub fn render_gui_all_tabs_snapshots(
    manifest: &SupervisorRunManifest,
    runs_root_dir: &Path,
    primary_png_out: &Path,
) -> Result<(), String> {
    let mut state = BioSuiteGuiState::new(
        manifest.clone(),
        runs_root_dir.to_path_buf(),
        PathBuf::from("cache"),
    );
    let mut fb = SoftwareFrameBuffer::new(1480, 880);

    if let Some(parent) = primary_png_out.parent() {
        let _ = fs::create_dir_all(parent);
    }

    // Render Tab 1 to primary_png_out
    state.active_tab = BioSuiteGuiTab::RunsAndPipelineBuilder;
    let _ = state.render_bio_suite_gui_frame(&mut fb);
    fs::write(primary_png_out, fb.encode_png()).map_err(|e| e.to_string())?;

    // Also save all 5 tab screenshots into `runs/<run_id>/figures/`
    let fig_dir = runs_root_dir.join(&manifest.run_id).join("figures");
    let _ = fs::create_dir_all(&fig_dir);

    for (tab, fname) in [
        (
            BioSuiteGuiTab::RunsAndPipelineBuilder,
            "gui-tab1-runs-and-pipeline-dag.png",
        ),
        (
            BioSuiteGuiTab::SuperBlastDataAudit,
            "gui-tab2-super-blast-hygiene-audit.png",
        ),
        (
            BioSuiteGuiTab::AlignViewInspector,
            "gui-tab3-align-view-inspector.png",
        ),
        (
            BioSuiteGuiTab::TreesAndSkylines,
            "gui-tab4-chronogram-and-skylines.png",
        ),
        (
            BioSuiteGuiTab::PdbAndFigureStudio,
            "gui-tab5-pdb-and-figure-studio.png",
        ),
    ] {
        state.active_tab = tab;
        let _ = state.render_bio_suite_gui_frame(&mut fb);
        fs::write(fig_dir.join(fname), fb.encode_png()).map_err(|e| e.to_string())?;
    }

    Ok(())
}

#[cfg(target_os = "windows")]
mod win32_suite_gui {
    use super::*;
    use std::ffi::OsStr;
    use std::os::windows::ffi::OsStrExt;
    use std::ptr::{null, null_mut};

    type HWND = *mut core::ffi::c_void;
    type HINSTANCE = *mut core::ffi::c_void;
    type HMENU = *mut core::ffi::c_void;
    type HICON = *mut core::ffi::c_void;
    type HCURSOR = *mut core::ffi::c_void;
    type HBRUSH = *mut core::ffi::c_void;
    type HDC = *mut core::ffi::c_void;
    type UINT = u32;
    type DWORD = u32;
    type LONG = i32;
    type WPARAM = usize;
    type LPARAM = isize;
    type LRESULT = isize;
    type BOOL = i32;

    const WS_OVERLAPPEDWINDOW: DWORD = 0x00CF_0000;
    const WS_VISIBLE: DWORD = 0x1000_0000;
    const CW_USEDEFAULT: i32 = -2147483648i32;
    const CS_HREDRAW: UINT = 0x0002;
    const CS_VREDRAW: UINT = 0x0001;
    const IDC_ARROW: *const u16 = 32512 as *const u16;
    const WM_DESTROY: UINT = 0x0002;
    const WM_SIZE: UINT = 0x0005;
    const WM_PAINT: UINT = 0x000F;
    const WM_CLOSE: UINT = 0x0010;
    const WM_ERASEBKGND: UINT = 0x0014;
    const WM_KEYDOWN: UINT = 0x0100;
    const WM_LBUTTONDOWN: UINT = 0x0201;
    const WM_MOUSEWHEEL: UINT = 0x020A;
    const VK_F1: usize = 0x70;
    const VK_F2: usize = 0x71;
    const VK_F3: usize = 0x72;
    const VK_F4: usize = 0x73;
    const VK_F5: usize = 0x74;
    const BI_RGB: DWORD = 0;
    const DIB_RGB_COLORS: UINT = 0;
    const SRCCOPY: DWORD = 0x00CC_0020;

    #[repr(C)]
    struct WNDCLASSEXW {
        cb_size: UINT,
        style: UINT,
        lpfn_wnd_proc: Option<unsafe extern "system" fn(HWND, UINT, WPARAM, LPARAM) -> LRESULT>,
        cb_cls_extra: i32,
        cb_wnd_extra: i32,
        h_instance: HINSTANCE,
        h_icon: HICON,
        h_cursor: HCURSOR,
        hbr_background: HBRUSH,
        lpsz_menu_name: *const u16,
        lpsz_class_name: *const u16,
        h_icon_sm: HICON,
    }

    #[repr(C)]
    struct POINT {
        x: LONG,
        y: LONG,
    }

    #[repr(C)]
    struct MSG {
        hwnd: HWND,
        message: UINT,
        wparam: WPARAM,
        lparam: LPARAM,
        time: DWORD,
        pt: POINT,
    }

    #[repr(C)]
    struct RECT {
        left: LONG,
        top: LONG,
        right: LONG,
        bottom: LONG,
    }

    #[repr(C)]
    struct PAINTSTRUCT {
        hdc: HDC,
        f_erase: BOOL,
        rc_paint: RECT,
        f_restore: BOOL,
        f_inc_update: BOOL,
        rgb_reserved: [u8; 32],
    }

    #[repr(C)]
    struct BITMAPINFOHEADER {
        bi_size: DWORD,
        bi_width: LONG,
        bi_height: LONG,
        bi_planes: u16,
        bi_bit_count: u16,
        bi_compression: DWORD,
        bi_size_image: DWORD,
        bi_x_pels_per_meter: LONG,
        bi_y_pels_per_meter: LONG,
        bi_clr_used: DWORD,
        bi_clr_important: DWORD,
    }

    #[repr(C)]
    struct BITMAPINFO {
        bmi_header: BITMAPINFOHEADER,
        bmi_colors: [u32; 1],
    }

    #[link(name = "kernel32")]
    extern "system" {
        fn GetModuleHandleW(lp_module_name: *const u16) -> HINSTANCE;
    }

    #[link(name = "user32")]
    extern "system" {
        fn RegisterClassExW(lpwcx: *const WNDCLASSEXW) -> u16;
        fn CreateWindowExW(
            dw_ex_style: DWORD,
            lp_class_name: *const u16,
            lp_window_name: *const u16,
            dw_style: DWORD,
            x: i32,
            y: i32,
            n_width: i32,
            n_height: i32,
            h_wnd_parent: HWND,
            h_menu: HMENU,
            h_instance: HINSTANCE,
            lp_param: *mut core::ffi::c_void,
        ) -> HWND;
        fn DefWindowProcW(hwnd: HWND, msg: UINT, wparam: WPARAM, lparam: LPARAM) -> LRESULT;
        fn GetMessageW(lp_msg: *mut MSG, hwnd: HWND, min: UINT, max: UINT) -> BOOL;
        fn TranslateMessage(lp_msg: *const MSG) -> BOOL;
        fn DispatchMessageW(lp_msg: *const MSG) -> LRESULT;
        fn PostQuitMessage(n_exit_code: i32);
        fn DestroyWindow(hwnd: HWND) -> BOOL;
        fn LoadCursorW(h_instance: HINSTANCE, lp_cursor_name: *const u16) -> HCURSOR;
        fn BeginPaint(hwnd: HWND, lp_paint: *mut PAINTSTRUCT) -> HDC;
        fn EndPaint(hwnd: HWND, lp_paint: *const PAINTSTRUCT) -> BOOL;
        fn InvalidateRect(hwnd: HWND, lp_rect: *const RECT, b_erase: BOOL) -> BOOL;
        fn GetClientRect(hwnd: HWND, lp_rect: *mut RECT) -> BOOL;
    }

    #[link(name = "gdi32")]
    extern "system" {
        fn StretchDIBits(
            hdc: HDC,
            x_dest: i32,
            y_dest: i32,
            dest_width: i32,
            dest_height: i32,
            x_src: i32,
            y_src: i32,
            src_width: i32,
            src_height: i32,
            lp_bits: *const core::ffi::c_void,
            lpbmi: *const BITMAPINFO,
            i_usage: UINT,
            rop: DWORD,
        ) -> i32;
    }

    fn to_wide(s: &str) -> Vec<u16> {
        OsStr::new(s)
            .encode_wide()
            .chain(std::iter::once(0))
            .collect()
    }

    struct Win32SuiteContext {
        state: BioSuiteGuiState,
        fb: SoftwareFrameBuffer,
        hits: Vec<UiHitRect>,
    }

    static mut GLOBAL_SUITE_CTX: *mut Win32SuiteContext = null_mut();

    unsafe extern "system" fn suite_wnd_proc(
        hwnd: HWND,
        msg: UINT,
        wparam: WPARAM,
        lparam: LPARAM,
    ) -> LRESULT {
        let ptr = GLOBAL_SUITE_CTX;
        if ptr.is_null() {
            return DefWindowProcW(hwnd, msg, wparam, lparam);
        }
        let ctx = &mut *ptr;

        match msg {
            WM_ERASEBKGND => 1,
            WM_SIZE => {
                let mut rc = RECT {
                    left: 0,
                    top: 0,
                    right: 1480,
                    bottom: 880,
                };
                GetClientRect(hwnd, &mut rc);
                let w = (rc.right - rc.left).max(640) as usize;
                let h = (rc.bottom - rc.top).max(480) as usize;
                ctx.fb.resize(w, h);
                ctx.hits = ctx.state.render_bio_suite_gui_frame(&mut ctx.fb);
                InvalidateRect(hwnd, null(), 0);
                0
            }
            WM_LBUTTONDOWN => {
                let mx = (lparam & 0xFFFF) as i16 as i32;
                let my = ((lparam >> 16) & 0xFFFF) as i16 as i32;
                let w = ctx.fb.width as i32;
                let mut matched = None;
                for h_rect in &ctx.hits {
                    if mx >= h_rect.x
                        && mx <= h_rect.x + h_rect.w
                        && my >= h_rect.y
                        && my <= h_rect.y + h_rect.h
                    {
                        matched = Some(h_rect.action_id);
                        break;
                    }
                }
                if let Some(act) = matched {
                    ctx.state.handle_gui_action(act, w);
                    ctx.hits = ctx.state.render_bio_suite_gui_frame(&mut ctx.fb);
                    InvalidateRect(hwnd, null(), 0);
                }
                0
            }
            WM_MOUSEWHEEL => {
                if ctx.state.active_tab == BioSuiteGuiTab::AlignViewInspector {
                    let delta = ((wparam >> 16) & 0xFFFF) as i16;
                    let w = ctx.fb.width as i32;
                    let step = (ctx.state.align_viewer.visible_cols_count(w) / 8).max(3);
                    if delta > 0 {
                        ctx.state.align_viewer.scroll_col =
                            ctx.state.align_viewer.scroll_col.saturating_sub(step);
                    } else {
                        ctx.state.align_viewer.scroll_col += step;
                    }
                    ctx.state
                        .align_viewer
                        .clamp_scroll(w, ctx.fb.height as i32);
                    ctx.hits = ctx.state.render_bio_suite_gui_frame(&mut ctx.fb);
                    InvalidateRect(hwnd, null(), 0);
                }
                0
            }
            WM_KEYDOWN => {
                let w = ctx.fb.width as i32;
                match wparam {
                    VK_F1 => ctx.state.handle_gui_action(GUI_ACT_TAB_1, w),
                    VK_F2 => ctx.state.handle_gui_action(GUI_ACT_TAB_2, w),
                    VK_F3 => ctx.state.handle_gui_action(GUI_ACT_TAB_3, w),
                    VK_F4 => ctx.state.handle_gui_action(GUI_ACT_TAB_4, w),
                    VK_F5 => ctx.state.handle_gui_action(GUI_ACT_TAB_5, w),
                    _ => {}
                }
                ctx.hits = ctx.state.render_bio_suite_gui_frame(&mut ctx.fb);
                InvalidateRect(hwnd, null(), 0);
                0
            }
            WM_PAINT => {
                let mut ps: PAINTSTRUCT = core::mem::zeroed();
                let hdc = BeginPaint(hwnd, &mut ps);
                let w = ctx.fb.width as i32;
                let h = ctx.fb.height as i32;
                let bmi = BITMAPINFO {
                    bmi_header: BITMAPINFOHEADER {
                        bi_size: core::mem::size_of::<BITMAPINFOHEADER>() as DWORD,
                        bi_width: w,
                        bi_height: -h,
                        bi_planes: 1,
                        bi_bit_count: 32,
                        bi_compression: BI_RGB,
                        bi_size_image: 0,
                        bi_x_pels_per_meter: 0,
                        bi_y_pels_per_meter: 0,
                        bi_clr_used: 0,
                        bi_clr_important: 0,
                    },
                    bmi_colors: [0],
                };
                StretchDIBits(
                    hdc,
                    0,
                    0,
                    w,
                    h,
                    0,
                    0,
                    w,
                    h,
                    ctx.fb.pixels.as_ptr() as *const core::ffi::c_void,
                    &bmi,
                    DIB_RGB_COLORS,
                    SRCCOPY,
                );
                EndPaint(hwnd, &ps);
                0
            }
            WM_CLOSE => {
                DestroyWindow(hwnd);
                0
            }
            WM_DESTROY => {
                PostQuitMessage(0);
                0
            }
            _ => DefWindowProcW(hwnd, msg, wparam, lparam),
        }
    }

    pub fn run_win32_suite_window(state: BioSuiteGuiState) -> Result<(), String> {
        unsafe {
            let h_instance = GetModuleHandleW(null());
            let class_name = to_wide("CrvBioSuiteSupervisorGuiWndClass");
            let wc = WNDCLASSEXW {
                cb_size: core::mem::size_of::<WNDCLASSEXW>() as UINT,
                style: CS_HREDRAW | CS_VREDRAW,
                lpfn_wnd_proc: Some(suite_wnd_proc),
                cb_cls_extra: 0,
                cb_wnd_extra: 0,
                h_instance,
                h_icon: null_mut(),
                h_cursor: LoadCursorW(null_mut(), IDC_ARROW),
                hbr_background: null_mut(),
                lpsz_menu_name: null(),
                lpsz_class_name: class_name.as_ptr(),
                h_icon_sm: null_mut(),
            };
            RegisterClassExW(&wc);

            let init_w = 1480usize;
            let init_h = 880usize;
            let mut fb = SoftwareFrameBuffer::new(init_w, init_h);
            let hits = state.render_bio_suite_gui_frame(&mut fb);

            let win_title = to_wide(&format!(
                "CRV Science Labs — bio-suite-supervisor & Multi-Tab Bioinformatics Suite | Active Run: {}",
                state.active_manifest.run_id
            ));

            let mut ctx = Box::new(Win32SuiteContext { state, fb, hits });
            GLOBAL_SUITE_CTX = &mut *ctx as *mut Win32SuiteContext;

            let hwnd = CreateWindowExW(
                0,
                class_name.as_ptr(),
                win_title.as_ptr(),
                WS_OVERLAPPEDWINDOW | WS_VISIBLE,
                CW_USEDEFAULT,
                CW_USEDEFAULT,
                init_w as i32,
                init_h as i32,
                null_mut(),
                null_mut(),
                h_instance,
                null_mut(),
            );
            if hwnd.is_null() {
                GLOBAL_SUITE_CTX = null_mut();
                return Err("Failed to create Win32 window for bio-suite GUI".to_string());
            }

            let mut msg: MSG = core::mem::zeroed();
            while GetMessageW(&mut msg, null_mut(), 0, 0) > 0 {
                TranslateMessage(&msg);
                DispatchMessageW(&msg);
            }
            GLOBAL_SUITE_CTX = null_mut();
            Ok(())
        }
    }
}

/// Launches the interactive multi-tab `bio-suite` GUI window.
pub fn launch_bio_suite_gui(
    manifest: SupervisorRunManifest,
    runs_root_dir: PathBuf,
    cache_dir: PathBuf,
) -> Result<(), String> {
    let state = BioSuiteGuiState::new(manifest, runs_root_dir, cache_dir);
    #[cfg(target_os = "windows")]
    {
        win32_suite_gui::run_win32_suite_window(state)
    }
    #[cfg(not(target_os = "windows"))]
    {
        let _ = state;
        Err("Interactive Win32 GUI is only supported on Windows targets".to_string())
    }
}
