//! Pure-Rust Desktop GUI Monitor (`--mode gui`) built with `eframe` / `egui`.
//!
//! Runs the 3-stage neural network schedule on a background worker thread while streaming
//! real-time `NeuronTelemetrySnapshot` frames over `crossbeam-channel` into an interactive
//! dashboard with 4 live views:
//! 1. **OTU x Locus 90% Coverage Heatmap** (single-voucher vs. chimeric multi-specimen cells)
//! 2. **Chimeric Species-per-Genus OTU Provenance Table**
//! 3. **Non-Destructive Soft-Quarantine Ledger** (`prior_quarantine_count` & `skepticism_prior`)
//! 4. **3-Stage Neuron Feedback & Convergence Log**

#[cfg(feature = "gui")]
mod monitor {
    use anyhow::Result;
    use crossbeam_channel::Receiver;
    use eframe::egui;
    use std::sync::atomic::Ordering;
    use std::sync::Arc;

    use crate::config::ConsensusConfig;
    use crate::db::SequenceRecord;
    use crate::neurons::{
        run_neural_network_schedule, EngineControlFlags, NeuronTelemetrySnapshot,
    };

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum ActiveTab {
        CoverageMatrix,
        ChimericProvenance,
        SoftQuarantineLedger,
        FeedbackLog,
    }

    pub struct SeqNetMonitorApp {
        rx: Receiver<NeuronTelemetrySnapshot>,
        controls: Arc<EngineControlFlags>,
        latest: Option<NeuronTelemetrySnapshot>,
        consensus_history: Vec<(usize, f64, f64)>,
        tab: ActiveTab,
        min_coverage_gate: f64,
    }

    impl SeqNetMonitorApp {
        pub fn new(
            rx: Receiver<NeuronTelemetrySnapshot>,
            controls: Arc<EngineControlFlags>,
            min_coverage_gate: f64,
        ) -> Self {
            Self {
                rx,
                controls,
                latest: None,
                consensus_history: Vec::new(),
                tab: ActiveTab::CoverageMatrix,
                min_coverage_gate,
            }
        }
    }

    impl eframe::App for SeqNetMonitorApp {
        fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
            while let Ok(snap) = self.rx.try_recv() {
                if self
                    .consensus_history
                    .last()
                    .map(|(e, _, _)| *e != snap.epoch)
                    .unwrap_or(true)
                {
                    self.consensus_history.push((
                        snap.epoch,
                        snap.global_consensus_score,
                        snap.matrix_occupancy,
                    ));
                }
                self.latest = Some(snap);
            }

            if self.latest.as_ref().map(|s| !s.finished).unwrap_or(true) {
                ctx.request_repaint_after(std::time::Duration::from_millis(80));
            }

            egui::TopBottomPanel::top("seqnet_header").show(ctx, |ui| {
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    ui.heading("bio-suite-seqnet | 3-Stage Sequence Analysis Neural Monitor");
                    ui.separator();
                    if let Some(ref s) = self.latest {
                        ui.label(format!("DB ID: {}", s.database_id));
                        ui.separator();
                        ui.label(format!("Taxon: {}", s.taxon_group));
                        ui.separator();
                        ui.label(format!("Epoch: {}/{}", s.epoch, s.max_epochs));
                        ui.separator();
                        ui.label(&s.active_stage);
                    }
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        if ui.button("Force Finalize & Export NEXUS").clicked() {
                            self.controls.force_finalize.store(true, Ordering::Relaxed);
                        }
                        let is_paused = self.controls.paused.load(Ordering::Relaxed);
                        let btn_label = if is_paused { "Resume Schedule" } else { "Pause Schedule" };
                        if ui.button(btn_label).clicked() {
                            self.controls.paused.store(!is_paused, Ordering::Relaxed);
                        }
                    });
                });
                ui.add_space(4.0);
            });

            egui::CentralPanel::default().show(ctx, |ui| {
                let Some(ref snap) = self.latest else {
                    ui.centered_and_justified(|ui| {
                        ui.label("Initializing SQLite cache & Stage 1 Input Neurons...");
                    });
                    return;
                };

                // 3-Stage Neuron Cluster Status Cards
                ui.columns(3, |cols| {
                    cols[0].group(|ui| {
                        ui.strong("Stage 1: Input Neuron Cluster");
                        ui.label(format!("SQLite Warm-Start Hits: {}", snap.cache_hits));
                        ui.label(format!(
                            "NCBI Entrez Delta Fetched: {}",
                            snap.live_delta_fetched
                        ));
                        ui.label(format!(
                            "Total Accessions Pooled: {}",
                            snap.total_sequences_evaluated
                        ));
                    });

                    cols[1].group(|ui| {
                        ui.strong("Stage 2: Hybrid Processing Cluster");
                        ui.label(format!("Clean Validated Accessions: {}", snap.clean_candidates_count));
                        ui.colored_label(
                            egui::Color32::from_rgb(225, 145, 55),
                            format!("Soft-Quarantined (Withheld): {}", snap.quarantined_count),
                        );
                        ui.label(format!(
                            "Rehabilitated from Prior DB Flags: {} | Loss: {:.4}",
                            snap.rehabilitated_count, snap.neural_bce_loss
                        ));
                    });

                    cols[2].group(|ui| {
                        ui.strong("Stage 3: Output & Valuation Cluster");
                        ui.label(format!(
                            "Retained Species Terminals: {} ({} Genera, {} Chimeric)",
                            snap.retained_otus_count,
                            snap.retained_genera_count,
                            snap.chimeric_otus_count
                        ));
                        ui.add(
                            egui::ProgressBar::new(snap.matrix_occupancy as f32).text(format!(
                                "Matrix Fill: {:.1}% (Gate >= {:.0}%)",
                                snap.matrix_occupancy * 100.0,
                                self.min_coverage_gate * 100.0
                            )),
                        );
                        ui.add(
                            egui::ProgressBar::new(snap.global_consensus_score as f32).text(
                                format!("Consensus: {:.4}", snap.global_consensus_score),
                            ),
                        );
                    });
                });

                ui.add_space(8.0);
                ui.horizontal(|ui| {
                    ui.selectable_value(
                        &mut self.tab,
                        ActiveTab::CoverageMatrix,
                        "OTU x Locus Occupancy Matrix",
                    );
                    ui.selectable_value(
                        &mut self.tab,
                        ActiveTab::ChimericProvenance,
                        "Chimeric OTU Accession Provenance",
                    );
                    ui.selectable_value(
                        &mut self.tab,
                        ActiveTab::SoftQuarantineLedger,
                        "Soft-Quarantine Ledger (DB Flagged)",
                    );
                    ui.selectable_value(
                        &mut self.tab,
                        ActiveTab::FeedbackLog,
                        "Epoch Feedback & Console Log",
                    );
                });
                ui.separator();

                egui::ScrollArea::both()
                    .auto_shrink([false, false])
                    .show(ui, |ui| match self.tab {
                        ActiveTab::CoverageMatrix => {
                            egui::Grid::new("otu_locus_heatmap")
                                .striped(true)
                                .spacing([4.0, 3.0])
                                .show(ui, |ui| {
                                    ui.strong("Species Terminal");
                                    ui.strong("Genus");
                                    ui.strong("Fill %");
                                    for loc in &snap.locus_names {
                                        ui.strong(loc);
                                    }
                                    ui.end_row();

                                    for otu in &snap.retained_otus_preview {
                                        ui.label(&otu.terminal_label);
                                        ui.label(&otu.genus);
                                        ui.label(format!(
                                            "{:.0}%",
                                            otu.locus_coverage_ratio * 100.0
                                        ));
                                        for loc in &snap.locus_names {
                                            if let Some(slot) = otu.locus_slots.get(loc) {
                                                let color = if otu.is_chimeric {
                                                    egui::Color32::from_rgb(56, 178, 172)
                                                } else {
                                                    egui::Color32::from_rgb(72, 187, 120)
                                                };
                                                ui.colored_label(color, "■").on_hover_text(
                                                    format!(
                                                        "{} | Voucher: {} | Conf: {:.3}",
                                                        slot.accession_version,
                                                        slot.voucher_id,
                                                        slot.neural_confidence
                                                    ),
                                                );
                                            } else {
                                                ui.colored_label(
                                                    egui::Color32::from_rgb(110, 115, 125),
                                                    "·",
                                                );
                                            }
                                        }
                                        ui.end_row();
                                    }
                                });
                        }
                        ActiveTab::ChimericProvenance => {
                            egui::Grid::new("chimera_table")
                                .striped(true)
                                .spacing([12.0, 5.0])
                                .show(ui, |ui| {
                                    ui.strong("Species Terminal");
                                    ui.strong("Genus");
                                    ui.strong("Chimeric?");
                                    ui.strong("Specimens");
                                    ui.strong("Loci Fill");
                                    ui.strong("Mean Conf");
                                    ui.strong("NCBI Accessions (Locus=Accession[Voucher])");
                                    ui.end_row();

                                    for otu in &snap.retained_otus_preview {
                                        ui.label(&otu.terminal_label);
                                        ui.label(&otu.genus);
                                        ui.label(if otu.is_chimeric { "YES" } else { "Single" });
                                        ui.label(otu.distinct_specimen_count.to_string());
                                        ui.label(format!(
                                            "{}/{} ({:.1}%)",
                                            otu.loci_filled,
                                            otu.loci_target,
                                            otu.locus_coverage_ratio * 100.0
                                        ));
                                        ui.label(format!("{:.3}", otu.mean_neural_confidence));
                                        let accs = otu
                                            .locus_slots
                                            .values()
                                            .map(|s| {
                                                format!(
                                                    "{}={}[{}]",
                                                    s.locus_canonical,
                                                    s.accession_version,
                                                    s.voucher_id
                                                )
                                            })
                                            .collect::<Vec<_>>()
                                            .join(", ");
                                        ui.label(accs);
                                        ui.end_row();
                                    }
                                });
                        }
                        ActiveTab::SoftQuarantineLedger => {
                            egui::Grid::new("quarantine_table")
                                .striped(true)
                                .spacing([12.0, 5.0])
                                .show(ui, |ui| {
                                    ui.strong("Accession");
                                    ui.strong("Species");
                                    ui.strong("Locus");
                                    ui.strong("Status");
                                    ui.strong("Reason Code");
                                    ui.strong("Prior Strikes");
                                    ui.strong("Skepticism Prior");
                                    ui.strong("K2P Dist");
                                    ui.strong("Neural Score");
                                    ui.end_row();

                                    for q in &snap.quarantine_preview {
                                        ui.label(&q.accession_version);
                                        ui.label(&q.species_binomial);
                                        ui.label(&q.locus_canonical);
                                        if q.currently_quarantined {
                                            ui.colored_label(
                                                egui::Color32::from_rgb(230, 120, 60),
                                                "QUARANTINED",
                                            );
                                        } else {
                                            ui.colored_label(
                                                egui::Color32::from_rgb(80, 200, 120),
                                                "REHABILITATED",
                                            );
                                        }
                                        ui.label(&q.last_reason_code);
                                        ui.label(q.prior_quarantine_count.to_string());
                                        ui.label(format!("{:.3}", q.skepticism_prior));
                                        ui.label(format!("{:.4}", q.within_locus_k2p_dist));
                                        ui.label(format!("{:.3}", q.neural_consensus_score));
                                        ui.end_row();
                                    }
                                });
                        }
                        ActiveTab::FeedbackLog => {
                            for line in &snap.feedback_log {
                                ui.monospace(line);
                            }
                        }
                    });
            });
        }
    }

    pub fn launch_gui_with_worker(
        cfg: ConsensusConfig,
        seed_records: Option<Vec<SequenceRecord>>,
    ) -> Result<()> {
        let (tx, rx) = crossbeam_channel::unbounded();
        let controls = EngineControlFlags::new();
        let worker_controls = Arc::clone(&controls);
        let min_cov = cfg.min_locus_coverage_ratio;

        std::thread::spawn(move || {
            if let Err(err) =
                run_neural_network_schedule(&cfg, seed_records, Some(tx), worker_controls)
            {
                eprintln!("[bio-suite-seqnet worker error] {:#}", err);
            }
        });

        let native_options = eframe::NativeOptions {
            viewport: egui::ViewportBuilder::default()
                .with_inner_size([1360.0, 820.0])
                .with_title("bio-suite-seqnet — Live 3-Stage Sequence Analysis Neural Network"),
            ..Default::default()
        };

        eframe::run_native(
            "bio-suite-seqnet",
            native_options,
            Box::new(move |_cc| Ok(Box::new(SeqNetMonitorApp::new(rx, controls, min_cov)))),
        )
        .map_err(|e| anyhow::anyhow!("eframe GUI terminated with error: {}", e))
    }
}

#[cfg(feature = "gui")]
pub use monitor::launch_gui_with_worker;

#[cfg(not(feature = "gui"))]
pub fn launch_gui_with_worker(
    _cfg: crate::config::ConsensusConfig,
    _seed_records: Option<Vec<crate::db::SequenceRecord>>,
) -> anyhow::Result<()> {
    anyhow::bail!(
        "Binary was compiled without the `gui` feature. Rebuild with `cargo build --release --features gui` or run with `--mode headless`."
    )
}
