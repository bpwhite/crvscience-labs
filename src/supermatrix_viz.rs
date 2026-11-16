//! `supermatrix_viz` — Pure-Rust Publication-Grade Graphical Supermatrix Positional Coverage,
//! Partition Architecture, Indel/Missing-Data Heatmap & Circular Mitogenome Origin-Phasing SVG Renderer.

use crate::prng_util::xml_escape;

#[derive(Debug, Clone)]
pub struct SupermatrixPartitionSpan {
    pub locus_name: String,
    pub locus_class_label: String,
    pub protein_coding: bool,
    pub is_nuclear: bool,
    pub start_col_1based: usize,
    pub end_col_1based: usize,
    pub taxa_present: usize,
}

#[derive(Debug, Clone)]
pub struct SupermatrixTaxonRow {
    pub tip_code: String,
    pub sector_code: String,
    pub species_label: String,
    /// Aligned string where `A/C/G/T/R/Y` = valid site, `-` = internal gap within present locus, `?` = missing locus block.
    pub aligned_row: String,
}

#[derive(Debug, Clone)]
pub struct CircularPhasingComparison {
    pub n_mitogenomes: usize,
    pub unphased_aligned_cols: usize,
    pub phased_aligned_cols: usize,
    pub unphased_mean_occ_pct: f64,
    pub phased_mean_occ_pct: f64,
    pub unphased_profile_bins: Vec<f64>,
    pub phased_profile_bins: Vec<f64>,
    pub raw_seqnet_cols: usize,
    pub trimmed_supermatrix_cols: usize,
    pub single_gene_cols: usize,
    pub single_gene_ufboot_pct: f64,
    pub supermatrix_ufboot_pct: f64,
    pub winning_model_label: String,
    pub unfiltered_all_otu_ufboot_pct: Option<f64>,
}

fn short_locus_label(name: &str) -> String {
    match name {
        "rrnL_16S" => "16S".to_string(),
        "rrnS_12S" => "12S".to_string(),
        "SSU_18S" => "18S".to_string(),
        "LSU_28S" => "28S".to_string(),
        other => other.to_string(),
    }
}

fn sector_accent_color(sector: &str) -> &'static str {
    if sector.contains('A') {
        "#475569"
    } else if sector.contains('B') {
        "#0284c7"
    } else if sector.contains('C') {
        "#0f766e"
    } else if sector.contains('D') {
        "#b45309"
    } else if sector.contains('E') {
        "#dc2626"
    } else if sector.contains('F') {
        "#7c3aed"
    } else if sector.contains('G') {
        "#be185d"
    } else {
        "#334155"
    }
}

/// Renders a 3-panel publication-grade vector SVG (`1480 x 1020`) depicting:
/// - Panel A: Partition architecture ribbon + stacked per-position base occupancy (`ACGT`),
///   internal indel burden (`-`), missing locus fraction (`?`), and parsimony-informative %
/// - Panel B: N-taxon × position 2D occupancy/compartment heatmap aligned column-for-column with Panel A
/// - Panel C: Circular mitogenome origin-phasing occupancy (`Raw GenBank Start` vs `COX1 Codon-1 Rotated`)
///   and alignment/UFBoot2 recovery summary.
pub fn render_supermatrix_positional_coverage_svg(
    study_subtitle: &str,
    partitions: &[SupermatrixPartitionSpan],
    rows: &[SupermatrixTaxonRow],
    phasing: &CircularPhasingComparison,
) -> String {
    let w = 1480.0f64;
    let h = 1020.0f64;
    let x_plot0 = 136.0f64;
    let x_plot1 = 1314.0f64;
    let plot_w = x_plot1 - x_plot0; // 1178.0 px

    let total_cols = rows.first().map(|r| r.aligned_row.len()).unwrap_or(0).max(1);
    let n_taxa = rows.len().max(1);
    let n_bins = 320usize.min(total_cols).max(24);

    // Precompute column compartment type (`0 = mito PCG, 1 = mito rRNA, 2 = nuc`)
    let mut col_comp = vec![0u8; total_cols];
    for p in partitions {
        let s = p.start_col_1based.saturating_sub(1).min(total_cols);
        let e = p.end_col_1based.min(total_cols);
        let code = if p.is_nuclear {
            2u8
        } else if !p.protein_coding {
            1u8
        } else {
            0u8
        };
        for c in s..e {
            col_comp[c] = code;
        }
    }

    // Compute per-column valid base count, gap count, missing count, and parsimony informativeness
    let row_bytes: Vec<&[u8]> = rows.iter().map(|r| r.aligned_row.as_bytes()).collect();
    let mut col_valid = vec![0usize; total_cols];
    let mut col_gap = vec![0usize; total_cols];
    let mut col_pi = vec![false; total_cols];
    let mut total_valid_cells = 0usize;
    let mut total_pi_sites = 0usize;

    for c in 0..total_cols {
        let mut cnt = [0usize; 4];
        let mut vc = 0usize;
        let mut gc = 0usize;
        for &rb in &row_bytes {
            let b = rb.get(c).copied().unwrap_or(b'?').to_ascii_uppercase();
            match b {
                b'A' => {
                    cnt[0] += 1;
                    vc += 1;
                }
                b'C' => {
                    cnt[1] += 1;
                    vc += 1;
                }
                b'G' => {
                    cnt[2] += 1;
                    vc += 1;
                }
                b'T' | b'U' => {
                    cnt[3] += 1;
                    vc += 1;
                }
                b'R' | b'Y' => {
                    cnt[if b == b'R' { 0 } else { 1 }] += 1;
                    vc += 1;
                }
                b'-' => gc += 1,
                _ => {}
            }
        }
        col_valid[c] = vc;
        col_gap[c] = gc;
        total_valid_cells += vc;
        let ge2 = cnt.iter().filter(|&&x| x >= 2).count();
        if ge2 >= 2 {
            col_pi[c] = true;
            total_pi_sites += 1;
        }
    }

    let mean_occ_pct = (total_valid_cells as f64) / ((n_taxa * total_cols) as f64) * 100.0;
    let pi_pct = (total_pi_sites as f64) / (total_cols as f64) * 100.0;

    // Bin positional profile curves across `n_bins`
    let mut bin_valid_pct = vec![0.0f64; n_bins];
    let mut bin_gap_pct = vec![0.0f64; n_bins];
    let mut bin_pi_pct = vec![0.0f64; n_bins];

    for bi in 0..n_bins {
        let s = (bi * total_cols) / n_bins;
        let e = (((bi + 1) * total_cols) / n_bins).max(s + 1).min(total_cols);
        let span = (e - s) as f64;
        let mut sv = 0usize;
        let mut sg = 0usize;
        let mut sp = 0usize;
        for c in s..e {
            sv += col_valid[c];
            sg += col_gap[c];
            if col_pi[c] {
                sp += 1;
            }
        }
        bin_valid_pct[bi] = (sv as f64 / (span * (n_taxa as f64))) * 100.0;
        bin_gap_pct[bi] = (sg as f64 / (span * (n_taxa as f64))) * 100.0;
        bin_pi_pct[bi] = (sp as f64 / span) * 100.0;
    }

    let mut out = String::with_capacity(180_000);
    out.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\" font-family=\"'Segoe UI', Arial, Helvetica, sans-serif\">\n",
        w as i32, h as i32
    ));
    out.push_str(&format!(
        "<rect x=\"0\" y=\"0\" width=\"{}\" height=\"{}\" fill=\"#ffffff\"/>\n",
        w as i32, h as i32
    ));

    // ================= PANEL A =================
    let sub_short: String = if study_subtitle.chars().count() > 64 {
        let mut s: String = study_subtitle.chars().take(61).collect();
        s.push_str("...");
        s
    } else {
        study_subtitle.to_string()
    };
    out.push_str("<rect x=\"16\" y=\"12\" width=\"1448\" height=\"306\" rx=\"8\" fill=\"#f8fafc\" stroke=\"#cbd5e1\" stroke-width=\"1.2\"/>\n");
    out.push_str(&format!(
        "<text x=\"28\" y=\"31\" font-size=\"12.5\" font-weight=\"700\" fill=\"#0f172a\">Panel A | Partition Architecture &amp; Positional Occupancy (L = {} bp)</text>\n",
        total_cols
    ));
    out.push_str(&format!(
        "<text x=\"1452\" y=\"31\" text-anchor=\"end\" font-size=\"10.5\" fill=\"#334155\">{} | {} Loci | Mean Occ: {:.1}% | PI: {} ({:.1}%)</text>\n",
        xml_escape(&sub_short),
        partitions.len(),
        mean_occ_pct,
        total_pi_sites,
        pi_pct
    ));

    // Partition Ribbon (`y = 66..88`) with staggered gene labels at `y = 49` / `y = 61`
    let y_rib = 66.0f64;
    let h_rib = 20.0f64;
    for (pi, p) in partitions.iter().enumerate() {
        let s0 = p.start_col_1based.saturating_sub(1) as f64;
        let e0 = p.end_col_1based.min(total_cols) as f64;
        let rx = x_plot0 + (s0 / total_cols as f64) * plot_w;
        let rw = (((e0 - s0) / total_cols as f64) * plot_w).max(2.0);
        let fill = if p.is_nuclear {
            "#7c3aed"
        } else if !p.protein_coding {
            "#0284c7"
        } else {
            "#0f766e"
        };
        let opacity = if pi % 2 == 0 { "0.92" } else { "0.74" };
        out.push_str(&format!(
            "<rect x=\"{:.2}\" y=\"{:.1}\" width=\"{:.2}\" height=\"{:.1}\" fill=\"{}\" fill-opacity=\"{}\" stroke=\"#ffffff\" stroke-width=\"0.9\"><title>{} ({}..{} bp, {}/{} taxa)</title></rect>\n",
            rx, y_rib, rw, h_rib, fill, opacity,
            xml_escape(&p.locus_name), p.start_col_1based, p.end_col_1based, p.taxa_present, n_taxa
        ));
        let lbl = short_locus_label(&p.locus_name);
        let cx = rx + rw * 0.5;
        if rw >= 34.0 {
            out.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-size=\"9.5\" font-weight=\"700\" fill=\"#ffffff\">{}</text>\n",
                cx, y_rib + 13.5, xml_escape(&lbl)
            ));
        } else {
            let ly = if pi % 2 == 0 { 50.0 } else { 61.0 };
            out.push_str(&format!(
                "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#64748b\" stroke-width=\"0.8\"/>\n",
                cx, ly + 2.0, cx, y_rib
            ));
            out.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-size=\"8.5\" font-weight=\"700\" fill=\"#1e293b\">{}</text>\n",
                cx, ly, xml_escape(&lbl)
            ));
        }
    }

    // Stacked Positional Profile Plot (`y = 96..280`, height = 184)
    let y_p0 = 96.0f64;
    let y_p1 = 280.0f64;
    let ph = y_p1 - y_p0;
    out.push_str(&format!(
        "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"#f1f5f9\" stroke=\"#94a3b8\" stroke-width=\"1\"/>\n",
        x_plot0, y_p0, plot_w, ph
    ));
    for pct in [0, 25, 50, 75, 100] {
        let yy = y_p1 - (pct as f64 / 100.0) * ph;
        out.push_str(&format!(
            "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#cbd5e1\" stroke-width=\"0.7\" stroke-dasharray=\"2,3\"/>\n",
            x_plot0, yy, x_plot1, yy
        ));
        out.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-size=\"9.5\" fill=\"#475569\">{}%</text>\n",
            x_plot0 - 6.0, yy + 3.2, pct
        ));
    }
    out.push_str(&format!(
        "<text x=\"36\" y=\"{:.1}\" text-anchor=\"middle\" font-size=\"10.5\" font-weight=\"700\" fill=\"#334155\" transform=\"rotate(-90 36 {:.1})\">Column Occupancy (%)</text>\n",
        (y_p0 + y_p1) * 0.5, (y_p0 + y_p1) * 0.5
    ));

    // Stacked bars: ACGT bottom (`#0d9488`), Internal Indel `-` middle (`#f97316`), Missing `?` = `#f1f5f9` background
    let bw = plot_w / (n_bins as f64);
    for bi in 0..n_bins {
        let bx = x_plot0 + (bi as f64) * bw;
        let vh = (bin_valid_pct[bi].clamp(0.0, 100.0) / 100.0) * ph;
        let gh = (bin_gap_pct[bi].clamp(0.0, 100.0 - bin_valid_pct[bi]) / 100.0) * ph;
        if vh > 0.3 {
            out.push_str(&format!(
                "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"#0d9488\" fill-opacity=\"0.88\"/>\n",
                bx, y_p1 - vh, bw + 0.35, vh
            ));
        }
        if gh > 0.3 {
            out.push_str(&format!(
                "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"#f97316\" fill-opacity=\"0.86\"/>\n",
                bx, y_p1 - vh - gh, bw + 0.35, gh
            ));
        }
    }

    // Partition divider lines across Panel A
    for p in partitions.iter().skip(1) {
        let px = x_plot0 + ((p.start_col_1based.saturating_sub(1) as f64) / total_cols as f64) * plot_w;
        out.push_str(&format!(
            "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#475569\" stroke-opacity=\"0.35\" stroke-width=\"0.7\" stroke-dasharray=\"2,2\"/>\n",
            px, y_p0, px, y_p1
        ));
    }

    // Parsimony-Informative % polyline (`#facc15` glow + `#1e1b4b` crisp line)
    let mut pi_pts = String::with_capacity(n_bins * 14);
    for bi in 0..n_bins {
        let cx = x_plot0 + ((bi as f64) + 0.5) * bw;
        let cy = y_p1 - (bin_pi_pct[bi].clamp(0.0, 100.0) / 100.0) * ph;
        if bi > 0 {
            pi_pts.push(' ');
        }
        pi_pts.push_str(&format!("{:.1},{:.1}", cx, cy));
    }
    out.push_str(&format!(
        "<polyline points=\"{}\" fill=\"none\" stroke=\"#fde047\" stroke-width=\"2.6\" stroke-opacity=\"0.92\"/>\n",
        pi_pts
    ));
    out.push_str(&format!(
        "<polyline points=\"{}\" fill=\"none\" stroke=\"#1e1b4b\" stroke-width=\"1.25\"/>\n",
        pi_pts
    ));

    // Mean occupancy horizontal reference
    let y_mean = y_p1 - (mean_occ_pct.clamp(0.0, 100.0) / 100.0) * ph;
    out.push_str(&format!(
        "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#0f172a\" stroke-width=\"1.0\" stroke-dasharray=\"5,3\"/>\n",
        x_plot0, y_mean, x_plot1, y_mean
    ));

    // X-axis bp ticks below Panel A
    let tick_step = if total_cols > 16000 { 2500 } else if total_cols > 8000 { 1500 } else if total_cols > 3000 { 500 } else { 150 };
    let mut tc = 0usize;
    while tc <= total_cols {
        let tx = x_plot0 + (tc as f64 / total_cols as f64) * plot_w;
        out.push_str(&format!(
            "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#475569\" stroke-width=\"0.9\"/>\n",
            tx, y_p1, tx, y_p1 + 4.5
        ));
        out.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-size=\"9.2\" fill=\"#334155\">{} bp</text>\n",
            tx, y_p1 + 15.0, tc.max(1)
        ));
        tc += tick_step;
    }

    // Panel A Right-Margin Legend (`x = 1324..1456`)
    out.push_str("<rect x=\"1324\" y=\"68\" width=\"132\" height=\"212\" rx=\"6\" fill=\"#ffffff\" stroke=\"#cbd5e1\" stroke-width=\"1\"/>\n");
    out.push_str("<text x=\"1332\" y=\"85\" font-size=\"10\" font-weight=\"700\" fill=\"#0f172a\">Locus &amp; State Legend</text>\n");
    let a_items: [(&str, &str); 7] = [
        ("#0f766e", "Mito PCG (RY3)"),
        ("#0284c7", "Mito rRNA (4-St)"),
        ("#7c3aed", "Nuclear (4-St)"),
        ("#0d9488", "Valid ACGT %"),
        ("#f97316", "Indel Gap (-) %"),
        ("#e2e8f0", "Missing Locus (?)"),
        ("#1e1b4b", "Parsim. Inf. (PI %)"),
    ];
    for (idx, &(col_hex, label)) in a_items.iter().enumerate() {
        let ly = 96.0 + (idx as f64) * 25.0;
        if idx == 6 {
            out.push_str(&format!(
                "<line x1=\"1332\" y1=\"{:.1}\" x2=\"1346\" y2=\"{:.1}\" stroke=\"#fde047\" stroke-width=\"3\"/>\n<line x1=\"1332\" y1=\"{:.1}\" x2=\"1346\" y2=\"{:.1}\" stroke=\"#1e1b4b\" stroke-width=\"1.2\"/>\n",
                ly + 5.0, ly + 5.0, ly + 5.0, ly + 5.0
            ));
        } else {
            out.push_str(&format!(
                "<rect x=\"1332\" y=\"{:.1}\" width=\"13\" height=\"11\" rx=\"2\" fill=\"{}\" stroke=\"#94a3b8\" stroke-width=\"0.5\"/>\n",
                ly, col_hex
            ));
        }
        out.push_str(&format!(
            "<text x=\"1350\" y=\"{:.1}\" font-size=\"9.5\" fill=\"#1e293b\">{}</text>\n",
            ly + 9.0, xml_escape(label)
        ));
    }

    // Per-row occupancy, distinct locus count, and >= 80% / >= 2-loci phylogram inclusion status
    let mut row_stats: Vec<(f64, usize, bool)> = Vec::with_capacity(n_taxa);
    for r in rows {
        let rb = r.aligned_row.as_bytes();
        let valid_t = rb
            .iter()
            .filter(|&&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'R' | b'Y'))
            .count();
        let mut loci_present = 0usize;
        for p in partitions {
            let s0 = p.start_col_1based.saturating_sub(1).min(rb.len());
            let e0 = p.end_col_1based.min(rb.len());
            if s0 < e0 {
                let valid_p = rb[s0..e0]
                    .iter()
                    .filter(|&&b| {
                        matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'R' | b'Y')
                    })
                    .count();
                if valid_p >= 10 {
                    loci_present += 1;
                }
            }
        }
        let occ_pct = (valid_t as f64 / total_cols as f64) * 100.0;
        let gated_in = crate::supervisor_study::passes_supermatrix_phylogram_coverage_gate(
            valid_t,
            total_cols,
            loci_present,
        );
        row_stats.push((occ_pct, loci_present, gated_in));
    }
    let n_gated_in = row_stats.iter().filter(|s| s.2).count();
    let n_popgen_only = n_taxa.saturating_sub(n_gated_in);

    // ================= PANEL B =================
    out.push_str("<rect x=\"16\" y=\"328\" width=\"1448\" height=\"416\" rx=\"8\" fill=\"#f8fafc\" stroke=\"#cbd5e1\" stroke-width=\"1.2\"/>\n");
    out.push_str(&format!(
        "<text x=\"28\" y=\"347\" font-size=\"12.5\" font-weight=\"700\" fill=\"#0f172a\">Panel B | {}-Taxon x Position Occupancy ({} Phylogram &gt;=80%; {} Coalescent/PopGen-Only)</text>\n",
        n_taxa, n_gated_in, n_popgen_only
    ));
    out.push_str("<text x=\"1452\" y=\"347\" text-anchor=\"end\" font-size=\"10.2\" fill=\"#334155\">Teal = Mito PCG | Blue = Mito rRNA | Purple = Nuclear | Orange = Gap (-) | Slate = Missing (?)</text>\n");

    let y_m0 = 363.0f64;
    let y_m1 = 724.0f64;
    let row_h = (y_m1 - y_m0) / (n_taxa as f64);

    // Track sector spans for left-margin brackets
    let mut sec_spans: Vec<(String, usize, usize)> = Vec::new();
    for (ri, r) in rows.iter().enumerate() {
        if let Some(last) = sec_spans.last_mut() {
            if last.0 == r.sector_code {
                last.2 = ri;
                continue;
            }
        }
        sec_spans.push((r.sector_code.clone(), ri, ri));
    }
    for (sec, r0, r1) in &sec_spans {
        let sy0 = y_m0 + (*r0 as f64) * row_h + 0.5;
        let sy1 = y_m0 + ((*r1 + 1) as f64) * row_h - 0.5;
        let sc = sector_accent_color(sec);
        out.push_str(&format!(
            "<rect x=\"24\" y=\"{:.1}\" width=\"5\" height=\"{:.1}\" rx=\"2\" fill=\"{}\"/>\n",
            sy0, (sy1 - sy0).max(2.0), sc
        ));
        let mid_y = 0.5 * (sy0 + sy1) + 3.5;
        out.push_str(&format!(
            "<text x=\"34\" y=\"{:.1}\" font-size=\"9.5\" font-weight=\"700\" fill=\"{}\">{} (n={})</text>\n",
            mid_y, sc, xml_escape(sec), r1 - r0 + 1
        ));
        if *r0 > 0 {
            out.push_str(&format!(
                "<line x1=\"96\" y1=\"{:.1}\" x2=\"1452\" y2=\"{:.1}\" stroke=\"#64748b\" stroke-opacity=\"0.35\" stroke-width=\"0.7\"/>\n",
                sy0 - 0.5, sy0 - 0.5
            ));
        }
    }

    // Render 2D occupancy heatmap using run-length merged `<rect>` spans
    out.push_str(&format!(
        "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"#e2e8f0\" stroke=\"#94a3b8\" stroke-width=\"0.8\"/>\n",
        x_plot0, y_m0, plot_w, y_m1 - y_m0
    ));
    let cell_h = (row_h - 0.65).max(2.0);
    let mut row_occ_pcts = Vec::with_capacity(n_taxa);

    for (ri, r) in rows.iter().enumerate() {
        let ry = y_m0 + (ri as f64) * row_h + 0.3;
        let rb = r.aligned_row.as_bytes();
        let valid_t = rb.iter().filter(|&&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'R' | b'Y')).count();
        let t_pct = (valid_t as f64 / total_cols as f64) * 100.0;
        row_occ_pcts.push(t_pct);

        if ri % 2 == 0 {
            out.push_str(&format!(
                "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-size=\"5.3\" fill=\"#475569\">{}</text>\n",
                x_plot0 - 4.0, ry + cell_h - 0.5, xml_escape(&r.tip_code)
            ));
        }

        let mut run_start = 0usize;
        let mut run_state: u8 = 255;
        let flush_run = |rs: usize, re: usize, state: u8, buf: &mut String| {
            // 0 = miss (`?` — skipped since background is #e2e8f0), 1 = mito PCG, 2 = mito rRNA, 3 = nuc, 4 = gap (`-`)
            let col_opt = match state {
                1 => Some("#0f766e"),
                2 => Some("#0284c7"),
                3 => Some("#7c3aed"),
                4 => Some("#fb923c"),
                _ => None,
            };
            if let Some(fill_c) = col_opt {
                let rx = x_plot0 + (rs as f64) * bw;
                let rw = ((re - rs) as f64) * bw + 0.25;
                buf.push_str(&format!(
                    "<rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"{}\"/>",
                    rx, ry, rw, cell_h, fill_c
                ));
            }
        };

        for bi in 0..n_bins {
            let s = (bi * total_cols) / n_bins;
            let e = (((bi + 1) * total_cols) / n_bins).max(s + 1).min(total_cols);
            let mut n_base = 0usize;
            let mut n_gap = 0usize;
            for c in s..e {
                let ch = rb.get(c).copied().unwrap_or(b'?').to_ascii_uppercase();
                if matches!(ch, b'A' | b'C' | b'G' | b'T' | b'R' | b'Y') {
                    n_base += 1;
                } else if ch == b'-' {
                    n_gap += 1;
                }
            }
            let mid_c = (s + e) / 2;
            let state = if n_base > 0 && n_base * 3 >= (e - s) {
                match col_comp.get(mid_c).copied().unwrap_or(0) {
                    2 => 3,
                    1 => 2,
                    _ => 1,
                }
            } else if n_gap > 0 && (n_base + n_gap) * 2 >= (e - s) {
                4
            } else {
                0
            };
            if state != run_state {
                if run_state != 255 {
                    flush_run(run_start, bi, run_state, &mut out);
                }
                run_start = bi;
                run_state = state;
            }
        }
        if run_state != 255 {
            flush_run(run_start, n_bins, run_state, &mut out);
        }
        out.push('\n');

        // Right-margin per-taxon completeness bar (`x = 1324..1428`) colored by 80% phylogram gate
        let bar_w = 104.0f64;
        let fill_w = (t_pct.clamp(0.0, 100.0) / 100.0) * bar_w;
        out.push_str(&format!(
            "<rect x=\"1324\" y=\"{:.2}\" width=\"{:.1}\" height=\"{:.2}\" fill=\"#e2e8f0\"/>\n",
            ry, bar_w, cell_h
        ));
        let gated_in = row_stats.get(ri).map(|s| s.2).unwrap_or(false);
        let bar_col = if gated_in { "#0f766e" } else { "#f59e0b" };
        out.push_str(&format!(
            "<rect x=\"1324\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"{}\"/>\n",
            ry, fill_w.max(1.5), cell_h, bar_col
        ));
        if ri % 4 == 0 {
            out.push_str(&format!(
                "<text x=\"1433\" y=\"{:.2}\" font-size=\"6.2\" fill=\"#334155\">{:.0}%</text>\n",
                ry + cell_h, t_pct
            ));
        }
    }
    let gate_x = 1324.0 + 0.80 * 104.0;
    out.push_str(&format!(
        "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#b91c1c\" stroke-width=\"1.4\" stroke-dasharray=\"4,3\"/>\n",
        gate_x, y_m0, gate_x, y_m1
    ));
    out.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-size=\"8.2\" font-weight=\"700\" fill=\"#b91c1c\">80% Gate</text>\n",
        gate_x, y_m0 - 2.5
    ));
    out.push_str("<text x=\"1376\" y=\"737\" text-anchor=\"middle\" font-size=\"9.2\" font-weight=\"700\" fill=\"#334155\">Completeness (Green &gt;=80%)</text>\n");

    // ================= PANEL C =================
    out.push_str("<rect x=\"16\" y=\"756\" width=\"1448\" height=\"250\" rx=\"8\" fill=\"#f8fafc\" stroke=\"#cbd5e1\" stroke-width=\"1.2\"/>\n");
    out.push_str("<text x=\"28\" y=\"776\" font-size=\"13.5\" font-weight=\"700\" fill=\"#0f172a\">Panel C | Circular Mitogenome Origin-Position Phasing (Table-5 COX1 Rotation) &amp; Partitioned UFBoot2 Recovery</text>\n");

    // Sub-panel C1 (`x = 70..830, y = 796..972`): Unphased vs COX1-Rotated Positional Occupancy Profile
    let cx0 = 76.0f64;
    let cx1 = 824.0f64;
    let cy0 = 804.0f64;
    let cy1 = 970.0f64;
    let cw = cx1 - cx0;
    let ch = cy1 - cy0;
    out.push_str(&format!(
        "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"#ffffff\" stroke=\"#94a3b8\" stroke-width=\"1\"/>\n",
        cx0, cy0, cw, ch
    ));
    out.push_str("<text x=\"82\" y=\"797\" font-size=\"10.5\" font-weight=\"700\" fill=\"#1e293b\">C1. Mitogenome gotoh-msa Positional Occupancy Before vs. After COX1 Codon-1 Rotation</text>\n");
    for pct in [0, 50, 100] {
        let yy = cy1 - (pct as f64 / 100.0) * ch;
        out.push_str(&format!(
            "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#e2e8f0\" stroke-width=\"0.8\"/>\n",
            cx0, yy, cx1, yy
        ));
        out.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-size=\"9\" fill=\"#475569\">{}%</text>\n",
            cx0 - 5.0, yy + 3.0, pct
        ));
    }

    let draw_curve = |vals: &[f64], color: &str, dashed: bool, buf: &mut String| {
        if vals.is_empty() {
            return;
        }
        let mut pts = String::with_capacity(vals.len() * 12);
        let m = (vals.len() - 1).max(1) as f64;
        for (i, &v) in vals.iter().enumerate() {
            let x = cx0 + (i as f64 / m) * cw;
            let y = cy1 - (v.clamp(0.0, 100.0) / 100.0) * ch;
            if i > 0 {
                pts.push(' ');
            }
            pts.push_str(&format!("{:.1},{:.1}", x, y));
        }
        let dash_attr = if dashed { " stroke-dasharray=\"6,4\"" } else { "" };
        buf.push_str(&format!(
            "<polyline points=\"{}\" fill=\"none\" stroke=\"{}\" stroke-width=\"2.1\"{}/>\n",
            pts, color, dash_attr
        ));
    };
    draw_curve(&phasing.unphased_profile_bins, "#dc2626", true, &mut out);
    draw_curve(&phasing.phased_profile_bins, "#0f766e", false, &mut out);

    let pct_drop = if phasing.unphased_aligned_cols > phasing.phased_aligned_cols && phasing.unphased_aligned_cols > 0 {
        (1.0 - (phasing.phased_aligned_cols as f64 / phasing.unphased_aligned_cols as f64)) * 100.0
    } else {
        0.0
    };
    out.push_str(&format!(
        "<text x=\"{:.1}\" y=\"988\" text-anchor=\"middle\" font-size=\"9.5\" fill=\"#334155\">Normalized Mitogenome Alignment Coordinate (0% -&gt; 100% of Aligned Columns; n = {} Complete Mitogenomes)</text>\n",
        0.5 * (cx0 + cx1), phasing.n_mitogenomes
    ));
    out.push_str("<rect x=\"514\" y=\"812\" width=\"302\" height=\"56\" rx=\"5\" fill=\"#ffffff\" fill-opacity=\"0.94\" stroke=\"#cbd5e1\" stroke-width=\"0.9\"/>\n");
    out.push_str(&format!(
        "<line x1=\"524\" y1=\"829\" x2=\"546\" y2=\"829\" stroke=\"#dc2626\" stroke-width=\"2.1\" stroke-dasharray=\"6,4\"/><text x=\"552\" y=\"833\" font-size=\"9.3\" font-weight=\"700\" fill=\"#991b1b\">Unphased Start: L = {} cols (occ {:.1}%)</text>\n",
        phasing.unphased_aligned_cols, phasing.unphased_mean_occ_pct
    ));
    out.push_str(&format!(
        "<line x1=\"524\" y1=\"851\" x2=\"546\" y2=\"851\" stroke=\"#0f766e\" stroke-width=\"2.1\"/><text x=\"552\" y=\"855\" font-size=\"9.3\" font-weight=\"700\" fill=\"#0f766e\">COX1-Rotated: L = {} cols ({:.1}% | -{:.0}% bloat)</text>\n",
        phasing.phased_aligned_cols, phasing.phased_mean_occ_pct, pct_drop
    ));

    // Sub-panel C2 (`x = 846..1452, y = 788..994`): Diagnostic Resolution Comparison
    out.push_str("<rect x=\"846\" y=\"788\" width=\"606\" height=\"206\" rx=\"6\" fill=\"#ffffff\" stroke=\"#cbd5e1\" stroke-width=\"1\"/>\n");
    out.push_str("<text x=\"858\" y=\"806\" font-size=\"10.8\" font-weight=\"700\" fill=\"#0f172a\">C2. Single-Locus &amp; Unfiltered vs. >=80% Gated Supermatrix Recovery</text>\n");

    let mut cards: Vec<(String, String, f64, &str)> = vec![
        (
            "Legacy 1-Locus Truncated COX1 Only".to_string(),
            format!(
                "{} bp aligned | UFBoot2 = {:.1}%",
                phasing.single_gene_cols, phasing.single_gene_ufboot_pct
            ),
            (phasing.single_gene_ufboot_pct / 100.0).clamp(0.05, 1.0),
            "#f59e0b",
        ),
        (
            "Untrimmed Center-Star seqnet NEXUS".to_string(),
            format!("{} cols (singlet flank bloat + chimeras)", phasing.raw_seqnet_cols),
            0.52,
            "#ef4444",
        ),
    ];
    if let Some(unfilt_uf) = phasing.unfiltered_all_otu_ufboot_pct {
        cards.push((
            format!("Unfiltered All-{}-OTU Supermatrix", n_taxa),
            format!(
                "{} sparse OTUs included | UFBoot2 = {:.1}%",
                n_popgen_only, unfilt_uf
            ),
            (unfilt_uf / 100.0).clamp(0.05, 1.0),
            "#ea580c",
        ));
    }
    cards.push((
        format!(">=80% Gated {}-OTU Supermatrix", n_gated_in),
        format!(
            "{} bp ({} loci) | {}: {:.1}% UFBoot2",
            phasing.trimmed_supermatrix_cols,
            partitions.len(),
            xml_escape(&phasing.winning_model_label),
            phasing.supermatrix_ufboot_pct
        ),
        (phasing.supermatrix_ufboot_pct / 100.0).clamp(0.05, 1.0),
        "#0f766e",
    ));
    let (cy0, dy, bar_h, fsz) = if cards.len() >= 4 {
        (823.0, 42.5, 12.0, 9.4)
    } else {
        (828.0, 54.0, 16.0, 10.0)
    };
    for (ci, (title, detail, frac, bar_c)) in cards.iter().enumerate() {
        let cy = cy0 + (ci as f64) * dy;
        out.push_str(&format!(
            "<text x=\"862\" y=\"{:.1}\" font-size=\"{:.1}\" font-weight=\"700\" fill=\"#1e293b\">{}</text>\n",
            cy, fsz, xml_escape(title)
        ));
        out.push_str(&format!(
            "<text x=\"1438\" y=\"{:.1}\" text-anchor=\"end\" font-size=\"{:.1}\" fill=\"#475569\">{}</text>\n",
            cy, fsz - 0.5, detail
        ));
        out.push_str(&format!(
            "<rect x=\"862\" y=\"{:.1}\" width=\"576\" height=\"{:.1}\" rx=\"4\" fill=\"#f1f5f9\" stroke=\"#cbd5e1\" stroke-width=\"0.7\"/>\n",
            cy + 5.5, bar_h
        ));
        out.push_str(&format!(
            "<rect x=\"862\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"4\" fill=\"{}\"/>\n",
            cy + 5.5, 576.0 * frac, bar_h, bar_c
        ));
    }

    out.push_str("</svg>\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::supervisor_study::passes_supermatrix_phylogram_coverage_gate;

    #[test]
    fn supermatrix_80pct_gate_and_svg_cards_render_cleanly() {
        assert!(!passes_supermatrix_phylogram_coverage_gate(79, 100, 3));
        assert!(!passes_supermatrix_phylogram_coverage_gate(95, 100, 1));
        assert!(passes_supermatrix_phylogram_coverage_gate(80, 100, 2));

        let parts = vec![
            SupermatrixPartitionSpan {
                locus_name: "COX1".to_string(),
                locus_class_label: "Mito PCG".to_string(),
                is_nuclear: false,
                protein_coding: true,
                start_col_1based: 1,
                end_col_1based: 120,
                taxa_present: 2,
            },
            SupermatrixPartitionSpan {
                locus_name: "CYTB".to_string(),
                locus_class_label: "Mito PCG".to_string(),
                is_nuclear: false,
                protein_coding: true,
                start_col_1based: 121,
                end_col_1based: 240,
                taxa_present: 1,
            },
        ];
        let full_row = format!("{}{}", "ACGT".repeat(30), "ACGT".repeat(30));
        let cox1_only = format!("{}{}", "ACGT".repeat(30), "?".repeat(120));
        let rows = vec![
            SupermatrixTaxonRow {
                tip_code: "A:T001".to_string(),
                sector_code: "[A]".to_string(),
                species_label: "Taxon One".to_string(),
                aligned_row: full_row,
            },
            SupermatrixTaxonRow {
                tip_code: "B:T002".to_string(),
                sector_code: "[B]".to_string(),
                species_label: "Taxon Two".to_string(),
                aligned_row: cox1_only,
            },
        ];
        let phasing = CircularPhasingComparison {
            n_mitogenomes: 1,
            unphased_aligned_cols: 500,
            phased_aligned_cols: 260,
            unphased_mean_occ_pct: 55.0,
            phased_mean_occ_pct: 96.0,
            unphased_profile_bins: vec![55.0; 16],
            phased_profile_bins: vec![96.0; 16],
            raw_seqnet_cols: 900,
            trimmed_supermatrix_cols: 240,
            single_gene_cols: 120,
            single_gene_ufboot_pct: 48.2,
            supermatrix_ufboot_pct: 99.1,
            winning_model_label: "TN+F+I+G4 (RY3)".to_string(),
            unfiltered_all_otu_ufboot_pct: Some(59.4),
        };
        let svg = render_supermatrix_positional_coverage_svg("Gate Unit Test", &parts, &rows, &phasing);
        assert!(svg.contains("1 Phylogram &gt;=80%; 1 Coalescent/PopGen-Only"));
        assert!(svg.contains("80% Gate"));
        assert!(svg.contains("Unfiltered All-2-OTU Supermatrix"));
        assert!(svg.contains("&gt;=80% Gated 1-OTU Supermatrix"));
    }
}
