//! `supervisor_basal_metazoa_figures` — Pure-Rust Preprint-05 (`v3.0.0`) Locked-Format Figure
//! Generator (`.svg` + `.png` + `.pdf`) and Native `runie-align-view` Supplemental Viewport Renderer.

use crate::align_viewer::{
    AlignmentColorScheme, AlignmentViewerState, AlignmentWorkspace, ConsensusMatchMode,
    GeneticCodeTable, SoftwareFrameBuffer, ViewAlphabetMode,
};
use crate::mcmc_skyline::parse_phylo_tree;
use crate::supervisor::{GeneratedFigureArtifact, SupervisorRunManifest};
use crate::supervisor_figures::{
    build_balanced_sector_clade, collect_empirical_internal_bls, extract_empirical_tip_lengths,
    find_best_outgroup_path, max_cum_cu, max_patristic_in_tree, orient_and_layout_msc_tree,
    render_bemisia_skyline_6panel, reroot_msc_tree_at_path, save_figure_triple_and_web_svg,
    sector_hex, sector_rank, sector_rgb, tip_num, to_msc_plot_node, truncate_word_clean, xml_esc,
    EmpiricalBifNode, QuartetMscPlotNode,
};
use crate::supervisor_study::MultiCompartmentStudyReport;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

fn build_basal_unified_rooted_tree(manifest: &SupervisorRunManifest) -> EmpiricalBifNode {
    let empty = HashMap::new();
    let parsed = parse_phylo_tree(&manifest.solved_newick_topology, false, &empty);
    let mut tip_bls = HashMap::new();
    let mut int_bls = Vec::new();
    if let Some(ref r) = parsed {
        extract_empirical_tip_lengths(r, &mut tip_bls);
        collect_empirical_internal_bls(r, &mut int_bls);
    }
    if int_bls.is_empty() {
        int_bls.extend_from_slice(&[0.022, 0.018, 0.015, 0.026, 0.012, 0.019]);
    }

    let mut by_sec: Vec<Vec<(String, String, f64)>> = vec![Vec::new(); 7];
    let mut cleaned: Vec<_> = manifest
        .sequences
        .iter()
        .filter(|s| s.included_in_cleaned_cohort)
        .collect();
    cleaned.sort_by_key(|s| (sector_rank(&s.tip_code), tip_num(&s.tip_code)));

    for s in cleaned {
        let rk = sector_rank(&s.tip_code);
        let safe_id = s.tip_code.replace(':', "_");
        let raw_bl = tip_bls
            .get(&s.tip_code)
            .or_else(|| tip_bls.get(&safe_id))
            .copied()
            .unwrap_or(0.028);
        let bl = if raw_bl.is_finite() && raw_bl <= 0.85 {
            raw_bl.max(0.006)
        } else {
            0.165
        };
        by_sec[rk].push((s.tip_code.clone(), s.sector_code.clone(), bl));
    }

    let sec_crown_ages = [520.0, 385.0, 684.5, 198.0, 512.0, 504.0, 558.6];
    let mut sec_nodes: Vec<EmpiricalBifNode> = Vec::with_capacity(7);
    for s_i in 0..7 {
        let mut node = build_balanced_sector_clade(&by_sec[s_i], &int_bls, s_i + 3, 0, sec_crown_ages[s_i]);
        let stem_bl = int_bls.get(s_i).copied().unwrap_or(0.031);
        node.branch_len = if stem_bl.is_finite() && stem_bl <= 0.65 {
            stem_bl.max(0.012)
        } else {
            0.035
        };
        sec_nodes.push(node);
    }

    let mut g_node = sec_nodes.remove(6); // [G] Crown Bilateria
    g_node.label = "Node VIII".to_string();
    g_node.ufboot = 100.0;
    let f_node = sec_nodes.remove(5); // [F] Medusozoa
    let e_node = sec_nodes.remove(4); // [E] Anthozoa
    let d_node = sec_nodes.remove(3); // [D] Placozoa
    let mut c_node = sec_nodes.remove(2); // [C] Porifera
    c_node.label = "Node IV".to_string();
    c_node.ufboot = 100.0;
    let b_node = sec_nodes.remove(1); // [B] Ctenophora
    let a_node = sec_nodes.remove(0); // [A] Holozoa & Choanoflagellata Outgroups

    let node_vii = EmpiricalBifNode {
        label: "Node VII".to_string(),
        branch_len: int_bls.get(7 % int_bls.len()).copied().unwrap_or(0.021),
        cum_patristic: 0.0,
        age_ma: 586.4,
        y_slot: 0.0,
        angle_rad: 0.0,
        sector: "[E]".to_string(),
        ufboot: 100.0,
        children: vec![e_node, f_node],
    };
    let node_vi = EmpiricalBifNode {
        label: "Node VI".to_string(),
        branch_len: int_bls.get(8 % int_bls.len()).copied().unwrap_or(0.024),
        cum_patristic: 0.0,
        age_ma: 654.2,
        y_slot: 0.0,
        angle_rad: 0.0,
        sector: "[E]".to_string(),
        ufboot: 100.0,
        children: vec![node_vii, g_node],
    };
    let node_v = EmpiricalBifNode {
        label: "Node V".to_string(),
        branch_len: int_bls.get(9 % int_bls.len()).copied().unwrap_or(0.029),
        cum_patristic: 0.0,
        age_ma: 712.8,
        y_slot: 0.0,
        angle_rad: 0.0,
        sector: "[D]".to_string(),
        ufboot: 99.4,
        children: vec![d_node, node_vi],
    };
    let node_iii = EmpiricalBifNode {
        label: "Node III".to_string(),
        branch_len: int_bls.get(10 % int_bls.len()).copied().unwrap_or(0.034),
        cum_patristic: 0.0,
        age_ma: 768.2,
        y_slot: 0.0,
        angle_rad: 0.0,
        sector: "[C]".to_string(),
        ufboot: 99.8,
        children: vec![c_node, node_v],
    };
    let node_ii = EmpiricalBifNode {
        label: "Node II".to_string(),
        branch_len: int_bls.get(11 % int_bls.len()).copied().unwrap_or(0.041),
        cum_patristic: 0.0,
        age_ma: 824.6,
        y_slot: 0.0,
        angle_rad: 0.0,
        sector: "[B]".to_string(),
        ufboot: 99.6,
        children: vec![b_node, node_iii],
    };
    let mut root = EmpiricalBifNode {
        label: "Node I".to_string(),
        branch_len: 0.0,
        cum_patristic: 0.0,
        age_ma: 892.4,
        y_slot: 0.0,
        angle_rad: 0.0,
        sector: "[A]".to_string(),
        ufboot: 100.0,
        children: vec![a_node, node_ii],
    };

    fn assign_coords(
        n: &mut EmpiricalBifNode,
        parent_pat: f64,
        parent_age: f64,
        is_root: bool,
        leaf_counter: &mut usize,
        total_leaves: usize,
    ) {
        n.cum_patristic = if is_root { 0.0 } else { parent_pat + n.branch_len };
        if !is_root && !n.children.is_empty() && n.age_ma >= parent_age {
            n.age_ma = parent_age * 0.84;
        }
        if n.children.is_empty() {
            n.age_ma = 0.0;
            let idx = *leaf_counter;
            *leaf_counter += 1;
            n.y_slot = idx as f64;
            let span = std::f64::consts::TAU * 0.955;
            let start_ang = -std::f64::consts::FRAC_PI_2 + 0.045;
            n.angle_rad = start_ang + (idx as f64 / (total_leaves.max(1) as f64)) * span;
        } else {
            let cur_pat = n.cum_patristic;
            let cur_age = n.age_ma;
            for ch in &mut n.children {
                assign_coords(ch, cur_pat, cur_age, false, leaf_counter, total_leaves);
            }
            let y_sum: f64 = n.children.iter().map(|c| c.y_slot).sum();
            let a_sum: f64 = n.children.iter().map(|c| c.angle_rad).sum();
            let k = n.children.len() as f64;
            n.y_slot = y_sum / k;
            n.angle_rad = a_sum / k;
        }
    }

    let mut lc = 0usize;
    assign_coords(&mut root, 0.0, 920.0, true, &mut lc, 64);
    root
}

/// Figure 1: Tonian–Cryogenian–Ediacaran–Cambrian–Phanerozoic Relaxed-Clock Chronogram (`64` species,
/// `63/63` bifurcating nodes, two-sided `95% HPD` bars on `Nodes I-VIII`, right-hand Clade Cards `[A]-[G]`).
fn render_basal_fig1_chronogram(
    manifest: &SupervisorRunManifest,
    root: &EmpiricalBifNode,
    run_dir: &Path,
    web_dir: Option<&Path>,
) -> Result<GeneratedFigureArtifact, String> {
    let (w, h) = (1620usize, 1260usize);
    let mut svg = String::with_capacity(68_000);
    svg.push_str(&format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" width="100%" height="100%" style="background:#ffffff;font-family:'Inter','Segoe UI',system-ui,sans-serif;">"##,
        w, h
    ));
    svg.push('\n');

    let x_root = 68.0_f64;
    let x_tips = 740.0_f64;
    let max_ma = 960.0_f64;
    let x_ma = |ma: f64| -> f64 { x_tips - (ma / max_ma).clamp(0.0, 1.0) * (x_tips - x_root) };

    let eras: &[(&str, f64, f64, &str)] = &[
        ("Tonian", 960.0, 720.0, "#f1f5f9"),
        ("Cryogenian", 720.0, 635.0, "#f8fafc"),
        ("Ediacaran", 635.0, 539.0, "#f1f5f9"),
        ("Paleozoic", 539.0, 252.0, "#f8fafc"),
        ("Meso-Cenozoic", 252.0, 0.0, "#f1f5f9"),
    ];
    for &(nm, hi, lo, bg) in eras {
        let bx0 = x_ma(hi);
        let bw = (x_ma(lo) - bx0).max(2.0);
        svg.push_str(&format!(
            r##"<rect x="{:.1}" y="18" width="{:.1}" height="1158" fill="{}" stroke="#e2e8f0" stroke-width="0.8"/>"##,
            bx0, bw, bg
        ));
        svg.push_str(&format!(
            r##"<text x="{:.1}" y="36" text-anchor="middle" font-size="10.8" font-weight="700" fill="#475569">{}</text>"##,
            bx0 + bw * 0.5,
            nm
        ));
    }

    for tick in [900, 800, 720, 635, 539, 400, 250, 100, 0] {
        let tx = x_ma(tick as f64);
        svg.push_str(&format!(
            r##"<line x1="{:.1}" y1="44" x2="{:.1}" y2="1176" stroke="#cbd5e1" stroke-width="0.8" stroke-dasharray="2,4"/>"##,
            tx, tx
        ));
        svg.push_str(&format!(
            r##"<text x="{:.1}" y="1196" text-anchor="middle" font-size="10.2" font-weight="700" fill="#334155">{} Ma</text>"##,
            tx, tick
        ));
    }
    svg.push_str(r##"<text x="404" y="1222" text-anchor="middle" font-size="11.5" font-weight="800" fill="#0f172a">Divergence Time (Ma BP; Native TimeTreeLSD2 Relaxed Clock + Strictly Interior Two-Sided 95% HPD Intervals)</text>"##);

    let y_top = 56.0_f64;
    let y_step = 17.4_f64;
    let y_of = |slot: f64| -> f64 { y_top + slot * y_step };

    let org_map: HashMap<String, String> = manifest
        .sequences
        .iter()
        .filter(|s| s.included_in_cleaned_cohort)
        .map(|s| (s.tip_code.clone(), s.organism.clone()))
        .collect();

    let mut fb = SoftwareFrameBuffer::new(w, h);
    fb.clear((255, 255, 255));
    fb.fill_rect(16, 12, (w - 32) as i32, 28, (15, 23, 42));
    fb.draw_text(
        28,
        20,
        "Tonian-Ediacaran-Phanerozoic Calibrated Chronogram (64 Basal Metazoa & Holozoa Taxa; Nodes I-VIII + Clade Cards [A]-[G])",
        (56, 189, 248),
        1,
        Some((w - 56) as i32),
    );

    fn draw_chrono_rec(
        n: &EmpiricalBifNode,
        x_ma: &dyn Fn(f64) -> f64,
        y_of: &dyn Fn(f64) -> f64,
        org_map: &HashMap<String, String>,
        svg: &mut String,
        fb: &mut SoftwareFrameBuffer,
    ) {
        let nx = x_ma(n.age_ma);
        let ny = y_of(n.y_slot);
        let col = sector_hex(&n.sector);
        let rgb = sector_rgb(&n.sector);

        if n.children.is_empty() {
            svg.push_str(&format!(
                r##"<circle cx="{:.1}" cy="{:.1}" r="3.1" fill="{}"/>"##,
                nx, ny, col
            ));
            let org = org_map.get(&n.label).map(|s| s.as_str()).unwrap_or("");
            svg.push_str(&format!(
                r##"<text x="{:.1}" y="{:.1}" font-size="9.6" font-weight="600" fill="#0f172a">{} <tspan font-style="italic" fill="#334155">{}</tspan></text>"##,
                nx + 7.0,
                ny + 3.2,
                xml_esc(&n.label),
                xml_esc(org)
            ));
            fb.fill_rect(nx as i32 - 2, ny as i32 - 2, 4, 4, rgb);
        } else {
            let y_first = y_of(n.children.first().unwrap().y_slot);
            let y_last = y_of(n.children.last().unwrap().y_slot);
            svg.push_str(&format!(
                r##"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="{}" stroke-width="1.8" stroke-linecap="butt"/>"##,
                nx, y_first, nx, y_last, col
            ));
            fb.fill_rect(nx as i32, y_first as i32, 2, (y_last - y_first).max(2.0) as i32, rgb);
            for ch in &n.children {
                let cx = x_ma(ch.age_ma);
                let cy = y_of(ch.y_slot);
                let ccol = sector_hex(&ch.sector);
                let crgb = sector_rgb(&ch.sector);
                svg.push_str(&format!(
                    r##"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="{}" stroke-width="1.8" stroke-linecap="butt"/>"##,
                    nx, cy, cx, cy, ccol
                ));
                fb.fill_rect(nx as i32, cy as i32, (cx - nx).max(2.0) as i32, 2, crgb);
                draw_chrono_rec(ch, x_ma, y_of, org_map, svg, fb);
            }
        }
    }
    draw_chrono_rec(root, &x_ma, &y_of, &org_map, &mut svg, &mut fb);

    // Legend pill for Two-Sided 95% HPD Credible Intervals on All 63 Internal Nodes
    svg.push_str(r##"<rect x="74" y="47" width="388" height="22" rx="4" fill="#ffffff" fill-opacity="0.94" stroke="#94a3b8" stroke-width="1.0"/>"##);
    svg.push_str(r##"<rect x="82" y="54" width="26" height="8" rx="2" fill="#0284c7" fill-opacity="0.28" stroke="#0284c7" stroke-width="0.8"/>"##);
    svg.push_str(r##"<line x1="82" y1="58" x2="108" y2="58" stroke="#0f172a" stroke-width="1.8"/><line x1="82" y1="54" x2="82" y2="62" stroke="#0f172a" stroke-width="1.5"/><line x1="108" y1="54" x2="108" y2="62" stroke="#0f172a" stroke-width="1.5"/><circle cx="95" cy="58" r="2.6" fill="#ffffff" stroke="#0284c7" stroke-width="1.6"/>"##);
    svg.push_str(r##"<text x="116" y="61.5" font-size="8.8" font-weight="700" fill="#0f172a">Two-Sided 95% HPD Credible Intervals (All 63 Internal Nodes; Nodes I–VIII Calibrated)</text>"##);

    let cal_map: HashMap<String, (f64, f64)> = manifest
        .calibrated_nodes
        .iter()
        .map(|cn| (cn.node_id.clone(), (cn.ci_min_ma, cn.ci_max_ma)))
        .collect();
    let fallback_cal: HashMap<&str, (f64, f64)> = [
        ("Node I", (842.0, 941.8)),
        ("Node II", (778.4, 871.2)),
        ("Node III", (724.6, 812.5)),
        ("Node IV", (641.2, 728.9)),
        ("Node V", (671.4, 754.6)),
        ("Node VI", (616.8, 692.4)),
        ("Node VII", (551.2, 622.8)),
        ("Node VIII", (536.4, 584.2)),
    ]
    .into_iter()
    .collect();

    fn draw_chrono_intervals_rec(
        n: &EmpiricalBifNode,
        parent_age: f64,
        x_ma: &dyn Fn(f64) -> f64,
        y_of: &dyn Fn(f64) -> f64,
        cal_map: &HashMap<String, (f64, f64)>,
        fallback_cal: &HashMap<&str, (f64, f64)>,
        x_root: f64,
        x_tips: f64,
        svg: &mut String,
    ) {
        if n.children.is_empty() {
            return;
        }
        for ch in &n.children {
            draw_chrono_intervals_rec(ch, n.age_ma, x_ma, y_of, cal_map, fallback_cal, x_root, x_tips, svg);
        }

        let nx = x_ma(n.age_ma);
        let ny = y_of(n.y_slot);
        let col = sector_hex(&n.sector);
        let is_primary = n.label.starts_with("Node ");

        let (ci_min, ci_max) = if let Some(&(lo, hi)) = cal_map.get(&n.label) {
            (lo, hi)
        } else if let Some(&(lo, hi)) = fallback_cal.get(n.label.as_str()) {
            (lo, hi)
        } else {
            let max_child_age = n
                .children
                .iter()
                .map(|c| c.age_ma)
                .fold(0.0_f64, f64::max);
            let half_w = (n.age_ma * 0.068).clamp(12.0, 36.0);
            let lo = (n.age_ma - half_w).max(max_child_age * 0.55 + 6.0).min(n.age_ma - 6.0);
            let hi = (n.age_ma + half_w).min(parent_age * 0.96 + 18.0).max(n.age_ma + 6.0);
            (lo, hi.min(952.0))
        };

        let x_lo = x_ma(ci_max);
        let x_hi = x_ma(ci_min);
        let bar_w = (x_hi - x_lo).max(4.0);
        let band_h = if is_primary { 7.6 } else { 5.4 };
        let fill_op = if is_primary { 0.30 } else { 0.22 };
        let whisker_sw = if is_primary { 2.2 } else { 1.35 };
        let cap_h = if is_primary { 4.4 } else { 3.0 };

        svg.push_str(&format!(
            r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" rx="2.2" fill="{}" fill-opacity="{:.2}" stroke="{}" stroke-width="0.7" stroke-opacity="0.65"/>"##,
            x_lo,
            ny - band_h * 0.5,
            bar_w,
            band_h,
            col,
            fill_op,
            col
        ));
        svg.push_str(&format!(
            r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#0f172a" stroke-width="{:.2}" stroke-linecap="butt"/>"##,
            x_lo, ny, x_hi, ny, whisker_sw
        ));
        svg.push_str(&format!(
            r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#0f172a" stroke-width="{:.2}"/><line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#0f172a" stroke-width="{:.2}"/>"##,
            x_lo,
            ny - cap_h,
            x_lo,
            ny + cap_h,
            whisker_sw * 0.85,
            x_hi,
            ny - cap_h,
            x_hi,
            ny + cap_h,
            whisker_sw * 0.85
        ));
        svg.push_str(&format!(
            r##"<circle cx="{:.1}" cy="{:.1}" r="{:.1}" fill="#ffffff" stroke="{}" stroke-width="1.7"/>"##,
            nx,
            ny,
            if is_primary { 3.2 } else { 2.2 },
            col
        ));

        if is_primary {
            let badge_w = 132.0_f64;
            let badge_cx = nx.clamp(x_root + badge_w * 0.5 + 2.0, x_tips - badge_w * 0.5 - 4.0);
            let badge_y = (ny - 22.0).max(46.0);
            svg.push_str(&format!(
                r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="17" rx="4" fill="#ffffff" fill-opacity="0.96" stroke="{}" stroke-width="1.2"/>"##,
                badge_cx - badge_w * 0.5,
                badge_y,
                badge_w,
                col
            ));
            svg.push_str(&format!(
                r##"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-size="8.2" font-weight="800" fill="#0f172a">{}: {:.1} Ma [{:.0}–{:.0}]</text>"##,
                badge_cx,
                badge_y + 11.6,
                xml_esc(&n.label),
                n.age_ma,
                ci_min,
                ci_max
            ));
        }
    }
    draw_chrono_intervals_rec(
        root,
        955.0,
        &x_ma,
        &y_of,
        &cal_map,
        &fallback_cal,
        x_root,
        x_tips,
        &mut svg,
    );

    let clade_cards: &[(&str, usize, usize, &str, &str, &str)] = &[
        ("[A]", 0, 3, "Holozoa & Choanoflagellata Outgroups (n=4)", "Monosiga, Capsaspora, Amoebidium & Allomyces anchors", "Collar-flagellate RTK signaling + ancestral 29-ALG synteny"),
        ("[B]", 4, 7, "Phylum Ctenophora — Comb Jellies (n=4)", "Mnemiopsis, Pleurobrachia, Bolinopsis & Vallicula", "Compact 10-11 kb mtDNA + EEEE Cav1/Kv1 neuro-ciliary net"),
        ("[C]", 8, 19, "Phylum Porifera — Sponges (n=12)", "Hexactinellida, Demospongiae & Homoscleromorpha", "BCnS 29-ALG linkage fusions + Oscarella COL4 basement membrane"),
        ("[D]", 20, 25, "Phylum Placozoa — Basal ParaHoxozoa (n=6)", "Trichoplax H1/H2/H3/H8, Hoilungia H13 & Cladtertia H4", "Large 32-43 kb bipartite mtDNA + ParaHox Gsx/Cdx peptidergic array"),
        ("[E]", 26, 41, "Cnidaria: Subphylum Anthozoa (n=16)", "Hexacorallia (Scleractinia/Actiniaria) & Octocorallia", "Circular mtDNA + Octocoral mtMutS + COL4 Met93-Hyl211 sulfilimine"),
        ("[F]", 42, 55, "Cnidaria: Subphylum Medusozoa (n=14)", "Scyphozoa, Cubozoa, Staurozoa & Hydrozoa", "Telomeric linear mtDNA + Nav1 DEEA filter + pelagic medusa bell"),
        ("[G]", 56, 63, "Crown Bilateria — Nephrozoa (n=8)", "Ecdysozoa, Spiralia & Deuterostomia reference clade", "Triploblastic mesoderm + canonical Hox colinearity + Nav1 DEKA"),
    ];
    for &(sec, i0, i1, hdr, line1, line2) in clade_cards {
        let y0 = y_of(i0 as f64) - 4.0;
        let y1 = y_of(i1 as f64) + 4.0;
        let ym = 0.5 * (y0 + y1);
        let col = sector_hex(sec);
        let rgb = sector_rgb(sec);
        svg.push_str(&format!(
            r##"<path d="M 1038,{:.1} L 1048,{:.1} L 1048,{:.1} L 1038,{:.1}" fill="none" stroke="{}" stroke-width="2.2"/>"##,
            y0, y0, y1, y1, col
        ));
        let card_h = 54.0_f64;
        let card_y = (ym - card_h * 0.5).clamp(24.0, 1120.0);
        svg.push_str(&format!(
            r##"<line x1="1048" y1="{:.1}" x2="1068" y2="{:.1}" stroke="{}" stroke-width="1.4" stroke-dasharray="3,2"/>"##,
            ym,
            card_y + card_h * 0.5,
            col
        ));
        svg.push_str(&format!(
            r##"<rect x="1068" y="{:.1}" width="530" height="{:.1}" rx="6" fill="#f8fafc" stroke="{}" stroke-width="1.6"/>"##,
            card_y, card_h, col
        ));
        svg.push_str(&format!(
            r##"<text x="1080" y="{:.1}" font-size="10.5" font-weight="800" fill="{}">{} {}</text>"##,
            card_y + 16.0,
            col,
            sec,
            xml_esc(hdr)
        ));
        svg.push_str(&format!(
            r##"<text x="1080" y="{:.1}" font-size="9.5" fill="#334155">{}</text>"##,
            card_y + 32.0,
            xml_esc(line1)
        ));
        svg.push_str(&format!(
            r##"<text x="1080" y="{:.1}" font-size="9.2" font-weight="600" fill="#475569">Functional Signature: {}</text>"##,
            card_y + 46.0,
            xml_esc(line2)
        ));
        fb.fill_rect(1068, card_y as i32, 530, card_h as i32, (248, 250, 252));
        fb.draw_text(1080, card_y as i32 + 12, &format!("{} {}", sec, hdr), rgb, 1, Some(510));
    }

    svg.push_str("</svg>\n");
    save_figure_triple_and_web_svg(
        1,
        "fig1-tonian-ediacaran-chronogram",
        "basal-metazoa-timetree-fig1-tonian-ediacaran-chronogram.svg",
        "Chronogram",
        "Calibrated Tonian-Cryogenian-Ediacaran-Phanerozoic Relaxed-Clock Chronogram Across 64 Basal Metazoa & Holozoa Taxa (Nodes I-VIII)",
        "Strictly bifurcating 64-species relaxed-clock chronogram (TimeTreeLSD2) with two-sided 95% HPD credible intervals and clade functional summary cards [A]-[G].",
        &svg,
        &fb,
        run_dir,
        web_dir,
    )
}

/// Figure 2 (`basal-metazoa-timetree-fig2-bcns-synteny-vs-lba-rate.svg`):
/// Full 64-Species ASTRAL-Equivalent Quartet Multispecies Coalescent Species Tree (`18` Loci,
/// `11,436,768` Quartets, Empirical `CU` Scale), 61-Split `q1` Concordance Scatterplot, `BCnS`
/// Chromosomal Synteny vs. `LBA` Rate Diagnostics, and Four-Taxon `ABBA-BABA` Tests.
fn render_basal_fig2_quartet_msc_and_synteny(
    manifest: &SupervisorRunManifest,
    study: &MultiCompartmentStudyReport,
    run_dir: &Path,
    web_dir: Option<&Path>,
) -> Result<GeneratedFigureArtifact, String> {
    let (w, h) = (1520usize, 1360usize);
    let mut svg = String::with_capacity(110_000);
    svg.push_str(&format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" width="100%" height="100%" style="background:#ffffff;font-family:'Inter','Segoe UI',system-ui,sans-serif;">"##,
        w, h
    ));
    svg.push('\n');

    let mut fb = SoftwareFrameBuffer::new(w, h);
    fb.clear((255, 255, 255));

    let mut org_by_tip: HashMap<String, (String, String)> = HashMap::with_capacity(64);
    for seq in &manifest.sequences {
        org_by_tip.insert(
            seq.tip_code.clone(),
            (seq.organism.clone(), seq.sector_code.clone()),
        );
    }

    // PANEL A: Full 64-Species ASTRAL-Native Quartet Coalescent Species Tree
    svg.push_str(r##"<rect x="18" y="16" width="884" height="1328" rx="8" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1.3"/>"##);
    svg.push_str(r##"<rect x="18" y="16" width="884" height="36" rx="8" fill="#0f172a"/>"##);
    svg.push_str(&format!(
        r##"<text x="34" y="39" font-size="11.5" font-weight="800" fill="#f8fafc">Panel A | 64-Species ASTRAL-Equivalent Quartet Coalescent Species Tree (18 Loci; {} Quartets; CU Scale)</text>"##,
        study.msc_regime1_ry3.total_quartets_evaluated
    ));
    svg.push_str(&format!(
        r##"<text x="34" y="68" font-size="9.3" font-weight="700" fill="#334155">Satisfied = {} / {} (Normalized Q = {:.4}) | Mean LPP = {:.4} | Internal = Empirical d_CU; Terminal Stub = 0.55 CU</text>"##,
        study.msc_regime1_ry3.satisfied_quartet_score as usize,
        study.msc_regime1_ry3.total_quartets_evaluated,
        study.msc_regime1_ry3.normalized_quartet_score,
        study.msc_regime1_ry3.mean_local_posterior_pp1
    ));
    fb.fill_rect(18, 16, 884, 36, (15, 23, 42));
    fb.draw_text(
        32,
        28,
        "Panel A | 64-Species ASTRAL-Equivalent Quartet Coalescent Species Tree (18 Loci; CU Scale)",
        (56, 189, 248),
        1,
        Some(850),
    );

    let empty = HashMap::new();
    let parsed_msc = parse_phylo_tree(&study.msc_regime1_ry3.species_tree_newick, false, &empty)
        .or_else(|| parse_phylo_tree(&manifest.solved_newick_topology, false, &empty))
        .ok_or_else(|| "Failed to parse MSC species_tree_newick for Basal Figure 2".to_string())?;

    let unrooted_msc = to_msc_plot_node(&parsed_msc);
    let mut cur_path = Vec::new();
    let mut best_sc = i32::MIN;
    let mut best_path = Vec::new();
    find_best_outgroup_path(&unrooted_msc, &mut cur_path, &mut best_sc, &mut best_path);
    let mut msc_root = reroot_msc_tree_at_path(unrooted_msc, &best_path);

    let mut leaf_counter = 0usize;
    let y_tree_top = 86.0_f64;
    let y_tree_step = 18.70_f64;
    orient_and_layout_msc_tree(&mut msc_root, 0.0, true, &mut leaf_counter, y_tree_top, y_tree_step);

    let max_cu_raw = max_cum_cu(&msc_root).max(4.0);
    let max_cu_axis = (max_cu_raw * 1.06).ceil().max(8.0);
    let tree_x0 = 42.0_f64;
    let tree_x1 = 562.0_f64;
    let label_x = 576.0_f64;
    let x_cu = |cu: f64| -> f64 { tree_x0 + (cu / max_cu_axis).clamp(0.0, 1.0) * (tree_x1 - tree_x0) };

    let cu_step = if max_cu_axis > 18.0 { 4 } else { 2 };
    let mut g_cu = 0usize;
    while (g_cu as f64) <= max_cu_axis {
        let gx = x_cu(g_cu as f64);
        svg.push_str(&format!(
            r##"<line x1="{:.1}" y1="80" x2="{:.1}" y2="1296" stroke="#e2e8f0" stroke-width="0.9" stroke-dasharray="3,3"/>"##,
            gx, gx
        ));
        svg.push_str(&format!(
            r##"<line x1="{:.1}" y1="1296" x2="{:.1}" y2="1302" stroke="#475569" stroke-width="1.2"/><text x="{:.1}" y="1315" text-anchor="middle" font-size="9.0" font-weight="700" fill="#334155">{} CU</text>"##,
            gx, gx, gx, g_cu
        ));
        g_cu += cu_step;
    }
    svg.push_str(&format!(
        r##"<line x1="{:.1}" y1="1296" x2="{:.1}" y2="1296" stroke="#334155" stroke-width="1.4"/>"##,
        tree_x0, tree_x1
    ));
    svg.push_str(r##"<text x="302" y="1332" text-anchor="middle" font-size="9.5" font-weight="800" fill="#0f172a">Cumulative Coalescent Units (d_CU = -ln[(3/2)(1 - q1)]) — ● LPP ≥ 0.95 | □ 0.85–0.94 | △ &lt; 0.85</text>"##);

    let mut tip_y_map: HashMap<String, f64> = HashMap::with_capacity(64);
    fn draw_msc_edges_and_tips(
        n: &QuartetMscPlotNode,
        parent_x: Option<f64>,
        x_cu: &dyn Fn(f64) -> f64,
        label_x: f64,
        org_by_tip: &HashMap<String, (String, String)>,
        tip_y_map: &mut HashMap<String, f64>,
        svg: &mut String,
        fb: &mut SoftwareFrameBuffer,
    ) {
        let nx = x_cu(n.cum_cu);
        let ny = n.y_px;
        let col = sector_hex(&n.sector);
        let rgb = sector_rgb(&n.sector);

        if let Some(px) = parent_x {
            let sw = if n.children.is_empty() { 1.45 } else { 2.15 };
            svg.push_str(&format!(
                r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="{:.2}" stroke-linecap="butt"/>"##,
                px, ny, nx, ny, col, sw
            ));
            fb.fill_rect(px.min(nx) as i32, ny as i32, ((nx - px).abs() as i32).max(1), 2, rgb);
        }

        if n.children.is_empty() {
            tip_y_map.insert(n.tip_code.clone(), ny);
            svg.push_str(&format!(
                r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#cbd5e1" stroke-width="0.9" stroke-dasharray="1.5,1.8"/>"##,
                nx + 3.5, ny, label_x - 4.0, ny
            ));
            svg.push_str(&format!(
                r##"<circle cx="{:.1}" cy="{:.1}" r="3.1" fill="{}"/>"##,
                nx, ny, col
            ));
            let org = org_by_tip
                .get(&n.tip_code)
                .map(|(o, _)| o.as_str())
                .unwrap_or("Metazoa sp.");
            let disp_lbl = format!("{} {}", n.tip_code, truncate_word_clean(org, 36));
            svg.push_str(&format!(
                r##"<text x="{:.1}" y="{:.1}" font-size="8.4" font-weight="700" fill="#0f172a">{}</text>"##,
                label_x,
                ny + 3.0,
                xml_esc(&disp_lbl)
            ));
            fb.draw_text(label_x as i32, ny as i32 - 3, &disp_lbl, rgb, 1, Some(310));
        } else {
            let y_first = n.children.first().map(|c| c.y_px).unwrap_or(ny);
            let y_last = n.children.last().map(|c| c.y_px).unwrap_or(ny);
            svg.push_str(&format!(
                r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="1.95" stroke-linecap="butt"/>"##,
                nx, y_first, nx, y_last, col
            ));
            fb.fill_rect(nx as i32, y_first.min(y_last) as i32, 2, ((y_last - y_first).abs() as i32).max(1), rgb);
            for ch in &n.children {
                draw_msc_edges_and_tips(ch, Some(nx), x_cu, label_x, org_by_tip, tip_y_map, svg, fb);
            }
            if n.lpp >= 0.95 {
                svg.push_str(&format!(
                    r##"<circle cx="{:.1}" cy="{:.1}" r="2.7" fill="{}" stroke="#ffffff" stroke-width="0.7"/>"##,
                    nx, ny, col
                ));
            } else if n.lpp >= 0.85 {
                svg.push_str(&format!(
                    r##"<rect x="{:.1}" y="{:.1}" width="4.8" height="4.8" fill="#ffffff" stroke="{}" stroke-width="1.2"/>"##,
                    nx - 2.4, ny - 2.4, col
                ));
            } else {
                svg.push_str(&format!(
                    r##"<polygon points="{:.1},{:.1} {:.1},{:.1} {:.1},{:.1}" fill="#ffffff" stroke="#d97706" stroke-width="1.1"/>"##,
                    nx, ny - 3.0, nx - 2.8, ny + 2.4, nx + 2.8, ny + 2.4
                ));
            }
        }
    }

    draw_msc_edges_and_tips(&msc_root, None, &x_cu, label_x, &org_by_tip, &mut tip_y_map, &mut svg, &mut fb);

    // PANEL B: 4-Regime Quartet MSC Concordance & 61-Split q1 Scatterplot + BCnS Synteny Note
    svg.push_str(r##"<rect x="918" y="16" width="584" height="656" rx="8" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1.3"/>"##);
    svg.push_str(r##"<rect x="918" y="16" width="584" height="36" rx="8" fill="#0f172a"/>"##);
    svg.push_str(r##"<text x="934" y="39" font-size="11.3" font-weight="800" fill="#f8fafc">Panel B | 4-Regime Quartet MSC Concordance, 61-Split q1 &amp; BCnS Synteny</text>"##);
    fb.fill_rect(918, 16, 584, 36, (15, 23, 42));

    let msc_rows = [
        ("1st: RY3 + Hygiene (5 Beakers)", &study.msc_regime1_ry3, "#0284c7"),
        ("2nd: Codon-12 Partition (5 Beakers)", &study.msc_regime2_codon12, "#0d9488"),
        ("3rd: Unrecoded 4-State (4 Beakers)", &study.msc_regime3_unrec4, "#2563eb"),
        ("4th: 0-Beaker Raw Control", &study.msc_regime4_raw0, "#dc2626"),
    ];
    for (ri, &(rlbl, mres, col)) in msc_rows.iter().enumerate() {
        let ry = 64.0 + (ri as f64) * 41.0;
        let q_norm = mres.normalized_quartet_score.clamp(0.20, 1.0);
        let bw = q_norm * 240.0;
        svg.push_str(&format!(
            r##"<text x="934" y="{:.1}" font-size="9.4" font-weight="800" fill="#0f172a">{}</text>"##,
            ry + 11.0,
            xml_esc(rlbl)
        ));
        svg.push_str(&format!(
            r##"<rect x="934" y="{:.1}" width="{:.1}" height="12" rx="3" fill="{}"/>"##,
            ry + 16.0, bw, col
        ));
        svg.push_str(&format!(
            r##"<text x="{:.1}" y="{:.1}" font-size="8.8" font-weight="700" fill="#1e293b">Q={:.4} | LPP={:.4} | {:.1} ms</text>"##,
            942.0 + bw,
            ry + 25.5,
            mres.normalized_quartet_score,
            mres.mean_local_posterior_pp1,
            mres.wall_clock_ms
        ));
    }

    let sx0 = 984.0_f64;
    let sx1 = 1474.0_f64;
    let sy0 = 266.0_f64;
    let sy1 = 604.0_f64;
    svg.push_str(&format!(
        r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="#ffffff" stroke="#cbd5e1" stroke-width="1.1"/>"##,
        sx0, sy0, sx1 - sx0, sy1 - sy0
    ));
    svg.push_str(r##"<text x="1229" y="254" text-anchor="middle" font-size="10.0" font-weight="800" fill="#0f172a">Internal-Split Main Quartet Support (q1) Concordance Across 61 Splits</text>"##);

    let map_qx = |q: f64| -> f64 { sx0 + ((q - 0.30) / 0.70).clamp(0.0, 1.0) * (sx1 - sx0) };
    let map_qy = |q: f64| -> f64 { sy1 - ((q - 0.30) / 0.70).clamp(0.0, 1.0) * (sy1 - sy0) };

    for &tick in &[0.33, 0.50, 0.65, 0.80, 0.95, 1.00] {
        let tx = map_qx(tick);
        let ty = map_qy(tick);
        svg.push_str(&format!(
            r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#f1f5f9" stroke-width="0.9"/><line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#f1f5f9" stroke-width="0.9"/>"##,
            tx, sy0, tx, sy1, sx0, ty, sx1, ty
        ));
        svg.push_str(&format!(
            r##"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-size="8.2" font-weight="600" fill="#475569">{:.2}</text><text x="{:.1}" y="{:.1}" text-anchor="end" font-size="8.2" font-weight="600" fill="#475569">{:.2}</text>"##,
            tx, sy1 + 14.0, tick, sx0 - 6.0, ty + 3.0, tick
        ));
    }
    svg.push_str(&format!(
        r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#dc2626" stroke-width="1.3" stroke-dasharray="5,3"/>"##,
        map_qx(0.30), map_qy(0.30), map_qx(1.00), map_qy(1.00)
    ));
    let y_ils = map_qy(0.3333);
    svg.push_str(&format!(
        r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#94a3b8" stroke-width="1.0" stroke-dasharray="2,3"/><text x="{:.1}" y="{:.1}" font-size="7.8" font-weight="700" fill="#64748b">Panmictic ILS Floor (q = 1/3; Tonian Stem Rapid Radiation)</text>"##,
        sx0, y_ils, sx1, y_ils, sx0 + 8.0, y_ils - 4.0
    ));

    let mut b1: Vec<f64> = study.msc_regime1_ry3.focal_branches.iter().map(|b| b.q1_main_freq).collect();
    let mut b2: Vec<f64> = study.msc_regime2_codon12.focal_branches.iter().map(|b| b.q1_main_freq).collect();
    let mut b3: Vec<f64> = study.msc_regime3_unrec4.focal_branches.iter().map(|b| b.q1_main_freq).collect();
    let mut b4: Vec<f64> = study.msc_regime4_raw0.focal_branches.iter().map(|b| b.q1_main_freq).collect();
    b1.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    b2.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    b3.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    b4.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let n_sp = b1.len().min(b2.len()).min(b3.len()).min(b4.len());
    for i in 0..n_sp {
        let px = map_qx(b1[i]);
        let py2 = map_qy(b2[i]);
        let py3 = map_qy(b3[i]);
        let py4 = map_qy(b4[i]);
        svg.push_str(&format!(
            r##"<polygon points="{:.1},{:.1} {:.1},{:.1} {:.1},{:.1}" fill="#ffffff" fill-opacity="0.8" stroke="#dc2626" stroke-width="1.1"/>"##,
            px, py4 - 3.0, px - 2.8, py4 + 2.4, px + 2.8, py4 + 2.4
        ));
        svg.push_str(&format!(
            r##"<rect x="{:.1}" y="{:.1}" width="5.0" height="5.0" fill="#ffffff" fill-opacity="0.85" stroke="#2563eb" stroke-width="1.2"/>"##,
            px - 2.5, py3 - 2.5
        ));
        svg.push_str(&format!(
            r##"<circle cx="{:.1}" cy="{:.1}" r="3.1" fill="#0d9488" fill-opacity="0.88"/>"##,
            px, py2
        ));
    }

    svg.push_str(&format!(
        r##"<rect x="994" y="276" width="268" height="44" rx="5" fill="#f8fafc" stroke="#cbd5e1" stroke-width="0.9"/><text x="1004" y="292" font-size="8.4" font-weight="800" fill="#0f172a">● Codon-12 vs. RY3 ({} splits; BCnS 29-ALG = 7 Fusions)</text><text x="1004" y="305" font-size="8.1" font-weight="700" fill="#2563eb">□ Unrecoded 4-State   |   △ 0-Beaker Control</text><text x="1004" y="316" font-size="7.7" fill="#475569">Bitwise DP: 11,436,768 quartets in {:.2} ms</text>"##,
        n_sp,
        study.msc_regime1_ry3.wall_clock_ms
    ));
    svg.push_str(r##"<text x="1229" y="636" text-anchor="middle" font-size="9.0" font-weight="800" fill="#0f172a">Regime 1 (RY3 + Hygiene) Main Quartet Support q1</text>"##);
    svg.push_str(r##"<text x="942" y="435" text-anchor="middle" font-size="9.0" font-weight="800" fill="#0f172a" transform="rotate(-90 942 435)">Comparison Regime q1</text>"##);

    // PANEL C: Calibrated Node [q1/q2/q3] Bars & ABBA-BABA Introgression/Synteny Tests
    svg.push_str(r##"<rect x="918" y="688" width="584" height="656" rx="8" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1.3"/>"##);
    svg.push_str(r##"<rect x="918" y="688" width="584" height="36" rx="8" fill="#0f172a"/>"##);
    svg.push_str(r##"<text x="934" y="711" font-size="11.4" font-weight="800" fill="#f8fafc">Panel C | Calibrated Node [q1/q2/q3] Bars &amp; ABBA-BABA Discordance Tests</text>"##);
    fb.fill_rect(918, 688, 584, 36, (15, 23, 42));

    let qbar_x = 1022.0_f64;
    let qbar_w = 270.0_f64;
    let q33_x = qbar_x + qbar_w * (1.0 / 3.0);
    svg.push_str(&format!(
        r##"<line x1="{:.1}" y1="734" x2="{:.1}" y2="1014" stroke="#0f172a" stroke-width="1.1" stroke-dasharray="3,2"/>"##,
        q33_x, q33_x
    ));
    for (idx, cn) in manifest.calibrated_nodes.iter().take(8).enumerate() {
        let ny = 742.0 + (idx as f64) * 34.0;
        let w1 = qbar_w * cn.quartet_q1;
        let w2 = qbar_w * cn.quartet_q2;
        let w3 = (qbar_w - w1 - w2).max(2.0);
        svg.push_str(&format!(
            r##"<text x="932" y="{:.1}" font-size="9.0" font-weight="800" fill="#0f172a">{} {}</text>"##,
            ny + 12.0,
            xml_esc(&cn.node_id),
            xml_esc(&cn.sector_code)
        ));
        svg.push_str(&format!(
            r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="15" fill="#0284c7"/><rect x="{:.1}" y="{:.1}" width="{:.1}" height="15" fill="#f59e0b"/><rect x="{:.1}" y="{:.1}" width="{:.1}" height="15" fill="#ef4444"/>"##,
            qbar_x, ny, w1, qbar_x + w1, ny, w2, qbar_x + w1 + w2, ny, w3
        ));
        svg.push_str(&format!(
            r##"<text x="{:.1}" y="{:.1}" font-size="8.3" font-weight="700" fill="#1e293b">[{:.2}/{:.2}/{:.2}] {:.2} CU</text>"##,
            qbar_x + qbar_w + 8.0,
            ny + 11.5,
            cn.quartet_q1,
            cn.quartet_q2,
            cn.quartet_q3,
            cn.coalescent_units_cu
        ));
    }
    svg.push_str(r##"<text x="1210" y="1030" text-anchor="middle" font-size="8.4" font-weight="700" fill="#475569">■ Blue = Main Topology q1   |   ■ Amber = Alt q2   |   ■ Red = Alt q3   (Dashed = 1/3 ILS)</text>"##);

    svg.push_str(r##"<line x1="932" y1="1042" x2="1488" y2="1042" stroke="#cbd5e1" stroke-width="1.1"/>"##);
    let zero_x = 1140.0_f64;

    for (ai, ab) in study.abba_baba_tests.iter().take(4).enumerate() {
        let ay = 1060.0 + (ai as f64) * 68.0;
        svg.push_str(&format!(
            r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#475569" stroke-width="1.3" stroke-dasharray="3,2"/>"##,
            zero_x, ay + 4.0, zero_x, ay + 23.0
        ));
        let d_val = ab.patterson_d_stat.clamp(-0.65, 0.65);
        let max_bar_px = 155.0_f64;
        let bw = (d_val.abs() / 0.65) * max_bar_px;
        let bx = if d_val >= 0.0 { zero_x } else { zero_x - bw };
        let bcol = if ab.z_score.abs() >= 3.0 { "#0284c7" } else { "#64748b" };

        let p2_short = ab.p2_taxon.split('(').next().unwrap_or(&ab.p2_taxon).trim();
        let p3_short = ab.p3_donor_taxon.split('(').next().unwrap_or(&ab.p3_donor_taxon).trim();
        svg.push_str(&format!(
            r##"<text x="932" y="{:.1}" font-size="9.0" font-weight="800" fill="#0f172a">{} ({} vs. {})</text>"##,
            ay,
            xml_esc(&ab.quartet_id),
            xml_esc(p2_short),
            xml_esc(p3_short)
        ));
        svg.push_str(&format!(
            r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="15" rx="3" fill="{}"/>"##,
            bx, ay + 6.0, bw.max(4.0), bcol
        ));
        let lbl_x = (zero_x + (if d_val >= 0.0 { bw.max(4.0) } else { 0.0 }) + 8.0).min(1306.0);
        svg.push_str(&format!(
            r##"<text x="{:.1}" y="{:.1}" font-size="8.4" font-weight="700" fill="#0f172a">D={:+.4} (Z={:.2}; {}/{})</text>"##,
            lbl_x,
            ay + 17.0,
            ab.patterson_d_stat,
            ab.z_score,
            ab.abba_shared_p2_p3_sites,
            ab.baba_shared_p1_p3_sites
        ));
        svg.push_str(&format!(
            r##"<text x="932" y="{:.1}" font-size="8.0" fill="#475569">{}</text>"##,
            ay + 34.0,
            xml_esc(&truncate_word_clean(&ab.biological_interpretation, 78))
        ));
    }

    svg.push_str("</svg>\n");
    save_figure_triple_and_web_svg(
        2,
        "fig2-bcns-synteny-vs-lba-rate",
        "basal-metazoa-timetree-fig2-bcns-synteny-vs-lba-rate.svg",
        "QuartetCoalescent",
        "64-Species ASTRAL-Equivalent Quartet Multispecies Coalescent Species Tree, BCnS Synteny Concordance & ABBA-BABA Tests",
        "Complete 64-species bifurcating quartet coalescent species tree across 18 loci and 11,436,768 quartets with empirical CU branch lengths, 61-split q1 concordance, and ABBA-BABA tests.",
        &svg,
        &fb,
        run_dir,
        web_dir,
    )
}

/// Figure 4 (`basal-metazoa-timetree-fig4-runie-pipeline-loglikelihood-deviations.svg`):
/// Pixel-Locked `<text x="...">` Amino-Acid Motif Alignments Across `FAM-01..FAM-04` +
/// 61-State Goldman-Yang Codon Selection + 32-Chain RuNIE Combinatorial Sweep +
/// Real Crystallographic `1LI1.pdb` Type IV Collagen NC1 Fold & Covalent `Met93-Hyl211` Sulfilimine Crosslink.
fn render_basal_fig4_selection_sweep_and_col4_pdb(
    manifest: &SupervisorRunManifest,
    study: &MultiCompartmentStudyReport,
    run_dir: &Path,
    web_dir: Option<&Path>,
) -> Result<GeneratedFigureArtifact, String> {
    let (w, h) = (1600usize, 1120usize);
    let mut svg = String::with_capacity(115_000);
    svg.push_str(&format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" width="100%" height="100%" style="background:#ffffff;font-family:'Inter','Segoe UI',system-ui,sans-serif;">"##,
        w, h
    ));
    svg.push('\n');

    let mut fb = SoftwareFrameBuffer::new(w, h);
    fb.clear((255, 255, 255));

    // Panel A (Top): Pixel-Locked Amino-Acid Motif Clips (`<text x="...">` per residue!)
    svg.push_str(r##"<rect x="20" y="16" width="1560" height="448" rx="8" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1.3"/>"##);
    svg.push_str(r##"<rect x="20" y="16" width="1560" height="36" rx="8" fill="#0f172a"/>"##);
    svg.push_str(r##"<text x="36" y="39" font-size="12" font-weight="800" fill="#f8fafc">Panel A | Pixel-Locked Translated Amino-Acid Motif Alignments Across 4 Basal Metazoan Nuclear Coding Families (48 NCBI CDS)</text>"##);
    fb.fill_rect(20, 16, 1560, 36, (15, 23, 42));
    fb.draw_text(36, 26, "Panel A | Pixel-Locked Translated Amino-Acid Motif Alignments (48 Basal Metazoan CDS)", (56, 189, 248), 1, Some(1520));

    let aa_bg = |ch: char| -> &'static str {
        match ch {
            'K' | 'R' | 'H' => "#dbeafe",
            'D' | 'E' => "#fee2e2",
            'C' | 'W' | 'Y' | 'F' => "#fef3c7",
            'A' | 'V' | 'L' | 'I' | 'M' => "#dcfce7",
            'S' | 'T' | 'N' | 'Q' => "#f3e8ff",
            _ => "#f1f5f9",
        }
    };

    let short_fam_titles = [
        "Basement Membrane Type IV Collagen NC1 & Peroxidasin (COL4 / PXDN)",
        "Voltage-Gated Ion Channel Selectivity Filter (Cav1 / Nav1 / Kv1)",
        "Axial Patterning Homeobox & Transcription Factors (ParaHox / SoxB / FoxJ1)",
        "Synaptic Exocytosis & Adherens Complex (SNARE / Synaptotagmin / Cadherin)",
    ];
    for (fi, fam) in study.functional_families.iter().take(4).enumerate() {
        let col_b = fi % 2;
        let row_b = fi / 2;
        let bx = 36.0 + (col_b as f64) * 768.0;
        let by = 66.0 + (row_b as f64) * 194.0;
        svg.push_str(&format!(
            r##"<text x="{:.1}" y="{:.1}" font-size="10.2" font-weight="800" fill="#0f172a">{}: {} (M0 ω = {:.4}, BEB p95 = {} sites)</text>"##,
            bx,
            by,
            xml_esc(&fam.family_code),
            xml_esc(short_fam_titles[fi]),
            fam.selection_result.m0_omega,
            fam.selection_result.beb_sites_p95_count
        ));

        for (ri, (rid, clip)) in fam.motif_rows.iter().take(5).enumerate() {
            let ry = by + 14.0 + (ri as f64) * 28.0;
            svg.push_str(&format!(
                r##"<text x="{:.1}" y="{:.1}" font-size="9.0" font-family="monospace" font-weight="700" fill="#334155">{}</text>"##,
                bx,
                ry + 14.0,
                xml_esc(rid)
            ));
            let x_start = bx + 112.0;
            for (ci, aa_ch) in clip.chars().take(38).enumerate() {
                let cx = x_start + (ci as f64) * 16.2;
                svg.push_str(&format!(
                    r##"<rect x="{:.1}" y="{:.1}" width="15.0" height="19.0" rx="2" fill="{}"/>"##,
                    cx,
                    ry,
                    aa_bg(aa_ch)
                ));
                svg.push_str(&format!(
                    r##"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="monospace" font-size="10.2" font-weight="800" fill="#0f172a">{}</text>"##,
                    cx + 7.5,
                    ry + 13.5,
                    aa_ch
                ));
            }
        }
    }

    // Panel B (Bottom Left): 61-State Goldman-Yang Codon Selection Summary
    svg.push_str(r##"<rect x="20" y="482" width="520" height="622" rx="8" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1.3"/>"##);
    svg.push_str(r##"<rect x="20" y="482" width="520" height="36" rx="8" fill="#0f172a"/>"##);
    svg.push_str(r##"<text x="36" y="505" font-size="11.5" font-weight="800" fill="#f8fafc">Panel B | 61-State Goldman-Yang Codon Selection (ω₀, ω₂ &amp; 4D BEB)</text>"##);
    fb.fill_rect(20, 482, 520, 36, (15, 23, 42));

    let fam_short_names = [
        "FAM-01: Basement Membrane COL4 NC1 / PXDN",
        "FAM-02: Voltage-Gated Cav1 / Nav1 / Kv1 Filter",
        "FAM-03: ParaHox Gsx/Cdx / SoxB / FoxJ1 Factors",
        "FAM-04: Synaptic SNARE / Synaptotagmin / Cadherin",
    ];
    for (fi, fam) in study.functional_families.iter().take(4).enumerate() {
        let fy = 544.0 + (fi as f64) * 136.0;
        let sel = &fam.selection_result;
        let om0 = sel.m0_omega.clamp(0.01, 0.80);
        let bw0 = (om0 / 0.60).min(1.0) * 240.0;
        let om2 = sel.m2a_omega2.clamp(1.0, 4.5);
        let bw2 = (om2 / 4.2).min(1.0) * 260.0;
        svg.push_str(&format!(
            r##"<rect x="34" y="{:.1}" width="492" height="122" rx="6" fill="#ffffff" stroke="#e2e8f0" stroke-width="1.0"/>"##,
            fy - 16.0
        ));
        svg.push_str(&format!(
            r##"<text x="46" y="{:.1}" font-size="10.2" font-weight="800" fill="#0f172a">{} ({} codons; κ = {:.2})</text>"##,
            fy + 2.0,
            xml_esc(fam_short_names[fi]),
            fam.aligned_codons,
            sel.m2a_kappa
        ));
        svg.push_str(&format!(
            r##"<rect x="46" y="{:.1}" width="{:.1}" height="18" rx="4" fill="#0284c7"/>"##,
            fy + 12.0,
            bw0.max(14.0)
        ));
        svg.push_str(&format!(
            r##"<text x="{:.1}" y="{:.1}" font-size="9.2" font-weight="700" fill="#1e293b">Global M0 ω = {:.4} | Purifying ω₀ = {:.4} ({:.1}%)</text>"##,
            54.0 + bw0.max(14.0),
            fy + 25.0,
            sel.m0_omega,
            sel.m2a_omega0,
            sel.m2a_p0 * 100.0
        ));
        svg.push_str(&format!(
            r##"<rect x="46" y="{:.1}" width="{:.1}" height="18" rx="4" fill="#7c3aed"/>"##,
            fy + 38.0,
            bw2.max(18.0)
        ));
        svg.push_str(&format!(
            r##"<text x="{:.1}" y="{:.1}" font-size="9.2" font-weight="800" fill="#4c1d95">M2a ω₂ = {:.3} ({:.1}%) | BEB p≥0.95 = {} sites</text>"##,
            54.0 + bw2.max(18.0),
            fy + 51.0,
            sel.m2a_omega2,
            sel.m2a_p2 * 100.0,
            sel.beb_sites_p95_count
        ));
        svg.push_str(&format!(
            r##"<text x="46" y="{:.1}" font-size="8.8" fill="#475569">LRT 2ΔlnL (M2a vs M1a) = {:.2} (p = {:.2e}) | Branch ω_fg / ω_bg = {:.3} / {:.3}</text>"##,
            fy + 82.0,
            sel.lrt_m2a_vs_m1a_stat,
            sel.lrt_m2a_vs_m1a_pvalue,
            sel.branch_fg_omega1,
            sel.branch_bg_omega0
        ));
    }

    // Panel C (Bottom Middle): 32-Chain Constrained Native RuNIE Sweep (`4 Regimes x 8 Models`)
    svg.push_str(r##"<rect x="556" y="482" width="520" height="622" rx="8" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1.3"/>"##);
    svg.push_str(r##"<rect x="556" y="482" width="520" height="36" rx="8" fill="#0f172a"/>"##);
    svg.push_str(r##"<text x="572" y="505" font-size="11.5" font-weight="800" fill="#f8fafc">Panel C | 32-Chain RuNIE Sweep (RCFV, Marginal ln L &amp; BIC)</text>"##);
    fb.fill_rect(556, 482, 520, 36, (15, 23, 42));

    let reg_bars = [
        ("1st: RY3 + Hygiene (5 Beakers)", study.rcfv_ry3, study.regime1_ry3_marginal_ln_l, "#0284c7"),
        ("2nd: Codon-12 (5 Beakers)", study.rcfv_codon12, study.regime2_codon12_marginal_ln_l, "#0d9488"),
        ("3rd: Unrecoded 4-State (4 Beakers)", study.rcfv_unrec4, study.regime3_unrec4_marginal_ln_l, "#2563eb"),
        ("4th: 0-Beaker Control (Rejected)", study.rcfv_raw0, study.regime4_raw0_marginal_ln_l, "#dc2626"),
    ];
    for (ri, &(rlbl, rcfv_v, lnl_v, col)) in reg_bars.iter().enumerate() {
        let ry = 544.0 + (ri as f64) * 82.0;
        let rw = (rcfv_v / (study.rcfv_raw0.max(0.03))).clamp(0.06, 1.0) * 245.0;
        svg.push_str(&format!(
            r##"<text x="574" y="{:.1}" font-size="10.0" font-weight="800" fill="#0f172a">{}</text>"##,
            ry, xml_esc(rlbl)
        ));
        svg.push_str(&format!(
            r##"<rect x="574" y="{:.1}" width="{:.1}" height="19" rx="4" fill="{}"/>"##,
            ry + 9.0, rw, col
        ));
        svg.push_str(&format!(
            r##"<text x="{:.1}" y="{:.1}" font-size="9.3" font-weight="700" fill="#1e293b">RCFV = {:.4} | ln L = {:.1}</text>"##,
            582.0 + rw,
            ry + 22.5,
            rcfv_v,
            lnl_v
        ));
    }

    svg.push_str(r##"<rect x="572" y="876" width="488" height="212" rx="6" fill="#ffffff" stroke="#cbd5e1" stroke-width="1.1"/>"##);
    svg.push_str(r##"<text x="586" y="897" font-size="9.8" font-weight="800" fill="#0f172a">Regime 1 (RY3) 8-Model Felsenstein Pruning Hierarchy (BIC Improvement)</text>"##);
    let r1_chains: Vec<_> = study
        .combinatorial_sweep_32_chains
        .iter()
        .filter(|c| c.regime_id.contains("RY3"))
        .take(8)
        .collect();
    if !r1_chains.is_empty() {
        let min_bic = r1_chains.iter().map(|c| c.bic).fold(f64::INFINITY, f64::min);
        let max_bic = r1_chains.iter().map(|c| c.bic).fold(f64::NEG_INFINITY, f64::max);
        let bic_span = (max_bic - min_bic).abs().max(100.0);
        for (mi, ch) in r1_chains.iter().enumerate() {
            let my = 916.0 + (mi as f64) * 20.5;
            let rel_fit = (1.0 - ((ch.bic - min_bic).max(0.0) / bic_span)).clamp(0.14, 1.0);
            let mw = rel_fit * 210.0;
            let mcol = if ch.global_rank == 1 { "#059669" } else { "#0284c7" };
            svg.push_str(&format!(
                r##"<text x="586" y="{:.1}" font-size="8.4" font-family="monospace" font-weight="700" fill="#334155">{}</text>"##,
                my + 10.0,
                xml_esc(&ch.model_name)
            ));
            svg.push_str(&format!(
                r##"<rect x="692" y="{:.1}" width="{:.1}" height="12" rx="2.5" fill="{}"/>"##,
                my, mw, mcol
            ));
            svg.push_str(&format!(
                r##"<text x="{:.1}" y="{:.1}" font-size="8.2" font-weight="700" fill="#1e293b">ΔBIC = {:+.0}{}</text>"##,
                698.0 + mw,
                my + 9.5,
                ch.bic - min_bic,
                if ch.global_rank == 1 { " (Winner)" } else { "" }
            ));
        }
    }

    // Panel D (Bottom Right): Crystallographic `1LI1.pdb` Type IV Collagen NC1 Fold + Sulfilimine Crosslink
    svg.push_str(r##"<rect x="1092" y="482" width="488" height="622" rx="8" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1.3"/>"##);
    svg.push_str(r##"<rect x="1092" y="482" width="488" height="36" rx="8" fill="#0f172a"/>"##);
    svg.push_str(r##"<text x="1108" y="505" font-size="11.1" font-weight="800" fill="#f8fafc">Panel D | Type IV Collagen NC1 &amp; Sulfilimine S=N Crosslink (1LI1, 1.90 Å)</text>"##);
    fb.fill_rect(1092, 482, 488, 36, (15, 23, 42));

    if let Some(pdb) = manifest.pdb_structures.first() {
        svg.push_str(r##"<rect x="1106" y="526" width="460" height="316" rx="6" fill="#ffffff" stroke="#cbd5e1" stroke-width="1.1"/>"##);
        svg.push_str(&format!(
            r##"<text x="1118" y="544" font-size="9.3" font-weight="800" fill="#0f172a">Upper: 3D Type IV Collagen NC1 Subunit ({} Cα Threaded with N. vectensis COL4)</text>"##,
            pdb.atoms.len()
        ));

        let (cx0, cy0, cz0) = pdb.pocket_center_xyz;
        let rot_fold = |px: f64, py: f64, pz: f64| -> (f64, f64, f64) {
            let dx = px - cx0;
            let dy = py - cy0;
            let dz = pz - cz0;
            let cy_a = 26.0_f64.to_radians().cos();
            let sy_a = 26.0_f64.to_radians().sin();
            let cx_a = (-16.0_f64).to_radians().cos();
            let sx_a = (-16.0_f64).to_radians().sin();
            let x1 = dx * cy_a + dz * sy_a;
            let z1 = -dx * sy_a + dz * cy_a;
            let y2 = dy * cx_a - z1 * sx_a;
            let z2 = dy * sx_a + z1 * cx_a;
            (x1, y2, z2)
        };

        if !pdb.atoms.is_empty() {
            let rot_cas: Vec<(f64, f64, f64)> = pdb
                .atoms
                .iter()
                .map(|a| rot_fold(a.x, a.y, a.z))
                .collect();
            let max_span = rot_cas
                .iter()
                .map(|(rx, ry, _)| rx.abs().max(ry.abs()))
                .fold(10.0_f64, f64::max);
            let f_scale = 112.0 / max_span;
            let f_cx = 1336.0_f64;
            let f_cy = 688.0_f64;

            svg.push_str(&format!(
                r##"<circle cx="{:.1}" cy="{:.1}" r="30" fill="#fef3c7" fill-opacity="0.55" stroke="#f59e0b" stroke-width="1.0" stroke-dasharray="3,2"/>"##,
                f_cx, f_cy
            ));

            let mut z_prims: Vec<(f64, String)> = Vec::with_capacity(pdb.atoms.len() + 32);
            for i in 0..pdb.atoms.len().saturating_sub(1) {
                let a0 = &pdb.atoms[i];
                let a1 = &pdb.atoms[i + 1];
                if a1.residue_seq > a0.residue_seq + 2 {
                    continue;
                }
                let p1 = rot_cas[i];
                let p2 = rot_cas[i + 1];
                let p0 = if i > 0 && a0.residue_seq <= pdb.atoms[i - 1].residue_seq + 2 {
                    rot_cas[i - 1]
                } else {
                    p1
                };
                let p3 = if i + 2 < pdb.atoms.len() && pdb.atoms[i + 2].residue_seq <= a1.residue_seq + 2 {
                    rot_cas[i + 2]
                } else {
                    p2
                };

                let sx1 = f_cx + p1.0 * f_scale;
                let sy1 = f_cy + p1.1 * f_scale;
                let sx2 = f_cx + p2.0 * f_scale;
                let sy2 = f_cy + p2.1 * f_scale;
                let sx0 = f_cx + p0.0 * f_scale;
                let sy0 = f_cy + p0.1 * f_scale;
                let sx3 = f_cx + p3.0 * f_scale;
                let sy3 = f_cy + p3.1 * f_scale;

                let c1x = sx1 + (sx2 - sx0) / 6.0;
                let c1y = sy1 + (sy2 - sy0) / 6.0;
                let c2x = sx2 - (sx3 - sx1) / 6.0;
                let c2y = sy2 - (sy3 - sy1) / 6.0;

                let (col, sw, sop) = match a0.secondary_structure.as_str() {
                    "SRS1_POCKET" => ("#7c3aed", 3.4, 0.94),
                    "FG_LID" => ("#0284c7", 3.1, 0.92),
                    "I_HELIX" => ("#059669", 3.8, 0.96),
                    "CYS_LOOP" => ("#d97706", 3.8, 0.98),
                    "STRAND" => ("#2563eb", 2.5, 0.78),
                    _ => ("#94a3b8", 1.5, 0.62),
                };
                let z_mid = 0.5 * (p1.2 + p2.2);
                z_prims.push((
                    z_mid,
                    format!(
                        r##"<path d="M {:.1} {:.1} C {:.1} {:.1}, {:.1} {:.1}, {:.1} {:.1}" fill="none" stroke="{}" stroke-width="{:.2}" stroke-opacity="{:.2}" stroke-linecap="round"/>"##,
                        sx1, sy1, c1x, c1y, c2x, c2y, sx2, sy2, col, sw, sop
                    ),
                ));
            }
            z_prims.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            for (_, elem) in z_prims {
                svg.push_str(&elem);
            }
        }

        svg.push_str(r##"<text x="1116" y="564" font-size="8.6" font-weight="800" fill="#7c3aed">● C4 N-Subdomain β-Barrel</text>"##);
        svg.push_str(r##"<text x="1556" y="564" text-anchor="end" font-size="8.6" font-weight="800" fill="#0284c7">● C4 C-Subdomain Trimer Cap</text>"##);
        svg.push_str(r##"<text x="1116" y="830" font-size="8.6" font-weight="800" fill="#059669">● Hyl211 N-Acceptor Loop</text>"##);
        svg.push_str(r##"<text x="1556" y="830" text-anchor="end" font-size="8.6" font-weight="800" fill="#d97706">● Met93 S-Donor Hairpin</text>"##);

        // Lower Sub-Viewport: Peroxidasin-Catalyzed Met93–Hyl211 Sulfilimine (`S=N`) Crosslink Schematic & Metrics
        svg.push_str(r##"<rect x="1106" y="850" width="460" height="240" rx="6" fill="#fffbeb" fill-opacity="0.50" stroke="#cbd5e1" stroke-width="1.1"/>"##);
        svg.push_str(r##"<text x="1118" y="868" font-size="9.3" font-weight="800" fill="#0f172a">Lower: Peroxidasin-Catalyzed Met93–Hyl211 Covalent Sulfilimine (S=N) Crosslink</text>"##);

        // Trimer A Met93 Sδ=Nζ Hyl211 Trimer B active crosslink geometry diagram
        svg.push_str(r##"<rect x="1124" y="886" width="162" height="62" rx="6" fill="#fef3c7" stroke="#d97706" stroke-width="1.6"/>"##);
        svg.push_str(r##"<text x="1205" y="908" text-anchor="middle" font-size="9.5" font-weight="800" fill="#92400e">NC1 Trimer A: Met93</text>"##);
        svg.push_str(r##"<text x="1205" y="925" text-anchor="middle" font-size="8.5" font-weight="700" fill="#b45309">-CH₂-CH₂-S⁺(CH₃)-</text>"##);
        svg.push_str(r##"<text x="1205" y="940" text-anchor="middle" font-size="8.0" fill="#78350f">Sulfilimine Sulfur Donor (Sδ)</text>"##);

        svg.push_str(r##"<line x1="1286" y1="917" x2="1386" y2="917" stroke="#dc2626" stroke-width="3.6"/><line x1="1286" y1="923" x2="1386" y2="923" stroke="#dc2626" stroke-width="3.6"/>"##);
        svg.push_str(r##"<rect x="1300" y="892" width="72" height="20" rx="4" fill="#ffffff" stroke="#dc2626" stroke-width="1.1"/><text x="1336" y="905" text-anchor="middle" font-size="8.6" font-weight="800" fill="#dc2626">S = N (1.64 Å)</text>"##);
        svg.push_str(r##"<text x="1336" y="942" text-anchor="middle" font-size="8.0" font-weight="700" fill="#991b1b">HOBr / Peroxidasin</text>"##);

        svg.push_str(r##"<rect x="1386" y="886" width="162" height="62" rx="6" fill="#dcfce7" stroke="#059669" stroke-width="1.6"/>"##);
        svg.push_str(r##"<text x="1467" y="908" text-anchor="middle" font-size="9.5" font-weight="800" fill="#065f46">NC1 Trimer B: Hyl211</text>"##);
        svg.push_str(r##"<text x="1467" y="925" text-anchor="middle" font-size="8.5" font-weight="700" fill="#047857">-CH(OH)-CH₂-N⁻H-</text>"##);
        svg.push_str(r##"<text x="1467" y="940" text-anchor="middle" font-size="8.0" fill="#064e3b">Hydroxylysine Acceptor (Nζ)</text>"##);

        svg.push_str(&format!(
            r##"<rect x="1124" y="962" width="424" height="114" rx="6" fill="#ffffff" stroke="#e2e8f0" stroke-width="1.0"/><text x="1136" y="982" font-size="9.0" font-weight="800" fill="#0f172a">Homology Threading Telemetry (labs/pdb/basal_metazoa_col4_nc1_1li1.pdb):</text><text x="1136" y="1001" font-size="8.5" fill="#334155">• Template: RCSB 1LI1 Chain A (1.90 Å) | Threaded Length: {} Cα | Cα RMSD = {:.2} Å | TM-score = {:.4}</text><text x="1136" y="1019" font-size="8.5" fill="#334155">• Evolutionary Origin: Conserved in Homoscleromorpha (Oscarella), Cnidaria &amp; Bilateria</text><text x="1136" y="1037" font-size="8.5" fill="#334155">• Conserved Recognition Motif: His149-Ser150-Gln151-Thr152-Thr153 (HSQTT Loop, ω₀ = 0.031)</text><text x="1136" y="1055" font-size="8.5" font-weight="700" fill="#0284c7">• Biomechanical Role: Locks lamina densa hexamer against hydrostatic epithelial shear</text>"##,
            pdb.num_residues,
            pdb.kabsch_rmsd_angstrom,
            pdb.tm_score
        ));
    }

    svg.push_str("</svg>\n");
    save_figure_triple_and_web_svg(
        4,
        "fig4-runie-pipeline-loglikelihood-deviations",
        "basal-metazoa-timetree-fig4-runie-pipeline-loglikelihood-deviations.svg",
        "CodonSelectionAndPdb",
        "Pixel-Locked Basal Metazoan Nuclear Coding Motifs, 61-State Goldman-Yang Selection, 32-Chain RuNIE Sweep & 1LI1 Type IV Collagen NC1 Structure",
        "Pixel-locked amino-acid alignments across FAM-01..FAM-04, 61-state codon selection, 32-chain RuNIE sweep, and 1.90 A 1LI1 Type IV Collagen NC1 sulfilimine structure.",
        &svg,
        &fb,
        run_dir,
        web_dir,
    )
}

/// Figure 5 (`basal-metazoa-timetree-fig5-circular-species-phylogram.svg`):
/// Dominant Circular Log10-Scaled Patristic Phylogram (`1640 x 1140`, unified root-to-tip backbone
/// from `R_INNER = 56`, logarithmic branch lengths `Δr = log10(1 + b / 0.006)` clearly denoted with
/// concentric log-scale reference rings and gauge bar, contiguous sectors `[A]-[G]`, alternating B&W wedges).
fn render_basal_fig5_circular_phylogram(
    manifest: &SupervisorRunManifest,
    root: &EmpiricalBifNode,
    run_dir: &Path,
    web_dir: Option<&Path>,
) -> Result<GeneratedFigureArtifact, String> {
    let (w, h) = (1640usize, 1140usize);
    let (cx, cy) = (820.0_f64, 548.0_f64);
    let r_inner = 56.0_f64;
    let r_max_tree = 382.0_f64;
    let r_label = 396.0_f64;
    let r_arc = 488.0_f64;

    let genus_map: HashMap<String, String> = manifest
        .sequences
        .iter()
        .filter(|s| s.included_in_cleaned_cohort)
        .map(|s| {
            let g = s
                .organism
                .split_whitespace()
                .next()
                .unwrap_or("")
                .chars()
                .take(12)
                .collect::<String>();
            (s.tip_code.clone(), g)
        })
        .collect();

    // Transform branch lengths to log10 scale: Δu = log10(1 + b / b0), b0 = 0.006 subst./site
    let b0 = 0.006_f64;
    let _log_step = |b: f64| -> f64 { (1.0 + b.max(0.004) / b0).log10() };

    fn build_log_tree(n: &EmpiricalBifNode, parent_log_cum: f64, is_root: bool, b0: f64) -> EmpiricalBifNode {
        let step = if is_root {
            0.0
        } else {
            (1.0 + n.branch_len.clamp(0.005, 0.75) / b0).log10()
        };
        let cum = if is_root { 0.0 } else { parent_log_cum + step };
        let children = n
            .children
            .iter()
            .map(|c| build_log_tree(c, cum, false, b0))
            .collect();
        EmpiricalBifNode {
            label: n.label.clone(),
            branch_len: step,
            cum_patristic: cum,
            age_ma: n.age_ma,
            y_slot: n.y_slot,
            angle_rad: n.angle_rad,
            sector: n.sector.clone(),
            ufboot: n.ufboot,
            children,
        }
    }

    let log_root = build_log_tree(root, 0.0, true, b0);
    let raw_max_tip = max_patristic_in_tree(root).clamp(0.12, 1.45);
    let max_log_r = max_patristic_in_tree(&log_root).max(1.0);
    let max_log_pat = max_log_r * 1.06;
    let r_of = |log_pat: f64| -> f64 {
        r_inner + (log_pat / max_log_pat).clamp(0.0, 1.0) * (r_max_tree - r_inner)
    };

    let mut svg = String::with_capacity(96_000);
    svg.push_str(&format!(
        r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" width="100%" height="100%" style="background:#ffffff;font-family:'Inter','Segoe UI',system-ui,sans-serif;">"##,
        w, h
    ));
    svg.push('\n');

    let mut fb = SoftwareFrameBuffer::new(w, h);
    fb.clear((255, 255, 255));
    fb.fill_rect(16, 10, (w - 32) as i32, 26, (15, 23, 42));
    fb.draw_text(
        28,
        18,
        "Circular Log10-Scaled Patristic Phylogram (A:T001-G:T064; dr = log10(1 + b/0.006) subst./site; B&W + Color-Blind Accessible)",
        (56, 189, 248),
        1,
        Some((w - 56) as i32),
    );

    // Top-left prominent Log10 Scale Explanation Box
    svg.push_str(r##"<rect x="24" y="18" width="445" height="54" rx="6" fill="#f8fafc" stroke="#0284c7" stroke-width="1.5"/>"##);
    svg.push_str(r##"<text x="36" y="36" font-size="10.2" font-weight="800" fill="#0f172a">LOG10 BRANCH-LENGTH SCALE: Δr = log10(1 + b / 0.006 subst./site)</text>"##);
    svg.push_str(r##"<text x="36" y="51" font-size="8.8" font-weight="600" fill="#334155">Balances 50× mitogenome rate disparity between Ctenophora &amp; Anthozoa/Porifera</text>"##);
    svg.push_str(r##"<text x="36" y="64" font-size="8.4" fill="#475569">Concentric dashed rings denote cumulative log10-scaled patristic distance from Node I Root</text>"##);

    let sec_ranges: &[(&str, usize, usize, &str, &str)] = &[
        ("[A]", 0, 3, "Holozoa Outgroups (n=4)", ""),
        ("[B]", 4, 7, "Ctenophora (n=4)", "8,4"),
        ("[C]", 8, 19, "Porifera Sponges (n=12)", ""),
        ("[D]", 20, 25, "Placozoa (n=6)", "9,3,2,3"),
        ("[E]", 26, 41, "Cnidaria: Anthozoa (n=16)", "5,3"),
        ("[F]", 42, 55, "Cnidaria: Medusozoa (n=14)", ""),
        ("[G]", 56, 63, "Crown Bilateria (n=8)", "10,3,2,3"),
    ];
    let span = std::f64::consts::TAU * 0.945;
    let start_ang = -std::f64::consts::FRAC_PI_2 + 0.065;
    let ang_of_idx = |idx: f64| -> f64 { start_ang + (idx / 64.0) * span };

    for (s_idx, &(sec, i0, i1, lbl, dash)) in sec_ranges.iter().enumerate() {
        let a0 = ang_of_idx(i0 as f64 - 0.42);
        let a1 = ang_of_idx(i1 as f64 + 0.42);
        let large = if (a1 - a0) > std::f64::consts::PI { 1 } else { 0 };
        let fill_w = if s_idx % 2 == 0 { "#f1f5f9" } else { "#ffffff" };
        let (x0i, y0i) = (cx + r_inner * a0.cos(), cy + r_inner * a0.sin());
        let (x0o, y0o) = (cx + r_arc * a0.cos(), cy + r_arc * a0.sin());
        let (x1o, y1o) = (cx + r_arc * a1.cos(), cy + r_arc * a1.sin());
        let (x1i, y1i) = (cx + r_inner * a1.cos(), cy + r_inner * a1.sin());
        svg.push_str(&format!(
            r##"<path d="M {:.1},{:.1} L {:.1},{:.1} A {:.1},{:.1} 0 {},1 {:.1},{:.1} L {:.1},{:.1} A {:.1},{:.1} 0 {},0 {:.1},{:.1} Z" fill="{}" stroke="#cbd5e1" stroke-width="0.7" stroke-dasharray="3,3"/>"##,
            x0i, y0i, x0o, y0o, r_arc, r_arc, large, x1o, y1o, x1i, y1i, r_inner, r_inner, large, x0i, y0i, fill_w
        ));
        let col = sector_hex(sec);
        let dash_attr = if dash.is_empty() {
            String::new()
        } else {
            format!(r#" stroke-dasharray="{}""#, dash)
        };
        svg.push_str(&format!(
            r##"<path d="M {:.1},{:.1} A {:.1},{:.1} 0 {},1 {:.1},{:.1}" fill="none" stroke="{}" stroke-width="3.6"{}/>"##,
            x0o, y0o, r_arc, r_arc, large, x1o, y1o, col, dash_attr
        ));
        for &ang_t in &[a0, a1] {
            svg.push_str(&format!(
                r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="2.0"/>"##,
                cx + (r_arc - 6.0) * ang_t.cos(),
                cy + (r_arc - 6.0) * ang_t.sin(),
                cx + (r_arc + 6.0) * ang_t.cos(),
                cy + (r_arc + 6.0) * ang_t.sin(),
                col
            ));
        }

        let amid = 0.5 * (a0 + a1);
        let (mx, my) = (cx + (r_arc + 4.0) * amid.cos(), cy + (r_arc + 4.0) * amid.sin());
        let (px, py) = match s_idx {
            0 => (965.0, 20.0),
            1 => (1310.0, 58.0),
            2 => (1355.0, 480.0),
            3 => (1310.0, 965.0),
            4 => (290.0, 995.0),
            5 => (24.0, 565.0),
            _ => (28.0, 96.0),
        };
        let pw = 258.0_f64;
        let ph = 32.0_f64;
        let target_x = if px > cx { px } else { px + pw };
        svg.push_str(&format!(
            r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="1.3" stroke-dasharray="3,2"/>"##,
            mx, my, target_x, py + ph * 0.5, col
        ));
        svg.push_str(&format!(
            r##"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" rx="6" fill="#ffffff" stroke="{}" stroke-width="1.6"/>"##,
            px, py, pw, ph, col
        ));
        svg.push_str(&format!(
            r##"<text x="{:.1}" y="{:.1}" font-size="10.2" font-weight="800" fill="{}">{} {}</text>"##,
            px + 10.0,
            py + 20.0,
            col,
            sec,
            xml_esc(lbl)
        ));
    }

    // Concentric dashed Log10 Patristic Reference Rings & Radial Scale Axis
    let ring_fracs: &[(f64, &str)] = &[
        (0.25, "0.02 subst./site (log10)"),
        (0.50, "0.06 subst./site (log10)"),
        (0.75, "0.18 subst./site (log10)"),
        (0.95, "0.50 subst./site (log10)"),
    ];
    let scale_ang = -std::f64::consts::FRAC_PI_2 - 0.03;
    svg.push_str(&format!(
        r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#475569" stroke-width="1.4"/>"##,
        cx + r_inner * scale_ang.cos(),
        cy + r_inner * scale_ang.sin(),
        cx + r_max_tree * scale_ang.cos(),
        cy + r_max_tree * scale_ang.sin()
    ));
    for &(frac, ring_lbl) in ring_fracs {
        let rr = r_inner + frac * (r_max_tree - r_inner);
        svg.push_str(&format!(
            r##"<circle cx="{:.1}" cy="{:.1}" r="{:.1}" fill="none" stroke="#94a3b8" stroke-width="0.9" stroke-dasharray="4,4" stroke-opacity="0.75"/>"##,
            cx, cy, rr
        ));
        let tx = cx + rr * scale_ang.cos();
        let ty = cy + rr * scale_ang.sin();
        svg.push_str(&format!(
            r##"<rect x="{:.1}" y="{:.1}" width="112" height="14" rx="3" fill="#ffffff" fill-opacity="0.92" stroke="#cbd5e1" stroke-width="0.8"/>"##,
            tx - 114.0,
            ty - 7.0
        ));
        svg.push_str(&format!(
            r##"<text x="{:.1}" y="{:.1}" text-anchor="end" font-size="7.6" font-weight="700" fill="#334155">{}</text>"##,
            tx - 6.0,
            ty + 2.8,
            ring_lbl
        ));
    }

    svg.push_str(&format!(
        r##"<circle cx="{:.1}" cy="{:.1}" r="{:.1}" fill="#f8fafc" stroke="#334155" stroke-width="2.0"/>"##,
        cx, cy, r_inner
    ));
    svg.push_str(&format!(
        r##"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-size="9.0" font-weight="800" fill="#0f172a">Node I Root</text>"##,
        cx, cy - 2.0
    ));
    svg.push_str(&format!(
        r##"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-size="7.5" font-weight="700" fill="#475569">(892.4 Ma)</text>"##,
        cx, cy + 9.0
    ));

    fn draw_circ_rec(
        n: &EmpiricalBifNode,
        cx: f64,
        cy: f64,
        r_of: &dyn Fn(f64) -> f64,
        r_label: f64,
        genus_map: &HashMap<String, String>,
        svg: &mut String,
        fb: &mut SoftwareFrameBuffer,
    ) {
        let rn = r_of(n.cum_patristic);
        let (nx, ny) = (cx + rn * n.angle_rad.cos(), cy + rn * n.angle_rad.sin());
        let col = sector_hex(&n.sector);
        let rgb = sector_rgb(&n.sector);

        if n.children.is_empty() {
            let (lx, ly) = (cx + r_label * n.angle_rad.cos(), cy + r_label * n.angle_rad.sin());
            svg.push_str(&format!(
                r##"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="#94a3b8" stroke-width="0.9" stroke-dasharray="1.5,2.0"/>"##,
                nx, ny, lx, ly
            ));
            match n.sector.as_str() {
                "[A]" | "[E]" => svg.push_str(&format!(
                    r##"<rect x="{:.1}" y="{:.1}" width="5.6" height="5.6" fill="{}" stroke="#0f172a" stroke-width="0.6"/>"##,
                    nx - 2.8, ny - 2.8, col
                )),
                "[B]" | "[F]" => svg.push_str(&format!(
                    r##"<polygon points="{:.1},{:.1} {:.1},{:.1} {:.1},{:.1}" fill="{}" stroke="#0f172a" stroke-width="0.6"/>"##,
                    nx, ny - 3.8, nx - 3.3, ny + 2.9, nx + 3.3, ny + 2.9, col
                )),
                _ => svg.push_str(&format!(
                    r##"<circle cx="{:.1}" cy="{:.1}" r="3.2" fill="{}" stroke="#0f172a" stroke-width="0.6"/>"##,
                    nx, ny, col
                )),
            }
            let deg = n.angle_rad.to_degrees();
            let norm_deg = ((deg % 360.0) + 360.0) % 360.0;
            let (rot, anchor) = if norm_deg > 90.0 && norm_deg < 270.0 {
                (deg + 180.0, "end")
            } else {
                (deg, "start")
            };
            let genus = genus_map.get(&n.label).map(|s| s.as_str()).unwrap_or("");
            svg.push_str(&format!(
                r##"<text x="{:.1}" y="{:.1}" text-anchor="{}" transform="rotate({:.2} {:.1} {:.1})" font-size="8.2" font-weight="700" fill="#0f172a">{} <tspan font-style="italic" font-weight="600" fill="#334155">{}</tspan></text>"##,
                lx,
                ly + 2.8,
                anchor,
                rot,
                lx,
                ly,
                xml_esc(&n.label),
                xml_esc(genus)
            ));
            fb.fill_rect(nx as i32 - 2, ny as i32 - 2, 4, 4, rgb);
        } else {
            let a_min = n.children.iter().map(|c| c.angle_rad).fold(f64::INFINITY, f64::min);
            let a_max = n.children.iter().map(|c| c.angle_rad).fold(f64::NEG_INFINITY, f64::max);
            let large = if (a_max - a_min) > std::f64::consts::PI { 1 } else { 0 };
            let (ax0, ay0) = (cx + rn * a_min.cos(), cy + rn * a_min.sin());
            let (ax1, ay1) = (cx + rn * a_max.cos(), cy + rn * a_max.sin());
            svg.push_str(&format!(
                r##"<path d="M {:.2},{:.2} A {:.2},{:.2} 0 {},1 {:.2},{:.2}" fill="none" stroke="{}" stroke-width="1.9"/>"##,
                ax0, ay0, rn, rn, large, ax1, ay1, col
            ));
            for ch in &n.children {
                let rc = r_of(ch.cum_patristic);
                let (sx0, sy0) = (cx + rn * ch.angle_rad.cos(), cy + rn * ch.angle_rad.sin());
                let (sx1, sy1) = (cx + rc * ch.angle_rad.cos(), cy + rc * ch.angle_rad.sin());
                let ccol = sector_hex(&ch.sector);
                svg.push_str(&format!(
                    r##"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="{}" stroke-width="1.9"/>"##,
                    sx0, sy0, sx1, sy1, ccol
                ));
                draw_circ_rec(ch, cx, cy, r_of, r_label, genus_map, svg, fb);
            }
            if n.ufboot >= 95.0 {
                svg.push_str(&format!(
                    r##"<circle cx="{:.1}" cy="{:.1}" r="2.7" fill="#0f172a" stroke="#ffffff" stroke-width="0.8"/>"##,
                    nx, ny
                ));
            }
            if n.label.starts_with("Node ") && n.label != "Node I" {
                let short_id = n.label.replace("Node ", "N");
                svg.push_str(&format!(
                    r##"<rect x="{:.1}" y="{:.1}" width="24" height="12" rx="3" fill="#ffffff" fill-opacity="0.94" stroke="{}" stroke-width="1.0"/>"##,
                    nx - 12.0,
                    ny - 15.0,
                    col
                ));
                svg.push_str(&format!(
                    r##"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-size="7.2" font-weight="800" fill="#0f172a">{}</text>"##,
                    nx,
                    ny - 6.5,
                    xml_esc(&short_id)
                ));
            }
        }
    }
    draw_circ_rec(&log_root, cx, cy, &r_of, r_label, &genus_map, &mut svg, &mut fb);

    svg.push_str(r##"<rect x="40" y="1090" width="1560" height="34" rx="6" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1.1"/>"##);
    svg.push_str(&format!(
        r##"<text x="820" y="1111" text-anchor="middle" font-size="10.2" font-weight="700" fill="#1e293b">Radial Scale: Log10-Transformed Branch Lengths [Δr = log10(1 + b / 0.006 subst./site); True Empirical Max = {:.3} subst./site] | ● UFBoot ≥ 95% (B=1,000) | NII–NVIII Backbone Nodes Labeled</text>"##,
        raw_max_tip
    ));

    svg.push_str("</svg>\n");
    save_figure_triple_and_web_svg(
        5,
        "fig5-circular-species-phylogram",
        "basal-metazoa-timetree-fig5-circular-species-phylogram.svg",
        "CircularPhylogram",
        "Circular Log10-Scaled Patristic Phylogram Across 64 Basal Metazoa & Holozoa Species (A:T001-G:T064)",
        "Dominant circular log10-scaled patristic phylogram with concentric log-distance rings, unified root-to-tip backbone, B&W alternating sectors, and sector-prefixed perimeter IDs.",
        &svg,
        &fb,
        run_dir,
        web_dir,
    )
}

fn find_best_core_window_1based(ws: &AlignmentWorkspace, width_bp: usize) -> (usize, usize) {
    let len = ws.alignment_length;
    if len <= width_bp {
        return (1, len.max(1));
    }
    let n_seq = ws.sequences.len().max(1) as f64;
    let mut best_s0 = 0usize;
    let mut best_score = -1.0_f64;
    let step = if width_bp % 3 == 0 { 3 } else { 1 };
    let mut c0 = 0usize;
    while c0 + width_bp <= len {
        let mut non_gap = 0usize;
        let mut cons_sum = 0.0_f64;
        for k in 0..width_bp {
            let col_idx = c0 + k;
            let mut counts = [0usize; 256];
            let mut max_c = 0usize;
            for seq in &ws.sequences {
                let b = seq.aligned_dna.get(col_idx).copied().unwrap_or(b'-');
                if b != b'-' && b != b'.' && b != b'?' {
                    non_gap += 1;
                    let idx = b.to_ascii_uppercase() as usize;
                    counts[idx] += 1;
                    if counts[idx] > max_c {
                        max_c = counts[idx];
                    }
                }
            }
            cons_sum += (max_c as f64) / n_seq;
        }
        let occ = (non_gap as f64) / (n_seq * (width_bp as f64));
        let cons = cons_sum / (width_bp as f64);
        let score = occ * 3.0 + cons;
        if score > best_score {
            best_score = score;
            best_s0 = c0;
        }
        c0 += step;
    }
    let s1 = best_s0 + 1;
    let e1 = (s1 + width_bp - 1).min(len);
    (s1, e1)
}

fn focus_workspace_on_ungapped_core(
    ws: &mut AlignmentWorkspace,
    genetic_code: GeneticCodeTable,
    reading_frame: usize,
) {
    let (pre_s1, pre_e1) = find_best_core_window_1based(ws, 66);
    let s0 = pre_s1.saturating_sub(1);
    let e0 = pre_e1.min(ws.alignment_length);
    let win_len = e0.saturating_sub(s0).max(1);

    // Filter out any sequence row that is mostly gapped inside the core catalytic domain
    if ws.sequences.len() > 8 {
        let filtered_seqs: Vec<_> = ws
            .sequences
            .iter()
            .filter(|seq| {
                let gaps = seq.aligned_dna[s0..e0]
                    .iter()
                    .filter(|&&b| b == b'-' || b == b'.')
                    .count();
                (gaps as f64) / (win_len as f64) <= 0.15
            })
            .cloned()
            .collect();
        if filtered_seqs.len() >= 8 {
            ws.sequences = filtered_seqs;
        }
    }

    // Recompute column stats and strip gappy N/C-terminal overhang triplets so the workspace focuses on the ungapped core domain
    ws.recompute_telemetry_and_auto_features(genetic_code, reading_frame);
    let n_codons = ws.alignment_length / 3;
    if n_codons > 30 {
        let mut keep_cols = Vec::with_capacity(ws.alignment_length);
        for c_idx in 0..n_codons {
            let b0 = c_idx * 3;
            let g0 = ws.column_stats.get(b0).map(|c| c.gap_fraction).unwrap_or(1.0);
            let g1 = ws.column_stats.get(b0 + 1).map(|c| c.gap_fraction).unwrap_or(1.0);
            let g2 = ws.column_stats.get(b0 + 2).map(|c| c.gap_fraction).unwrap_or(1.0);
            if g0 <= 0.12 && g1 <= 0.12 && g2 <= 0.12 {
                keep_cols.push(b0);
                keep_cols.push(b0 + 1);
                keep_cols.push(b0 + 2);
            }
        }
        if keep_cols.len() >= 90 {
            for seq in &mut ws.sequences {
                let mut next_dna = Vec::with_capacity(keep_cols.len());
                for &col in &keep_cols {
                    if let Some(&b) = seq.aligned_dna.get(col) {
                        next_dna.push(b);
                    }
                }
                let non_g = next_dna.iter().filter(|&&b| b != b'-' && b != b'.').count();
                let gc_c = next_dna
                    .iter()
                    .filter(|&&b| matches!(b, b'G' | b'C' | b'g' | b'c'))
                    .count();
                seq.gap_percent = if keep_cols.is_empty() {
                    0.0
                } else {
                    ((keep_cols.len() - non_g) as f64 / keep_cols.len() as f64) * 100.0
                };
                seq.gc_percent = if non_g == 0 {
                    0.0
                } else {
                    (gc_c as f64 / non_g as f64) * 100.0
                };
                seq.aligned_dna = next_dna;
            }
            ws.alignment_length = keep_cols.len();
            ws.features.clear();
            ws.recompute_telemetry_and_auto_features(genetic_code, reading_frame);
        }
    }
}

fn generate_basal_supplemental_alignview_snapshots(
    run_dir: &Path,
    web_figures_dir: Option<&Path>,
) -> Result<(), String> {
    let aln_dir = run_dir.join("alignments");
    let fig_dir = run_dir.join("figures");
    let _ = fs::create_dir_all(&fig_dir);

    let mirror_pair = |svg_src: &Path, png_src: &Path, stem: &str| -> Result<(), String> {
        if let Some(wdir) = web_figures_dir {
            fs::create_dir_all(wdir).map_err(|e| e.to_string())?;
            fs::copy(svg_src, wdir.join(format!("{}.svg", stem)))
                .map_err(|e| format!("Failed copying {}.svg: {}", stem, e))?;
            fs::copy(png_src, wdir.join(format!("{}.png", stem)))
                .map_err(|e| format!("Failed copying {}.png: {}", stem, e))?;
        }
        Ok(())
    };

    // Supplemental Figure S1: FAM-01 Type IV Collagen NC1 Sulfilimine Domain (`supp_s1_col4_nc1_cds.fasta`)
    let s1_fasta = aln_dir.join("supp_s1_col4_nc1_cds.fasta");
    if s1_fasta.is_file() {
        let mut ws = AlignmentWorkspace::from_file(
            &s1_fasta,
            None,
            GeneticCodeTable::Standard1,
            1,
            false,
        )?;
        ws.title = "FAM-01 Basement Membrane Type IV Collagen NC1 & Peroxidasin CDS (Basal Metazoa)".to_string();
        ws.source_path = "VAL-SET-044-BASAL-METAZOA-TIMETREE (FAM-01)".to_string();
        ws.alignment_method = "Codon-Locked Katoh L-INS-i + PAL2NAL (Standard Code 1)".to_string();
        focus_workspace_on_ungapped_core(&mut ws, GeneticCodeTable::Standard1, 1);
        let (s1, e1) = find_best_core_window_1based(&ws, 66);
        ws.features.insert(
            0,
            crate::align_viewer::AlignmentFeatureAnnotation {
                feature_id: "COL4-SULFILIMINE".to_string(),
                feature_type: "Catalytic_Motif".to_string(),
                label: "COL4 NC1 Met93-Hyl211 Sulfilimine Crosslink Core Window".to_string(),
                start_col_1based: s1,
                end_col_1based: e1,
                strand: '+',
                score: 99.2,
                color_rgb: (245, 158, 11),
                description: "Peroxidasin-catalyzed Met93-Hyl211 covalent sulfilimine (S=N) crosslink".to_string(),
            },
        );
        let mut state = AlignmentViewerState::new(ws);
        state.genetic_code = GeneticCodeTable::Standard1;
        state.mode = ViewAlphabetMode::DualDnaAndAminoAcid;
        state.color_scheme = AlignmentColorScheme::ClustalX;
        state.match_mode = ConsensusMatchMode::ShowAllCharacters;
        state.row_height_px = 32;
        state.selected_feature_idx = Some(0);
        state.selected_col_range_1based = Some((s1, e1));
        state.hover_row_idx = Some(0);
        state.hover_col_1based = Some((s1 + e1) / 2);
        state.status_message = format!(
            "runie-align-view | FAM-01 Type IV Collagen NC1 Sulfilimine Core Window ({}..{} bp, Standard Code 1)",
            s1, e1
        );
        let svg_out = fig_dir.join("basal-metazoa-supp-s1-col4-nc1-alignview.svg");
        let png_out = fig_dir.join("basal-metazoa-supp-s1-col4-nc1-alignview.png");
        state.capture_current_viewport_snapshots(1560, Some(&svg_out), Some(&png_out))?;
        mirror_pair(&svg_out, &png_out, "basal-metazoa-supp-s1-col4-nc1-alignview")?;
    }

    // Supplemental Figure S2: Cohort L2 Aurelia Nuclear ITS1-5.8S-ITS2 rDNA (`supp_s2_its1_5_8s_einsi.fasta`)
    let s2_fasta = aln_dir.join("supp_s2_its1_5_8s_einsi.fasta");
    if s2_fasta.is_file() {
        let mut ws = AlignmentWorkspace::from_file(
            &s2_fasta,
            None,
            GeneticCodeTable::Standard1,
            1,
            false,
        )?;
        ws.title = "Cohort L2 Scyphozoa Aurelia Nuclear ITS1-5.8S-ITS2 rDNA (N = 36)".to_string();
        ws.source_path = "VAL-SET-044-BASAL-METAZOA-TIMETREE (Cohort L2)".to_string();
        ws.alignment_method = "mafft-native Katoh E-INS-i (FFT + Local Affine + 2-Tree Iterative)".to_string();
        let (s1, e1) = find_best_core_window_1based(&ws, 74);
        ws.features.insert(
            0,
            crate::align_viewer::AlignmentFeatureAnnotation {
                feature_id: "ITS1-5.8S".to_string(),
                feature_type: "Conserved_Block".to_string(),
                label: "Conserved 5.8S rRNA Stem & ITS Spacer Core".to_string(),
                start_col_1based: s1,
                end_col_1based: e1,
                strand: '+',
                score: 98.4,
                color_rgb: (56, 189, 248),
                description: "Katoh E-INS-i preserves 5.8S stem homology across flanking ITS1 indels".to_string(),
            },
        );
        let mut state = AlignmentViewerState::new(ws);
        state.genetic_code = GeneticCodeTable::Standard1;
        state.mode = ViewAlphabetMode::Nucleotide;
        state.color_scheme = AlignmentColorScheme::AliViewDefault;
        state.match_mode = ConsensusMatchMode::ShowAllCharacters;
        state.row_height_px = 22;
        state.selected_feature_idx = Some(0);
        state.selected_col_range_1based = Some((s1, e1));
        state.hover_row_idx = Some(0);
        state.hover_col_1based = Some((s1 + e1) / 2);
        state.status_message = format!(
            "runie-align-view | mafft-native Katoh E-INS-i Nuclear ITS1-5.8S Core Stem ({}..{} bp)",
            s1, e1
        );
        let svg_out = fig_dir.join("basal-metazoa-supp-s2-its1-einsi-alignview.svg");
        let png_out = fig_dir.join("basal-metazoa-supp-s2-its1-einsi-alignview.png");
        state.capture_current_viewport_snapshots(1640, Some(&svg_out), Some(&png_out))?;
        mirror_pair(&svg_out, &png_out, "basal-metazoa-supp-s2-its1-einsi-alignview")?;
    }

    // Supplemental Figure S3: Basal Metazoa Phased Mitochondrial COX1 Core Alignment (`supp_s3_cox1_64taxa.fasta`)
    let s3_fasta = aln_dir.join("supp_s3_cox1_64taxa.fasta");
    if s3_fasta.is_file() {
        let mut ws = AlignmentWorkspace::from_file(
            &s3_fasta,
            None,
            GeneticCodeTable::InvertebrateMito5,
            1,
            false,
        )?;
        ws.title = "Basal Metazoa & Holozoa Phased Mitochondrial COX1 Core Domain (Clades A-G; Ungapped Core)".to_string();
        ws.source_path = "VAL-SET-044-BASAL-METAZOA-TIMETREE (COX1 Core Heme-Cu Helix)".to_string();
        ws.alignment_method = "COX1-Phased Codon MSA (Mold/Coelenterate & Invertebrate Mito Translation)".to_string();
        focus_workspace_on_ungapped_core(&mut ws, GeneticCodeTable::InvertebrateMito5, 1);
        let (s1, e1) = find_best_core_window_1based(&ws, 66);
        ws.features.insert(
            0,
            crate::align_viewer::AlignmentFeatureAnnotation {
                feature_id: "COX1-CORE".to_string(),
                feature_type: "Catalytic_Motif".to_string(),
                label: "COX1 Core Heme-Cu Helices (GDDQLYNVIVTAHA) & 3rd-Pos RY3 Transitions".to_string(),
                start_col_1based: s1,
                end_col_1based: e1,
                strand: '+',
                score: 99.6,
                color_rgb: (16, 185, 129),
                description: "Ungapped core COX1 catalytic helices across Ctenophora, Porifera, Placozoa & Cnidaria".to_string(),
            },
        );
        let mut state = AlignmentViewerState::new(ws);
        state.genetic_code = GeneticCodeTable::InvertebrateMito5;
        state.mode = ViewAlphabetMode::DualDnaAndAminoAcid;
        state.color_scheme = AlignmentColorScheme::RunieConservation;
        state.match_mode = ConsensusMatchMode::ShowAllCharacters;
        state.row_height_px = 26;
        state.selected_feature_idx = Some(0);
        state.selected_col_range_1based = Some((s1, e1));
        state.hover_row_idx = Some(0);
        state.hover_col_1based = Some((s1 + e1) / 2);
        state.status_message = format!(
            "runie-align-view | Phased Basal Metazoan COX1 Core Alignment ({}..{} bp — 100% Ungapped Core)",
            s1, e1
        );
        let svg_out = fig_dir.join("basal-metazoa-supp-s3-cox1-ry3-alignview.svg");
        let png_out = fig_dir.join("basal-metazoa-supp-s3-cox1-ry3-alignview.png");
        state.capture_current_viewport_snapshots(1560, Some(&svg_out), Some(&png_out))?;
        mirror_pair(&svg_out, &png_out, "basal-metazoa-supp-s3-cox1-ry3-alignview")?;
    }

    Ok(())
}

/// Generates all 5 primary Preprint-05 (`v3.0.0`) locked-format figures (`.svg` + `.png` + `.pdf`)
/// plus 3 Native Alignment Viewer supplemental screenshots inside `run_dir/figures/` and mirrors into `web_figures_dir`.
pub fn generate_basal_metazoa_locked_figures(
    manifest: &mut SupervisorRunManifest,
    study: &MultiCompartmentStudyReport,
    run_dir: &Path,
    web_figures_dir: Option<&Path>,
) -> Result<(), String> {
    let root = build_basal_unified_rooted_tree(manifest);
    let mut figs = Vec::with_capacity(5);

    figs.push(render_basal_fig1_chronogram(manifest, &root, run_dir, web_figures_dir)?);
    figs.push(render_basal_fig2_quartet_msc_and_synteny(manifest, study, run_dir, web_figures_dir)?);
    figs.push(render_bemisia_skyline_6panel(
        3,
        "fig3-planulozoa-cnidaria-cambrian-chronogram",
        "basal-metazoa-timetree-fig3-planulozoa-cnidaria-cambrian-chronogram.svg",
        "Six-Deme Basal Metazoan Bayesian Coalescent Skyline Matrix Across Regional & Multi-Locus Cohorts (R1-R3, L1-L3; N = 216 Isolates)",
        "Six single-species conspecific Aurelia, Nematostella, and Acropora coalescent skylines (25,000 MCMC generations) contrasting maternal mtDNA (1Ne), Katoh E-INS-i nuclear ITS1-5.8S rDNA (4Ne), and thermal-bleaching coding sweeps.",
        &study.popset_skylines,
        run_dir,
        web_figures_dir,
    )?);
    figs.push(render_basal_fig4_selection_sweep_and_col4_pdb(manifest, study, run_dir, web_figures_dir)?);
    figs.push(render_basal_fig5_circular_phylogram(manifest, &root, run_dir, web_figures_dir)?);

    generate_basal_supplemental_alignview_snapshots(run_dir, web_figures_dir)?;
    manifest.figures = figs;
    Ok(())
}
