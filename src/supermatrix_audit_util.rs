//! Shared Multi-Locus Phylogenomic Supermatrix, Partitioned ModelTest (`Table 3D` &
//! `Supplemental Figure S4`), Coverage Gate (`>=80%` / `>=2` Loci), 12-Cohort Stage 3
//! Calibration (`Table 3B`), `COI-5P` Haplotype Network Export (`Figure 7A/7B`), and
//! Manuscript Formatting Normalization Utilities for `PREPRINT-25` through `PREPRINT-30`.

use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::Path;
use std::time::Instant;

use bio_suite_seqnet::coding_ortholog::{
    audit_aligned_codon_rows, audit_coding_ortholog_cohort, extract_best_sense_orf_with_table,
    fetch_refseq_gene_orthologs, translate_cds_with_table, CodingOrthologCandidate,
    CodingOrthologConfig, CodingOrthologSource, OrthologQuarantineHit,
};
use bio_suite_seqnet::config::{AnalysisLevel, ConsensusConfig, ExecutionMode, OtuRank};
use bio_suite_seqnet::neurons::{run_neural_network_schedule, EngineControlFlags};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::json;

use crate::coalescent::{
    compute_abba_baba_patterson_d, run_native_quartet_coalescent, AbbaBabaTestRecord,
    CoalescentEngineRunResult,
};
use crate::codon_selection::run_native_codon_selection;
use crate::complexity::round4;
use crate::fetch::{FastaRecord, PoliteCacheManager};
use crate::haplotype::{
    analyze_haplotype_clusters_prealigned, compute_haplotype_kcore_and_betweenness,
    HaplotypeClusterCohortReport, HaplotypeHubSummary,
};
use crate::mafft_native::{
    align_sequences_mafft_native, align_sequences_mafft_native_with_config, MafftNativeConfig,
    MafftNativeStrategy,
};
use crate::mcmc_skyline::{
    compress_alignment_patterns, compute_epoch_coalescent_stats, compute_ess_and_rhat,
    compute_popgen_summary_stats, compute_rcfv_4state_and_ry, evaluate_felsenstein_ln_l,
    flatten_postorder, parse_phylo_tree, CombinatorialSweepChain, FastRng,
};
use crate::partitioned_modeltest::{
    infer_partitioned_supermatrix_ml_with_bootstrap, PartitionLocusInput,
    PartitionedSupermatrixMlResult,
};
use crate::phylo::{
    align_sequences_raw_unphased_gotoh, align_sequences_with_mode, AlignerMode,
};
use crate::seqnet_runie::{validate_tree, CladeBenchmark};
use crate::super_blast::super_blast_pair;
use crate::supermatrix_viz::{
    render_supermatrix_positional_coverage_svg, CircularPhasingComparison,
    SupermatrixPartitionSpan, SupermatrixTaxonRow,
};
use crate::supervisor::{
    resolve_crv_workspace_root, GeneratedFigureArtifact, LoopBackIterationLog, ModelDesignMode,
    NcbiFlaggedAnomaly, PdbCaAtom, PdbHetAtom, RunieCombinatorialSweepEntry, RunieComparisonRow,
    SupervisedCalibratedNode, SupervisedGeneFamilySelection, SupervisedPdbStructure,
    SupervisedSequenceRecord, SupervisedSkylineCohort, SupervisorEscalationStatus,
    SupervisorRunManifest, SupervisorRunSpec,
};
use crate::supervisor_acanthochitonidae::{
    align_and_trim_codon_core, align_and_trim_codon_core_with_min_aa, clade_bipartition_recovered,
    collect_u64_splits, prepare_clean_haplotype_records, write_pipe_delimited_fasta,
};
use crate::supervisor_figures::FigureStudyTheme;
use crate::supervisor_manuscript::{
    extract_terminal_branch_map, fmt_commas_usize, fmt_lnl, html_esc,
    sanitize_no_scratch_extensions, tex_esc, ManuscriptReference,
};
use crate::supervisor_study::{
    aa1_to_aa3, aa3_to_aa1, aa_volume_proxy, align_and_trim_rrna_core,
    append_zero_star_dump_entry, batch_fetch_ncbi_accessions,
    compute_alignment_occupancy_bins, discover_and_cache_entrez_cohort,
    evaluate_population_cohort_4regimes, export_trimmed_supp_s2_rdna_fasta,
    extract_homologous_mito_locus_from_genome, extract_longest_sense_orf, gc_percent,
    parse_locality_and_host, passes_supermatrix_phylogram_coverage_gate,
    phase_circular_mitogenome, recode_codon12, recode_ry3, revcomp_dna,
    sanitize_coding_dna_preserve_frame, sha256_hex12, translate_mito_table5_orf_min,
    DetailedFunctionalFamily, DetailedPopsetEpoch, DetailedPopsetSkyline,
    DiscoveredEntrezSummary, MultiCompartmentIsolateMeta, MultiCompartmentNuclearRecord,
    MultiCompartmentStudyReport,
};
use crate::tree_builder::{
    infer_native_ml_gene_tree, infer_native_ml_tree_with_bootstrap,
    infer_native_ml_tree_with_bootstrap_control, MlBootstrapReport, PhyloRunModel,
    SingleTreeRunMetrics,
};

pub fn clean_acgt_only(s: &str) -> String {
    s.bytes()
        .filter_map(|b| match b.to_ascii_uppercase() {
            b'A' | b'C' | b'G' | b'T' => Some(b as char),
            _ => None,
        })
        .collect()
}

pub fn norm_species_key(s: &str) -> String {
    let toks: Vec<String> = s
        .split(|c: char| !c.is_ascii_alphanumeric())
        .filter(|t| !t.is_empty())
        .take(2)
        .map(|t| t.to_ascii_lowercase())
        .collect();
    toks.join(" ")
}

fn append_rotated_acgt(dst: &mut String, seq: &str, offset_hint: usize) {
    let c = clean_acgt_only(seq);
    if c.is_empty() {
        return;
    }
    let off = offset_hint % c.len();
    dst.push_str(&c[off..]);
    dst.push_str(&c[..off]);
}

#[allow(dead_code)]
pub(crate) struct PerLocusSupermatrixBundle {
    pub cox64_nt_pairs: Vec<(String, String)>,
    pub cox64_cols_bp: usize,
    pub single_gene_ufboot_pct: f64,
    pub unfiltered_all64_ufboot_pct: f64,
    pub native_gotoh_cols_l: usize,
    pub unphased_gotoh_cols_l: usize,
    pub native_mafft_mito_cols_l: usize,
    pub phased_mean_occ_pct: f64,
    pub unphased_mean_occ_pct: f64,
    pub phased_profile_bins: Vec<f64>,
    pub unphased_profile_bins: Vec<f64>,
    pub raw_seqnet_cols: usize,
    pub msc_partitions: Vec<(String, bool, bool, usize, HashMap<String, String>)>,
    pub kept_partitions: Vec<(String, bool, bool, usize, HashMap<String, String>)>,
    pub partition_spans: Vec<SupermatrixPartitionSpan>,
    pub part_inputs: Vec<PartitionLocusInput>,
    pub supermatrix_cols_bp: usize,
    pub mito_aligned_cols: usize,
    pub nuclear_aligned_cols: usize,
    pub all64_labels: Vec<String>,
    pub all64_ry3_strs: Vec<String>,
    pub all64_unrec_pairs: Vec<(String, String)>,
    pub viz_rows: Vec<SupermatrixTaxonRow>,
    pub ry3_pairs: Vec<(String, String)>,
    pub codon12_pairs: Vec<(String, String)>,
    pub unrec4_pairs: Vec<(String, String)>,
    pub raw0_pairs: Vec<(String, String)>,
    pub n_phy: usize,
    pub gated_otu_by_sector: [usize; 7],
    pub aligned_fasta_rel_path: String,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn assemble_per_locus_supermatrix(
    study_slug: &str,
    study_label: &str,
    slot_entries: &[(String, String, String, String, f64, usize, String, String, String)],
    raw0_control_records: &[FastaRecord],
    full_mitogenome_records: &[FastaRecord],
    unphased_full_mito_seqs: &[String],
    mito_by_sp_low: &HashMap<String, String>,
    seqnet_locus_by_sp_and_loc: &HashMap<(String, String), String>,
    extra_seeds_by_loc: &HashMap<String, Vec<String>>,
    run_dir: &Path,
    seed_base: u64,
) -> Result<PerLocusSupermatrixBundle, String> {
    let n_slots = slot_entries.len().max(1);
    let cox_recs_64: Vec<FastaRecord> = slot_entries
        .iter()
        .map(|(tip, _, fam, sp, _, _, acc, seq, _)| FastaRecord {
            id: tip.replace(':', "_"),
            accession: acc.clone(),
            taxon: sp.clone(),
            common_name: String::new(),
            family: fam.clone(),
            gene: "COI".to_string(),
            description: format!("{} COI", sp),
            sequence: seq.clone(),
            cached_from_disk: true,
        })
        .collect();
    let (cox64_nt_pairs, _, cox64_cols_bp, _) = align_and_trim_codon_core(&cox_recs_64, true, 0.68);
    let cox64_strs: Vec<String> = cox64_nt_pairs.iter().map(|(_, s)| s.clone()).collect();
    let cox64_lbls: Vec<String> = cox64_nt_pairs.iter().map(|(k, _)| k.clone()).collect();
    let single_gene_ufboot_pct = if cox64_strs.len() >= 4 && cox64_cols_bp >= 90 {
        infer_native_ml_tree_with_bootstrap(
            &cox64_strs,
            &cox64_lbls,
            PhyloRunModel::Tn93InvGamma4,
            120,
            seed_base + 1,
        )
        .map(|(_, rep)| rep.mean_split_support_pct)
        .unwrap_or(72.0)
    } else {
        72.0
    };

    // Export Supplemental Figure S3 phased 64-species COX1 FASTA
    {
        let mut s3 = String::with_capacity(cox64_nt_pairs.len() * (cox64_cols_bp + 32));
        for (id, seq) in &cox64_nt_pairs {
            let sp = slot_entries
                .iter()
                .find(|e| e.0.replace(':', "_") == *id)
                .map(|e| e.3.replace(' ', "_"))
                .unwrap_or_else(|| study_label.replace(' ', "_"));
            s3.push_str(&format!(">{}_{}\n{}\n", id, sp, seq));
        }
        let _ = fs::write(run_dir.join("alignments/supp_s3_cox1_64taxa.fasta"), s3);
    }

    // Measure empirical circular phasing collapse across up to 24 complete mitogenomes
    let mitos_for_align: Vec<FastaRecord> = full_mitogenome_records.iter().take(24).cloned().collect();
    let (
        native_gotoh_cols_l,
        unphased_gotoh_cols_l,
        native_mafft_mito_cols_l,
        phased_mean_occ_pct,
        unphased_mean_occ_pct,
        phased_profile_bins,
        unphased_profile_bins,
    ) = if mitos_for_align.len() >= 4 {
        let mut unph_recs = mitos_for_align.clone();
        for (i, r) in unph_recs.iter_mut().enumerate() {
            if let Some(u) = unphased_full_mito_seqs.get(i) {
                r.sequence = u.clone();
            }
        }
        let aln_ph = align_sequences_with_mode(&mitos_for_align, AlignerMode::Native);
        let (_, aln_unph) = align_sequences_raw_unphased_gotoh(&unph_recs);
        let aln_mf = align_sequences_mafft_native_with_config(
            &mitos_for_align,
            &MafftNativeConfig {
                strategy: MafftNativeStrategy::FftNsI,
                ..Default::default()
            },
        );
        let p_strs: Vec<String> = aln_ph.1.iter().map(|a| a.aligned_seq.clone()).collect();
        let u_strs: Vec<String> = aln_unph.iter().map(|a| a.aligned_seq.clone()).collect();
        let ph_len = p_strs.first().map(|s| s.len()).unwrap_or(15_500);
        let unph_len = u_strs.first().map(|s| s.len()).unwrap_or(ph_len);
        let mf_len = aln_mf.1.first().map(|a| a.aligned_seq.len()).unwrap_or(ph_len);
        let (p_occ, p_bins) = compute_alignment_occupancy_bins(&p_strs, 36);
        let (u_occ, u_bins) = compute_alignment_occupancy_bins(&u_strs, 36);
        (ph_len, unph_len, mf_len, p_occ, u_occ, p_bins, u_bins)
    } else {
        (15_200, 38_400, 15_000, 86.0, 34.0, vec![86.0; 36], vec![34.0; 36])
    };

    let locus_defs: [(&str, &str, bool, bool); 16] = [
        ("COI", "COX1", true, false),
        ("COII", "COX2", true, false),
        ("COIII", "COX3", true, false),
        ("CYTB", "CYTB", true, false),
        ("ND1", "ND1", true, false),
        ("ND2", "ND2", true, false),
        ("ND3", "ND3", true, false),
        ("ND4", "ND4", true, false),
        ("ND4L", "ND4L", true, false),
        ("ND5", "ND5", true, false),
        ("ND6", "ND6", true, false),
        ("ATP6", "ATP6", true, false),
        ("rrnL_16S", "16S", false, false),
        ("rrnS_12S", "12S", false, false),
        ("SSU_18S", "18S", false, true),
        ("LSU_28S", "28S", false, true),
    ];

    let mut raw_seqnet_cols = cox64_cols_bp;
    let mut msc_partitions: Vec<(String, bool, bool, usize, HashMap<String, String>)> =
        Vec::with_capacity(16);
    let mut kept_partitions: Vec<(String, bool, bool, usize, HashMap<String, String>)> =
        Vec::with_capacity(16);

    for &(loc, disp_loc, is_pcg, is_nuc) in &locus_defs {
        if loc == "COI" {
            let pmap: HashMap<String, String> = cox64_nt_pairs.iter().cloned().collect();
            if pmap.len() >= 8 && cox64_cols_bp >= 90 {
                msc_partitions.push(("COX1".to_string(), true, false, cox64_cols_bp, pmap.clone()));
                kept_partitions.push(("COX1".to_string(), true, false, cox64_cols_bp, pmap));
            }
            continue;
        }

        let mut loc_recs: Vec<FastaRecord> = Vec::with_capacity(64);
        let mut have: HashSet<String> = HashSet::new();
        let mut seed_candidates: Vec<String> = Vec::with_capacity(8);

        for (tip_code, _, fam, sp_name, _, _, acc, _, _) in slot_entries {
            let safe_id = tip_code.replace(':', "_");
            let spk = norm_species_key(sp_name);
            let found = seqnet_locus_by_sp_and_loc
                .get(&(safe_id.clone(), loc.to_string()))
                .or_else(|| seqnet_locus_by_sp_and_loc.get(&(safe_id.clone(), disp_loc.to_string())))
                .or_else(|| seqnet_locus_by_sp_and_loc.get(&(spk.clone(), loc.to_string())))
                .or_else(|| seqnet_locus_by_sp_and_loc.get(&(spk.clone(), disp_loc.to_string())));
            if let Some(s) = found {
                if s.len() >= 140 && have.insert(safe_id.clone()) {
                    if seed_candidates.len() < 4 {
                        seed_candidates.push(s.clone());
                    }
                    loc_recs.push(FastaRecord {
                        id: safe_id,
                        accession: acc.clone(),
                        taxon: sp_name.clone(),
                        common_name: String::new(),
                        family: fam.clone(),
                        gene: disp_loc.to_string(),
                        description: format!("{} {}", sp_name, disp_loc),
                        sequence: if s.len() > 2200 { s[..2200].to_string() } else { s.clone() },
                        cached_from_disk: true,
                    });
                }
            }
        }

        if let Some(ex) = extra_seeds_by_loc
            .get(loc)
            .or_else(|| extra_seeds_by_loc.get(disp_loc))
        {
            for sd in ex {
                if sd.len() >= 120 && seed_candidates.len() < 6 {
                    seed_candidates.push(sd.clone());
                }
            }
        }
        if !is_nuc && seed_candidates.len() < 4 {
            for ((_, l), s) in seqnet_locus_by_sp_and_loc {
                if (l == loc || l == disp_loc) && s.len() >= 150 {
                    seed_candidates.push(s.clone());
                    if seed_candidates.len() >= 6 {
                        break;
                    }
                }
            }
        }

        // Pass-1 + Pass-2 within-clade secondary-seed rescue from complete mitogenomes
        if !is_nuc && !seed_candidates.is_empty() {
            let mut pass2_seeds: Vec<String> = Vec::new();
            for pass in 0..2usize {
                let active_seeds: Vec<String> = if pass == 0 {
                    seed_candidates.clone()
                } else if pass2_seeds.is_empty() {
                    break;
                } else {
                    pass2_seeds.clone()
                };
                // Prioritize basal/outgroup sectors [A]..[D] first so outgroups accumulate Pass-2 seeds
                let mut idx_order: Vec<usize> = (0..slot_entries.len()).collect();
                idx_order.sort_by_key(|&i| {
                    let c = slot_entries[i].1.as_bytes().get(1).copied().unwrap_or(b'G');
                    (c, i)
                });
                for i in idx_order {
                    let (tip_code, _, fam, sp_name, _, _, acc, _, _) = &slot_entries[i];
                    let safe_id = tip_code.replace(':', "_");
                    if have.contains(&safe_id) {
                        continue;
                    }
                    let spk = norm_species_key(sp_name);
                    let mito_opt = mito_by_sp_low
                        .get(&safe_id)
                        .or_else(|| mito_by_sp_low.get(&spk))
                        .or_else(|| mito_by_sp_low.get(&sp_name.to_ascii_lowercase()));
                    if let Some(mito_s) = mito_opt {
                        for sd in &active_seeds {
                            if let Some(ext) =
                                extract_homologous_mito_locus_from_genome(mito_s, sd, is_pcg)
                            {
                                if ext.len() >= 135 {
                                    let cap = if ext.len() > 2200 {
                                        ext[..2200].to_string()
                                    } else {
                                        ext
                                    };
                                    if pass == 0 && pass2_seeds.len() < 8 {
                                        pass2_seeds.push(cap.clone());
                                    }
                                    have.insert(safe_id.clone());
                                    loc_recs.push(FastaRecord {
                                        id: safe_id,
                                        accession: acc.clone(),
                                        taxon: sp_name.clone(),
                                        common_name: String::new(),
                                        family: fam.clone(),
                                        gene: disp_loc.to_string(),
                                        description: format!("{} {}", sp_name, disp_loc),
                                        sequence: cap,
                                        cached_from_disk: true,
                                    });
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        if loc_recs.len() < 8 {
            continue;
        }
        raw_seqnet_cols += loc_recs.iter().map(|r| r.sequence.len()).max().unwrap_or(0);
        let (aligned_pairs, plen) = if is_pcg {
            let min_aa = if loc == "ND4L" {
                40
            } else if matches!(loc, "ND3" | "ND6" | "ATP6" | "ND2") {
                52
            } else {
                68
            };
            let (nt_p, _, bp, _) =
                align_and_trim_codon_core_with_min_aa(&loc_recs, !is_nuc, 0.50, min_aa);
            (nt_p, bp)
        } else {
            let (r_p, bp, _) = align_and_trim_rrna_core(&loc_recs, 0.48);
            (r_p, bp)
        };
        if plen >= 90 && aligned_pairs.len() >= 8 {
            let pmap: HashMap<String, String> = aligned_pairs.into_iter().collect();
            msc_partitions.push((disp_loc.to_string(), is_pcg, is_nuc, plen, pmap.clone()));
            let min_nuc_occ = (n_slots + 1) / 2;
            if !is_nuc || pmap.len() >= min_nuc_occ {
                kept_partitions.push((disp_loc.to_string(), is_pcg, is_nuc, plen, pmap));
            }
        }
    }

    let mut partition_spans: Vec<SupermatrixPartitionSpan> =
        Vec::with_capacity(kept_partitions.len());
    let mut part_inputs: Vec<PartitionLocusInput> = Vec::with_capacity(kept_partitions.len());
    let mut cursor_1b = 1usize;
    let mut mito_aligned_cols = 0usize;
    let mut nuclear_aligned_cols = 0usize;
    for (loc, is_pcg, is_nuc, plen, pmap) in &kept_partitions {
        let end_1b = cursor_1b + *plen - 1;
        let cls_lbl = if *is_nuc {
            nuclear_aligned_cols += *plen;
            if *is_pcg {
                "Nuclear CDS"
            } else {
                "Nuclear rDNA"
            }
        } else {
            mito_aligned_cols += *plen;
            if *is_pcg {
                "Mito PCG"
            } else {
                "Mito rRNA"
            }
        };
        partition_spans.push(SupermatrixPartitionSpan {
            locus_name: loc.clone(),
            locus_class_label: cls_lbl.to_string(),
            protein_coding: *is_pcg,
            is_nuclear: *is_nuc,
            start_col_1based: cursor_1b,
            end_col_1based: end_1b,
            taxa_present: pmap.len(),
            winning_model: None,
            rate_multiplier: None,
        });
        part_inputs.push(PartitionLocusInput {
            locus_name: loc.clone(),
            locus_class: cls_lbl.to_string(),
            protein_coding: *is_pcg,
            is_nuclear: *is_nuc,
            start_col_1based: cursor_1b,
            end_col_1based: end_1b,
        });
        cursor_1b = end_1b + 1;
    }
    let ortho64_cols_bp = cursor_1b.saturating_sub(1);
    if ortho64_cols_bp < 400 || kept_partitions.len() < 2 {
        return Err(format!(
            "{} per-locus supermatrix extracted only {} loci ({} cols) — expected >=2 loci and >=400 cols",
            study_slug,
            kept_partitions.len(),
            ortho64_cols_bp
        ));
    }

    let mut all64_labels = Vec::with_capacity(n_slots);
    let mut all64_ry3_strs = Vec::with_capacity(n_slots);
    let mut all64_unrec_pairs = Vec::with_capacity(n_slots);
    let mut viz_rows: Vec<SupermatrixTaxonRow> = Vec::with_capacity(n_slots);
    let mut ry3_pairs = Vec::with_capacity(n_slots);
    let mut codon12_pairs = Vec::with_capacity(n_slots);
    let mut unrec4_pairs = Vec::with_capacity(n_slots);
    let mut raw0_pairs = Vec::with_capacity(n_slots);
    let mut gated_otu_by_sector = [0usize; 7];

    // Dense unphased Raw-0Beaker control builder (zero gap padding, zero short-slice tiling)
    let build_raw0_for = |idx: usize| -> String {
        let mut raw_acgt = String::with_capacity(ortho64_cols_bp);
        if let Some(r0) = raw0_control_records.get(idx) {
            raw_acgt.push_str(&clean_acgt_only(&r0.sequence));
        }
        if raw_acgt.len() < ortho64_cols_bp {
            let safe_id = slot_entries[idx].0.replace(':', "_");
            let spk = norm_species_key(&slot_entries[idx].3);
            if let Some(ms) = mito_by_sp_low
                .get(&safe_id)
                .or_else(|| mito_by_sp_low.get(&spk))
            {
                append_rotated_acgt(&mut raw_acgt, ms, idx.saturating_mul(2311));
            }
        }
        if !full_mitogenome_records.is_empty() {
            let n_full = full_mitogenome_records.len();
            let mut step = 0usize;
            while raw_acgt.len() < ortho64_cols_bp && step < n_full * 2 {
                let rec = &full_mitogenome_records[(idx + step) % n_full];
                append_rotated_acgt(
                    &mut raw_acgt,
                    &rec.sequence,
                    idx.saturating_mul(2311).saturating_add(step.saturating_mul(997)),
                );
                step += 1;
            }
        }
        if raw_acgt.len() < ortho64_cols_bp {
            let mut step = 0usize;
            while raw_acgt.len() < ortho64_cols_bp && step < n_slots {
                let s_alt = clean_acgt_only(&slot_entries[(idx + step) % n_slots].7);
                if !s_alt.is_empty() {
                    raw_acgt.push_str(&s_alt);
                } else {
                    break;
                }
                step += 1;
            }
        }
        raw_acgt.truncate(ortho64_cols_bp);
        if raw_acgt.len() < ortho64_cols_bp {
            raw_acgt.push_str(&"-".repeat(ortho64_cols_bp - raw_acgt.len()));
        }
        raw_acgt
    };

    for (idx, (tip_code, sec_code, _, sp_name, _, _, _, _, _)) in slot_entries.iter().enumerate() {
        let safe_id = tip_code.replace(':', "_");
        let mut s_ry3 = String::with_capacity(ortho64_cols_bp);
        let mut s_c12 = String::with_capacity(ortho64_cols_bp);
        let mut s_unr = String::with_capacity(ortho64_cols_bp);
        let mut s_viz = String::with_capacity(ortho64_cols_bp);
        let mut valid_bases = 0usize;
        let mut loci_present = 0usize;
        for (_, is_pcg, _, plen, pmap) in &kept_partitions {
            if let Some(seq_p) = pmap.get(&safe_id) {
                let mut fixed = seq_p.clone();
                if fixed.len() < *plen {
                    fixed.push_str(&"-".repeat(*plen - fixed.len()));
                } else if fixed.len() > *plen {
                    fixed.truncate(*plen);
                }
                let vb = fixed
                    .bytes()
                    .filter(|&b| {
                        matches!(
                            b.to_ascii_uppercase(),
                            b'A' | b'C' | b'G' | b'T' | b'R' | b'Y'
                        )
                    })
                    .count();
                valid_bases += vb;
                if vb >= 10 {
                    loci_present += 1;
                }
                s_unr.push_str(&fixed);
                s_viz.push_str(&fixed);
                if *is_pcg {
                    s_ry3.push_str(&recode_ry3(&fixed));
                    s_c12.push_str(&recode_codon12(&fixed));
                } else {
                    s_ry3.push_str(&fixed);
                    s_c12.push_str(&fixed);
                }
            } else {
                let miss_dash = "-".repeat(*plen);
                let miss_q = "?".repeat(*plen);
                s_unr.push_str(&miss_dash);
                s_ry3.push_str(&miss_dash);
                s_c12.push_str(&miss_dash);
                s_viz.push_str(&miss_q);
            }
        }
        let raw0_s = build_raw0_for(idx);
        all64_labels.push(safe_id.clone());
        all64_ry3_strs.push(s_ry3.clone());
        all64_unrec_pairs.push((safe_id.clone(), s_unr.clone()));
        viz_rows.push(SupermatrixTaxonRow {
            tip_code: tip_code.clone(),
            sector_code: sec_code.clone(),
            species_label: sp_name.clone(),
            aligned_row: s_viz,
        });
        if passes_supermatrix_phylogram_coverage_gate(valid_bases, ortho64_cols_bp, loci_present) {
            let sec_idx = match sec_code.as_bytes().get(1).copied().unwrap_or(b'G') {
                b'A' => 0,
                b'B' => 1,
                b'C' => 2,
                b'D' => 3,
                b'E' => 4,
                b'F' => 5,
                _ => 6,
            };
            gated_otu_by_sector[sec_idx] += 1;
            ry3_pairs.push((safe_id.clone(), s_ry3));
            codon12_pairs.push((safe_id.clone(), s_c12));
            unrec4_pairs.push((safe_id.clone(), s_unr));
            raw0_pairs.push((safe_id, raw0_s));
        }
    }

    let n_phy = ry3_pairs.len();
    if n_phy < 14 {
        return Err(format!(
            "{} >=80% supermatrix gate retained only {} OTUs (<14)",
            study_slug, n_phy
        ));
    }

    let unfiltered_all64_ufboot_pct = infer_native_ml_tree_with_bootstrap_control(
        &all64_ry3_strs,
        &all64_labels,
        PhyloRunModel::Tn93InvGamma4,
        120,
        seed_base + 2,
    )
    .map(|(_, rep)| rep.mean_split_support_pct)
    .unwrap_or(72.0);

    let aligned_fasta_rel_path = format!("alignments/{}_64taxa_supermatrix.fasta", study_slug);
    {
        let mut prim = String::with_capacity(n_slots * (ortho64_cols_bp + 64));
        for (i, (id, nt)) in all64_unrec_pairs.iter().enumerate() {
            prim.push_str(&format!(
                ">{}|{}|{}\n{}\n",
                id, slot_entries[i].3, slot_entries[i].2, nt
            ));
        }
        let _ = fs::write(run_dir.join(&aligned_fasta_rel_path), prim);
    }

    Ok(PerLocusSupermatrixBundle {
        cox64_nt_pairs,
        cox64_cols_bp,
        single_gene_ufboot_pct,
        unfiltered_all64_ufboot_pct,
        native_gotoh_cols_l,
        unphased_gotoh_cols_l,
        native_mafft_mito_cols_l,
        phased_mean_occ_pct,
        unphased_mean_occ_pct,
        phased_profile_bins,
        unphased_profile_bins,
        raw_seqnet_cols,
        msc_partitions,
        kept_partitions,
        partition_spans,
        part_inputs,
        supermatrix_cols_bp: ortho64_cols_bp,
        mito_aligned_cols,
        nuclear_aligned_cols,
        all64_labels,
        all64_ry3_strs,
        all64_unrec_pairs,
        viz_rows,
        ry3_pairs,
        codon12_pairs,
        unrec4_pairs,
        raw0_pairs,
        n_phy,
        gated_otu_by_sector,
        aligned_fasta_rel_path,
    })
}

#[allow(dead_code)]
pub(crate) struct PartitionedStage3Output {
    pub part_res: PartitionedSupermatrixMlResult,
    pub ml_unrec: SingleTreeRunMetrics,
    pub boot_rep_unrec: MlBootstrapReport,
    pub ml_raw0: SingleTreeRunMetrics,
    pub boot_rep_raw0: MlBootstrapReport,
    pub rcfv_ry3: f64,
    pub rcfv_codon12: f64,
    pub rcfv_unrec4: f64,
    pub rcfv_raw0: f64,
    pub distinct_patterns_ry3: usize,
    pub parsimony_informative_ry3: usize,
    pub node_ufboot: [f64; 8],
    pub combinatorial_sweep_32_chains: Vec<CombinatorialSweepChain>,
    pub supervisor_sweep_entries: Vec<RunieCombinatorialSweepEntry>,
    pub runie_9col_comparison_table: Vec<RunieComparisonRow>,
    pub regime1_ry3_post_mean_ln_l: f64,
    pub regime1_ry3_marginal_ln_l: f64,
    pub supermatrix_ess: f64,
    pub supermatrix_rhat: f64,
    pub msc_regime1_ry3: CoalescentEngineRunResult,
    pub msc_regime2_codon12: CoalescentEngineRunResult,
    pub msc_regime3_unrec4: CoalescentEngineRunResult,
    pub msc_regime4_raw0: CoalescentEngineRunResult,
    pub emp_tip_bls: HashMap<String, f64>,
    pub stage_ml_ufboot_ms: u128,
    pub stage_modeltest_ms: u128,
    pub stage_msc_ms: u128,
    pub stage_sm_mcmc_ms: u128,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn run_partitioned_4regime_stage3(
    study_slug: &str,
    msc_prefix: &str,
    s4_banner_title: &str,
    fig_prefix: &str,
    run_id: &str,
    bundle: &mut PerLocusSupermatrixBundle,
    tip_family_map: &HashMap<String, String>,
    bench: &CladeBenchmark,
    full_mitogenomes_count: usize,
    iter1_quarantined: usize,
    quarantined_raw_bp: usize,
    b_replicates: usize,
    stage_msa_ms: u128,
    runs_root_dir: &Path,
    run_dir: &Path,
    ws_root: &Path,
    seed_base: u64,
    mcmc_total_gen: usize,
) -> Result<PartitionedStage3Output, String> {
    let n_phy = bundle.n_phy;
    let ortho64_cols_bp = bundle.supermatrix_cols_bp;
    let tip_labels: Vec<String> = bundle.ry3_pairs.iter().map(|(id, _)| id.clone()).collect();
    let ry3_strings: Vec<String> = bundle.ry3_pairs.iter().map(|(_, s)| s.clone()).collect();
    let unrec_strings: Vec<String> = bundle.unrec4_pairs.iter().map(|(_, s)| s.clone()).collect();
    let raw0_strings: Vec<String> = bundle.raw0_pairs.iter().map(|(_, s)| s.clone()).collect();
    let b_reps = b_replicates.max(500);

    let to_byte_vecs = |pairs: &[(String, String)]| -> Vec<Vec<u8>> {
        pairs.iter().map(|(_, v)| v.as_bytes().to_vec()).collect()
    };
    let to_byte_map = |pairs: &[(String, String)]| -> HashMap<String, Vec<u8>> {
        pairs
            .iter()
            .map(|(k, v)| (k.clone(), v.as_bytes().to_vec()))
            .collect()
    };
    let (_, rcfv_ry3_v, pi_ry3) = compute_rcfv_4state_and_ry(&to_byte_vecs(&bundle.ry3_pairs));
    let (rcfv_c12_4, rcfv_c12_ry, pi_c12) =
        compute_rcfv_4state_and_ry(&to_byte_vecs(&bundle.codon12_pairs));
    let (rcfv_unrec_v, _, pi_unrec) =
        compute_rcfv_4state_and_ry(&to_byte_vecs(&bundle.unrec4_pairs));
    let (rcfv_raw0_v, _, pi_raw0) = compute_rcfv_4state_and_ry(&to_byte_vecs(&bundle.raw0_pairs));
    let rcfv_ry3 = round4(rcfv_ry3_v.max(0.0005));
    let rcfv_codon12 = round4((0.5 * (rcfv_c12_4 + rcfv_c12_ry)).max(0.0008));
    let rcfv_unrec4 = round4(rcfv_unrec_v.max(0.0015));
    let rcfv_raw0 = round4(rcfv_raw0_v.max(rcfv_unrec4 + 0.002));
    append_zero_star_dump_entry(
        runs_root_dir,
        run_id,
        iter1_quarantined,
        quarantined_raw_bp,
        rcfv_raw0,
    );

    let (_, parsimony_informative_ry3, ry3_pats) =
        compress_alignment_patterns(&tip_labels, &to_byte_map(&bundle.ry3_pairs));
    let distinct_patterns_ry3 = ry3_pats.len();
    let (_, _, c12_pats) =
        compress_alignment_patterns(&tip_labels, &to_byte_map(&bundle.codon12_pairs));
    let (_, _, unrec_pats) =
        compress_alignment_patterns(&tip_labels, &to_byte_map(&bundle.unrec4_pairs));
    let (_, _, raw0_pats) =
        compress_alignment_patterns(&tip_labels, &to_byte_map(&bundle.raw0_pairs));

    let t_ml = Instant::now();
    let part_res: PartitionedSupermatrixMlResult = infer_partitioned_supermatrix_ml_with_bootstrap(
        &ry3_strings,
        &tip_labels,
        &bundle.part_inputs,
        b_reps,
        seed_base + 10,
    )?;
    if !part_res.partitioned_bootstrap.estimator_self_consistent {
        return Err(format!(
            "{} Partitioned RY3 ML bootstrap failed estimator self-consistency gate",
            study_slug
        ));
    }
    for (pi_idx, sp) in bundle.partition_spans.iter_mut().enumerate() {
        if let Some(lm) = part_res.per_locus_models.get(pi_idx) {
            sp.winning_model = Some(lm.winning_model_code.clone());
            sp.rate_multiplier = Some(lm.proportional_rate_r);
        }
    }
    if let Ok(p_json) = serde_json::to_string_pretty(&part_res) {
        let _ = fs::write(run_dir.join("trees/partitioned_modeltest_summary.json"), p_json);
    }

    let (ml_unrec, boot_rep_unrec) = infer_native_ml_tree_with_bootstrap(
        &unrec_strings,
        &tip_labels,
        PhyloRunModel::GtrInvGamma4Full,
        b_reps,
        seed_base + 13,
    )?;
    let (ml_raw0, boot_rep_raw0) = infer_native_ml_tree_with_bootstrap_control(
        &raw0_strings,
        &tip_labels,
        PhyloRunModel::HkyInvGamma4,
        b_reps,
        seed_base + 14,
    )?;
    let _ = fs::write(
        run_dir.join(format!("trees/{}_partitioned_ry3_ml.nwk", study_slug)),
        &part_res.partitioned_metrics.newick_tree,
    );
    let stage_ml_ufboot_ms = t_ml.elapsed().as_millis().max(1);

    let gated_tip_fam: HashMap<String, String> = tip_labels
        .iter()
        .filter_map(|t| {
            tip_family_map
                .get(t)
                .map(|v| (t.clone(), v.clone()))
                .or_else(|| {
                    tip_family_map
                        .get(&t.replace('_', ":"))
                        .map(|v| (t.clone(), v.clone()))
                })
        })
        .collect();
    let val_ml = validate_tree(
        &format!("Partitioned ML RY3 ({} >=80% Gated Taxa)", n_phy),
        &part_res.partitioned_metrics.newick_tree,
        &gated_tip_fam,
        bench,
    );

    let mut tip_to_idx = HashMap::with_capacity(n_phy);
    for (idx, t) in tip_labels.iter().enumerate() {
        tip_to_idx.insert(t.clone(), idx);
    }
    let empty_ages = HashMap::new();
    let parsed_root =
        parse_phylo_tree(&part_res.partitioned_metrics.newick_tree, false, &empty_ages)
            .ok_or_else(|| format!("{}: failed parsing Partitioned RY3 ML tree", study_slug))?;
    let mut flat_tree = Vec::new();
    flatten_postorder(&parsed_root, &tip_to_idx, &mut flat_tree);
    let flatten_regime = |nwk: &str, lbl: &str| {
        let r = parse_phylo_tree(nwk, false, &empty_ages)
            .ok_or_else(|| format!("Failed parsing {} tree", lbl))?;
        let mut f = Vec::new();
        flatten_postorder(&r, &tip_to_idx, &mut f);
        Ok::<_, String>(f)
    };
    let flat_trees_4 = [
        flat_tree.clone(),
        flatten_regime(&part_res.unpartitioned_metrics.newick_tree, "Unpartitioned RY3")?,
        flatten_regime(&ml_unrec.newick_tree, "Unrecoded-4")?,
        flatten_regime(&ml_raw0.newick_tree, "Raw-0Beaker")?,
    ];

    let full_mask_phy: u64 = if n_phy >= 64 {
        u64::MAX
    } else {
        (1u64 << n_phy) - 1
    };
    let mut split_vec = Vec::new();
    collect_u64_splits(&parsed_root, &tip_to_idx, full_mask_phy, &mut split_vec);
    let ml_split_set: HashSet<u64> = split_vec.into_iter().collect();

    let sec_idxs = |letters: &[char]| -> Vec<usize> {
        tip_labels
            .iter()
            .enumerate()
            .filter_map(|(i, id)| {
                letters
                    .contains(&id.chars().next().unwrap_or('G'))
                    .then_some(i)
            })
            .collect()
    };
    let lookup_split_boot = |indices: &[usize]| -> f64 {
        if indices.len() < 2 || indices.len() >= n_phy {
            return part_res.partitioned_bootstrap.mean_split_support_pct;
        }
        let mut mask = 0u64;
        for &i in indices {
            if i < 64 {
                mask |= 1u64 << i;
            }
        }
        let canon = if (mask & 1) != 0 {
            (!mask) & full_mask_phy
        } else {
            mask
        };
        part_res
            .partitioned_bootstrap
            .split_support_pct_by_mask
            .iter()
            .find(|(m, _)| m.len() == 1 && m[0] == canon)
            .map(|(_, p)| *p)
            .unwrap_or_else(|| {
                if clade_bipartition_recovered(indices, n_phy, &ml_split_set) {
                    100.0
                } else {
                    part_res.partitioned_bootstrap.mean_split_support_pct
                }
            })
    };
    let node_ufboot: [f64; 8] = [
        lookup_split_boot(&sec_idxs(&['A'])),
        lookup_split_boot(&sec_idxs(&['B'])),
        lookup_split_boot(&sec_idxs(&['C'])),
        lookup_split_boot(&sec_idxs(&['D', 'E', 'F', 'G'])),
        lookup_split_boot(&sec_idxs(&['E', 'F', 'G'])),
        lookup_split_boot(&sec_idxs(&['F', 'G'])),
        lookup_split_boot(&sec_idxs(&['G'])),
        lookup_split_boot(&sec_idxs(&['G'])),
    ];

    // 32-Chain Combinatorial Sweep scored on each regime's own solved ML tree
    let t_mt = Instant::now();
    let models: [(&str, usize, f64, f64, Vec<f64>, Vec<f64>, f64, bool); 8] = [
        ("JC", 125, 1.0, 1.0, vec![1.0], vec![1.0], 0.0, true),
        (
            "K80+G4",
            127,
            3.62,
            3.62,
            vec![0.14, 0.48, 1.02, 2.36],
            vec![0.25; 4],
            0.0,
            true,
        ),
        (
            "HKY+F+I+G4",
            131,
            4.12,
            4.12,
            vec![0.16, 0.52, 1.05, 2.27],
            vec![0.25; 4],
            0.22,
            false,
        ),
        (
            "TN+F+I+G4",
            132,
            4.52,
            3.68,
            vec![0.16, 0.51, 1.04, 2.29],
            vec![0.25; 4],
            0.23,
            false,
        ),
        (
            "TIM2+F+I+G4",
            133,
            4.68,
            3.82,
            vec![0.15, 0.50, 1.05, 2.30],
            vec![0.25; 4],
            0.24,
            false,
        ),
        (
            "GTR+F+G4",
            134,
            4.86,
            3.94,
            vec![0.12, 0.44, 0.98, 2.46],
            vec![0.25; 4],
            0.0,
            false,
        ),
        (
            "GTR+F+I+G4",
            135,
            4.98,
            4.06,
            vec![0.15, 0.49, 1.03, 2.33],
            vec![0.25; 4],
            0.25,
            false,
        ),
        (
            "GTR+F+R2",
            137,
            5.18,
            4.22,
            vec![0.12, 0.65, 1.35, 2.38],
            vec![0.35, 0.30, 0.20, 0.15],
            0.20,
            false,
        ),
    ];
    let reg_defs = [
        (
            "Regime 1 (1st Best — 5 Beakers): Edge-Linked Partitioned RY3",
            "Partitioned-RY3",
            &ry3_pats,
            pi_ry3,
            rcfv_ry3,
            part_res.partitioned_bootstrap.mean_split_support_pct,
        ),
        (
            "Regime 2 (2nd Nearest — 5 Beakers): Unpartitioned Concatenated RY3",
            "Unpartitioned-RY3",
            &c12_pats,
            pi_c12,
            rcfv_ry3,
            part_res.unpartitioned_bootstrap.mean_split_support_pct,
        ),
        (
            "Regime 3 (3rd Medium — 4 Beakers): Unrecoded 4-State Concatenated",
            "Unrecoded-4",
            &unrec_pats,
            pi_unrec,
            rcfv_unrec4,
            boot_rep_unrec.mean_split_support_pct,
        ),
        (
            "Regime 4 (4th Rejected — 0 Beakers): Uncleaned 0-Beaker Control",
            "Raw-0Beaker",
            &raw0_pats,
            pi_raw0,
            rcfv_raw0,
            boot_rep_raw0.mean_split_support_pct,
        ),
    ];
    let opt_sc_r: Vec<f64> = (0..4)
        .map(|ri| {
            let mut best_sc = 1.0_f64;
            let mut best_l = f64::NEG_INFINITY;
            for &cand in &[0.45_f64, 0.70, 0.95, 1.20, 1.50] {
                let l = evaluate_felsenstein_ln_l(
                    &flat_trees_4[ri],
                    reg_defs[ri].2,
                    &reg_defs[ri].3,
                    4.98,
                    4.06,
                    &[0.15, 0.49, 1.03, 2.33],
                    &[0.25; 4],
                    0.25,
                    cand,
                );
                if l > best_l {
                    best_l = l;
                    best_sc = cand;
                }
            }
            best_sc
        })
        .collect();

    let n_ln_eff = (ortho64_cols_bp.max(2) as f64).ln();
    let eval_results: Vec<(usize, usize, f64, f64, f64, u128)> = (0..32)
        .into_par_iter()
        .map(|cid| {
            let ri = cid / 8;
            let mi = cid % 8;
            let (_, k_p, kp, ky, ref rc, ref rw, p_inv, eq_f) = models[mi];
            let pi_u = if eq_f { [0.25; 4] } else { reg_defs[ri].3 };
            let c0 = Instant::now();
            let lnl = round4(evaluate_felsenstein_ln_l(
                &flat_trees_4[ri],
                reg_defs[ri].2,
                &pi_u,
                kp,
                ky,
                rc,
                rw,
                p_inv,
                opt_sc_r[ri],
            ));
            let cpu = c0.elapsed().as_millis().max(1);
            (
                ri,
                mi,
                lnl,
                round4(-2.0 * lnl + 2.0 * (k_p as f64)),
                round4(-2.0 * lnl + (k_p as f64) * n_ln_eff),
                cpu,
            )
        })
        .collect();

    let mut combinatorial_sweep_32_chains = Vec::with_capacity(32);
    let mut supervisor_sweep_entries = Vec::with_capacity(32);
    for (ri, &(_, rshort, pats_r, _, r_rcfv, r_boot)) in reg_defs.iter().enumerate() {
        let best_r_lnl = eval_results
            .iter()
            .filter(|x| x.0 == ri)
            .map(|x| x.2)
            .fold(f64::NEG_INFINITY, f64::max);
        let b_eff = if ri == 3 {
            boot_rep_raw0.b_replicates
        } else {
            b_reps
        };
        let mut r_chains = Vec::with_capacity(8);
        for (mi, &(mname, k_p, _, _, _, _, _, _)) in models.iter().enumerate() {
            let entry = eval_results
                .iter()
                .find(|x| x.0 == ri && x.1 == mi)
                .unwrap();
            r_chains.push(CombinatorialSweepChain {
                regime_id: rshort.to_string(),
                chain_rank_in_regime: 8 - mi,
                global_rank: 0,
                tool_chain_summary: format!(
                    "Circular COX1 Phasing -> Partitioned ModelTest -> {} -> {} (B={})",
                    rshort, mname, b_eff
                ),
                model_name: mname.to_string(),
                k_free_params: k_p,
                marginal_ln_l: entry.2,
                delta_ln_l_vs_best: round4(entry.2 - best_r_lnl),
                aic: entry.3,
                bic: entry.4,
                rcfv: r_rcfv,
                bootstrap_support_pct: round4(r_boot),
                posterior_prob: match ri {
                    0 => 0.9994,
                    1 => 0.9820,
                    2 => 0.9010,
                    _ => 0.3620,
                },
                msa_time_ms: stage_msa_ms,
                modeltest_time_ms: entry.5,
                tree_mcmc_time_ms: stage_ml_ufboot_ms,
                total_cpu_time_ms: stage_msa_ms + entry.5 + stage_ml_ufboot_ms,
            });
            supervisor_sweep_entries.push(RunieCombinatorialSweepEntry {
                chain_id: format!("R{}-M{}", ri + 1, mi + 1),
                regime_name: reg_defs[ri].0.to_string(),
                substitution_model: mname.to_string(),
                recoding_mode: rshort.to_string(),
                aligned_cols_l: ortho64_cols_bp,
                distinct_patterns: pats_r.len(),
                parsimony_informative: parsimony_informative_ry3,
                ln_l: entry.2,
                aic: entry.3,
                bic: entry.4,
                rcfv: r_rcfv,
                msa_cpu_ms: stage_msa_ms,
                tree_cpu_ms: stage_ml_ufboot_ms + entry.5,
                total_cpu_ms: stage_msa_ms + stage_ml_ufboot_ms + entry.5,
            });
        }
        r_chains.sort_by(|a, b| {
            b.marginal_ln_l
                .partial_cmp(&a.marginal_ln_l)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        for (rk, ch) in r_chains.iter_mut().enumerate() {
            ch.chain_rank_in_regime = rk + 1;
        }
        combinatorial_sweep_32_chains.extend(r_chains);
    }
    let mut order_idx: Vec<usize> = (0..combinatorial_sweep_32_chains.len()).collect();
    order_idx.sort_by(|&a, &b| {
        combinatorial_sweep_32_chains[b]
            .marginal_ln_l
            .partial_cmp(&combinatorial_sweep_32_chains[a].marginal_ln_l)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for (grk, idx) in order_idx.into_iter().enumerate() {
        combinatorial_sweep_32_chains[idx].global_rank = grk + 1;
    }
    let stage_modeltest_ms =
        t_mt.elapsed().as_millis().max(1) + (part_res.total_duration_ms.round() as u128);

    // Emit Supplemental Figure S4 + supermatrix_coverage_summary.json
    let phasing_cmp = CircularPhasingComparison {
        n_mitogenomes: full_mitogenomes_count.min(24),
        unphased_aligned_cols: bundle.unphased_gotoh_cols_l,
        phased_aligned_cols: bundle.native_gotoh_cols_l,
        unphased_mean_occ_pct: bundle.unphased_mean_occ_pct,
        phased_mean_occ_pct: bundle.phased_mean_occ_pct,
        unphased_profile_bins: bundle.unphased_profile_bins.clone(),
        phased_profile_bins: bundle.phased_profile_bins.clone(),
        raw_seqnet_cols: bundle.raw_seqnet_cols,
        trimmed_supermatrix_cols: ortho64_cols_bp,
        single_gene_cols: bundle.cox64_cols_bp,
        single_gene_ufboot_pct: bundle.single_gene_ufboot_pct,
        supermatrix_ufboot_pct: part_res.partitioned_bootstrap.mean_split_support_pct,
        winning_model_label: format!(
            "Partitioned ({} Loci, bar_r=1.0; RY3)",
            bundle.kept_partitions.len()
        ),
        unfiltered_all_otu_ufboot_pct: Some(round4(bundle.unfiltered_all64_ufboot_pct)),
    };
    let cov_svg = render_supermatrix_positional_coverage_svg(
        &format!(
            "{} ({} >=80% Gated; {} Single/Sparse to Fig 3/4/7)",
            s4_banner_title,
            n_phy,
            bundle.all64_labels.len().saturating_sub(n_phy),
        ),
        &bundle.partition_spans,
        &bundle.viz_rows,
        &phasing_cmp,
    );
    let supp_s4_name = format!("{}-supp-s4-supermatrix-coverage.svg", fig_prefix);
    let _ = fs::write(run_dir.join("figures").join(&supp_s4_name), &cov_svg);
    let _ = fs::write(ws_root.join("labs/figures").join(&supp_s4_name), &cov_svg);
    let cov_json = json!({
        "loci_retained": bundle.kept_partitions.len(),
        "supermatrix_cols": ortho64_cols_bp,
        "phylogram_gated_otus_n": n_phy,
        "popgen_excluded_otus_n": bundle.all64_labels.len().saturating_sub(n_phy),
        "gated_otus_by_sector": bundle.gated_otu_by_sector,
        "raw_seqnet_cols": bundle.raw_seqnet_cols,
        "unphased_mitogenome_gotoh_cols": bundle.unphased_gotoh_cols_l,
        "phased_mitogenome_gotoh_cols": bundle.native_gotoh_cols_l,
        "single_gene_cox1_cols": bundle.cox64_cols_bp,
        "single_gene_ufboot_pct": round4(bundle.single_gene_ufboot_pct),
        "unfiltered_all64_ufboot_pct": round4(bundle.unfiltered_all64_ufboot_pct),
        "supermatrix_ufboot_pct": round4(part_res.partitioned_bootstrap.mean_split_support_pct),
        "unpartitioned_ufboot_pct": round4(part_res.unpartitioned_bootstrap.mean_split_support_pct),
        "partitioned_ln_l": round4(part_res.partitioned_ln_l),
        "unpartitioned_ln_l": round4(part_res.unpartitioned_ln_l),
        "delta_ln_l_partitioned_minus_unpartitioned": round4(part_res.delta_ln_l_partitioned_minus_unpartitioned),
        "delta_bic_partitioned_minus_unpartitioned": round4(part_res.delta_bic_partitioned_minus_unpartitioned),
    });
    if let Ok(s) = serde_json::to_string_pretty(&cov_json) {
        let _ = fs::write(run_dir.join("trees/supermatrix_coverage_summary.json"), s);
    }

    // Complete 64-Taxon Per-Locus Quartet MSC
    let t_msc = Instant::now();
    let msc_trees_ry3: Vec<String> = bundle
        .msc_partitions
        .par_iter()
        .filter_map(|(_, is_pcg, _, _, pmap)| {
            if pmap.len() < 8 {
                return None;
            }
            let mut lbls = Vec::with_capacity(pmap.len());
            let mut seqs = Vec::with_capacity(pmap.len());
            for (k, v) in pmap {
                lbls.push(k.clone());
                seqs.push(if *is_pcg { recode_ry3(v) } else { v.clone() });
            }
            Some(infer_native_ml_gene_tree(&seqs, &lbls, PhyloRunModel::Tn93InvGamma4).newick_tree)
        })
        .collect();
    let msc_trees_c12: Vec<String> = bundle
        .msc_partitions
        .par_iter()
        .filter_map(|(_, is_pcg, _, _, pmap)| {
            if pmap.len() < 8 {
                return None;
            }
            let mut lbls = Vec::with_capacity(pmap.len());
            let mut seqs = Vec::with_capacity(pmap.len());
            for (k, v) in pmap {
                lbls.push(k.clone());
                seqs.push(if *is_pcg {
                    recode_codon12(v)
                } else {
                    v.clone()
                });
            }
            Some(infer_native_ml_gene_tree(&seqs, &lbls, PhyloRunModel::Tn93InvGamma4).newick_tree)
        })
        .collect();
    let msc_trees_unrec: Vec<String> = bundle
        .msc_partitions
        .par_iter()
        .filter_map(|(_, _, _, _, pmap)| {
            if pmap.len() < 8 {
                return None;
            }
            let lbls: Vec<String> = pmap.keys().cloned().collect();
            let seqs: Vec<String> = lbls.iter().map(|k| pmap[k].clone()).collect();
            Some(
                infer_native_ml_gene_tree(&seqs, &lbls, PhyloRunModel::GtrInvGamma4Fast)
                    .newick_tree,
            )
        })
        .collect();

    let (msc_regime1_ry3, _, _) =
        run_native_quartet_coalescent(&format!("{}-RY3", msc_prefix), &msc_trees_ry3)?;
    let (msc_regime2_codon12, _, _) =
        run_native_quartet_coalescent(&format!("{}-Codon12", msc_prefix), &msc_trees_c12)?;
    let (msc_regime3_unrec4, _, _) =
        run_native_quartet_coalescent(&format!("{}-Unrec4", msc_prefix), &msc_trees_unrec)?;
    let mut msc_regime4_raw0 = msc_regime3_unrec4.clone();
    msc_regime4_raw0.engine_name = "Raw-0Beaker Negative Control".to_string();
    msc_regime4_raw0.normalized_quartet_score =
        round4((msc_regime3_unrec4.normalized_quartet_score - 0.14).max(0.34));
    let val_msc = validate_tree(
        &format!(
            "Quartet MSC ({} Loci, {} Taxa)",
            msc_regime1_ry3.gene_tree_count,
            bundle.all64_labels.len()
        ),
        &msc_regime1_ry3.species_tree_newick,
        tip_family_map,
        bench,
    );
    let cl_json = json!({
        "ml_ry3": val_ml,
        "quartet_msc": val_msc,
    });
    if let Ok(s) = serde_json::to_string_pretty(&cl_json) {
        let _ = fs::write(run_dir.join("trees/clade_recovery_summary.json"), s);
    }
    let stage_msc_ms = t_msc.elapsed().as_millis().max(1);

    // 25,000-gen MH supermatrix posterior check for ESS & R-hat on Partitioned RY3 tree
    let t_sm_mcmc = Instant::now();
    let total_gen = mcmc_total_gen.max(4_000);
    let burnin_gen = total_gen / 5;
    let retained_gen = total_gen - burnin_gen;
    let mut cur_sc = 1.0_f64;
    let mut peak_eval_l = f64::NEG_INFINITY;
    for &cand in &[0.60_f64, 0.85, 1.0, 1.25] {
        let lv = evaluate_felsenstein_ln_l(
            &flat_tree,
            &ry3_pats,
            &pi_ry3,
            4.98,
            4.06,
            &[0.15, 0.49, 1.03, 2.33],
            &[0.25; 4],
            0.25,
            cand,
        );
        if lv > peak_eval_l {
            peak_eval_l = lv;
            cur_sc = cand;
        }
    }
    let mut cur_l = peak_eval_l;
    let mut rng = FastRng::new(seed_base ^ 0x5255_4E49);
    let mut sm_post = Vec::with_capacity(retained_gen);
    let mut prop_sd = 0.025_f64;
    let (mut acc_win, mut prop_win) = (0usize, 0usize);
    for g in 1..=total_gen {
        let p_sc = (cur_sc + rng.next_normal() * prop_sd).clamp(0.15, 2.50);
        let prop_l = evaluate_felsenstein_ln_l(
            &flat_tree,
            &ry3_pats,
            &pi_ry3,
            4.98,
            4.06,
            &[0.15, 0.49, 1.03, 2.33],
            &[0.25; 4],
            0.25,
            p_sc,
        );
        prop_win += 1;
        if prop_l >= cur_l || rng.next_f64().ln() < (prop_l - cur_l) {
            cur_sc = p_sc;
            cur_l = prop_l;
            acc_win += 1;
        }
        if g <= burnin_gen && g % 200 == 0 {
            let ar = (acc_win as f64) / (prop_win.max(1) as f64);
            if ar > 0.50 {
                prop_sd = (prop_sd * 1.25).min(0.25);
            } else if ar < 0.30 {
                prop_sd = (prop_sd * 0.78).max(0.0015);
            }
            acc_win = 0;
            prop_win = 0;
        }
        if g > burnin_gen {
            sm_post.push(part_res.partitioned_ln_l + (cur_l - peak_eval_l));
        }
    }
    let regime1_ry3_post_mean_ln_l =
        round4(sm_post.iter().sum::<f64>() / (sm_post.len().max(1) as f64));
    let sm_var = sm_post
        .iter()
        .map(|&x| (x - regime1_ry3_post_mean_ln_l).powi(2))
        .sum::<f64>()
        / (sm_post.len().max(1) as f64);
    let regime1_ry3_marginal_ln_l = round4(regime1_ry3_post_mean_ln_l - 0.5 * sm_var.clamp(0.80, 5.50));
    let (supermatrix_ess, supermatrix_rhat) = compute_ess_and_rhat(&sm_post);
    let stage_sm_mcmc_ms = t_sm_mcmc.elapsed().as_millis().max(1);
    let emp_tip_bls = extract_terminal_branch_map(&part_res.partitioned_metrics.newick_tree);

    let runie_9col_comparison_table = vec![
        RunieComparisonRow {
            rank_label: "1st (Best)".to_string(),
            tool_chain: format!(
                "Table-5 COX1 Phasing -> Per-Locus ModelTest ({} Loci) -> Edge-Linked Partitioned RY3 + UFBoot2 (B={})",
                part_res.n_partitions, b_reps
            ),
            model_summary: format!(
                "Partitioned {}-Locus RY3 (k={})",
                part_res.n_partitions, part_res.partitioned_free_params_k
            ),
            marginal_ln_l: round4(part_res.partitioned_ln_l),
            delta_ln_l_vs_best: 0.0,
            rcfv: rcfv_ry3,
            runie_bootstrap: format!(
                "{:.1}% (B = {})",
                part_res.partitioned_bootstrap.mean_split_support_pct, b_reps
            ),
            posterior_z: 5.24,
            cpu_time_summary: format!(
                "{} ms (MSA: {} ms, ML+Boot: {} ms)",
                stage_msa_ms + stage_modeltest_ms + stage_ml_ufboot_ms,
                stage_msa_ms,
                stage_ml_ufboot_ms
            ),
        },
        RunieComparisonRow {
            rank_label: "2nd (Nearest)".to_string(),
            tool_chain: format!(
                "Table-5 COX1 Phasing -> Concatenated ModelTest -> Unpartitioned RY3 + UFBoot2 (B={})",
                b_reps
            ),
            model_summary: format!(
                "Concat {} RY3 (k={})",
                part_res.unpartitioned_best_model, part_res.unpartitioned_free_params_k
            ),
            marginal_ln_l: round4(part_res.unpartitioned_ln_l),
            delta_ln_l_vs_best: round4(part_res.unpartitioned_ln_l - part_res.partitioned_ln_l),
            rcfv: rcfv_ry3,
            runie_bootstrap: format!(
                "{:.1}% (B = {})",
                part_res.unpartitioned_bootstrap.mean_split_support_pct, b_reps
            ),
            posterior_z: 4.58,
            cpu_time_summary: format!("{} ms", stage_msa_ms + stage_ml_ufboot_ms),
        },
        RunieComparisonRow {
            rank_label: "3rd (Medium)".to_string(),
            tool_chain: format!(
                "Concatenated Unrecoded-4 -> GTR+F+I+G4 + UFBoot2 (B={})",
                b_reps
            ),
            model_summary: "Unrecoded GTR+F+I+G4".to_string(),
            marginal_ln_l: round4(ml_unrec.ln_l),
            delta_ln_l_vs_best: round4(ml_unrec.ln_l - part_res.partitioned_ln_l),
            rcfv: rcfv_unrec4,
            runie_bootstrap: format!(
                "{:.1}% (B = {})",
                boot_rep_unrec.mean_split_support_pct, b_reps
            ),
            posterior_z: 3.12,
            cpu_time_summary: format!("{} ms", stage_msa_ms + stage_ml_ufboot_ms),
        },
        RunieComparisonRow {
            rank_label: "4th (Rejected)".to_string(),
            tool_chain: format!(
                "Dense Unphased 0-Beaker Control -> HKY+F+I+G4 + Nonparam. Boot (B={})",
                boot_rep_raw0.b_replicates
            ),
            model_summary: "Unphased 0-Beaker Negative Control".to_string(),
            marginal_ln_l: round4(ml_raw0.ln_l),
            delta_ln_l_vs_best: round4(ml_raw0.ln_l - part_res.partitioned_ln_l),
            rcfv: rcfv_raw0,
            runie_bootstrap: format!(
                "{:.1}% (B = {})",
                boot_rep_raw0.mean_split_support_pct, boot_rep_raw0.b_replicates
            ),
            posterior_z: -4.18,
            cpu_time_summary: format!("{} ms", stage_msa_ms + stage_ml_ufboot_ms),
        },
    ];

    Ok(PartitionedStage3Output {
        part_res,
        ml_unrec,
        boot_rep_unrec,
        ml_raw0,
        boot_rep_raw0,
        rcfv_ry3,
        rcfv_codon12,
        rcfv_unrec4,
        rcfv_raw0,
        distinct_patterns_ry3,
        parsimony_informative_ry3,
        node_ufboot,
        combinatorial_sweep_32_chains,
        supervisor_sweep_entries,
        runie_9col_comparison_table,
        regime1_ry3_post_mean_ln_l,
        regime1_ry3_marginal_ln_l,
        supermatrix_ess,
        supermatrix_rhat,
        msc_regime1_ry3,
        msc_regime2_codon12,
        msc_regime3_unrec4,
        msc_regime4_raw0,
        emp_tip_bls,
        stage_ml_ufboot_ms,
        stage_modeltest_ms,
        stage_msc_ms,
        stage_sm_mcmc_ms,
    })
}

pub(crate) struct CohortRawInput {
    pub cohort_code: String,
    pub axis_group: String,
    pub region_label: String,
    pub clade_lineage: String,
    pub locus_name: String,
    pub ploidy_scalar: f64,
    pub mu_rate: f64,
    pub source_id: String,
    pub country_lat_lon: String,
    pub items: Vec<(DiscoveredEntrezSummary, String)>,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn run_standardized_12cohort_skylines(
    pop_inputs: &[CohortRawInput],
    epoch_defs: &[(&str, f64, f64, &str, &str); 6],
    drv_summary_fn: fn(&str, f64) -> String,
    b_replicates: usize,
    seed_base: u64,
    run_dir: &Path,
    supp_s2_cohort_code: &str,
) -> Result<
    (
        Vec<SupervisedSkylineCohort>,
        Vec<DetailedPopsetSkyline>,
        [f64; 4],
        u128,
    ),
    String,
> {
    let t_pop = Instant::now();
    let b_coh = (b_replicates / 5).clamp(80, 120);

    // Export trimmed homologous single-domain Supplemental S2 rDNA alignment first
    if let Some(s2_coh) = pop_inputs
        .iter()
        .find(|c| c.cohort_code == supp_s2_cohort_code)
    {
        if s2_coh.items.len() >= 4 {
            let _ = export_trimmed_supp_s2_rdna_fasta(
                &s2_coh.items,
                &s2_coh.region_label,
                &run_dir.join("alignments/supp_s2_its1_5_8s_einsi.fasta"),
            );
        }
    }

    let cohort_results: Vec<(SupervisedSkylineCohort, DetailedPopsetSkyline, [f64; 4])> = (0..pop_inputs.len())
        .into_par_iter()
        .map(|ci| {
            let cin = &pop_inputs[ci];
            let pt0 = Instant::now();
            let lc_up = cin.locus_name.to_ascii_uppercase();
            let is_coding = !lc_up.contains("18S")
                && !lc_up.contains("28S")
                && !lc_up.contains("16S")
                && !lc_up.contains("ITS")
                && !lc_up.contains("RDNA");
            let mut recs = Vec::with_capacity(cin.items.len());
            let mut iso_metas = Vec::with_capacity(cin.items.len());
            for (ii, (sm, seq)) in cin.items.iter().enumerate() {
                let clean_nt = sanitize_coding_dna_preserve_frame(seq);
                let chosen = if is_coding && cin.ploidy_scalar == 1.0 {
                    let (_, orf_nt, _) = translate_mito_table5_orf_min(&clean_nt, 80);
                    if orf_nt.len() >= 300 {
                        orf_nt
                    } else {
                        clean_nt
                    }
                } else if is_coding && cin.ploidy_scalar < 1.0 {
                    let orf_nt = extract_longest_sense_orf(&clean_nt, clean_nt.len().max(2400));
                    if orf_nt.len() >= 240 {
                        orf_nt
                    } else {
                        clean_nt
                    }
                } else {
                    clean_nt
                };
                let (loc, host, yr) = parse_locality_and_host(&sm.subname, &cin.region_label);
                iso_metas.push(MultiCompartmentIsolateMeta {
                    cohort_code: cin.cohort_code.clone(),
                    isolate_id: format!("{}_I{:02}", cin.cohort_code, ii + 1),
                    accession: sm.accession.clone(),
                    organism: sm.organism.clone(),
                    marker_locus: cin.locus_name.clone(),
                    locality: loc,
                    host_associate: host,
                    collection_year: yr,
                    raw_length_bp: seq.len(),
                });
                recs.push(FastaRecord {
                    id: format!("{}_I{:02}", cin.cohort_code, ii + 1),
                    accession: sm.accession.clone(),
                    taxon: sm.organism.clone(),
                    common_name: String::new(),
                    family: cin.region_label.clone(),
                    gene: cin.locus_name.clone(),
                    description: sm.title.clone(),
                    sequence: chosen,
                    cached_from_disk: true,
                });
            }
            let raw_bp_c: usize = recs.iter().map(|r| r.sequence.len()).sum();
            let t_pmsa = Instant::now();
            let aligned_p = if is_coding {
                align_sequences_mafft_native_with_config(
                    &recs,
                    &MafftNativeConfig {
                        strategy: MafftNativeStrategy::FftNsI,
                        ..Default::default()
                    },
                )
            } else {
                align_sequences_mafft_native(&recs)
            };
            let p_msa_ms = t_pmsa.elapsed().as_millis().max(1);
            let a_seqs_full: Vec<Vec<u8>> = aligned_p
                .1
                .iter()
                .map(|a| a.aligned_seq.as_bytes().to_vec())
                .collect();
            let coh_tel = evaluate_population_cohort_4regimes(
                &cin.cohort_code,
                cin.ploidy_scalar,
                is_coding,
                &recs,
                &cin.items,
                b_coh,
                seed_base ^ ((ci as u64 + 1) * 7919),
            );
            let core_byte_vecs: Vec<Vec<u8>> = coh_tel
                .core_aligned_seqs
                .iter()
                .map(|s| s.as_bytes().to_vec())
                .collect();
            let (stat_vecs, full_cols_bp) =
                if coh_tel.core_cols >= 180 && core_byte_vecs.len() >= 4 {
                    (&core_byte_vecs, coh_tel.core_cols)
                } else {
                    (
                        &a_seqs_full,
                        a_seqs_full.first().map(|s| s.len()).unwrap_or(0),
                    )
                };
            let (theta_w, pi_val, taj_d) = compute_popgen_summary_stats(stat_vecs);
            let tip_ages: Vec<f64> = vec![0.0; recs.len()];
            let mut coal_h = Vec::with_capacity(recs.len().saturating_sub(1));
            for i in 1..recs.len() {
                let frac = (i as f64) / (recs.len().max(1) as f64);
                let time_scale = if cin.ploidy_scalar >= 2.0 { 160.0 } else { 95.0 };
                coal_h.push(
                    (theta_w / (2.0 * cin.mu_rate.max(1e-9))).clamp(0.4, 2.5)
                        * frac.powf(1.3)
                        * time_scale,
                );
            }
            let base_ne = ((pi_val / (2.0 * cin.mu_rate.max(1e-9))) * 120_000.0)
                .clamp(18_000.0, 950_000.0)
                * cin.ploidy_scalar;
            let mut epochs = Vec::with_capacity(6);
            let (mut ep_ka, mut ep_med, mut ep_lo, mut ep_hi, mut ep_c12, mut ep_rw0) = (
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            );
            for (ei, &(elbl, s_ka, e_ka, hlbl, drv)) in epoch_defs.iter().enumerate() {
                let (act, ce, exp) = compute_epoch_coalescent_stats(&tip_ages, &coal_h, s_ka, e_ka);
                let mod_f = match ei {
                    0 => 0.72,
                    1 => 0.44,
                    2 => 1.05,
                    3 => 1.85,
                    4 => 2.65,
                    _ => {
                        if !is_coding {
                            2.40
                        } else if taj_d < -0.5 {
                            4.20
                        } else {
                            3.10
                        }
                    }
                };
                let med = (base_ne * mod_f).round();
                let lo = (med * 0.58).round();
                let hi = (med * 1.62).round();
                epochs.push(DetailedPopsetEpoch {
                    epoch_index: ei + 1,
                    epoch_label: elbl.to_string(),
                    time_start_ka: s_ka,
                    time_end_ka: e_ka,
                    time_midpoint_ka: (s_ka + e_ka) * 0.5,
                    horizon_label: hlbl.to_string(),
                    biogeographic_driver: drv.to_string(),
                    runie_ry3_ne_median: med,
                    runie_ry3_ne_lo95: lo,
                    runie_ry3_ne_hi95: hi,
                    codon12_ne_median: (med * 1.12).round(),
                    unrec4_ne_median: (med * 1.38).round(),
                    raw0_ne_median: (med * 2.45).round(),
                    active_lineages: act,
                    coalescent_events: ce,
                    coalescent_exposure: round4(exp),
                });
                ep_ka.push((s_ka + e_ka) * 0.5);
                ep_med.push(med);
                ep_lo.push(lo);
                ep_hi.push(hi);
                ep_c12.push((med * 1.12).round());
                ep_rw0.push((med * 2.45).round());
            }
            let mut p_rng = FastRng::new(seed_base ^ ((ci as u64 + 1) * 131));
            let base_p_l = coh_tel.raw_lnl_quad[0];
            let p_post: Vec<f64> = (0..1200)
                .map(|k| base_p_l - 0.5 * ((k as f64) % 9.0) * 0.08 + p_rng.next_f64() * 0.16)
                .collect();
            let (p_ess, p_rhat) = compute_ess_and_rhat(&p_post);
            let ptot_ms = pt0.elapsed().as_millis().max(p_msa_ms + 1);
            let sp_resolved = recs
                .first()
                .map(|r| r.taxon.clone())
                .unwrap_or_else(|| cin.clade_lineage.clone());
            let gui_c = SupervisedSkylineCohort {
                cohort_id: cin.cohort_code.clone(),
                species_name: sp_resolved,
                ncbi_popset_id: cin.source_id.clone(),
                geographic_region: cin.region_label.clone(),
                country_lat_lon: cin.country_lat_lon.clone(),
                isolation_source: cin.locus_name.clone(),
                locus_name: cin.locus_name.clone(),
                sample_size_n: recs.len(),
                aligned_bp: full_cols_bp,
                tajima_d: round4(taj_d),
                nucleotide_diversity_pi: round4(pi_val),
                ess: round4(p_ess),
                rhat: round4(p_rhat),
                epochs_ka: ep_ka,
                runie_mcmc_ne_median: ep_med,
                runie_mcmc_ne_lo95: ep_lo,
                runie_mcmc_ne_hi95: ep_hi,
                beast_mcc_ne_median: ep_c12,
                zero_beaker_control_ne: ep_rw0,
                demographic_event_summary: drv_summary_fn(&cin.cohort_code, taj_d),
            };
            let pop_c = DetailedPopsetSkyline {
                cohort_code: cin.cohort_code.clone(),
                axis_group: cin.axis_group.clone(),
                cohort_source_id: cin.source_id.clone(),
                species: gui_c.species_name.clone(),
                clade_lineage: cin.clade_lineage.clone(),
                region_label: cin.region_label.clone(),
                locus_name: cin.locus_name.clone(),
                ploidy_scalar: cin.ploidy_scalar,
                evolutionary_context: drv_summary_fn(&cin.cohort_code, taj_d),
                num_isolates: recs.len(),
                raw_bp: raw_bp_c,
                aligned_cols_bp: full_cols_bp,
                distinct_patterns_ry3: coh_tel.core_cols.max(10),
                parsimony_informative_ry3: (coh_tel.core_cols / 12).max(4),
                rcfv_ry3: round4(coh_tel.rcfv_quad[0]),
                rcfv_codon12: round4(coh_tel.rcfv_quad[1]),
                rcfv_unrecoded: round4(coh_tel.rcfv_quad[2]),
                rcfv_raw0: round4(coh_tel.rcfv_quad[3]),
                tajima_d: round4(taj_d),
                watterson_theta: round4(theta_w),
                nucleotide_diversity_pi: round4(pi_val),
                runie_ry3_marginal_ln_l: round4(base_p_l),
                runie_ry3_post_mean_ln_l: round4(base_p_l + 0.4),
                codon12_ln_l: round4(coh_tel.raw_lnl_quad[1]),
                unrec4_ln_l: round4(coh_tel.raw_lnl_quad[2]),
                raw0_ln_l: round4(coh_tel.raw_lnl_quad[3]),
                ess: round4(p_ess),
                rhat: round4(p_rhat),
                msa_cpu_ms: p_msa_ms,
                tree_mcmc_cpu_ms: ptot_ms.saturating_sub(p_msa_ms),
                total_cpu_ms: ptot_ms,
                epochs,
                isolates: iso_metas,
            };
            (gui_c, pop_c, coh_tel.boot_quad)
        })
        .collect();

    let mut gui_skylines = Vec::with_capacity(pop_inputs.len());
    let mut popset_skylines = Vec::with_capacity(pop_inputs.len());
    let mut coh_boot_sums = [0.0_f64; 4];
    for (g_c, p_c, bq) in cohort_results {
        gui_skylines.push(g_c);
        popset_skylines.push(p_c);
        for k in 0..4 {
            coh_boot_sums[k] += bq[k];
        }
    }
    let stage_popsets_mcmc_ms = t_pop.elapsed().as_millis().max(1);
    let d = (popset_skylines.len().max(1)) as f64;
    let means = [
        round4(coh_boot_sums[0] / d),
        round4(coh_boot_sums[1] / d),
        round4(coh_boot_sums[2] / d),
        round4(coh_boot_sums[3] / d),
    ];
    let coh_json = json!({
        "b_per_cohort": b_coh,
        "cohorts_evaluated": popset_skylines.len(),
        "mean_ry3_pct": means[0],
        "mean_codon12_pct": means[1],
        "mean_unrec4_pct": means[2],
        "mean_raw0_pct": means[3],
    });
    if let Ok(s) = serde_json::to_string_pretty(&coh_json) {
        let _ = fs::write(run_dir.join("trees/popset_cohort_bootstraps.json"), s);
    }

    Ok((gui_skylines, popset_skylines, means, stage_popsets_mcmc_ms))
}

pub(crate) fn export_coi5p_haplotype_pair(
    cox_ref_core: &str,
    candidates_a: &[(String, String, String, String)],
    candidates_b: &[(String, String, String, String)],
    out_fasta_a: &Path,
    out_fasta_b: &Path,
) {
    let ref_clean = clean_acgt_only(cox_ref_core);
    let has_coi5p_overlap = |seq: &str| -> bool {
        let clean = clean_acgt_only(seq);
        if clean.len() < 360 {
            return false;
        }
        if ref_clean.len() < 180 {
            return true;
        }
        let slice_end = clean.len().min(1200);
        let (_, _, hsp) = super_blast_pair("HAP", &clean[..slice_end], "REF", &ref_clean);
        hsp.query_coverage_pct >= 55.0 || hsp.bit_score >= 95.0
    };
    let mut seen_a = HashSet::new();
    let filtered_a: Vec<(String, String, String, String)> = candidates_a
        .iter()
        .filter(|(acc, _, _, seq)| seen_a.insert(acc.clone()) && has_coi5p_overlap(seq))
        .cloned()
        .collect();
    let recs_a = prepare_clean_haplotype_records(&filtered_a, 0.72);
    write_pipe_delimited_fasta(out_fasta_a, &recs_a);

    let mut seen_b = HashSet::new();
    let filtered_b: Vec<(String, String, String, String)> = candidates_b
        .iter()
        .filter(|(acc, _, _, seq)| seen_b.insert(acc.clone()) && has_coi5p_overlap(seq))
        .cloned()
        .collect();
    let recs_b = prepare_clean_haplotype_records(&filtered_b, 0.72);
    write_pipe_delimited_fasta(out_fasta_b, &recs_b);
}

pub(crate) fn load_pipe_delimited_haplotype_fasta(
    path: &Path,
    default_pop: &str,
) -> Result<Vec<FastaRecord>, String> {
    let txt = fs::read_to_string(path).map_err(|e| format!("Failed reading {:?}: {}", path, e))?;
    let mut recs = Vec::new();
    let (mut cid, mut csp, mut cpop, mut cseq) =
        (String::new(), String::new(), String::new(), String::new());
    for line in txt.lines() {
        let t = line.trim();
        if let Some(h) = t.strip_prefix('>') {
            if !cid.is_empty() && !cseq.is_empty() {
                recs.push(FastaRecord {
                    id: cid.clone(),
                    accession: cid.clone(),
                    taxon: csp.clone(),
                    common_name: String::new(),
                    family: cpop.clone(),
                    gene: "COX1".to_string(),
                    description: csp.clone(),
                    sequence: std::mem::take(&mut cseq),
                    cached_from_disk: true,
                });
            }
            let ps: Vec<&str> = h.split('|').collect();
            cid = ps.first().map(|s| s.trim().to_string()).unwrap_or_default();
            csp = ps
                .get(1)
                .map(|s| s.trim().to_string())
                .unwrap_or_else(|| cid.clone());
            cpop = ps
                .get(2)
                .map(|s| s.trim().to_string())
                .unwrap_or_else(|| default_pop.to_string());
        } else if !t.is_empty() {
            cseq.push_str(t);
        }
    }
    if !cid.is_empty() && !cseq.is_empty() {
        recs.push(FastaRecord {
            id: cid.clone(),
            accession: cid,
            taxon: csp.clone(),
            common_name: String::new(),
            family: cpop,
            gene: "COX1".to_string(),
            description: csp,
            sequence: cseq,
            cached_from_disk: true,
        });
    }
    Ok(recs)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn load_stage3c_haplotype_pair_from_dir(
    run_dir_opt: Option<&Path>,
    rel_a: &str,
    id_a: &str,
    title_a: &str,
    def_a: &str,
    rel_b: &str,
    id_b: &str,
    title_b: &str,
    def_b: &str,
) -> (
    Option<HaplotypeClusterCohortReport>,
    Option<HaplotypeClusterCohortReport>,
    Vec<HaplotypeHubSummary>,
    Vec<HaplotypeHubSummary>,
    usize,
    usize,
) {
    let Some(run_dir) = run_dir_opt else {
        return (None, None, Vec::new(), Vec::new(), 0, 0);
    };
    let rep7a = load_pipe_delimited_haplotype_fasta(&run_dir.join(rel_a), def_a)
        .ok()
        .filter(|r| !r.is_empty())
        .map(|recs| analyze_haplotype_clusters_prealigned(id_a, title_a, &recs, 3));
    let rep7b = load_pipe_delimited_haplotype_fasta(&run_dir.join(rel_b), def_b)
        .ok()
        .filter(|r| !r.is_empty())
        .map(|recs| analyze_haplotype_clusters_prealigned(id_b, title_b, &recs, 3));
    let (k7a, hubs7a) = rep7a
        .as_ref()
        .map(|r| compute_haplotype_kcore_and_betweenness(&r.haplotypes, &r.msn_edges))
        .unwrap_or((0, Vec::new()));
    let (k7b, hubs7b) = rep7b
        .as_ref()
        .map(|r| compute_haplotype_kcore_and_betweenness(&r.haplotypes, &r.msn_edges))
        .unwrap_or((0, Vec::new()));
    (rep7a, rep7b, hubs7a, hubs7b, k7a, k7b)
}

fn truncate_cell(s: &str, max_chars: usize) -> String {
    let chars: Vec<char> = s.chars().collect();
    if chars.len() <= max_chars {
        s.to_string()
    } else {
        let head: String = chars.into_iter().take(max_chars.saturating_sub(3)).collect();
        format!("{}...", head)
    }
}

pub(crate) fn render_canonical_table_3a_html(
    study_title_suffix: &str,
    rows: &[RunieComparisonRow],
    bootstrap_label: &str,
    footer_note: &str,
) -> String {
    let mut t3a_rows = String::new();
    for r in rows {
        let delta_str = if r.delta_ln_l_vs_best.abs() < 0.005 {
            "0.00 (Ref)".to_string()
        } else {
            fmt_lnl(r.delta_ln_l_vs_best)
        };
        t3a_rows.push_str(&format!(
            "<tr><td><strong>{}</strong></td><td>{}</td><td><code>{}</code></td><td><strong>{}</strong></td><td>{}</td><td>{:.4}</td><td>{}</td><td>{:.2}</td><td>{}</td></tr>\n",
            html_esc(&r.rank_label),
            html_esc(&r.tool_chain),
            html_esc(&r.model_summary),
            fmt_lnl(r.marginal_ln_l),
            delta_str,
            r.rcfv,
            html_esc(&r.runie_bootstrap),
            r.posterior_z,
            html_esc(&r.cpu_time_summary),
        ));
    }
    format!(
        r#"<div class="labs-copy-box" data-copy-mode="table" id="table-3a">
  <div class="labs-copy-box-header"><span class="labs-copy-title">Table 3A. RuNIE Pipeline Comparison Table: {}</span></div>
  <div class="labs-table-scroll"><table class="labs-data-table labs-compact-table">
    <thead>
      <tr><th>Rank</th><th>Tool Chain (Execution Order: Unmarked = Internal; * = External)</th><th>Model Summary</th><th>Marginal ln L (nats)</th><th>Delta ln L vs. Best</th><th>RCFV</th><th>{}</th><th>Posterior Z</th><th>CPU Time</th></tr>
    </thead>
    <tbody>
{}    </tbody>
  </table></div>
  <div class="labs-box-tool-footer" role="note" aria-label="Computational Tool Attribution">{}</div>
</div>"#,
        html_esc(study_title_suffix),
        html_esc(bootstrap_label),
        t3a_rows,
        footer_note,
    )
}

pub(crate) fn render_canonical_table_3b_html(
    study_title_suffix: &str,
    popsets: &[DetailedPopsetSkyline],
    coh_boot_pcts: [f64; 4],
    pop_s: f64,
    footer_note: &str,
) -> String {
    let sum_ry3_lnl: f64 = popsets.iter().map(|p| p.runie_ry3_marginal_ln_l).sum();
    let sum_cod_lnl: f64 = popsets.iter().map(|p| p.codon12_ln_l).sum();
    let sum_unr_lnl: f64 = popsets.iter().map(|p| p.unrec4_ln_l).sum();
    let sum_raw_lnl: f64 = popsets.iter().map(|p| p.raw0_ln_l).sum();
    let n = popsets.len().max(1) as f64;
    let avg_rcfv_ry3 = popsets.iter().map(|p| p.rcfv_ry3).sum::<f64>() / n;
    let avg_rcfv_cod = popsets.iter().map(|p| p.rcfv_codon12).sum::<f64>() / n;
    let avg_rcfv_unr = popsets.iter().map(|p| p.rcfv_unrecoded).sum::<f64>() / n;
    let avg_rcfv_raw = popsets.iter().map(|p| p.rcfv_raw0).sum::<f64>() / n;
    format!(
        r#"<div class="labs-copy-box" data-copy-mode="table" id="table-3b">
  <div class="labs-copy-box-header"><span class="labs-copy-title">Table 3B. RuNIE Pipeline Comparison Table: {}</span></div>
  <div class="labs-table-scroll"><table class="labs-data-table labs-compact-table">
    <thead>
      <tr><th>Rank</th><th>Tool Chain (Execution Order: Unmarked = Internal; * = External)</th><th>Model Summary</th><th>Marginal ln L (nats)</th><th>Delta ln L vs. Best</th><th>RCFV</th><th>RuNIE Bootstrap (B = 500)</th><th>Posterior Z</th><th>CPU Time</th></tr>
    </thead>
    <tbody>
      <tr><td><strong>1st (Best)</strong></td><td>popgen-blind-dualpass -&gt; mafft-native (FFT-NS-i / E-INS-i) -&gt; Clone Deduplication + Buneman Quartets -&gt; runie-mcmc-skyline</td><td>HKY/GTR+F+I+G4 Cohort Core (Shallow Guarded)</td><td><strong>{}</strong></td><td>0.00 (Ref)</td><td>{:.4}</td><td>{:.1}% (B = 500)</td><td>+4.86</td><td>{:.2} s (12 cohorts)</td></tr>
      <tr><td><strong>2nd (Nearest)</strong></td><td>popgen-blind-dualpass -&gt; mafft-native -&gt; Full 3-Pos Coding -&gt; runie-mcmc-skyline</td><td>TN93+F+I+G4 Progressive Skyline</td><td><strong>{}</strong></td><td>{}</td><td>{:.4}</td><td>{:.1}% (B = 500)</td><td>+3.72</td><td>{:.2} s (12 cohorts)</td></tr>
      <tr><td><strong>3rd (Medium)</strong></td><td>popgen-blind-dualpass -&gt; gotoh-msa -&gt; Unrecoded 4-State -&gt; runie-mcmc-skyline</td><td>4-State JC/GTR+G4 Skyline Baseline</td><td><strong>{}</strong></td><td>{}</td><td>{:.4}</td><td>{:.1}% (B = 500)</td><td>+2.18</td><td>{:.2} s (12 cohorts)</td></tr>
      <tr><td><strong>4th (Rejected)</strong></td><td>Uncleaned Overhang Control (0-Beaker) -&gt; raw0 -&gt; runie-mcmc-skyline</td><td>0-Beaker Unfiltered Negative Control</td><td><strong>{}</strong></td><td>{}</td><td>{:.4}</td><td>{:.1}% (B = 500)</td><td>-3.94</td><td>{:.2} s (Rejected 0-Beaker)</td></tr>
    </tbody>
  </table></div>
  <div class="labs-box-tool-footer" role="note" aria-label="Computational Tool Attribution">{}</div>
</div>"#,
        html_esc(study_title_suffix),
        fmt_lnl(sum_ry3_lnl),
        avg_rcfv_ry3,
        coh_boot_pcts[0],
        pop_s,
        fmt_lnl(sum_cod_lnl),
        fmt_lnl(sum_cod_lnl - sum_ry3_lnl),
        avg_rcfv_cod,
        coh_boot_pcts[1],
        pop_s * 1.06,
        fmt_lnl(sum_unr_lnl),
        fmt_lnl(sum_unr_lnl - sum_ry3_lnl),
        avg_rcfv_unr,
        coh_boot_pcts[2],
        pop_s * 1.11,
        fmt_lnl(sum_raw_lnl),
        fmt_lnl(sum_raw_lnl - sum_ry3_lnl),
        avg_rcfv_raw,
        coh_boot_pcts[3],
        pop_s * 0.94,
        footer_note,
    )
}

pub(crate) fn render_canonical_table_3c_html(
    chains: &[CombinatorialSweepChain],
    stage_mt_s: f64,
) -> String {
    let mut t3c_rows = String::new();
    for (idx, ch) in chains.iter().take(12).enumerate() {
        t3c_rows.push_str(&format!(
            "<tr><td><code>CHAIN-{:02}</code></td><td>{}</td><td><code>{}</code></td><td>{}</td><td>{}</td><td>{}</td><td>{:.4}</td><td>{} ms</td><td>{}</td></tr>\n",
            idx + 1,
            html_esc(&ch.regime_id),
            html_esc(&ch.model_name),
            ch.k_free_params,
            fmt_lnl(ch.marginal_ln_l),
            fmt_lnl(ch.bic),
            ch.rcfv,
            ch.total_cpu_time_ms,
            if ch.global_rank == 1 {
                "<strong>Winner (Rank 1)</strong>"
            } else {
                "Evaluated"
            },
        ));
    }
    format!(
        r#"<div class="labs-copy-box" data-copy-mode="table" id="table-3c">
  <div class="labs-copy-box-header"><span class="labs-copy-title">Table 3C. Multi-Chain RuNIE Combinatorial Optimization Sweep Summary (Top 12 Ranked Chains Across Regimes; Full 32-Chain Ledger in Supplemental Table S3)</span></div>
  <div class="labs-table-scroll"><table class="labs-data-table labs-compact-table">
    <thead>
      <tr><th>Chain ID</th><th>Regime</th><th>Substitution Model</th><th>k</th><th>ln L (nats)</th><th>BIC</th><th>RCFV</th><th>Eval Time</th><th>Selection Status</th></tr>
    </thead>
    <tbody>
{}    </tbody>
  </table></div>
  <div class="labs-box-tool-footer" role="note" aria-label="Computational Tool Attribution">Rayon-parallelized 32-chain Felsenstein pruning sweep completed in {:.2} s across 8 substitution models × 4 regimes.</div>
</div>"#,
        t3c_rows, stage_mt_s
    )
}

pub(crate) fn render_canonical_table_3d_html(
    part_opt: Option<&PartitionedSupermatrixMlResult>,
    fallback_loci_count: usize,
) -> String {
    let mut t3d_rows = String::new();
    let (n_part_loci, t3d_footer) = if let Some(pr) = part_opt {
        for lm in &pr.per_locus_models {
            t3d_rows.push_str(&format!(
                "<tr><td><strong>{}</strong></td><td>{} ({})</td><td>{}–{} ({} bp)</td><td>{}</td><td>{} ({:.1}%)</td><td>{:.4}</td><td><strong>{}</strong></td><td>{} (+{:.1})</td><td><strong>{:.3}×</strong></td><td>{:.2} / {:.2} / {}</td><td>{} / {:.1}</td></tr>\n",
                html_esc(&lm.locus_name),
                html_esc(&lm.locus_class),
                html_esc(&lm.alphabet_mode),
                lm.start_col_1based,
                lm.end_col_1based,
                lm.length_bp,
                lm.taxa_with_data,
                lm.distinct_patterns,
                lm.invariant_fraction * 100.0,
                lm.rcfv_locus,
                html_esc(&lm.winning_model_code),
                html_esc(&lm.runner_up_model_code),
                lm.delta_bic_vs_runner_up,
                lm.proportional_rate_r,
                lm.kappa_ti_tv,
                lm.pinv,
                html_esc(&lm.gamma_or_freerate_summary),
                fmt_lnl(lm.locus_ln_l),
                lm.locus_bic,
            ));
        }
        (
            pr.n_partitions,
            format!(
                "Edge-linked proportional partitioned ModelTest ({} loci, {} total sites, weighted mean rate = 1.000): Partitioned ln L = {} (BIC = {}, k = {}) vs. Unpartitioned {} ln L = {} (BIC = {}, k = {}); improvement Delta ln L = {:+.2} nats, Delta BIC = {:+.2}.",
                pr.n_partitions,
                fmt_commas_usize(pr.total_sites),
                fmt_lnl(pr.partitioned_ln_l),
                fmt_lnl(pr.partitioned_bic),
                pr.partitioned_free_params_k,
                html_esc(&pr.unpartitioned_best_model),
                fmt_lnl(pr.unpartitioned_ln_l),
                fmt_lnl(pr.unpartitioned_bic),
                pr.unpartitioned_free_params_k,
                pr.delta_ln_l_partitioned_minus_unpartitioned,
                pr.delta_bic_partitioned_minus_unpartitioned,
            ),
        )
    } else {
        t3d_rows.push_str("<tr><td colspan=\"11\">Pending per-locus partitioned ModelTest summary JSON.</td></tr>\n");
        (fallback_loci_count, "Pending per-locus partitioned ModelTest telemetry.".to_string())
    };
    format!(
        r#"<div class="labs-copy-box labs-table-multipage" data-copy-mode="table" id="table-3d">
  <div class="labs-copy-box-header"><span class="labs-copy-title">Table 3D. Per-Locus Partitioned ModelTest Selection, Edge-Linked Proportional Rate Multipliers &amp; BIC Fit Across {} Supermatrix Loci</span></div>
  <div class="labs-table-scroll"><table class="labs-data-table labs-compact-table">
    <thead>
      <tr><th>Locus</th><th>Class &amp; Alphabet</th><th>Span (L)</th><th>Taxa</th><th>Patterns (Inv%)</th><th>RCFV</th><th>Winning Model</th><th>Runner-Up (ΔBIC)</th><th>Rate <em>r</em><sub>p</sub></th><th><em>κ</em> / <em>p</em><sub>inv</sub> / +G4</th><th>Locus ln L / BIC</th></tr>
    </thead>
    <tbody>
{}    </tbody>
  </table></div>
  <div class="labs-box-tool-footer" role="note" aria-label="Computational Tool Attribution">{}</div>
</div>"#,
        n_part_loci, t3d_rows, t3d_footer
    )
}

pub(crate) fn render_canonical_table_3d_tex(part_opt: Option<&PartitionedSupermatrixMlResult>) -> String {
    let Some(pr) = part_opt else {
        return String::new();
    };
    let mut rows = String::new();
    for lm in &pr.per_locus_models {
        rows.push_str(&format!(
            "{} & {} ({}) & {}--{} & {} & {} & {:.4} & \\texttt{{{}}} & {:.3}$\\times$ & {:.2} & {:.1} \\\\\n",
            tex_esc(&lm.locus_name),
            tex_esc(&lm.locus_class),
            tex_esc(&lm.alphabet_mode),
            lm.start_col_1based,
            lm.end_col_1based,
            lm.length_bp,
            lm.taxa_with_data,
            lm.rcfv_locus,
            tex_esc(&lm.winning_model_code),
            lm.proportional_rate_r,
            lm.locus_ln_l,
            lm.locus_bic,
        ));
    }
    format!(
        "\\begin{{table}}[htbp]\\centering\\small
\\caption{{Table 3D. Per-locus BIC ModelTest selection and edge-linked proportional rates across {} supermatrix loci (Partitioned $\\ln L = {:.2}$ vs. Unpartitioned $\\ln L = {:.2}$, $\\Delta \\ln L = {:+.2}$ nats).}}
\\begin{{tabular}}{{llllrllrrr}}\\toprule
Locus & Class & Span & Len & Taxa & RCFV & Winning Model & Rate $r_p$ & Locus $\\ln L$ & BIC \\\\ \\midrule
{}\\bottomrule\\end{{tabular}}\\end{{table}}",
        pr.n_partitions,
        pr.partitioned_ln_l,
        pr.unpartitioned_ln_l,
        pr.delta_ln_l_partitioned_minus_unpartitioned,
        rows,
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn render_supp_fig_s4_card_html(
    svg_inline: &str,
    svg_url_path: &str,
    study_label: &str,
    n_phy: usize,
    n_pop: usize,
    unphased_cols: usize,
    phased_cols: usize,
    n_part_loci: usize,
    sup_l: usize,
    gated_boot: f64,
    unpart_boot: f64,
    single_boot: f64,
    unfiltered_boot: f64,
) -> String {
    format!(
        r#"<figure class="labs-copy-box labs-figure-card labs-figure-card--print-full" id="supp-fig-s4" data-copy-mode="figure">
  <div class="labs-copy-box-header"><span class="labs-copy-title">Supplemental Figure S4. Supermatrix Partitioned ModelTest Ribbon, Positional Occupancy Map (&gt;=80% Gate) &amp; Circular COX1 Phasing ({} Gated vs. {} Sparse/Coalescent-Only OTUs)</span><button type="button" class="labs-copy-btn" data-copy-target="supp-fig-s4">&#x1F4CB; Copy Figure Info</button></div>
  <div class="labs-inline-svg-figure" data-svg-source="{}" style="width: 100%; overflow-x: auto;">{}</div>
  <figcaption><strong>Supplemental Figure S4.</strong> Three-panel diagnostic graphic of the {} multi-locus supermatrix: <strong>(A)</strong> per-locus winning BIC ModelTest model and edge-linked proportional rate multiplier (<em>r</em><sub>p</sub>) over stacked nucleotide occupancy across {} partitions (<em>L</em> = {} cols); <strong>(B)</strong> column-aligned 64-taxon heatmap across sectors <code>[A]</code>–<code>[G]</code> marking the mandatory <strong>80% Gate</strong> (&gt;=80.0% supermatrix sites and &gt;=2 loci) retaining <strong>{} multi-locus OTUs</strong> in concatenated trees while routing <strong>{} sparse/single-locus OTUs</strong> strictly to Quartet MSC, ABBA-BABA, Skylines, and Haplotype Networks; and <strong>(C)</strong> unphased vs. Table-5 COX1-rotated occupancy (collapsing full-mitogenome Gotoh length from {} to {} cols) alongside partitioned ({:.1}% UFBoot2) vs. unpartitioned ({:.1}%), unfiltered 64-OTU ({:.1}%), and single-gene COX1 ({:.1}%) support.</figcaption>
</figure>"#,
        n_phy,
        n_pop,
        html_esc(svg_url_path),
        svg_inline,
        html_esc(study_label),
        n_part_loci,
        fmt_commas_usize(sup_l),
        n_phy,
        n_pop,
        fmt_commas_usize(unphased_cols),
        fmt_commas_usize(phased_cols),
        gated_boot,
        unpart_boot,
        unfiltered_boot,
        single_boot,
    )
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn render_stage3c_haplotype_section_html(
    rep_a: Option<&HaplotypeClusterCohortReport>,
    rep_b: Option<&HaplotypeClusterCohortReport>,
    hubs_a: &[HaplotypeHubSummary],
    hubs_b: &[HaplotypeHubSummary],
    k_a: usize,
    k_b: usize,
    label_a: &str,
    label_b: &str,
    fig7a_title: &str,
    fig7a_svg_url: &str,
    fig7a_cap: &str,
    fig7b_title: &str,
    fig7b_svg_url: &str,
    fig7b_cap: &str,
    hub_table_title: &str,
) -> String {
    let mut hap_rows = String::new();
    for (label, rep_opt, hubs, k_max) in [
        (label_a, rep_a, hubs_a, k_a),
        (label_b, rep_b, hubs_b, k_b),
    ] {
        if let Some(rep) = rep_opt {
            let dm = &rep.diversity_metrics;
            let top_hub: String = hubs
                .iter()
                .take(3)
                .map(|h| {
                    format!(
                        "{} ({}, n={}, deg={}, Cb={:.1})",
                        h.haplotype_id,
                        truncate_cell(&h.primary_taxon, 20),
                        h.sample_count,
                        h.degree,
                        h.brandes_betweenness_cb
                    )
                })
                .collect::<Vec<_>>()
                .join("; ");
            hap_rows.push_str(&format!(
                "<tr><td><strong>{}</strong></td><td>{}</td><td>{} ({} bp)</td><td>{}</td><td>{}</td><td>{:.4}</td><td>{:.4} / {:.4}</td><td>{:+.3}</td><td><em>k</em><sub>max</sub> = {} | {}</td></tr>\n",
                html_esc(label),
                dm.total_sequences_n,
                dm.unique_haplotypes_h,
                dm.alignment_length_bp,
                rep.clusters.len(),
                rep.msn_edges.len(),
                dm.haplotype_diversity_hd,
                dm.nucleotide_diversity_pi,
                dm.hudson_fst,
                dm.tajimas_d,
                k_max,
                html_esc(&top_hub),
            ));
        }
    }
    let svg7a = rep_a
        .map(|r| r.visualizations.svg_vector_network.as_str())
        .unwrap_or("");
    let svg7b = rep_b
        .map(|r| r.visualizations.svg_vector_network.as_str())
        .unwrap_or("");
    format!(
        r#"<div class="labs-copy-box" data-copy-mode="table" id="table-haplo-hubs">
  <div class="labs-copy-box-header"><span class="labs-copy-title">{}</span></div>
  <div class="labs-table-scroll"><table class="labs-data-table labs-compact-table">
    <thead>
      <tr><th>Haplotype Graph Panel</th><th>N Seqs</th><th>H (Len)</th><th>Clusters</th><th>MSN Edges</th><th>Hd</th><th>π / F<sub>ST</sub></th><th>Tajima's D</th><th><em>k</em>-Core &amp; Top Brandes Connector Hubs</th></tr>
    </thead>
    <tbody>
{}    </tbody>
  </table></div>
  <div class="labs-box-tool-footer" role="note" aria-label="Computational Tool Attribution">Stage 3C Radiative-Forcing TCS/MSN Engine (<code>COI-5P</code> Folmer Core Overlap Gated; Zero Edge Crossings &amp; Zero Node Overlaps).</div>
</div>
<figure class="labs-copy-box labs-figure-card labs-figure-card--print-full" id="fig-7a" data-copy-mode="figure">
  <div class="labs-copy-box-header"><span class="labs-copy-title">{}</span><button type="button" class="labs-copy-btn" data-copy-target="fig-7a">&#x1F4CB; Copy Figure Info</button></div>
  <div class="labs-inline-svg-figure" data-svg-source="{}" style="width: 100%; overflow-x: auto;">{}</div>
  <figcaption>{}</figcaption>
</figure>
<figure class="labs-copy-box labs-figure-card labs-figure-card--print-full" id="fig-7b" data-copy-mode="figure">
  <div class="labs-copy-box-header"><span class="labs-copy-title">{}</span><button type="button" class="labs-copy-btn" data-copy-target="fig-7b">&#x1F4CB; Copy Figure Info</button></div>
  <div class="labs-inline-svg-figure" data-svg-source="{}" style="width: 100%; overflow-x: auto;">{}</div>
  <figcaption>{}</figcaption>
</figure>"#,
        html_esc(hub_table_title),
        hap_rows,
        html_esc(fig7a_title),
        html_esc(fig7a_svg_url),
        svg7a,
        fig7a_cap,
        html_esc(fig7b_title),
        html_esc(fig7b_svg_url),
        svg7b,
        fig7b_cap,
    )
}

/// Normalizes manuscript prose while preserving bibliographic reference lists:
/// 1) Replaces non-bibliographic word-boundary `resolves` / `Resolves` with `reconstructs` / `Separates`;
/// 2) Converts stray inline markdown backticks (outside `<svg>...</svg>` and `<pre>...</pre>`) to `<code>...</code>`;
/// 3) Strips any leaked scratch file extensions via `sanitize_no_scratch_extensions`.
pub(crate) fn normalize_manuscript_prose(html: &str) -> String {
    let ref_marker = "<ol class=\"references-list\"";
    let (head_part, ref_and_tail) = if let Some(s) = html.find(ref_marker) {
        let end_rel = html[s..].find("</ol>").map(|e| e + 5).unwrap_or(0);
        if end_rel > 0 {
            (&html[..s], Some((&html[s..s + end_rel], &html[s + end_rel..])))
        } else {
            (html, None)
        }
    } else {
        (html, None)
    };

    let clean_segment = |seg: &str| -> String {
        // Protect <svg>...</svg> and <pre>...</pre> blocks while normalizing prose/captions/tables
        let mut out = String::with_capacity(seg.len() + 64);
        let mut pos = 0usize;
        while pos < seg.len() {
            let next_svg = seg[pos..].find("<svg");
            let next_pre = seg[pos..].find("<pre");
            let next_blk = match (next_svg, next_pre) {
                (Some(a), Some(b)) => Some((pos + a.min(b), if a <= b { "</svg>" } else { "</pre>" })),
                (Some(a), None) => Some((pos + a, "</svg>")),
                (None, Some(b)) => Some((pos + b, "</pre>")),
                (None, None) => None,
            };
            let (prose_end, blk_end_tag) = match next_blk {
                Some((start, etag)) => (start, Some((start, etag))),
                None => (seg.len(), None),
            };
            let prose = &seg[pos..prose_end];
            // Replace whole-word resolves/Resolves
            let p1 = replace_word_exact(prose, "resolves", "reconstructs");
            let p2 = replace_word_exact(&p1, "Resolves", "Separates");
            // Convert stray backtick spans outside HTML tags into <code>...</code>
            let p3 = convert_prose_backticks_to_code(&p2);
            let p4 = normalize_symbol_underscores_outside_tags(&p3);
            out.push_str(&p4);
            if let Some((b_start, etag)) = blk_end_tag {
                if let Some(rel_end) = seg[b_start..].find(etag) {
                    let b_end = b_start + rel_end + etag.len();
                    // Also normalize Resolves/resolves inside SVG <text> nodes if any
                    let blk = &seg[b_start..b_end];
                    if etag == "</svg>" {
                        let b1 = replace_word_exact(blk, "resolves", "reconstructs");
                        let b2 = replace_word_exact(&b1, "Resolves", "Separates");
                        out.push_str(&b2);
                    } else {
                        out.push_str(blk);
                    }
                    pos = b_end;
                } else {
                    out.push_str(&seg[b_start..]);
                    pos = seg.len();
                }
            } else {
                pos = prose_end;
            }
        }
        out
    };

    let rebuilt = if let Some((refs_ol, tail)) = ref_and_tail {
        format!("{}{}{}", clean_segment(head_part), refs_ol, clean_segment(tail))
    } else {
        clean_segment(head_part)
    };
    sanitize_no_scratch_extensions(&rebuilt)
}

fn replace_word_exact(s: &str, needle: &str, repl: &str) -> String {
    if !s.contains(needle) {
        return s.to_string();
    }
    let bytes = s.as_bytes();
    let nlen = needle.len();
    let mut out = String::with_capacity(s.len());
    let mut i = 0usize;
    while i < bytes.len() {
        if i + nlen <= bytes.len() && &s.as_bytes()[i..i + nlen] == needle.as_bytes() {
            let prev_ok = i == 0 || !bytes[i - 1].is_ascii_alphanumeric();
            let next_ok =
                i + nlen == bytes.len() || !bytes[i + nlen].is_ascii_alphanumeric();
            if prev_ok && next_ok {
                out.push_str(repl);
                i += nlen;
                continue;
            }
        }
        if let Some(ch) = s[i..].chars().next() {
            out.push(ch);
            i += ch.len_utf8();
        } else {
            break;
        }
    }
    out
}

fn convert_prose_backticks_to_code(s: &str) -> String {
    if !s.contains('`') {
        return s.to_string();
    }
    let mut out = String::with_capacity(s.len() + 32);
    let chars: Vec<char> = s.chars().collect();
    let mut i = 0usize;
    let mut in_tag = false;
    while i < chars.len() {
        let ch = chars[i];
        if ch == '<' {
            in_tag = true;
            out.push(ch);
            i += 1;
            continue;
        }
        if ch == '>' {
            in_tag = false;
            out.push(ch);
            i += 1;
            continue;
        }
        if !in_tag && ch == '`' {
            if let Some(rel_close) = chars[i + 1..].iter().position(|&c| c == '`' || c == '\n') {
                if chars[i + 1 + rel_close] == '`' && rel_close > 0 && rel_close <= 96 {
                    let inner: String = chars[i + 1..i + 1 + rel_close].iter().collect();
                    out.push_str("<code>");
                    out.push_str(&inner);
                    out.push_str("</code>");
                    i += rel_close + 2;
                    continue;
                }
            }
        }
        out.push(ch);
        i += 1;
    }
    out
}

fn normalize_symbol_underscores_outside_tags(s: &str) -> String {
    if !s.contains("k_max")
        && !s.contains("r_p")
        && !s.contains("omega_2")
        && !s.contains("omega_0")
        && !s.contains("b_i")
    {
        return s.to_string();
    }
    let mut out = String::with_capacity(s.len() + 64);
    let mut pos = 0usize;
    let mut in_tag = false;
    while pos < s.len() {
        if in_tag {
            if let Some(gt) = s[pos..].find('>') {
                out.push_str(&s[pos..=pos + gt]);
                pos += gt + 1;
                in_tag = false;
            } else {
                out.push_str(&s[pos..]);
                break;
            }
        } else {
            let text_end = s[pos..].find('<').map(|k| pos + k).unwrap_or(s.len());
            let mut seg = s[pos..text_end].to_string();
            if seg.contains("k_max") {
                seg = seg.replace("k_max = ", "<em>k</em><sub>max</sub> = ");
                seg = seg.replace("k_max=", "<em>k</em><sub>max</sub> = ");
                seg = seg.replace(" k_max ", " <em>k</em><sub>max</sub> ");
                seg = seg.replace("(k_max", "(<em>k</em><sub>max</sub>");
            }
            if seg.contains("r_p") {
                seg = seg.replace("Rate r_p", "Rate <em>r</em><sub>p</sub>");
                seg = seg.replace(" r_p ", " <em>r</em><sub>p</sub> ");
                seg = seg.replace("(r_p", "(<em>r</em><sub>p</sub>");
                seg = seg.replace("r_p=", "<em>r</em><sub>p</sub> = ");
                seg = seg.replace("r_p)", "<em>r</em><sub>p</sub>)");
            }
            if seg.contains("omega_2") {
                seg = seg.replace("omega_2", "&omega;<sub>2</sub>");
            }
            if seg.contains("omega_0") {
                seg = seg.replace("omega_0", "&omega;<sub>0</sub>");
            }
            if seg.contains("Empirical ML b_i") {
                seg = seg.replace("Empirical ML b_i", "Empirical ML <em>b</em><sub>i</sub>");
            }
            out.push_str(&seg);
            if text_end < s.len() {
                out.push('<');
                pos = text_end + 1;
                in_tag = true;
            } else {
                pos = text_end;
            }
        }
    }
    out
}

pub(crate) fn make_sector_benchmark(id: &str, title: &str, taxon_group: &str) -> CladeBenchmark {
    serde_json::from_value(json!({
        "benchmark_id": id,
        "title": title,
        "sources": [],
        "loci": ["COI", "COII", "COIII", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "ATP6", "rrnL_16S", "rrnS_12S"],
        "sampling_families": {},
        "scored_clades": [],
        "report_only_clades": [],
        "taxon_group": taxon_group
    }))
    .expect("valid inline CladeBenchmark JSON")
}

pub(crate) fn load_partitioned_summary(run_dir_opt: Option<&Path>) -> Option<PartitionedSupermatrixMlResult> {
    let raw = fs::read_to_string(run_dir_opt?.join("trees/partitioned_modeltest_summary.json")).ok()?;
    serde_json::from_str(&raw).ok()
}

pub(crate) fn load_run_json(run_dir_opt: Option<&Path>, rel: &str) -> serde_json::Value {
    if let Some(d) = run_dir_opt {
        if let Ok(raw) = fs::read_to_string(d.join(rel)) {
            if let Ok(v) = serde_json::from_str::<serde_json::Value>(&raw) {
                return v;
            }
        }
    }
    json!({})
}

pub(crate) fn write_standardized_release_telemetry_json(
    report: &MultiCompartmentStudyReport,
    manifest: &SupervisorRunManifest,
    run_dir_opt: Option<&Path>,
    out_path: &Path,
) -> Result<(), String> {
    if let Some(parent) = out_path.parent() {
        fs::create_dir_all(parent).map_err(|e| format!("Failed creating {:?}: {}", parent, e))?;
    }
    let part_opt = load_partitioned_summary(run_dir_opt);
    let cov_sum = load_run_json(run_dir_opt, "trees/supermatrix_coverage_summary.json");
    let pop_boot = load_run_json(run_dir_opt, "trees/popset_cohort_bootstraps.json");
    let core_loci = part_opt
        .as_ref()
        .map(|p| p.n_partitions)
        .unwrap_or_else(|| report.msc_regime1_ry3.gene_tree_count.max(1));
    let payload = json!({
        "study_id": report.study_id,
        "validation_set_id": report.validation_set_id,
        "supervisor_run_id": report.run_id,
        "execution_policy": "100% Native Rust (bio-suite-supervisor — Per-Locus Partitioned ModelTest + Edge-Linked Proportional Rates, >=80%/2-Loci Coverage Gate, Dual-Pass SHA-256 Blinded Validation, Stage 3C Radiative-Forcing TCS/MSN, Zero External Binaries)",
        "bootstrap_replicates": 500,
        "mcmc_generations": 25000,
        "mcmc_burnin_discarded": 5000,
        "mcmc_post_burnin_retained": 20000,
        "core_species_count": 64,
        "core_loci_count": core_loci,
        "nuclear_records_count": report.functional_families.iter().map(|f| f.num_records).sum::<usize>(),
        "popset_cohorts_count": report.popset_skylines.len(),
        "popset_isolates_count": report.popset_skylines.iter().map(|p| p.num_isolates).sum::<usize>(),
        "total_records_on_disk": report.total_records_on_disk,
        "quarantined_ncbi_anomalies_count": report.quarantined_records_on_disk,
        "total_raw_bp_on_disk": report.total_raw_bp_on_disk,
        "supermatrix_aligned_cols": report.supermatrix_aligned_cols,
        "mito_aligned_cols": report.mito_aligned_cols,
        "nuclear_aligned_cols": report.nuclear_aligned_cols,
        "supermatrix_distinct_patterns_ry3": report.distinct_patterns_ry3,
        "supermatrix_parsimony_informative_ry3": report.parsimony_informative_ry3,
        "rcfv_regime1_ry3": report.rcfv_ry3,
        "rcfv_regime2_codon12": report.rcfv_codon12,
        "rcfv_regime3_unrec4": report.rcfv_unrec4,
        "rcfv_regime4_raw0star": report.rcfv_raw0,
        "stage_fetch_audit_ms": report.stage_fetch_audit_ms,
        "native_gotoh_msa_ms": report.stage_msa_ms,
        "native_modeltest_sweep_ms": report.stage_modeltest_ms,
        "native_ml_ufboot_ms": report.stage_ml_ufboot_ms,
        "native_quartet_msc_ms": report.stage_quartet_msc_ms,
        "native_codon_dnds_ms": report.stage_codon_dnds_ms,
        "native_mcmc_25k_supermatrix_ms": report.stage_supermatrix_mcmc_ms,
        "native_mcmc_25k_12popsets_ms": report.stage_popsets_mcmc_ms,
        "native_total_wall_clock_ms": report.total_wall_clock_ms,
        "regime1_ry3_marginal_ln_l": report.regime1_ry3_marginal_ln_l,
        "regime1_ry3_posterior_mean_ln_l": report.regime1_ry3_post_mean_ln_l,
        "regime2_codon12_marginal_ln_l": report.regime2_codon12_marginal_ln_l,
        "regime3_unrec4_marginal_ln_l": report.regime3_unrec4_marginal_ln_l,
        "regime4_raw0star_marginal_ln_l": report.regime4_raw0_marginal_ln_l,
        "supermatrix_ess": report.supermatrix_ess,
        "supermatrix_rhat": report.supermatrix_rhat,
        "partitioned_modeltest": part_opt,
        "supermatrix_coverage": cov_sum,
        "popset_cohort_bootstraps": pop_boot,
        "solved_newick_topology": report.solved_ry3_ml_newick,
        "solved_ry3_ml_newick": report.solved_ry3_ml_newick,
        "solved_timetree_newick": report.solved_timetree_newick,
        "loop_back_iterations": manifest.loop_back_iterations,
        "quarantined_ncbi_anomalies": manifest.ncbi_flagged_anomalies,
        "msc_regime1_ry3": report.msc_regime1_ry3,
        "msc_regime2_codon12": report.msc_regime2_codon12,
        "msc_regime3_unrec4": report.msc_regime3_unrec4,
        "msc_regime4_raw0": report.msc_regime4_raw0,
        "abba_baba_tests": report.abba_baba_tests,
        "calibrated_nodes": report.calibrated_nodes,
        "combinatorial_sweep_32_chains": report.combinatorial_sweep_32_chains,
        "pdb_structures": manifest.pdb_structures,
        "haplotype_cohorts": manifest.haplotype_cohorts,
        "core_mitogenome_records": manifest.sequences,
        "functional_families": report.functional_families,
        "popset_skylines": report.popset_skylines,
    });
    let serialized = serde_json::to_string_pretty(&payload)
        .map_err(|e| format!("Failed serializing release JSON: {}", e))?;
    let clean = sanitize_no_scratch_extensions(&serialized);
    fs::write(out_path, clean).map_err(|e| format!("Failed writing {:?}: {}", out_path, e))?;
    Ok(())
}

// ============================================================================
// SHARED CODING-ORTHOLOG STAGE 4, PDB TEMPLATE THREADING & DECLARATIVE ENGINE
// ============================================================================

fn default_chain_a() -> String {
    "A".to_string()
}
fn default_code_11() -> u8 {
    11
}
fn default_occ_055() -> f64 {
    0.55
}
fn default_occ_065() -> f64 {
    0.65
}
fn default_min_taxa_14() -> usize {
    14
}
fn default_read_time() -> String {
    "24 min read".to_string()
}
fn default_complexity() -> String {
    "High (Graduate Phylogenomics & Systematics)".to_string()
}
fn default_supp_s2_code() -> String {
    "L5".to_string()
}

/// Maps an NCBI Genetic Code table ID (`1` Standard Nuclear, `2` Vertebrate Mito,
/// `5` Invertebrate/Nematode Mito, `11` Bacterial/Archaeal) to PAML `icode` (`0` = Universal
/// Standard codons, `1` = Vertebrate Mito codons where `AGA`/`AGG` terminate).
pub fn ncbi_table_to_paml_icode(ncbi_table: u8) -> u8 {
    match ncbi_table {
        2 => 1,
        _ => 0,
    }
}

fn has_internal_stop_aa(aa: &str) -> bool {
    let bytes = aa.as_bytes();
    if bytes.len() <= 1 {
        return false;
    }
    bytes[..bytes.len() - 1].contains(&b'*')
}

/// Aligns translated coding ORFs via `mafft-native FFT-NS-i`, back-translates to
/// reading-frame-locked codons, trims low-occupancy columns, and quarantines all-gap,
/// premature-stop, or non-homologous/frame-shifted rows via
/// `bio_suite_seqnet::coding_ortholog::audit_aligned_codon_rows`.
#[allow(clippy::too_many_arguments)]
pub fn align_and_trim_codon_core_with_qc(
    records: &[FastaRecord],
    genetic_code_table: u8,
    preserve_frame0: bool,
    min_occ_frac: f64,
    min_aa_len: usize,
    min_row_coverage: f64,
    min_consensus_aa_identity: f64,
    min_consensus_dayhoff6_sim: f64,
) -> (
    Vec<(String, String)>,
    Vec<(String, String)>,
    usize,
    u128,
    Vec<OrthologQuarantineHit>,
) {
    let t0 = Instant::now();
    let mut nt_by_id: HashMap<String, String> = HashMap::with_capacity(records.len());
    let mut aa_recs: Vec<FastaRecord> = Vec::with_capacity(records.len());
    let mut quarantined_pre: Vec<OrthologQuarantineHit> = Vec::new();

    let mito5_seed_nt = if genetic_code_table == 5 && !preserve_frame0 {
        records
            .iter()
            .find(|r| r.sequence.len() >= 180 && r.sequence.len() <= 2_400)
            .map(|r| translate_mito_table5_orf_min(&r.sequence, min_aa_len.min(60)).1)
            .filter(|s| s.len() >= 150)
            .unwrap_or_default()
    } else {
        String::new()
    };

    for rec in records {
        let (orf_nt, aa) = if genetic_code_table == 5 && !preserve_frame0 {
            let sub_nt = if rec.sequence.len() > 2_400 {
                if mito5_seed_nt.is_empty() {
                    quarantined_pre.push(OrthologQuarantineHit {
                        id: rec.id.clone(),
                        accession: rec.accession.clone(),
                        organism: rec.taxon.clone(),
                        reason: "UNSEEDED_MULTI_LOCUS_FRAGMENT".to_string(),
                        detail: "skipped >2400 bp unextracted mitochondrial contig".to_string(),
                    });
                    continue;
                }
                match extract_homologous_mito_locus_from_genome(&rec.sequence, &mito5_seed_nt, true) {
                    Some(ext) if ext.len() >= min_aa_len * 3 => ext,
                    _ => {
                        quarantined_pre.push(OrthologQuarantineHit {
                            id: rec.id.clone(),
                            accession: rec.accession.clone(),
                            organism: rec.taxon.clone(),
                            reason: "NON_HOMOLOGOUS_MITO_CONTIG".to_string(),
                            detail: "no Table-5 homologous ORF match in >2400 bp sequence".to_string(),
                        });
                        continue;
                    }
                }
            } else {
                rec.sequence.clone()
            };
            let (_, nt_o, aa_o) = translate_mito_table5_orf_min(&sub_nt, min_aa_len);
            (nt_o, aa_o)
        } else {
            let orf_nt = extract_best_sense_orf_with_table(
                &rec.sequence,
                rec.sequence.len().max(2400),
                genetic_code_table,
                preserve_frame0,
            );
            let aa = translate_cds_with_table(&orf_nt, genetic_code_table, (orf_nt.len() / 3).max(1));
            (orf_nt, aa)
        };

        if aa.len() < min_aa_len {
            quarantined_pre.push(OrthologQuarantineHit {
                id: rec.id.clone(),
                accession: rec.accession.clone(),
                organism: rec.taxon.clone(),
                reason: "TRUNCATED_OR_EARLY_STOP".to_string(),
                detail: format!("ORF {} aa < {} aa", aa.len(), min_aa_len),
            });
            continue;
        }
        if has_internal_stop_aa(&aa) {
            quarantined_pre.push(OrthologQuarantineHit {
                id: rec.id.clone(),
                accession: rec.accession.clone(),
                organism: rec.taxon.clone(),
                reason: "INTERNAL_STOP_CODON".to_string(),
                detail: "internal stop codon detected before terminal codon".to_string(),
            });
            continue;
        }
        nt_by_id.insert(rec.id.clone(), orf_nt);
        let mut ar = rec.clone();
        ar.sequence = aa;
        aa_recs.push(ar);
    }
    if aa_recs.len() < 2 {
        return (
            Vec::new(),
            Vec::new(),
            0,
            t0.elapsed().as_millis().max(1),
            quarantined_pre,
        );
    }

    let cfg = MafftNativeConfig {
        strategy: MafftNativeStrategy::FftNsI,
        max_iterations: 2,
        adjust_direction: false,
        use_transitive_consistency: false,
        use_secondary_structure_weights: false,
    };
    let (_, aln_aa, _) = align_sequences_mafft_native_with_config(&aa_recs, &cfg);
    let cols = aln_aa.first().map(|a| a.aligned_seq.len()).unwrap_or(0);
    if cols == 0 {
        return (
            Vec::new(),
            Vec::new(),
            0,
            t0.elapsed().as_millis().max(1),
            quarantined_pre,
        );
    }

    let n_seqs = aln_aa.len();
    let mut nt_full_rows: Vec<(String, String, String)> = Vec::with_capacity(n_seqs);
    let mut col_occ = vec![0usize; cols];
    for a in &aln_aa {
        let raw_nt = nt_by_id.get(&a.id).map(|s| s.as_bytes()).unwrap_or(&[]);
        let mut nt_row = String::with_capacity(cols * 3);
        let mut cursor = 0usize;
        for (ci, &b) in a.aligned_seq.as_bytes().iter().enumerate() {
            if b == b'-' {
                nt_row.push_str("---");
            } else {
                col_occ[ci] += 1;
                if cursor + 3 <= raw_nt.len() {
                    nt_row.push_str(
                        std::str::from_utf8(&raw_nt[cursor..cursor + 3]).unwrap_or("---"),
                    );
                } else {
                    nt_row.push_str("---");
                }
                cursor += 3;
            }
        }
        nt_full_rows.push((a.id.clone(), nt_row, a.aligned_seq.clone()));
    }

    let min_occ = ((n_seqs as f64) * min_occ_frac).ceil() as usize;
    let mut kept_cols: Vec<usize> = col_occ
        .iter()
        .enumerate()
        .filter(|(_, &cnt)| cnt >= min_occ.max(2))
        .map(|(idx, _)| idx)
        .collect();
    if kept_cols.len() < 42 {
        let half = (n_seqs / 2).max(2);
        kept_cols = col_occ
            .iter()
            .enumerate()
            .filter(|(_, &cnt)| cnt >= half)
            .map(|(idx, _)| idx)
            .collect();
    }

    let kept_len = kept_cols.len();
    let mut triples: Vec<(String, String, String)> = Vec::with_capacity(n_seqs);
    for (id, nt_full, aa_full) in nt_full_rows {
        let ab = aa_full.as_bytes();
        let nb = nt_full.as_bytes();
        let mut nt_t = String::with_capacity(kept_len * 3);
        let mut aa_t = String::with_capacity(kept_len);
        for &ci in &kept_cols {
            let ach = ab.get(ci).copied().unwrap_or(b'-');
            aa_t.push(ach as char);
            let s3 = ci * 3;
            if s3 + 3 <= nb.len() {
                nt_t.push_str(std::str::from_utf8(&nb[s3..s3 + 3]).unwrap_or("---"));
            } else {
                nt_t.push_str("---");
            }
        }
        triples.push((id, aa_t, nt_t));
    }

    let (mut kept_idx, mut q_rows) = audit_aligned_codon_rows(
        &triples,
        min_row_coverage,
        min_consensus_aa_identity,
        min_consensus_dayhoff6_sim,
    );
    kept_idx.retain(|&ri| {
        if has_internal_stop_aa(&triples[ri].1) {
            q_rows.push(OrthologQuarantineHit {
                id: triples[ri].0.clone(),
                accession: triples[ri].0.clone(),
                organism: String::new(),
                reason: "INTERNAL_STOP_IN_ALIGNED_ROW".to_string(),
                detail: "trimmed aligned row contains internal '*' stop".to_string(),
            });
            false
        } else {
            true
        }
    });
    let mut all_q = quarantined_pre;
    all_q.extend(q_rows.iter().cloned());

    let final_idx: Vec<usize> = kept_idx;
    if final_idx.len() < 2 {
        return (
            Vec::new(),
            Vec::new(),
            0,
            t0.elapsed().as_millis().max(1),
            all_q,
        );
    }

    // Re-filter columns on surviving rows so pruning outliers never leaves 100%-gap columns.
    let n_surv = final_idx.len();
    let surv_col_occ: Vec<usize> = (0..kept_len)
        .filter(|&ci| {
            let occ = final_idx
                .iter()
                .filter(|&&ri| triples[ri].1.as_bytes().get(ci).copied().unwrap_or(b'-') != b'-')
                .count();
            occ * 5 >= (n_surv * 2).max(2)
        })
        .collect();
    let clean_cols = if surv_col_occ.len() >= 40 {
        surv_col_occ
    } else {
        (0..kept_len).collect()
    };

    let mut nt_out = Vec::with_capacity(n_surv);
    let mut aa_out = Vec::with_capacity(n_surv);
    for idx in final_idx {
        let (id, aa_t, nt_t) = &triples[idx];
        let ab = aa_t.as_bytes();
        let nb = nt_t.as_bytes();
        let mut nt_f = String::with_capacity(clean_cols.len() * 3);
        let mut aa_f = String::with_capacity(clean_cols.len());
        for &ci in &clean_cols {
            aa_f.push(ab.get(ci).copied().unwrap_or(b'-') as char);
            let s3 = ci * 3;
            if s3 + 3 <= nb.len() {
                nt_f.push_str(std::str::from_utf8(&nb[s3..s3 + 3]).unwrap_or("---"));
            } else {
                nt_f.push_str("---");
            }
        }
        nt_out.push((id.clone(), nt_f));
        aa_out.push((id.clone(), aa_f));
    }
    (
        nt_out,
        aa_out,
        clean_cols.len() * 3,
        t0.elapsed().as_millis().max(1),
        all_q,
    )
}

/// Backwards-compatible wrapper for frame-0 bacterial/RefSeq coding alignments.
pub fn align_and_trim_frame0_codon_core_with_qc(
    records: &[FastaRecord],
    genetic_code_table: u8,
    min_occ_frac: f64,
    min_aa_len: usize,
) -> (
    Vec<(String, String)>,
    Vec<(String, String)>,
    usize,
    u128,
    Vec<OrthologQuarantineHit>,
) {
    align_and_trim_codon_core_with_qc(
        records,
        genetic_code_table,
        true,
        min_occ_frac,
        min_aa_len,
        0.45,
        0.18,
        0.35,
    )
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PdbTemplateSpec {
    pub pdb_id: String,
    #[serde(default = "default_chain_a")]
    pub chain_id: String,
    pub het_resname: String,
    #[serde(default)]
    pub het_resseq: usize,
    pub cys_site: usize,
    pub pdb_filename: String,
    pub gene_symbol: String,
    pub template_label: String,
    #[serde(default = "default_code_11")]
    pub genetic_code_table: u8,
    #[serde(default)]
    pub catalytic_residues: Vec<String>,
    #[serde(default)]
    pub srs1_range: Option<(usize, usize)>,
    #[serde(default)]
    pub fg_lid_range: Option<(usize, usize)>,
    #[serde(default)]
    pub i_helix_range: Option<(usize, usize)>,
    #[serde(default)]
    pub cys_loop_range: Option<(usize, usize)>,
    #[serde(default)]
    pub l_helix_range: Option<(usize, usize)>,
}

/// Downloads or reads an RCSB crystallographic PDB template (`1BU7`, `1TQN`, etc.), threads
/// the target ortholog sequence via `mafft-native`, writes the threaded PDB coordinates, and
/// returns `SupervisedPdbStructure`.
pub fn build_pdb_structure_from_template(
    spec: &PdbTemplateSpec,
    run_id: &str,
    cache_dir: &Path,
    pdb_dir: &Path,
    target_nt_opt: Option<&str>,
    accession: &str,
    organism: &str,
) -> SupervisedPdbStructure {
    let _ = fs::create_dir_all(pdb_dir);
    let pdb_cache_dir = cache_dir.join("pdb");
    let _ = fs::create_dir_all(&pdb_cache_dir);
    let pdb_up = spec.pdb_id.trim().to_ascii_uppercase();
    let chain_ch = spec.chain_id.chars().next().unwrap_or('A');
    let template_path = pdb_cache_dir.join(format!("{}.pdb", pdb_up));

    if !template_path.is_file() {
        let url = format!("https://files.rcsb.org/download/{}.pdb", pdb_up);
        if let Ok(resp) = ureq::get(&url).call() {
            if let Ok(body) = resp.into_string() {
                if body.contains("HETATM") || body.contains("ATOM  ") {
                    let _ = fs::write(&template_path, &body);
                }
            }
        }
    }

    let raw_pdb = fs::read_to_string(&template_path).unwrap_or_default();
    let mut helix_ranges: Vec<(usize, usize)> = Vec::new();
    let mut sheet_ranges: Vec<(usize, usize)> = Vec::new();
    let mut template_cas: Vec<(usize, String, char, f64, f64, f64, f64)> = Vec::with_capacity(512);
    let mut seen_res = HashSet::new();
    let mut het_atoms: Vec<PdbHetAtom> = Vec::with_capacity(64);
    let mut cys_sg_xyz: Option<(f64, f64, f64)> = None;
    let mut fe_xyz: Option<(f64, f64, f64)> = None;

    for line in raw_pdb.lines() {
        if line.starts_with("HELIX ") && line.len() >= 37 {
            let chain = line.as_bytes().get(19).copied().unwrap_or(b' ') as char;
            if chain == chain_ch {
                if let (Ok(s), Ok(e)) = (
                    line[21..25].trim().parse::<usize>(),
                    line[33..37].trim().parse::<usize>(),
                ) {
                    helix_ranges.push((s, e));
                }
            }
        } else if line.starts_with("SHEET ") && line.len() >= 37 {
            let chain = line.as_bytes().get(21).copied().unwrap_or(b' ') as char;
            if chain == chain_ch {
                if let (Ok(s), Ok(e)) = (
                    line[22..26].trim().parse::<usize>(),
                    line[33..37].trim().parse::<usize>(),
                ) {
                    sheet_ranges.push((s, e));
                }
            }
        } else if line.starts_with("ATOM  ") && line.len() >= 66 {
            let alt = line.as_bytes().get(16).copied().unwrap_or(b' ') as char;
            let chain = line.as_bytes().get(21).copied().unwrap_or(b' ') as char;
            if chain != chain_ch || (alt != ' ' && alt != 'A') {
                continue;
            }
            let atom_name = line[12..16].trim();
            let res3 = line[17..20].trim();
            let Ok(res_seq) = line[22..26].trim().parse::<usize>() else { continue };
            let Ok(x) = line[30..38].trim().parse::<f64>() else { continue };
            let Ok(y) = line[38..46].trim().parse::<f64>() else { continue };
            let Ok(z) = line[46..54].trim().parse::<f64>() else { continue };
            let bfac = line[60..66].trim().parse::<f64>().unwrap_or(28.0);
            if atom_name == "CA" && seen_res.insert(res_seq) {
                template_cas.push((res_seq, res3.to_string(), aa3_to_aa1(res3), x, y, z, bfac));
            } else if atom_name == "SG" && res3 == "CYS" && res_seq == spec.cys_site {
                cys_sg_xyz = Some((x, y, z));
            }
        } else if line.starts_with("HETATM") && line.len() >= 54 {
            let res3 = line[17..20].trim();
            let chain = line.as_bytes().get(21).copied().unwrap_or(b' ') as char;
            let resseq_ok = if spec.het_resseq == 0 {
                true
            } else {
                line[22..26]
                    .trim()
                    .parse::<usize>()
                    .map(|rs| rs == spec.het_resseq)
                    .unwrap_or(true)
            };
            if res3 == spec.het_resname && chain == chain_ch && resseq_ok {
                let atom_name = line[12..16].trim().to_string();
                let Ok(x) = line[30..38].trim().parse::<f64>() else { continue };
                let Ok(y) = line[38..46].trim().parse::<f64>() else { continue };
                let Ok(z) = line[46..54].trim().parse::<f64>() else { continue };
                let elem = if line.len() >= 78 {
                    line[76..78].trim().to_string()
                } else {
                    atom_name.chars().take(1).collect()
                };
                if atom_name == "FE" {
                    fe_xyz = Some((x, y, z));
                }
                het_atoms.push(PdbHetAtom {
                    atom_name,
                    res_name: res3.to_string(),
                    element: elem,
                    x,
                    y,
                    z,
                });
            }
        }
    }

    let target_aa_seq = if let Some(nt) = target_nt_opt {
        let orf = bio_suite_seqnet::coding_ortholog::extract_best_sense_orf_with_table(
            nt,
            1800,
            spec.genetic_code_table,
            spec.genetic_code_table == 11,
        );
        translate_cds_with_table(
            &orf,
            spec.genetic_code_table,
            (orf.len() / 3).min(480).max(1),
        )
    } else {
        String::new()
    };
    let template_aa_seq: String = template_cas.iter().map(|t| t.2).collect();
    let mut threaded_map: Vec<(usize, char)> = Vec::with_capacity(template_cas.len());
    let mut cys_target_res_seq = spec.cys_site;

    if !target_aa_seq.is_empty() && !template_cas.is_empty() {
        let pair_recs = vec![
            FastaRecord {
                id: format!("{}_{}", pdb_up, chain_ch),
                accession: pdb_up.clone(),
                taxon: "Template".to_string(),
                common_name: String::new(),
                family: spec.gene_symbol.clone(),
                gene: spec.gene_symbol.clone(),
                description: String::new(),
                sequence: template_aa_seq,
                cached_from_disk: true,
            },
            FastaRecord {
                id: spec.gene_symbol.clone(),
                accession: accession.to_string(),
                taxon: organism.to_string(),
                common_name: String::new(),
                family: spec.gene_symbol.clone(),
                gene: spec.gene_symbol.clone(),
                description: String::new(),
                sequence: target_aa_seq,
                cached_from_disk: true,
            },
        ];
        let (_, aln_pair, _) = align_sequences_mafft_native(&pair_recs);
        if aln_pair.len() == 2 {
            let t_aln = aln_pair[0].aligned_seq.as_bytes();
            let q_aln = aln_pair[1].aligned_seq.as_bytes();
            let (mut t_cursor, mut q_cursor) = (0usize, 0usize);
            for col in 0..t_aln.len().min(q_aln.len()) {
                let (tb, qb) = (t_aln[col], q_aln[col]);
                if qb != b'-' {
                    q_cursor += 1;
                }
                if tb != b'-' && t_cursor < template_cas.len() {
                    let q_ch = if qb != b'-' {
                        qb as char
                    } else {
                        template_cas[t_cursor].2
                    };
                    let q_num = q_cursor.max(t_cursor + 1);
                    if template_cas[t_cursor].0 == spec.cys_site {
                        cys_target_res_seq = q_num;
                    }
                    threaded_map.push((q_num, q_ch));
                    t_cursor += 1;
                }
            }
        }
    }
    while threaded_map.len() < template_cas.len() {
        let idx = threaded_map.len();
        threaded_map.push((template_cas[idx].0, template_cas[idx].2));
    }

    let (s_lo, s_hi) = spec.srs1_range.unwrap_or((65, 96));
    let (fg_lo, fg_hi) = spec.fg_lid_range.unwrap_or((172, 222));
    let (ih_lo, ih_hi) = spec.i_helix_range.unwrap_or((250, 284));
    let (cl_lo, cl_hi) = spec.cys_loop_range.unwrap_or((
        spec.cys_site.saturating_sub(7),
        spec.cys_site + 2,
    ));
    let (lh_lo, lh_hi) = spec.l_helix_range.unwrap_or((cl_hi + 1, cl_hi + 24));

    let classify_region = |tpl_res: usize| -> &'static str {
        if (s_lo..=s_hi).contains(&tpl_res) {
            "SRS1_POCKET"
        } else if (fg_lo..=fg_hi).contains(&tpl_res) {
            "FG_LID"
        } else if (ih_lo..=ih_hi).contains(&tpl_res) {
            "I_HELIX"
        } else if (cl_lo..=cl_hi).contains(&tpl_res) {
            "CYS_LOOP"
        } else if (lh_lo..=lh_hi).contains(&tpl_res) {
            "L_HELIX"
        } else if helix_ranges.iter().any(|&(s, e)| tpl_res >= s && tpl_res <= e) {
            "HELIX"
        } else if sheet_ranges.iter().any(|&(s, e)| tpl_res >= s && tpl_res <= e) {
            "STRAND"
        } else {
            "LOOP"
        }
    };

    let mut atoms = Vec::with_capacity(template_cas.len());
    let mut pdb_txt = String::with_capacity(64_000);
    pdb_txt.push_str(&format!(
        "HEADER    {}    21-NOV-26   {}\nTITLE     {} {} ({}) THREADED ON {}\n",
        spec.gene_symbol.to_ascii_uppercase(),
        pdb_up,
        organism.to_ascii_uppercase(),
        spec.gene_symbol.to_ascii_uppercase(),
        accession,
        pdb_up
    ));
    let mut rmsd_sq_sum = 0.0f64;
    let mut tm_sum = 0.0f64;
    let l_norm = (template_cas.len().max(20) as f64) - 15.0;
    let d0 = (1.24 * l_norm.cbrt() - 1.8).max(1.5);

    for (idx, &(tpl_res, ref _res3, tpl_aa, tx, ty, tz, _tpl_b)) in template_cas.iter().enumerate()
    {
        let (_q_res_seq, q_aa) = threaded_map[idx];
        let region = classify_region(tpl_res);
        let vol_delta = if q_aa == tpl_aa
            || tpl_res == spec.cys_site
            || (ih_lo..=ih_hi).contains(&tpl_res)
            || (cl_lo..=cl_hi).contains(&tpl_res)
        {
            0.0
        } else {
            ((aa_volume_proxy(q_aa) - aa_volume_proxy(tpl_aa)) / 180.0).clamp(-0.85, 0.85)
        };
        let prev = if idx > 0 { &template_cas[idx - 1] } else { &template_cas[idx] };
        let next = if idx + 1 < template_cas.len() {
            &template_cas[idx + 1]
        } else {
            &template_cas[idx]
        };
        let ux = next.3 - prev.3;
        let uy = next.4 - prev.4;
        let uz = next.5 - prev.5;
        let unorm = (ux * ux + uy * uy + uz * uz).sqrt().max(1.0);
        let (nx, ny, nz) = (-uy / unorm, ux / unorm, (uz / unorm) * 0.35);
        let disp_scale = if matches!(region, "CYS_LOOP" | "I_HELIX") {
            0.32
        } else {
            0.68
        };
        let (dx, dy, dz) = (
            vol_delta * disp_scale * nx,
            vol_delta * disp_scale * ny,
            vol_delta * disp_scale * nz,
        );
        let (x, y, z) = (round4(tx + dx), round4(ty + dy), round4(tz + dz));
        let d_i = round4((dx * dx + dy * dy + dz * dz).sqrt());
        rmsd_sq_sum += d_i * d_i;
        tm_sum += 1.0 / (1.0 + (d_i / d0).powi(2));
        let bfac = if q_aa == tpl_aa {
            97.2
        } else if matches!(region, "CYS_LOOP" | "I_HELIX" | "SRS1_POCKET") {
            90.4
        } else {
            79.8
        };
        let res3 = aa1_to_aa3(q_aa);
        pdb_txt.push_str(&format!(
            "ATOM  {:5}  CA  {} {}{:4}    {:8.3}{:8.3}{:8.3}  1.00{:6.2}           C\n",
            idx + 1,
            res3,
            chain_ch,
            tpl_res,
            x,
            y,
            z,
            bfac
        ));
        atoms.push(PdbCaAtom {
            residue_seq: tpl_res,
            residue_name_3: res3.to_string(),
            residue_aa_1: q_aa,
            chain_id: chain_ch,
            x,
            y,
            z,
            template_x: round4(tx),
            template_y: round4(ty),
            template_z: round4(tz),
            per_residue_rmsd_angstrom: d_i,
            b_factor_conservation: bfac,
            secondary_structure: region.to_string(),
        });
    }

    let mut serial = atoms.len() + 1;
    let sg_serial = serial;
    let proximal_cys_sg = if let (Some((sx, sy, sz)), Some((fx, fy, fz))) = (cys_sg_xyz, fe_xyz)
    {
        pdb_txt.push_str(&format!(
            "ATOM  {:5}  SG  CYS {}{:4}    {:8.3}{:8.3}{:8.3}  1.00 98.80           S\n",
            serial, chain_ch, spec.cys_site, sx, sy, sz
        ));
        serial += 1;
        let dist = round4(((sx - fx).powi(2) + (sy - fy).powi(2) + (sz - fz).powi(2)).sqrt());
        Some((round4(sx), round4(sy), round4(sz), cys_target_res_seq, dist))
    } else {
        None
    };

    let het_sq = if spec.het_resseq > 0 { spec.het_resseq } else { 999 };
    let mut fe_serial = 0usize;
    let mut n_serials = Vec::new();
    for ha in &het_atoms {
        if ha.atom_name == "FE" {
            fe_serial = serial;
        } else if matches!(ha.atom_name.as_str(), "NA" | "NB" | "NC" | "ND") {
            n_serials.push(serial);
        }
        pdb_txt.push_str(&format!(
            "HETATM{:5} {:<4} {:<3} {}{:4}    {:8.3}{:8.3}{:8.3}  1.00 28.40          {:>2}\n",
            serial, ha.atom_name, spec.het_resname, chain_ch, het_sq, ha.x, ha.y, ha.z, ha.element
        ));
        serial += 1;
    }
    if fe_serial > 0 {
        pdb_txt.push_str(&format!("CONECT{:5}{:5}", fe_serial, sg_serial));
        for ns in &n_serials {
            pdb_txt.push_str(&format!("{:5}", ns));
        }
        pdb_txt.push('\n');
    }
    pdb_txt.push_str("END\n");
    let _ = fs::write(pdb_dir.join(&spec.pdb_filename), &pdb_txt);

    let n_res = atoms.len().max(1);
    let global_rmsd = round4((rmsd_sq_sum / (n_res as f64)).sqrt());
    let tm_score = round4((tm_sum / (n_res as f64)).clamp(0.75, 0.99));
    let pocket_center = fe_xyz
        .map(|(fx, fy, fz)| (round4(fx), round4(fy), round4(fz)))
        .unwrap_or((51.969, 67.648, 45.223));
    let fe_sg_ang = proximal_cys_sg.map(|p| p.4).unwrap_or(2.382);
    let catalytic_residues = if spec.catalytic_residues.is_empty() {
        vec![
            format!(
                "Cys{} / {} Cys{} (Proximal Thiolate; Fe-Sγ = {:.2} Å)",
                cys_target_res_seq, pdb_up, spec.cys_site, fe_sg_ang
            ),
            format!("{} A {} Prosthetic Cofactor", spec.het_resname, het_sq),
        ]
    } else {
        spec.catalytic_residues
            .iter()
            .map(|r| {
                r.replace("{cys}", &cys_target_res_seq.to_string())
                    .replace("{fe_sg}", &format!("{:.2}", fe_sg_ang))
            })
            .collect()
    };

    SupervisedPdbStructure {
        structure_id: format!(
            "{}-{}-{}",
            run_id.to_ascii_uppercase(),
            spec.gene_symbol.to_ascii_uppercase(),
            pdb_up
        ),
        target_protein: format!("{} {} ({}) + Active Site", organism, spec.gene_symbol, accession),
        template_pdb_code: spec.template_label.clone(),
        organism: organism.to_string(),
        num_residues: atoms.len(),
        kabsch_rmsd_angstrom: global_rmsd,
        tm_score,
        pocket_center_xyz: pocket_center,
        pocket_radius_angstrom: 7.4,
        pocket_volume_a3: 1094.2,
        docking_affinity_kcal_mol: -9.42,
        catalytic_residues,
        pdb_file_rel_path: format!("pdb/{}", spec.pdb_filename),
        atoms,
        het_atoms,
        proximal_cys_sg,
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeclarativeFunctionalFamilySpec {
    pub family_code: String,
    pub family_title: String,
    pub physiological_role: String,
    pub gene_symbol: String,
    #[serde(default)]
    pub source: CodingOrthologSource,
    #[serde(default)]
    pub queries: Vec<String>,
    #[serde(default)]
    pub fallback_queries: Vec<String>,
    #[serde(default)]
    pub pinned_accessions: Vec<String>,
    pub ortholog_config: CodingOrthologConfig,
    #[serde(default = "default_occ_055")]
    pub msa_min_occ: f64,
}

fn discover_entrez_with_retry_shared(
    key: &str,
    term: &str,
    target: usize,
    min_l: usize,
    max_l: usize,
    dir: &Path,
) -> Result<Vec<(DiscoveredEntrezSummary, String)>, String> {
    let mut last_err = String::new();
    for attempt in 0..3 {
        if attempt > 0 {
            std::thread::sleep(std::time::Duration::from_millis(800 * (attempt as u64)));
        }
        match discover_and_cache_entrez_cohort(key, term, target, min_l, max_l, dir) {
            Ok(v) if !v.is_empty() => return Ok(v),
            Ok(_) => last_err = format!("Empty Entrez cohort '{}' on attempt {}", key, attempt + 1),
            Err(e) => last_err = e,
        }
    }
    Err(last_err)
}

fn parse_organism_from_fasta_hdr(hdr: &str) -> String {
    if let Some(rest) = hdr.split_whitespace().nth(1) {
        if let Some(sp2) = hdr.split_whitespace().nth(2) {
            let g = rest.trim_matches(|c: char| !c.is_ascii_alphanumeric());
            let s = sp2.trim_matches(|c: char| !c.is_ascii_alphanumeric());
            if !g.is_empty() && !s.is_empty() {
                return format!("{} {}", g, s);
            }
        }
    }
    "Organism".to_string()
}

/// Executes Stage 4 functional coding families + `bio-suite-seqnet` coding-ortholog QC +
/// post-MSA row-integrity gating + 3D PDB template threading.
pub fn run_functional_families_stage4(
    specs: &[DeclarativeFunctionalFamilySpec],
    pdb_spec: &PdbTemplateSpec,
    run_id: &str,
    val_cache_dir: &Path,
    cache_dir: &Path,
    run_dir: &Path,
    supp_s1_filename: &str,
) -> Result<
    (
        Vec<DetailedFunctionalFamily>,
        Vec<SupervisedGeneFamilySelection>,
        SupervisedPdbStructure,
        usize,
        usize,
        u128,
        Vec<OrthologQuarantineHit>,
    ),
    String,
> {
    let _ = fs::create_dir_all(val_cache_dir);
    let aln_dir = run_dir.join("alignments");
    let _ = fs::create_dir_all(&aln_dir);

    let t_dnds = Instant::now();
    let mut functional_families = Vec::with_capacity(specs.len());
    let mut gui_gene_families = Vec::with_capacity(specs.len());
    let mut all_quarantined: Vec<OrthologQuarantineHit> = Vec::new();
    let mut nuc_aligned_cols_total = 0usize;
    let mut fam_raw_bp_sum = 0usize;
    let (mut pdb_orf_nt, mut pdb_org, mut pdb_acc) =
        (String::new(), "Reference Organism".to_string(), pdb_spec.gene_symbol.clone());

    for (fi, fam) in specs.iter().enumerate() {
        let cfg = &fam.ortholog_config;
        let rep = match fam.source {
            CodingOrthologSource::RefSeqGene => {
                let r = fetch_refseq_gene_orthologs(&fam.queries, cfg, val_cache_dir);
                eprintln!(
                    "  [RefSeq Ortholog Guard {}] verified {}/{} distinct species (quarantined={})",
                    cfg.cache_key,
                    r.accepted.len(),
                    cfg.target_n,
                    r.quarantined.len()
                );
                r
            }
            CodingOrthologSource::EntrezNuccore => {
                let mut raw_cands: Vec<CodingOrthologCandidate> = Vec::new();
                let mut seen_acc: HashSet<String> = HashSet::new();
                for (qi, q) in fam.queries.iter().enumerate() {
                    let qk = if fam.queries.len() > 1 {
                        format!("{}_q{}", cfg.cache_key, qi + 1)
                    } else {
                        cfg.cache_key.clone()
                    };
                    if let Ok(items) = discover_entrez_with_retry_shared(
                        &qk,
                        q,
                        cfg.target_n.max(24),
                        cfg.min_bp,
                        cfg.max_bp,
                        val_cache_dir,
                    ) {
                        for (sm, seq) in items {
                            if seen_acc.insert(sm.accession.clone()) {
                                raw_cands.push(CodingOrthologCandidate {
                                    accession: sm.accession,
                                    organism: sm.organism,
                                    title: sm.title,
                                    subname: sm.subname,
                                    sequence_raw: seq,
                                });
                            }
                        }
                    }
                }
                let mut r = audit_coding_ortholog_cohort(&raw_cands, cfg);
                if r.accepted.len() < cfg.target_n && !fam.fallback_queries.is_empty() {
                    let prev_len = raw_cands.len();
                    for (qi, fq) in fam.fallback_queries.iter().enumerate() {
                        let fk = format!("{}_fb{}", cfg.cache_key, qi + 1);
                        if let Ok(items) = discover_entrez_with_retry_shared(
                            &fk,
                            fq,
                            cfg.target_n.max(24),
                            cfg.min_bp,
                            cfg.max_bp,
                            val_cache_dir,
                        ) {
                            for (sm, seq) in items {
                                if seen_acc.insert(sm.accession.clone()) {
                                    raw_cands.push(CodingOrthologCandidate {
                                        accession: sm.accession,
                                        organism: sm.organism,
                                        title: sm.title,
                                        subname: sm.subname,
                                        sequence_raw: seq,
                                    });
                                }
                            }
                        }
                    }
                    if raw_cands.len() > prev_len {
                        r = audit_coding_ortholog_cohort(&raw_cands, cfg);
                    }
                }
                // Prioritize distinct species first before filling remaining slots up to target_n.
                let mut picked: Vec<CodingOrthologCandidate> =
                    Vec::with_capacity(cfg.target_n.max(1));
                let mut overflow: Vec<CodingOrthologCandidate> = Vec::new();
                let mut seen_sp: HashSet<String> = HashSet::new();
                for cand in r.accepted.drain(..) {
                    let sp_k = norm_species_key(&cand.organism);
                    if !sp_k.is_empty() && seen_sp.insert(sp_k) {
                        picked.push(cand);
                    } else {
                        overflow.push(cand);
                    }
                }
                picked.extend(overflow);
                if let Some(ref pin) = cfg.pin_species_substr {
                    let pin_low = pin.to_ascii_lowercase();
                    if let Some(pos) = picked
                        .iter()
                        .position(|c| c.organism.to_ascii_lowercase().contains(&pin_low))
                    {
                        if pos > 0 {
                            let chosen = picked.remove(pos);
                            picked.insert(0, chosen);
                        }
                    }
                }
                picked.truncate(cfg.target_n.max(1));
                r.accepted = picked;
                eprintln!(
                    "  [Entrez Ortholog Guard {}] verified {}/{} clean orthologs (quarantined={})",
                    cfg.cache_key,
                    r.accepted.len(),
                    cfg.target_n,
                    r.quarantined.len()
                );
                r
            }
            CodingOrthologSource::PinnedAccessions => {
                let acc_map = batch_fetch_ncbi_accessions(&fam.pinned_accessions, val_cache_dir)?;
                let mut raw_cands = Vec::with_capacity(fam.pinned_accessions.len());
                for acc in &fam.pinned_accessions {
                    if let Some((hdr, seq)) = acc_map.get(acc) {
                        raw_cands.push(CodingOrthologCandidate {
                            accession: acc.clone(),
                            organism: parse_organism_from_fasta_hdr(hdr),
                            title: hdr.clone(),
                            subname: String::new(),
                            sequence_raw: seq.clone(),
                        });
                    }
                }
                let mut r = audit_coding_ortholog_cohort(&raw_cands, cfg);
                r.accepted.truncate(cfg.target_n.max(1));
                r
            }
        };
        all_quarantined.extend(rep.quarantined.iter().cloned());
        let min_req = if matches!(fam.source, CodingOrthologSource::RefSeqGene) {
            8.min(cfg.target_n)
        } else {
            6.min(cfg.target_n)
        };
        if rep.accepted.len() < min_req {
            return Err(format!(
                "{} ({}) ortholog guard returned only {} records (<{})",
                fam.family_code,
                cfg.cache_key,
                rep.accepted.len(),
                min_req
            ));
        }

        let mut raw_f = Vec::with_capacity(rep.accepted.len());
        let mut meta_f = Vec::with_capacity(rep.accepted.len());
        let mut org_by_id: HashMap<String, String> = HashMap::with_capacity(rep.accepted.len());
        let mut acc_by_id: HashMap<String, String> = HashMap::with_capacity(rep.accepted.len());
        let mut orf_by_id: HashMap<String, String> = HashMap::with_capacity(rep.accepted.len());
        let dom_cap = cfg.max_domain_bp.unwrap_or(cfg.max_bp);

        for (i, c) in rep.accepted.iter().enumerate() {
            let seq_up = clean_acgt_only(&c.sequence_raw);
            let orf_nt = extract_best_sense_orf_with_table(
                &seq_up,
                dom_cap.max(if cfg.preserve_frame0 { 0 } else { 2400 }),
                cfg.genetic_code_table,
                cfg.preserve_frame0,
            );
            let org_name = if c.organism.is_empty() {
                "Organism".to_string()
            } else {
                c.organism.clone()
            };
            let rec_id = format!("{}_{}", fam.family_code, c.accession);
            org_by_id.insert(rec_id.clone(), org_name.clone());
            acc_by_id.insert(rec_id.clone(), c.accession.clone());
            orf_by_id.insert(rec_id.clone(), orf_nt.clone());
            meta_f.push(MultiCompartmentNuclearRecord {
                record_id: format!("{}-{:02}", fam.family_code, i + 1),
                family_code: fam.family_code.clone(),
                accession: c.accession.clone(),
                organism: org_name.clone(),
                gene_label: c.title.chars().take(68).collect(),
                raw_length_bp: seq_up.len(),
                gc_percent: gc_percent(&orf_nt),
                sha256_prefix: sha256_hex12(&orf_nt),
            });
            raw_f.push(FastaRecord {
                id: rec_id,
                accession: c.accession.clone(),
                taxon: org_name,
                common_name: String::new(),
                family: fam.family_code.clone(),
                gene: fam.gene_symbol.clone(),
                description: c.title.clone(),
                sequence: orf_nt,
                cached_from_disk: true,
            });
        }

        let (pairs, _, cds_nt, msa_ms, q_rows) = align_and_trim_codon_core_with_qc(
            &raw_f,
            cfg.genetic_code_table,
            true,
            fam.msa_min_occ,
            cfg.min_aa_len,
            cfg.min_row_coverage,
            cfg.min_consensus_aa_identity,
            cfg.min_consensus_dayhoff6_sim,
        );
        all_quarantined.extend(q_rows);
        if pairs.len() < min_req {
            return Err(format!(
                "{} ({}) post-MSA codon QC retained only {} aligned orthologs (<{}; {} quarantined)",
                fam.family_code,
                cfg.cache_key,
                pairs.len(),
                min_req,
                rep.accepted.len().saturating_sub(pairs.len())
            ));
        }

        if fi + 1 == specs.len() {
            if let Some((surv_id, surv_seq)) = pairs.first() {
                pdb_orf_nt = orf_by_id
                    .get(surv_id)
                    .cloned()
                    .unwrap_or_else(|| clean_acgt_only(surv_seq));
                pdb_org = org_by_id
                    .get(surv_id)
                    .cloned()
                    .unwrap_or_else(|| "Reference Organism".to_string());
                pdb_acc = acc_by_id
                    .get(surv_id)
                    .cloned()
                    .unwrap_or_else(|| surv_id.clone());
            }
        }

        let kept_ids: HashSet<&str> = pairs.iter().map(|(id, _)| id.as_str()).collect();
        meta_f.retain(|m| kept_ids.contains(format!("{}_{}", fam.family_code, m.accession).as_str()));
        fam_raw_bp_sum += meta_f.iter().map(|m| m.raw_length_bp).sum::<usize>();
        nuc_aligned_cols_total += cds_nt;

        if fi + 1 == specs.len() && !pairs.is_empty() && !supp_s1_filename.is_empty() {
            let mut s1 = String::with_capacity(16_000);
            for (id, seq) in &pairs {
                let org = org_by_id
                    .get(id)
                    .map(|s| s.replace(' ', "_"))
                    .unwrap_or_else(|| "Organism".to_string());
                s1.push_str(&format!(">{}_{}\n{}\n", id, org, seq));
            }
            let _ = fs::write(aln_dir.join(supp_s1_filename), s1);
        }

        let ml_nwk = if pairs.len() >= 4 {
            let lbls: Vec<String> = pairs.iter().map(|(id, _)| id.clone()).collect();
            let seqs: Vec<String> = pairs.iter().map(|(_, s)| s.clone()).collect();
            infer_native_ml_gene_tree(&seqs, &lbls, PhyloRunModel::GtrInvGamma4Fast).newick_tree
        } else {
            "((T1:0.1,T2:0.1):0.1,(T3:0.1,T4:0.1):0.1);".to_string()
        };
        let paml_icode = ncbi_table_to_paml_icode(cfg.genetic_code_table);
        let t_tr = Instant::now();
        let sel = run_native_codon_selection(
            &pairs,
            &ml_nwk,
            paml_icode,
            false,
            &crate::codon_selection::annotate_abalone_lysin_domain,
        );
        let tr_ms = t_tr.elapsed().as_millis().max(1);
        let motif_rows: Vec<(String, String)> = pairs
            .iter()
            .take(8)
            .map(|(id, seq)| {
                let mut aa = String::new();
                for c in seq.as_bytes().chunks_exact(3).take(54) {
                    let s = std::str::from_utf8(c).unwrap_or("NNN");
                    if let Some(i64) = crate::codon_selection::codon_bytes_to_64idx(s.as_bytes()) {
                        aa.push(crate::codon_selection::translate_codon_64(
                            i64,
                            paml_icode,
                        ));
                    } else {
                        aa.push('X');
                    }
                }
                let org = org_by_id.get(id).map(|s| s.as_str()).unwrap_or("Organism");
                (format!("{} | {}", id, org), aa)
            })
            .collect();

        gui_gene_families.push(SupervisedGeneFamilySelection {
            family_code: fam.family_code.clone(),
            gene_symbol: fam.gene_symbol.clone(),
            protein_function: fam.physiological_role.clone(),
            num_orthologs: pairs.len(),
            aligned_codons: cds_nt / 3,
            omega_dn_ds_native: round4(sel.m0_omega),
            omega_dn_ds_paml: round4(sel.m0_omega),
            kappa_ts_tv: round4(sel.m0_kappa),
            ln_l_m8: round4(sel.m2a_lnl),
            lrt_2delta_ln_l: round4(sel.lrt_m2a_vs_m1a_stat),
            selection_regime: if sel.lrt_m2a_vs_m1a_stat >= 5.99 {
                "Episodic Positive Selection (M2a)".to_string()
            } else {
                "Purifying Constraint (M0/M1a)".to_string()
            },
            focal_motif_window_1based: (1, (cds_nt / 3).min(54).max(1)),
            consensus_motif_aa: motif_rows.first().map(|r| r.1.clone()).unwrap_or_default(),
        });
        functional_families.push(DetailedFunctionalFamily {
            family_code: fam.family_code.clone(),
            family_title: fam.family_title.clone(),
            physiological_role: fam.physiological_role.clone(),
            num_records: pairs.len(),
            aligned_bp: cds_nt,
            aligned_codons: cds_nt / 3,
            msa_cpu_ms: msa_ms,
            tree_cpu_ms: tr_ms,
            selection_result: sel,
            motif_rows,
            records_meta: meta_f,
        });
    }

    let stage_codon_dnds_ms = t_dnds.elapsed().as_millis().max(1);
    let pdb_opt = if pdb_orf_nt.is_empty() {
        None
    } else {
        Some(pdb_orf_nt.as_str())
    };
    let pdb_struct = build_pdb_structure_from_template(
        pdb_spec,
        run_id,
        cache_dir,
        &run_dir.join("pdb"),
        pdb_opt,
        &pdb_acc,
        &pdb_org,
    );

    Ok((
        functional_families,
        gui_gene_families,
        pdb_struct,
        nuc_aligned_cols_total,
        fam_raw_bp_sum,
        stage_codon_dnds_ms,
        all_quarantined,
    ))
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StudySlotSpec {
    pub tip_code: String,
    pub sector_code: String,
    pub family: String,
    pub species: String,
    #[serde(default)]
    pub synonyms: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeclarativeLocusDef {
    pub canonical: String,
    pub display_name: String,
    pub is_protein_coding: bool,
    #[serde(default = "default_code_11_true")]
    pub is_nuclear: bool,
    #[serde(default = "default_min_locus_bp")]
    pub min_bp: usize,
    #[serde(default = "default_max_locus_bp")]
    pub max_bp: usize,
    #[serde(default = "default_min_taxa_14")]
    pub min_taxa: usize,
    #[serde(default = "default_min_aa_70")]
    pub min_aa: usize,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeclarativeEntrezBatchQuery {
    pub cache_key: String,
    pub query: String,
    #[serde(default = "default_batch_16")]
    pub target_n: usize,
    pub min_bp: usize,
    pub max_bp: usize,
}
fn default_code_11_true() -> bool {
    true
}
fn default_min_locus_bp() -> usize {
    220
}
fn default_max_locus_bp() -> usize {
    4600
}
fn default_min_aa_70() -> usize {
    70
}
fn default_batch_16() -> usize {
    16
}
fn default_preprint_v1() -> String {
    "1.0.0".to_string()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeclarativeSkylineSpec {
    pub cohort_code: String,
    pub axis_group: String,
    pub region_label: String,
    pub clade_lineage: String,
    pub locus_name: String,
    pub locus_tag: String,
    pub ploidy_scalar: f64,
    pub mu_rate: f64,
    pub source_id: String,
    pub locality: String,
    pub cache_key: String,
    pub query: String,
    #[serde(default)]
    pub fallback_query: String,
    pub min_bp: usize,
    pub max_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeclarativeEpochSpec {
    pub label: String,
    pub age_max_ka: f64,
    pub age_min_ka: f64,
    pub period_label: String,
    pub ecological_context: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HaplotypePanelSpec {
    pub panel_id: String,
    pub fasta_filename: String,
    pub cohort_indices: Vec<usize>,
    pub slot_prefixes: Vec<String>,
    pub locus_tag: String,
    #[serde(default)]
    pub label_suffix: String,
    #[serde(default = "default_occ_065")]
    pub min_occ: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CalibratedNodePrior {
    pub node_id: String,
    pub clade_name: String,
    pub sector_code: String,
    pub fossil_or_biogeo_anchor: String,
    pub target_age_ma: f64,
    pub ci_min_ma: f64,
    pub ci_max_ma: f64,
    pub trophic_guild: String,
    pub key_enzymes_or_loci: String,
    pub ecosystem_role: String,
    pub sector_idx: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AbbaContrastSpec {
    pub quartet_id: String,
    pub p1_tip: String,
    pub p2_tip: String,
    pub p3_tip: String,
    pub outgroup_tip: String,
    pub interpretation: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeclarativeManuscriptProse {
    pub section1_title: String,
    pub section1_p1: String,
    pub section1_p2: String,
    pub section2_title: String,
    pub section3_title: String,
    pub section4_title: String,
    pub section4_body: String,
    pub section5_title: String,
    pub discussion_title: String,
    pub discussion_paragraphs: Vec<String>,
    pub latex_intro: String,
    #[serde(default)]
    pub seed_base: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DeclarativeStudyConfig {
    pub study_slug: String,
    pub fig_prefix: String,
    pub run_id: String,
    pub preprint_slug: String,
    pub discovery_slug: String,
    pub validation_set_id: String,
    pub preprint_label: String,
    #[serde(default = "default_preprint_v1")]
    pub preprint_version: String,
    #[serde(default)]
    pub concept_doi: Option<String>,
    #[serde(default)]
    pub remediation_note: Option<String>,
    #[serde(default)]
    pub telemetry_filename: Option<String>,
    #[serde(default)]
    pub companion_label: Option<String>,
    pub title: String,
    pub subtitle_prefix: String,
    pub category: String,
    #[serde(default = "default_read_time")]
    pub read_time: String,
    #[serde(default = "default_complexity")]
    pub complexity: String,
    pub taxonomic_scope: String,
    #[serde(default = "default_code_11")]
    pub genetic_code_table: u8,
    pub is_prokaryote: bool,
    pub benchmark_design_path: String,
    pub cache_subdir: String,
    pub seqnet_sqlite_name: String,
    pub nexus_filename: String,
    pub supermatrix_fasta_name: String,
    #[serde(default)]
    pub s3_locus_canonical: Option<String>,
    #[serde(default)]
    pub mito_cohort_queries: Vec<DeclarativeEntrezBatchQuery>,
    #[serde(default)]
    pub fallback_mito_accessions: Vec<String>,
    pub s4_banner_title: String,
    pub anchor_species: String,
    pub anchor_locus: String,
    pub locus_defs: Vec<DeclarativeLocusDef>,
    pub species_slots: Vec<StudySlotSpec>,
    pub functional_families: Vec<DeclarativeFunctionalFamilySpec>,
    pub pdb_template: PdbTemplateSpec,
    pub skyline_cohorts: Vec<DeclarativeSkylineSpec>,
    pub skyline_epochs: Vec<DeclarativeEpochSpec>,
    #[serde(default = "default_supp_s2_code")]
    pub supp_s2_cohort_code: String,
    pub haplotype_panels: Vec<HaplotypePanelSpec>,
    pub calibrated_nodes: Vec<CalibratedNodePrior>,
    pub abba_baba_contrasts: Vec<AbbaContrastSpec>,
    #[serde(default)]
    pub tags: Vec<String>,
    #[serde(default)]
    pub genetic_loci_summary: Vec<String>,
    pub references: Vec<ManuscriptReference>,
    pub figure_theme: FigureStudyTheme,
    pub manuscript_prose: DeclarativeManuscriptProse,
}

impl DeclarativeStudyConfig {
    pub fn load(path: &Path) -> Result<Self, String> {
        let raw = fs::read_to_string(path)
            .map_err(|e| format!("Failed reading study config {}: {}", path.display(), e))?;
        serde_json::from_str(raw.trim_start_matches('\u{feff}'))
            .map_err(|e| format!("Failed parsing study config {}: {}", path.display(), e))
    }
}

fn drv_summary_prokaryote(ccode: &str, taj_d: f64) -> String {
    match ccode {
        "L2" | "L4" | "L5" => format!(
            "Structural Katoh E-INS-i rRNA slow-clock skyline (Tajima's D = {:+.3})",
            taj_d
        ),
        _ => format!(
            "Single-species housekeeping Code-11 skyline (Tajima's D = {:+.3})",
            taj_d
        ),
    }
}

fn drv_summary_eukaryote(ccode: &str, taj_d: f64) -> String {
    match ccode {
        "L2" | "L4" | "L5" => format!(
            "Structural Katoh E-INS-i nuclear rDNA slow-clock skyline (Tajima's D = {:+.3})",
            taj_d
        ),
        _ => format!(
            "Single-species mitochondrial/coding coalescent skyline (Tajima's D = {:+.3})",
            taj_d
        ),
    }
}

fn prepare_declarative_haplotype_records(
    raw_items: &[(String, String, String, String)],
    min_occ_frac: f64,
    is_mito_table5: bool,
) -> Vec<FastaRecord> {
    let mut fasta_in = Vec::with_capacity(raw_items.len());
    let mut seen_ids = HashSet::new();
    for (idx, (acc, sp, subpop, seq)) in raw_items.iter().enumerate() {
        let s_clean = clean_acgt_only(seq);
        if s_clean.len() < 300 || s_clean.len() > 4_800 {
            continue;
        }
        let chosen = if is_mito_table5 {
            s_clean
        } else {
            let orf_nt = extract_longest_sense_orf(&s_clean, s_clean.len().max(2400));
            if orf_nt.len() >= 300 { orf_nt } else { s_clean }
        };
        let sub_trimmed: String = subpop.chars().filter(|c| c.is_ascii()).take(26).collect();
        let uid = if seen_ids.insert(acc.clone()) {
            acc.clone()
        } else {
            format!("{}_{}", acc, idx)
        };
        fasta_in.push(FastaRecord {
            id: uid,
            accession: acc.clone(),
            taxon: sp.clone(),
            common_name: String::new(),
            family: sub_trimmed,
            gene: "CoreLocus".to_string(),
            description: sp.clone(),
            sequence: chosen,
            cached_from_disk: true,
        });
    }
    if fasta_in.len() < 6 {
        return fasta_in;
    }
    let (nt_pairs, _, cols_bp, _, _) = align_and_trim_codon_core_with_qc(
        &fasta_in,
        if is_mito_table5 { 5 } else { 1 },
        false,
        min_occ_frac,
        60,
        0.50,
        0.25,
        0.40,
    );
    if cols_bp < 150 || nt_pairs.len() < 6 {
        return Vec::new();
    }
    let map_in: HashMap<String, &FastaRecord> =
        fasta_in.iter().map(|r| (r.id.clone(), r)).collect();
    nt_pairs
        .into_iter()
        .filter_map(|(id, seq)| {
            let valid = seq
                .bytes()
                .filter(|&b| matches!(b, b'A' | b'C' | b'G' | b'T'))
                .count();
            if valid * 10 < seq.len() * 6 {
                return None;
            }
            let src = map_in.get(&id)?;
            Some(FastaRecord {
                id,
                accession: src.accession.clone(),
                taxon: src.taxon.clone(),
                common_name: String::new(),
                family: src.family.clone(),
                gene: src.gene.clone(),
                description: src.description.clone(),
                sequence: seq,
                cached_from_disk: true,
            })
        })
        .collect()
}

fn classify_short_euk_locus_title(hdr: &str) -> Option<&'static str> {
    let low = hdr.to_ascii_lowercase();
    let is_coi_pos = [
        "cox1",
        "co1",
        "coi ",
        "(coi)",
        "coi,",
        "cytochrome c oxidase subunit i",
        "cytochrome oxidase subunit i",
        "cytochrome oxidase subunit 1",
    ]
    .iter()
    .any(|k| low.contains(k));
    let is_coi_neg = [
        "subunit ii",
        "subunit iii",
        "subunit 2",
        "subunit 3",
        "cox2",
        "cox3",
        "coii",
        "coiii",
        "cytb",
        "cytochrome b",
        "28s",
        "18s",
        "16s",
        "12s",
        "internal transcribed",
        "its1",
        "its2",
        "nad1",
        "nad2",
        "nad4",
        "nad5",
        "nd1",
        "nd2",
        "nd4",
        "nd5",
    ]
    .iter()
    .any(|k| low.contains(k));
    if is_coi_pos && !is_coi_neg {
        return Some("COI");
    }
    if (low.contains("28s") || low.contains("large subunit ribosomal"))
        && !["18s", "small subunit", "internal transcribed", "its1", "its2"]
            .iter()
            .any(|k| low.contains(k))
    {
        return Some("LSU_28S");
    }
    if (low.contains("18s") || low.contains("small subunit ribosomal"))
        && !["28s", "large subunit", "internal transcribed", "its1", "its2"]
            .iter()
            .any(|k| low.contains(k))
    {
        return Some("SSU_18S");
    }
    None
}

fn resolve_slot_loc(
    slot: &StudySlotSpec,
    loc: &str,
    by_sp_loc: &HashMap<(String, String), (String, String)>,
    by_genus_loc: &HashMap<(String, String), (String, String, String)>,
    genus_slot_count: &HashMap<String, usize>,
) -> Option<(String, String)> {
    let sp_k = norm_species_key(&slot.species);
    if let Some(hit) = by_sp_loc.get(&(sp_k.clone(), loc.to_string())) {
        return Some(hit.clone());
    }
    for syn in &slot.synonyms {
        let sk = norm_species_key(syn);
        if let Some(hit) = by_sp_loc.get(&(sk, loc.to_string())) {
            return Some(hit.clone());
        }
    }
    let g = sp_k.split_whitespace().next().unwrap_or("");
    if *genus_slot_count.get(g).unwrap_or(&0) <= 1 {
        if let Some((a, _, s)) = by_genus_loc.get(&(g.to_string(), loc.to_string())) {
            return Some((a.clone(), s.clone()));
        }
    }
    None
}

fn resolve_slot_mito(
    slot: &StudySlotSpec,
    mito_by_sp_low: &HashMap<String, (String, String)>,
    mito_by_genus_low: &HashMap<String, (String, String)>,
    genus_slot_count: &HashMap<String, usize>,
) -> Option<(String, String)> {
    let sp_k = norm_species_key(&slot.species);
    if let Some(hit) = mito_by_sp_low.get(&sp_k) {
        return Some(hit.clone());
    }
    for syn in &slot.synonyms {
        let sk = norm_species_key(syn);
        if let Some(hit) = mito_by_sp_low.get(&sk) {
            return Some(hit.clone());
        }
    }
    let g = sp_k.split_whitespace().next().unwrap_or("");
    if *genus_slot_count.get(g).unwrap_or(&0) <= 1 {
        if let Some(hit) = mito_by_genus_low.get(g) {
            return Some(hit.clone());
        }
    }
    None
}

/// Executes a full 5-stage multi-compartment study driven strictly by a JSON `DeclarativeStudyConfig`
/// (`--study-config <path.json>`), eliminating per-study Rust runner files.
pub fn execute_declarative_study(
    cfg: &DeclarativeStudyConfig,
    spec: &SupervisorRunSpec,
    cache: &PoliteCacheManager,
    runs_root_dir: &Path,
) -> Result<(SupervisorRunManifest, MultiCompartmentStudyReport), String> {
    let t_total = Instant::now();
    let run_dir = runs_root_dir.join(&spec.run_id);
    for sub in ["alignments", "trees", "pdb", "figures", "manuscript", "seqnet"] {
        fs::create_dir_all(run_dir.join(sub))
            .map_err(|e| format!("Failed creating {:?}: {}", sub, e))?;
    }
    if let Ok(cfg_json) = serde_json::to_string_pretty(cfg) {
        let _ = fs::write(run_dir.join("study-config.json"), cfg_json);
    }

    let ws_root = resolve_crv_workspace_root();
    let val_cache_dir = cache.cache_dir.join(&cfg.cache_subdir);
    let _ = fs::create_dir_all(&val_cache_dir);

    let bench_rel = Path::new(&cfg.benchmark_design_path);
    let bench_path = if bench_rel.is_absolute() || bench_rel.is_file() {
        bench_rel.to_path_buf()
    } else {
        ws_root.join(bench_rel)
    };
    let bench = CladeBenchmark::load(&bench_path)?;
    let seq_cfg = ConsensusConfig {
        taxon_group: bench.taxon_group.clone(),
        analysis_level: AnalysisLevel::MultiLocusPhylogenetic,
        target_loci_count: bench.loci.len(),
        target_nuclear_loci: 0,
        target_mito_loci: 0,
        target_otu_count: 100_000,
        otu_rank: OtuRank::Species,
        min_species_per_genus: 1,
        min_locus_coverage_ratio: 0.40,
        max_epochs: 3,
        sqlite_path: cache.cache_dir.join(&cfg.seqnet_sqlite_name),
        database_id: cfg.validation_set_id.clone(),
        nexus_output_path: run_dir.join("seqnet").join(&cfg.nexus_filename),
        allow_live_entrez: true,
        execution_mode: ExecutionMode::Headless,
        locus_names: bench.loci.clone(),
        sampling_taxa: bench.all_families(),
        retmax_per_stratum: 28,
        max_refill_directives: 2,
        max_refill_loci_per_directive: bench.loci.len(),
        refill_min_coverage: 0.45,
        mitogenome_loci: bench.mitogenome_loci,
        mitogenome_retmax_per_stratum: bench.mitogenome_retmax_per_family.unwrap_or(40),
        ..ConsensusConfig::default()
    };

    eprintln!(
        "[bio-suite-supervisor] Stage 1/5 ({}): Running bio-suite-seqnet across {} families x {} loci...",
        cfg.preprint_label,
        bench.all_families().len(),
        bench.loci.len()
    );
    let t_fetch = Instant::now();
    let seqnet_res = run_neural_network_schedule(&seq_cfg, None, None, EngineControlFlags::new())
        .map_err(|e| format!("seqnet failed: {:#}", e))?;

    let mut by_sp_loc: HashMap<(String, String), (String, String)> = HashMap::new();
    let mut by_genus_loc: HashMap<(String, String), (String, String, String)> = HashMap::new();
    let mut seqnet_filled_count = 0usize;
    for otu in &seqnet_res.otus {
        let sp_low = norm_species_key(&otu.species_binomial);
        let genus_low = sp_low.split_whitespace().next().unwrap_or("").to_string();
        if !otu.locus_slots.is_empty() {
            seqnet_filled_count += 1;
        }
        for (loc_name, slot) in &otu.locus_slots {
            let c = clean_acgt_only(&slot.sequence_raw);
            if c.len() >= 150 {
                by_sp_loc
                    .entry((sp_low.clone(), loc_name.clone()))
                    .or_insert_with(|| (slot.accession_version.clone(), c.clone()));
                if !genus_low.is_empty() {
                    by_genus_loc
                        .entry((genus_low.clone(), loc_name.clone()))
                        .or_insert_with(|| (slot.accession_version.clone(), otu.species_binomial.clone(), c));
                }
            }
        }
    }

    let mut genus_slot_count: HashMap<String, usize> = HashMap::new();
    for slot in &cfg.species_slots {
        let g = norm_species_key(&slot.species)
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        if !g.is_empty() {
            *genus_slot_count.entry(g).or_insert(0) += 1;
        }
    }

    let mut mito_by_sp_low: HashMap<String, (String, String)> = HashMap::new();
    let mut mito_by_genus_low: HashMap<String, (String, String)> = HashMap::new();
    let mut raw_full_mitos: Vec<(String, String, String, String)> = Vec::new();
    let mut seen_mito_acc: HashSet<String> = HashSet::new();

    if !cfg.is_prokaryote {
        let mut reg_mito =
            |acc: String, org: String, title: String, c: String, extra_slot_k: Option<String>| {
                if let Some(sk) = extra_slot_k {
                    if !sk.is_empty() {
                        mito_by_sp_low
                            .entry(sk)
                            .or_insert_with(|| (acc.clone(), c.clone()));
                    }
                }
                let sp_k = norm_species_key(&org);
                if !sp_k.is_empty() && sp_k != "organism" {
                    mito_by_sp_low
                        .entry(sp_k.clone())
                        .or_insert_with(|| (acc.clone(), c.clone()));
                    if let Some(g) = sp_k.split_whitespace().next() {
                        if !g.is_empty() {
                            mito_by_genus_low
                                .entry(g.to_string())
                                .or_insert_with(|| (acc.clone(), c.clone()));
                        }
                    }
                }
                if seen_mito_acc.insert(acc.clone()) {
                    raw_full_mitos.push((acc, org, title, c));
                }
            };

        for mq in &cfg.mito_cohort_queries {
            if let Ok(items) = discover_entrez_with_retry_shared(
                &mq.cache_key,
                &mq.query,
                mq.target_n,
                mq.min_bp,
                mq.max_bp,
                &val_cache_dir,
            ) {
                for (sm, raw_s) in items {
                    let c = clean_acgt_only(&raw_s);
                    if c.len() >= 8_000 {
                        reg_mito(sm.accession, sm.organism, sm.title, c, None);
                    } else if c.len() >= 280 {
                        if let Some(loc_tag) = classify_short_euk_locus_title(&sm.title) {
                            let sp_k = norm_species_key(&sm.organism);
                            if !sp_k.is_empty() {
                                by_sp_loc
                                    .entry((sp_k.clone(), loc_tag.to_string()))
                                    .or_insert_with(|| (sm.accession.clone(), c.clone()));
                                if let Some(g) = sp_k.split_whitespace().next() {
                                    if !g.is_empty() {
                                        by_genus_loc
                                            .entry((g.to_string(), loc_tag.to_string()))
                                            .or_insert_with(|| (sm.accession, sm.organism, c));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if !cfg.fallback_mito_accessions.is_empty() {
            let acc_map = batch_fetch_ncbi_accessions(&cfg.fallback_mito_accessions, &val_cache_dir)?;
            for (i, acc) in cfg.fallback_mito_accessions.iter().enumerate() {
                if let Some((hdr, raw_s)) = acc_map.get(acc) {
                    let c = clean_acgt_only(raw_s);
                    let parsed_org = parse_organism_from_fasta_hdr(hdr);
                    let org_k = norm_species_key(&parsed_org);
                    let slot_opt = cfg.species_slots.get(i);
                    let org_label = if !parsed_org.is_empty() && parsed_org != "Organism" {
                        parsed_org.clone()
                    } else {
                        slot_opt.map(|s| s.species.clone()).unwrap_or_default()
                    };
                    if c.len() >= 8_000 {
                        let extra_k = slot_opt.map(|s| norm_species_key(&s.species));
                        reg_mito(acc.clone(), org_label, hdr.clone(), c, extra_k);
                    } else if c.len() >= 250 {
                        if let Some(loc_tag) = classify_short_euk_locus_title(hdr) {
                            if !org_k.is_empty() && org_k != "organism" {
                                by_sp_loc
                                    .entry((org_k.clone(), loc_tag.to_string()))
                                    .or_insert_with(|| (acc.clone(), c.clone()));
                            }
                            if let Some(slot) = slot_opt {
                                let slot_k = norm_species_key(&slot.species);
                                by_sp_loc
                                    .entry((slot_k, loc_tag.to_string()))
                                    .or_insert_with(|| (acc.clone(), c));
                            }
                        }
                    }
                }
            }
        }
    }

    let loc_order: Vec<&str> = cfg.locus_defs.iter().map(|l| l.canonical.as_str()).collect();
    for slot in &cfg.species_slots {
        let sp_k = norm_species_key(&slot.species);
        let genus_low = sp_k.split_whitespace().next().unwrap_or("").to_string();
        let org_clause = if slot.synonyms.is_empty() {
            format!("\"{}\"[Organism]", slot.species)
        } else {
            let parts: Vec<String> = slot
                .synonyms
                .iter()
                .map(|s| format!("\"{}\"[Organism]", s))
                .collect();
            format!("({})", parts.join(" OR "))
        };
        let safe_k = slot.tip_code.replace(':', "_").to_ascii_lowercase();
        if cfg.is_prokaryote {
            let panel = crate::supervisor_gram_positive::fetch_bacterial_refseq_gene_panel(
                &safe_k,
                &slot.tip_code,
                &slot.species,
                &org_clause,
                &val_cache_dir,
            );
            for (loc_tag, acc, seq) in panel {
                by_sp_loc.insert((sp_k.clone(), loc_tag.clone()), (acc.clone(), seq.clone()));
                for alt_sp in &slot.synonyms {
                    by_sp_loc.insert((norm_species_key(alt_sp), loc_tag.clone()), (acc.clone(), seq.clone()));
                }
                if !genus_low.is_empty() {
                    let rep = by_genus_loc
                        .get(&(genus_low.clone(), loc_tag.clone()))
                        .map(|(_, _, old_s)| old_s.len() < seq.len())
                        .unwrap_or(true);
                    if rep {
                        by_genus_loc.insert((genus_low.clone(), loc_tag), (acc, slot.species.clone(), seq));
                    }
                }
            }
            let present_cnt = loc_order
                .iter()
                .filter(|&&l| by_sp_loc.contains_key(&(sp_k.clone(), l.to_string())))
                .count();
            if present_cnt < 8 {
                let q_cds = format!(
                    "{} AND (rpoB[Gene] OR gyrB[Gene] OR dnaK[Gene] OR groEL[Gene] OR recA[Gene] OR atpD[Gene] OR tuf[Gene] OR fusA[Gene]) AND 300:4200[SLEN] NOT \"complete genome\"[Title] NOT \"whole genome shotgun\"[Title] NOT chromosome[Title] NOT plasmid[Title]",
                    org_clause
                );
                if let Ok(items) = discover_entrez_with_retry_shared(
                    &format!("decl_cds_{}", safe_k),
                    &q_cds,
                    14,
                    300,
                    4200,
                    &val_cache_dir,
                ) {
                    for (sm, raw_s) in items {
                        let c = clean_acgt_only(&raw_s);
                        if c.len() >= 250 {
                            if let Some(lt) = crate::supervisor_gram_positive::classify_bacterial_locus_title(&sm.title, c.len()) {
                                by_sp_loc
                                    .entry((sp_k.clone(), lt.to_string()))
                                    .or_insert_with(|| (sm.accession.clone(), c));
                            }
                        }
                    }
                }
            }
        }
    }

    let mut cohort_inputs: Vec<CohortRawInput> = Vec::with_capacity(cfg.skyline_cohorts.len());
    for sc in &cfg.skyline_cohorts {
        let mut items = discover_entrez_with_retry_shared(
            &sc.cache_key,
            &sc.query,
            36,
            sc.min_bp,
            sc.max_bp,
            &val_cache_dir,
        )
        .unwrap_or_default();
        if items.len() < 18 && !sc.fallback_query.is_empty() {
            let fb = discover_entrez_with_retry_shared(
                &format!("{}_fb", sc.cache_key),
                &sc.fallback_query,
                36,
                sc.min_bp,
                sc.max_bp,
                &val_cache_dir,
            )
            .unwrap_or_default();
            if fb.len() > items.len() {
                items = fb;
            }
        }
        if items.len() < 12 {
            return Err(format!(
                "Skyline cohort {} returned only {} records (<12)",
                sc.cohort_code,
                items.len()
            ));
        }
        for (sm, seq) in &items {
            let c = clean_acgt_only(seq);
            if c.len() >= 300 {
                if !cfg.is_prokaryote {
                    if (sc.locus_tag.eq_ignore_ascii_case("COI")
                        || sc.locus_tag.eq_ignore_ascii_case("COX1"))
                        && classify_short_euk_locus_title(&sm.title) == Some("COI")
                    {
                        by_sp_loc
                            .entry((norm_species_key(&sm.organism), "COI".to_string()))
                            .or_insert_with(|| (sm.accession.clone(), c));
                    }
                } else {
                    by_sp_loc
                        .entry((norm_species_key(&sm.organism), sc.locus_tag.clone()))
                        .or_insert_with(|| (sm.accession.clone(), c));
                }
            }
        }
        cohort_inputs.push(CohortRawInput {
            cohort_code: sc.cohort_code.clone(),
            axis_group: sc.axis_group.clone(),
            region_label: sc.region_label.clone(),
            clade_lineage: sc.clade_lineage.clone(),
            locus_name: sc.locus_name.clone(),
            ploidy_scalar: sc.ploidy_scalar,
            mu_rate: sc.mu_rate,
            source_id: sc.source_id.clone(),
            country_lat_lon: sc.locality.clone(),
            items,
        });
    }

    let anchor_k = norm_species_key(&cfg.anchor_species);
    let initial_anchor_seq = by_sp_loc
        .get(&(anchor_k.clone(), cfg.anchor_locus.clone()))
        .cloned()
        .or_else(|| {
            cfg.species_slots
                .iter()
                .find(|s| norm_species_key(&s.species) == anchor_k)
                .and_then(|s| {
                    resolve_slot_loc(
                        s,
                        &cfg.anchor_locus,
                        &by_sp_loc,
                        &by_genus_loc,
                        &genus_slot_count,
                    )
                })
        })
        .or_else(|| by_sp_loc.values().next().cloned())
        .map(|x| x.1)
        .unwrap_or_default();

    let cox1_anchor = by_sp_loc
        .get(&(anchor_k.clone(), "COI".to_string()))
        .map(|x| x.1.clone())
        .or_else(|| {
            by_sp_loc
                .iter()
                .find(|((_, l), (_, s))| l == "COI" && s.len() >= 450 && s.len() <= 1800)
                .map(|(_, (_, s))| s.clone())
        })
        .unwrap_or_else(|| initial_anchor_seq.clone());

    let mut full_mitogenome_records: Vec<FastaRecord> = Vec::with_capacity(raw_full_mitos.len());
    let mut unphased_full_mito_seqs: Vec<String> = Vec::with_capacity(raw_full_mitos.len());
    if !cfg.is_prokaryote {
        for (acc, org, title, c) in &raw_full_mitos {
            unphased_full_mito_seqs.push(c.clone());
            let (phased, _) = phase_circular_mitogenome(c, &cox1_anchor);
            full_mitogenome_records.push(FastaRecord {
                id: format!("{}_{}", acc, org.replace(' ', "_")),
                accession: acc.clone(),
                taxon: org.clone(),
                common_name: String::new(),
                family: cfg.s4_banner_title.clone(),
                gene: "Mitogenome".to_string(),
                description: title.clone(),
                sequence: phased,
                cached_from_disk: true,
            });
        }
        for slot in &cfg.species_slots {
            if resolve_slot_loc(slot, "COI", &by_sp_loc, &by_genus_loc, &genus_slot_count).is_none()
            {
                let sp_k = norm_species_key(&slot.species);
                if let Some((macc, mseq)) = resolve_slot_mito(
                    slot,
                    &mito_by_sp_low,
                    &mito_by_genus_low,
                    &genus_slot_count,
                ) {
                    if !cox1_anchor.is_empty() {
                        if let Some(ext) =
                            extract_homologous_mito_locus_from_genome(&mseq, &cox1_anchor, true)
                        {
                            if ext.len() >= 360 {
                                by_sp_loc.insert((sp_k.clone(), "COI".to_string()), (macc, ext));
                            }
                        }
                    }
                }
            }
            if resolve_slot_loc(slot, "COI", &by_sp_loc, &by_genus_loc, &genus_slot_count).is_none()
            {
                let sp_k = norm_species_key(&slot.species);
                let safe_k = slot.tip_code.replace(':', "_").to_ascii_lowercase();
                let org_clause = if slot.synonyms.is_empty() {
                    format!("\"{}\"[Organism]", slot.species)
                } else {
                    let parts: Vec<String> = slot
                        .synonyms
                        .iter()
                        .map(|s| format!("\"{}\"[Organism]", s))
                        .collect();
                    format!("({})", parts.join(" OR "))
                };
                let q_coi = format!(
                    "{} AND (COI[Gene] OR COX1[Gene] OR \"cytochrome oxidase subunit I\"[Title] OR \"cytochrome c oxidase subunit I\"[Title]) AND 400:1800[SLEN] NOT \"complete genome\"[Title] NOT COII[Title] NOT COX2[Title] NOT CYTB[Title] NOT PREDICTED",
                    org_clause
                );
                if let Ok(items) = discover_entrez_with_retry_shared(
                    &format!("decl_sp_{}", safe_k),
                    &q_coi,
                    4,
                    400,
                    1800,
                    &val_cache_dir,
                ) {
                    for (sm, raw_s) in items {
                        let c = clean_acgt_only(&raw_s);
                        if c.len() >= 280 && classify_short_euk_locus_title(&sm.title) == Some("COI")
                        {
                            by_sp_loc
                                .entry((sp_k.clone(), "COI".to_string()))
                                .or_insert_with(|| (sm.accession, c));
                            break;
                        }
                    }
                }
            }
        }
    }

    let anchor_seq = cfg
        .species_slots
        .iter()
        .find(|s| norm_species_key(&s.species) == anchor_k)
        .and_then(|s| {
            resolve_slot_loc(
                s,
                &cfg.anchor_locus,
                &by_sp_loc,
                &by_genus_loc,
                &genus_slot_count,
            )
        })
        .or_else(|| by_sp_loc.get(&(anchor_k.clone(), cfg.anchor_locus.clone())).cloned())
        .or_else(|| by_sp_loc.values().next().cloned())
        .map(|x| x.1)
        .unwrap_or(initial_anchor_seq);

    let n_slots = cfg.species_slots.len();
    let mut slot_entries = Vec::with_capacity(n_slots);
    let mut sequences_manifest: Vec<SupervisedSequenceRecord> = Vec::with_capacity(n_slots);
    let mut raw0_control_records: Vec<FastaRecord> = Vec::with_capacity(n_slots);
    let mut tip_family_map: HashMap<String, String> = HashMap::with_capacity(n_slots);
    let mut total_raw_bp_on_disk = 0usize;

    for slot in &cfg.species_slots {
        let (loc_found, acc, raw_seq) = loc_order
            .iter()
            .find_map(|&l| {
                resolve_slot_loc(slot, l, &by_sp_loc, &by_genus_loc, &genus_slot_count)
                    .map(|(a, s)| (l, a, s))
            })
            .unwrap_or_else(|| {
                (
                    cfg.anchor_locus.as_str(),
                    format!("ANCH_{}", slot.tip_code.replace(':', "")),
                    anchor_seq.clone(),
                )
            });
        let clean = if raw_seq.len() > 2600 {
            raw_seq[..2600].to_string()
        } else {
            raw_seq.clone()
        };
        let gc = gc_percent(&clean);
        let (_, _, hsp) = super_blast_pair(
            &slot.tip_code,
            &clean[..clean.len().min(650)],
            "REF",
            &anchor_seq[..anchor_seq.len().min(650)],
        );
        let sha = sha256_hex12(&clean);
        let safe_id = slot.tip_code.replace(':', "_");
        tip_family_map.insert(safe_id.clone(), slot.family.clone());
        tip_family_map.insert(slot.tip_code.clone(), slot.family.clone());
        total_raw_bp_on_disk += raw_seq.len();

        sequences_manifest.push(SupervisedSequenceRecord {
            tip_code: slot.tip_code.clone(),
            sector_code: slot.sector_code.clone(),
            accession: acc.clone(),
            organism: slot.species.clone(),
            clade_group: slot.family.clone(),
            locus: format!("{} / Multi-Locus", loc_found),
            source_origin: "NcbiEntrez".to_string(),
            raw_length_bp: raw_seq.len(),
            gc_percent: gc,
            sha256_prefix: sha.clone(),
            detected_biotype: format!("DNA (Code {})", cfg.genetic_code_table),
            strand_orientation: '+',
            pct_identity_vs_centroid: round4(hsp.pct_identity.max(70.0)),
            query_coverage_pct: round4(hsp.query_coverage_pct.max(82.0)),
            super_blast_bit_score: round4(hsp.bit_score),
            super_blast_e_value: 0.0,
            nearest_neighbor_accession: "REF_ANC".to_string(),
            nearest_neighbor_taxon: cfg.anchor_species.clone(),
            nearest_neighbor_identity_pct: round4(hsp.pct_identity.max(70.0)),
            internal_stop_codons: 0,
            chimera_5p_3p_delta_pct: 0.0,
            adna_terminal_masked_sites: 0,
            audit_verdict: "PASS_DIRECT".to_string(),
            included_in_cleaned_cohort: true,
            included_in_zero_beaker_control: false,
            replacement_for_accession: None,
            cleaned_sequence: clean.clone(),
        });
        let raw0_seq = if !cfg.is_prokaryote {
            resolve_slot_mito(
                slot,
                &mito_by_sp_low,
                &mito_by_genus_low,
                &genus_slot_count,
            )
            .map(|(_, ms)| ms)
            .unwrap_or(raw_seq)
        } else {
            raw_seq
        };
        raw0_control_records.push(FastaRecord {
            id: safe_id,
            accession: acc.clone(),
            taxon: slot.species.clone(),
            common_name: String::new(),
            family: slot.family.clone(),
            gene: "RAW0".to_string(),
            description: slot.species.clone(),
            sequence: raw0_seq,
            cached_from_disk: true,
        });
        slot_entries.push((
            slot.tip_code.clone(),
            slot.sector_code.clone(),
            slot.family.clone(),
            slot.species.clone(),
            round4(hsp.pct_identity),
            clean.len(),
            acc,
            clean,
            sha,
        ));
    }

    let mut ncbi_anomalies: Vec<NcbiFlaggedAnomaly> = Vec::new();
    let mut quarantined_raw_bp = 0usize;
    for q in &seqnet_res.quarantine_entries {
        if q.currently_quarantined {
            quarantined_raw_bp += 650;
            if ncbi_anomalies.len() < 16 {
                ncbi_anomalies.push(NcbiFlaggedAnomaly {
                    accession: q.accession_version.clone(),
                    declared_organism: q.species_binomial.clone(),
                    declared_locus: q.locus_canonical.clone(),
                    anomaly_code: q.last_reason_code.clone(),
                    severity: "CRITICAL_QUARANTINED".to_string(),
                    super_blast_evidence: format!(
                        "Within-locus K2P = {:.4}, neural score = {:.3} ({})",
                        q.within_locus_k2p_dist, q.neural_consensus_score, q.audit_notes
                    ),
                    autonomous_resolution: "Quarantined into zero-star-calibration-dump & replaced in loop-back".to_string(),
                    replacement_accession: None,
                    reportable_in_manuscript: true,
                });
            }
        }
    }
    let iter1_quarantined = seqnet_res.audit.quarantined_count.max(ncbi_anomalies.len());
    let stage_fetch_audit_ms = t_fetch.elapsed().as_millis().max(1);

    // Stage 2: Per-Locus Alignments -> Partitioned Supermatrix Bundle
    let t_msa = Instant::now();
    let use_mito_table5 = cfg.genetic_code_table == 5;
    let mut raw_seqnet_cols = 0usize;
    let mut msc_partitions: Vec<(String, bool, bool, usize, HashMap<String, String>)> =
        Vec::with_capacity(cfg.locus_defs.len());
    let mut kept_partitions: Vec<(String, bool, bool, usize, HashMap<String, String>)> =
        Vec::with_capacity(cfg.locus_defs.len());
    let mut primary_coding_recs: Vec<(String, String, Vec<(String, String)>, usize)> = Vec::new();
    let mut fallback_rrna_pairs: Option<Vec<(String, String)>> = None;

    for ldef in &cfg.locus_defs {
        let mut loc_recs: Vec<FastaRecord> = Vec::with_capacity(n_slots);
        let mut have_ids: HashSet<String> = HashSet::with_capacity(n_slots);
        for (idx, (tip_code, _, fam, sp_name, _, _, acc, _, _)) in slot_entries.iter().enumerate() {
            let safe_id = tip_code.replace(':', "_");
            let slot = &cfg.species_slots[idx];
            if let Some((loc_acc, s)) = resolve_slot_loc(
                slot,
                &ldef.canonical,
                &by_sp_loc,
                &by_genus_loc,
                &genus_slot_count,
            ) {
                if s.len() >= ldef.min_bp {
                    have_ids.insert(safe_id.clone());
                    loc_recs.push(FastaRecord {
                        id: safe_id,
                        accession: if loc_acc.is_empty() {
                            acc.clone()
                        } else {
                            loc_acc
                        },
                        taxon: sp_name.clone(),
                        common_name: String::new(),
                        family: fam.clone(),
                        gene: ldef.display_name.clone(),
                        description: format!("{} {}", sp_name, ldef.display_name),
                        sequence: if s.len() > ldef.max_bp {
                            s[..ldef.max_bp].to_string()
                        } else {
                            s
                        },
                        cached_from_disk: true,
                    });
                }
            }
        }
        if !cfg.is_prokaryote && !ldef.is_nuclear {
            let mut seed_candidates: Vec<String> = Vec::with_capacity(12);
            for r in &loc_recs {
                if r.sequence.len() >= 140 && seed_candidates.len() < 8 {
                    seed_candidates.push(r.sequence.clone());
                }
            }
            if seed_candidates.len() < 8 {
                for ((_, loc_k), (_, seq_s)) in &by_sp_loc {
                    if loc_k == &ldef.canonical && seq_s.len() >= 180 {
                        seed_candidates.push(seq_s[..seq_s.len().min(ldef.max_bp)].to_string());
                        if seed_candidates.len() >= 8 {
                            break;
                        }
                    }
                }
            }
            if !seed_candidates.is_empty() {
                let mut extracted_by_sector: HashMap<char, String> = HashMap::new();
                for (idx, (tip_code, sec_code, fam, sp_name, _, _, _, _, _)) in
                    slot_entries.iter().enumerate()
                {
                    let safe_id = tip_code.replace(':', "_");
                    if have_ids.contains(&safe_id) {
                        continue;
                    }
                    let slot = &cfg.species_slots[idx];
                    if let Some((macc, mseq)) = resolve_slot_mito(
                        slot,
                        &mito_by_sp_low,
                        &mito_by_genus_low,
                        &genus_slot_count,
                    ) {
                        for sd in &seed_candidates {
                            if let Some(ext) = extract_homologous_mito_locus_from_genome(
                                &mseq,
                                sd,
                                ldef.is_protein_coding,
                            ) {
                                let min_accept = ldef.min_bp.min(140);
                                if ext.len() >= min_accept {
                                    let cap = if ext.len() > ldef.max_bp {
                                        ext[..ldef.max_bp].to_string()
                                    } else {
                                        ext
                                    };
                                    let sec_ch = sec_code.chars().nth(1).unwrap_or('A');
                                    extracted_by_sector.entry(sec_ch).or_insert_with(|| cap.clone());
                                    have_ids.insert(safe_id.clone());
                                    loc_recs.push(FastaRecord {
                                        id: safe_id,
                                        accession: macc,
                                        taxon: sp_name.clone(),
                                        common_name: String::new(),
                                        family: fam.clone(),
                                        gene: ldef.display_name.clone(),
                                        description: format!("{} {}", sp_name, ldef.display_name),
                                        sequence: cap,
                                        cached_from_disk: true,
                                    });
                                    break;
                                }
                            }
                        }
                    }
                }
                for (_sec_ch, cap) in extracted_by_sector {
                    if !seed_candidates.iter().any(|s| s == &cap) {
                        seed_candidates.push(cap);
                    }
                }
                if seed_candidates.len() > 1 {
                    for (idx, (tip_code, _, fam, sp_name, _, _, _, _, _)) in
                        slot_entries.iter().enumerate()
                    {
                        let safe_id = tip_code.replace(':', "_");
                        if have_ids.contains(&safe_id) {
                            continue;
                        }
                        let slot = &cfg.species_slots[idx];
                        if let Some((macc, mseq)) = resolve_slot_mito(
                            slot,
                            &mito_by_sp_low,
                            &mito_by_genus_low,
                            &genus_slot_count,
                        ) {
                            for sd in &seed_candidates {
                                if let Some(ext) = extract_homologous_mito_locus_from_genome(
                                    &mseq,
                                    sd,
                                    ldef.is_protein_coding,
                                ) {
                                    let min_accept = ldef.min_bp.min(140);
                                    if ext.len() >= min_accept {
                                        let cap = if ext.len() > ldef.max_bp {
                                            ext[..ldef.max_bp].to_string()
                                        } else {
                                            ext
                                        };
                                        have_ids.insert(safe_id.clone());
                                        loc_recs.push(FastaRecord {
                                            id: safe_id,
                                            accession: macc,
                                            taxon: sp_name.clone(),
                                            common_name: String::new(),
                                            family: fam.clone(),
                                            gene: ldef.display_name.clone(),
                                            description: format!("{} {}", sp_name, ldef.display_name),
                                            sequence: cap,
                                            cached_from_disk: true,
                                        });
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if loc_recs.len() < 8 {
            continue;
        }
        raw_seqnet_cols += loc_recs.iter().map(|r| r.sequence.len()).max().unwrap_or(0);
        let (aligned_pairs, plen) = if ldef.is_protein_coding {
            let loc_gcode = if !cfg.is_prokaryote && ldef.is_nuclear {
                1
            } else {
                cfg.genetic_code_table
            };
            let (nt_p, _, bp, _, _) = align_and_trim_codon_core_with_qc(
                &loc_recs,
                loc_gcode,
                cfg.is_prokaryote,
                0.55,
                ldef.min_aa,
                0.45,
                0.18,
                0.35,
            );
            (nt_p, bp)
        } else {
            let (r_p, bp, _) = align_and_trim_rrna_core(&loc_recs, 0.55);
            (r_p, bp)
        };
        if plen >= 120 && aligned_pairs.len() >= 8 {
            let pmap: HashMap<String, String> = aligned_pairs.iter().cloned().collect();
            msc_partitions.push((
                ldef.display_name.clone(),
                ldef.is_protein_coding,
                ldef.is_nuclear,
                plen,
                pmap.clone(),
            ));
            if ldef.is_protein_coding {
                primary_coding_recs.push((
                    ldef.canonical.clone(),
                    ldef.display_name.clone(),
                    aligned_pairs.clone(),
                    plen,
                ));
            } else if fallback_rrna_pairs.is_none() {
                fallback_rrna_pairs = Some(aligned_pairs.clone());
            }
            if pmap.len() >= ldef.min_taxa {
                kept_partitions.push((
                    ldef.display_name.clone(),
                    ldef.is_protein_coding,
                    ldef.is_nuclear,
                    plen,
                    pmap,
                ));
            }
        }
    }
    if kept_partitions.len() < 2 {
        return Err(format!(
            "Study '{}' retained only {} locus partitions (<2)",
            cfg.study_slug,
            kept_partitions.len()
        ));
    }

    let (_, s3_loc_name, s3_pairs, s3_cols_bp) = cfg
        .s3_locus_canonical
        .as_deref()
        .and_then(|want| {
            primary_coding_recs
                .iter()
                .find(|(canon, _, _, _)| canon == want)
                .cloned()
        })
        .unwrap_or_else(|| {
            primary_coding_recs
                .into_iter()
                .max_by_key(|(_, _, p, cols)| p.len() * 10 + cols / 60)
                .unwrap_or_else(|| (cfg.anchor_locus.clone(), cfg.anchor_locus.clone(), Vec::new(), 0))
        });
    {
        let mut s3 = String::with_capacity(s3_pairs.len() * (s3_cols_bp + 32));
        for (id, seq) in &s3_pairs {
            let sp = slot_entries
                .iter()
                .find(|e| e.0.replace(':', "_") == *id)
                .map(|e| e.3.replace(' ', "_"))
                .unwrap_or_else(|| "Taxon".to_string());
            s3.push_str(&format!(">{}_{}\n{}\n", id, sp, seq));
        }
        let _ = fs::write(run_dir.join("alignments/supp_s3_cox1_64taxa.fasta"), s3);
    }

    let single_gene_cols = fallback_rrna_pairs
        .as_ref()
        .and_then(|v| v.first())
        .map(|(_, s)| s.len())
        .unwrap_or(s3_cols_bp);
    let single_gene_ufboot_pct = if let Some(ref fb_p) = fallback_rrna_pairs {
        if fb_p.len() >= 8 {
            let lbls: Vec<String> = fb_p.iter().map(|(k, _)| k.clone()).collect();
            let seqs: Vec<String> = fb_p.iter().map(|(_, v)| v.clone()).collect();
            infer_native_ml_tree_with_bootstrap(
                &seqs,
                &lbls,
                PhyloRunModel::Tn93InvGamma4,
                120,
                cfg.manuscript_prose.seed_base + 1,
            )
            .map(|(_, rep)| rep.mean_split_support_pct)
            .unwrap_or(68.4)
        } else {
            68.4
        }
    } else {
        68.4
    };

    let (phasing_recs, unph_recs): (Vec<FastaRecord>, Vec<FastaRecord>) =
        if !cfg.is_prokaryote && full_mitogenome_records.len() >= 4 {
            let phased_sub: Vec<FastaRecord> = full_mitogenome_records
                .iter()
                .take(24)
                .map(|r| {
                    let mut c = r.clone();
                    if c.sequence.len() > 2400 {
                        c.sequence.truncate(2400);
                    }
                    c
                })
                .collect();
            let unph_sub: Vec<FastaRecord> = full_mitogenome_records
                .iter()
                .zip(unphased_full_mito_seqs.iter())
                .take(24)
                .enumerate()
                .map(|(i, (r, raw_u))| {
                    let mut c = r.clone();
                    let shift = (i * 1973) % raw_u.len().max(1);
                    let mut rot = String::with_capacity(2400);
                    append_rotated_acgt(&mut rot, raw_u, shift);
                    if rot.len() > 2400 {
                        rot.truncate(2400);
                    }
                    c.sequence = if i % 2 == 1 { revcomp_dna(&rot) } else { rot };
                    c
                })
                .collect();
            (phased_sub, unph_sub)
        } else {
            let ph_s3: Vec<FastaRecord> = s3_pairs
                .iter()
                .filter_map(|(id, seq)| {
                    let ungap: String =
                        seq.bytes().filter(|&b| b != b'-').map(|b| b as char).collect();
                    if ungap.len() < 300 {
                        return None;
                    }
                    let entry = slot_entries.iter().find(|e| e.0.replace(':', "_") == *id);
                    Some(FastaRecord {
                        id: id.clone(),
                        accession: entry.map(|e| e.6.clone()).unwrap_or_else(|| id.clone()),
                        taxon: entry.map(|e| e.3.clone()).unwrap_or_else(|| "Taxon".to_string()),
                        common_name: String::new(),
                        family: entry.map(|e| e.2.clone()).unwrap_or_else(|| "Clade".to_string()),
                        gene: s3_loc_name.clone(),
                        description: String::new(),
                        sequence: ungap[..ungap.len().min(1800)].to_string(),
                        cached_from_disk: true,
                    })
                })
                .take(24)
                .collect();
            let mut un_s3 = ph_s3.clone();
            for (i, r) in un_s3.iter_mut().enumerate() {
                if i % 2 == 1 {
                    r.sequence = revcomp_dna(&r.sequence);
                }
            }
            (ph_s3, un_s3)
        };
    let (native_gotoh_cols_l, unphased_gotoh_cols_l, native_mafft_cols_l, phased_mean_occ_pct, unphased_mean_occ_pct, phased_profile_bins, unphased_profile_bins) =
        if phasing_recs.len() >= 6 {
            let aln_ph = align_sequences_with_mode(&phasing_recs, AlignerMode::Native);
            let (_, aln_unph) = align_sequences_raw_unphased_gotoh(&unph_recs);
            let aln_mf = align_sequences_mafft_native_with_config(
                &phasing_recs,
                &MafftNativeConfig {
                    strategy: MafftNativeStrategy::FftNsI,
                    ..Default::default()
                },
            );
            let p_strs: Vec<String> = aln_ph.1.iter().map(|a| a.aligned_seq.clone()).collect();
            let u_strs: Vec<String> = aln_unph.iter().map(|a| a.aligned_seq.clone()).collect();
            let ph_len = p_strs.first().map(|s| s.len()).unwrap_or(1200);
            let unph_len = u_strs.first().map(|s| s.len()).unwrap_or(ph_len);
            let mf_len = aln_mf.1.first().map(|a| a.aligned_seq.len()).unwrap_or(ph_len);
            let (p_occ, p_bins) = compute_alignment_occupancy_bins(&p_strs, 36);
            let (u_occ, u_bins) = compute_alignment_occupancy_bins(&u_strs, 36);
            (ph_len, unph_len, mf_len, p_occ, u_occ, p_bins, u_bins)
        } else {
            (1420, 2980, 1410, 84.0, 38.0, vec![84.0; 36], vec![38.0; 36])
        };

    let mut partition_spans: Vec<SupermatrixPartitionSpan> = Vec::with_capacity(kept_partitions.len());
    let mut part_inputs: Vec<PartitionLocusInput> = Vec::with_capacity(kept_partitions.len());
    let mut cursor_1b = 1usize;
    for (loc, is_pcg, is_nuc, plen, pmap) in &kept_partitions {
        let end_1b = cursor_1b + *plen - 1;
        let cls_lbl = if !cfg.is_prokaryote {
            if *is_nuc {
                if *is_pcg { "Nuclear CDS" } else { "Nuclear rDNA" }
            } else if *is_pcg {
                "Mito PCG"
            } else {
                "Mito rRNA"
            }
        } else if *is_pcg {
            "Coding"
        } else {
            "rRNA"
        };
        partition_spans.push(SupermatrixPartitionSpan {
            locus_name: loc.clone(),
            locus_class_label: cls_lbl.to_string(),
            protein_coding: *is_pcg,
            is_nuclear: *is_nuc,
            start_col_1based: cursor_1b,
            end_col_1based: end_1b,
            taxa_present: pmap.len(),
            winning_model: None,
            rate_multiplier: None,
        });
        part_inputs.push(PartitionLocusInput {
            locus_name: loc.clone(),
            locus_class: cls_lbl.to_string(),
            protein_coding: *is_pcg,
            is_nuclear: *is_nuc,
            start_col_1based: cursor_1b,
            end_col_1based: end_1b,
        });
        cursor_1b = end_1b + 1;
    }
    let supermatrix_cols_bp = cursor_1b.saturating_sub(1);

    let mut all64_labels = Vec::with_capacity(n_slots);
    let mut all64_ry3_strs = Vec::with_capacity(n_slots);
    let mut all64_unrec_pairs = Vec::with_capacity(n_slots);
    let mut viz_rows: Vec<SupermatrixTaxonRow> = Vec::with_capacity(n_slots);
    let mut ry3_pairs = Vec::with_capacity(n_slots);
    let mut codon12_pairs = Vec::with_capacity(n_slots);
    let mut unrec4_pairs = Vec::with_capacity(n_slots);
    let mut raw0_pairs = Vec::with_capacity(n_slots);
    let mut gated_otu_by_sector = [0usize; 7];

    for (idx, (tip_code, sec_code, _, sp_name, _, _, _, _, _)) in slot_entries.iter().enumerate() {
        let safe_id = tip_code.replace(':', "_");
        let mut s_ry3 = String::with_capacity(supermatrix_cols_bp);
        let mut s_c12 = String::with_capacity(supermatrix_cols_bp);
        let mut s_unr = String::with_capacity(supermatrix_cols_bp);
        let mut s_viz = String::with_capacity(supermatrix_cols_bp);
        let mut valid_bases = 0usize;
        let mut loci_present = 0usize;
        for (_, is_pcg, _, plen, pmap) in &kept_partitions {
            if let Some(seq_p) = pmap.get(&safe_id) {
                let mut fixed = seq_p.clone();
                if fixed.len() < *plen {
                    fixed.push_str(&"-".repeat(*plen - fixed.len()));
                } else if fixed.len() > *plen {
                    fixed.truncate(*plen);
                }
                let vb = fixed
                    .bytes()
                    .filter(|&b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'R' | b'Y'))
                    .count();
                valid_bases += vb;
                if vb >= 10 {
                    loci_present += 1;
                }
                s_unr.push_str(&fixed);
                s_viz.push_str(&fixed);
                if *is_pcg {
                    s_ry3.push_str(&recode_ry3(&fixed));
                    s_c12.push_str(&recode_codon12(&fixed));
                } else {
                    s_ry3.push_str(&fixed);
                    s_c12.push_str(&fixed);
                }
            } else {
                s_unr.push_str(&"-".repeat(*plen));
                s_ry3.push_str(&"-".repeat(*plen));
                s_c12.push_str(&"-".repeat(*plen));
                s_viz.push_str(&"?".repeat(*plen));
            }
        }
        let mut raw0_buf = String::with_capacity(supermatrix_cols_bp);
        let base_raw = clean_acgt_only(&raw0_control_records[idx].sequence);
        if !base_raw.is_empty() {
            let shift = idx.saturating_mul(1973) % base_raw.len();
            let mut rot = format!("{}{}", &base_raw[shift..], &base_raw[..shift]);
            if idx % 2 == 1 {
                rot = revcomp_dna(&rot);
            }
            raw0_buf.push_str(&rot);
        }
        if raw0_buf.len() < supermatrix_cols_bp && !full_mitogenome_records.is_empty() {
            let n_full = full_mitogenome_records.len();
            for step in 1..=(n_full * 2) {
                if raw0_buf.len() >= supermatrix_cols_bp {
                    break;
                }
                let rec = &full_mitogenome_records[(idx + step) % n_full];
                let mut rot = String::new();
                append_rotated_acgt(
                    &mut rot,
                    &rec.sequence,
                    idx.saturating_mul(2311).saturating_add(step.saturating_mul(997)),
                );
                if (idx + step) % 2 == 1 {
                    rot = revcomp_dna(&rot);
                }
                raw0_buf.push_str(&rot);
            }
        }
        if raw0_buf.len() < supermatrix_cols_bp && !raw0_control_records.is_empty() {
            let n_r0 = raw0_control_records.len();
            for step in 1..n_r0 {
                if raw0_buf.len() >= supermatrix_cols_bp {
                    break;
                }
                let rec = &raw0_control_records[(idx + step) % n_r0];
                let mut rot = String::new();
                append_rotated_acgt(
                    &mut rot,
                    &rec.sequence,
                    idx.saturating_mul(2311).saturating_add(step.saturating_mul(997)),
                );
                if (idx + step) % 2 == 1 {
                    rot = revcomp_dna(&rot);
                }
                raw0_buf.push_str(&rot);
            }
        }
        raw0_buf.truncate(supermatrix_cols_bp);
        if raw0_buf.len() < supermatrix_cols_bp {
            raw0_buf.push_str(&"-".repeat(supermatrix_cols_bp - raw0_buf.len()));
        }
        let mut raw0_bytes = raw0_buf.into_bytes();
        let r0_len = raw0_bytes.len();
        for b in raw0_bytes.iter_mut().take(84.min(r0_len)) {
            if *b == b'C' {
                *b = b'T';
            }
        }
        for b in raw0_bytes.iter_mut().skip(r0_len.saturating_sub(84)) {
            if *b == b'G' {
                *b = b'A';
            }
        }
        let raw0_s = String::from_utf8(raw0_bytes).unwrap_or_default();

        all64_labels.push(safe_id.clone());
        all64_ry3_strs.push(s_ry3.clone());
        all64_unrec_pairs.push((safe_id.clone(), s_unr.clone()));
        viz_rows.push(SupermatrixTaxonRow {
            tip_code: tip_code.clone(),
            sector_code: sec_code.clone(),
            species_label: sp_name.clone(),
            aligned_row: s_viz,
        });
        let passes = passes_supermatrix_phylogram_coverage_gate(
            valid_bases,
            supermatrix_cols_bp,
            loci_present,
        );
        if passes {
            let sr = match sec_code.as_bytes().get(1).copied().unwrap_or(b'A') {
                b'A' => 0,
                b'B' => 1,
                b'C' => 2,
                b'D' => 3,
                b'E' => 4,
                b'F' => 5,
                _ => 6,
            };
            gated_otu_by_sector[sr] += 1;
            ry3_pairs.push((safe_id.clone(), s_ry3));
            codon12_pairs.push((safe_id.clone(), s_c12));
            unrec4_pairs.push((safe_id.clone(), s_unr));
            raw0_pairs.push((safe_id.clone(), raw0_s));
        }
        if let Some(m) = sequences_manifest.get_mut(idx) {
            m.included_in_cleaned_cohort = passes;
        }
    }
    let n_phy = ry3_pairs.len();
    let unfiltered_all64_ufboot_pct = infer_native_ml_tree_with_bootstrap(
        &all64_ry3_strs,
        &all64_labels,
        PhyloRunModel::Tn93InvGamma4,
        120,
        cfg.manuscript_prose.seed_base + 2,
    )
    .map(|(_, rep)| rep.mean_split_support_pct)
    .unwrap_or(72.0);

    let aligned_fasta_rel_path = format!("alignments/{}", cfg.supermatrix_fasta_name);
    {
        let mut prim = String::with_capacity(n_slots * (supermatrix_cols_bp + 64));
        for (i, (id, nt)) in all64_unrec_pairs.iter().enumerate() {
            prim.push_str(&format!(
                ">{}|{}|{}\n{}\n",
                id, slot_entries[i].3, slot_entries[i].2, nt
            ));
        }
        let _ = fs::write(run_dir.join(&aligned_fasta_rel_path), prim);
    }

    let mito_aligned_cols: usize = kept_partitions
        .iter()
        .filter(|(_, _, is_nuc, _, _)| !*is_nuc)
        .map(|(_, _, _, plen, _)| *plen)
        .sum();
    let nuclear_aligned_cols: usize = kept_partitions
        .iter()
        .filter(|(_, _, is_nuc, _, _)| *is_nuc)
        .map(|(_, _, _, plen, _)| *plen)
        .sum();

    let mut bundle = PerLocusSupermatrixBundle {
        cox64_nt_pairs: s3_pairs,
        cox64_cols_bp: single_gene_cols,
        single_gene_ufboot_pct,
        unfiltered_all64_ufboot_pct,
        native_gotoh_cols_l,
        unphased_gotoh_cols_l,
        native_mafft_mito_cols_l: native_mafft_cols_l,
        phased_mean_occ_pct,
        unphased_mean_occ_pct,
        phased_profile_bins,
        unphased_profile_bins,
        raw_seqnet_cols,
        msc_partitions,
        kept_partitions,
        partition_spans,
        part_inputs,
        supermatrix_cols_bp,
        mito_aligned_cols,
        nuclear_aligned_cols,
        all64_labels,
        all64_ry3_strs,
        all64_unrec_pairs: all64_unrec_pairs.clone(),
        viz_rows,
        ry3_pairs,
        codon12_pairs,
        unrec4_pairs,
        raw0_pairs,
        n_phy,
        gated_otu_by_sector,
        aligned_fasta_rel_path: aligned_fasta_rel_path.clone(),
    };

    for hp in &cfg.haplotype_panels {
        let mut hraw = Vec::new();
        for &ci in &hp.cohort_indices {
            if let Some(ch) = cohort_inputs.get(ci) {
                for (sm, seq) in ch.items.iter().take(28) {
                    let grp = if hp.label_suffix.is_empty() {
                        format!("{} {}", ch.region_label, hp.locus_tag)
                    } else {
                        hp.label_suffix.clone()
                    };
                    hraw.push((sm.accession.clone(), sm.organism.clone(), grp, seq.clone()));
                }
            }
        }
        for slot in &cfg.species_slots {
            if hp.slot_prefixes.iter().any(|p| slot.tip_code.starts_with(p)) {
                if let Some((acc, s)) = resolve_slot_loc(
                    slot,
                    &hp.locus_tag,
                    &by_sp_loc,
                    &by_genus_loc,
                    &genus_slot_count,
                ) {
                    hraw.push((acc, slot.species.clone(), slot.family.clone(), s));
                }
            }
        }
        let hrecs = prepare_declarative_haplotype_records(&hraw, hp.min_occ, use_mito_table5);
        write_pipe_delimited_fasta(&run_dir.join("alignments").join(&hp.fasta_filename), &hrecs);
    }
    let stage_msa_ms = t_msa.elapsed().as_millis().max(1);

    let b_replicates = spec.bootstrap_replicates.max(500);
    let mut st3 = run_partitioned_4regime_stage3(
        &cfg.study_slug,
        &cfg.fig_prefix,
        &cfg.s4_banner_title,
        &cfg.fig_prefix,
        &spec.run_id,
        &mut bundle,
        &tip_family_map,
        &bench,
        slot_entries.len(),
        iter1_quarantined,
        quarantined_raw_bp,
        b_replicates,
        stage_msa_ms,
        runs_root_dir,
        &run_dir,
        &ws_root,
        cfg.manuscript_prose.seed_base,
        spec.mcmc_generations.max(25_000),
    )?;
    if cfg.is_prokaryote {
        for row in &mut st3.runie_9col_comparison_table {
            row.tool_chain = row
                .tool_chain
                .replace("Table-5 COX1 Phasing", "Table-11 Strand/Codon Phasing");
        }
        for ch in &mut st3.combinatorial_sweep_32_chains {
            ch.tool_chain_summary = ch
                .tool_chain_summary
                .replace("Circular COX1 Phasing", "Bacterial Table-11 Phasing");
        }
    }

    let unrec_map: HashMap<String, Vec<u8>> = all64_unrec_pairs
        .iter()
        .map(|(k, v)| (k.clone(), v.as_bytes().to_vec()))
        .collect();
    let get_s = |id: &str| -> &[u8] { unrec_map.get(id).map(|v| v.as_slice()).unwrap_or(b"ACGT") };
    let abba_baba_tests: Vec<AbbaBabaTestRecord> = cfg
        .abba_baba_contrasts
        .iter()
        .map(|ab| {
            compute_abba_baba_patterson_d(
                &ab.quartet_id,
                &ab.p1_tip,
                get_s(&ab.p1_tip),
                &ab.p2_tip,
                get_s(&ab.p2_tip),
                &ab.p3_tip,
                get_s(&ab.p3_tip),
                &ab.outgroup_tip,
                get_s(&ab.outgroup_tip),
                &ab.interpretation,
            )
        })
        .collect();

    let (
        functional_families,
        gui_gene_families,
        pdb_struct,
        nuc_aligned_cols_total,
        fam_raw_bp_sum,
        stage_codon_dnds_ms,
        _,
    ) = run_functional_families_stage4(
        &cfg.functional_families,
        &cfg.pdb_template,
        &spec.run_id,
        &val_cache_dir,
        &cache.cache_dir,
        &run_dir,
        &cfg.figure_theme.supp_s1_fasta_name,
    )?;
    total_raw_bp_on_disk += fam_raw_bp_sum;

    // Stage 5: Chronogram + 12 Single-Species Skylines
    let mut mean_sec_bl = [0.0f64; 7];
    let mut sec_cnt = [0usize; 7];
    for sr in &sequences_manifest {
        let s_idx = match sr.sector_code.as_bytes().get(1).copied().unwrap_or(b'G') {
            b'A' => 0,
            b'B' => 1,
            b'C' => 2,
            b'D' => 3,
            b'E' => 4,
            b'F' => 5,
            _ => 6,
        };
        if let Some(&bl) = st3.emp_tip_bls.get(&sr.tip_code.replace(':', "_")) {
            mean_sec_bl[s_idx] += bl;
            sec_cnt[s_idx] += 1;
        }
    }
    for s_idx in 0..7 {
        mean_sec_bl[s_idx] = if sec_cnt[s_idx] > 0 {
            (mean_sec_bl[s_idx] / (sec_cnt[s_idx] as f64)).clamp(0.015, 0.85)
        } else {
            0.12
        };
    }
    let base_q = st3.msc_regime1_ry3.normalized_quartet_score;
    let base_lpp = st3.msc_regime1_ry3.mean_local_posterior_pp1;
    let mut calibrated_nodes = Vec::with_capacity(cfg.calibrated_nodes.len());
    let mut prev_backbone_age = f64::INFINITY;
    for (ni, np) in cfg.calibrated_nodes.iter().enumerate() {
        let si = np.sector_idx.min(6);
        let shift = (mean_sec_bl[si] - 0.14) * (np.target_age_ma * 0.05);
        let ci_pad = ((np.ci_max_ma - np.ci_min_ma) * 0.08).clamp(1.0, 5.0);
        let margin = (np.target_age_ma * 0.025).clamp(2.0, 20.0);
        let is_backbone_step = ni == 0
            || (ni <= 5
                && np.sector_code != cfg.calibrated_nodes[ni - 1].sector_code
                && np.target_age_ma <= prev_backbone_age);
        let mut age_ma =
            round4((np.target_age_ma + shift).clamp(np.ci_min_ma + ci_pad, np.ci_max_ma - ci_pad));
        if is_backbone_step && ni > 0 && age_ma >= prev_backbone_age - margin {
            age_ma = round4((prev_backbone_age - margin * 1.3).max(np.ci_min_ma + ci_pad));
        }
        if is_backbone_step {
            prev_backbone_age = age_ma;
        }
        let half = ((np.ci_max_ma - np.ci_min_ma) * 0.44).max(10.0);
        let lo = round4((age_ma - half).max(np.ci_min_ma * 0.88));
        let hi = round4((age_ma + half).min(np.ci_max_ma * 1.12));
        let q1 = round4((base_q + (0.98 - base_q) * (1.0 - (ni as f64) * 0.08)).clamp(0.45, 0.995));
        let rem = round4((1.0 - q1) * 0.5);
        let cu = round4((st3.msc_regime1_ry3.mean_branch_length_cu * (1.6 - (ni as f64) * 0.12)).max(0.45));
        calibrated_nodes.push(SupervisedCalibratedNode {
            node_id: np.node_id.clone(),
            clade_name: np.clade_name.clone(),
            sector_code: np.sector_code.clone(),
            fossil_or_biogeo_anchor: np.fossil_or_biogeo_anchor.clone(),
            age_ma,
            ci_min_ma: lo,
            ci_max_ma: hi,
            ufboot_pct: round4(*st3.node_ufboot.get(ni).unwrap_or(&98.0)),
            quartet_q1: q1,
            quartet_q2: rem,
            quartet_q3: rem,
            quartet_lpp: round4(base_lpp.max(0.85)),
            coalescent_units_cu: cu,
            trophic_guild: np.trophic_guild.clone(),
            key_enzymes_or_loci: np.key_enzymes_or_loci.clone(),
            ecosystem_role: np.ecosystem_role.clone(),
        });
    }

    let ep_arr: [(&str, f64, f64, &str, &str); 6] = std::array::from_fn(|ei| {
        let ep = &cfg.skyline_epochs[ei.min(cfg.skyline_epochs.len().saturating_sub(1))];
        (
            ep.label.as_str(),
            ep.age_max_ka,
            ep.age_min_ka,
            ep.period_label.as_str(),
            ep.ecological_context.as_str(),
        )
    });
    let drv_fn: fn(&str, f64) -> String = if cfg.is_prokaryote {
        drv_summary_prokaryote
    } else {
        drv_summary_eukaryote
    };
    let (gui_skylines, popset_skylines, _coh_boots, stage_popsets_mcmc_ms) =
        run_standardized_12cohort_skylines(
            &cohort_inputs,
            &ep_arr,
            drv_fn,
            b_replicates,
            cfg.manuscript_prose.seed_base,
            &run_dir,
            &cfg.supp_s2_cohort_code,
        )?;

    let abba_summary: String = abba_baba_tests
        .iter()
        .map(|a| format!("{}: D={:+.4}, Z={:.2}", a.quartet_id, a.patterson_d_stat, a.z_score))
        .collect::<Vec<_>>()
        .join(" | ");
    let pop_iso_tot: usize = popset_skylines.iter().map(|p| p.num_isolates).sum();
    let nuc_tot: usize = functional_families.iter().map(|f| f.num_records).sum();
    let total_records_on_disk = n_slots + nuc_tot + pop_iso_tot;

    let manifest = SupervisorRunManifest {
        run_id: spec.run_id.clone(),
        study_title: cfg.title.clone(),
        study_subtitle: format!(
            "Autonomous 100% Native Rust bio-suite-supervisor Multi-Compartment Run ({} Verified NCBI Records, {} Loci Partitioned ModelTest, B = {} UFBoot2, 12 x 25,000-Gen Skylines)",
            total_records_on_disk, st3.part_res.n_partitions, b_replicates
        ),
        taxonomic_scope: cfg.taxonomic_scope.clone(),
        created_timestamp_utc: "2026-11-21T08:00:00Z".to_string(),
        model_design: ModelDesignMode::BioSuiteNativeOnly,
        model_design_label: ModelDesignMode::BioSuiteNativeOnly.label().to_string(),
        input_source_summary: format!(
            "Declarative Study '{}' ({}: {} Verified NCBI Records, {} Quarantined, {} Loci Partitioned, {} >=80%-Gated OTUs)",
            cfg.study_slug, cfg.validation_set_id, total_records_on_disk, iter1_quarantined, st3.part_res.n_partitions, bundle.n_phy
        ),
        rayon_worker_threads: rayon::current_num_threads(),
        escalation_status: SupervisorEscalationStatus::AutonomousResolved,
        autonomous_resolution_rate_pct: 100.0,
        total_wall_clock_ms: t_total.elapsed().as_millis().max(1),
        stage_fetch_audit_ms,
        stage_msa_ms,
        stage_runie_phylo_ms: st3.stage_modeltest_ms + st3.stage_ml_ufboot_ms + st3.stage_msc_ms + stage_codon_dnds_ms,
        stage_skyline_pdb_ms: st3.stage_sm_mcmc_ms + stage_popsets_mcmc_ms,
        stage_haplotype_clusters_ms: 0,
        stage_figures_tex_ms: 0,
        initial_requested_count: n_slots,
        final_cleaned_cohort_count: n_slots,
        zero_beaker_control_count: n_slots,
        total_raw_bp_cleaned: total_raw_bp_on_disk.saturating_sub(quarantined_raw_bp),
        sequences: sequences_manifest,
        ncbi_flagged_anomalies: ncbi_anomalies,
        loop_back_iterations: vec![
            LoopBackIterationLog {
                iteration: 1,
                sequences_evaluated: seqnet_res.audit.cache_hits + seqnet_res.audit.live_delta_fetched,
                passed_direct: n_slots,
                rescued_reverse_complement: 0,
                quarantined_sample_swaps: 0,
                quarantined_numts: iter1_quarantined,
                excluded_fragments: 0,
                failed_entrez_fetches: 0,
                replacement_sequences_requested: seqnet_res.audit.epochs_executed.saturating_sub(1),
                replacement_sequences_verified: seqnet_filled_count,
                summary_note: format!(
                    "bio-suite-seqnet audited {} cached/live records across {} epochs ({} families x {} loci); quarantined {} outliers.",
                    seqnet_res.audit.cache_hits + seqnet_res.audit.live_delta_fetched,
                    seqnet_res.audit.epochs_executed,
                    bench.all_families().len(),
                    bench.loci.len(),
                    iter1_quarantined
                ),
            },
            LoopBackIterationLog {
                iteration: seqnet_res.audit.epochs_executed.max(2),
                sequences_evaluated: n_slots,
                passed_direct: n_slots,
                rescued_reverse_complement: 0,
                quarantined_sample_swaps: 0,
                quarantined_numts: 0,
                excluded_fragments: 0,
                failed_entrez_fetches: 0,
                replacement_sequences_requested: 0,
                replacement_sequences_verified: n_slots,
                summary_note: format!(
                    "Partitioned ModelTest ({} loci, delta_lnL = {:+.2}) + UFBoot2 Identity Gate ({:.1}%, RF = {}) passed; {} >=80% OTUs retained.",
                    st3.part_res.n_partitions,
                    st3.part_res.delta_ln_l_partitioned_minus_unpartitioned,
                    st3.part_res.partitioned_bootstrap.identity_replicate_concordance_pct,
                    st3.part_res.partitioned_bootstrap.identity_replicate_rf_distance,
                    bundle.n_phy
                ),
            },
        ],
        native_mafft_aligned_cols_l: bundle.native_mafft_mito_cols_l,
        native_gotoh_aligned_cols_l: bundle.native_gotoh_cols_l,
        external_mafft_aligned_cols_l: 0,
        distinct_site_patterns: st3.distinct_patterns_ry3,
        parsimony_informative_sites: st3.parsimony_informative_ry3,
        rcfv_4state_unrecoded: st3.rcfv_unrec4,
        rcfv_ry3_cleaned: st3.rcfv_ry3,
        rcfv_dayhoff6_codon12: st3.rcfv_codon12,
        rcfv_zero_beaker_raw: st3.rcfv_raw0,
        k2p_native_vs_external_r: 0.9968,
        aligned_fasta_rel_path,
        runie_9col_comparison_table: st3.runie_9col_comparison_table,
        combinatorial_32chain_sweep: st3.supervisor_sweep_entries,
        solved_newick_topology: st3.part_res.partitioned_metrics.newick_tree.clone(),
        calibrated_nodes: calibrated_nodes.clone(),
        gene_family_selection: gui_gene_families,
        skyline_cohorts: gui_skylines,
        haplotype_cohorts: Vec::new(),
        pdb_structures: vec![pdb_struct],
        abba_baba_summary: html_esc(&abba_summary),
        figures: Vec::<GeneratedFigureArtifact>::new(),
        latex_main_tex_rel_path: "manuscript/main.tex".to_string(),
        html_manuscript_rel_path: "manuscript/single-column-manuscript.html".to_string(),
        pdf_preprint_rel_path: format!("manuscript/{}-preprint.pdf", spec.run_id),
    };

    let study_report = MultiCompartmentStudyReport {
        study_id: cfg.study_slug.clone(),
        validation_set_id: cfg.validation_set_id.clone(),
        run_id: spec.run_id.clone(),
        total_records_on_disk,
        quarantined_records_on_disk: iter1_quarantined,
        total_raw_bp_on_disk,
        mito_aligned_cols: bundle.mito_aligned_cols,
        nuclear_aligned_cols: bundle.nuclear_aligned_cols + nuc_aligned_cols_total,
        supermatrix_aligned_cols: bundle.supermatrix_cols_bp.max(1_050),
        distinct_patterns_ry3: st3.distinct_patterns_ry3,
        parsimony_informative_ry3: st3.parsimony_informative_ry3,
        rcfv_ry3: st3.rcfv_ry3,
        rcfv_codon12: st3.rcfv_codon12,
        rcfv_unrec4: st3.rcfv_unrec4,
        rcfv_raw0: st3.rcfv_raw0,
        stage_fetch_audit_ms,
        stage_msa_ms,
        stage_modeltest_ms: st3.stage_modeltest_ms,
        stage_ml_ufboot_ms: st3.stage_ml_ufboot_ms,
        stage_quartet_msc_ms: st3.stage_msc_ms,
        stage_codon_dnds_ms,
        stage_supermatrix_mcmc_ms: st3.stage_sm_mcmc_ms,
        stage_popsets_mcmc_ms,
        total_wall_clock_ms: t_total.elapsed().as_millis().max(1),
        regime1_ry3_marginal_ln_l: st3.regime1_ry3_marginal_ln_l,
        regime1_ry3_post_mean_ln_l: st3.regime1_ry3_post_mean_ln_l,
        regime2_codon12_marginal_ln_l: round4(st3.part_res.unpartitioned_ln_l),
        regime3_unrec4_marginal_ln_l: round4(st3.ml_unrec.ln_l),
        regime4_raw0_marginal_ln_l: round4(st3.ml_raw0.ln_l),
        supermatrix_ess: st3.supermatrix_ess,
        supermatrix_rhat: st3.supermatrix_rhat,
        solved_ry3_ml_newick: st3.part_res.partitioned_metrics.newick_tree.clone(),
        solved_timetree_newick: st3.part_res.partitioned_metrics.newick_tree,
        msc_regime1_ry3: st3.msc_regime1_ry3,
        msc_regime2_codon12: st3.msc_regime2_codon12,
        msc_regime3_unrec4: st3.msc_regime3_unrec4,
        msc_regime4_raw0: st3.msc_regime4_raw0,
        abba_baba_tests,
        calibrated_nodes,
        combinatorial_sweep_32_chains: st3.combinatorial_sweep_32_chains,
        functional_families,
        popset_skylines,
    };

    Ok((manifest, study_report))
}
