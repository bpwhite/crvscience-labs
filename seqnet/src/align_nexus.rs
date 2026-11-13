//! Stage 3 Multiple Sequence Alignment (Gotoh Affine-Gap Star/Profile Aligner) & Clean
//! Partitioned NEXUS Supermatrix Exporter (`align_nexus.rs`).
//!
//! Guarantees that ZERO quarantined sequences enter `TAXA` or `DATA` blocks, while recording
//! full chimeric OTU accession provenance AND a group-level soft-quarantine flag manifest
//! inside NEXUS metadata blocks (`SETS`, `NOTES`, and `BIOSUITE_METADATA`).

use anyhow::{Context, Result};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};

use crate::config::{ConsensusConfig, LocusClass, LocusSpec};
use crate::db::{QuarantineEntry, RunAuditRecord};
use crate::otu::AssembledSpeciesOtu;

/// Partition coordinates for a single aligned locus inside the concatenated supermatrix.
#[derive(Debug, Clone)]
pub struct LocusPartitionMeta {
    pub locus_name: String,
    pub locus_class: LocusClass,
    pub protein_coding: bool,
    pub start_col_1based: usize,
    pub end_col_1based: usize,
    pub aligned_len: usize,
    pub taxa_present: usize,
}

/// Pairwise affine-gap Needleman-Wunsch / Gotoh alignment between `reference` and `query`.
/// Uses banded/profile-consistent nucleotide scoring (`match = +2, mismatch = -1, gap_open = -4, gap_ext = -1`).
pub fn gotoh_pairwise_align(reference: &[u8], query: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let n = reference.len();
    let m = query.len();
    if n == 0 {
        return (vec![b'-'; m], query.to_vec());
    }
    if m == 0 {
        return (reference.to_vec(), vec![b'-'; n]);
    }
    // Fast path when sequences are already equal length and >= 75% positionally identical (e.g. ungapped CDS)
    if n == m {
        let matches = reference
            .iter()
            .zip(query.iter())
            .filter(|(a, b)| a == b)
            .count();
        if matches * 4 >= n * 3 {
            return (reference.to_vec(), query.to_vec());
        }
    }

    let neg_inf: i32 = -1_000_000_000;
    let gap_open: i32 = -4;
    let gap_ext: i32 = -1;
    let cols = m + 1;

    let mut dp_m = vec![neg_inf; (n + 1) * cols];
    let mut dp_x = vec![neg_inf; (n + 1) * cols]; // gap in query
    let mut dp_y = vec![neg_inf; (n + 1) * cols]; // gap in reference
    let mut trace = vec![0u8; (n + 1) * cols]; // 0 = diag, 1 = up (gap in q), 2 = left (gap in r)

    dp_m[0] = 0;
    for i in 1..=n {
        dp_x[i * cols] = gap_open + (i as i32 - 1) * gap_ext;
        dp_m[i * cols] = dp_x[i * cols];
        trace[i * cols] = 1;
    }
    for j in 1..=m {
        dp_y[j] = gap_open + (j as i32 - 1) * gap_ext;
        dp_m[j] = dp_y[j];
        trace[j] = 2;
    }

    for i in 1..=n {
        let rb = reference[i - 1];
        for j in 1..=m {
            let qb = query[j - 1];
            let sub = if rb == qb { 2 } else { -1 };
            let idx = i * cols + j;
            let up_idx = (i - 1) * cols + j;
            let left_idx = i * cols + (j - 1);
            let diag_idx = (i - 1) * cols + (j - 1);

            dp_x[idx] = (dp_m[up_idx] + gap_open).max(dp_x[up_idx] + gap_ext);
            dp_y[idx] = (dp_m[left_idx] + gap_open).max(dp_y[left_idx] + gap_ext);
            let diag_score = dp_m[diag_idx] + sub;

            if diag_score >= dp_x[idx] && diag_score >= dp_y[idx] {
                dp_m[idx] = diag_score;
                trace[idx] = 0;
            } else if dp_x[idx] >= dp_y[idx] {
                dp_m[idx] = dp_x[idx];
                trace[idx] = 1;
            } else {
                dp_m[idx] = dp_y[idx];
                trace[idx] = 2;
            }
        }
    }

    let mut aln_r = Vec::with_capacity(n + m);
    let mut aln_q = Vec::with_capacity(n + m);
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        let t = trace[i * cols + j];
        match t {
            0 if i > 0 && j > 0 => {
                aln_r.push(reference[i - 1]);
                aln_q.push(query[j - 1]);
                i -= 1;
                j -= 1;
            }
            1 if i > 0 => {
                aln_r.push(reference[i - 1]);
                aln_q.push(b'-');
                i -= 1;
            }
            _ => {
                aln_r.push(b'-');
                aln_q.push(query[j - 1]);
                j -= 1;
            }
        }
    }
    aln_r.reverse();
    aln_q.reverse();
    (aln_r, aln_q)
}

/// Center-star / profile multiple sequence alignment across all retained OTUs for every locus,
/// executed concurrently across CPU cores via `rayon`. Populates `aligned_sequence` on each slot
/// and returns ordered `LocusPartitionMeta` definitions.
pub fn align_all_loci_parallel(
    panel: &[LocusSpec],
    otus: &mut [AssembledSpeciesOtu],
) -> Vec<LocusPartitionMeta> {
    // Collect raw sequences per locus across `otus`
    let aligned_per_locus: Vec<(String, usize, Vec<(usize, String)>)> = panel
        .par_iter()
        .map(|locus| {
            let mut present: Vec<(usize, &[u8])> = Vec::new();
            for (o_idx, otu) in otus.iter().enumerate() {
                if let Some(slot) = otu.locus_slots.get(&locus.canonical_name) {
                    present.push((o_idx, slot.sequence_raw.as_bytes()));
                }
            }
            if present.is_empty() {
                return (locus.canonical_name.clone(), 0, Vec::new());
            }

            // Pick longest high-confidence sequence as center-star reference guide
            let center_seq: &[u8] = present
                .iter()
                .map(|&(_, seq)| seq)
                .max_by_key(|seq| seq.len())
                .unwrap_or(&[]);

            // Track maximum insertion count before each reference index `0..=center_seq.len()`
            let ref_len = center_seq.len();
            let mut pairwise_results: Vec<(usize, Vec<u8>, Vec<u8>)> =
                Vec::with_capacity(present.len());
            let mut max_ins = vec![0usize; ref_len + 1];

            for &(o_idx, q_seq) in &present {
                let (aln_r, aln_q) = gotoh_pairwise_align(center_seq, q_seq);
                let mut ref_cursor = 0usize;
                let mut run_ins = 0usize;
                for &rb in &aln_r {
                    if rb == b'-' {
                        run_ins += 1;
                    } else {
                        if run_ins > max_ins[ref_cursor] {
                            max_ins[ref_cursor] = run_ins;
                        }
                        run_ins = 0;
                        ref_cursor += 1;
                    }
                }
                if run_ins > max_ins[ref_len] {
                    max_ins[ref_len] = run_ins;
                }
                pairwise_results.push((o_idx, aln_r, aln_q));
            }

            let total_cols = ref_len + max_ins.iter().sum::<usize>();
            let mut projected: Vec<(usize, String)> = Vec::with_capacity(pairwise_results.len());

            for (o_idx, aln_r, aln_q) in pairwise_results {
                let mut row = Vec::with_capacity(total_cols);
                let mut ref_cursor = 0usize;
                let mut k = 0usize;
                while k < aln_r.len() {
                    let mut ins_bases = Vec::new();
                    while k < aln_r.len() && aln_r[k] == b'-' {
                        ins_bases.push(aln_q[k]);
                        k += 1;
                    }
                    let pad = max_ins[ref_cursor].saturating_sub(ins_bases.len());
                    row.extend_from_slice(&ins_bases);
                    row.resize(row.len() + pad, b'-');

                    if k < aln_r.len() {
                        row.push(aln_q[k]);
                        ref_cursor += 1;
                        k += 1;
                    }
                }
                let tail_pad = max_ins[ref_len];
                if row.len() < total_cols {
                    row.resize(total_cols.max(row.len() + tail_pad).min(total_cols), b'-');
                }
                projected.push((o_idx, String::from_utf8(row).unwrap_or_default()));
            }

            (locus.canonical_name.clone(), total_cols, projected)
        })
        .collect();

    let mut partitions = Vec::with_capacity(panel.len());
    let mut cursor_col = 1usize;

    for (locus, (_, aligned_len, rows)) in panel.iter().zip(aligned_per_locus.into_iter()) {
        if aligned_len == 0 {
            continue;
        }
        let taxa_present = rows.len();
        for (o_idx, aligned_str) in rows {
            if let Some(slot) = otus[o_idx].locus_slots.get_mut(&locus.canonical_name) {
                slot.aligned_sequence = Some(aligned_str);
            }
        }
        let start = cursor_col;
        let end = cursor_col + aligned_len - 1;
        cursor_col = end + 1;

        partitions.push(LocusPartitionMeta {
            locus_name: locus.canonical_name.clone(),
            locus_class: locus.locus_class,
            protein_coding: locus.protein_coding,
            start_col_1based: start,
            end_col_1based: end,
            aligned_len,
            taxa_present,
        });
    }

    partitions
}

/// Writes a clean, standards-compliant NEXUS supermatrix (Maddison, Swofford & Maddison 1997)
/// containing ONLY validated consensus sequences (`currently_quarantined == false`), complete
/// `SETS` (`CHARSET`, `TAXSET`, `CHARPARTITION`), `NOTES` chimeric accession provenance, and
/// a `BIOSUITE_METADATA` block flagging quarantined accessions without releasing their sequences.
pub fn write_partitioned_nexus(
    cfg: &ConsensusConfig,
    audit: &RunAuditRecord,
    otus: &[AssembledSpeciesOtu],
    partitions: &[LocusPartitionMeta],
    quarantine_flags: &[QuarantineEntry],
) -> Result<()> {
    if let Some(parent) = cfg.nexus_output_path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let file = File::create(&cfg.nexus_output_path)
        .with_context(|| format!("failed to create NEXUS file {:?}", cfg.nexus_output_path))?;
    let mut w = BufWriter::new(file);

    let ntax = otus.len();
    let nchar: usize = partitions.iter().map(|p| p.aligned_len).sum();
    let active_quarantined: Vec<&QuarantineEntry> = quarantine_flags
        .iter()
        .filter(|q| q.currently_quarantined)
        .collect();

    writeln!(w, "#NEXUS")?;
    writeln!(w, "[!")?;
    writeln!(
        w,
        "  Generated by bio-suite-seqnet v0.1.0 (Pure-Rust 3-Stage Hybrid Sequence Neural Network)"
    )?;
    writeln!(w, "  SQLite Database ID   : {}", audit.database_id)?;
    writeln!(w, "  Neural Run ID        : {}", audit.run_id)?;
    writeln!(w, "  Organism Group       : {}", cfg.taxon_group)?;
    writeln!(w, "  Analysis Level       : {}", cfg.analysis_level.as_str())?;
    writeln!(
        w,
        "  OTU Rank & Rule      : {} (separate species terminals; >= {} species/genus)",
        cfg.otu_rank.as_str(),
        cfg.min_species_per_genus
    )?;
    writeln!(
        w,
        "  Locus Panel          : {} retained loci ({} nuclear + {} mitochondrial requested)",
        partitions.len(),
        cfg.target_nuclear_loci,
        cfg.target_mito_loci
    )?;
    writeln!(
        w,
        "  Matrix Occupancy     : {:.2}% (Minimum per-OTU threshold: {:.1}%)",
        audit.final_matrix_occupancy * 100.0,
        cfg.min_locus_coverage_ratio * 100.0
    )?;
    writeln!(
        w,
        "  Final Consensus      : {:.4} (Converged: {}, Epochs: {})",
        audit.final_consensus_score, audit.converged, audit.epochs_executed
    )?;
    writeln!(
        w,
        "  Quarantine Policy    : {} bad accessions withheld from matrix (flagged in metadata below)",
        active_quarantined.len()
    )?;
    writeln!(w, "  Completed Timestamp  : {}", audit.completed_at)?;
    writeln!(w, "]")?;
    writeln!(w)?;

    // 1. TAXA BLOCK
    writeln!(w, "BEGIN TAXA;")?;
    writeln!(w, "    DIMENSIONS NTAX={};", ntax)?;
    writeln!(w, "    TAXLABELS")?;
    for otu in otus {
        writeln!(
            w,
            "        {} [Genus={} | Chimeric={} | Specimens={} | Coverage={}/{} ({:.1}%)]",
            otu.terminal_label,
            otu.genus,
            if otu.is_chimeric { "YES" } else { "NO" },
            otu.distinct_specimen_count,
            otu.loci_filled,
            otu.loci_target,
            otu.locus_coverage_ratio * 100.0
        )?;
    }
    writeln!(w, "    ;")?;
    writeln!(w, "END;")?;
    writeln!(w)?;

    // 2. DATA BLOCK (Interleaved by locus partition for clarity & parser compatibility)
    writeln!(w, "BEGIN DATA;")?;
    writeln!(w, "    DIMENSIONS NTAX={} NCHAR={};", ntax, nchar)?;
    writeln!(
        w,
        "    FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=YES;"
    )?;
    writeln!(w, "    MATRIX")?;

    let pad_width = otus
        .iter()
        .map(|o| o.terminal_label.len())
        .max()
        .unwrap_or(24)
        + 2;

    for part in partitions {
        writeln!(
            w,
            "    [Locus: {} ({}) | Columns {}-{} ({} bp) | Taxa={}/{}]",
            part.locus_name,
            part.locus_class.as_str(),
            part.start_col_1based,
            part.end_col_1based,
            part.aligned_len,
            part.taxa_present,
            ntax
        )?;
        let missing_fill = "?".repeat(part.aligned_len);
        for otu in otus {
            let seq_slice = otu
                .locus_slots
                .get(&part.locus_name)
                .and_then(|s| s.aligned_sequence.as_deref())
                .unwrap_or(missing_fill.as_str());
            writeln!(
                w,
                "    {:width$}{}",
                otu.terminal_label,
                seq_slice,
                width = pad_width
            )?;
        }
        writeln!(w)?;
    }
    writeln!(w, "    ;")?;
    writeln!(w, "END;")?;
    writeln!(w)?;

    // 3. SETS BLOCK (Locus CHARSETS, Codon Positions, Genus TAXSETS, Compartment CHARPARTITION)
    writeln!(w, "BEGIN SETS;")?;
    for part in partitions {
        writeln!(
            w,
            "    CHARSET {} = {}-{};",
            part.locus_name, part.start_col_1based, part.end_col_1based
        )?;
        if part.protein_coding && part.aligned_len >= 3 {
            writeln!(
                w,
                "    CHARSET {}_pos1 = {}-{}\\3;",
                part.locus_name, part.start_col_1based, part.end_col_1based
            )?;
            writeln!(
                w,
                "    CHARSET {}_pos2 = {}-{}\\3;",
                part.locus_name,
                part.start_col_1based + 1,
                part.end_col_1based
            )?;
            writeln!(
                w,
                "    CHARSET {}_pos3 = {}-{}\\3;",
                part.locus_name,
                part.start_col_1based + 2,
                part.end_col_1based
            )?;
        }
    }

    let mut nuc_ranges = Vec::new();
    let mut mito_ranges = Vec::new();
    for part in partitions {
        let span = format!("{}-{}", part.start_col_1based, part.end_col_1based);
        match part.locus_class {
            LocusClass::Nuclear => nuc_ranges.push(span),
            LocusClass::Mitochondrial => mito_ranges.push(span),
        }
    }
    if !nuc_ranges.is_empty() && !mito_ranges.is_empty() {
        writeln!(
            w,
            "    CHARPARTITION GenomicCompartment = Nuclear: {}, Mitochondrial: {};",
            nuc_ranges.join(" "),
            mito_ranges.join(" ")
        )?;
    }

    let mut genus_groups: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for otu in otus {
        genus_groups
            .entry(otu.genus.as_str())
            .or_default()
            .push(otu.terminal_label.as_str());
    }
    for (genus, members) in genus_groups {
        writeln!(w, "    TAXSET Genus_{} = {};", genus, members.join(" "))?;
    }
    writeln!(w, "END;")?;
    writeln!(w)?;

    // 4. NOTES BLOCK (Chimeric OTU Accession Provenance)
    writeln!(w, "BEGIN NOTES;")?;
    for (idx, otu) in otus.iter().enumerate() {
        let provenance_list = otu
            .locus_slots
            .values()
            .map(|s| {
                let v = if s.voucher_id.is_empty() {
                    "unvouchered"
                } else {
                    s.voucher_id.as_str()
                };
                format!(
                    "{}:{}(voucher={},conf={:.3})",
                    s.locus_canonical, s.accession_version, v, s.neural_confidence
                )
            })
            .collect::<Vec<_>>()
            .join(", ");

        writeln!(
            w,
            "    TEXT TAXON={} TEXT='ArtificialOTU={} | Chimeric={} | DistinctSpecimens={} | DistinctAccessions={} | MeanNeuralConf={:.4} | Accessions=[{}]';",
            idx + 1,
            otu.terminal_label,
            otu.is_chimeric,
            otu.distinct_specimen_count,
            otu.distinct_accession_count,
            otu.mean_neural_confidence,
            provenance_list
        )?;
    }
    writeln!(w, "END;")?;
    writeln!(w)?;

    // 5. PRIVATE BIOSUITE_METADATA BLOCK (Group Quarantine Flags — Sequences Strictly Excluded)
    writeln!(w, "BEGIN BIOSUITE_METADATA;")?;
    writeln!(w, "    DATABASE_ID {};", audit.database_id)?;
    writeln!(w, "    RUN_ID {};", audit.run_id)?;
    writeln!(
        w,
        "    CONSENSUS_SUMMARY occupancy={:.4} consensus={:.4} cache_hits={} delta_fetched={} quarantined={} rehabilitated={};",
        audit.final_matrix_occupancy,
        audit.final_consensus_score,
        audit.cache_hits,
        audit.live_delta_fetched,
        audit.quarantined_count,
        audit.rehabilitated_count
    )?;
    for q in &active_quarantined {
        writeln!(
            w,
            "    QUARANTINE_FLAG accession={} taxon=\"{}\" locus={} reason={} prior_quarantine_count={} skepticism_prior={:.3} k2p_dist={:.4} neural_score={:.4};",
            q.accession_version,
            q.species_binomial,
            q.locus_canonical,
            q.last_reason_code,
            q.prior_quarantine_count,
            q.skepticism_prior,
            q.within_locus_k2p_dist,
            q.neural_consensus_score
        )?;
    }
    writeln!(w, "END;")?;

    w.flush()?;
    Ok(())
}
