//! `seqnet_runie` — RuNIE pipeline fed by `bio-suite-seqnet` (`bio-suite seqnet-runie`).
//!
//! Stage 1 is `bio-suite-seqnet`: query-driven NCBI Entrez acquisition into a SQLite cache,
//! sequence QC (NUMT stop codons, ambiguity, within-locus outliers, neural consensus) with a
//! non-destructive quarantine ledger, and species-OTU assembly with gap-refill feedback. It
//! supersedes `fetch::PoliteCacheManager` and hard-coded accession lists for new pipelines.
//!
//! 1. seqnet acquisition, QC and OTU assembly
//! 2. per-locus native MAFFT-mimetic alignment
//! 3. per-locus native model selection (BIC) + native ML gene trees + supermatrix ML tree
//! 4. native quartet MSC species tree from the gene trees
//! 5. validation of both trees against clades reported in published studies (benchmark JSON)

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use bio_suite_seqnet::config::{AnalysisLevel, ConsensusConfig, ExecutionMode, OtuRank};
use bio_suite_seqnet::neurons::{run_neural_network_schedule, EngineControlFlags};
use bio_suite_seqnet::otu::AssembledSpeciesOtu;

use crate::coalescent::run_native_quartet_coalescent;
use crate::fetch::FastaRecord;
use crate::mafft_native::align_sequences_mafft_native;
use crate::modeltest::run_native_modeltest;
use crate::phylo::AlignedSequence;
use crate::tree_builder::{infer_native_ml_gene_tree, PhyloRunModel};

// ---------------------------------------------------------------------------
// Benchmark definition (published reference clades)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct BenchmarkSource {
    pub key: String,
    pub citation: String,
    pub doi: String,
    #[serde(default)]
    pub pmcid: String,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ReferenceClade {
    pub id: String,
    pub name: String,
    pub families: Vec<String>,
    #[serde(default)]
    pub why_not_scored: Option<String>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct CladeBenchmark {
    pub benchmark_id: String,
    pub title: String,
    pub sources: Vec<BenchmarkSource>,
    pub loci: Vec<String>,
    pub sampling_families: BTreeMap<String, Vec<String>>,
    pub scored_clades: Vec<ReferenceClade>,
    #[serde(default)]
    pub report_only_clades: Vec<ReferenceClade>,
}

impl CladeBenchmark {
    pub fn load(path: &Path) -> Result<Self, String> {
        let raw = fs::read_to_string(path).map_err(|e| format!("read {}: {}", path.display(), e))?;
        serde_json::from_str(raw.trim_start_matches('\u{feff}'))
            .map_err(|e| format!("parse {}: {}", path.display(), e))
    }

    pub fn all_families(&self) -> Vec<String> {
        self.sampling_families.values().flatten().cloned().collect()
    }
}

// ---------------------------------------------------------------------------
// Pipeline specification & report
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
pub struct SeqnetRunieSpec {
    pub benchmark_path: PathBuf,
    pub run_dir: PathBuf,
    pub sqlite_path: PathBuf,
    pub database_id: String,
    pub retmax_per_stratum: usize,
    pub refill_directives: usize,
    pub min_coverage: f64,
    pub offline_cache_only: bool,
    pub max_epochs: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct LocusStageReport {
    pub locus: String,
    pub taxa: usize,
    pub aligned_length: usize,
    pub best_model_bic: String,
    pub ml_model_used: String,
    pub gene_tree_ln_l: f64,
    pub align_ms: f64,
    pub tree_ms: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct CladeValidation {
    pub id: String,
    pub name: String,
    pub scored: bool,
    pub testable: bool,
    pub tips_inside: usize,
    pub tips_outside: usize,
    pub recovered: bool,
    /// Support label on the recovering branch when the tree carries one.
    pub support: Option<f64>,
    /// Best-matching split when not recovered: (Jaccard, intruders, missing).
    pub best_jaccard: f64,
    pub intruders: Vec<String>,
    pub missing: Vec<String>,
    pub note: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct TreeValidation {
    pub tree: String,
    pub tips: usize,
    pub scored_testable: usize,
    pub scored_recovered: usize,
    pub recovery_fraction: f64,
    pub clades: Vec<CladeValidation>,
}

#[derive(Debug, Clone, Serialize)]
pub struct SeqnetRunieReport {
    pub benchmark_id: String,
    pub run_dir: String,
    pub seqnet_run_id: String,
    pub seqnet_converged: bool,
    pub seqnet_epochs: usize,
    pub sequences_cached: usize,
    pub sequences_fetched_live: usize,
    pub sequences_quarantined: usize,
    pub otus: usize,
    pub families_sampled: BTreeMap<String, usize>,
    pub families_missing: Vec<String>,
    pub matrix_occupancy: f64,
    pub loci: Vec<LocusStageReport>,
    pub supermatrix_length: usize,
    pub supermatrix_ln_l: f64,
    pub msc_normalized_quartet_score: f64,
    pub validation_supermatrix: TreeValidation,
    pub validation_msc: Option<TreeValidation>,
    pub stage_seconds: BTreeMap<String, f64>,
}

// ---------------------------------------------------------------------------
// Pipeline
// ---------------------------------------------------------------------------

pub fn run_seqnet_runie_pipeline(spec: &SeqnetRunieSpec) -> Result<SeqnetRunieReport, String> {
    let bench = CladeBenchmark::load(&spec.benchmark_path)?;
    let families = bench.all_families();
    let aln_dir = spec.run_dir.join("alignments");
    let tree_dir = spec.run_dir.join("trees");
    fs::create_dir_all(&aln_dir).map_err(|e| e.to_string())?;
    fs::create_dir_all(&tree_dir).map_err(|e| e.to_string())?;
    let mut stage_seconds = BTreeMap::new();

    // ---- Stage 1: seqnet ------------------------------------------------------------------
    eprintln!(
        "[seqnet-runie] Stage 1/5: seqnet acquisition + QC + OTU assembly ({} families x {} loci)",
        families.len(),
        bench.loci.len()
    );
    let t = Instant::now();
    let cfg = ConsensusConfig {
        taxon_group: "Coleoptera".to_string(),
        analysis_level: AnalysisLevel::MultiLocusPhylogenetic,
        target_loci_count: bench.loci.len(),
        target_nuclear_loci: 0,
        target_mito_loci: 0,
        target_otu_count: 100_000,
        otu_rank: OtuRank::Species,
        min_species_per_genus: 1,
        min_locus_coverage_ratio: spec.min_coverage,
        max_epochs: spec.max_epochs,
        sqlite_path: spec.sqlite_path.clone(),
        database_id: spec.database_id.clone(),
        nexus_output_path: spec.run_dir.join("seqnet_supermatrix.nex"),
        allow_live_entrez: !spec.offline_cache_only,
        execution_mode: ExecutionMode::Headless,
        locus_names: bench.loci.clone(),
        sampling_taxa: families.clone(),
        retmax_per_stratum: spec.retmax_per_stratum,
        max_refill_directives: spec.refill_directives,
        max_refill_loci_per_directive: bench.loci.len(),
        refill_min_coverage: 0.5,
        ..ConsensusConfig::default()
    };
    let seq = run_neural_network_schedule(&cfg, None, None, EngineControlFlags::new())
        .map_err(|e| format!("seqnet stage failed: {e:#}"))?;
    stage_seconds.insert("1_seqnet".to_string(), t.elapsed().as_secs_f64());

    // Keep only OTUs from the benchmark's families (taxonomy strings name the family directly)
    let fam_set: BTreeSet<&str> = families.iter().map(|s| s.as_str()).collect();
    let otus: Vec<AssembledSpeciesOtu> = seq
        .otus
        .iter()
        .filter(|o| fam_set.contains(o.family.as_str()))
        .cloned()
        .collect();
    if otus.len() < 8 {
        return Err(format!("only {} OTUs after seqnet; too few for a tree", otus.len()));
    }
    let tip_family: HashMap<String, String> =
        otus.iter().map(|o| (o.terminal_label.clone(), o.family.clone())).collect();
    let mut families_sampled: BTreeMap<String, usize> = BTreeMap::new();
    for o in &otus {
        *families_sampled.entry(o.family.clone()).or_default() += 1;
    }
    let families_missing: Vec<String> = families
        .iter()
        .filter(|f| !families_sampled.contains_key(*f))
        .cloned()
        .collect();
    eprintln!(
        "[seqnet-runie]   {} OTUs from {} of {} families (occupancy {:.1}%, {} quarantined)",
        otus.len(),
        families_sampled.len(),
        families.len(),
        seq.audit.final_matrix_occupancy * 100.0,
        seq.audit.quarantined_count
    );

    // ---- Stage 2-3: per-locus alignment, model selection, gene trees -----------------------
    eprintln!("[seqnet-runie] Stage 2/5 + 3/5: native MSA, model selection and ML gene trees per locus");
    let labels_all: Vec<String> = otus.iter().map(|o| o.terminal_label.clone()).collect();
    let mut locus_reports = Vec::new();
    let mut gene_trees: Vec<String> = Vec::new();
    let mut aligned_by_locus: Vec<(String, usize, HashMap<String, String>)> = Vec::new();
    let t = Instant::now();
    for locus in &bench.loci {
        let records: Vec<FastaRecord> = otus
            .iter()
            .filter_map(|o| {
                o.locus_slots.get(locus).map(|slot| FastaRecord {
                    id: o.terminal_label.clone(),
                    accession: slot.accession_version.clone(),
                    taxon: o.species_binomial.clone(),
                    common_name: String::new(),
                    family: o.family.clone(),
                    gene: locus.clone(),
                    description: format!("{} {} {}", slot.accession_version, o.species_binomial, locus),
                    sequence: slot.sequence_raw.clone(),
                    cached_from_disk: true,
                })
            })
            .collect();
        if records.len() < 4 {
            eprintln!("[seqnet-runie]   {}: {} taxa, skipped", locus, records.len());
            continue;
        }
        let ta = Instant::now();
        let (_, aligned, _) = align_sequences_mafft_native(&records);
        let align_ms = ta.elapsed().as_secs_f64() * 1000.0;
        write_fasta(&aln_dir.join(format!("{locus}.fasta")), &aligned)?;
        let aln_len = aligned.first().map(|a| a.aligned_seq.len()).unwrap_or(0);

        let mt = run_native_modeltest(&aligned);
        let model = map_model(&mt.best_model_bic);
        let tt = Instant::now();
        let labels: Vec<String> = aligned.iter().map(|a| a.id.clone()).collect();
        let seqs: Vec<String> = aligned.iter().map(|a| a.aligned_seq.clone()).collect();
        let gt = infer_native_ml_gene_tree(&seqs, &labels, model);
        let tree_ms = tt.elapsed().as_secs_f64() * 1000.0;
        fs::write(tree_dir.join(format!("{locus}.nwk")), &gt.newick_tree).map_err(|e| e.to_string())?;
        gene_trees.push(gt.newick_tree.clone());
        eprintln!(
            "[seqnet-runie]   {:<9} {:>4} taxa x {:>5} bp | BIC {:<10} -> {:<8} | lnL {:.1}",
            locus, aligned.len(), aln_len, mt.best_model_bic, model.code(), gt.ln_l
        );
        locus_reports.push(LocusStageReport {
            locus: locus.clone(),
            taxa: aligned.len(),
            aligned_length: aln_len,
            best_model_bic: mt.best_model_bic.clone(),
            ml_model_used: model.code().to_string(),
            gene_tree_ln_l: gt.ln_l,
            align_ms,
            tree_ms,
        });
        aligned_by_locus.push((
            locus.clone(),
            aln_len,
            aligned.into_iter().map(|a| (a.id, a.aligned_seq)).collect(),
        ));
    }
    stage_seconds.insert("2_3_loci".to_string(), t.elapsed().as_secs_f64());

    // Supermatrix (absent loci padded as gaps = missing data in the ML likelihood)
    let t = Instant::now();
    let mut super_rows: Vec<String> = Vec::with_capacity(labels_all.len());
    for label in &labels_all {
        let mut row = String::new();
        for (_, len, rows) in &aligned_by_locus {
            match rows.get(label) {
                Some(s) => row.push_str(s),
                None => row.push_str(&"-".repeat(*len)),
            }
        }
        super_rows.push(row);
    }
    let super_len = super_rows.first().map(|r| r.len()).unwrap_or(0);
    let super_aligned: Vec<AlignedSequence> = labels_all
        .iter()
        .zip(super_rows.iter())
        .map(|(l, s)| AlignedSequence {
            id: l.clone(),
            taxon: l.clone(),
            common_name: String::new(),
            family: tip_family.get(l).cloned().unwrap_or_default(),
            aligned_seq: s.clone(),
        })
        .collect();
    write_fasta(&aln_dir.join("supermatrix.fasta"), &super_aligned)?;
    eprintln!(
        "[seqnet-runie]   supermatrix {} taxa x {} bp, ML tree (GTR+F+I+G4 fast NNI)...",
        labels_all.len(),
        super_len
    );
    let st = infer_native_ml_gene_tree(&super_rows, &labels_all, PhyloRunModel::GtrInvGamma4Fast);
    fs::write(tree_dir.join("supermatrix.nwk"), &st.newick_tree).map_err(|e| e.to_string())?;
    stage_seconds.insert("3_supermatrix_ml".to_string(), t.elapsed().as_secs_f64());

    // ---- Stage 4: quartet MSC ---------------------------------------------------------------
    eprintln!("[seqnet-runie] Stage 4/5: native quartet MSC species tree from {} gene trees", gene_trees.len());
    let t = Instant::now();
    let msc = if gene_trees.len() >= 2 {
        match run_native_quartet_coalescent(&bench.benchmark_id, &gene_trees) {
            Ok((res, _, _)) => {
                fs::write(tree_dir.join("msc_species_tree.nwk"), &res.species_tree_newick)
                    .map_err(|e| e.to_string())?;
                Some(res)
            }
            Err(e) => {
                eprintln!("[seqnet-runie]   quartet MSC failed: {e}");
                None
            }
        }
    } else {
        None
    };
    stage_seconds.insert("4_msc".to_string(), t.elapsed().as_secs_f64());

    // ---- Stage 5: validation against published clades --------------------------------------
    eprintln!("[seqnet-runie] Stage 5/5: validation against {} published clades", bench.scored_clades.len());
    let validation_supermatrix = validate_tree("supermatrix_ml", &st.newick_tree, &tip_family, &bench);
    let validation_msc = msc
        .as_ref()
        .map(|m| validate_tree("quartet_msc", &m.species_tree_newick, &tip_family, &bench));

    let report = SeqnetRunieReport {
        benchmark_id: bench.benchmark_id.clone(),
        run_dir: spec.run_dir.display().to_string(),
        seqnet_run_id: seq.audit.run_id.clone(),
        seqnet_converged: seq.audit.converged,
        seqnet_epochs: seq.audit.epochs_executed,
        sequences_cached: seq.audit.cache_hits,
        sequences_fetched_live: seq.audit.live_delta_fetched,
        sequences_quarantined: seq.audit.quarantined_count,
        otus: otus.len(),
        families_sampled,
        families_missing,
        matrix_occupancy: seq.audit.final_matrix_occupancy,
        loci: locus_reports,
        supermatrix_length: super_len,
        supermatrix_ln_l: st.ln_l,
        msc_normalized_quartet_score: msc.as_ref().map(|m| m.normalized_quartet_score).unwrap_or(0.0),
        validation_supermatrix,
        validation_msc,
        stage_seconds,
    };
    fs::write(
        spec.run_dir.join("seqnet_runie_report.json"),
        serde_json::to_string_pretty(&report).unwrap_or_default(),
    )
    .map_err(|e| e.to_string())?;
    Ok(report)
}

fn map_model(best_bic: &str) -> PhyloRunModel {
    let base = best_bic.split('+').next().unwrap_or("").trim().to_ascii_uppercase();
    match base.as_str() {
        "JC" | "F81" => PhyloRunModel::JcUniform,
        "K80" | "K2P" => PhyloRunModel::K80Gamma4,
        "HKY" => PhyloRunModel::HkyInvGamma4,
        "TN93" | "TRN" | "TN" => PhyloRunModel::Tn93InvGamma4,
        _ => PhyloRunModel::GtrInvGamma4Fast,
    }
}

fn write_fasta(path: &Path, rows: &[AlignedSequence]) -> Result<(), String> {
    let mut out = String::new();
    for r in rows {
        out.push('>');
        out.push_str(&r.id);
        out.push('\n');
        for chunk in r.aligned_seq.as_bytes().chunks(80) {
            out.push_str(&String::from_utf8_lossy(chunk));
            out.push('\n');
        }
    }
    fs::write(path, out).map_err(|e| format!("write {}: {}", path.display(), e))
}

// ---------------------------------------------------------------------------
// Newick splits & clade validation
// ---------------------------------------------------------------------------

/// Parsed Newick: tip labels and, for every internal node, its descendant tip set + support label.
#[derive(Debug, Default)]
pub struct NewickSplits {
    pub tips: Vec<String>,
    pub clades: Vec<(BTreeSet<String>, Option<f64>)>,
}

pub fn parse_newick_splits(newick: &str) -> Result<NewickSplits, String> {
    let s: Vec<char> = newick.trim().chars().collect();
    let mut pos = 0usize;
    let mut out = NewickSplits::default();
    parse_node(&s, &mut pos, &mut out)?;
    Ok(out)
}

fn read_label(s: &[char], pos: &mut usize) -> String {
    let mut label = String::new();
    if *pos < s.len() && s[*pos] == '\'' {
        *pos += 1;
        while *pos < s.len() && s[*pos] != '\'' {
            label.push(s[*pos]);
            *pos += 1;
        }
        *pos += 1;
        return label;
    }
    while *pos < s.len() && !matches!(s[*pos], ':' | ',' | ')' | '(' | ';' | '[') {
        label.push(s[*pos]);
        *pos += 1;
    }
    label.trim().to_string()
}

fn skip_branch_length_and_comments(s: &[char], pos: &mut usize) {
    loop {
        if *pos < s.len() && s[*pos] == '[' {
            while *pos < s.len() && s[*pos] != ']' {
                *pos += 1;
            }
            *pos += 1;
        } else if *pos < s.len() && s[*pos] == ':' {
            *pos += 1;
            while *pos < s.len() && !matches!(s[*pos], ',' | ')' | ';' | '[') {
                *pos += 1;
            }
        } else {
            break;
        }
    }
}

fn parse_node(s: &[char], pos: &mut usize, out: &mut NewickSplits) -> Result<BTreeSet<String>, String> {
    if *pos < s.len() && s[*pos] == '(' {
        *pos += 1;
        let mut set = BTreeSet::new();
        loop {
            let child = parse_node(s, pos, out)?;
            set.extend(child);
            match s.get(*pos) {
                Some(',') => *pos += 1,
                Some(')') => {
                    *pos += 1;
                    break;
                }
                other => return Err(format!("unexpected {:?} at {}", other, *pos)),
            }
        }
        let label = read_label(s, pos);
        let support = label.parse::<f64>().ok();
        skip_branch_length_and_comments(s, pos);
        out.clades.push((set.clone(), support));
        Ok(set)
    } else {
        let label = read_label(s, pos);
        skip_branch_length_and_comments(s, pos);
        if label.is_empty() {
            return Err(format!("empty tip label at {}", *pos));
        }
        out.tips.push(label.clone());
        Ok(BTreeSet::from([label]))
    }
}

pub fn validate_tree(
    tree_name: &str,
    newick: &str,
    tip_family: &HashMap<String, String>,
    bench: &CladeBenchmark,
) -> TreeValidation {
    let splits = match parse_newick_splits(newick) {
        Ok(s) => s,
        Err(e) => {
            return TreeValidation {
                tree: format!("{tree_name} (unparseable: {e})"),
                tips: 0,
                scored_testable: 0,
                scored_recovered: 0,
                recovery_fraction: 0.0,
                clades: Vec::new(),
            }
        }
    };
    let all: BTreeSet<String> = splits.tips.iter().cloned().collect();
    let mut clades = Vec::new();
    let defs = bench
        .scored_clades
        .iter()
        .map(|c| (c, true))
        .chain(bench.report_only_clades.iter().map(|c| (c, false)));
    for (def, scored) in defs {
        let fams: BTreeSet<&str> = def.families.iter().map(|f| f.as_str()).collect();
        let inside: BTreeSet<String> = all
            .iter()
            .filter(|t| tip_family.get(*t).map(|f| fams.contains(f.as_str())).unwrap_or(false))
            .cloned()
            .collect();
        let outside = all.len() - inside.len();
        let testable = inside.len() >= 2 && outside >= 2;
        let complement: BTreeSet<String> = all.difference(&inside).cloned().collect();
        let mut recovered = false;
        let mut support = None;
        let mut best = (0.0f64, BTreeSet::new());
        if testable {
            for (set, sup) in &splits.clades {
                if *set == inside || *set == complement {
                    recovered = true;
                    support = *sup;
                    break;
                }
                for side in [set.clone(), all.difference(set).cloned().collect::<BTreeSet<_>>()] {
                    let inter = side.intersection(&inside).count() as f64;
                    let uni = side.union(&inside).count() as f64;
                    let j = if uni > 0.0 { inter / uni } else { 0.0 };
                    if j > best.0 {
                        best = (j, side);
                    }
                }
            }
        }
        let (intruders, missing) = if recovered || !testable {
            (Vec::new(), Vec::new())
        } else {
            (
                best.1.difference(&inside).cloned().collect(),
                inside.difference(&best.1).cloned().collect(),
            )
        };
        clades.push(CladeValidation {
            id: def.id.clone(),
            name: def.name.clone(),
            scored,
            testable,
            tips_inside: inside.len(),
            tips_outside: outside,
            recovered,
            support,
            best_jaccard: if recovered { 1.0 } else { best.0 },
            intruders,
            missing,
            note: def.why_not_scored.clone(),
        });
    }
    let scored_testable = clades.iter().filter(|c| c.scored && c.testable).count();
    let scored_recovered = clades.iter().filter(|c| c.scored && c.testable && c.recovered).count();
    TreeValidation {
        tree: tree_name.to_string(),
        tips: all.len(),
        scored_testable,
        scored_recovered,
        recovery_fraction: if scored_testable > 0 {
            scored_recovered as f64 / scored_testable as f64
        } else {
            0.0
        },
        clades,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bench() -> CladeBenchmark {
        CladeBenchmark {
            benchmark_id: "t".into(),
            title: "t".into(),
            sources: vec![],
            loci: vec![],
            sampling_families: BTreeMap::new(),
            scored_clades: vec![
                ReferenceClade { id: "ab".into(), name: "AB".into(), families: vec!["A".into(), "B".into()], why_not_scored: None },
                ReferenceClade { id: "ac".into(), name: "AC".into(), families: vec!["A".into(), "C".into()], why_not_scored: None },
            ],
            report_only_clades: vec![],
        }
    }

    #[test]
    fn splits_are_rooting_independent() {
        let fam: HashMap<String, String> = [("a1", "A"), ("a2", "A"), ("b1", "B"), ("c1", "C"), ("c2", "C"), ("d1", "D")]
            .iter()
            .map(|(t, f)| (t.to_string(), f.to_string()))
            .collect();
        // ((a1,a2),b1) clade is AB; rerooting inside AB must still recover it via the complement split
        for nwk in [
            "(((a1:0.1,a2:0.1)90:0.2,b1:0.3)85:0.1,(c1,c2),d1);",
            "(a1,(a2,(b1,((c1,c2),d1)77)));",
        ] {
            let v = validate_tree("t", nwk, &fam, &bench());
            let ab = v.clades.iter().find(|c| c.id == "ab").unwrap();
            assert!(ab.testable && ab.recovered, "AB should be recovered in {nwk}");
            let ac = v.clades.iter().find(|c| c.id == "ac").unwrap();
            assert!(ac.testable && !ac.recovered, "AC should not be recovered in {nwk}");
            assert_eq!(v.scored_testable, 2);
            assert_eq!(v.scored_recovered, 1);
        }
    }

    #[test]
    fn model_mapping() {
        assert_eq!(map_model("GTR+I+G").code(), PhyloRunModel::GtrInvGamma4Fast.code());
        assert_eq!(map_model("HKY+G").code(), PhyloRunModel::HkyInvGamma4.code());
        assert_eq!(map_model("JC").code(), PhyloRunModel::JcUniform.code());
    }
}
