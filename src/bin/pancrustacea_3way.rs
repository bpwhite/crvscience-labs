use bio_suite::fetch::FastaRecord;
use bio_suite::modeltest::run_native_modeltest;
use bio_suite::phylo::{align_sequences_with_mode, AlignedSequence, AlignerMode};
use bio_suite::tree_builder::{infer_native_ml_gene_tree, PhyloRunModel, SingleTreeRunMetrics};
use rayon::prelude::*;
use serde::Serialize;
use std::env;
use std::fs;
use std::time::Instant;

#[derive(Serialize, Clone)]
struct CandidateChainTrace {
    chain_id: String,
    model_code: String,
    engine_name: String,
    duration_ms: f64,
    ln_l: f64,
    aic: f64,
    bic: f64,
    total_tree_length: f64,
    mean_support: f64,
    nni_evaluated: usize,
    nni_accepted: usize,
}

#[derive(Serialize)]
struct NativeRust3WayOutput {
    msa_method: String,
    msa_duration_ms: f64,
    alignment_length_bp: usize,
    modeltest_best_bic: String,
    modeltest_best_aic: String,
    modeltest_duration_ms: f64,
    modeltest_best_ln_l: f64,
    modeltest_best_bic_score: f64,
    modeltest_candidates_evaluated: usize,
    candidate_chains_evaluated: usize,
    candidate_chains: Vec<CandidateChainTrace>,
    ml_tree_engine: String,
    ml_tree_model: String,
    ml_tree_duration_ms: f64,
    ml_tree_ln_l: f64,
    ml_tree_aic: f64,
    ml_tree_bic: f64,
    ml_tree_total_length: f64,
    ml_tree_mean_support: f64,
    ml_tree_nni_evaluated: usize,
    ml_tree_nni_accepted: usize,
    ml_tree_newick: String,
    total_native_duration_ms: f64,
}

fn load_fasta_records(path: &str) -> Vec<FastaRecord> {
    let content = fs::read_to_string(path).expect("Failed reading input FASTA");
    let mut records = Vec::new();
    let mut cur_id = String::new();
    let mut cur_seq = String::new();
    for line in content.lines() {
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        if let Some(rest) = t.strip_prefix('>') {
            if !cur_id.is_empty() {
                records.push(FastaRecord {
                    id: cur_id.clone(),
                    accession: cur_id.clone(),
                    taxon: cur_id.clone(),
                    common_name: cur_id.clone(),
                    family: "Panarthropoda".to_string(),
                    gene: "COX1".to_string(),
                    description: cur_id.clone(),
                    sequence: cur_seq.clone(),
                    cached_from_disk: true,
                });
            }
            cur_id = rest.split_whitespace().next().unwrap_or("UNKNOWN").to_string();
            cur_seq.clear();
        } else {
            cur_seq.push_str(t);
        }
    }
    if !cur_id.is_empty() {
        records.push(FastaRecord {
            id: cur_id.clone(),
            accession: cur_id.clone(),
            taxon: cur_id.clone(),
            common_name: cur_id.clone(),
            family: "Panarthropoda".to_string(),
            gene: "COX1".to_string(),
            description: cur_id,
            sequence: cur_seq,
            cached_from_disk: true,
        });
    }
    records
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let in_fasta = args
        .get(1)
        .expect("Usage: pancrustacea_3way <input.fasta> <output.json> [--ry3]");
    let out_json = args
        .get(2)
        .expect("Usage: pancrustacea_3way <input.fasta> <output.json> [--ry3]");
    let apply_ry3 = args.iter().any(|a| a == "--ry3");

    let t_total = Instant::now();
    let records = load_fasta_records(in_fasta);
    let total_raw_bp: usize = records.iter().map(|r| r.sequence.len()).sum();
    eprintln!(
        "[pancrustacea_3way] Loaded {} FASTA records ({} total bp) from {} (apply_ry3={})",
        records.len(),
        total_raw_bp,
        in_fasta,
        apply_ry3
    );

    // Stage 1: Native Rust Gotoh + Codon-Anchored MSA on the full downloaded sequences
    let (msa_method, mut aligned_seqs, msa_tel) =
        align_sequences_with_mode(&records, AlignerMode::Native);

    // Stage 1b: Optional Native Rust RY3 post-MSA compositional recoding (every 3rd aligned column: G->A, C->T)
    if apply_ry3 {
        aligned_seqs = aligned_seqs
            .into_iter()
            .map(|mut s| {
                let mut bytes = s.aligned_seq.into_bytes();
                for (i, b) in bytes.iter_mut().enumerate() {
                    if i % 3 == 2 {
                        if *b == b'G' {
                            *b = b'A';
                        } else if *b == b'C' {
                            *b = b'T';
                        }
                    }
                }
                s.aligned_seq = String::from_utf8(bytes).unwrap_or_default();
                s
            })
            .collect::<Vec<AlignedSequence>>();
    }

    let raw_aligned_strings: Vec<String> =
        aligned_seqs.iter().map(|s| s.aligned_seq.clone()).collect();
    let labels: Vec<String> = aligned_seqs.iter().map(|s| s.id.clone()).collect();
    let align_len = raw_aligned_strings.first().map(|s| s.len()).unwrap_or(0);
    eprintln!(
        "[pancrustacea_3way] Native MSA ({}) completed in {:.2} ms ({} aligned columns across {} taxa)",
        msa_method, msa_tel.duration_ms, align_len, labels.len()
    );

    // Stage 2: Native 24-Model Felsenstein Pruning + Yang G4 ModelTest
    let mt = run_native_modeltest(&aligned_seqs);
    eprintln!(
        "[pancrustacea_3way] Native ModelTest evaluated {} models in {:.2} ms (Best BIC: {}, ln L = {:.2})",
        mt.models_evaluated.len(),
        mt.duration_ms,
        mt.best_model_bic,
        mt.best_bic_ln_l
    );

    // Stage 3: RuNIE Multi-Chain Native Rust Pipeline Optimization across all 8 PhyloRunModel candidates in parallel
    let candidate_models = [
        PhyloRunModel::JcUniform,
        PhyloRunModel::K80Gamma4,
        PhyloRunModel::HkyInvGamma4,
        PhyloRunModel::Tn93InvGamma4,
        PhyloRunModel::Tim2InvGamma4,
        PhyloRunModel::GtrInvGamma4Fast,
        PhyloRunModel::GtrInvGamma4Full,
        PhyloRunModel::GtrFreeRate2,
    ];

    let evaluated: Vec<(usize, SingleTreeRunMetrics)> = candidate_models
        .par_iter()
        .enumerate()
        .map(|(idx, &m)| {
            let res = infer_native_ml_gene_tree(&raw_aligned_strings, &labels, m);
            eprintln!(
                "  [Native Chain {}/{}] Model {:<18} -> ln L = {:.2}, BIC = {:.2}, len = {:.4}, NNI = {}/{}, time = {:.2} ms",
                idx + 1,
                candidate_models.len(),
                m.code(),
                res.ln_l,
                res.bic,
                res.total_tree_length,
                res.nni_moves_accepted,
                res.nni_moves_evaluated,
                res.duration_ms
            );
            (idx, res)
        })
        .collect();

    let mut candidate_chains = Vec::new();
    let mut best_tree: Option<SingleTreeRunMetrics> = None;

    for (idx, res) in evaluated {
        candidate_chains.push(CandidateChainTrace {
            chain_id: format!("NATIVE-CHAIN-{:02}", idx + 1),
            model_code: res.run_model_code.clone(),
            engine_name: res.engine_name.clone(),
            duration_ms: res.duration_ms,
            ln_l: res.ln_l,
            aic: res.aic,
            bic: res.bic,
            total_tree_length: res.total_tree_length,
            mean_support: res.mean_branch_support_pct,
            nni_evaluated: res.nni_moves_evaluated,
            nni_accepted: res.nni_moves_accepted,
        });

        if best_tree
            .as_ref()
            .map(|b| res.ln_l > b.ln_l)
            .unwrap_or(true)
        {
            best_tree = Some(res);
        }
    }

    // Sort candidate chains by ln_l descending
    candidate_chains.sort_by(|a, b| b.ln_l.partial_cmp(&a.ln_l).unwrap_or(std::cmp::Ordering::Equal));
    let winner = best_tree.expect("At least one candidate chain must succeed");
    let total_ms = t_total.elapsed().as_secs_f64() * 1000.0;

    eprintln!(
        "[pancrustacea_3way] RuNIE Native Winner: {} (ln L = {:.2}, BIC = {:.2}, NNI = {}/{}, winner tree time = {:.2} ms, total 8-chain sweep = {:.2} ms)",
        winner.run_model_code,
        winner.ln_l,
        winner.bic,
        winner.nni_moves_accepted,
        winner.nni_moves_evaluated,
        winner.duration_ms,
        total_ms
    );

    let out = NativeRust3WayOutput {
        msa_method,
        msa_duration_ms: msa_tel.duration_ms,
        alignment_length_bp: align_len,
        modeltest_best_bic: mt.best_model_bic,
        modeltest_best_aic: mt.best_model_aic,
        modeltest_duration_ms: mt.duration_ms,
        modeltest_best_ln_l: mt.best_bic_ln_l,
        modeltest_best_bic_score: mt.best_bic_score,
        modeltest_candidates_evaluated: mt.models_evaluated.len(),
        candidate_chains_evaluated: candidate_chains.len(),
        candidate_chains,
        ml_tree_engine: winner.engine_name,
        ml_tree_model: winner.run_model_code,
        ml_tree_duration_ms: winner.duration_ms,
        ml_tree_ln_l: winner.ln_l,
        ml_tree_aic: winner.aic,
        ml_tree_bic: winner.bic,
        ml_tree_total_length: winner.total_tree_length,
        ml_tree_mean_support: winner.mean_branch_support_pct,
        ml_tree_nni_evaluated: winner.nni_moves_evaluated,
        ml_tree_nni_accepted: winner.nni_moves_accepted,
        ml_tree_newick: winner.newick_tree,
        total_native_duration_ms: total_ms,
    };

    fs::write(
        out_json,
        serde_json::to_string_pretty(&out).expect("Failed serializing output JSON"),
    )
    .expect("Failed writing output JSON");
}
