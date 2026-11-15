use bio_suite::fetch::FastaRecord;
use bio_suite::haplotype::analyze_haplotype_clusters_prealigned;
use std::fs;
use std::path::PathBuf;

fn load_fasta(path: &str) -> Vec<FastaRecord> {
    let text = fs::read_to_string(path).expect("read fasta");
    let mut out = Vec::new();
    let mut id = String::new();
    let mut taxon = String::new();
    let mut family = String::new();
    let mut seq = String::new();
    for line in text.lines() {
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        if let Some(h) = t.strip_prefix('>') {
            if !id.is_empty() && !seq.is_empty() {
                out.push(FastaRecord {
                    id: id.clone(),
                    accession: id.clone(),
                    taxon: taxon.clone(),
                    common_name: String::new(),
                    family: family.clone(),
                    gene: "locus".to_string(),
                    sequence: seq.clone(),
                    description: id.clone(),
                    cached_from_disk: true,
                });
            }
            let parts: Vec<&str> = h.split('|').collect();
            id = parts[0].split_whitespace().next().unwrap_or("S").to_string();
            taxon = parts
                .get(1)
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .unwrap_or_else(|| id.clone());
            family = parts
                .get(2)
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .unwrap_or_else(|| taxon.clone());
            seq.clear();
        } else {
            seq.push_str(t);
        }
    }
    if !id.is_empty() && !seq.is_empty() {
        out.push(FastaRecord {
            id: id.clone(),
            accession: id,
            taxon,
            common_name: String::new(),
            family,
            gene: "locus".to_string(),
            sequence: seq,
            description: String::new(),
            cached_from_disk: true,
        });
    }
    out
}

fn mutate_seq(base: &[u8], positions: &[usize], seed: u64) -> String {
    let mut s = base.to_vec();
    let bases = [b'A', b'C', b'G', b'T'];
    let n = s.len();
    for (idx, &pos) in positions.iter().enumerate() {
        let cur = s[pos % n];
        let next = bases[(cur as usize + idx + seed as usize + 1) % 4];
        s[pos % n] = if next == cur { bases[(idx + 2) % 4] } else { next };
    }
    String::from_utf8(s).unwrap()
}

fn build_multi_cluster_benchmark(
    n_clusters: usize,
    haps_per_cluster: usize,
) -> Vec<FastaRecord> {
    let len = 480;
    let bases = [b'A', b'C', b'G', b'T'];
    let root: Vec<u8> = (0..len).map(|i| bases[(i * 7 + 3) % 4]).collect();
    let pop_names = [
        "Pop_A_Coastal",
        "Pop_B_Highland",
        "Pop_C_Basin",
        "Pop_D_Island",
        "Pop_E_Estuary",
        "Pop_F_Steppe",
        "Pop_G_Montane",
        "Pop_H_Savanna",
        "Pop_I_Archipelago",
    ];

    let mut records = Vec::new();
    let mut seq_counter = 0usize;

    for c in 0..n_clusters {
        // Cluster hub separated by 7..13 SNPs from root
        let mut hub_muts = Vec::new();
        for m in 0..(7 + (c % 5)) {
            hub_muts.push((c * 37 + m * 3) % len);
        }
        let hub_seq = mutate_seq(&root, &hub_muts, (c + 1) as u64);
        let hub_bytes = hub_seq.as_bytes();

        // Dominant hub has 14..32 samples, some shared across 2-3 populations for pie charts
        let hub_samples = 14 + ((c * 5) % 19);
        for s in 0..hub_samples {
            seq_counter += 1;
            let pop = if s % 5 == 0 {
                pop_names[(c + 1) % pop_names.len()]
            } else if s % 9 == 0 {
                pop_names[(c + 2) % pop_names.len()]
            } else {
                pop_names[c % pop_names.len()]
            };
            records.push(FastaRecord {
                id: format!("S{:04}", seq_counter),
                accession: format!("ACC{:04}", seq_counter),
                taxon: pop.to_string(),
                common_name: String::new(),
                family: pop.to_string(),
                gene: "COX1".to_string(),
                sequence: hub_seq.clone(),
                description: pop.to_string(),
                cached_from_disk: true,
            });
        }

        // Satellites: 1..3 SNPs from hub
        for h in 1..haps_per_cluster {
            let n_muts = 1 + (h % 3);
            let mut sat_muts = Vec::new();
            for m in 0..n_muts {
                sat_muts.push((200 + c * 29 + h * 5 + m) % len);
            }
            let sat_seq = mutate_seq(hub_bytes, &sat_muts, (c * 100 + h) as u64);
            let sat_samples = if h % 4 == 0 { 4 + (h % 6) } else { 1 + (h % 2) };
            for s in 0..sat_samples {
                seq_counter += 1;
                let pop = if h % 3 == 0 && s % 2 == 1 {
                    pop_names[(c + 1) % pop_names.len()]
                } else {
                    pop_names[c % pop_names.len()]
                };
                records.push(FastaRecord {
                    id: format!("S{:04}", seq_counter),
                    accession: format!("ACC{:04}", seq_counter),
                    taxon: pop.to_string(),
                    common_name: String::new(),
                    family: pop.to_string(),
                    gene: "COX1".to_string(),
                    sequence: sat_seq.clone(),
                    description: pop.to_string(),
                    cached_from_disk: true,
                });
            }
        }
    }
    records
}

fn main() {
    let out_dir = PathBuf::from(
        r"C:\Users\bpcwh\.gemini\antigravity\brain\d72b3afe-a160-42e3-9e7e-58491c65c79e\scratch\haplo-calibration",
    );
    fs::create_dir_all(&out_dir).expect("create out_dir");

    // 1. Real 64-taxon Nematoda COX1
    let nem_cox1 = load_fasta("runs/runie-native-nematoda-pests-27/alignments/supp_s3_cox1_64taxa.fasta");
    let r_nem_cox1 = analyze_haplotype_clusters_prealigned(
        "calib-01-nematoda-cox1",
        "Calibration 1: Nematoda COX1 (64 Taxa)",
        &nem_cox1,
        15,
    );
    fs::write(out_dir.join("calib-01-nematoda-cox1.svg"), &r_nem_cox1.visualizations.svg_vector_network).unwrap();
    println!(
        "Calib 1 (Nematoda COX1): N={}, h={}, clusters={}, edges={}",
        r_nem_cox1.diversity_metrics.total_sequences_n,
        r_nem_cox1.diversity_metrics.unique_haplotypes_h,
        r_nem_cox1.clusters.len(),
        r_nem_cox1.msn_edges.len()
    );

    // 2. Medium 4-Cluster Pie-Chart Benchmark (~40 haplotypes)
    let med_recs = build_multi_cluster_benchmark(4, 10);
    let r_med = analyze_haplotype_clusters_prealigned(
        "calib-02-4cluster-pie",
        "Calibration 2: 4-Cluster Multi-Pop Pie Network (~40 Haplotypes)",
        &med_recs,
        3,
    );
    fs::write(out_dir.join("calib-02-4cluster-pie.svg"), &r_med.visualizations.svg_vector_network).unwrap();
    println!(
        "Calib 2 (4-Cluster Pie): N={}, h={}, clusters={}, edges={}",
        r_med.diversity_metrics.total_sequences_n,
        r_med.diversity_metrics.unique_haplotypes_h,
        r_med.clusters.len(),
        r_med.msn_edges.len()
    );

    // 3. Dense 9-Cluster Pie-Chart Stress Benchmark (~81 haplotypes)
    let dense_recs = build_multi_cluster_benchmark(9, 9);
    let r_dense = analyze_haplotype_clusters_prealigned(
        "calib-03-9cluster-dense",
        "Calibration 3: 9-Cluster Dense Pie Network (~81 Haplotypes)",
        &dense_recs,
        3,
    );
    fs::write(out_dir.join("calib-03-9cluster-dense.svg"), &r_dense.visualizations.svg_vector_network).unwrap();
    println!(
        "Calib 3 (9-Cluster Dense): N={}, h={}, clusters={}, edges={}",
        r_dense.diversity_metrics.total_sequences_n,
        r_dense.diversity_metrics.unique_haplotypes_h,
        r_dense.clusters.len(),
        r_dense.msn_edges.len()
    );
}
