//! Tier 1: Polite Data Acquisition & Content-Addressable Local Cache (`bio-suite fetch`)
//!
//! **Superseded for new pipelines.** Stage-1 data acquisition for RuNIE pipelines is
//! `bio-suite-seqnet` (workspace crate `seqnet/`, wired in by `seqnet_runie.rs` /
//! `bio-suite seqnet-runie`): query-driven Entrez discovery, SQLite cache with `mdat` delta
//! sync, sequence QC with a quarantine ledger, and species-OTU assembly. `PoliteCacheManager`
//! and the hard-coded accession lists built on it stay only so existing studies and published
//! lab reports remain reproducible; do not use them in new pipelines. `FastaRecord` remains
//! the shared sequence record type.
//!
//! Enforces strict NCBI E-utilities / Datasets rate-limiting (<= 2.5 req/s without key)
//! and SHA-256 content-addressable local disk caching so repeated pipelines never
//! re-hit NCBI servers unnecessarily.

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::fs;
use std::path::{Path, PathBuf};
use std::thread;
use std::time::{Duration, SystemTime, UNIX_EPOCH};

const NCBI_EUTILS_BASE: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";
const MIN_REQUEST_INTERVAL_MS: u64 = 400; // <= 2.5 req/sec polite limit (NCBI limit is 3 req/sec)

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FastaRecord {
    pub id: String,
    pub accession: String,
    pub taxon: String,
    pub common_name: String,
    pub family: String,
    pub gene: String,
    pub description: String,
    pub sequence: String,
    pub cached_from_disk: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenomicAnnotation {
    pub feature_id: String,
    pub feature_type: String, // "CDS", "gene", "repeat_region", "regulatory"
    pub symbol: String,
    pub start: usize,
    pub end: usize,
    pub strand: char,
    pub description: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoiSpecimenTarget {
    pub accession: &'static str,
    pub taxon: &'static str,
    pub common_name: &'static str,
    pub family: &'static str,
}

/// Curated 16-specimen NCBI RefSeq / GenBank Cytochrome c Oxidase Subunit I (COX1 / COI)
/// barcode reference panel across 4 Carnivora families (Ursidae, Felidae, Canidae, Mustelidae),
/// including multiple intraspecific pairs (Ursus arctos, Ursus maritimus, Ursus americanus,
/// Panthera tigris, Canis lupus) to rigorously test the Automatic Barcode Gap Discovery (ABGD)
/// threshold and phylogenetic clade resolution.
pub const CARNIVORA_COI_PANEL: &[CoiSpecimenTarget] = &[
    // Family Ursidae (Bears — 7 specimens across 4 species)
    CoiSpecimenTarget {
        accession: "NC_003427.1:6484-8028",
        taxon: "Ursus arctos",
        common_name: "Brown bear (Eurasian RefSeq)",
        family: "Ursidae",
    },
    CoiSpecimenTarget {
        accession: "AP012578.1:6392-7936",
        taxon: "Ursus arctos",
        common_name: "Brown bear (Hokkaido isolate)",
        family: "Ursidae",
    },
    CoiSpecimenTarget {
        accession: "NC_003428.1:6479-8023",
        taxon: "Ursus maritimus",
        common_name: "Polar bear (RefSeq)",
        family: "Ursidae",
    },
    CoiSpecimenTarget {
        accession: "GU573490.1:6272-7816",
        taxon: "Ursus maritimus",
        common_name: "Polar bear (Isolate 495)",
        family: "Ursidae",
    },
    CoiSpecimenTarget {
        accession: "NC_003426.1:6291-7835",
        taxon: "Ursus americanus",
        common_name: "American black bear (RefSeq)",
        family: "Ursidae",
    },
    CoiSpecimenTarget {
        accession: "JX196366.1:5883-7427",
        taxon: "Ursus americanus",
        common_name: "American black bear (Isolate BLK014)",
        family: "Ursidae",
    },
    CoiSpecimenTarget {
        accession: "NC_009971.1:6243-7787",
        taxon: "Ursus thibetanus",
        common_name: "Asian black bear (RefSeq)",
        family: "Ursidae",
    },
    // Family Felidae (Cats — 4 specimens across 3 species)
    CoiSpecimenTarget {
        accession: "NC_010642.1:6280-7824",
        taxon: "Panthera tigris",
        common_name: "Tiger (RefSeq)",
        family: "Felidae",
    },
    CoiSpecimenTarget {
        accession: "EF551003.1:6280-7824",
        taxon: "Panthera tigris",
        common_name: "Amur tiger (Isolate)",
        family: "Felidae",
    },
    CoiSpecimenTarget {
        accession: "NC_010641.1:6305-7849",
        taxon: "Panthera pardus",
        common_name: "Leopard (RefSeq)",
        family: "Felidae",
    },
    CoiSpecimenTarget {
        accession: "NC_001700.1:6216-7760",
        taxon: "Felis catus",
        common_name: "Domestic cat (RefSeq)",
        family: "Felidae",
    },
    // Family Canidae (Canids — 3 specimens across 2 species)
    CoiSpecimenTarget {
        accession: "NC_008092.1:5350-6894",
        taxon: "Canis lupus",
        common_name: "Gray wolf (RefSeq)",
        family: "Canidae",
    },
    CoiSpecimenTarget {
        accession: "NC_002008.4:5349-6893",
        taxon: "Canis lupus familiaris",
        common_name: "Domestic dog (RefSeq)",
        family: "Canidae",
    },
    CoiSpecimenTarget {
        accession: "NC_008434.1:5350-6894",
        taxon: "Vulpes vulpes",
        common_name: "Red fox (RefSeq)",
        family: "Canidae",
    },
    // Family Mustelidae (Mustelids — 2 specimens across 2 species)
    CoiSpecimenTarget {
        accession: "NC_009677.1:5339-6883",
        taxon: "Meles anakuma",
        common_name: "Japanese badger (RefSeq)",
        family: "Mustelidae",
    },
    CoiSpecimenTarget {
        accession: "NC_020638.1:5343-6887",
        taxon: "Mustela putorius",
        common_name: "European polecat (RefSeq)",
        family: "Mustelidae",
    },
];

pub struct PoliteCacheManager {
    pub cache_dir: PathBuf,
}

impl PoliteCacheManager {
    pub fn new(cache_dir: impl AsRef<Path>) -> Self {
        let dir = cache_dir.as_ref().to_path_buf();
        let _ = fs::create_dir_all(&dir);
        Self { cache_dir: dir }
    }

    /// Resolve the canonical `crvscience-labs/cache` directory whether invoked from workspace root or `crvscience-labs/`.
    pub fn default_workspace_cache() -> Self {
        let ws_cache = Path::new("crvscience-labs").join("cache");
        if ws_cache.exists() || Path::new("crvscience-labs").join("Cargo.toml").exists() {
            Self::new(ws_cache)
        } else {
            Self::new("cache")
        }
    }

    fn cache_key_path(&self, key: &str, ext: &str) -> PathBuf {
        let mut hasher = Sha256::new();
        hasher.update(key.as_bytes());
        let hash_hex = format!("{:x}", hasher.finalize());
        let safe_prefix: String = key
            .chars()
            .map(|c| if c.is_ascii_alphanumeric() || c == '_' || c == '-' { c } else { '_' })
            .take(36)
            .collect();
        self.cache_dir
            .join(format!("{}_{}.{}", safe_prefix, &hash_hex[..12], ext))
    }

    /// Enforce polite NCBI rate limiting across processes using a timestamp lockfile.
    fn wait_for_polite_slot(&self) {
        let lock_path = self.cache_dir.join(".ncbi_rate_limit.stamp");
        let now_ms = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_millis() as u64;

        if let Ok(contents) = fs::read_to_string(&lock_path) {
            if let Ok(last_ms) = contents.trim().parse::<u64>() {
                let elapsed = now_ms.saturating_sub(last_ms);
                if elapsed < MIN_REQUEST_INTERVAL_MS {
                    thread::sleep(Duration::from_millis(MIN_REQUEST_INTERVAL_MS - elapsed));
                }
            }
        }

        let updated_ms = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_millis() as u64;
        let _ = fs::write(&lock_path, updated_ms.to_string());
    }

    /// Fetch a nucleotide slice from NCBI Entrez `efetch.fcgi` with cache-first lookup.
    pub fn fetch_ncbi_fasta_slice(
        &self,
        accession_spec: &str,
        taxon: &str,
        common_name: &str,
        family: &str,
        gene: &str,
    ) -> Result<FastaRecord, String> {
        let cache_file = self.cache_key_path(accession_spec, "fasta");
        if cache_file.exists() {
            if let Ok(raw) = fs::read_to_string(&cache_file) {
                if let Some((desc, seq)) = parse_single_fasta(&raw) {
                    if !seq.is_empty() {
                        return Ok(FastaRecord {
                            id: accession_spec.to_string(),
                            accession: accession_spec.to_string(),
                            taxon: taxon.to_string(),
                            common_name: common_name.to_string(),
                            family: family.to_string(),
                            gene: gene.to_string(),
                            description: desc,
                            sequence: seq,
                            cached_from_disk: true,
                        });
                    }
                }
            }
        }

        // Parse "ACCESSION:START-STOP" or plain "ACCESSION"
        let (acc, seq_start, seq_stop) = if let Some((a, range)) = accession_spec.split_once(':') {
            if let Some((s, e)) = range.split_once('-') {
                (
                    a.trim(),
                    s.trim().parse::<usize>().ok(),
                    e.trim().parse::<usize>().ok(),
                )
            } else {
                (a.trim(), None, None)
            }
        } else {
            (accession_spec.trim(), None, None)
        };

        self.wait_for_polite_slot();

        let mut url = format!(
            "{}/efetch.fcgi?db=nuccore&id={}&rettype=fasta&retmode=text&tool=crvscience_bio_suite&email=labs@crvscience.com",
            NCBI_EUTILS_BASE, acc
        );
        if let (Some(s), Some(e)) = (seq_start, seq_stop) {
            url.push_str(&format!("&seq_start={}&seq_stop={}", s, e));
        }

        // Retry with exponential backoff on transient 429/5xx
        let mut last_err = String::new();
        for attempt in 0..3 {
            if attempt > 0 {
                thread::sleep(Duration::from_millis(750 * (1 << attempt)));
            }
            match ureq::get(&url)
                .set("User-Agent", "crvscience-bio-suite/0.1.0 (+https://crvscience.com/labs; GPL-3.0)")
                .timeout(Duration::from_secs(20))
                .call()
            {
                Ok(resp) => {
                    if let Ok(body) = resp.into_string() {
                        if let Some((desc, seq)) = parse_single_fasta(&body) {
                            if seq.len() >= 100 {
                                let _ = fs::write(&cache_file, &body);
                                return Ok(FastaRecord {
                                    id: accession_spec.to_string(),
                                    accession: accession_spec.to_string(),
                                    taxon: taxon.to_string(),
                                    common_name: common_name.to_string(),
                                    family: family.to_string(),
                                    gene: gene.to_string(),
                                    description: desc,
                                    sequence: seq,
                                    cached_from_disk: false,
                                });
                            }
                        }
                        last_err = format!("Empty or invalid FASTA body for {}", accession_spec);
                    }
                }
                Err(e) => {
                    last_err = format!("HTTP error fetching {}: {}", accession_spec, e);
                }
            }
        }

        Err(last_err)
    }

    /// Fetch the entire Carnivora COI Barcode panel (cache-first, polite NCBI fallback).
    pub fn fetch_carnivora_coi_panel(&self) -> Result<Vec<FastaRecord>, String> {
        let mut records = Vec::with_capacity(CARNIVORA_COI_PANEL.len());
        for item in CARNIVORA_COI_PANEL {
            let rec = self.fetch_ncbi_fasta_slice(
                item.accession,
                item.taxon,
                item.common_name,
                item.family,
                "COX1",
            )?;
            records.push(rec);
        }
        Ok(records)
    }

    /// Fetch Human Chromosome 1 Section A (GRCh38 NC_000001.11:1000000-1060000, 60 kb Telomeric/Subtelomeric 1p36.33 locus)
    /// along with curated RefSeq gene/repeat annotations for that locus.
    pub fn fetch_human_chr1_section_a(&self) -> Result<(FastaRecord, Vec<GenomicAnnotation>), String> {
        let rec = self.fetch_ncbi_fasta_slice(
            "NC_000001.11:1000001-1060000",
            "Homo sapiens",
            "Human GRCh38 chr1p36.33 Section A (1.00–1.06 Mb)",
            "Hominidae",
            "CHR1_1P36_33_SEC_A",
        )?;

        // Curated GRCh38.p14 RefSeq & RepeatMasker/Dfam annotations within NC_000001.11:1,000,001..1,060,000 (relative 0..60,000 bp)
        let annotations = vec![
            GenomicAnnotation {
                feature_id: "GENE_HES4".to_string(),
                feature_type: "CDS".to_string(),
                symbol: "HES4".to_string(),
                start: 1100,
                end: 2350,
                strand: '-',
                description: "hes family bHLH transcription factor 4 (CpG-dense Notch target)".to_string(),
            },
            GenomicAnnotation {
                feature_id: "GENE_ISG15".to_string(),
                feature_type: "CDS".to_string(),
                symbol: "ISG15".to_string(),
                start: 13400,
                end: 14550,
                strand: '+',
                description: "ISG15 ubiquitin-like modifier (interferon-stimulated innate immune effector)".to_string(),
            },
            GenomicAnnotation {
                feature_id: "REPEAT_VNTR_1P36".to_string(),
                feature_type: "repeat_region".to_string(),
                symbol: "VNTR-1p36-Subtelomeric".to_string(),
                start: 19500,
                end: 24800,
                strand: '+',
                description: "GC-rich subtelomeric tandem repeat array and G-quadruplex cluster".to_string(),
            },
            GenomicAnnotation {
                feature_id: "GENE_AGRN_N".to_string(),
                feature_type: "CDS".to_string(),
                symbol: "AGRN (5' exons 1-8)".to_string(),
                start: 20100,
                end: 36500,
                strand: '+',
                description: "Agrin neuromuscular junction proteoglycan (exons 1-8 overlapping GC-rich VNTR)".to_string(),
            },
            GenomicAnnotation {
                feature_id: "GENE_AGRN_C".to_string(),
                feature_type: "CDS".to_string(),
                symbol: "AGRN (3' exons 9-36)".to_string(),
                start: 36501,
                end: 56100,
                strand: '+',
                description: "Agrin laminin-G and EGF-like catalytic C-terminal domains".to_string(),
            },
        ];

        Ok((rec, annotations))
    }
}

pub fn parse_single_fasta(raw: &str) -> Option<(String, String)> {
    let mut header = String::new();
    let mut seq = String::with_capacity(raw.len());
    for line in raw.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('>') {
            if header.is_empty() {
                header = rest.trim().to_string();
            } else {
                break; // only take first record
            }
        } else {
            for ch in trimmed.chars() {
                let up = ch.to_ascii_uppercase();
                if matches!(up, 'A' | 'C' | 'G' | 'T') {
                    seq.push(up);
                } else if up.is_ascii_alphabetic() {
                    seq.push('N'); // Preserve IUPAC degenerate bases in-frame without frameshifting
                }
            }
        }
    }
    if seq.is_empty() {
        None
    } else {
        Some((header, seq))
    }
}

/// Sequence character normalization mode for shared multi-FASTA parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FastaSeqMode {
    /// Keep each trimmed non-empty sequence line verbatim (case and gap/symbol characters preserved).
    RawTrimmed,
    /// Uppercase each trimmed sequence line while preserving non-whitespace symbols.
    UppercaseTrimmed,
    /// Strip embedded whitespace inside sequence lines and convert characters to ASCII uppercase.
    UppercaseNonWhitespace,
    /// Retain only ASCII alphabetic characters and `-` gaps, converted to ASCII uppercase.
    UppercaseAlphaGap,
    /// Retain only ASCII alphabetic characters, converted to ASCII uppercase.
    UppercaseAlphaOnly,
}

/// Parses a single-record (or last-header concatenated) cached FASTA slice without altering sequence characters.
pub fn parse_single_fasta_raw(raw: &str) -> Option<(String, String)> {
    let mut desc = String::new();
    let mut seq = String::new();
    for line in raw.lines() {
        let t = line.trim();
        if let Some(rest) = t.strip_prefix('>') {
            desc = rest.to_string();
        } else if !t.is_empty() {
            seq.push_str(t);
        }
    }
    if seq.is_empty() {
        None
    } else {
        Some((desc, seq))
    }
}

/// Parses a multi-FASTA buffer into `(record_id, sequence)` pairs using the requested normalization mode.
///
/// `record_id` extracts the first whitespace-delimited token after `>` (falling back to `default_id`).
pub fn parse_multi_fasta_pairs(
    raw: &str,
    default_id: &str,
    mode: FastaSeqMode,
) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut cur_id = String::new();
    let mut cur_seq = String::new();
    for line in raw.lines() {
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        if let Some(rest) = t.strip_prefix('>') {
            if !cur_id.is_empty() {
                out.push((cur_id.clone(), cur_seq.clone()));
            }
            cur_id = rest.split_whitespace().next().unwrap_or(default_id).to_string();
            cur_seq.clear();
        } else {
            match mode {
                FastaSeqMode::RawTrimmed => cur_seq.push_str(t),
                FastaSeqMode::UppercaseTrimmed => {
                    for ch in t.chars() {
                        cur_seq.push(ch.to_ascii_uppercase());
                    }
                }
                FastaSeqMode::UppercaseNonWhitespace => {
                    for ch in t.chars() {
                        if !ch.is_whitespace() {
                            cur_seq.push(ch.to_ascii_uppercase());
                        }
                    }
                }
                FastaSeqMode::UppercaseAlphaGap => {
                    for b in t.bytes() {
                        let u = b.to_ascii_uppercase();
                        if u.is_ascii_alphabetic() || u == b'-' {
                            cur_seq.push(u as char);
                        }
                    }
                }
                FastaSeqMode::UppercaseAlphaOnly => {
                    for b in t.bytes() {
                        let u = b.to_ascii_uppercase();
                        if u.is_ascii_alphabetic() {
                            cur_seq.push(u as char);
                        }
                    }
                }
            }
        }
    }
    if !cur_id.is_empty() {
        out.push((cur_id, cur_seq));
    }
    out
}

/// Reads an aligned or unaligned multi-FASTA file keeping uppercase letters and `-` gap characters.
///
/// Panics with path context if the file cannot be read, matching the standalone audit contract.
pub fn read_multi_fasta_alpha_gap(path: impl AsRef<Path>) -> Vec<(String, String)> {
    let p = path.as_ref();
    let content = fs::read_to_string(p)
        .unwrap_or_else(|e| panic!("Failed reading multi-FASTA {}: {}", p.display(), e));
    parse_multi_fasta_pairs(&content, "SEQ", FastaSeqMode::UppercaseAlphaGap)
}

/// Reads a multi-FASTA file gracefully (returning empty on I/O error), uppercasing and stripping whitespace.
pub fn read_multi_fasta_lossy_upper(path: impl AsRef<Path>) -> Vec<(String, String)> {
    let raw = fs::read_to_string(path.as_ref()).unwrap_or_default();
    parse_multi_fasta_pairs(&raw, "", FastaSeqMode::UppercaseNonWhitespace)
}

/// Reads a multi-FASTA file into [`FastaRecord`] values using raw trimmed sequence lines.
pub fn load_fasta_records_raw(path: impl AsRef<Path>, family: &str, gene: &str) -> Vec<FastaRecord> {
    let content = fs::read_to_string(path.as_ref()).expect("Failed reading input FASTA");
    parse_multi_fasta_pairs(&content, "UNKNOWN", FastaSeqMode::RawTrimmed)
        .into_iter()
        .map(|(id, sequence)| FastaRecord {
            id: id.clone(),
            accession: id.clone(),
            taxon: id.clone(),
            common_name: id.clone(),
            family: family.to_string(),
            gene: gene.to_string(),
            description: id,
            sequence,
            cached_from_disk: true,
        })
        .collect()
}

#[cfg(test)]
mod fasta_parser_tests {
    use super::*;

    #[test]
    fn parses_single_and_raw_fasta() {
        let raw = ">NC_001234.1 Example taxon\nacgt ry\n>SECOND\ntttt\n";
        assert_eq!(
            parse_single_fasta(raw),
            Some(("NC_001234.1 Example taxon".to_string(), "ACGTNN".to_string()))
        );
        let single_cached = ">hdr one\nacgt\n>hdr two\nTGCA\n";
        assert_eq!(
            parse_single_fasta_raw(single_cached),
            Some(("hdr two".to_string(), "acgtTGCA".to_string()))
        );
        assert_eq!(parse_single_fasta_raw(">empty_only\n"), None);
    }

    #[test]
    fn parses_multi_fasta_all_modes() {
        let fasta = ">T001 first taxon\nacgt-?\n  tg ca \n>T002\nGgAa\n";
        assert_eq!(
            parse_multi_fasta_pairs(fasta, "SEQ", FastaSeqMode::RawTrimmed),
            vec![
                ("T001".to_string(), "acgt-?tg ca".to_string()),
                ("T002".to_string(), "GgAa".to_string()),
            ]
        );
        assert_eq!(
            parse_multi_fasta_pairs(fasta, "SEQ", FastaSeqMode::UppercaseTrimmed),
            vec![
                ("T001".to_string(), "ACGT-?TG CA".to_string()),
                ("T002".to_string(), "GGAA".to_string()),
            ]
        );
        assert_eq!(
            parse_multi_fasta_pairs(fasta, "SEQ", FastaSeqMode::UppercaseNonWhitespace),
            vec![
                ("T001".to_string(), "ACGT-?TGCA".to_string()),
                ("T002".to_string(), "GGAA".to_string()),
            ]
        );
        assert_eq!(
            parse_multi_fasta_pairs(fasta, "SEQ", FastaSeqMode::UppercaseAlphaGap),
            vec![
                ("T001".to_string(), "ACGT-TGCA".to_string()),
                ("T002".to_string(), "GGAA".to_string()),
            ]
        );
        assert_eq!(
            parse_multi_fasta_pairs(fasta, "SEQ", FastaSeqMode::UppercaseAlphaOnly),
            vec![
                ("T001".to_string(), "ACGTTGCA".to_string()),
                ("T002".to_string(), "GGAA".to_string()),
            ]
        );
    }
}
