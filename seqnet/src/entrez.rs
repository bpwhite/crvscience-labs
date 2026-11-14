//! Stage 1 Input Support: Rate-limited NCBI Entrez E-utilities client (`esearch.fcgi` +
//! `efetch.fcgi` INSDSeq XML) with `datetype=mdat` incremental delta watermarks, plus
//! local GenBank/FASTA parsing and deterministic offline Coleoptera benchmark generation.

use anyhow::{Context, Result};
use chrono::Utc;
use std::sync::Mutex;
use std::time::{Duration, Instant};

use crate::config::{ConsensusConfig, LocusClass, LocusSpec};
use crate::db::SequenceRecord;

const ESEARCH_URL: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi";
const EFETCH_URL: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";

/// Computes nucleotide composition statistics `(gc_fraction, ambig_fraction, shannon_entropy)`.
pub fn compute_seq_metrics(seq: &str) -> (f64, f64, f64) {
    if seq.is_empty() {
        return (0.0, 1.0, 0.0);
    }
    let mut counts = [0usize; 4]; // A, C, G, T
    let mut ambig = 0usize;
    for b in seq.bytes() {
        match b.to_ascii_uppercase() {
            b'A' => counts[0] += 1,
            b'C' => counts[1] += 1,
            b'G' => counts[2] += 1,
            b'T' | b'U' => counts[3] += 1,
            b'-' => ambig += 1,
            _ => ambig += 1,
        }
    }
    let total = seq.len() as f64;
    let unambig = (counts[0] + counts[1] + counts[2] + counts[3]).max(1) as f64;
    let gc = (counts[1] + counts[2]) as f64 / unambig;
    let ambig_frac = ambig as f64 / total;

    let mut entropy = 0.0;
    for &c in &counts {
        if c > 0 {
            let p = c as f64 / unambig;
            entropy -= p * p.log2();
        }
    }
    (gc, ambig_frac, entropy)
}

/// Checks whether a protein-coding sequence translates cleanly without internal stop codons
/// in at least one forward reading frame (`frame` in `0..3`) under NCBI genetic code 1 or 5.
pub fn check_reading_frame_clean(seq: &str, genetic_code_table: u8) -> bool {
    let clean: Vec<u8> = seq
        .bytes()
        .map(|b| b.to_ascii_uppercase())
        .filter(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
        .collect();
    if clean.len() < 30 {
        return false;
    }

    for frame in 0..3 {
        let mut internal_stops = 0usize;
        let slice = &clean[frame..];
        let n_codons = slice.len() / 3;
        if n_codons < 8 {
            continue;
        }
        // Ignore final codon (legitimate terminal stop codon)
        for c_idx in 0..(n_codons.saturating_sub(1)) {
            let codon = &slice[c_idx * 3..c_idx * 3 + 3];
            let is_stop = match genetic_code_table {
                // Table 5: Invertebrate Mitochondrial (TAA, TAG are stops; TGA = Trp, AGA/AGG = Ser)
                5 => codon == b"TAA" || codon == b"TAG",
                // Table 2: Vertebrate Mitochondrial (TAA, TAG, AGA, AGG)
                2 => codon == b"TAA" || codon == b"TAG" || codon == b"AGA" || codon == b"AGG",
                // Table 1: Standard Nuclear (TAA, TAG, TGA)
                _ => codon == b"TAA" || codon == b"TAG" || codon == b"TGA",
            };
            if is_stop {
                internal_stops += 1;
                break;
            }
        }
        if internal_stops == 0 {
            return true;
        }
    }
    false
}

/// Token-bucket polite NCBI Entrez HTTP client (`3 req/s` unauthenticated, `10 req/s` with API key).
pub struct EntrezClient {
    agent: ureq::Agent,
    api_key: Option<String>,
    email: String,
    min_interval: Duration,
    last_call: Mutex<Option<Instant>>,
}

impl EntrezClient {
    pub fn new(cfg: &ConsensusConfig) -> Self {
        let min_interval = if cfg.ncbi_api_key.is_some() {
            Duration::from_millis(105) // <= 9.5 req/s polite ceiling
        } else {
            Duration::from_millis(345) // <= 2.9 req/s polite ceiling
        };
        let agent = ureq::AgentBuilder::new()
            .timeout_connect(Duration::from_secs(10))
            .timeout_read(Duration::from_secs(30))
            .user_agent("bio-suite-seqnet/0.1 (https://github.com/bpwhite/crvscience-labs)")
            .build();

        Self {
            agent,
            api_key: cfg.ncbi_api_key.clone(),
            email: cfg.ncbi_email.clone(),
            min_interval,
            last_call: Mutex::new(None),
        }
    }

    fn wait_rate_limit(&self) {
        if let Ok(mut guard) = self.last_call.lock() {
            if let Some(prev) = *guard {
                let elapsed = prev.elapsed();
                if elapsed < self.min_interval {
                    std::thread::sleep(self.min_interval - elapsed);
                }
            }
            *guard = Some(Instant::now());
        }
    }

    /// Queries NCBI `nuccore` for a `(taxon_query, locus)` pair. If `since_mdat` (`YYYY/MM/DD`)
    /// is supplied from SQLite `entrez_sync_watermarks`, issues an incremental delta query
    /// (`datetype=mdat&mindate=...&maxdate=...`) so cached datasets only pull new accessions.
    /// Retries transient failures (timeouts, HTTP 429/5xx) up to 3 times with backoff.
    pub fn fetch_locus_delta(
        &self,
        taxon_query: &str,
        locus: &LocusSpec,
        since_mdat: Option<&str>,
        retmax: usize,
    ) -> Result<Vec<SequenceRecord>> {
        let term = build_locus_query(taxon_query, locus);
        self.with_retries(taxon_query, &locus.canonical_name, || {
            let accessions = self.esearch_accessions(&term, since_mdat, retmax)?;
            let mut records = Vec::new();
            for xml_body in self.efetch_insdseq(&accessions, 100)? {
                records.extend(parse_insdseq_xml(&xml_body, locus));
            }
            Ok(records)
        })
    }

    /// Mitogenome mode: one query per stratum for complete/partial mitochondrial genomes; every
    /// panel locus found in a record's CDS/rRNA feature table becomes its own `SequenceRecord`
    /// (same accession, different `locus_canonical`).
    pub fn fetch_mitogenome_delta(
        &self,
        taxon_query: &str,
        mito_loci: &[LocusSpec],
        since_mdat: Option<&str>,
        retmax: usize,
    ) -> Result<Vec<SequenceRecord>> {
        let term = build_mitogenome_query(taxon_query);
        self.with_retries(taxon_query, "mitogenome", || {
            let accessions = self.esearch_accessions(&term, since_mdat, retmax)?;
            let mut records = Vec::new();
            // Mitogenome XML is ~40 kB per record with features; keep POST bodies modest.
            for xml_body in self.efetch_insdseq(&accessions, 40)? {
                records.extend(parse_mitogenome_insdseq(&xml_body, mito_loci));
            }
            Ok(records)
        })
    }

    /// Retries transient failures (timeouts, HTTP 429/5xx) up to 3 times with backoff.
    fn with_retries<T>(&self, taxon_query: &str, what: &str, f: impl Fn() -> Result<T>) -> Result<T> {
        let mut last_err = None;
        for (attempt, wait_s) in [0u64, 2, 5, 12].into_iter().enumerate() {
            if wait_s > 0 {
                std::thread::sleep(Duration::from_secs(wait_s));
            }
            match f() {
                Ok(r) => return Ok(r),
                Err(e) => {
                    eprintln!(
                        "  [Stage 1] Entrez attempt {} for {} / {} failed: {:#}",
                        attempt + 1,
                        taxon_query,
                        what,
                        e
                    );
                    last_err = Some(e);
                }
            }
        }
        Err(last_err.unwrap_or_else(|| anyhow::anyhow!("Entrez request failed")))
    }

    fn esearch_accessions(&self, term: &str, since_mdat: Option<&str>, retmax: usize) -> Result<Vec<String>> {
        let today = Utc::now().format("%Y/%m/%d").to_string();
        self.wait_rate_limit();

        let mut req = self
            .agent
            .get(ESEARCH_URL)
            .query("db", "nuccore")
            .query("retmode", "json")
            .query("idtype", "acc")
            .query("retmax", &retmax.to_string())
            .query("tool", "bio-suite-seqnet")
            .query("email", &self.email)
            .query("term", term);

        if let Some(ref key) = self.api_key {
            req = req.query("api_key", key);
        }
        if let Some(min_d) = since_mdat {
            req = req
                .query("datetype", "mdat")
                .query("mindate", min_d)
                .query("maxdate", &today);
        }

        let esearch_json: serde_json::Value = req
            .call()
            .context("NCBI ESearch request failed")?
            .into_json()
            .context("failed to parse NCBI ESearch JSON")?;

        let id_list = esearch_json["esearchresult"]["idlist"]
            .as_array()
            .cloned()
            .unwrap_or_default();
        let accessions: Vec<String> = id_list
            .into_iter()
            .filter_map(|v| v.as_str().map(|s| s.to_string()))
            .collect();
        Ok(accessions)
    }

    /// EFetch INSDSeq XML bodies for `accessions`, `chunk` records per POST.
    fn efetch_insdseq(&self, accessions: &[String], chunk: usize) -> Result<Vec<String>> {
        let mut bodies = Vec::new();
        for chunk in accessions.chunks(chunk.max(1)) {
            self.wait_rate_limit();
            let id_csv = chunk.join(",");
            let freq = self.agent.post(EFETCH_URL);
            // `rettype=gbc` returns INSDSeq XML, which `parse_insdseq_xml` reads.
            // (`rettype=gb` + XML returns GBSeq tags and would parse to nothing.)
            let mut form_pairs = vec![
                ("db", "nuccore"),
                ("id", id_csv.as_str()),
                ("rettype", "gbc"),
                ("retmode", "xml"),
                ("tool", "bio-suite-seqnet"),
                ("email", self.email.as_str()),
            ];
            if let Some(ref key) = self.api_key {
                form_pairs.push(("api_key", key.as_str()));
            }

            let xml_body = freq
                .send_form(&form_pairs)
                .context("NCBI EFetch INSDSeq XML request failed")?
                .into_string()
                .context("failed reading EFetch XML body")?;

            if !xml_body.contains("<INSDSeq>") {
                anyhow::bail!(
                    "EFetch returned no INSDSeq records for {} requested accessions (unexpected response format)",
                    chunk.len()
                );
            }
            bodies.push(xml_body);
        }
        Ok(bodies)
    }
}

/// Entrez `nuccore` query for complete or partial mitochondrial genomes of one stratum.
/// GenBank records only: RefSeq `NC_` mitogenomes duplicate a GenBank submission.
pub fn build_mitogenome_query(taxon_query: &str) -> String {
    format!(
        "\"{}\"[Organism] AND mitochondrion[Filter] AND (\"complete genome\"[Title] OR \"partial genome\"[Title]) AND 10000:30000[SLEN] NOT srcdb_refseq[PROP]",
        taxon_query
    )
}

/// Maps a mitochondrial CDS/rRNA `gene` or `product` qualifier to a canonical panel locus.
pub fn mito_feature_canonical(name: &str) -> Option<&'static str> {
    let n: String = name
        .chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .collect::<String>()
        .to_ascii_lowercase();
    let n = n.strip_suffix("gene").unwrap_or(&n);
    let out = match n {
        "cox1" | "coi" | "co1" | "coxi" | "cytochromecoxidasesubuniti" | "cytochromecoxidasesubunit1"
        | "cytochromeoxidasesubuniti" | "cytochromeoxidasesubunit1" => "COI",
        "cox2" | "coii" | "co2" | "coxii" | "cytochromecoxidasesubunitii" | "cytochromecoxidasesubunit2"
        | "cytochromeoxidasesubunitii" | "cytochromeoxidasesubunit2" => "COII",
        "cox3" | "coiii" | "co3" | "coxiii" | "cytochromecoxidasesubunitiii" | "cytochromecoxidasesubunit3"
        | "cytochromeoxidasesubunitiii" | "cytochromeoxidasesubunit3" => "COIII",
        "cytb" | "cob" | "cyb" | "cytochromeb" | "apocytochromeb" => "CYTB",
        "nd1" | "nad1" | "nadhdehydrogenasesubunit1" | "nadhdehydrogenasesubuniti" => "ND1",
        "nd2" | "nad2" | "nadhdehydrogenasesubunit2" | "nadhdehydrogenasesubunitii" => "ND2",
        "nd3" | "nad3" | "nadhdehydrogenasesubunit3" | "nadhdehydrogenasesubunitiii" => "ND3",
        "nd4" | "nad4" | "nadhdehydrogenasesubunit4" | "nadhdehydrogenasesubunitiv" => "ND4",
        "nd4l" | "nad4l" | "nadhdehydrogenasesubunit4l" | "nadhdehydrogenasesubunitivl" => "ND4L",
        "nd5" | "nad5" | "nadhdehydrogenasesubunit5" | "nadhdehydrogenasesubunitv" => "ND5",
        "nd6" | "nad6" | "nadhdehydrogenasesubunit6" | "nadhdehydrogenasesubunitvi" => "ND6",
        "atp6" | "atpase6" | "atpsynthasef0subunit6" | "atpsynthasefosubunit6" | "atpsynthasesubunit6"
        | "atpasesubunit6" => "ATP6",
        "atp8" | "atpase8" | "atpsynthasef0subunit8" | "atpsynthasefosubunit8" | "atpsynthasesubunit8"
        | "atpasesubunit8" => "ATP8",
        "rrnl" | "16s" | "lrrna" | "16srrna" | "16sribosomalrna" | "largesubunitribosomalrna" | "lsurrna"
        | "rnr2" | "16sribosomalrnalargesubunit" => "rrnL_16S",
        "rrns" | "12s" | "srrna" | "12srrna" | "12sribosomalrna" | "smallsubunitribosomalrna" | "ssurrna"
        | "rnr1" | "12sribosomalrnasmallsubunit" => "rrnS_12S",
        _ => return None,
    };
    Some(out)
}

/// Slices CDS/rRNA features of mitogenome INSDSeq records into per-locus records (first feature
/// per locus per record; length bounds from the locus spec).
pub fn parse_mitogenome_insdseq(xml: &str, mito_loci: &[LocusSpec]) -> Vec<SequenceRecord> {
    let mut out = Vec::new();
    for block in xml.split("<INSDSeq>").skip(1) {
        let Some(end_idx) = block.find("</INSDSeq>") else {
            continue;
        };
        let rec_xml = &block[..end_idx];
        let Some(seq) = extract_xml_tag(rec_xml, "INSDSeq_sequence") else {
            continue;
        };
        let seq = seq.to_ascii_uppercase();
        let table = match (rec_xml.find("<INSDSeq_feature-table>"), rec_xml.find("</INSDSeq_feature-table>")) {
            (Some(a), Some(b)) if b > a => &rec_xml[a..b],
            _ => continue,
        };
        let mut seen: std::collections::HashSet<&'static str> = std::collections::HashSet::new();
        for feat in table.split("<INSDFeature>").skip(1) {
            let key = extract_xml_tag(feat, "INSDFeature_key").unwrap_or_default();
            if key != "CDS" && key != "rRNA" {
                continue;
            }
            let canon = ["gene", "product"]
                .iter()
                .filter_map(|q| extract_insd_qualifier(feat, q))
                .find_map(|v| mito_feature_canonical(&v));
            let Some(canon) = canon else { continue };
            let Some(spec) = mito_loci.iter().find(|l| l.canonical_name == canon) else {
                continue;
            };
            if (key == "CDS") != spec.protein_coding || seen.contains(canon) {
                continue;
            }
            let Some(loc) = extract_xml_tag(feat, "INSDFeature_location") else { continue };
            let Some(sub) = extract_location(seq.as_bytes(), &loc) else { continue };
            if sub.len() < spec.min_expected_bp || sub.len() > spec.max_expected_bp {
                continue;
            }
            if let Some(rec) = build_record(rec_xml, spec, sub) {
                seen.insert(canon);
                out.push(rec);
            }
        }
    }
    out
}

/// Extracts an INSDFeature location (`a..b`, `complement(..)`, `join(..)`/`order(..)`, partial
/// `<`/`>` ends) from a record sequence. Remote locations (`ACC:a..b`) return `None`.
pub fn extract_location(seq: &[u8], loc: &str) -> Option<String> {
    fn wrapper<'a>(s: &'a str, name: &str) -> Option<&'a str> {
        let rest = s.strip_prefix(name)?.strip_prefix('(')?;
        rest.strip_suffix(')')
    }
    fn split_top(s: &str) -> Vec<&str> {
        let (mut depth, mut start, mut parts) = (0i32, 0usize, Vec::new());
        for (i, c) in s.char_indices() {
            match c {
                '(' => depth += 1,
                ')' => depth -= 1,
                ',' if depth == 0 => {
                    parts.push(&s[start..i]);
                    start = i + 1;
                }
                _ => {}
            }
        }
        parts.push(&s[start..]);
        parts
    }
    fn inner(seq: &[u8], s: &str) -> Option<Vec<u8>> {
        let s = s.trim();
        if let Some(x) = wrapper(s, "complement") {
            return inner(seq, x).map(|v| revcomp(&v));
        }
        if let Some(x) = wrapper(s, "join").or_else(|| wrapper(s, "order")) {
            let mut out = Vec::new();
            for part in split_top(x) {
                out.extend(inner(seq, part)?);
            }
            return Some(out);
        }
        let clean: String = s.chars().filter(|c| !matches!(c, '<' | '>')).collect();
        let (a, b) = match clean.split_once("..") {
            Some((a, b)) => (a.trim().parse::<usize>().ok()?, b.trim().parse::<usize>().ok()?),
            None => {
                let p = clean.trim().parse::<usize>().ok()?;
                (p, p)
            }
        };
        if a == 0 || a > b || b > seq.len() {
            return None;
        }
        Some(seq[a - 1..b].to_vec())
    }
    if loc.contains(':') {
        return None;
    }
    inner(seq, &loc.replace(char::is_whitespace, "")).map(|v| String::from_utf8_lossy(&v).into_owned())
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' | b'U' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'R' => b'Y',
            b'Y' => b'R',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            other => other,
        })
        .collect()
}

/// Entrez `nuccore` query for one `(taxon, locus)` pair.
///
/// Synonyms are matched in record **titles** only: `[All Fields]` also hits feature notes and
/// annotations, which pulled RefSeq predicted mRNAs and mitochondrial rRNA into nuclear rRNA
/// loci. Predicted RefSeq models (`XM_`) are excluded, the genomic compartment is enforced with
/// the `mitochondrion` filter, and rRNA loci exclude mRNA records.
pub fn build_locus_query(taxon_query: &str, locus: &LocusSpec) -> String {
    let syn_clause = locus
        .entrez_synonyms
        .iter()
        .map(|s| format!("\"{}\"[Title]", s))
        .collect::<Vec<_>>()
        .join(" OR ");
    let compartment = match locus.locus_class {
        LocusClass::Mitochondrial => " AND mitochondrion[Filter]",
        LocusClass::Nuclear => " NOT mitochondrion[Filter]",
    };
    let molecule = if locus.protein_coding { "" } else { " NOT biomol_mrna[PROP]" };
    format!(
        "\"{}\"[Organism] AND ({}) AND {}:{}[SLEN]{}{} NOT srcdb_refseq_model[PROP] NOT \"whole genome shotgun\"[Title]",
        taxon_query, syn_clause, locus.min_expected_bp, locus.max_expected_bp, compartment, molecule
    )
}

/// Lightweight, allocation-conscious INSDSeq XML extractor for NCBI `efetch` GenBank XML.
pub fn parse_insdseq_xml(xml: &str, locus: &LocusSpec) -> Vec<SequenceRecord> {
    let mut out = Vec::new();
    for block in xml.split("<INSDSeq>").skip(1) {
        let Some(end_idx) = block.find("</INSDSeq>") else {
            continue;
        };
        let rec_xml = &block[..end_idx];
        let seq_raw = extract_xml_tag(rec_xml, "INSDSeq_sequence")
            .unwrap_or_default()
            .to_ascii_uppercase();
        if let Some(rec) = build_record(rec_xml, locus, seq_raw) {
            out.push(rec);
        }
    }
    out
}

/// Builds a `SequenceRecord` for `locus` from one INSDSeq record's metadata and `seq_raw`.
/// Returns `None` for missing metadata or unresolved (`sp.`, `cf.`, ...) species names.
fn build_record(rec_xml: &str, locus: &LocusSpec, seq_raw: String) -> Option<SequenceRecord> {
    {
        let acc = extract_xml_tag(rec_xml, "INSDSeq_accession-version")
            .or_else(|| extract_xml_tag(rec_xml, "INSDSeq_primary-accession"))
            .unwrap_or_default();
        let organism = extract_xml_tag(rec_xml, "INSDSeq_organism").unwrap_or_default();
        let taxonomy = extract_xml_tag(rec_xml, "INSDSeq_taxonomy").unwrap_or_default();
        if acc.is_empty() || organism.is_empty() || seq_raw.is_empty() {
            return None;
        }

        let parts: Vec<&str> = organism.split_whitespace().collect();
        if parts.len() < 2 {
            return None;
        }
        // Reject unverified `sp.`, `cf.`, `aff.`, or environmental placeholders at species level
        if matches!(
            parts[1].to_ascii_lowercase().as_str(),
            "sp." | "sp" | "cf." | "aff." | "gen." | "gen" | "nr." | "environmental" | "uncultured"
        ) {
            return None;
        }
        let genus = parts[0].to_string();
        let species_binomial = format!("{} {}", parts[0], parts[1]);

        let family = taxonomy
            .split(';')
            .map(|s| s.trim())
            .find(|t| t.ends_with("idae"))
            .unwrap_or("")
            .to_string();

        let voucher_id = extract_insd_qualifier(rec_xml, "specimen_voucher").unwrap_or_default();
        let isolate_code = extract_insd_qualifier(rec_xml, "isolate").unwrap_or_default();
        let country_locality = extract_insd_qualifier(rec_xml, "country")
            .or_else(|| extract_insd_qualifier(rec_xml, "geo_loc_name"))
            .unwrap_or_default();
        // A record can carry several db_xref qualifiers (e.g. BOLD:..., taxon:...); use the taxon one.
        let taxid = extract_all_insd_qualifiers(rec_xml, "db_xref")
            .into_iter()
            .find_map(|x| x.strip_prefix("taxon:").and_then(|n| n.parse::<u64>().ok()))
            .unwrap_or(0);

        let update_raw = extract_xml_tag(rec_xml, "INSDSeq_update-date").unwrap_or_default();
        let ncbi_mod_date = normalize_genbank_date(&update_raw);

        let (gc, ambig, entropy) = compute_seq_metrics(&seq_raw);
        let rf_ok = if locus.protein_coding {
            check_reading_frame_clean(&seq_raw, locus.genetic_code_table)
        } else {
            true
        };

        Some(SequenceRecord {
            accession_version: acc,
            ncbi_taxid: taxid,
            organism_name: organism,
            genus,
            species_binomial,
            family,
            voucher_id,
            isolate_code,
            country_locality,
            locus_canonical: locus.canonical_name.clone(),
            locus_class: locus.locus_class,
            seq_length: seq_raw.len(),
            sequence_raw: seq_raw,
            gc_fraction: gc,
            ambig_fraction: ambig,
            shannon_entropy: entropy,
            reading_frame_ok: rf_ok,
            ncbi_mod_date,
        })
    }
}

fn extract_xml_tag(xml: &str, tag: &str) -> Option<String> {
    let open = format!("<{}>", tag);
    let close = format!("</{}>", tag);
    let start = xml.find(&open)? + open.len();
    let end = xml[start..].find(&close)? + start;
    let raw = xml[start..end].trim();
    if !raw.contains('&') {
        return Some(raw.to_string());
    }
    // XML entities: partial feature ends arrive as `&lt;1487..3019`, names may carry `&amp;`.
    Some(
        raw.replace("&lt;", "<")
            .replace("&gt;", ">")
            .replace("&quot;", "\"")
            .replace("&apos;", "'")
            .replace("&amp;", "&"),
    )
}

fn extract_insd_qualifier(xml: &str, qual_name: &str) -> Option<String> {
    let needle = format!("<INSDQualifier_name>{}</INSDQualifier_name>", qual_name);
    let pos = xml.find(&needle)? + needle.len();
    extract_xml_tag(&xml[pos..], "INSDQualifier_value")
}

fn extract_all_insd_qualifiers(xml: &str, qual_name: &str) -> Vec<String> {
    let needle = format!("<INSDQualifier_name>{}</INSDQualifier_name>", qual_name);
    let mut out = Vec::new();
    let mut rest = xml;
    while let Some(pos) = rest.find(&needle) {
        rest = &rest[pos + needle.len()..];
        if let Some(v) = extract_xml_tag(rest, "INSDQualifier_value") {
            out.push(v);
        }
    }
    out
}

fn normalize_genbank_date(gb_date: &str) -> String {
    // Converts `14-OCT-2026` -> `2026/10/14`; falls back to current UTC date.
    let segs: Vec<&str> = gb_date.split('-').collect();
    if segs.len() == 3 {
        let month = match segs[1].to_ascii_uppercase().as_str() {
            "JAN" => "01",
            "FEB" => "02",
            "MAR" => "03",
            "APR" => "04",
            "MAY" => "05",
            "JUN" => "06",
            "JUL" => "07",
            "AUG" => "08",
            "SEP" => "09",
            "OCT" => "10",
            "NOV" => "11",
            "DEC" => "12",
            _ => "01",
        };
        return format!("{}/{}/{}", segs[2], month, segs[0]);
    }
    Utc::now().format("%Y/%m/%d").to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::build_named_locus_panel;

    /// INSDSeq as returned by `efetch rettype=gbc retmode=xml` (trimmed real-record shape).
    const FIXTURE: &str = r#"<INSDSet><INSDSeq>
<INSDSeq_locus>AF398722</INSDSeq_locus><INSDSeq_length>40</INSDSeq_length>
<INSDSeq_update-date>14-OCT-2026</INSDSeq_update-date>
<INSDSeq_primary-accession>AF398722</INSDSeq_primary-accession>
<INSDSeq_accession-version>AF398722.1</INSDSeq_accession-version>
<INSDSeq_organism>Carabus nemoralis</INSDSeq_organism>
<INSDSeq_taxonomy>Eukaryota; Metazoa; Arthropoda; Insecta; Coleoptera; Adephaga; Caraboidea; Carabidae; Carabinae; Carabus</INSDSeq_taxonomy>
<INSDSeq_feature-table><INSDFeature><INSDFeature_quals>
<INSDQualifier><INSDQualifier_name>specimen_voucher</INSDQualifier_name><INSDQualifier_value>CRV 12</INSDQualifier_value></INSDQualifier>
<INSDQualifier><INSDQualifier_name>db_xref</INSDQualifier_name><INSDQualifier_value>BOLD:AAA1234.12S</INSDQualifier_value></INSDQualifier>
<INSDQualifier><INSDQualifier_name>db_xref</INSDQualifier_name><INSDQualifier_value>taxon:41073</INSDQualifier_value></INSDQualifier>
</INSDFeature_quals></INSDFeature></INSDSeq_feature-table>
<INSDSeq_sequence>acgtacgtacgtacgtacgtacgtacgtacgtacgtacgt</INSDSeq_sequence>
</INSDSeq>
<INSDSeq><INSDSeq_accession-version>XX000001.1</INSDSeq_accession-version>
<INSDSeq_organism>Carabus sp. ABC</INSDSeq_organism><INSDSeq_sequence>acgt</INSDSeq_sequence></INSDSeq></INSDSet>"#;

    #[test]
    fn parses_insdseq_and_finds_taxon_db_xref() {
        let panel = build_named_locus_panel(&["SSU_18S".to_string()], true);
        let recs = parse_insdseq_xml(FIXTURE, &panel[0]);
        assert_eq!(recs.len(), 1, "unresolved `sp.` names are rejected");
        let r = &recs[0];
        assert_eq!(r.accession_version, "AF398722.1");
        assert_eq!(r.species_binomial, "Carabus nemoralis");
        assert_eq!(r.family, "Carabidae");
        assert_eq!(r.voucher_id, "CRV 12");
        assert_eq!(r.ncbi_taxid, 41073, "taxon: db_xref found even after a BOLD db_xref");
        assert_eq!(r.ncbi_mod_date, "2026/10/14");
        assert_eq!(r.sequence_raw, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
    }

    #[test]
    fn locus_queries_are_title_restricted_and_compartment_filtered() {
        let panel = build_named_locus_panel(&["SSU_18S".into(), "COI".into()], true);
        let q18 = build_locus_query("Carabidae", &panel[0]);
        assert!(q18.contains("\"18S\"[Title]") && !q18.contains("[All Fields]"));
        assert!(!q18.contains("\"SSU\""), "bare SSU also names mitochondrial 12S");
        assert!(q18.contains("NOT mitochondrion[Filter]") && q18.contains("NOT biomol_mrna[PROP]"));
        assert!(q18.contains("NOT srcdb_refseq_model[PROP]"));
        let qcoi = build_locus_query("Carabidae", &panel[1]);
        assert!(qcoi.contains("AND mitochondrion[Filter]") && !qcoi.contains("biomol_mrna"));
    }

    #[test]
    fn feature_locations_cover_complement_join_and_partials() {
        let seq = b"AACCGGTTAC";
        assert_eq!(extract_location(seq, "3..6").as_deref(), Some("CCGG"));
        assert_eq!(extract_location(seq, "<1..>4").as_deref(), Some("AACC"));
        // complement(7..10) = revcomp("TTAC") = "GTAA"
        assert_eq!(extract_location(seq, "complement(7..10)").as_deref(), Some("GTAA"));
        assert_eq!(extract_location(seq, "join(9..10,1..2)").as_deref(), Some("ACAA"));
        assert_eq!(extract_location(seq, "complement(join(1..2,9..10))").as_deref(), Some("GTTT"));
        assert_eq!(extract_location(seq, "X12345.1:1..4"), None);
        assert_eq!(extract_location(seq, "5..11"), None, "out of range");
    }

    #[test]
    fn mito_feature_names_map_to_panel_loci() {
        assert_eq!(mito_feature_canonical("COX1"), Some("COI"));
        assert_eq!(mito_feature_canonical("cytochrome c oxidase subunit I"), Some("COI"));
        assert_eq!(mito_feature_canonical("nad4L"), Some("ND4L"));
        assert_eq!(mito_feature_canonical("ND4"), Some("ND4"));
        assert_eq!(mito_feature_canonical("ATP synthase F0 subunit 8"), Some("ATP8"));
        assert_eq!(mito_feature_canonical("l-rRNA"), Some("rrnL_16S"));
        assert_eq!(mito_feature_canonical("12S ribosomal RNA"), Some("rrnS_12S"));
        assert_eq!(mito_feature_canonical("tRNA-Leu"), None);
    }

    #[test]
    fn mitogenome_record_is_sliced_into_panel_loci() {
        // 400 bp toy "mitogenome": a forward CDS (ATP8, 150 bp) and a reverse-strand rRNA (12S, 330 bp
        // would exceed the record, so use a short 12S span and widen bounds via a custom spec).
        let codon_body: String = "ATGAAA".repeat(25); // 150 bp, no TAA/TAG
        let filler: String = "C".repeat(250);
        let seq = format!("{}{}", codon_body, filler).to_ascii_lowercase();
        let xml = format!(
            r#"<INSDSet><INSDSeq><INSDSeq_accession-version>MT000001.1</INSDSeq_accession-version>
<INSDSeq_organism>Tenomerga mucida</INSDSeq_organism>
<INSDSeq_taxonomy>Eukaryota; Metazoa; Arthropoda; Insecta; Coleoptera; Archostemata; Cupedidae; Tenomerga</INSDSeq_taxonomy>
<INSDSeq_feature-table>
<INSDFeature><INSDFeature_key>source</INSDFeature_key><INSDFeature_location>1..400</INSDFeature_location><INSDFeature_quals>
<INSDQualifier><INSDQualifier_name>db_xref</INSDQualifier_name><INSDQualifier_value>taxon:999</INSDQualifier_value></INSDQualifier></INSDFeature_quals></INSDFeature>
<INSDFeature><INSDFeature_key>gene</INSDFeature_key><INSDFeature_location>1..150</INSDFeature_location><INSDFeature_quals>
<INSDQualifier><INSDQualifier_name>gene</INSDQualifier_name><INSDQualifier_value>ATP8</INSDQualifier_value></INSDQualifier></INSDFeature_quals></INSDFeature>
<INSDFeature><INSDFeature_key>CDS</INSDFeature_key><INSDFeature_location>&lt;1..150</INSDFeature_location><INSDFeature_quals>
<INSDQualifier><INSDQualifier_name>gene</INSDQualifier_name><INSDQualifier_value>ATP8</INSDQualifier_value></INSDQualifier></INSDFeature_quals></INSDFeature>
<INSDFeature><INSDFeature_key>rRNA</INSDFeature_key><INSDFeature_location>complement(51..400)</INSDFeature_location><INSDFeature_quals>
<INSDQualifier><INSDQualifier_name>product</INSDQualifier_name><INSDQualifier_value>12S ribosomal RNA</INSDQualifier_value></INSDQualifier></INSDFeature_quals></INSDFeature>
<INSDFeature><INSDFeature_key>tRNA</INSDFeature_key><INSDFeature_location>151..170</INSDFeature_location><INSDFeature_quals>
<INSDQualifier><INSDQualifier_name>product</INSDQualifier_name><INSDQualifier_value>tRNA-Leu</INSDQualifier_value></INSDQualifier></INSDFeature_quals></INSDFeature>
</INSDSeq_feature-table>
<INSDSeq_sequence>{seq}</INSDSeq_sequence></INSDSeq></INSDSet>"#
        );
        let panel = build_named_locus_panel(&["ATP8".into(), "rrnS_12S".into(), "COI".into()], true);
        let recs = parse_mitogenome_insdseq(&xml, &panel);
        assert_eq!(recs.len(), 2, "ATP8 CDS + 12S rRNA; `gene` and tRNA features ignored");
        let atp8 = recs.iter().find(|r| r.locus_canonical == "ATP8").unwrap();
        assert_eq!(atp8.seq_length, 150);
        assert!(atp8.reading_frame_ok);
        assert_eq!(atp8.family, "Cupedidae");
        assert_eq!(atp8.ncbi_taxid, 999);
        let s12 = recs.iter().find(|r| r.locus_canonical == "rrnS_12S").unwrap();
        assert_eq!(s12.seq_length, 350);
        assert!(s12.sequence_raw.starts_with("GGGG"), "reverse strand is reverse-complemented");
        assert!(build_mitogenome_query("Cupedidae").contains("NOT srcdb_refseq[PROP]"));
    }

    #[test]
    fn gbseq_xml_parses_to_nothing() {
        // Guards the rettype choice: GBSeq tags (rettype=gb) are not INSDSeq.
        let panel = build_named_locus_panel(&["SSU_18S".to_string()], true);
        let gbseq = FIXTURE.replace("INSDSeq", "GBSeq").replace("INSDSet", "GBSet");
        assert!(parse_insdseq_xml(&gbseq, &panel[0]).is_empty());
    }
}

/// Deterministic offline Coleoptera multi-locus sequence benchmark generator for `demo-sim`
/// and integration testing (simulates realistic vouchers, chimeric coverage across specimens,
/// NUMT internal-stop pseudogenes, Diptera sample-swaps, high-N low-quality reads, and
/// a singleton genus that fails the `>= 2 species/genus` constraint).
pub fn generate_coleoptera_benchmark_cohort(
    panel: &[LocusSpec],
    include_delta_wave: bool,
) -> Vec<SequenceRecord> {
    // 6 major Coleoptera families, genera with >=2 species + 1 singleton genus (`Cupes`)
    let taxa: &[(&str, &str, &str, u64)] = &[
        ("Carabidae", "Carabus", "Carabus nemoralis", 41073),
        ("Carabidae", "Carabus", "Carabus granulatus", 41074),
        ("Carabidae", "Carabus", "Carabus auratus", 41075),
        ("Carabidae", "Calosoma", "Calosoma scrutator", 49412),
        ("Carabidae", "Calosoma", "Calosoma sycophanta", 49413),
        ("Scarabaeidae", "Onthophagus", "Onthophagus taurus", 166361),
        ("Scarabaeidae", "Onthophagus", "Onthophagus nuchicornis", 166362),
        ("Coccinellidae", "Coccinella", "Coccinella septempunctata", 41139),
        ("Coccinellidae", "Coccinella", "Coccinella transversoguttata", 41140),
        ("Chrysomelidae", "Leptinotarsa", "Leptinotarsa decemlineata", 7539),
        ("Chrysomelidae", "Leptinotarsa", "Leptinotarsa juncta", 7540),
        ("Curculionidae", "Dendroctonus", "Dendroctonus ponderosae", 77166),
        ("Curculionidae", "Dendroctonus", "Dendroctonus frontalis", 77167),
        ("Cerambycidae", "Anoplophora", "Anoplophora glabripennis", 217634),
        ("Cerambycidae", "Anoplophora", "Anoplophora chinensis", 217635),
        // Singleton genus (only 1 species -> should be filtered out by `min_species_per_genus = 2`!)
        ("Cupedidae", "Cupes", "Cupes capitatus", 110815),
    ];

    let mut records = Vec::new();
    let bases = [b'A', b'C', b'G', b'T'];

    for (l_idx, locus) in panel.iter().enumerate() {
        let len_codons = ((locus.min_expected_bp + 90) / 3).max(90);
        // Deterministic ancestral template without stop codons (`CAG`, `CTG`, `ATC`, `GTC`, `AAC`, `GGT`)
        let safe_codons: [&[u8; 3]; 8] = [
            b"CAG", b"CTG", b"ATC", b"GTC", b"AAC", b"GGT", b"TCC", b"GAC",
        ];
        let mut ancestral = Vec::with_capacity(len_codons * 3);
        for c in 0..len_codons {
            let pick = safe_codons[(c * 7 + l_idx * 13) % safe_codons.len()];
            ancestral.extend_from_slice(pick);
        }

        for (t_idx, (family, genus, sp_bin, taxid)) in taxa.iter().enumerate() {
            // Leave 1 missing locus on `Anoplophora chinensis` until `include_delta_wave` if l_idx < 3
            if !include_delta_wave && *sp_bin == "Anoplophora chinensis" && l_idx >= panel.len().saturating_sub(4) {
                continue;
            }
            // Drop at most 1 locus per species so clean chimeras exceed 90% (24/25 = 96%)
            if (t_idx + l_idx * 5) % 29 == 0 && l_idx > 0 {
                continue;
            }

            // Split coverage across two museum specimens (`VOUCHER-A` covers ~65%, `VOUCHER-B` ~45%)
            // so the neural network constructs multi-specimen chimeric species terminals!
            let specimen_slots: &[(&str, &str)] = if l_idx % 3 == 0 {
                &[("CRV-COL-A", "2026/08/15"), ("CRV-COL-B", "2026/09/01")]
            } else if l_idx % 2 == 0 {
                &[("CRV-COL-A", "2026/08/15")]
            } else {
                &[("CRV-COL-B", "2026/09/01")]
            };

            for (s_idx, (v_prefix, mod_date)) in specimen_slots.iter().enumerate() {
                let mut seq = ancestral.clone();
                // Introduce realistic hierarchical phylogenetic substitutions (3rd-codon synonymous changes)
                for c in 0..len_codons {
                    if (c + t_idx / 2) % 19 == 0 {
                        // Genus-shared third-position transition (C <-> T)
                        let pos = c * 3 + 2;
                        seq[pos] = if seq[pos] == b'C' { b'T' } else { b'C' };
                    }
                    if (c * 3 + t_idx) % 41 == 0 {
                        let pos = c * 3 + 2;
                        seq[pos] = if seq[pos] == b'G' { b'A' } else { b'C' };
                    }
                    if s_idx == 1 && c == 12 {
                        // Tiny intraspecific 1-bp polymorphism between VOUCHER-A and VOUCHER-B
                        let pos = c * 3 + 2;
                        seq[pos] = bases[(t_idx + l_idx) % 4];
                    }
                }

                let seq_str = String::from_utf8(seq).unwrap_or_default();
                let (gc, ambig, entropy) = compute_seq_metrics(&seq_str);
                let rf_ok = if locus.protein_coding {
                    check_reading_frame_clean(&seq_str, locus.genetic_code_table)
                } else {
                    true
                };

                let acc = format!("MN{:02}{:02}{:02}.1", l_idx + 1, t_idx + 1, s_idx + 1);
                let voucher = format!("{}-{:03}", v_prefix, t_idx + 1);
                records.push(SequenceRecord {
                    accession_version: acc,
                    ncbi_taxid: *taxid,
                    organism_name: (*sp_bin).to_string(),
                    genus: (*genus).to_string(),
                    species_binomial: (*sp_bin).to_string(),
                    family: (*family).to_string(),
                    voucher_id: voucher,
                    isolate_code: format!("ISO-{:02}", s_idx + 1),
                    country_locality: "USA: Oregon, Corvallis".to_string(),
                    locus_canonical: locus.canonical_name.clone(),
                    locus_class: locus.locus_class,
                    seq_length: seq_str.len(),
                    sequence_raw: seq_str,
                    gc_fraction: gc,
                    ambig_fraction: ambig,
                    shannon_entropy: entropy,
                    reading_frame_ok: rf_ok,
                    ncbi_mod_date: (*mod_date).to_string(),
                });
            }
        }

        // Inject realistic bad GenBank accessions for Stage 2 Processing Neurons to quarantine:
        // 1. NUMT pseudogene with internal `TAA` stop codons in `COI` (`Carabus nemoralis`)
        if locus.canonical_name == "COI" {
            let mut numt = ancestral.clone();
            if numt.len() > 90 {
                numt[30..33].copy_from_slice(b"TAA");
                numt[60..63].copy_from_slice(b"TAG");
                numt[75..78].copy_from_slice(b"TAA");
            }
            let s = String::from_utf8(numt).unwrap_or_default();
            let (gc, ambig, entropy) = compute_seq_metrics(&s);
            records.push(SequenceRecord {
                accession_version: "KX990001.1".to_string(),
                ncbi_taxid: 41073,
                organism_name: "Carabus nemoralis".to_string(),
                genus: "Carabus".to_string(),
                species_binomial: "Carabus nemoralis".to_string(),
                family: "Carabidae".to_string(),
                voucher_id: "BAD-NUMT-01".to_string(),
                isolate_code: "NUMT1".to_string(),
                country_locality: "Unknown".to_string(),
                locus_canonical: locus.canonical_name.clone(),
                locus_class: LocusClass::Mitochondrial,
                seq_length: s.len(),
                sequence_raw: s,
                gc_fraction: gc,
                ambig_fraction: ambig,
                shannon_entropy: entropy,
                reading_frame_ok: false,
                ncbi_mod_date: "2026/05/10".to_string(),
            });
        }

        // 2. Cross-phylum sample swap / extreme K2P outlier in `EF1a` (`Coccinella septempunctata`)
        if locus.canonical_name == "EF1a" || l_idx == 0 {
            let mut swapped = Vec::with_capacity(ancestral.len());
            for i in 0..ancestral.len() {
                swapped.push(bases[(i * 11 + 3) % 4]);
            }
            // Strip accidental stops so specifically the homology/K2P neuron catches the swap!
            for chunk in swapped.chunks_exact_mut(3) {
                if chunk == b"TAA" || chunk == b"TAG" || chunk == b"TGA" {
                    chunk.copy_from_slice(b"CAG");
                }
            }
            let s = String::from_utf8(swapped).unwrap_or_default();
            let (gc, ambig, entropy) = compute_seq_metrics(&s);
            records.push(SequenceRecord {
                accession_version: format!("JQ8800{:02}.1", l_idx + 1),
                ncbi_taxid: 41139,
                organism_name: "Coccinella septempunctata".to_string(),
                genus: "Coccinella".to_string(),
                species_binomial: "Coccinella septempunctata".to_string(),
                family: "Coccinellidae".to_string(),
                voucher_id: "SWAP-CONTAM-99".to_string(),
                isolate_code: "ERR99".to_string(),
                country_locality: "Lab Contaminant".to_string(),
                locus_canonical: locus.canonical_name.clone(),
                locus_class: locus.locus_class,
                seq_length: s.len(),
                sequence_raw: s,
                gc_fraction: gc,
                ambig_fraction: ambig,
                shannon_entropy: entropy,
                reading_frame_ok: true,
                ncbi_mod_date: "2026/06/12".to_string(),
            });
        }

        // 3. High-ambiguity degraded read (>18% 'N's) in `LSU_28S`
        if l_idx == 1 {
            let mut degraded = ancestral.clone();
            for i in (0..degraded.len()).step_by(5) {
                degraded[i] = b'N';
            }
            let s = String::from_utf8(degraded).unwrap_or_default();
            let (gc, ambig, entropy) = compute_seq_metrics(&s);
            records.push(SequenceRecord {
                accession_version: "AY770012.1".to_string(),
                ncbi_taxid: 7539,
                organism_name: "Leptinotarsa decemlineata".to_string(),
                genus: "Leptinotarsa".to_string(),
                species_binomial: "Leptinotarsa decemlineata".to_string(),
                family: "Chrysomelidae".to_string(),
                voucher_id: "LOWQUAL-N12".to_string(),
                isolate_code: "LQ12".to_string(),
                country_locality: "USA".to_string(),
                locus_canonical: locus.canonical_name.clone(),
                locus_class: locus.locus_class,
                seq_length: s.len(),
                sequence_raw: s,
                gc_fraction: gc,
                ambig_fraction: ambig,
                shannon_entropy: entropy,
                reading_frame_ok: true,
                ncbi_mod_date: "2026/04/02".to_string(),
            });
        }
    }

    records
}
