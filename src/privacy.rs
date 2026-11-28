//! `privacy` — keep local machine details out of published artifacts.
//!
//! Release JSON and telemetry are published on crvscience.com and GitHub, so captured tool
//! output must not carry absolute local paths (which include the operator's user name).

use std::path::PathBuf;

/// Default working directory for audit/sweep binaries when none is passed on the command line:
/// `<system temp>/bio-suite/<sub>`. Keeps operator-specific paths out of source code.
pub fn default_scratch_dir(sub: &str) -> PathBuf {
    let base = std::env::temp_dir().join("bio-suite");
    if sub.is_empty() { base } else { base.join(sub) }
}

/// Root of this crate (the crvscience-labs checkout).
pub fn labs_repo_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Root of the CRV Science web checkout: `CRV_WEB_ROOT` if set, else the parent of crvscience-labs.
pub fn web_root_dir() -> PathBuf {
    std::env::var_os("CRV_WEB_ROOT")
        .map(PathBuf::from)
        .unwrap_or_else(|| labs_repo_dir().parent().map(PathBuf::from).unwrap_or_else(labs_repo_dir))
}

/// Repository folder names; an absolute path inside one of them is rewritten repo-relative.
const REPO_MARKERS: [&str; 2] = ["crvscience-web", "crvscience-labs"];

fn is_path_char(c: char) -> bool {
    !(c.is_whitespace() || matches!(c, '"' | '\'' | '`' | '<' | '>' | '|' | '(' | ')' | ',' | ';'))
}

/// Rewrite one absolute path token: repo-relative when it sits inside the repository,
/// otherwise just its final component.
fn rewrite_path(token: &str) -> String {
    let parts: Vec<&str> = token.split(['\\', '/']).filter(|p| !p.is_empty()).collect();
    if let Some(i) = parts
        .iter()
        .rposition(|p| REPO_MARKERS.iter().any(|m| p.eq_ignore_ascii_case(m)))
    {
        let rest = &parts[i + 1..];
        return if rest.is_empty() { ".".to_string() } else { rest.join("/") };
    }
    parts.last().map(|s| s.to_string()).unwrap_or_default()
}

/// Replace absolute local paths (`C:\Users\...`, `C:/...`, `/Users/...`, `/home/...`) in free text,
/// such as captured stdout/stderr of external tools.
pub fn redact_local_paths(text: &str) -> String {
    let chars: Vec<char> = text.chars().collect();
    let mut out = String::with_capacity(text.len());
    let mut i = 0;
    while i < chars.len() {
        let at_boundary = i == 0 || !is_path_char(chars[i - 1]) || chars[i - 1] == '=';
        // Drive letters can follow progress counters directly ("…8C:\Users"), so only a preceding
        // letter (as in "http:") rules a Windows path out.
        let windows_abs = (i == 0 || !chars[i - 1].is_ascii_alphabetic())
            && i + 2 < chars.len()
            && chars[i].is_ascii_alphabetic()
            && chars[i + 1] == ':'
            && (chars[i + 2] == '\\' || chars[i + 2] == '/');
        let unix_home = at_boundary
            && chars[i] == '/'
            && (text[char_offset(&chars, i)..].starts_with("/Users/") || text[char_offset(&chars, i)..].starts_with("/home/"));
        if windows_abs || unix_home {
            let start = i;
            while i < chars.len() && is_path_char(chars[i]) {
                i += 1;
            }
            let token: String = chars[start..i].iter().collect();
            out.push_str(&rewrite_path(&token));
            continue;
        }
        out.push(chars[i]);
        i += 1;
    }
    out
}

fn char_offset(chars: &[char], idx: usize) -> usize {
    chars[..idx].iter().map(|c| c.len_utf8()).sum()
}

/// Classification for terms found in downloaded FASTA/GenBank headers or metadata logs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MetadataContextClass {
    BenignScientificContext,
    ReviewRawHeader,
}

/// One detected metadata term with its benign-context classification.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MetadataContextFlag {
    pub term: String,
    pub class: MetadataContextClass,
}

/// Redact query-string credentials (`api_key=`, `token=`, `password=`, `email=`, `%40`) from log lines.
pub fn redact_url_credentials(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for (idx, line) in text.split_inclusive('\n').enumerate() {
        let _ = idx;
        let mut first = true;
        let has_nl = line.ends_with('\n');
        let trimmed_nl = line.strip_suffix('\n').unwrap_or(line);
        for tok in trimmed_nl.split(' ') {
            if !first {
                out.push(' ');
            }
            first = false;
            if let Some(q_idx) = tok.find(".fcgi?") {
                out.push_str(&tok[..q_idx + 5]);
                out.push_str("?...[query-redacted]");
                if tok.ends_with(':') {
                    out.push(':');
                }
            } else if tok.contains("api_key=")
                || tok.contains("token=")
                || tok.contains("password=")
                || tok.contains("email=")
                || tok.contains("%40")
            {
                out.push_str("[redacted-cred]");
            } else {
                out.push_str(tok);
            }
        }
        if has_nl {
            out.push('\n');
        }
    }
    out
}

/// Abbreviate binomial species labels in CLI stderr logs (e.g., `"Genus species"` -> `"G. species"`,
/// `"Candidatus Taxon"` -> `"Ca. Taxon"`) to prevent high-density raw taxon lists in terminal telemetry.
pub fn abbreviate_binomial_for_cli_log(taxon: &str) -> String {
    let parts: Vec<&str> = taxon.split_whitespace().collect();
    match parts.as_slice() {
        [] => String::new(),
        [single] => (*single).to_string(),
        [first, second, rest @ ..] if first.eq_ignore_ascii_case("candidatus") => {
            if rest.is_empty() {
                format!("Ca. {}", second)
            } else {
                format!("Ca. {} {}", second, rest.join(" "))
            }
        }
        [genus, epithet, rest @ ..] => {
            let g_init = genus.chars().next().unwrap_or('?');
            if rest.is_empty() {
                format!("{}. {}", g_init, epithet)
            } else {
                format!("{}. {} {}", g_init, epithet, rest.join(" "))
            }
        }
    }
}

const TRIGGER_TERMS: &[&str] = &[
    "toxin",
    "antitoxin",
    "hemolysin",
    "virulence",
    "pathogen",
    "lethal",
    "attack",
    "agent",
    "blast",
    "shell",
    "kill",
    "mycotoxin",
    "amatoxin",
    "amanitin",
    "aflatoxin",
    "puccinia",
    "magnaporthe",
    "ralstonia",
    "synchytrium",
    "batrachochytrium",
    "coccidioides",
    "cryptococcus",
    "conidiobolus",
    "mucormycosis",
];

const TRIGGER_SUBSTRINGS: &[&str] = &["toxicaria", "trichothec"];

const BENIGN_CONTEXT_CUES: &[&str] = &[
    "cds",
    "mrna",
    "gene",
    "protein",
    "subunit",
    "domain",
    "homolog",
    "phylogen",
    "plasmid",
    "replicon",
    "blastp",
    "blastn",
    "blastoderm",
    "super-blast",
    "super_blast",
    "reciprocal blast",
    "bit-score",
    "e-value",
    "shell matrix",
    "shell loss",
    "shell-less",
    "nacre",
    "mollus",
    "nudibranch",
    "gastropod",
    "polyplacophora",
    "8-valve",
    "dorsal",
    "valve",
    "chiton",
    "killifish",
    "biocontrol",
    "mitochondrial",
    "mitochondrion",
    "ribosomal",
    "rrna",
    "its1",
    "its2",
    "complete sequence",
    "complete genome",
    "insecticide",
    "pyrethroid",
    "diamide",
    "xenobiotic",
    "detoxification",
    "cyp",
    "cultivar",
    "quarantine",
    "tylenchomorpha",
    "heteroderidae",
    "meloidogyne",
    "entomopathogenic",
    "symbiont",
    "mycorrhiz",
    "ectomycorrhiz",
    "arbuscular",
    "glomeromycot",
    "agaricomycot",
    "basidiomycot",
    "ascomycot",
    "saprotroph",
    "endophyte",
    "hartig",
    "rhizosphere",
    "ortholog",
    "supermatrix",
    "phylo",
    "taxon",
    "clade",
    "cohort",
];

fn contains_ascii_word(haystack_lower: &str, needle: &str) -> bool {
    let bytes = haystack_lower.as_bytes();
    let nbytes = needle.as_bytes();
    if nbytes.is_empty() || bytes.len() < nbytes.len() {
        return false;
    }
    let mut i = 0;
    while i + nbytes.len() <= bytes.len() {
        if &bytes[i..i + nbytes.len()] == nbytes {
            let left_ok = i == 0 || !bytes[i - 1].is_ascii_alphanumeric();
            let right_ok = i + nbytes.len() == bytes.len() || !bytes[i + nbytes.len()].is_ascii_alphanumeric();
            if left_ok && right_ok {
                return true;
            }
        }
        i += 1;
    }
    false
}

/// Scan a FASTA header or sequence annotation string and classify trigger keywords as
/// `BenignScientificContext` when accompanied by standard molecular/taxonomic cues, or
/// `ReviewRawHeader` when raw free text lacks structured biological context.
pub fn classify_sequence_metadata_false_positives(header: &str) -> Vec<MetadataContextFlag> {
    let lower = header.to_ascii_lowercase();
    let has_benign_cue = BENIGN_CONTEXT_CUES.iter().any(|cue| lower.contains(cue));
    let class = if has_benign_cue {
        MetadataContextClass::BenignScientificContext
    } else {
        MetadataContextClass::ReviewRawHeader
    };
    let mut flags = Vec::new();
    for &term in TRIGGER_TERMS {
        if contains_ascii_word(&lower, term) {
            flags.push(MetadataContextFlag {
                term: term.to_string(),
                class,
            });
        }
    }
    for &sub in TRIGGER_SUBSTRINGS {
        if lower.contains(sub) {
            flags.push(MetadataContextFlag {
                term: sub.to_string(),
                class,
            });
        }
    }
    flags
}

/// Normalize a raw NCBI FASTA header into a concise, structured `accession | taxon | gene`
/// summary free of noisy free-text descriptions before emitting to logs or model prompts.
pub fn sanitize_fasta_header_for_log(accession: &str, taxon: &str, gene: &str, raw_desc: &str) -> String {
    let acc = accession.trim();
    let clean_acc = if !acc.is_empty() {
        acc.to_string()
    } else {
        raw_desc
            .trim_start_matches('>')
            .split_whitespace()
            .next()
            .unwrap_or("UNKNOWN")
            .to_string()
    };
    let clean_taxon: String = taxon
        .chars()
        .map(|c| if matches!(c, '\r' | '\n' | '|' | '<' | '>' | '`') { ' ' } else { c })
        .collect();
    let clean_gene: String = gene
        .chars()
        .map(|c| if matches!(c, '\r' | '\n' | '|' | '<' | '>' | '`') { ' ' } else { c })
        .collect();
    let t = clean_taxon.trim();
    let g = clean_gene.trim();
    if !t.is_empty() || !g.is_empty() {
        return format!(
            "{} | {} | {}",
            clean_acc,
            if t.is_empty() { "unclassified" } else { t },
            if g.is_empty() { "locus" } else { g }
        );
    }
    clean_acc
}

#[cfg(test)]
mod tests {
    use super::{
        abbreviate_binomial_for_cli_log, classify_sequence_metadata_false_positives,
        redact_local_paths, redact_url_credentials, sanitize_fasta_header_for_log,
        MetadataContextClass,
    };

    #[test]
    fn rewrites_repo_paths_and_hides_user_dirs() {
        assert_eq!(
            redact_local_paths(r"in C:\Users\someone\Documents\Code\CRVScience-web\labs/figures\a.svg done"),
            "in labs/figures/a.svg done"
        );
        assert_eq!(redact_local_paths(r"tmp=C:\Users\someone\AppData\Local\Temp\x.fa"), "tmp=x.fa");
        assert_eq!(redact_local_paths("/home/someone/crvscience-labs/runs/r1.log ok"), "runs/r1.log ok");
        assert_eq!(redact_local_paths("ratio 3:1 and http://a/b stays"), "ratio 3:1 and http://a/b stays");
        assert_eq!(
            redact_local_paths("/   8C:\\Users\\someone\\x\\crvscience-labs\\tools\\mafft.exe (nuc)"),
            "/   8tools/mafft.exe (nuc)"
        );
    }

    #[test]
    fn classifies_and_sanitizes_sequence_metadata() {
        let benign = classify_sequence_metadata_false_positives(
            ">NC_012345.1 Acanthochitona shell matrix protein mRNA, complete cds",
        );
        assert_eq!(benign.len(), 1);
        assert_eq!(benign[0].term, "shell");
        assert_eq!(benign[0].class, MetadataContextClass::BenignScientificContext);

        let expanded = classify_sequence_metadata_false_positives(
            "super-blast reciprocal blast bit-score nudibranch shell loss meloidogyne virulence",
        );
        assert_eq!(expanded.len(), 3);
        assert!(expanded
            .iter()
            .all(|f| f.class == MetadataContextClass::BenignScientificContext));

        let review = classify_sequence_metadata_false_positives(">XX_999999.1 unknown agent isolate");
        assert_eq!(review.len(), 1);
        assert_eq!(review[0].term, "agent");
        assert_eq!(review[0].class, MetadataContextClass::ReviewRawHeader);

        assert_eq!(
            sanitize_fasta_header_for_log("AB576780.1", "Aedes albopictus", "apyrase", ">AB576780.1 raw"),
            "AB576780.1 | Aedes albopictus | apyrase"
        );

        let fungal_benign = classify_sequence_metadata_false_positives(
            ">BK072141.1 TPA_asm: Amanita brunneitoxicaria mitochondrion, complete sequence",
        );
        assert_eq!(fungal_benign.len(), 1);
        assert_eq!(fungal_benign[0].term, "toxicaria");
        assert_eq!(fungal_benign[0].class, MetadataContextClass::BenignScientificContext);

        let rust_benign = classify_sequence_metadata_false_positives(
            ">AY125409.1 Puccinia graminis f. sp. tritici 18S ribosomal RNA gene, partial sequence",
        );
        assert_eq!(rust_benign.len(), 1);
        assert_eq!(rust_benign[0].term, "puccinia");
        assert_eq!(rust_benign[0].class, MetadataContextClass::BenignScientificContext);

        assert_eq!(abbreviate_binomial_for_cli_log("Puccinia graminis"), "P. graminis");
        assert_eq!(abbreviate_binomial_for_cli_log("Candidatus Heimdallarchaeota"), "Ca. Heimdallarchaeota");

        assert_eq!(
            redact_url_credentials("GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&email=labs%40crvscience.com: status 429"),
            "GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?...[query-redacted]: status 429"
        );
    }
}
