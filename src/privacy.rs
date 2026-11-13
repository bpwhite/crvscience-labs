//! `privacy` — keep local machine details out of published artifacts.
//!
//! Release JSON and telemetry are published on crvscience.com and GitHub, so captured tool
//! output must not carry absolute local paths (which include the operator's user name).

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

#[cfg(test)]
mod tests {
    use super::redact_local_paths;

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
}
