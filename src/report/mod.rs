//! Report output. v0.1 writes JSON only.

pub mod html;
pub mod json;
pub mod multiqc;
pub mod svg;
pub mod thresholds;

use std::io::Write;
use std::path::Path;

use crate::error::{IoContext, Result};

/// Writes `bytes` to `path` via a temporary file and a rename, so a reader
/// never sees a partial file and a failed run leaves nothing behind.
pub fn write_atomic(path: &Path, bytes: &[u8]) -> Result<()> {
    let mut tmp = path.as_os_str().to_owned();
    tmp.push(".tmp");
    let tmp = std::path::PathBuf::from(tmp);
    let res = (|| -> std::io::Result<()> {
        let mut f = std::fs::File::create(&tmp)?;
        f.write_all(bytes)?;
        f.sync_all()?;
        std::fs::rename(&tmp, path)
    })();
    if res.is_err() {
        let _ = std::fs::remove_file(&tmp);
    }
    res.ctx(format!("writing {}", path.display()))
}

pub fn with_suffix(prefix: &Path, suffix: &str) -> std::path::PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(suffix);
    std::path::PathBuf::from(s)
}

/// Writes the files rendered from a report: `<prefix>.qc.html` and the
/// MultiQC custom content. Returns the paths written.
pub fn write_rendered(r: &json::Report, prefix: &Path) -> Result<Vec<std::path::PathBuf>> {
    let html = with_suffix(prefix, ".qc.html");
    write_atomic(&html, html::render_html(r).as_bytes())?;
    let mut written = vec![html];
    for (suffix, v) in multiqc::files(r) {
        let p = with_suffix(prefix, suffix);
        let mut b = serde_json::to_vec_pretty(&v).expect("MultiQC JSON serialises");
        b.push(b'\n');
        write_atomic(&p, &b)?;
        written.push(p);
    }
    Ok(written)
}

/// Writes `<prefix>.qc.json` first, then everything rendered from it.
pub fn write_all(r: &json::Report, prefix: &Path) -> Result<Vec<std::path::PathBuf>> {
    let json = with_suffix(prefix, ".qc.json");
    if let Some(dir) = json.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir).ctx(format!("creating {}", dir.display()))?;
    }
    write_atomic(&json, &r.to_json())?;
    let mut written = vec![json];
    written.extend(write_rendered(r, prefix)?);
    Ok(written)
}

/// The output prefix for `render`: `--out` if given (a trailing `.qc.html`
/// or `.qc.json` is dropped), else the JSON path without `.qc.json`.
pub fn render_prefix(json: &Path, out: Option<&Path>) -> std::path::PathBuf {
    let strip = |p: &Path| {
        let s = p.to_string_lossy();
        for suffix in [".qc.json", ".qc.html", ".json", ".html"] {
            if let Some(stem) = s.strip_suffix(suffix) {
                return std::path::PathBuf::from(stem);
            }
        }
        p.to_path_buf()
    };
    strip(out.unwrap_or(json))
}

/// `warpqc render`: rebuilds the HTML and MultiQC files from a report. With
/// thresholds, re-evaluates `qc_status` and writes `<prefix>.qc.json` too.
pub fn render(
    json: &Path,
    thresholds: Option<&Path>,
    out: Option<&Path>,
) -> Result<Vec<std::path::PathBuf>> {
    let bytes = std::fs::read(json).ctx(format!("reading {}", json.display()))?;
    let mut report = json::Report::from_json(&bytes)?;
    let prefix = render_prefix(json, out);
    match thresholds {
        Some(t) => {
            let t = thresholds::Thresholds::load(t)?;
            let v = serde_json::to_value(&report).expect("report serialises");
            t.validate(&v)?;
            report.apply_thresholds(&t);
            write_all(&report, &prefix)
        }
        None => write_rendered(&report, &prefix),
    }
}
