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
