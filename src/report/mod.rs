//! Report output. v0.1 writes JSON only.

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
