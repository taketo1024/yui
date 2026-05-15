//! Runtime resolution of the user-data directory for yui.
//!
//! Lookup order:
//!   1. `$YUI_DATA_DIR` if set.
//!   2. Platform user-data dir via the `directories` crate
//!      (e.g. `~/Library/Application Support/yui/` on macOS).
//!
//! Subdirectories under the data dir partition the data by kind, e.g.
//! `<data_dir>/links/3_1.json`, `<data_dir>/braid/3_1.json`. The set of
//! kinds is open — callers pass the subdirectory name they want.

use std::path::PathBuf;

use directories::ProjectDirs;

pub const ENV_VAR: &str = "YUI_DATA_DIR";

pub fn resolve_data_dir() -> Result<PathBuf, String> {
    if let Some(dir) = std::env::var_os(ENV_VAR) {
        return Ok(PathBuf::from(dir));
    }
    let proj = ProjectDirs::from("", "", "yui")
        .ok_or_else(|| "no platform user-data directory available".to_string())?;
    Ok(proj.data_dir().to_path_buf())
}

pub fn load_json(kind: &str, name: &str) -> Result<String, Box<dyn std::error::Error>> {
    let dir = resolve_data_dir()?;
    let path = dir.join(kind).join(format!("{name}.json"));
    if !path.exists() {
        return Err(format!(
            "no `{kind}/{name}.json` under {}. \
             Set ${ENV_VAR} or run scripts/fetch-knotdata.sh to populate it.",
            dir.display()
        ).into());
    }
    Ok(std::fs::read_to_string(&path)?)
}
