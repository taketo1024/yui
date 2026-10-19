//! `setup`: download the knot data into the data directory.

use std::error::Error;
use std::fs;
use std::io::{Cursor, Read};
use std::path::{Path, PathBuf};

use calamine::{Reader, Xls, open_workbook_from_rs};
use serde_json::Value;
use yui_core::util::data_dir::resolve_data_dir;
use zip::ZipArchive;

use crate::app::args::LogArg;

const KNOTINFO_URL: &str = "https://raw.githubusercontent.com/knotinfo-org/assets/main/knotinfo_data_complete.xls.zip";
const KNOTINFO_XLS: &str = "knotinfo_data_complete.xls";
const INV_LINK_URL: &str = "https://raw.githubusercontent.com/taketo1024/yui/setup/resources/inv_link";
const INV_LINK_ZIPS: [&str; 3] = ["A.zip", "B.zip", "two_bridge.zip"];
const MAX_DOWNLOAD: u64 = 256 << 20;

#[derive(Clone, Debug, clap::Args)]
pub struct Args {
    /// Output directory (default: $YUI_DATA_DIR, else the platform user-data dir)
    #[arg(long)]
    pub out: Option<PathBuf>,

    /// Read this local KnotInfo file (.xls or .xls.zip) instead of downloading
    #[arg(long = "in", value_name = "FILE")]
    pub input: Option<PathBuf>,

    /// Overwrite existing files
    #[arg(long)]
    pub force: bool,

    /// Remove `inv_link/` before extracting, so diagrams dropped from the bundle do not linger
    #[arg(long)]
    pub clean: bool,

    #[arg(long, default_value = "0")]
    pub log: u8,
}

impl LogArg for Args {
    fn log_value(&self) -> u8 { self.log }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Outcome {
    Written,
    Kept,
    Skipped,
}

pub fn dispatch(args: &Args) -> Result<String, Box<dyn Error>> {
    let out_dir = match &args.out {
        Some(dir) => dir.clone(),
        None => resolve_data_dir()?,
    };
    let xls = match &args.input {
        Some(path) => read_local(path)?,
        None => unzip_xls(&download(KNOTINFO_URL)?)?,
    };

    let knots = save_knots(xls, &out_dir, args.force)?;
    let inv_links = save_inv_links(&out_dir, args.force, args.clean)?;

    Ok(format!("{knots}\n{inv_links}"))
}

fn download(url: &str) -> Result<Vec<u8>, Box<dyn Error>> {
    eprintln!("==> downloading {url}");
    let bytes = ureq::get(url).call()?
        .into_body()
        .into_with_config()
        .limit(MAX_DOWNLOAD)
        .read_to_vec()?;
    Ok(bytes)
}

fn read_local(path: &Path) -> Result<Vec<u8>, Box<dyn Error>> {
    eprintln!("==> reading {}", path.display());
    let bytes = fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    if path.extension().is_some_and(|ext| ext == "zip") {
        unzip_xls(&bytes)
    } else {
        Ok(bytes)
    }
}

fn unzip_xls(zip: &[u8]) -> Result<Vec<u8>, Box<dyn Error>> {
    let mut archive = ZipArchive::new(Cursor::new(zip))?;
    let mut file = archive.by_name(KNOTINFO_XLS)?;
    let mut xls = Vec::with_capacity(file.size() as usize);
    file.read_to_end(&mut xls)?;
    Ok(xls)
}

/// Saves the PD codes to `links/` and the braid words to `braid/`.
fn save_knots(xls: Vec<u8>, out_dir: &Path, force: bool) -> Result<String, Box<dyn Error>> {
    let links_dir = out_dir.join("links");
    let braid_dir = out_dir.join("braid");
    fs::create_dir_all(&links_dir)?;
    fs::create_dir_all(&braid_dir)?;

    eprintln!("==> converting PD codes and braid words");
    let rows = read_rows(xls)?;
    let links = rows.iter().map(|[name, pd, _]|
        save(&links_dir, name, pd, force)
    ).collect::<Result<Vec<_>, _>>()?;
    let braids = rows.iter().map(|[name, _, braid]|
        save(&braid_dir, name, braid, force)
    ).collect::<Result<Vec<_>, _>>()?;

    Ok(format!(
        "{}\n{}",
        summary("PD codes", &links_dir, &links),
        summary("braid words", &braid_dir, &braids),
    ))
}

/// The `[name, pd_notation, braid_notation]` of all rows in the KnotInfo table.
fn read_rows(xls: Vec<u8>) -> Result<Vec<[String; 3]>, Box<dyn Error>> {
    let mut workbook: Xls<_> = open_workbook_from_rs(Cursor::new(xls))?;
    let sheet = workbook.worksheet_range_at(0).ok_or("the workbook has no sheet")??;

    let mut rows = sheet.rows();
    let header = rows.next().ok_or("the sheet is empty")?;
    let col = |key: &str| {
        header.iter().position(|h| h == key).ok_or_else(|| format!("column `{key}` not found"))
    };
    let cols = [col("name")?, col("pd_notation")?, col("braid_notation")?];

    // note: the row after the header holds human-readable column titles
    let rows = rows.skip(1).map(|row|
        cols.map(|i| row[i].to_string().trim().to_string())
    ).collect();

    Ok(rows)
}

fn save(dir: &Path, name: &str, notation: &str, force: bool) -> Result<Outcome, Box<dyn Error>> {
    if name.is_empty() || notation.is_empty() {
        return Ok(Outcome::Skipped);
    }
    if name.contains(['/', '\\']) {
        eprintln!("  skip {name}: invalid name");
        return Ok(Outcome::Skipped);
    }

    let value: Value = match serde_json::from_str(notation) {
        Ok(value) => value,
        Err(e) => {
            eprintln!("  skip {name}: {e}");
            return Ok(Outcome::Skipped);
        }
    };

    let path = dir.join(format!("{name}.json"));
    if path.exists() && !force {
        return Ok(Outcome::Kept);
    }
    fs::write(path, value.to_string())?;
    Ok(Outcome::Written)
}

/// Downloads Lamm's symmetric diagrams and extracts them, flattened, into `inv_link/`.
fn save_inv_links(out_dir: &Path, force: bool, clean: bool) -> Result<String, Box<dyn Error>> {
    let dir = out_dir.join("inv_link");
    if clean && dir.is_dir() {
        eprintln!("==> removing {}", dir.display());
        fs::remove_dir_all(&dir)?;
    }
    fs::create_dir_all(&dir)?;

    let outcomes = INV_LINK_ZIPS.iter().map(|zip| {
        let bytes = download(&format!("{INV_LINK_URL}/{zip}"))?;
        extract_all(&bytes, &dir, force)
    }).collect::<Result<Vec<_>, _>>()?.concat();

    Ok(summary("symmetric diagrams", &dir, &outcomes))
}

fn extract_all(zip: &[u8], dir: &Path, force: bool) -> Result<Vec<Outcome>, Box<dyn Error>> {
    let mut archive = ZipArchive::new(Cursor::new(zip))?;
    (0..archive.len()).map(|i| {
        let mut file = archive.by_index(i)?;
        let path = file.enclosed_name()
            .and_then(|p| p.file_name().map(|name| dir.join(name)))
            .ok_or_else(|| format!("invalid entry `{}`", file.name()))?;

        if path.exists() && !force {
            return Ok(Outcome::Kept);
        }
        let mut json = Vec::with_capacity(file.size() as usize);
        file.read_to_end(&mut json)?;
        fs::write(path, json)?;
        Ok(Outcome::Written)
    }).collect()
}

fn summary(kind: &str, dir: &Path, outcomes: &[Outcome]) -> String {
    let count = |o: Outcome| outcomes.iter().filter(|&&x| x == o).count();
    format!(
        "wrote {} {kind} to {} (kept {} existing, skipped {})",
        count(Outcome::Written),
        dir.display(),
        count(Outcome::Kept),
        count(Outcome::Skipped),
    )
}
