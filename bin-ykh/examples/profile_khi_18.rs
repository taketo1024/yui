//! Profile `khi` on an 18-crossing strongly invertible knot via the CLI `App`,
//! with the PD baked in — pass any `khi` flags after `--`.
//!
//! ```
//! cargo run -r --example profile_khi_18 -- --no-preprocess --log 2
//! ```

use clap::Parser;
use ykh::{App, CliArgs};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

const PD: &str = "[[1,27,2,26],[5,16,6,17],[6,32,7,31],[10,27,11,28],[11,1,12,36],[13,8,14,9],[14,20,15,19],[17,4,18,5],[18,24,19,23],[21,32,22,33],[22,16,23,15],[25,3,26,2],[28,9,29,10],[29,24,30,25],[30,4,31,3],[33,20,34,21],[34,8,35,7],[35,13,36,12]]";

fn main() {
    let mut argv = vec!["ykh".to_string(), "khi".to_string(), PD.to_string()];
    argv.extend(std::env::args().skip(1)); // forward flags passed after `--`

    let app = App { args: CliArgs::parse_from(argv) };
    match app.run() {
        Ok(out) => println!("{out}"),
        Err(e)  => { eprintln!("error: {e}"); std::process::exit(1) }
    }
}
