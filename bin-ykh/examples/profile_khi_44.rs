//! Profile `khi` on a 44-crossing strongly invertible knot via the CLI `App`,
//! with the PD baked in — pass any `khi` flags after `--`.
//!
//! ```
//! cargo run -r --example profile_khi_44 -- -c H --node min-cut --mode min-fill --chunk 28 --h-range "..1" --log 2
//! ```

use clap::Parser;
use ykh::{App, CliArgs};

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

const PD: &str = "[[2,64,3,63],[4,62,5,61],[7,31,8,30],[8,79,9,80],[11,19,12,18],[12,47,13,48],[13,75,14,74],[14,35,15,36],[19,7,20,6],[20,59,21,60],[24,86,25,85],[26,84,27,83],[27,43,28,42],[28,67,29,68],[31,11,32,10],[32,55,33,56],[37,51,38,50],[38,15,39,16],[39,35,40,34],[40,75,41,76],[43,23,44,22],[45,61,46,60],[46,5,47,6],[51,37,52,36],[52,73,53,74],[53,49,54,48],[54,17,55,18],[57,81,58,80],[58,29,59,30],[62,4,63,3],[64,2,65,1],[65,45,66,44],[66,21,67,22],[69,77,70,76],[70,33,71,34],[71,17,72,16],[72,49,73,50],[77,57,78,56],[78,9,79,10],[81,69,82,68],[82,41,83,42],[84,26,85,25],[86,24,1,23]]";

fn main() {
    let mut argv = vec!["ykh".to_string(), "khi".to_string(), PD.to_string()];
    argv.extend(std::env::args().skip(1)); // forward flags passed after `--`

    let app = App { args: CliArgs::parse_from(argv) };
    match app.run() {
        Ok(out) => println!("{out}"),
        Err(e)  => { eprintln!("error: {e}"); std::process::exit(1) }
    }
}
