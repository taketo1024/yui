//! Profile `khi` on a 32-crossing strongly invertible knot via the CLI `App`,
//! with the PD baked in — pass any `khi` flags after `--`.
//!
//! ```
//! cargo run -r --example profile_khi_32 -- --node min-cut --mode min-fill --chunk 20 --log 2
//! ```

use clap::Parser;
use ykh::{App, CliArgs};

#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

const PD: &str = "[[48,2,49,1],[2,48,3,47],[40,3,41,4],[39,47,40,46],[9,5,10,4],[10,45,11,46],[8,23,9,24],[41,25,42,24],[7,59,8,58],[42,57,43,58],[59,7,60,6],[60,43,61,44],[22,5,23,6],[21,45,22,44],[61,57,62,56],[62,25,63,26],[19,27,20,26],[20,55,21,56],[18,64,19,63],[64,18,1,17],[11,55,12,54],[12,27,13,28],[38,53,39,54],[37,29,38,28],[29,37,30,36],[30,13,31,14],[52,35,53,36],[51,15,52,14],[15,51,16,50],[16,31,17,32],[33,33,34,32],[34,49,35,50]]";

fn main() {
    let mut argv = vec!["ykh".to_string(), "khi".to_string(), PD.to_string()];
    argv.extend(std::env::args().skip(1)); // forward flags passed after `--`

    let app = App { args: CliArgs::parse_from(argv) };
    match app.run() {
        Ok(out) => println!("{out}"),
        Err(e)  => { eprintln!("error: {e}"); std::process::exit(1) }
    }
}
