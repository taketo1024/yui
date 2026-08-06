//! The `ykh` command-line entry point; see the crate docs for the commands.

use ykh::App;

// Opt-in (`--features mimalloc`): mimalloc showed 10–20× slowdowns and one unreproduced
// crash on the Mac Studio's cobordism-elimination workload — see wh-p5-h1-success-2026-07.md.
#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

fn main() {
    let app = App::new();
    let res = app.run();

    match res {
        Ok(output) => println!("{output}"),
        Err(e) => {
            log::error!("{}", e);
            eprintln!("\x1b[0;31merror\x1b[0m: {e}");
            std::process::exit(1)
        }
    }
}