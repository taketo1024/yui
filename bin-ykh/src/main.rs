//! The `ykh` command-line entry point; see the crate docs for the commands.

use ykh::App;

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