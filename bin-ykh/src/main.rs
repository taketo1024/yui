use ykh::App;

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