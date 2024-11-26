pub mod app;
use app::App;

fn main() {
    let app = App::new();
    let res = app.run();
    
    match &res { 
        Ok(s)  => println!("{s}"),
        Err(e) => eprintln!("\x1b[0;31merror\x1b[0m: {e}")
    }

    if res.is_err() { 
        std::process::exit(1)
    }
}