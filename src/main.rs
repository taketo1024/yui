mod app;
mod utils;
use app::App;

fn main() {
    let app = App::new();
    let res = app.run();
    
    match res { 
        Ok(output) => println!("{output}"),
        Err(code)  => std::process::exit(code)
    }
}