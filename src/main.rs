mod app;
mod utils;

fn main() {
    use app::App;
    
    let mut app = App::new();
    let result = app.run();

    if let Ok(output) = result { 
        println!("{output}");
    } else if let Err(code) = result { 
        std::process::exit(code);
    }
}