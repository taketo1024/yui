use clap::{Parser, Subcommand, ValueEnum};
use log::{info, error};
use simplelog::*;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    #[arg(long, default_value_t = false)]
    debug: bool
}

#[derive(Subcommand, Debug)]
enum Commands {
    SS {
        name: String,
        c_value: String,

        #[arg(short = 't', long, default_value_t = ss::CType::i64)]
        c_type: ss::CType,

        #[arg(short, long)]
        data: Option<String>,

        #[arg(short, long)]
        output: Option<String>,
    },
}

fn main() {
    let args = Cli::parse();

    if args.debug { 
        init_logger();
    }

    let (res, time) = measure(|| 
        guard_panic(||
            match args.command { 
            Commands::SS { name, c_value, c_type, data, output } 
                => ss::run(name, c_value, c_type, data, output)
            }
        )
    );

    if res.is_ok() { 
        info!("time: {:?}", time);
        std::process::exit(0)

    } else if let Err(e) = res { 
        error!("{}", e);
        eprintln!("{e}");
        std::process::exit(1)
    }
}

mod ss { 
    use super::*;
    use yui::links::Link;
    use yui::khovanov::invariants::ss::ss_invariant;
    use derive_more::Display;

    #[allow(non_camel_case_types)]
    #[derive(Clone, ValueEnum, Debug, Display)]
    pub enum CType { 
        i32, i64, i128, bigint, gauss, eisen
    }

    pub fn run(name: String, c_value: String, c_type: CType, data: Option<String>, output: Option<String>) -> Result<(), Box<dyn std::error::Error>> {
        info!("compute ss: {name}, c = {c_value}");

        let l = load_link(&name, &data)?;
        let res = guard_panic(|| 
            compute_ss(&l, &c_value, &c_type)
        );

        if let Ok(s) = res { 
            info!("result: {name}, c = {c_value}: s = {s}");
        }

        if let Some(output) = output { 
            write_res(&name, &c_value, &res, &output)?;
        }

        res.map(|_| ())
    }

    macro_rules! compute_ss_int {
        ($l:expr, $c_value: expr, $t:ty) => {{
            if let Ok(c) = $c_value.parse::<$t>() { 
                Ok( ss_invariant($l, &c, true) )
            } else { 
                err!("invalid c: {} as type: {}.", $c_value, std::any::type_name::<$t>())
            }
        }};
    }

    fn compute_ss(l: &Link, c_value: &String, c_type: &CType) -> Result<i32, Box<dyn std::error::Error>> { 
        use num_bigint::BigInt;
        // use yui::math::types::quad_int::{GaussInt, EisenInt};

        match c_type { 
            CType::i32    => compute_ss_int!(l, c_value, i32),
            CType::i64    => compute_ss_int!(l, c_value, i64),
            CType::i128   => compute_ss_int!(l, c_value, i128),
            CType::bigint => compute_ss_int!(l, c_value, BigInt),
            _ => err!("invalid c-type: {c_type}.")
        }
    }

    fn write_res(name: &String, c_value: &String, res: &Result<i32, Box<dyn std::error::Error>>, output: &String) -> Result<(), Box<dyn std::error::Error>> { 
        use std::fs::OpenOptions;
        use std::path::Path;

        let file_exists = Path::new(output).exists();
        let file = OpenOptions::new()
            .write(true)
            .create(true)
            .append(true)
            .open(output)?;
        let mut wtr = csv::Writer::from_writer(file);

        if !file_exists { 
            wtr.write_record(vec!["name", "c", "ss"])?;
        }

        let res = if let Ok(s) = res { 
            s.to_string()
        } else { 
            "!".to_string()
        };

        wtr.write_record(vec![name, c_value, &res])?;
        wtr.flush()?;

        info!("write: {output}");

        Ok(())
    }
}

fn load_link(name: &String, path: &Option<String>) -> Result<yui::links::Link, Box<dyn std::error::Error>> { 
    use indexmap::IndexMap;
    use yui::links::{Link, links::Edge};
    type Data = IndexMap<String, Vec<[Edge; 4]>>;

    if let Ok(l) = Link::load(&name) { 
        Ok(l)
    } else if let Some(path) = path {
        let json = std::fs::read_to_string(path)?;
        let data: Data = serde_json::from_str(&json)?;
        let Some(code) = data.get(name) else { 
            return err!("{name} not found in file: {path}")
        };
        let l = Link::from(code);
        Ok(l)
    } else { 
        err!("cannot load {name}")
    }
}

fn init_logger() {
    TermLogger::init(
        LevelFilter::Info,
        Config::default(),
        TerminalMode::Mixed,
        ColorChoice::Auto
    ).unwrap()
}

fn measure<F, Res>(proc: F) -> (Res, std::time::Duration) 
where F: FnOnce() -> Res { 
    let start = std::time::Instant::now();
    let res = proc();
    let time = start.elapsed();
    (res, time)
}

fn guard_panic<F, R>(f: F) -> Result<R, Box<dyn std::error::Error>>
where F: FnOnce() -> Result<R, Box<dyn std::error::Error>> + std::panic::UnwindSafe {
    std::panic::catch_unwind(|| {
        f()
    }).unwrap_or_else(|e| {
        let info = match e.downcast::<String>() {
            Ok(v) => *v,
            Err(e) => match e.downcast::<&str>() {
                Ok(v) => v.to_string(),
                _ => "Unknown Source of Error".to_owned()
            }
        };
        err!("panic: {info}")
    })
}

#[derive(Debug, derive_more::Display)]
struct Error(String);
impl std::error::Error for Error {}

macro_rules! err {
    ($($arg:tt)*) => {{
        let msg = format!($($arg)*);
        Err( Error(msg).into() )
    }}
}
pub(crate) use err;