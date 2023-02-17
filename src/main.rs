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

    SSBatch {
        targets: String,
        c_value: String,
        
        #[arg(short = 't', long, default_value_t = ss::CType::i64)]
        c_type: ss::CType,

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
                => ss::run_single(name, c_value, c_type, data, output),
            Commands::SSBatch { targets, c_value, c_type, output }
                => ss::run_batch(targets, c_value, c_type, output)
            }
        )
    );

    if let Ok(res) = res { 
        info!("time: {:?}", time);
        println!("{res}");
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
        i64, i128, bigint, 
        gauss, gauss_128, gauss_big,
        eisen, eisen_128, eisen_big,
    }

    pub fn run_single(name: String, c_value: String, c_type: CType, data: Option<String>, output: Option<String>) -> Result<String, Box<dyn std::error::Error>> {
        let l = load_link(&name, &data)?;
        run(&name, &l, &c_value, &c_type, &output)
    }

    pub fn run_batch(targets: String, c_value: String, c_type: CType, output: Option<String>) -> Result<String, Box<dyn std::error::Error>> {
        let data = load_data(&targets)?;
        let mut all_res = String::from("");

        for (name, code) in data { 
            let l = Link::from(&code);
            let res = run(&name, &l, &c_value, &c_type, &output)?;

            if !all_res.is_empty() { 
                all_res += "\n";
            }
            all_res += &format!("{name}: {res}");
        }
        
        Ok(all_res)
    }

    fn run(name: &String, link: &Link, c_value: &String, c_type: &CType, output: &Option<String>) -> Result<String, Box<dyn std::error::Error>> {
        info!("compute ss: {name}, c = {c_value}");

        let res = guard_panic(|| 
            compute_ss(name, link, c_value, c_type)
        );

        if let Ok(s) = res { 
            info!("{name}: s = {s} (c = {c_value})")
        }

        if let Some(output) = output { 
            write_res(&name, &c_value, &res, &output)?;
        }

        res.map(|s| s.to_string() )
    }

    macro_rules! ss_parse {
        ($l:expr, $c_value: expr, $t:ty) => {{
            if let Ok(c) = $c_value.parse::<$t>() { 
                Ok( ss_invariant($l, &c, true) )
            } else { 
                err!("invalid c: {} as type: {}.", $c_value, std::any::type_name::<$t>())
            }
        }};
    }

    fn compute_ss(name: &String, l: &Link, c_value: &String, c_type: &CType) -> Result<i32, Box<dyn std::error::Error>> { 
        use num_bigint::BigInt;
        use yui::math::types::quad_int::{GaussInt, EisenInt};

        info!("compute ss: {name}, c = {c_value}");

        match c_type { 
            CType::i64        => ss_parse!(l, c_value, i64),
            CType::i128       => ss_parse!(l, c_value, i128),
            CType::bigint     => ss_parse!(l, c_value, BigInt),
            CType::gauss      => ss_parse!(l, c_value, GaussInt<i64>), 
            CType::gauss_128  => ss_parse!(l, c_value, GaussInt<i128>), 
            CType::gauss_big  => ss_parse!(l, c_value, GaussInt<BigInt>), 
            CType::eisen      => ss_parse!(l, c_value, EisenInt<i64>), 
            CType::eisen_128  => ss_parse!(l, c_value, EisenInt<i128>), 
            CType::eisen_big  => ss_parse!(l, c_value, EisenInt<BigInt>), 
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

use indexmap::IndexMap;
use yui::links::{Link, links::Edge};

type Data = IndexMap<String, Vec<[Edge; 4]>>;
fn load_data(path: &String) -> Result<Data, Box<dyn std::error::Error>> { 
    let json = std::fs::read_to_string(path)?;
    let data: Data = serde_json::from_str(&json)?;
    Ok(data)
}

fn load_link(name: &String, path: &Option<String>) -> Result<Link, Box<dyn std::error::Error>> { 
    if let Ok(l) = Link::load(&name) { 
        Ok(l)
    } else if let Some(path) = path {
        let data = load_data(path)?;
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