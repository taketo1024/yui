use clap::{Parser, Subcommand, ValueEnum};
use log::{info, error};
use simplelog::*;

#[allow(non_camel_case_types)]
#[derive(Clone, ValueEnum, Debug, derive_more::Display)]
pub enum CType { 
    Z, Z_128, Z_big, 
    Q, Q_128, Q_big, 
    F2, F3, F5, F7,
    PolyQ, PolyQ_128, PolyQ_big, 
    PolyF2, PolyF3, PolyF5, PolyF7,
    Gauss, Gauss_128, Gauss_big,
    Eisen, Eisen_128, Eisen_big,
    None
}

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
    Kh {
        name: String,

        #[arg(short, long)]
        link: Option<String>,

        #[arg(short, long, default_value = "0")]
        c_value: String,

        #[arg(short = 't', long, default_value = "z")]
        c_type: CType,

        #[arg(short, long)]
        mirror: bool,

        #[arg(short, long)]
        reduced: bool,

        #[arg(short, long)]
        bigraded: bool,
    },

    SS {
        name: String,

        #[arg(short, long)]
        link: Option<String>,

        #[arg(short, long)]
        c_value: String,

        #[arg(short = 't', long, default_value = "z")]
        c_type: CType,

        #[arg(short, long)]
        output: Option<String>,
    },

    SSBatch {
        targets: String,
        c_value: String,
        
        #[arg(short = 't', long, default_value = "z")]
        c_type: CType,

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
            Commands::Kh { name, link, c_value, c_type, mirror, reduced, bigraded }
                => kh::run(name, link, c_value, c_type, mirror, reduced, bigraded),
            Commands::SS { name, link, c_value, c_type, output } 
                => ss::run_single(name, link, c_value, c_type, output),
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

mod kh { 
    use std::str::FromStr;
    use yui::khovanov::homology::{KhHomologyBigraded, KhHomology};
    use yui::links::Link;
    use yui::math::homology::base::PrintTable;
    use yui::math::traits::{EucRing, EucRingOps};
    use super::*;

    pub fn run(name: String, link: Option<String>, c_value: String, c_type: CType, mirror: bool, reduced: bool, bigraded: bool) -> Result<String, Box<dyn std::error::Error>> {
        let mut l = load_link(&name, &link)?;
        if mirror { 
            l = l.mirror();
        }

        if bigraded { 
            dispatch!(compute_bigraded, &c_type, l, &c_value, reduced)
        } else { 
            dispatch!(compute_homology, &c_type, l, &c_value, reduced)
        }
    }

    fn compute_bigraded<R>(l: Link, c_value: &String, reduced: bool) -> Result<String, Box<dyn std::error::Error>>
    where R: EucRing + FromStr, for<'x> &'x R: EucRingOps<R> { 
        if let Ok(c) = c_value.parse::<R>() { 
            if !c.is_zero() { 
                return err!("--bigraded only supported for `c = 0`.")
            }
            let h = KhHomologyBigraded::new(l, reduced);
            let table = h.table();
            Ok(table)
        } else { 
            err!("cannot parse c: '{}' as type: {}.", c_value, std::any::type_name::<R>())
        }
    }

    fn compute_homology<R>(l: Link, c_value: &String, reduced: bool) -> Result<String, Box<dyn std::error::Error>>
    where R: EucRing + FromStr, for<'x> &'x R: EucRingOps<R> { 
        if let Ok(c) = c_value.parse::<R>() { 
            let h = KhHomology::new(l, c, R::zero(), reduced);
            let res = h.to_string();
            Ok(res)
        } else { 
            err!("cannot parse c: '{}' as type: {}.", c_value, std::any::type_name::<R>())
        }
    }
}

mod ss { 
    use std::str::FromStr;
    use yui::links::Link;
    use yui::math::traits::{EucRing, EucRingOps};
    use yui::khovanov::invariants::ss::ss_invariant;
    use super::*;

    pub fn run_single(name: String, link: Option<String>, c_value: String, c_type: CType, output: Option<String>) -> Result<String, Box<dyn std::error::Error>> {
        let l = load_link(&name, &link)?;
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
            dispatch!(compute_ss, c_type, link, c_value)
        );

        if let Ok(s) = res { 
            info!("{name}: s = {s} (c = {c_value})")
        }

        if let Some(output) = output { 
            write_res(&name, &c_value, &res, &output)?;
        }

        res.map(|s| s.to_string() )
    }

    fn compute_ss<R>(l: &Link, c_value: &String) -> Result<i32, Box<dyn std::error::Error>>
    where R: EucRing + FromStr, for<'x> &'x R: EucRingOps<R> { 
        if let Ok(c) = c_value.parse::<R>() { 
            Ok( ss_invariant(l, &c, true) )
        } else { 
            err!("cannot parse c: '{}' as type: {}.", c_value, std::any::type_name::<R>())
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

type PDCode = Vec<[Edge; 4]>;
type Data = IndexMap<String, PDCode>;

fn load_data(path: &String) -> Result<Data, Box<dyn std::error::Error>> { 
    let json = std::fs::read_to_string(path)?;
    let data: Data = serde_json::from_str(&json)?;
    Ok(data)
}

fn load_link(name: &String, pd_code: &Option<String>) -> Result<Link, Box<dyn std::error::Error>> { 
    if let Some(pd_code) = pd_code { 
        let pd_code: PDCode = serde_json::from_str(&pd_code)?;
        let l = Link::from(&pd_code);
        Ok(l)
    } else if let Ok(l) = Link::load(&name) { 
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

use num_bigint::BigInt;
use yui::math::types::quad_int::{GaussInt, EisenInt};
use yui::math::types::{ratio::Ratio, fin_field::FF};
use yui::math::types::polynomial::Poly;

macro_rules! dispatch {
    ($method:ident, $c_type:expr $(, $args:expr)*) => {
        match $c_type { 
            CType::Z          => $method::<i64>($($args),*),
            CType::Z_128      => $method::<i128>($($args),*),
            CType::Z_big      => $method::<BigInt>($($args),*),
            CType::Q          => $method::<Ratio<i64>>($($args),*),
            CType::Q_128      => $method::<Ratio<i128>>($($args),*),
            CType::Q_big      => $method::<Ratio<BigInt>>($($args),*),
            CType::F2         => $method::<FF<2>>($($args),*),
            CType::F3         => $method::<FF<3>>($($args),*),
            CType::F5         => $method::<FF<5>>($($args),*),
            CType::F7         => $method::<FF<7>>($($args),*),
            CType::PolyQ      => $method::<Poly<'H', Ratio<i64>>>($($args),*),
            CType::PolyQ_128  => $method::<Poly<'H', Ratio<i128>>>($($args),*),
            CType::PolyQ_big  => $method::<Poly<'H', Ratio<BigInt>>>($($args),*),
            CType::PolyF2     => $method::<Poly<'H', FF<2>>>($($args),*),
            CType::PolyF3     => $method::<Poly<'H', FF<3>>>($($args),*),
            CType::PolyF5     => $method::<Poly<'H', FF<5>>>($($args),*),
            CType::PolyF7     => $method::<Poly<'H', FF<7>>>($($args),*),
            CType::Gauss      => $method::<GaussInt<i64>>($($args),*), 
            CType::Gauss_128  => $method::<GaussInt<i128>>($($args),*), 
            CType::Gauss_big  => $method::<GaussInt<BigInt>>($($args),*), 
            CType::Eisen      => $method::<EisenInt<i64>>($($args),*), 
            CType::Eisen_128  => $method::<EisenInt<i128>>($($args),*), 
            CType::Eisen_big  => $method::<EisenInt<BigInt>>($($args),*), 
            _ => err!("cannot dispatch {} for c-type: {}", stringify!($method), $c_type)
        }
    };
}

pub(crate) use dispatch;

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