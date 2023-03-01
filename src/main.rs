use clap::{Parser, Subcommand, ValueEnum};
use log::{info, error};
use simplelog::*;

#[allow(non_camel_case_types)]
#[derive(Clone, ValueEnum, Debug, derive_more::Display)]
pub enum CType { 
    Z, Z_128, Z_big, 
    Q, Q_128, Q_big, 
    F2, F3, F5, F7,
    ZPoly, ZPoly_128, ZPoly_big, 
    QPoly, QPoly_128, QPoly_big, 
    F2Poly, F3Poly, F5Poly, F7Poly,
    ZMpoly, ZMpoly_128, ZMpoly_big,
    QMpoly, QMpoly_128, QMpoly_big, 
    F2Mpoly, F3Mpoly, F5Mpoly, F7Mpoly,
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
    
    Ckh {
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
            Commands::Ckh { name, link, c_value, c_type, mirror, reduced }
                => ckh::run(name, link, c_value, c_type, mirror, reduced),
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
    use super::*;
    use std::str::FromStr;
    use yui::khovanov::homology::{KhHomologyBigraded, KhHomology};
    use yui::links::Link;
    use yui::math::homology::base::PrintTable;
    use yui_core::{EucRing, EucRingOps};

    pub fn run(name: String, link: Option<String>, c_value: String, c_type: CType, mirror: bool, reduced: bool, bigraded: bool) -> Result<String, Box<dyn std::error::Error>> {
        let mut l = load_link(&name, &link)?;
        if mirror { 
            l = l.mirror();
        }

        if bigraded { 
            dispatch_eucring!(compute_bigraded, &c_type, &l, &c_value, reduced)
        } else { 
            dispatch_eucring!(compute_homology, &c_type, &l, &c_value, reduced)
        }
    }

    fn compute_bigraded<R>(l: &Link, c_value: &String, reduced: bool) -> Result<String, Box<dyn std::error::Error>>
    where R: EucRing + FromStr, for<'x> &'x R: EucRingOps<R> { 
        if c_value.as_str() != "0" { 
            return err!("--bigraded only supported for `c = 0`.")
        }
        let h = KhHomologyBigraded::new(l.clone(), reduced);
        let table = h.table();
        Ok(table)
    }

    fn compute_homology<R>(l: &Link, c_value: &String, reduced: bool) -> Result<String, Box<dyn std::error::Error>>
    where R: EucRing + FromStr, for<'x> &'x R: EucRingOps<R> { 
        let (h, t) = parse_pair::<R>(c_value)?;
        if reduced && !t.is_zero() { 
            return err!("{t} != 0 is not allowed for reduced.");
        }

        let kh = KhHomology::new(l.clone(), h, t, reduced);
        let res = kh.to_string();
        Ok(res)
    }
}

mod ckh { 
    use super::*;
    use std::str::FromStr;
    use yui::khovanov::complex::KhComplex;
    use yui::links::Link;
    use yui::math::homology::reduced::Reduced;
    use yui::math::homology::complex::*;
    use yui::math::homology::base::*;
    use yui_core::{Ring, RingOps};
    use yui::math::matrix::dense::*;
    
    pub fn run(name: String, link: Option<String>, c_value: String, c_type: CType, mirror: bool, reduced: bool) -> Result<String, Box<dyn std::error::Error>> {
        let mut l = load_link(&name, &link)?;
        if mirror { 
            l = l.mirror();
        }

        dispatch_ring!(describe_ckh, &c_type, &l, &c_value, reduced)
    }

    fn describe_ckh<R>(l: &Link, c_value: &String, reduced: bool) -> Result<String, Box<dyn std::error::Error>>
    where R: Ring + FromStr, for<'x> &'x R: RingOps<R> { 
        use string_builder::Builder;

        let (h, t) = parse_pair::<R>(c_value)?;
        if reduced && !t.is_zero() { 
            return err!("{t} != 0 is not allowed for reduced.");
        }
        
        let ckh = KhComplex::new(l.clone(), h, t, reduced);
        let ckh = Reduced::from(ckh);

        let mut b = Builder::new(1024);
        for i in ckh.range() {
            b.append(format!("C[{}]: {} -> {}\n", i, ckh[i], ckh[i+1]));

            let d = ckh.d_matrix(i);
            if d.rows() > 0 && d.cols() > 0 {
                b.append(d.to_dense().to_string());
                b.append("\n");
            }
            b.append("\n");
        }

        let res = b.string()?;
        Ok(res)
    }
}

mod ss { 
    use super::*;
    use std::str::FromStr;
    use yui::links::Link;
    use yui_core::{EucRing, EucRingOps};
    use yui::khovanov::invariants::ss::ss_invariant;

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
            dispatch_eucring!(compute_ss, c_type, link, c_value)
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
        let Ok(c) = R::from_str(c_value) else { 
            return err!("cannot parse c: '{}' as type: {}.", c_value, std::any::type_name::<R>())
        };

        let ss = ss_invariant(l, &c, true);
        Ok(ss)
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

use std::str::FromStr;
use num_traits::Zero;

fn parse_pair<R: FromStr + Zero>(s: &String) -> Result<(R, R), Box<dyn std::error::Error>> { 
    if let Ok(c) = R::from_str(s) { 
        return Ok((c, R::zero()))
    }

    let r = regex::Regex::new(r"^(.+),(.+)$").unwrap();
    if let Some(m) = r.captures(&s) { 
        let (s1, s2) = (&m[1], &m[2]);
        if let (Ok(a), Ok(b)) = (R::from_str(s1), R::from_str(s2)) {
            return Ok((a, b))
        }
    }

    err!("cannot parse '{}' as {}.", s, std::any::type_name::<R>())
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
use yui_ratio::Ratio;
use yui_ff::FF;
use yui_quad_int::{GaussInt, EisenInt};
use yui_polynomial::{Poly, MPoly};

macro_rules! dispatch_ring {
    ($method:ident, $c_type:expr $(, $args:expr)*) => {{
        match $c_type { 
            CType::ZPoly      => $method::<Poly<'H', i64>>($($args),*),
            CType::ZPoly_128  => $method::<Poly<'H', i128>>($($args),*),
            CType::ZPoly_big  => $method::<Poly<'H', BigInt>>($($args),*),
            CType::ZMpoly     => $method::<MPoly<'H', i64>>($($args),*),
            CType::ZMpoly_128 => $method::<MPoly<'H', i128>>($($args),*),
            CType::ZMpoly_big => $method::<MPoly<'H', BigInt>>($($args),*),
            CType::QMpoly     => $method::<MPoly<'H', Ratio<i64>>>($($args),*),
            CType::QMpoly_128 => $method::<MPoly<'H', Ratio<i128>>>($($args),*),
            CType::QMpoly_big => $method::<MPoly<'H', Ratio<BigInt>>>($($args),*),
            CType::F2Mpoly    => $method::<MPoly<'H', FF<2>>>($($args),*),
            CType::F3Mpoly    => $method::<MPoly<'H', FF<3>>>($($args),*),
            CType::F5Mpoly    => $method::<MPoly<'H', FF<5>>>($($args),*),
            CType::F7Mpoly    => $method::<MPoly<'H', FF<7>>>($($args),*),
            _ => dispatch_eucring!($method, $c_type, $($args),*)
        }
    }};
}

macro_rules! dispatch_eucring {
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
            CType::QPoly      => $method::<Poly<'H', Ratio<i64>>>($($args),*),
            CType::QPoly_128  => $method::<Poly<'H', Ratio<i128>>>($($args),*),
            CType::QPoly_big  => $method::<Poly<'H', Ratio<BigInt>>>($($args),*),
            CType::F2Poly     => $method::<Poly<'H', FF<2>>>($($args),*),
            CType::F3Poly     => $method::<Poly<'H', FF<3>>>($($args),*),
            CType::F5Poly     => $method::<Poly<'H', FF<5>>>($($args),*),
            CType::F7Poly     => $method::<Poly<'H', FF<7>>>($($args),*),
            CType::Gauss      => $method::<GaussInt<i64>>($($args),*), 
            CType::Gauss_128  => $method::<GaussInt<i128>>($($args),*), 
            CType::Gauss_big  => $method::<GaussInt<BigInt>>($($args),*), 
            CType::Eisen      => $method::<EisenInt<i64>>($($args),*), 
            CType::Eisen_128  => $method::<EisenInt<i128>>($($args),*), 
            CType::Eisen_big  => $method::<EisenInt<BigInt>>($($args),*), 
            _ => err!("{} is not supported for: {}", stringify!($method), $c_type)
        }
    };
}

pub(crate) use {dispatch_ring, dispatch_eucring};

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