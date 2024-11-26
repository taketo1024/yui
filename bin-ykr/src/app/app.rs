use itertools::Itertools;
use log::{info, error};
use clap::{Parser, ValueEnum};
use derive_more::Display;
use num_bigint::BigInt;
use yui::{Ratio, EucRing, EucRingOps};
use yui::tex::TeX;
use yui_kr::KRHomologyStr;
use yui_link::{Braid, Link};
use super::calc::{KRCalc, KRCalcMode};
use super::utils::*;

#[derive(Parser, Debug, Default)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct CliArgs {
    pub target: String,

    #[arg(short, long)]
    pub braid: Option<String>,

    #[arg(short, long, default_value = "i64")]
    pub int_type: IntType,

    #[arg(short, long)]
    pub mirror: bool,

    #[arg(short, long, default_value = "poincare")]
    pub format: Output,

    #[arg(short = 'F', long)]
    pub force_compute: bool,

    #[arg(short, long, default_value = "default")]
    pub compute_mode: KRCalcMode,

    #[arg(short, long)]
    pub limit: Option<usize>,

    #[arg(short = 'p', long)]
    pub save_progress: bool,

    #[arg(long)]
    pub debug: bool
}

#[derive(ValueEnum, Clone, Copy, Default, Debug)]
#[clap(rename_all="lower")]
pub enum IntType { 
    #[default]
    I64, 
    I128, 
    BigInt
}

#[derive(ValueEnum, Clone, Copy, Default, Debug)]
pub enum Output { 
    #[default]
    Poincare,
    PoincareTex,
    Delta,
    Table, 
}

#[derive(Debug, Display)]
pub struct AppErr(pub(crate) String);
impl std::error::Error for AppErr {}

macro_rules! err {
    ($($arg:tt)*) => {{
        let msg = format!($($arg)*);
        let e = AppErr(msg);
        Result::Err( Box::new(e) )
    }}
}

pub(crate) use err;

pub struct App {
    args: CliArgs
}

impl App { 
    pub fn new() -> Self { 
        let args = CliArgs::parse();
        Self::new_with(args)
    }

    pub fn new_with(args: CliArgs) -> Self { 
        if args.debug { 
            std::env::set_var("RUST_BACKTRACE", "1");
            Self::init_logger();
        }
        App { args }
    }

    fn init_logger() {
        let l = log::LevelFilter::Info;
        yui::util::log::init_simple_logger(l).unwrap()
    }

    pub fn run(&self) -> Result<String, Box<dyn std::error::Error>> { 
        info!("args: {:?}", self.args);

        let n_threads = std::thread::available_parallelism().map(|x| x.get()).unwrap_or(1);
        info!("n-threads: {n_threads}");

        let (res, time) = measure(|| guard_panic(||
            self.compute()
        ));

        if let Some(e) = res.as_ref().err() { 
            error!("{e}");
        }

        let res = res?;
        let output = self.format(&res);

        info!("time: {:?}", time);

        Ok(output)
    }

    pub fn compute(&self) -> Result<KRHomologyStr, Box<dyn std::error::Error>> { 
        let target = &self.args.target;
        let file = File::Result(target);

        if !self.args.force_compute && self.args.braid.is_none() && file.exists() { 
            info!("result exists: {}", target);
            let res: KRHomologyStr = file.read()?;
            let res = if self.args.mirror { res.mirror() } else { res };
            return Ok(res)
        }
        
        let link = self.make_link()?;

        let res = if link.writhe() > 0 { 
            info!("compute from mirror.");
            self.compute_kr_dispatch(&link.mirror())?.mirror()
        } else { 
            self.compute_kr_dispatch(&link)?
        };

        Ok(res)
    }

    fn compute_kr_dispatch(&self, link: &Link) -> Result<KRHomologyStr, Box<dyn std::error::Error>> { 
        match self.args.int_type { 
            IntType::I64    => self.compute_kr::<Ratio<i64>>(link),
            IntType::I128   => self.compute_kr::<Ratio<i128>>(link),
            IntType::BigInt => self.compute_kr::<Ratio<BigInt>>(link),
        }
    }

    fn compute_kr<R>(&self, link: &Link) -> Result<KRHomologyStr, Box<dyn std::error::Error>>
    where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
        let mut calc = KRCalc::init(&self.args.target, link);
        
        calc.mode = self.args.compute_mode;

        if let Some(limit) = self.args.limit { 
            calc.size_limit = limit;
        }

        if self.args.save_progress { 
            if self.args.force_compute { 
                calc.clear()?;
            } else { 
                calc.load_if_exists()?;
            }
            calc.save_progress = true;
        }

        calc.compute()?;

        Ok(calc.into_result())
    }

    fn format(&self, res: &KRHomologyStr) -> String { 
        let str = match &self.args.format {
            Output::Poincare    => res.poincare_poly().to_string(),
            Output::PoincareTex => res.poincare_poly().tex_string(),
            Output::Delta => res.delta_table(),
            Output::Table => res.qpoly_table(),
        };

        if res.is_determined() { 
            str
        } else { 
            let non_det = res.non_determined().map(|idx| format!("{idx:?}")).join(", ");
            str + &format!("\n\x1b[0;31mNon-determined\x1b[0m: {non_det}")
        }
    }

    fn make_link(&self) -> Result<Link, Box<dyn std::error::Error>> { 
        let braid = if let Some(b) = &self.args.braid { 
            let seq: Vec<i32> = serde_json::from_str(&b)?;
            Braid::from_iter(seq)
        } else { 
            let target = &self.args.target;
            Braid::load(target)?
        };

        info!("braid: {}", braid);
        info!("braid-len: {}", braid.len());

        let link = braid.closure();
        let link = if self.args.mirror { link.mirror() } else { link };

        Ok(link)
    }
}