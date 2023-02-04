#![allow(unused_imports)]

use std::panic::UnwindSafe;
use derive_more::{Display, Error};
use log::{info, error};
use simplelog::*;
use num_bigint::BigInt;
use indexmap::IndexMap;
use yui::math::types::quad_int::{GaussInt, EisenInt};
use yui::math::traits::{EucRing, EucRingOps};
use yui::links::{Link, links::Edge};
use yui::khovanov::invariants::ss::ss_invariant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    set_logger()?;
    run_all("computation/targets/18.json")
}

fn set_logger() -> Result<(), log::SetLoggerError> {
    TermLogger::init(
        LevelFilter::Info,
        Config::default(),
        TerminalMode::Mixed,
        ColorChoice::Auto
    )
}

fn run_all(target: &str) -> Result<(), Box<dyn std::error::Error>> {
    let c2: BigInt = 2.into();
    let c3: BigInt = 3.into();
    let cx: GaussInt<BigInt> = GaussInt::new(1.into(), 1.into());
    let cy: EisenInt<BigInt> = EisenInt::new(1.into(), 1.into());
    
    let file = format!("result.csv");
    let mut csv = csv::Writer::from_path(file)?;

    csv.write_record(vec![
        "name", 
        "(Z, 2)", 
        "(Z, 3)", 
        "(Z[i], 1+i)", 
        "(Z[j], 1+j)"
    ])?;
    csv.flush()?;

    let data = load_data(target)?;

    for (name, code) in data { 
        let l = Link::from(&code);
        let list = run_for!(&l, &name, &c2, &c3, &cx, &cy);
        let record = [vec![name], list].concat();

        csv.write_record(record)?;
        csv.flush()?;
    }

    Ok(())
}

macro_rules! run_for {
    ($l:expr, $name:expr) => { vec![] };
    ($l:expr, $name:expr, $c:expr $(,$next:expr)*) => {{
        let res = run($l, $name, $c);
        let str = if let Ok(res) = res { 
            res.to_string()
        } else { 
            "!".to_string()
        };
        let mut list = vec![str];
        let mut other = run_for!($l, $name $(,$next)*);
        list.append(&mut other);

        list
    }};
}

use run_for;

fn run<R>(l: &Link, name: &str, c: &R) -> Result<i32, Box<dyn std::error::Error>> 
where R: EucRing, for<'x> &'x R: EucRingOps<R> + UnwindSafe { 
    info!("start: {name}");

    let (res, time) = measure(|| {
        std::panic::catch_unwind(|| {
            let s = ss_invariant(l, c, true);
            Ok(s)
        }).unwrap_or_else(|_| {
            error!("panic!");
            Err(Panic().into())
        })
    });

    if let Ok(s) = res {
        info!("{name}: ss = {s}");
        info!("time: {:?}", time);
        info!("");
    }

    res
}

#[derive(Debug, Display, Error)]
struct Panic();

fn measure<F, Res>(proc: F) -> (Res, std::time::Duration) 
where F: FnOnce() -> Res { 
    let start = std::time::Instant::now();
    let res = proc();
    let time = start.elapsed();
    (res, time)
}

type Data = IndexMap<String, Vec<[Edge; 4]>>;

fn load_data(target: &str) -> Result<Data, Box<dyn std::error::Error>> {
    let json = std::fs::read_to_string(target)?;
    let data: Data = serde_json::from_str(&json)?;
    Ok(data)
}