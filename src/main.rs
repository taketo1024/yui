#![allow(dead_code)]

use std::collections::BTreeMap;
use std::panic::UnwindSafe;
use log::{info, error};
use num_bigint::{BigInt, ToBigInt};
use simplelog::*;
use yui::math::ext::int_ext::{Integer, IntOps};
use yui::math::ext::quad_int::{GaussInt, EisenInt};
use yui::math::traits::{EucRing, EucRingOps};
use yui::links::{Link, links::Edge};
use yui::khovanov::invariants::ss::ss_invariant;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    set_logger()?;
    run_all("computation/targets/18.json")
}

fn set_logger() -> Result<(), log::SetLoggerError> {
    TermLogger::init(
        LevelFilter::Trace,
        Config::default(),
        TerminalMode::Mixed,
        ColorChoice::Auto
    )
}

fn run_all(target: &str) -> Result<(), Box<dyn std::error::Error>> {
    let c2: i64 = 2;
    let c3: i64 = 3;
    let cx: GaussInt<i64> = GaussInt::new(1, 1);
    let cy: EisenInt<i64> = EisenInt::new(1, 1);
    
    let file = format!("result.csv");
    let mut csv = csv::Writer::from_path(file)?;

    csv.write_record(vec![
        "name", 
        "s_2", 
        "s_3", 
        "s_{1+i}", 
        "s_{1+ω}"
    ])?;

    let data = load_data(target)?;

    for (name, code) in data { 
        let l = Link::from(&code);
        let s2 = run(&l, &name, &c2);
        let s3 = run(&l, &name, &c3);
        let sx = run(&l, &name, &cx);
        let sy = run(&l, &name, &cy);

        info!("{name}: s2 = {s2}, s3 = {s3}, s_(1+i) = {sx}, s_(1+ω) = {sy}");
        
        csv.write_record(vec![
            name, 
            s2.to_string(), 
            s3.to_string(), 
            sx.to_string(),
            sy.to_string()
        ])?;
        csv.flush()?;
    }

    Ok(())
}

fn run<R>(l: &Link, name: &str, c: &R) -> i32
where 
    R: EucRing, for<'x> &'x R: EucRingOps<R> + UnwindSafe,
    R: ToBig, for<'x> &'x <R as ToBig>::BigR: EucRingOps<<R as ToBig>::BigR>
{ 
    info!("compute {name}, c = {c}");

    let (s, time) = measure(|| {
        std::panic::catch_unwind(|| {
            ss_invariant(l, c.clone(), true)
        }).unwrap_or_else(|_| {
            error!("");
            info!("retry with c = {} ({})", c.to_big(), std::any::type_name::<R::BigR>());
            ss_invariant(l, c.to_big(), true)
        })
    });

    info!("l = {name}, c = {c}");
    info!("s = {s}");
    info!("time: {:?}", time);
    info!("");

    s
}

fn measure<F, Res>(proc: F) -> (Res, std::time::Duration) 
where F: FnOnce() -> Res { 
    let start = std::time::Instant::now();
    let res = proc();
    let time = start.elapsed();
    (res, time)
}

type Data = BTreeMap<String, Vec<[Edge; 4]>>;

fn load_data(target: &str) -> Result<Data, Box<dyn std::error::Error>> {
    let json = std::fs::read_to_string(target)?;
    let data: Data = serde_json::from_str(&json)?;
    Ok(data)
}

trait ToBig
where Self::BigR: EucRing, for<'x> &'x Self::BigR: EucRingOps<Self::BigR> { 
    type BigR;
    fn to_big(&self) -> Self::BigR;
}

impl ToBig for i32 { 
    type BigR = BigInt;
    fn to_big(&self) -> Self::BigR {
        self.to_bigint().unwrap()
    }
}

impl ToBig for i64 { 
    type BigR = BigInt;
    fn to_big(&self) -> Self::BigR {
        self.to_bigint().unwrap()
    }
}

impl<I> ToBig for GaussInt<I>
where I: Integer + ToBigInt, for<'x> &'x I: IntOps<I> {
    type BigR = GaussInt<BigInt>;

    fn to_big(&self) -> Self::BigR {
        let (a, b) = self.pair();
        Self::BigR::new(a.to_bigint().unwrap(), b.to_bigint().unwrap())
    }
}

impl<I> ToBig for EisenInt<I>
where I: Integer + ToBigInt, for<'x> &'x I: IntOps<I> {
    type BigR = EisenInt<BigInt>;

    fn to_big(&self) -> Self::BigR {
        let (a, b) = self.pair();
        Self::BigR::new(a.to_bigint().unwrap(), b.to_bigint().unwrap())
    }
}