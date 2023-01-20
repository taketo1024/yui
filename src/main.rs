#![allow(dead_code)]

use std::collections::BTreeMap;
use yui::math::traits::{EucRing, EucRingOps};
use yui::links::{Link, links::Edge};
use yui::khovanov::invariants::ss::ss_invariant;

fn main() {
    let target = "14n19265";
    let c: i64 = 2;
    
    let data = load_data().expect("");
    let l = Link::from(data[target].clone());
    run(&l, &target, c);
}

fn run<R>(l: &Link, name: &str, c: R) 
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    let ((h, s), time) = measure(|| {
        ss_invariant(l, c.clone(), true)
    });

    println!("l = {name}, c = {c}");
    println!("h = {h}");
    println!("s = {s}");

    println!("\ntime: {:?}\n", time);
}

fn measure<F, Res>(proc: F) -> (Res, std::time::Duration) 
where F: FnOnce() -> Res { 
    let start = std::time::Instant::now();
    let res = proc();
    let time = start.elapsed();
    (res, time)
}

type Data = BTreeMap<String, Vec<[Edge; 4]>>;
fn load_data() -> Result<Data, Box<dyn std::error::Error>> {
    let path = "resources/targets/targets.json";
    let json = std::fs::read_to_string(path)?;
    let data: Data = serde_json::from_str(&json)?;
    Ok(data)
}

