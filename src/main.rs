#![allow(dead_code)]

use std::collections::BTreeMap;
use yui::math::traits::{EucRing, EucRingOps};
use yui::links::{Link, links::Edge};
use yui::khovanov::invariants::ss::ss_invariant;

fn main() {
    let name = "3_1";
    let l = Link::load(name).unwrap();
    run(&l, name, 2);
}

fn run_all() { 
    let data = load_data().expect("");
    for (name, pd) in data.into_iter() { 
        let l = Link::from(pd);
        for c in [2 as i64, 3] { 
            run(&l, &name, c)
        }
    }
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

