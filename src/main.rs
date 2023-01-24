#![allow(dead_code)]

use std::collections::BTreeMap;
use yui::math::ext::quad_int::{GaussInt, EisenInt};
use yui::math::traits::{EucRing, EucRingOps};
use yui::links::{Link, links::Edge};
use yui::khovanov::invariants::ss::ss_invariant;

fn main() {
    run_all().expect("");
}

fn run_all() -> Result<(), Box<dyn std::error::Error>> {
    let c2: i64 = 2;
    let c3: i64 = 3;
    let cx: GaussInt<i64> = GaussInt::new(1, 1);
    let cy: EisenInt<i64> = EisenInt::new(1, 1);
    
    let mut csv = csv::Writer::from_path("result.csv")?;
    csv.write_record(vec![
        "name", 
        "s_2", 
        "s_3", 
        "s_{1+i}", 
        "s_{1+ω}"
    ])?;

    let data = load_data()?;

    for (target, code) in data { 
        let l = Link::from(&code);
        let s2 = run(&l, &target, c2);
        let s3 = run(&l, &target, c3);
        let sx = run(&l, &target, cx.clone());
        let sy = run(&l, &target, cy.clone());

        println!("{target}: s2 = {s2}, s3 = {s3}, s_(1+i) = {sx}, s_(1+ω) = {sy}\n");
        
        csv.write_record(vec![
            target, 
            s2.to_string(), 
            s3.to_string(), 
            sx.to_string(),
            sy.to_string()
        ])?;
        csv.flush()?;
    }

    Ok(())
}

fn run<R>(l: &Link, name: &str, c: R) -> i32
where R: EucRing, for<'x> &'x R: EucRingOps<R> { 
    let ((h, s), time) = measure(|| {
        ss_invariant(l, c.clone(), true)
    });

    println!("l = {name}, c = {c}");
    println!("h = {h}");
    println!("s = {s}");

    println!("\ntime: {:?}\n", time);

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
fn load_data() -> Result<Data, Box<dyn std::error::Error>> {
    let path = "resources/targets/targets.json";
    let json = std::fs::read_to_string(path)?;
    let data: Data = serde_json::from_str(&json)?;
    Ok(data)
}

