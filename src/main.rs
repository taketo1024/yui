use std::collections::BTreeMap;
use yui::links::{Link, links::Edge};
use yui::khovanov::homology::KhHomology;

fn main() {
    let data = load_data().expect("");
    let pd = data["14n19265"].clone();
    let l = Link::from(pd);

    measure(|| {
        let h = KhHomology::<i32>::new(l);
        println!("{h}");
    });
}

type Data = BTreeMap<String, Vec<[Edge; 4]>>;
fn load_data() -> Result<Data, Box<dyn std::error::Error>> {
    let path = "resources/targets/targets.json";
    let json = std::fs::read_to_string(path)?;
    let data: Data = serde_json::from_str(&json)?;
    Ok(data)
}

fn measure<F>(proc: F) where F: FnOnce() -> () { 
    let start = std::time::Instant::now();
    proc();
    let dur = start.elapsed();
    println!("time: {:?}", dur);
}