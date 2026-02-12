use std::collections::HashMap;

use itertools::Itertools;
use yui_core::num::Ratio;
use yui_homology::{DisplaySeq, DisplayTable, SummandTrait};
use yui_kh::kh::{KhComplex, KhHomology};
use yui_link::Link;

type Dict = HashMap<String, Vec<String>>;
fn main() {
    let n = 11;
    let dict = make_dict(n);
    let targets = dupl_list(dict);
    let res = targets.iter().map(|list| { 
        println!("distinguish: {list:?}");
        distinguish(&list)
    }).collect_vec();

    for (i, list) in res.iter().enumerate() { 
        println!("group: {i}");
        for (e, l) in list { 
            println!("\t{}\n\t{e}", l.join(" = "));
        }
        println!();
    }
}

fn compute_kh(name: &String, l: &Link) -> String { 
    println!("compute: {name}");

    let c = Ratio::from(0_i128);
    let kh = KhHomology::new(l, &c, &c, true);
    let grid = kh.gen_grid();

    kh.print_seq("i");

    let str = grid.iter().flat_map(|(idx, summand)| { 
        let r = summand.rank();
        if r > 0 { 
            Some(format!("{idx}:{r}"))
        } else { 
            None
        }
    }).join(",");
    str
}

fn make_dict(l: usize) -> Dict { 
    let mut dict = Dict::new();

    let p = usize::max(l, 10);
    for i in 3..=p { 
        for j in 1..=200 { 
            let name = format!("{i}_{j}");
            let Ok(l) = Link::load(&name) else { 
                break;
            };
            update_dict(&name, &l, &mut dict);
        }
    }

    for i in 11..=p { 
        for t in ["a", "n"] { 
            for j in 1..=500 { 
                let name = format!("K{i}{t}{j}");
                let Ok(l) = Link::load(&name) else { 
                    break;
                };
                update_dict(&name, &l, &mut dict);
            }
        }
    }

    dict.retain(|_, v| v.len() > 1);
    dict
}

fn update_dict(name: &String, l: &Link, dict: &mut Dict) { 
    let m_name = format!("m{name}");

    let a = compute_kh(name, &l);
    let b = compute_kh(&m_name, &l.mirror());

    if dict.contains_key(&a) { 
        dict.get_mut(&a).unwrap().push(name.clone());
    } else if dict.contains_key(&b){ 
        dict.get_mut(&b).unwrap().push(m_name.clone());
    } else {
        dict.insert(a.clone(), vec![name.clone()]);
    }

    if a == b { 
        dict.get_mut(&a).unwrap().push(m_name);
    }
}

fn dupl_list(dict: Dict) -> Vec<Vec<String>> { 
    dict.into_values().sorted_by_key(|v| name_to_code(&v[0])).collect()
}

fn distinguish(targets: &Vec<String>) -> Vec<(String, Vec<String>)> { 
    let mut res: HashMap<String, Vec<String>> = HashMap::new();
    
    for name in targets { 
        let l = load_link(name).unwrap();
        let e_str = compute_e_str(name, &l);
        if res.contains_key(&e_str) { 
            res.get_mut(&e_str).unwrap().push(name.clone());
        } else { 
            res.insert(e_str, vec![name.clone()]);
        }
    }

    res.into_iter().sorted_by_key(|(_, list)| name_to_code(&list[0])).collect()
}

fn compute_e_str(name: &String, l: &Link) -> String { 
    println!("compute: {name}");

    let c = Ratio::from(0_i128);
    let ckh = KhComplex::new_no_simplify(l, &c, &c, true);
    let kh = ckh.homology();
    let e_str = ckh.sl2_map(l).string_decomp(&kh);

    println!("{name}");
    kh.gen_grid().print_table("i", "j");
    println!("{e_str}\n");

    e_str.to_string()
}

fn load_link(name: &str) -> Result<Link, Box<dyn std::error::Error>> { 
    if name.starts_with("m") { 
        let name = String::from(&name[1..]);
        Link::load(&name).map(|l| l.mirror())
    } else {
        Link::load(name)
    }
}

fn name_to_code(name: &str) -> [usize; 4] { // [crossing-num, a/n, index, mirror] 
    let re = regex::Regex::new(r"^(m?)(\d+)_(\d+)$").unwrap();
    if let Some(caps) = re.captures(name) {
        let m = if &caps[1] == "m" { 1 } else { 0 };
        let n: usize = caps[2].parse().unwrap();
        let i: usize = caps[3].parse().unwrap();
        return [n, 0, i, m];
    }

    let re = regex::Regex::new(r"^(m)?K(\d+)([a|n])(\d+)$").unwrap();
    if let Some(caps) = re.captures(name) {
        let m = if caps.get(1).is_some() { 1 } else { 0 };
        let n: usize = caps[2].parse().unwrap();
        let t = if &caps[3] == "a" { 0 } else { 1 };
        let i: usize = caps[4].parse().unwrap();
        return [n, t, i, m];
    }

    panic!("Invalid name format: {}", name);
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_name_to_code() {
        let cases = vec![
            ("3_1",   [3, 0, 1, 0]),
            ("m4_2",  [4, 0, 2, 1]),
            ("K5a3",  [5, 0, 3, 0]),
            ("K6n4",  [6, 1, 4, 0]),
            ("mK7a5", [7, 0, 5, 1]),
            ("mK8n6", [8, 1, 6, 1]),
        ];

        for (input, expected) in cases {
            assert_eq!(name_to_code(input), expected);
        }
    }
}
