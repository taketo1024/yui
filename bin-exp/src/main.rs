use std::collections::HashMap;

use itertools::Itertools;
use yui_core::num::Ratio;
use yui_homology::{DisplaySeq, SummandTrait};
use yui_kh::kh::KhHomology;
use yui_link::Link;

type Dict = HashMap<String, Vec<String>>;
fn main() {
    let dict = make_dict(10);
    let list = dupl_list(dict);
    println!("{list:?}");
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
    dict.into_values().sorted_by_key(|v| v.first().unwrap().clone()).collect()
}