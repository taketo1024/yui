use itertools::Itertools;
use yui_core::{Ring, RingOps};

pub fn r_mod_str<'a, R, Itr>(rank: usize, tors: Itr) -> String
where R: Ring, for<'x> &'x R: RingOps<R>, Itr: Iterator<Item = &'a R> {
    use yui_utils::superscript;

    let tors = tors
        .into_group_map_by(|r| r.to_string())
        .into_iter().map(|(k, list)| (k, list.len()))
        .collect_vec();

    if rank == 0 && tors.is_empty() { 
        return ".".to_string()
    }

    let mut res = vec![];
    let symbol = R::math_symbol();

    if rank > 1 {
        let str = format!("{}{}", symbol, superscript(rank as isize));
        res.push(str);
    } else if rank == 1 { 
        let str = format!("{}", symbol);
        res.push(str);
    }
    
    for (t, r) in tors.iter() { 
        let str = if r > &1 { 
            format!("({}/{}){}", symbol, t, superscript(*r as isize))
        } else { 
            format!("({}/{})", symbol, t)
        };
        res.push(str);
    }

    let str = res.join(" ⊕ ");
    str
}