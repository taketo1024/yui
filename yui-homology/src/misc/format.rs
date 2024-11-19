use itertools::Itertools;
use yui::{Ring, RingOps};

pub fn rmod_str_symbol<R>(rank: usize, tors: &[R], dflt: &str) -> String
where R: Ring, for<'x> &'x R: RingOps<R> {
    use yui::util::format::superscript;
    make_string(
        R::math_symbol(), 
        rank, 
        &tors.iter().map(|t| t.to_string()).collect_vec(), 
        dflt, 
        superscript, 
        "⊕"
    )
}

#[cfg(feature = "tex")]
pub fn tex_rmod_str_symbol<R>(rank: usize, tors: &[R], dflt: &str) -> String
where R: Ring + yui::TeX, for<'x> &'x R: RingOps<R> {
    make_string(
        R::tex_math_symbol(), 
        rank, 
        &tors.iter().map(|t| t.tex_string()).collect_vec(), 
        dflt, 
        |a| format!("^{a}"), 
        "\\oplus"
    )
}

fn make_string<F>(symbol: String, rank: usize, tors: &[String], dflt: &str, superscript: F, oplus: &str) -> String
where F: Fn(usize) -> String {
    use std::collections::BTreeMap;

    if rank == 0 && tors.is_empty() { 
        return dflt.to_string()
    }

    let mut res = vec![];

    if rank > 1 {
        let str = format!("{}{}", symbol, superscript(rank));
        res.push(str);
    } else if rank == 1 { 
        res.push(symbol.clone());
    }

    let mut tors_acc = BTreeMap::<String, usize>::new();
    for t in tors { 
        if let Some(v) = tors_acc.get_mut(t) { 
            *v += 1;
        } else { 
            tors_acc.insert(t.clone(), 1);
        }
    }
    
    for (t, r) in tors_acc.iter() { 
        let str = if r > &1 { 
            format!("({}/{}){}", symbol, t, superscript(*r))
        } else { 
            format!("({}/{})", symbol, t)
        };
        res.push(str);
    }

    res.join(&format!(" {oplus} "))
}

#[cfg(test)]
mod tests {
    use super::*;
 
    #[test]
    fn zero() { 
        let s = rmod_str_symbol::<i32>(0, &[], "0");
        assert_eq!(s, "0");
    }

    #[test]
    fn rank1() { 
        let s = rmod_str_symbol::<i32>(1, &[], "0");
        assert_eq!(s, "Z");
    }

    #[test]
    fn rank2() { 
        let s = rmod_str_symbol::<i32>(2, &[], "0");
        assert_eq!(s, "Z²");
    }

    #[test]
    fn tor() { 
        let s = rmod_str_symbol(0, &[2], "0");
        assert_eq!(s, "(Z/2)");
    }

    #[test]
    fn tor2() { 
        let s = rmod_str_symbol(0, &[2,2,3], "0");
        assert_eq!(s, "(Z/2)² ⊕ (Z/3)");
    }

    #[test]
    fn mix() { 
        let s = rmod_str_symbol(2, &[2,2,3], "0");
        assert_eq!(s, "Z² ⊕ (Z/2)² ⊕ (Z/3)");
    }

    #[cfg(feature = "tex")]
    #[test]
    fn tex() { 
        let s = tex_rmod_str_symbol(2, &[2,2,3], "0");
        assert_eq!(s, "\\mathbb{Z}^2 \\oplus (\\mathbb{Z}/2)^2 \\oplus (\\mathbb{Z}/3)");
    }
}