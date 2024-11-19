use yui::{Ring, RingOps};

pub fn rmod_str_symbol<R>(rank: usize, tors: &[R], dflt: &str) -> String
where R: Ring, for<'x> &'x R: RingOps<R> {
    use std::collections::BTreeMap;
    use yui::util::format::superscript;

    if rank == 0 && tors.is_empty() { 
        return dflt.to_string()
    }

    let mut res = vec![];
    let symbol = R::math_symbol();

    if rank > 1 {
        let str = format!("{}{}", symbol, superscript(rank as isize));
        res.push(str);
    } else if rank == 1 { 
        let str = symbol.to_string();
        res.push(str);
    }

    let mut tors_acc = BTreeMap::<String, usize>::new();
    for t in tors { 
        let t_str = t.to_string();
        if let Some(v) = tors_acc.get_mut(&t_str) { 
            *v += 1;
        } else { 
            tors_acc.insert(t_str, 1);
        }
    }
    
    for (t, r) in tors_acc.iter() { 
        let str = if r > &1 { 
            format!("({}/{}){}", symbol, t, superscript(*r as isize))
        } else { 
            format!("({}/{})", symbol, t)
        };
        res.push(str);
    }

    res.join(" ⊕ ")
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
}