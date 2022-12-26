use itertools::Itertools;
use crate::math::traits::{Ring, RingOps};

use super::digits::Digits;

pub fn subscript(i: isize) -> String {
    if i == 0 { return String::from('\u{2080}') }

    let (mut res, i) = if i > 0 { 
        (vec![], i as usize)
    } else { 
        (vec!['\u{208B}'], -i as usize)
    };

    let mut s = i.digits().into_iter().map(|d| {
        char::from_u32(('\u{2080}' as u32) + (d as u32)).unwrap()
    }).collect_vec();

    res.append(&mut s);
    res.iter().collect()
}

pub fn superscript(i: isize) -> String {
    if i == 0 { return String::from('\u{2070}') }

    let (mut res, i) = if i > 0 { 
        (vec![], i as usize)
    } else { 
        (vec!['\u{207B}'], -i as usize)
    };

    let mut s = i.digits().into_iter().map(|d| {
        match d { 
            1 => '\u{00B9}',
            2 => '\u{00B2}',
            3 => '\u{00B3}',
            _ => char::from_u32(('\u{2070}' as u32) + (d as u32)).unwrap()
        }
    }).collect_vec();

    res.append(&mut s);
    res.iter().collect()
}

pub fn module_descr<R>(rank: usize, tors: &Vec<R>) -> String 
where R: Ring, for<'x> &'x R: RingOps<R>  
{ 
    use grouping_by::GroupingBy;
    if rank == 0 && tors.is_empty() { 
        return String::from("0")
    }

    let mut res = vec![];
    let symbol = R::symbol();

    if rank > 1 {
        let str = format!("{}{}", symbol, superscript(rank as isize));
        res.push(str);
    } else if rank == 1 { 
        let str = format!("{}", symbol);
        res.push(str);
    }

    for (t, r) in tors.iter().counter(|&t| t) { 
        let str = if r > 1 { 
            format!("({}/{}){}", symbol, t, superscript(r as isize))
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
    fn test_subscript() { 
        assert_eq!(subscript(0), "₀");
        assert_eq!(subscript(1234567890), "₁₂₃₄₅₆₇₈₉₀");
        assert_eq!(subscript(-1234567890), "₋₁₂₃₄₅₆₇₈₉₀");
    }

    #[test]
    fn test_superscript() { 
        assert_eq!(superscript(0), "⁰");
        assert_eq!(superscript(1234567890), "¹²³⁴⁵⁶⁷⁸⁹⁰");
        assert_eq!(superscript(-1234567890), "⁻¹²³⁴⁵⁶⁷⁸⁹⁰");
    }
}
