use std::str::FromStr;
use yui_core::{EucRing, EucRingOps};
use yui_link::Link;
use yui_homology::PrintTable;
use yui_khovanov::{KhHomology, KhHomologyBigraded};
use crate::utils::*;

pub fn run(name: String, link: Option<String>, c_value: String, c_type: CType, mirror: bool, reduced: bool, bigraded: bool) -> Result<String, Box<dyn std::error::Error>> {
    let mut l = load_link(&name, &link)?;
    if mirror { 
        l = l.mirror();
    }

    if bigraded { 
        dispatch_eucring!(c_type, compute_bigraded, &l, &c_value, reduced)
    } else { 
        dispatch_eucring!(c_type, compute_homology, &l, &c_value, reduced)
    }
}

fn compute_bigraded<R>(l: &Link, c_value: &String, reduced: bool) -> Result<String, Box<dyn std::error::Error>>
where R: EucRing + FromStr, for<'x> &'x R: EucRingOps<R> { 
    if c_value.as_str() != "0" { 
        return err!("--bigraded only supported for `c = 0`.")
    }
    let h = KhHomologyBigraded::new(l.clone(), reduced);
    let table = h.table();
    Ok(table)
}

fn compute_homology<R>(l: &Link, c_value: &String, reduced: bool) -> Result<String, Box<dyn std::error::Error>>
where R: EucRing + FromStr, for<'x> &'x R: EucRingOps<R> { 
    let (h, t) = parse_pair::<R>(c_value)?;
    if reduced && !t.is_zero() { 
        return err!("{t} != 0 is not allowed for reduced.");
    }

    let kh = KhHomology::new(l.clone(), h, t, reduced);
    let res = kh.to_string();
    Ok(res)
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn test1() { 
        let name = "3_1".to_string();
        let c_value = "0".to_string();
        let c_type = CType::Z;
        let mirror = false;
        let reduced = false;
        let bigraded = false;

        let res = run(name, None, c_value, c_type, mirror, reduced, bigraded);
        assert!(res.is_ok());
    }

    #[test]
    fn test2() { 
        let name = "4_1".to_string();
        let c_value = "0".to_string();
        let c_type = CType::Z;
        let mirror = true;
        let reduced = true;
        let bigraded = true;

        let res = run(name, None, c_value, c_type, mirror, reduced, bigraded);
        assert!(res.is_ok());
    }
}