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