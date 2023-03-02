use std::str::FromStr;
use yui_core::{Ring, RingOps};
use yui_link::Link;
use yui_matrix::dense::*;
use yui_homology::{RModGrid, ChainComplex, Reduced};
use yui_khovanov::KhComplex;
use crate::utils::*;

pub fn run(name: String, link: Option<String>, c_value: String, c_type: CType, mirror: bool, reduced: bool) -> Result<String, Box<dyn std::error::Error>> {
    let mut l = load_link(&name, &link)?;
    if mirror { 
        l = l.mirror();
    }

    dispatch_ring!(describe_ckh, &c_type, &l, &c_value, reduced)
}

fn describe_ckh<R>(l: &Link, c_value: &String, reduced: bool) -> Result<String, Box<dyn std::error::Error>>
where R: Ring + FromStr, for<'x> &'x R: RingOps<R> { 
    use string_builder::Builder;

    let (h, t) = parse_pair::<R>(c_value)?;
    if reduced && !t.is_zero() { 
        return err!("{t} != 0 is not allowed for reduced.");
    }
    
    let ckh = KhComplex::new(l.clone(), h, t, reduced);
    let ckh = Reduced::from(ckh);

    let mut b = Builder::new(1024);
    for i in ckh.indices() {
        b.append(format!("C[{}]: {} -> {}\n", i, ckh[i], ckh[i+1]));

        let d = ckh.d_matrix(i);
        if d.rows() > 0 && d.cols() > 0 {
            b.append(d.to_dense().to_string());
            b.append("\n");
        }
        b.append("\n");
    }

    let res = b.string()?;
    Ok(res)
}