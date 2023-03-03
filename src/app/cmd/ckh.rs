use std::str::FromStr;
use itertools::Itertools;
use yui_core::{Ring, RingOps};
use yui_homology::utils::ChainReducer;
use yui_link::Link;
use yui_matrix::dense::*;
use yui_homology::RModGrid;
use yui_khovanov::KhComplex;
use crate::utils::*;

pub fn run(name: String, link: Option<String>, c_value: String, c_type: CType, mirror: bool, reduced: bool, with_alpha: bool) -> Result<String, Box<dyn std::error::Error>> {
    let mut l = load_link(&name, &link)?;
    if mirror { 
        l = l.mirror();
    }

    dispatch_ring!(c_type, describe_ckh, &l, &c_value, reduced, with_alpha)
}

fn describe_ckh<R>(l: &Link, c_value: &String, reduced: bool, with_alpha: bool) -> Result<String, Box<dyn std::error::Error>>
where R: Ring + FromStr, for<'x> &'x R: RingOps<R> { 
    use string_builder::Builder;

    let (h, t) = parse_pair::<R>(c_value)?;
    if reduced && !t.is_zero() { 
        return err!("{t} != 0 is not allowed for reduced.");
    }
    
    let ckh = KhComplex::new(l.clone(), h, t, reduced);
    let mut red = ChainReducer::new(&ckh);

    if with_alpha { 
        for z in ckh.canon_cycles() {
            let v = ckh[0].vectorize(&z);
            red.set_vec(0, v);
        }
    }

    red.process();

    let mut b = Builder::new(1024);
    let symbol = R::set_symbol();
    let superscript = |n: usize| yui_utils::superscript(n as isize);

    b.append(format!("{}\n", ckh.to_string()));

    for i in ckh.indices() {
        let d = red.matrix(i).unwrap();
        let (m, n) = d.shape();
        b.append(format!("d[{}]: {}{} -> {}{}\n\n", i, symbol, superscript(m), symbol, superscript(n)));
        b.append(d.to_dense().to_string());
        b.append("\n\n");
    }

    if with_alpha { 
        let vs = red.vecs(0).unwrap();
        for (i, v) in vs.iter().enumerate() {
            b.append(format!("a[{i}] = [{}]\n", v.to_dense().iter().map(|r| r.to_string()).collect_vec().join(", ")));
        }
        b.append("\n");
    }

    let res = b.string()?;
    Ok(res)
}