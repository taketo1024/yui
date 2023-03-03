use std::str::FromStr;
use itertools::Itertools;
use yui_core::{Ring, RingOps};
use yui_homology::{ChainComplex, GenericChainComplex, utils::ChainReducer};
use yui_link::Link;
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
    
    let c = KhComplex::new(l.clone(), h, t, reduced);
    let mut red = ChainReducer::new(&c);

    if with_alpha { 
        for z in c.canon_cycles() {
            let v = c[0].vectorize(&z);
            red.set_vec(0, v);
        }
    }

    red.process();

    let c = GenericChainComplex::generate(
        c.indices(), 
        c.d_degree(), 
        |i| Some(red.take_matrix(i))
    );
    let vs = red.take_vecs(0);

    let mut b = Builder::new(1024);

    b.append(format!("-------\nComplex\n-------\n\n"));
    b.append(format!("{}\n", c.to_string()));

    b.append(format!("-------------\nDifferentials\n-------------\n\n"));
    for i in c.indices() {
        let d = c.d_matrix(i);
        b.append(format!("d[{}]: {} -> {}\n\n", i, c[i], c[i+1]));
        b.append(d.to_dense().to_string());
        b.append("\n\n");
    }

    if with_alpha { 
        b.append(format!("----------\nLee cycles\n----------\n\n"));
        for (i, v) in vs.iter().enumerate() {
            b.append(format!("a[{i}] = [{}]\n", v.to_dense().iter().map(|r| r.to_string()).collect_vec().join(", ")));
        }
    }

    let res = b.string()?;
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
        let with_alpha = false;

        let res = run(name, None, c_value, c_type, mirror, reduced, with_alpha);
        assert!(res.is_ok());
    }

    #[test]
    fn test2() { 
        let name = "4_1".to_string();
        let c_value = "2".to_string();
        let c_type = CType::Z;
        let mirror = true;
        let reduced = true;
        let with_alpha = true;

        let res = run(name, None, c_value, c_type, mirror, reduced, with_alpha);
        assert!(res.is_ok());
    }
}