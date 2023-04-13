use std::str::FromStr;
use itertools::Itertools;
use yui_core::{Ring, RingOps};
use yui_homology::ChainComplex;
use yui_homology::Grid;
use yui_khovanov::KhComplex;
use crate::utils::*;

#[derive(Debug, clap::Args)]
pub struct Args { 
    name: String,

    #[arg(short, long)]
    link: Option<String>,

    #[arg(short, long, default_value = "0")]
    c_value: String,

    #[arg(short = 't', long, default_value = "z")]
    c_type: CType,

    #[arg(short, long)]
    mirror: bool,

    #[arg(short, long)]
    reduced: bool,

    #[arg(short = 'a', long)]
    with_alpha: bool,

    #[arg(short = 'n', long)]
    no_simplify: bool,
}

pub fn run(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    dispatch_ring!(&args.c_value, &args.c_type, describe_ckh, args)
}

fn describe_ckh<R>(args: &Args) -> Result<String, Box<dyn std::error::Error>>
where R: Ring + FromStr, for<'x> &'x R: RingOps<R> { 
    use string_builder::Builder;
    
    let (h, t) = parse_pair::<R>(&args.c_value)?;
    if args.reduced && !t.is_zero() { 
        return err!("{t} != 0 is not allowed for reduced.");
    }
    
    let l = load_link(&args.name, &args.link, args.mirror)?;
    let c = KhComplex::new(l, h, t, args.reduced);
    
    let vs = if args.with_alpha { 
        c.canon_cycles().into_iter().map(|z| 
            (0, c[0].vectorize(&z))
        ).collect_vec()
    } else { 
        vec![]
    };

    let (c, vs) = if args.no_simplify { 
        (c.as_generic(), vs)
    } else { 
        c.as_generic().simplify_with(vs)
    };

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

    if args.with_alpha { 
        b.append(format!("----------\nLee cycles\n----------\n\n"));
        for (i, (_, v)) in vs.iter().enumerate() {
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
        let args = Args {
        	name: "3_1".to_string(),
            link: None,
        	c_value: "0".to_string(),
        	c_type: CType::Z,
        	mirror: false,
        	reduced: false,
        	with_alpha: false,
            no_simplify: false
        };
        let res = run(&args);
        assert!(res.is_ok());
    }

    #[test]
    fn test2() { 
        let args = Args {
        	name: "".to_string(),
            link: Some("[[1,4,2,5],[3,6,4,1],[5,2,6,3]]".to_string()),
        	c_value: "2".to_string(),
        	c_type: CType::Z,
        	mirror: true,
        	reduced: true,
        	with_alpha: true,
            no_simplify: false
        };
        let res = run(&args);
        assert!(res.is_ok());
    }

    #[cfg(feature = "poly")]
    mod poly_tests { 
        use super::*;
        
        #[test]
        fn test_zpoly_h() { 
            let args = Args {
                name: "3_1".to_string(),
                link: None,
                c_value: "H".to_string(),
                c_type: CType::Z,
                mirror: false,
                reduced: false,
                with_alpha: false,
                no_simplify: false
            };
            let res = run(&args);
            assert!(res.is_ok());
        }

        #[test]
        fn test_zpoly_t() { 
            let args = Args {
                name: "3_1".to_string(),
                link: None,
                c_value: "T".to_string(),
                c_type: CType::Z,
                mirror: false,
                reduced: false,
                with_alpha: false,
                no_simplify: false
            };
            let res = run(&args);
            assert!(res.is_ok());
        }

        #[test]
        fn test_zpoly_ht() { 
            let args = Args {
                name: "3_1".to_string(),
                link: None,
                c_value: "H,T".to_string(),
                c_type: CType::Z,
                mirror: false,
                reduced: false,
                with_alpha: false,
                no_simplify: false
            };
            let res = run(&args);
            assert!(res.is_ok());
        }
    }
}