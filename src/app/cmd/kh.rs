use std::str::FromStr;
use yui_core::{EucRing, EucRingOps};
use yui_homology::fmt::FmtTable;
use yui_khovanov::{KhHomology, KhHomologyBigraded};
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

    #[arg(short, long)]
    bigraded: bool,
}

pub fn run(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    if args.bigraded { 
        dispatch_eucring!(&args.c_value, &args.c_type, compute_bigraded, args)
    } else { 
        dispatch_eucring!(&args.c_value, &args.c_type, compute_homology, args)
    }
}

fn compute_bigraded<R>(args: &Args) -> Result<String, Box<dyn std::error::Error>>
where R: EucRing + FromStr, for<'x> &'x R: EucRingOps<R> { 
    let l = load_link(&args.name, &args.link, args.mirror)?;
    if args.c_value != "0" { 
        return err!("--bigraded only supported for `c = 0`.")
    }

    let h = KhHomologyBigraded::new(l, args.reduced);
    let table = h.table_string();
    Ok(table)
}

fn compute_homology<R>(args: &Args) -> Result<String, Box<dyn std::error::Error>>
where R: EucRing + FromStr, for<'x> &'x R: EucRingOps<R> { 
    let l = load_link(&args.name, &args.link, args.mirror)?;
    let (h, t) = parse_pair::<R>(&args.c_value)?;

    if args.reduced && !t.is_zero() { 
        return err!("{t} != 0 is not allowed for reduced.");
    }

    let kh = KhHomology::new(l.clone(), &h, &t, args.reduced);
    let res = kh.to_string();
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
            bigraded: false
        };
        let res = run(&args);
        assert!(res.is_ok());
    }

    #[test]
    fn test2() { 
        let args = Args { 
            name: "".to_string(),
            link: Some("[[1,4,2,5],[3,6,4,1],[5,2,6,3]]".to_string()),
            c_value: "0".to_string(),
            c_type: CType::Z,
            mirror: true,
            reduced: true,
            bigraded: true
        };
        let res = run(&args);
        assert!(res.is_ok());
    }

    #[cfg(feature = "poly")]
    mod poly_tests { 
        use super::*;
        
        #[test]
        fn test_qpoly_h() { 
            let args = Args {
                name: "3_1".to_string(),
                link: None,
                c_value: "H".to_string(),
                c_type: CType::Q,
                mirror: false,
                reduced: false,
                bigraded: false
            };
            let res = run(&args);
            assert!(res.is_ok());
        }

        #[test]
        fn test_qpoly_t() { 
            let args = Args {
                name: "3_1".to_string(),
                link: None,
                c_value: "0,T".to_string(),
                c_type: CType::Q,
                mirror: false,
                reduced: false,
                bigraded: false
            };
            let res = run(&args);
            assert!(res.is_ok());
        }
    }
}