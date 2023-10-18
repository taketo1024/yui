use std::str::FromStr;
use yui_core::{EucRing, EucRingOps};
use yui_homology::v2::{PrintTable, PrintSeq};
use yui_khovanov::{KhHomology, KhHomologyBigraded};
use crate::utils::*;

#[derive(Debug, clap::Args)]
pub struct Args { 
    link: String,

    #[arg(short, long, default_value = "0")]
    c_value: String,

    #[arg(short = 't', long, default_value = "Z")]
    c_type: CType,

    #[arg(short, long)]
    mirror: bool,

    #[arg(short, long)]
    reduced: bool,

    #[arg(short, long)]
    bigraded: bool,

    #[arg(long)]
    old: bool,

    #[arg(long)]
    pub debug: bool
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
    let l = load_link(&args.link, args.mirror)?;
    if args.c_value != "0" { 
        return err!("--bigraded only supported for `c = 0`.")
    }

    let kh = if args.old { 
        KhHomologyBigraded::new_v1(l, args.reduced)
    } else {
        KhHomologyBigraded::new(l, args.reduced)
    };

    let table = kh.display_table();
    Ok(table)
}

fn compute_homology<R>(args: &Args) -> Result<String, Box<dyn std::error::Error>>
where R: EucRing + FromStr, for<'x> &'x R: EucRingOps<R> { 
    let l = load_link(&args.link, args.mirror)?;
    let (h, t) = parse_pair::<R>(&args.c_value)?;

    if args.reduced && !t.is_zero() { 
        return err!("{t} != 0 is not allowed for reduced.");
    }

    let kh = if args.old { 
        KhHomology::new_v1(&l, &h, &t, args.reduced)
    } else {
        KhHomology::new(&l, &h, &t, args.reduced)
    };

    let res = kh.display_seq();
    Ok(res)
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn test1() { 
        let args = Args { 
            link: "3_1".to_string(), 
            c_value: "0".to_string(), 
            c_type: CType::Z, 
            mirror: false, 
            reduced: false, 
            bigraded: false,
            old: false,
            debug: false
        };
        let res = run(&args);
        assert!(res.is_ok());
    }

    #[test]
    fn test2() { 
        let args = Args { 
            link: "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]".to_string(),
            c_value: "0".to_string(),
            c_type: CType::Z,
            mirror: true,
            reduced: true,
            bigraded: true,
            old: false,
            debug: false
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
                link: "3_1".to_string(),
                c_value: "H".to_string(),
                c_type: CType::Q,
                mirror: false,
                reduced: false,
                bigraded: false,
                old: false,
                debug: false
            };
            let res = run(&args);
            assert!(res.is_ok());
        }

        #[test]
        fn test_qpoly_t() { 
            let args = Args {
                link: "3_1".to_string(),
                c_value: "0,T".to_string(),
                c_type: CType::Q,
                mirror: false,
                reduced: false,
                bigraded: false,
                old: false,
                debug: false
            };
            let res = run(&args);
            assert!(res.is_ok());
        }
    }
}