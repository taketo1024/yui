use std::str::FromStr;
use log::info;
use yui_core::{EucRing, EucRingOps};
use yui_khovanov::{ss_invariant, ss_invariant_v1};
use crate::utils::*;

#[derive(Debug, clap::Args)]
pub struct Args { 
    pub link: String,

    #[arg(short, long)]
    pub c_value: String,

    #[arg(short = 't', long, default_value = "z")]
    pub c_type: CType,

    #[arg(short, long)]
    pub mirror: bool,

    #[arg(short = 'n', long, default_value = "1")]
    pub order: usize,

    #[arg(short, long)]
    pub reduced: bool,

    #[arg(long)]
    pub old: bool,

    #[arg(long)]
    pub debug: bool
}

pub fn run(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    info!("compute ss: {}, c = {}", args.link, args.c_value);

    guard_panic(|| 
        dispatch_eucring!(&args.c_value, &args.c_type, compute_ss, args)
        .map(|ss| ss.to_string())
    )
}

fn compute_ss<R>(args: &Args) -> Result<i32, Box<dyn std::error::Error>>
where R: EucRing + FromStr, for<'x> &'x R: EucRingOps<R> { 
    let l = load_link(&args.link, args.mirror)?;
    let Ok(c) = R::from_str(&args.c_value) else { 
        return err!("cannot parse c: '{}' as type: {}.", &args.c_value, std::any::type_name::<R>())
    };

    let ss = if args.old { 
        ss_invariant_v1(&l, &c, args.order, args.reduced)
    } else { 
        ss_invariant(&l, &c, args.order, args.reduced)
    };

    Ok(ss)
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn test1() { 
        let args = Args {
            link: "3_1".to_string(),
        	c_value: "2".to_string(),
        	c_type: CType::Z,
            mirror: false,
            order: 1,
            reduced: true,
            old: false,
            debug: false
        };
        let res = run(&args);

        assert!(res.is_ok());
        assert_eq!(res.unwrap(), "-2".to_string());
    }

    #[test]
    fn test2() { 
        let args = Args {
        	link: "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]".to_string(),
        	c_value: "3".to_string(),
        	c_type: CType::Z,
            mirror: false,
            order: 1,
            reduced: true,
            old: false,
            debug: false
        };
        let res = run(&args);

        assert!(res.is_ok());
        assert_eq!(res.unwrap(), "-2".to_string());
    }

    #[cfg(feature = "poly")]
    mod poly_tests { 
        use super::*;
        #[test]
        fn test1() { 
            let args = Args {
                link: "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]".to_string(),
                c_value: "H".to_string(),
                c_type: CType::Q,
                mirror: false,
                order: 1,
                reduced: true,
                old: false,
                debug: false
            };
            let res = run(&args);
            
            assert!(res.is_ok());
            assert_eq!(res.unwrap(), "-2".to_string());
        }    
    }
}