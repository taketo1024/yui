use crate::app::utils::*;
use crate::app::err::*;
use std::marker::PhantomData;
use std::str::FromStr;
use yui::tex::TeX;
use yui::{EucRing, EucRingOps};
use yui_homology::DisplaySeq;
use yui_homology::{DisplayTable, GridTrait, SummandTrait};
use yui_kh::kh::ext::cc::KhChainMap;
use yui_kh::kh::KhComplex;
use yui_kh::kh::KhHomology;

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    dispatch_eucring!(App, args)
}

#[derive(Clone, Default, Debug, clap::Args)]
pub struct Args {
    pub link: String,

    #[arg(short = 'i', long, required = true)]
    pub cc_index: usize,

    #[arg(short = 'f', long, required = true)]
    pub map_type: usize,

    #[arg(short = 't', long, default_value = "Z")]
    pub c_type: CType,

    #[arg(short, long, default_value = "0")]
    pub c_value: String,

    #[arg(short, long)]
    pub mirror: bool,

    #[arg(short, long)]
    pub reduced: bool,

    #[arg(short = 'b', long)]
    pub reverse: bool,

    #[arg(short = 'g', long)]
    pub show_gens: bool,

    #[arg(long, default_value = "0")]
    pub log: u8,
}

pub struct App<R>
where
    R: EucRing + FromStr + TeX,
    for<'x> &'x R: EucRingOps<R>,
{
    args: Args,
    buff: String,
    _ring: PhantomData<R>
}

impl<R> App<R>
where
    R: EucRing + FromStr + TeX,
    for<'x> &'x R: EucRingOps<R>,
{
    pub fn new(args: Args) -> Self { 
        let buff = String::with_capacity(1024);
        App { args, buff, _ring: PhantomData }
    }

    pub fn run(&mut self) -> Result<String, Box<dyn std::error::Error>> {
        assert!([0, 1].contains(&self.args.map_type));

        let (h, t) = parse_pair::<R>(&self.args.c_value)?;

        if self.args.reduced {
            ensure!(t.is_zero(), "`t` must be zero for reduced.");
        }
    
        let r = self.args.reduced;
        let bigraded = (h.is_zero() && t.is_zero()) || 
            ["H", "0,T"].contains(&self.args.c_value.as_str());
        let gens = self.args.show_gens;
        
        let l = load_link(&self.args.link, self.args.mirror)?;
        let i = self.args.cc_index;

        let l = if self.args.reverse { 
            l.crossing_change(i)
        } else { 
            l
        };

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, r, i);
        let (f_name, f) = if self.args.map_type == 0 { 
            ("f0", KhComplex::cc_map0(&c1, &c2, i)) 
        } else { 
            ("f1", KhComplex::cc_map1(&c1, &c2, &l, i))
        };

        let (h1, h2) = (c1.homology(), c2.homology());

        self.show_table("from:", &h1, bigraded, gens);
        self.show_table("to:",   &h2, bigraded, gens);
        self.show_map(f_name, &h1, &h2, &f);
        
        let res = self.flush();
        Ok(res)
    }

    fn show_table(&mut self, label: &str, h: &KhHomology<R>, bigraded: bool, with_gens: bool) { 
        let table = if bigraded { 
            h.gen_grid().display_table("i", "j")
        } else { 
            h.display_seq("i")
        };

        self.out(label);
        self.out(&table);

        if with_gens { 
            self.show_gens(h);
        }
    }

    fn show_gens(&mut self, h: &KhHomology<R>) { 
        for i in h.support() {
            if h[i].is_zero() { continue }

            self.out(&format!("({i}): {}", h[i]));

            for (k, z) in h[i].gens().enumerate() { 
                self.out(&format!("  {k}: {z}"));
            }
            self.out("");
        }
    }

    fn show_map(&mut self, f_name: &str, h1: &KhHomology<R>, h2: &KhHomology<R>, f: &KhChainMap<R>) { 
        self.out(&format!("{f_name}: deg {}\n", f.deg()));

        for i in h1.h_range() { 
            let j = i + f.deg();
            self.out(&format!("({i}) {} -> ({j}) {}", h1[i], h2[j]));

            for z in h1[i].gens() { 
                let w = f.apply(i, &z);
                let x = h1[i].vectorize_euc(&z).into_vec();
                let y = h2[j].vectorize_euc(&w).into_vec();
                self.out(&format!("\t{:?} -> {:?}", x, y));
            }
            self.out("");
        }
    }

    fn out(&mut self, str: &str) { 
        self.buff.push_str(str);
        self.buff.push('\n');
    }

    fn flush(&mut self) -> String { 
        let res = std::mem::take(&mut self.buff);
        res.trim_end().to_string()
    }
}

// #[cfg(_test)]
// mod tests {
//     use super::*;

//     #[test]
//     fn test1() {
//         let args = Args {
//             link: "3_1".to_string(),
//             c_value: "0".to_string(),
//             ..Default::default()
//         };
//         let res = dispatch(&args);
//         assert!(res.is_ok());
//     }

//     #[test]
//     fn test2() {
//         let args = Args {
//             link: "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]".to_string(),
//             c_value: "2".to_string(),
//             mirror: true,
//             reduced: true,
//             show_alpha: true,
//             ..Default::default()
//         };
//         let res = dispatch(&args);
//         assert!(res.is_ok());
//     }

//     #[cfg(feature = "poly")]
//     mod poly_tests {
//         use super::*;

//         #[test]
//         fn test_zpoly_h() {
//             let args = Args {
//                 link: "3_1".to_string(),
//                 c_value: "H".to_string(),
//                 c_type: CType::Z,
//                 ..Default::default()
//             };
//             let res = dispatch(&args);
//             assert!(res.is_ok());
//         }

//         #[test]
//         fn test_zpoly_t() {
//             let args = Args {
//                 link: "3_1".to_string(),
//                 c_value: "0,T".to_string(),
//                 c_type: CType::Z,
//                 ..Default::default()
//             };
//             let res = dispatch(&args);
//             assert!(res.is_ok());
//         }

//         #[test]
//         fn test_zpoly_ht() {
//             let args = Args {
//                 link: "3_1".to_string(),
//                 c_value: "H,T".to_string(),
//                 c_type: CType::Z,
//                 ..Default::default()
//             };
//             let res = dispatch(&args);
//             assert!(res.is_ok());
//         }
//     }
// }
