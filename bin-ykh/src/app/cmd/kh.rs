use std::marker::PhantomData;
use std::str::FromStr;
use yui::tex::TeX;
use yui::{EucRing, EucRingOps};
use yui_homology::{DisplaySeq, DisplayTable, GridTrait, SummandTrait, tex::TeXTable};
use yui_kh::kh::KhHomology;
use yui_kh::kh::KhChainExt;
use yui_link::Link;
use crate::app::utils::*;
use crate::app::err::*;

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    dispatch_eucring!(App, args)
}

#[derive(Clone, Default, Debug, clap::Args)]
pub struct Args { 
    pub link: String,

    #[arg(short = 't', long, default_value = "Z")]
    pub c_type: CType,

    #[arg(short, long, default_value = "0")]
    pub c_value: String,

    #[arg(short, long)]
    pub mirror: bool,

    #[arg(short, long)]
    pub reduced: bool,

    #[arg(short = 'g', long)]
    pub show_gens: bool,

    #[arg(short = 'a', long)]
    pub show_alpha: bool,

    #[arg(short = 's', long)]
    pub show_ss: bool,

    #[arg(short = 'n', long)]
    pub no_simplify: bool,

    #[arg(short, long, default_value = "unicode")]
    pub format: Format,

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
        let (h, t) = parse_pair::<R>(&self.args.c_value)?;
    
        if self.args.reduced { 
            ensure!(t.is_zero(), "`t` must be zero for reduced.");
        }
        if self.args.show_alpha { 
            ensure!(t.is_zero(), "`t` must be zero to have alpha.");
        }
        if self.args.show_ss { 
            ensure!(!h.is_zero() && !h.is_unit(), "`h` must be non-zero, non-invertible to compute ss.");
            ensure!(t.is_zero(), "`t` must be zero to compute ss.");
        }
    
        let bigraded = (h.is_zero() && t.is_zero()) || 
            ["H", "0,T"].contains(&self.args.c_value.as_str());
    
        let l = load_link(&self.args.link, self.args.mirror)?;
        
        let kh = if self.args.no_simplify {
            KhHomology::new_no_simplify(&l, &h, &t, self.args.reduced)
        } else { 
            KhHomology::new(&l, &h, &t, self.args.reduced)
        } ;

        // print Kh
        let table = if bigraded { 
            let grid = kh.gen_grid();
            match self.args.format {
                Format::Unicode => grid.display_table("i", "j"),
                Format::TeX     => grid.tex_table("$\\mathit{Kh}$", " ")
            }
        } else { 
            kh.display_seq("i")
        };
        self.out(&table);

        if self.args.show_gens { 
            self.show_gens(&kh);
        }

        if self.args.show_alpha { 
            self.show_alpha(&kh);
        }

        if self.args.show_ss { 
            self.show_ss(&l, &h, &kh)?;
        }
    
        Ok(self.flush())
    }

    fn show_gens(&mut self, kh: &KhHomology<R>) { 
        for i in kh.support() {
            let h = &kh[i];
            if h.is_zero() { continue }

            self.out(&format!("Kh[{i}]: {}", h));

            let r = h.rank() + h.tors().len();
            for i in 0..r { 
                let z = h.gen(i);
                self.out(&format!("  {i}: {z}"));
            }
            self.out("");
        }
    }

    fn show_alpha(&mut self, kh: &KhHomology<R>) { 
        let zs = kh.canon_cycles();
        for (i, z) in zs.iter().enumerate() { 
            let h = z.h_deg();
            let v = kh[h].vectorize_euc(z);
            self.out(&format!("a[{i}] in Kh[{h}]: {}", vec2str(&v)));
            self.out(&format!("  {z}\n"));
        }
    }

    fn show_ss(&mut self, l: &Link, c: &R, kh: &KhHomology<R>) -> Result<(), Box<dyn std::error::Error>> { 
        assert!(!c.is_unit() && !c.is_unit());

        use yui_kh::misc::div_vec;

        let w = l.writhe();
        let r = l.seifert_circles().len() as i32;
        let zs = kh.canon_cycles();

        for (i, z) in zs.iter().enumerate() { 
            let h = &kh[z.h_deg()];
            let v = h.vectorize(z).subvec(0..h.rank());
            let d = div_vec(&v, c);

            ensure!(d.is_some(), "invalid div for a = {z}.");

            let d = d.unwrap();
            let s = 2 * d + w - r + 1;

            self.out(&format!("ss[{i}] = {s} (d = {d}, w = {w}, r = {r})"));
        }

        Ok(())
    }

    fn out(&mut self, str: &str) { 
        self.buff.push_str(str);
        self.buff.push('\n');
    }

    fn flush(&mut self) -> String { 
        let res = std::mem::take(&mut self.buff);
        res.trim().to_string()
    }
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn test1() { 
        let args = Args { 
            link: "3_1".to_string(), 
            c_value: "0".to_string(), 
            ..Default::default()
        };
        let res = dispatch(&args);
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
            ..Default::default()
        };
        let res = dispatch(&args);
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
                ..Default::default()
            };
            let res = dispatch(&args);
            assert!(res.is_ok());
        }

        #[test]
        fn test_qpoly_t() { 
            let args = Args {
                link: "3_1".to_string(),
                c_value: "0,T".to_string(),
                c_type: CType::Q,
                ..Default::default()
            };
            let res = dispatch(&args);
            assert!(res.is_ok());
        }
    }
}