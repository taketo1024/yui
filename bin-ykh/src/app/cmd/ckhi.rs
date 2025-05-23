use crate::app::utils::*;
use crate::app::err::*;
use std::marker::PhantomData;
use std::str::FromStr;
use yui::tex::TeX;
use yui::{Ring, RingOps};
use yui_homology::DisplayTable;
use yui_homology::{ChainComplexTrait, GridTrait, SummandTrait, tex::TeXTable};
use yui_kh::kh::KhChainExt;
use yui_kh::khi::KhIComplex;

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    dispatch_ring!(App, args)
}

#[derive(Clone, Default, Debug, clap::Args)]
pub struct Args {
    pub link: String,

    #[arg(short = 't', long, default_value = "F2")]
    pub c_type: CType,

    #[arg(short, long, default_value = "0")]
    pub c_value: String,

    #[arg(short, long)]
    pub mirror: bool,

    #[arg(short, long)]
    pub reduced: bool,

    #[arg(short = 'g', long)]
    pub show_gens: bool,

    #[arg(short = 'd', long)]
    pub show_diff: bool,

    #[arg(short = 'a', long)]
    pub show_alpha: bool,

    #[arg(short, long, default_value = "unicode")]
    pub format: Format,

    #[arg(long, default_value = "0")]
    pub log: u8,
}

pub struct App<R>
where
    R: Ring + FromStr + TeX,
    for<'x> &'x R: RingOps<R>,
{
    args: Args,
    buff: String,
    _ring: PhantomData<R>
}

impl<R> App<R>
where
    R: Ring + FromStr + TeX,
    for<'x> &'x R: RingOps<R>,
{
    pub fn new(args: Args) -> Self { 
        let buff = String::with_capacity(1024);
        App { args, buff, _ring: PhantomData }
    }

    pub fn run(&mut self) -> Result<String, Box<dyn std::error::Error>> {
        let (h, t) = parse_pair::<R>(&self.args.c_value)?;
    
        ensure!(self.args.c_type == CType::F2, "Only `-t F2` is supported.");

        if self.args.reduced {
            ensure!(t.is_zero(), "`t` must be zero for reduced.");
        }
        if self.args.show_alpha { 
            ensure!(t.is_zero(), "`t` must be zero to have alpha.");
        }
    
        let l = load_sinv_knot(&self.args.link, self.args.mirror)?;
        let ckhi = KhIComplex::new(&l, &h, &t, self.args.reduced);
        
        // CKh generators
        let grid = ckhi.gen_grid();
        let table = match self.args.format {
            Format::Unicode => grid.display_table("i", "j"),
            Format::TeX     => grid.tex_table("$\\mathit{CKhI}$", " ")
        };
        self.out(&table);

        // Generators
        if self.args.show_gens { 
            self.show_gens(&ckhi);
        }

        // Diff
        if self.args.show_diff { 
            let bigraded = h.is_zero() && t.is_zero();
            self.show_diff(&ckhi, bigraded);
        }
    
        // Alpha
        if self.args.show_alpha { 
            self.show_alpha(&ckhi);
        }

        let res = self.flush();
        Ok(res)
    }

    fn show_gens(&mut self, ckh: &KhIComplex<R>) { 
        for i in ckh.support() {
            let c = &ckh[i];
            if c.is_zero() { continue }
            
            self.out(&format!("C[{i}]: {}", c));
    
            let r = c.rank() + c.tors().len();
            for i in 0..r { 
                let z = c.gen(i);
                self.out(&format!("  {i}: {z}"));
            }
            self.out("");
        }
    }

    fn show_diff(&mut self, ckh: &KhIComplex<R>, bigraded: bool) { 
        if bigraded { 
            let ckh = ckh.clone().into_bigraded();
            self.out(&ckh.display_d());
        } else { 
            self.out(&ckh.display_d());
        }
    }

    fn show_alpha(&mut self, ckh: &KhIComplex<R>) { 
        for (i, z) in ckh.canon_cycles().iter().enumerate() { 
            let h = z.h_deg();
            let v = ckh[h].vectorize(z);
            self.out(&format!("a[{i}] in CKhI[{h}]: {}", vec2str(&v)));
            self.out(&format!("  {z}\n"));
        }
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
            c_type: CType::F2,
            ..Default::default()
        };
        let res = dispatch(&args);
        assert!(res.is_ok());
    }

    #[test]
    fn test2() {
        let args = Args {
            link: "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]".to_string(),
            c_value: "1".to_string(),
            c_type: CType::F2,
            mirror: true,
            reduced: true,
            show_alpha: true,
            ..Default::default()
        };
        let res = dispatch(&args);
        assert!(res.is_ok());
    }

    #[cfg(feature = "poly")]
    mod poly_tests {
        use super::*;

        #[test]
        fn test_poly_h() {
            let args = Args {
                link: "3_1".to_string(),
                c_value: "H".to_string(),
                c_type: CType::F2,
                ..Default::default()
            };
            let res = dispatch(&args);
            assert!(res.is_ok());
        }

        #[test]
        fn test_poly_t() {
            let args = Args {
                link: "3_1".to_string(),
                c_value: "0,T".to_string(),
                c_type: CType::F2,
                ..Default::default()
            };
            let res = dispatch(&args);
            assert!(res.is_ok());
        }

        #[test]
        fn test_poly_ht() {
            let args = Args {
                link: "3_1".to_string(),
                c_value: "H,T".to_string(),
                c_type: CType::F2,
                ..Default::default()
            };
            let res = dispatch(&args);
            assert!(res.is_ok());
        }
    }
}
