use crate::app::args::*;
use crate::app::utils::*;
use crate::app::err::*;
use std::marker::PhantomData;
use std::ops::RangeInclusive;
use std::str::FromStr;
use yui_core::TeX;
use yui_core::{Ring, RingOps};
use yui_homology::ToTableString;
use yui_kh::khi::KhIComplex;
use yui_kh::tng::builder::{SymBuildConfig, BuildMode, NodeOrder};

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    dispatch_ring!(App, boot, args)
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

    #[arg(short = 'n', long)]
    pub no_simplify: bool,

    #[arg(long, value_parser = parse_h_range)]
    pub h_range: Option<RangeInclusive<isize>>,

    #[arg(long)]
    pub chunk: Option<usize>,

    #[arg(long, value_parser = parse_build_mode, default_value = "greedy")]
    pub mode: BuildMode,

    // crossing order: min-cut (default; bounds cutwidth) or given (PD order, debug).
    #[arg(long, value_parser = parse_node_order, default_value = "min-cut")]
    pub node_order: NodeOrder,

    // skip the half-build/τ-mirror preprocess (which materializes the unbridged off-axis product).
    #[arg(long)]
    pub no_preprocess: bool,

    #[arg(long, default_value = "0")]
    pub log: u8,
}

impl AppArgs for Args {
    fn c_type(&self) -> CType { self.c_type }
    fn c_value(&self) -> &String { &self.c_value }
    fn log(&self) -> u8 { self.log }
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
    pub fn boot(args: &Args) -> Result<String, Box<dyn std::error::Error>> { 
        let mut app = Self::new(args.clone());
        app.run()
    }

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

        let ckhi = if self.args.no_simplify {
            KhIComplex::new_no_simplify(&l, &h, &t, self.args.reduced)
        } else {
            let config = SymBuildConfig {
                h_range: self.args.h_range.clone(), // open ends are clamped inside the build
                chunk_bound: self.args.chunk,
                mode: self.args.mode,
                node_order: self.args.node_order,
                preprocess: !self.args.no_preprocess,
                ..Default::default()
            };
            KhIComplex::new_with_config(&l, &h, &t, self.args.reduced, config)
        };
        
        // CKh generators
        let table = ckhi.to_table_string();
        self.out(&table);

        // Generators
        if self.args.show_gens { 
            self.show_gens(&ckhi);
        }

        // Diff
        if self.args.show_diff { 
            self.show_diff(&ckhi);
        }
    
        // Alpha
        if self.args.show_alpha { 
            self.show_alpha(&ckhi);
        }

        let res = self.flush();
        Ok(res)
    }

    fn show_gens(&mut self, ckh: &KhIComplex<R>) { 
        for &i in ckh.support() {
            let c = &ckh[i];
            if c.is_zero() { continue }
            
            self.out(&format!("C[{i}]: {}", c));
    
            let r = c.rank() + c.tors().len();
            for i in 0..r { 
                let z = c.generator(i);
                self.out(&format!("  {i}: {z}"));
            }
            self.out("");
        }
    }

    fn show_diff(&mut self, ckh: &KhIComplex<R>) { 
        self.out(&ckh.describe_d());
    }

    fn show_alpha(&mut self, ckh: &KhIComplex<R>) {
        for (i, z) in ckh.canon_cycles().iter().enumerate() {
            let h = ckh.h_deg_of_chain(z);
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
        res.trim_end().to_string()
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
