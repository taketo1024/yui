use crate::app::utils::*;
use crate::app::err::*;
use std::marker::PhantomData;
use std::str::FromStr;
use yui_core::tex::TeX;
use yui_core::{EucRing, EucRingOps};
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
    
        let r = self.args.reduced;
        let bigraded = (h.is_zero() && t.is_zero()) || 
            ["H", "0,T"].contains(&self.args.c_value.as_str());
        
        let l = load_link(&self.args.link, self.args.mirror)?;
        assert!(l.is_knot());

        let c = KhComplex::new_no_simplify(&l, &h, &t, r);
        let e = c.e_map(&l);
        let h = c.homology();

        self.show_table(&h, bigraded);
        self.show_map("e", &h, &e);
        
        let res = self.flush();
        Ok(res)
    }

    fn show_table(&mut self, h: &KhHomology<R>, bigraded: bool) { 
        let table = if bigraded { 
            h.gen_grid().display_table("i", "j")
        } else { 
            h.display_seq("i")
        };

        self.out(&table);

        if self.args.show_gens { 
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

    fn show_map(&mut self, f_name: &str, h: &KhHomology<R>, f: &KhChainMap<R>) { 
        self.out(&format!("{f_name}: deg {}\n", f.deg()));

        for i in h.h_range().rev() { 
            let j = i + f.deg();
            if h[i].is_zero() || h[j].is_zero() { continue; }

            self.out(&format!("({i}) {} -> ({j}) {}", h[i], h[j]));

            for z in h[i].gens() { 
                let w = f.apply(i, &z);
                let x = h[i].vectorize_euc(&z).into_vec();
                let y = h[j].vectorize_euc(&w).into_vec();
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