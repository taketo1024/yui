use crate::app::utils::*;
use crate::app::err::*;
use std::collections::HashMap;
use std::marker::PhantomData;
use std::str::FromStr;
use itertools::Itertools;
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

    #[arg(short = 'M', long)]
    pub show_matrix: bool,

    #[arg(long)]
    pub verify: bool,

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

        if self.args.verify {
            e.check_all(c.inner(), c.inner());
        }

        let h = c.homology();

        self.show_table(&h, bigraded);
        self.show_map("e", &h, &e, (-2, -4));
        
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

    fn show_map(&mut self, f_name: &str, h: &KhHomology<R>, f: &KhChainMap<R>, deg: (isize, isize)) { 
        self.out(&format!("{f_name}: deg {deg:?}\n"));

        let mut ranks: HashMap<isize, HashMap<isize, usize>> = HashMap::new();

        let grid = h.gen_grid();
        for d in h.delta_range().step_by(2) { 
            if self.args.show_matrix { 
                self.out(&format!("delta: {d}"));
            }

            let mut ranks_d = HashMap::new();

            for i1 in h.h_range() { 
                let j1 = 2 * i1 - d;
                let i2 = i1 + deg.0;
                let j2 = j1 + deg.1;

                if grid[(i1, j1)].is_zero() || grid[(i2, j2)].is_zero() { continue; }

                let mat = f.as_matrix(i1, &grid[(i1, j1)], &grid[(i2, j2)]).into_dense();
                let r = mat.rank();
                
                if self.args.show_matrix { 
                    self.out(&format!("  ({i1}, {j1}) {} -> ({i2}, {j2}) {}; rank: {}", grid[(i1, j1)], grid[(i2, j2)], r));
                    self.out(&format!("{}", mat));
                }

                ranks_d.insert(i1, r);
            }

            ranks.insert(d, ranks_d);
        }

        if self.args.show_matrix { 
            self.out("");
        }

        for d in h.delta_range().step_by(2) { 
            let ranks_d = &ranks[&d];

            if ranks_d.is_empty() { continue; }

            let i0 = ranks_d.keys().min().unwrap();
            let rank_str = ranks_d.iter().sorted_by_key(|v| v.0).map(|(_, r)| r).join(", ");            

            self.out(&format!("rank({d}):\t[{i0}; {rank_str}]"));
        }
        self.out("");
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