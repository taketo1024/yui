use crate::app::args::*;
use crate::app::utils::dispatch::dispatch_field;
use crate::app::utils::*;
use crate::app::err::*;
use std::marker::PhantomData;
use std::str::FromStr;
use yui_core::Field;
use yui_core::FieldOps;
use yui_core::TeX;
use yui_core::{EucRing, EucRingOps};
use yui_homology::DisplaySeq;
use yui_homology::{DisplayTable, GridTrait, SummandTrait};
use yui_kh::kh::ext::cc::KhChainMap;
use yui_kh::kh::KhComplex;
use yui_kh::kh::KhHomology;
use yui_kh::kh::ext::sl2::KhSl2Map;

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    if args.is_field() { 
        dispatch_field!(App, boot_field, args)
    } else {
        dispatch_eucring!(App, boot, args)
    }
}

#[derive(Clone, Default, Debug, clap::Args)]
pub struct Args {
    pub link: String,

    #[arg(short = 't', long, default_value = "Q")]
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

impl AppArgs for Args {
    fn c_type(&self) -> CType { self.c_type }
    fn c_value(&self) -> &String { &self.c_value }
    fn log(&self) -> u8 { self.log }
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
    pub fn boot_field(args: &Args) -> Result<String, Box<dyn std::error::Error>>
    where R: Field, for<'x> &'x R: FieldOps<R> {
        let mut app = Self::new(args.clone());
        app.run_field()
    } 

    pub fn boot(args: &Args) -> Result<String, Box<dyn std::error::Error>> { 
        let mut app = Self::new(args.clone());
        app.run()
    }

    pub fn new(args: Args) -> Self { 
        let buff = String::with_capacity(1024);
        App { args, buff, _ring: PhantomData }
    }

    pub fn run_field(&mut self) -> Result<String, Box<dyn std::error::Error>>
    where R: Field, for<'x> &'x R: FieldOps<R> {
        let (kh, e, ht) = self.compute()?;

        let e_str = e.string_decomp(&kh);
        let e = e.into_chain_map();

        self.show_results(&kh, &e, &ht);
        self.out(&e_str.to_string());
        self.out("");
        
        let res = self.flush();
        Ok(res)
    }

    pub fn run(&mut self) -> Result<String, Box<dyn std::error::Error>> {
        let (kh, e, ht) = self.compute()?;
        let e = e.into_chain_map();

        self.show_results(&kh, &e, &ht);

        let res = self.flush();
        Ok(res)
    }

    pub fn show_results(&mut self, kh: &KhHomology<R>, e: &KhChainMap<R>, ht: &(R, R)) {
        let (h, t) = ht;
        let is_zero = h.is_zero() && t.is_zero();
        let bigraded = is_zero || ["H", "0,T"].contains(&self.args.c_value.as_str());

        self.show_table(kh, bigraded);

        if self.args.show_matrix { 
            if is_zero {
                self.show_matrix_bigr(kh, e, (-2, -4));
            } else { 
                self.show_matrix(kh, e, -2);
            }
        }
    }

    fn compute(&self) -> Result<(KhHomology<R>, KhSl2Map<R>, (R, R)), Box<dyn std::error::Error>> { 
        let (h, t) = parse_pair::<R>(&self.args.c_value)?;

        if self.args.reduced {
            ensure!(t.is_zero(), "`t` must be zero for reduced.");
        }
    
        let r = self.args.reduced;
        let l = load_link(&self.args.link, self.args.mirror)?;
        assert!(l.is_knot());

        let c = KhComplex::new_no_simplify(&l, &h, &t, r);
        let e = c.sl2_map(&l);
        let kh = c.homology();

        Ok((kh, e, (h, t)))
    }

    fn show_table(&mut self, h: &KhHomology<R>, bigraded: bool) { 
        let table = if bigraded { 
            h.display_table()
        } else { 
            h.display_seq()
        };

        self.out(&table);

        if self.args.show_gens { 
            self.show_gens(h);
        }
    }

    fn show_gens(&mut self, h: &KhHomology<R>) { 
        for &i in h.support() {
            if h[i].is_zero() { continue }

            self.out(&format!("({i}): {}", h[i]));

            for (k, z) in h[i].generators().enumerate() { 
                self.out(&format!("  {k}: {z}"));
            }
            self.out("");
        }
    }

    fn show_matrix_bigr(&mut self, h: &KhHomology<R>, f: &KhChainMap<R>, deg: (isize, isize)) { 
        for d in h.delta_range().step_by(2) { 
            self.out(&format!("delta: {d}\n"));

            for i1 in h.h_range().rev() { 
                let j1 = 2 * i1 - d;
                let i2 = i1 + deg.0;
                let j2 = j1 + deg.1;

                let h1 = &h[(i1, j1)];
                let h2 = &h[(i2, j2)];

                if h1.is_zero() || h2.is_zero() { continue; }

                let mat = h1.make_matrix_euc(h2, |z| f.apply(i1, z)).into_dense();
                let r = mat.rank();
                
                self.out(&format!("  ({i1}, {j1}): {} -> ({i2}, {j2}): {}; rank: {}", h[(i1, j1)], h[(i2, j2)], r));
                self.out(&format!("{}\n", mat.to_string().trim_end()));
            }
        }
        self.out("");
    }

    fn show_matrix(&mut self, h: &KhHomology<R>, f: &KhChainMap<R>, deg: isize) { 
        for i1 in h.h_range().rev() { 
            let i2 = i1 + deg;
            let h1 = &h[i1];
            let h2 = &h[i2];

            if h1.is_zero() || h2.is_zero() { continue; }

            let mat = h1.make_matrix_euc(h2, |z| f.apply(i1, z)).into_dense();
            let r = mat.rank();
            
            if self.args.show_matrix { 
                self.out(&format!("  {i1}: {} -> {i2}: {}; rank: {}", h[i1], h[i2], r));
                self.out(&format!("{}\n", mat.to_string().trim_end()));
            }
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