//! `cc`: the map on Khovanov homology induced by a crossing change.

use smart_default::SmartDefault;
use crate::app::args::*;
use crate::app::utils::*;
use crate::app::err::*;
use std::marker::PhantomData;
use std::str::FromStr;
use yui_core::util::tex::TeX;
use yui_core::abst::{EucRing, EucRingOps};
use yui_homology::{ToSeqString, ToTableString};
use yui_kh::kh::ext::cc::KhChainMap;
use yui_kh::kh::KhComplex;
use yui_kh::kh::KhHomology;

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    dispatch_eucring!(App, boot, args)
}

#[derive(Clone, SmartDefault, PartialEq, Debug, clap::Args)]
pub struct Args {
    pub link: String,

    #[arg(short = 'i', long, required = true)]
    pub cc_index: usize,

    #[arg(short = 'f', long, required = true)]
    pub map_type: usize,

    #[arg(short = 't', long, default_value = "Z")]
    #[default(CType::Z)]
    pub c_type: CType,

    #[arg(short, long, default_value = "0")]
    #[default("0".to_string())]
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
    pub fn boot(args: &Args) -> Result<String, Box<dyn std::error::Error>> { 
        let mut app = Self::new(args.clone());
        app.run()
    }

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
            l.cc_at(i)
        } else { 
            l
        };

        let (c1, c2) = KhComplex::cc_pair(&l, &h, &t, r, i);

        // MEMO: f0, f1 are named after the paper. 
        let (f_name, f) = if self.args.map_type == 0 { 
            ("f0", KhComplex::cc_map1(&c1, &c2, &l, i))
        } else { 
            ("f1", KhComplex::cc_map0(&c1, &c2, i)) 
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
            h.to_table_string()
        } else { 
            h.to_seq_string()
        };

        self.out(label);
        self.out(&table);

        if with_gens { 
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

    fn show_map(&mut self, f_name: &str, h1: &KhHomology<R>, h2: &KhHomology<R>, f: &KhChainMap<R>) { 
        self.out(&format!("{f_name}: deg {}\n", f.deg()));

        for i in h1.h_range() { 
            let j = i + f.deg();
            let (s1, s2) = (&h1[i], &h2[j]);

            if s1.is_zero() || s2.is_zero() { 
                self.out(&format!("({i}) {s1} -> ({j}) {s2}\n"));
                continue;
            }

            let mat = s1.make_matrix_euc(s2, |z| f.apply(i, z)).into_dense();
            self.out(&format!("({i}) {s1} -> ({j}) {s2}; rank: {}", mat.rank()));
            self.out(&format!("{}\n", mat.to_string().trim_end()));
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
    use clap::Parser;
    use crate::app::app::{CliArgs, Cmd};
    use crate::app::cmd::test_utils::{pd, assert_out, assert_cli_default};

    #[test]
    fn cli_defaults() {
        let link = pd("3_1");
        let Cmd::CC(a) = CliArgs::parse_from(["ykh", "cc", &link, "-i", "0", "-f", "0"]).command else {
            panic!("`cc` routed to the wrong subcommand")
        };
        assert_cli_default(&a, &Args { link, cc_index: 0, map_type: 0, ..Default::default() });
    }

    // f0 lowers the h-degree by 2: the negative-to-positive crossing change on the trefoil.
    #[test]
    fn cc_trefoil_map0() {
        let args = Args {
            link: pd("3_1"),
            cc_index: 0,
            map_type: 0,
            ..Default::default()
        };
        assert_out(dispatch(&args), r"
            from:
             j\i  0  1  2  3
             9    .  .  .  Z
             7    .  .  .  (Z/2)
             5    .  .  Z  .
             3    Z  .  .  .
             1    Z  .  .  .

            to:
             j\i  0
             1    Z
             -1   Z

            f0: deg -2

            (0) Z² -> (-2) 0

            (1) 0 -> (-1) 0

            (2) Z -> (0) Z²; rank: 1

              ┌    ┐
              │  0 │
              │ -1 │
              └    ┘

            (3) Z ⊕ (Z/2) -> (1) 0
        ");
    }

    // f1 preserves the h-degree.
    #[test]
    fn cc_trefoil_map1() {
        let args = Args {
            link: pd("3_1"),
            cc_index: 0,
            map_type: 1,
            ..Default::default()
        };
        assert_out(dispatch(&args), r"
            from:
             j\i  0  1  2  3
             9    .  .  .  Z
             7    .  .  .  (Z/2)
             5    .  .  Z  .
             3    Z  .  .  .
             1    Z  .  .  .

            to:
             j\i  0
             1    Z
             -1   Z

            f1: deg 0

            (0) Z² -> (0) Z²; rank: 2

              ┌       ┐
              │  1  0 │
              │  0 -1 │
              └       ┘

            (1) 0 -> (1) 0

            (2) Z -> (2) 0

            (3) Z ⊕ (Z/2) -> (3) 0
        ");
    }
}
