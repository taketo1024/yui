//! `ckh`: print the Khovanov chain complex of a link.

use smart_default::SmartDefault;
use crate::app::args::*;
use crate::app::utils::*;
use crate::app::err::*;
use std::marker::PhantomData;
use std::ops::RangeInclusive;
use std::str::FromStr;
use yui_core::util::tex::TeX;
use yui_homology::tex::ToTexTable;
use yui_core::abst::{Ring, RingOps};
use yui_homology::ToTableString;
use yui_kh::kh::KhComplex;
use yui_kh::tng::builder::{BuildConfig, CutOption, NodeOrder, Strategy};

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    dispatch_ring!(App, boot, args)
}

#[derive(Clone, SmartDefault, PartialEq, Debug, clap::Args)]
pub struct Args {
    pub link: String,

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

    #[arg(short = 'g', long)]
    pub show_gens: bool,

    #[arg(short = 'd', long)]
    pub show_diff: bool,

    #[arg(short = 'a', long)]
    pub show_alpha: bool,

    #[arg(short = 'n', long)]
    pub no_simplify: bool,

    #[arg(long, value_parser = parse_h_range, allow_hyphen_values = true)]
    pub h_range: Option<RangeInclusive<isize>>,

    // chunking: `N` (auto, by cutwidth) or `at(c,..)` (manual crossing counts).
    #[arg(long, value_parser = parse_cut)]
    pub cut: Option<CutOption>,

    // cap the per-elimination fill cost; survivors defer to the matrix reduction.
    #[arg(long)]
    pub max_elim_cost: Option<usize>,

    #[arg(long, value_parser = parse_strategy, default_value = "greedy")]
    pub strategy: Strategy,

    // crossing order: min-cut (default; bounds cutwidth) or given (PD order, debug).
    #[arg(long, value_parser = parse_node_order, default_value = "min-cut")]
    pub node_order: NodeOrder,

    // skip the final deloop/eliminate; remaining circles defer to the matrix reducer.
    #[arg(long)]
    pub no_full_deloop: bool,

    #[arg(short, long, default_value = "unicode")]
    #[default(Format::Unicode)]
    pub format: Format,

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

        if self.args.reduced {
            ensure!(t.is_zero(), "`t` must be zero for reduced.");
        }
        if self.args.show_alpha {
            ensure!(t.is_zero(), "`t` must be zero to have alpha.");
        }

        let l = load_link(&self.args.link, self.args.mirror)?;

        let ckh = if self.args.no_simplify {
            KhComplex::new_no_simplify(&l, &h, &t, self.args.reduced)
        } else {
            let config = BuildConfig {
                strategy: self.args.strategy,
                node_order: self.args.node_order,
                cut: self.args.cut.clone().unwrap_or_default(),
                h_range: self.args.h_range.clone(),
                max_elim_cost: self.args.max_elim_cost,
                no_full_deloop: self.args.no_full_deloop,
                ..Default::default()
            };
            KhComplex::new_with_config(&l, &h, &t, self.args.reduced, config)
        };

        // CKh generators
        let table = match self.args.format {
            Format::TeX => ckh.tex_table("CKh"),
            _           => ckh.to_table_string(),
        };
        self.out(&table);

        // Generators
        if self.args.show_gens {
            self.show_gens(&ckh);
        }

        // Diff
        if self.args.show_diff {
            self.show_diff(&ckh);
        }

        // Alpha
        if self.args.show_alpha {
            self.show_alpha(&ckh);
        }

        let res = self.flush();
        Ok(res)
    }

    fn show_gens(&mut self, ckh: &KhComplex<R>) {
        for &i in ckh.support() {
            let c = &ckh[i];
            if c.is_zero() { continue }

            self.out(&format!("C[{i}]: {}", c));

            let r = c.n_generators();
            for i in 0..r {
                let z = c.generator(i);
                self.out(&format!("  {i}: {z}"));
            }
            self.out("");
        }
    }

    fn show_diff(&mut self, ckh: &KhComplex<R>) {
        self.out(&ckh.describe_d());
    }

    fn show_alpha(&mut self, ckh: &KhComplex<R>) {
        for (i, z) in ckh.canon_cycles().iter().enumerate() {
            let h = ckh.h_deg_of_chain(z);
            let v = ckh[h].vectorize(z);
            self.out(&format!("a[{i}] in CKh[{h}]: {}", vec2str(&v)));
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
    use clap::Parser;
    use crate::app::app::{CliArgs, Cmd};
    use crate::app::cmd::test_utils::{pd, assert_out, assert_cli_default};

    #[test]
    fn cli_defaults() {
        let link = pd("3_1");
        let Cmd::CKh(a) = CliArgs::parse_from(["ykh", "ckh", &link]).command else {
            panic!("`ckh` routed to the wrong subcommand")
        };
        assert_cli_default(&a, &Args { link, ..Default::default() });
    }

    #[test]
    fn ckh_trefoil_z() {
        let args = Args {
            link: pd("3_1"),
            ..Default::default()
        };
        assert_out(dispatch(&args), r"
             j\i  0  1  2  3
             9    .  .  .  Z
             7    .  .  Z  Z
             5    .  .  Z  .
             3    Z  .  .  .
             1    Z  .  .  .
        ");
    }

    #[test]
    fn ckh_trefoil_mirror_reduced_alpha() {
        let args = Args {
            link: pd("3_1"),
            c_value: "2".to_string(),
            mirror: true,
            reduced: true,
            show_alpha: true,
            ..Default::default()
        };
        assert_out(dispatch(&args), r"
             j\i  -3  -2  -1  0
             -2   .   .   .   Z
             -4   .   .   .   .
             -6   .   Z   .   .
             -8   Z   .   .   .

             a[0] in CKh[0]: (-2)
               -2(1X)₁₁₁
        ");
    }

    #[test]
    fn ckh_trefoil_zpoly_h() {
        let args = Args {
            link: pd("3_1"),
            c_value: "H".to_string(),
            ..Default::default()
        };
        assert_out(dispatch(&args), r"
             j\i  0     1  2     3
             9    .     .  .     Z[H]
             7    .     .  Z[H]  Z[H]
             5    .     .  Z[H]  .
             3    Z[H]  .  .     .
             1    Z[H]  .  .     .
        ");
    }

    #[test]
    fn ckh_trefoil_zpoly_t() {
        let args = Args {
            link: pd("3_1"),
            c_value: "0,T".to_string(),
            ..Default::default()
        };
        assert_out(dispatch(&args), r"
             j\i  0     1  2     3
             9    .     .  .     Z[T]
             7    .     .  Z[T]  Z[T]
             5    .     .  Z[T]  .
             3    Z[T]  .  .     .
             1    Z[T]  .  .     .
        ");
    }

    #[test]
    fn ckh_trefoil_zpoly_ht() {
        let args = Args {
            link: pd("3_1"),
            c_value: "H,T".to_string(),
            ..Default::default()
        };
        assert_out(dispatch(&args), r"
             j\i  0        1  2        3
             9    .        .  .        Z[H, T]
             7    .        .  Z[H, T]  Z[H, T]
             5    .        .  Z[H, T]  .
             3    Z[H, T]  .  .        .
             1    Z[H, T]  .  .        .
        ");
    }
}
