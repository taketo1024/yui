//! `ss`: the Rasmussen `s` invariant of a knot, over a chosen ring.

use smart_default::SmartDefault;
use std::marker::PhantomData;
use std::ops::RangeInclusive;
use std::str::FromStr;
use yui_core::abst::{EucRing, EucRingOps, Field, FieldOps};
use yui_kh::ss::{s_invariant_with, ss_invariant_with, SsVersion};
use yui_kh::tng::builder::{BuildConfig, CutOption, NodeOrder, Strategy};

use crate::app::args::*;
use crate::app::utils::*;
use crate::app::utils::dispatch::dispatch_field;
use crate::app::err::*;

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    // over `F[H]` the invariant is the Rasmussen `s`, which has a specialized pipeline.
    if args.c_type.is_field() && args.c_value == "H" {
        dispatch_field!(SApp, boot, args)
    } else {
        dispatch_eucring!(App, boot, args)
    }
}

/// The slice-torus invariant `ss̃_c`. Over `F[H]` (a field `F`, `c = H`) this is the Rasmussen
/// invariant `s` over `F`, and is computed by the `H = 1` pipeline.
#[derive(Clone, SmartDefault, PartialEq, Debug, clap::Args)]
pub struct Args {
    pub link: String,

    #[arg(short = 't', long, default_value = "Q")]
    #[default(CType::Q)]
    pub c_type: CType,

    #[arg(short, long, default_value = "H")]
    #[default(String::from("H"))]
    pub c_value: String,

    #[arg(short, long)]
    pub mirror: bool,

    #[arg(short, long)]
    pub reduced: bool,

    // pipeline: v2 (default, `H = 1` q-truncated solves) or v1 (divisibility in the homology).
    // Only applies over `F[H]`; the general `(R, c)` route is always v1.
    #[arg(long, value_parser = parse_ss_version, default_value = "v2")]
    #[default(SsVersion::V2)]
    pub ver: SsVersion,

    #[arg(long, value_parser = parse_h_range, allow_hyphen_values = true)]
    pub h_range: Option<RangeInclusive<isize>>,

    #[arg(long, value_parser = parse_strategy, default_value = "greedy")]
    pub strategy: Strategy,

    // crossing order: min-cut (default; bounds cutwidth) or given (PD order, debug).
    #[arg(long, value_parser = parse_node_order, default_value = "min-cut")]
    pub node_order: NodeOrder,

    // cap the per-elimination fill cost; survivors defer to the matrix reduction.
    #[arg(long)]
    pub max_elim_cost: Option<usize>,

    // skip the final deloop/eliminate; remaining circles defer to the matrix reducer.
    #[arg(long)]
    pub no_full_deloop: bool,

    // chunking: `N` (cutwidth, N pieces) or `at(c,..)` (cut after the given crossing counts).
    #[arg(long, value_parser = parse_cut)]
    pub cut: Option<CutOption>,

    #[arg(long, default_value = "0")]
    pub log: u8,
}

impl AppArgs for Args {
    fn c_type(&self) -> CType { self.c_type }
    fn c_value(&self) -> &String { &self.c_value }
    fn log(&self) -> u8 { self.log }
}

impl Args {
    fn build_config(&self) -> BuildConfig {
        BuildConfig {
            h_range: self.h_range.clone(),
            strategy: self.strategy,
            node_order: self.node_order,
            cut: self.cut.clone().unwrap_or_default(),
            max_elim_cost: self.max_elim_cost,
            no_full_deloop: self.no_full_deloop,
            ..Default::default()
        }
    }
}

/// `s` over a field `F`, i.e. `ss̃_H` over `F[H]`.
pub struct SApp<F> {
    _field: PhantomData<F>,
}

impl<F> SApp<F>
where F: Field + FromStr + FieldOps<F>, for<'x> &'x F: FieldOps<F> {
    pub fn boot(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
        let l = load_link(&args.link, args.mirror)?;
        let s = s_invariant_with::<F>(&l, args.reduced, args.build_config(), args.ver);
        Ok(format!("s = {s}"))
    }
}

/// `ss̃_c` for a general Euclidean ring `R` and non-unit `c`.
pub struct App<R> {
    _ring: PhantomData<R>,
}

impl<R> App<R>
where R: EucRing + FromStr, for<'x> &'x R: EucRingOps<R> {
    pub fn boot(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
        let Ok(c) = R::from_str(&args.c_value) else {
            return err!("cannot parse '{}' as {}.", args.c_value, std::any::type_name::<R>())
        };
        ensure!(!c.is_zero() && !c.is_unit(), "`c` must be non-zero, non-invertible.");

        let l = load_link(&args.link, args.mirror)?;
        let ss = ss_invariant_with(&l, &c, args.reduced, args.build_config());
        Ok(format!("ss = {ss}"))
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
        let Cmd::Ss(a) = CliArgs::parse_from(["ykh", "ss", &link]).command else {
            panic!("`ss` routed to the wrong subcommand")
        };
        assert_cli_default(&a, &Args { link, ..Default::default() });
        // the default is `s` over Q, not the general `ss̃_c`
        assert_eq!(a.c_type, CType::Q);
        assert_eq!(a.c_value, "H");
    }

    #[test]
    fn s_over_q_by_default() {
        let args = Args { link: pd("3_1"), ..Default::default() };
        assert_out(dispatch(&args), "s = 2");
    }

    #[test]
    fn s_of_the_mirror() {
        let args = Args { link: pd("3_1"), mirror: true, ..Default::default() };
        assert_out(dispatch(&args), "s = -2");
    }

    #[test]
    fn ss_falls_back_for_a_general_ring() {
        // c = 2 over Z is not of the form F[H], so this takes the `ss̃_c` route.
        let args = Args {
            link: pd("3_1"),
            c_type: CType::Z,
            c_value: "2".to_string(),
            ..Default::default()
        };
        assert_out(dispatch(&args), "ss = 2");
    }

    #[test]
    fn rejects_an_invertible_c() {
        let args = Args {
            link: pd("3_1"),
            c_type: CType::Z,
            c_value: "1".to_string(),
            ..Default::default()
        };
        assert!(dispatch(&args).is_err());
    }

    // 14n_19265: `s` over F2 differs from `s` over Q, and `ss̃_c` differs for c = 2 vs 3.
    #[test]
    fn ring_dependence() {
        let link = pd("14n_19265");
        let s = |c_type, c_value: &str| {
            let args = Args { link: link.clone(), c_type, c_value: c_value.to_string(), ..Default::default() };
            dispatch(&args).unwrap()
        };

        assert_eq!(s(CType::Q,  "H"), "s = 0");
        assert_eq!(s(CType::F2, "H"), "s = -2");
        assert_eq!(s(CType::F3, "H"), "s = 0");
        assert_eq!(s(CType::Z,  "2"), "ss = -2");
        assert_eq!(s(CType::Z,  "3"), "ss = 0");
    }
}
