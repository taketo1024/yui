use smart_default::SmartDefault;
use std::ops::RangeInclusive;
use yui_kh::khi::{ssi_invariant_ver, SsiVersion};
use yui_kh::tng::builder::{CutOption, NodeOrder, Strategy, SymBuildConfig};

use crate::app::args::*;
use crate::app::utils::*;

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    App::new(args.clone()).run()
}

/// The equivariant Rasmussen invariant. Always computed over `F2[H]`, so unlike `khi` this
/// command takes no coefficient ring and prints no table.
#[derive(Clone, SmartDefault, PartialEq, Debug, clap::Args)]
pub struct Args {
    pub link: String,

    // fixed: ssi is an F2[H] computation, so there is no ring to choose.
    #[arg(skip = CType::F2)]
    #[default(CType::F2)]
    pub c_type: CType,

    #[arg(skip = String::from("H"))]
    #[default("H".to_string())]
    pub c_value: String,

    #[arg(short, long)]
    pub mirror: bool,

    #[arg(short, long)]
    pub reduced: bool,

    // the guessed s-value seeding the high-q build cut; `None` = full build.
    #[arg(short, long)]
    pub expected: Option<isize>,

    // pipeline: v2 (default, `H = 1` q-truncated solves) or v1 (homology of the cone).
    #[arg(long, value_parser = parse_ssi_version, default_value = "v2")]
    #[default(SsiVersion::V2)]
    pub ver: SsiVersion,

    #[arg(long, value_parser = parse_h_range)]
    pub h_range: Option<RangeInclusive<isize>>,

    #[arg(long, value_parser = parse_strategy, default_value = "greedy")]
    pub strategy: Strategy,

    // crossing order: min-cut (default; bounds cutwidth) or given (PD order, debug).
    #[arg(long, value_parser = parse_node_order, default_value = "min-cut")]
    pub node_order: NodeOrder,

    // skip the half-build/τ-mirror preprocess (which materializes the unbridged off-axis product).
    #[arg(long)]
    pub no_preprocess: bool,

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

pub struct App {
    args: Args,
}

impl App {
    pub fn new(args: Args) -> Self {
        App { args }
    }

    pub fn run(&mut self) -> Result<String, Box<dyn std::error::Error>> {
        let l = load_sinv_knot(&self.args.link, self.args.mirror)?;

        let config = SymBuildConfig {
            h_range: self.args.h_range.clone(), // open ends are clamped inside the build
            strategy: self.args.strategy,
            node_order: self.args.node_order,
            preprocess: !self.args.no_preprocess,
            cut: self.args.cut.clone().unwrap_or_default(),
            max_elim_cost: self.args.max_elim_cost,
            no_full_deloop: self.args.no_full_deloop,
            ..Default::default()
        };

        let (s0, s1) = ssi_invariant_ver(&l, self.args.reduced, config, self.args.expected, self.args.ver);

        Ok(format!("ssi = ({s0}, {s1})"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use crate::app::app::{CliArgs, Cmd};
    use crate::app::cmd::test_utils::{inv_pd, assert_out, assert_cli_default};

    #[test]
    fn cli_defaults() {
        let link = inv_pd("3_1");
        let Cmd::Ssi(a) = CliArgs::parse_from(["ykh", "ssi", &link]).command else {
            panic!("`ssi` routed to the wrong subcommand")
        };
        assert_cli_default(&a, &Args { link, ..Default::default() });
        // the ring is fixed, not selectable
        assert_eq!(a.c_type, CType::F2);
        assert_eq!(a.c_value, "H");
    }

    #[test]
    fn ssi_trefoil() {
        let args = Args { link: inv_pd("3_1"), ..Default::default() };
        assert_out(dispatch(&args), "ssi = (2, 2)");
    }

    #[test]
    fn ssi_trefoil_mirror() {
        // (s_lo, s_hi)(-K) = (-s_hi, -s_lo)
        let args = Args { link: inv_pd("3_1"), mirror: true, ..Default::default() };
        assert_out(dispatch(&args), "ssi = (-2, -2)");
    }

    // 8_21b is one of the knots with s_lo != s_hi ([Sano, InvKh II, Prop. 1.3]).
    #[test]
    fn ssi_8_21b() {
        let args = Args { link: inv_pd("8_21b"), ..Default::default() };
        assert_out(dispatch(&args), "ssi = (2, 4)");
    }

    #[test]
    fn ssi_v1_agrees_with_v2() {
        let v2 = Args { link: inv_pd("8_21b"), ..Default::default() };
        let v1 = Args { ver: SsiVersion::V1, ..v2.clone() };
        assert_eq!(dispatch(&v1).unwrap(), dispatch(&v2).unwrap());
    }
}
