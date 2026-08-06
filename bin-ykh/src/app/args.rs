//! Argument types shared by the subcommands — the coefficient ring, output
//! format, build strategy and chunking — with their `clap` value parsers.

use std::ops::RangeInclusive;
use clap::ValueEnum;
use derive_more::Display;
use yui_link::Link;
use yui_kh::ss::SsVersion;
use yui_kh::tng::builder::{Strategy, NodeOrder, CutOption};

pub trait AppArgs {
    fn c_type(&self) -> CType;
    fn c_value(&self) -> &String;
    fn log(&self) -> u8;

    fn poly_vars(&self) -> PolyVars {
        parse_poly_vars(self.c_value())
    }

    fn is_poly(&self) -> bool {
        self.poly_vars() != PolyVars::None
    }

    fn is_field(&self) -> bool {
        self.c_type().is_field() && !self.is_poly()
    }

    fn is_euc_ring(&self) -> bool {
        self.c_type() == CType::Z && !self.is_poly() ||
        self.c_type().is_field() && self.poly_vars().nvars() == 1
    }

    fn log_level(&self) -> log::LevelFilter {
        use log::LevelFilter::*;
        match self.log() {
            1 => Info,
            2 => Debug,
            3 => Trace,
            _ => Off,
        }
    }
}

// no `Default`: each command's `Args` states its own, matching its clap `default_value`.
#[derive(Clone, Copy, PartialEq, Eq, ValueEnum, Display, Debug)]
#[clap(rename_all="verbatim")]
pub enum CType {
    Z, Q, F2, F3
}

impl CType {
    pub fn is_field(&self) -> bool {
        use CType::*;
        matches!(self, Q | F2 | F3)
    }
}

#[derive(PartialEq, Eq)]
pub(crate) enum PolyVars {
    H, T, HT, None
}

impl PolyVars {
    pub fn nvars(&self) -> usize {
        match self {
            PolyVars::H | PolyVars::T  => 1,
            PolyVars::HT => 2,
            _ => 0,
        }
    }
}

pub(crate) fn parse_poly_vars(c_value: &String) -> PolyVars {
    use std::collections::HashSet;

    let s: HashSet<_> = c_value.split(',').collect();
    match (s.contains("H"), s.contains("T")) {
        (true,  true)  => PolyVars::HT,
        (true,  false) => PolyVars::H,
        (false, true)  => PolyVars::T,
        (false, false) => PolyVars::None
    }
}

#[derive(Clone, Copy, PartialEq, Eq, ValueEnum, Display, Debug, Default)]
#[clap(rename_all="lower")]
pub enum Format {
    #[default] Unicode,
    TeX
}

// Parse an inclusive degree range like `-1..=2` (also accepts `-1..2`). An open bottom `..=2`
// parses with a near-min sentinel start, an open top `-1..=` with a near-max end; the command
// clamps these to the complex's degree span (so the bottom resolves to deg_shift.0). The `±1`
// margin keeps a stray `start-1` / `end+1` widening from overflowing before the clamp.
pub fn parse_h_range(s: &str) -> Result<RangeInclusive<isize>, String> {
    let (lo, hi) = s.split_once("..=")
        .or_else(|| s.split_once(".."))
        .ok_or_else(|| format!("invalid range `{s}`, expected e.g. `-1..=2` or `..=2`"))?;
    let parse = |x: &str, open: isize| {
        let x = x.trim();
        if x.is_empty() { Ok(open) } else { x.parse::<isize>().map_err(|e| format!("`{x}`: {e}")) }
    };
    Ok(parse(lo, -(Link::MAX_CROSSING as isize))? ..= parse(hi, Link::MAX_CROSSING as isize + 2)?)
}

// parse a build strategy: greedy | min-fill | no-elim | none.
pub fn parse_strategy(s: &str) -> Result<Strategy, String> {
    match s.to_lowercase().as_str() {
        "greedy"               => Ok(Strategy::Greedy),
        "min-fill" | "minfill" => Ok(Strategy::MinFill),
        "no-elim" | "noelim"   => Ok(Strategy::NoElim),
        "none"                 => Ok(Strategy::None),
        _ => Err(format!("invalid strategy `{s}`, expected greedy|min-fill|no-elim|none")),
    }
}

// parse the ss/ssi pipeline version: v1 | v2.
pub fn parse_ss_version(s: &str) -> Result<SsVersion, String> {
    match s.to_lowercase().as_str() {
        "v1" | "1" => Ok(SsVersion::V1),
        "v2" | "2" => Ok(SsVersion::V2),
        _ => Err(format!("invalid version `{s}`, expected v1|v2")),
    }
}

// parse a node order: min-cut | given.
pub fn parse_node_order(s: &str) -> Result<NodeOrder, String> {
    match s.to_lowercase().as_str() {
        "min-cut" | "mincut" | "min" => Ok(NodeOrder::MinCut),
        "given"                      => Ok(NodeOrder::Given),
        _ => Err(format!("invalid node order `{s}`, expected min-cut|given")),
    }
}

// parse `--cut`: `N` for cutwidth chunking into `N` pieces, `at(c,..)` to cut after the given
// cumulative crossing counts.
pub fn parse_cut(s: &str) -> Result<CutOption, String> {
    if let Some(inner) = s.trim().strip_prefix("at(").and_then(|x| x.strip_suffix(')')) {
        let counts = inner.split(',')
            .map(|c| c.trim().parse::<usize>().map_err(|e| format!("at(c,..): `{c}`: {e}")))
            .collect::<Result<Vec<usize>, _>>()?;
        return Ok(CutOption::At(counts));
    }
    if let Ok(k) = s.trim().parse::<usize>() {
        return Ok(CutOption::Auto(k));
    }
    Err(format!("invalid cut `{s}`, expected N|at(c,..)"))
}