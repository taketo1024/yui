use std::ops::RangeInclusive;
use clap::ValueEnum;
use derive_more::Display;
use yui_kh::tng::builder::{BuildMode, NodeOrder, ChunkStrategy};
use yui_link::Edge;

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

#[derive(Clone, Copy, PartialEq, Eq, ValueEnum, Display, Debug, Default)]
#[clap(rename_all="verbatim")]
pub enum CType { 
    #[default] Z, 
    Q, F2, F3
}

impl CType { 
    pub fn is_field(&self) -> bool { 
        use CType::*;
        match self {
            Q | F2 | F3 => true,
            _ => false
        }
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
    Ok(parse(lo, isize::MIN + 1)? ..= parse(hi, isize::MAX - 1)?)
}

// parse a build mode: greedy | selective | min-fill | none.
pub fn parse_build_mode(s: &str) -> Result<BuildMode, String> {
    match s.to_lowercase().as_str() {
        "greedy"               => Ok(BuildMode::Greedy),
        "selective"            => Ok(BuildMode::Selective),
        "min-fill" | "minfill" => Ok(BuildMode::MinFill),
        "no-elim" | "noelim"   => Ok(BuildMode::NoElim),
        "none"                 => Ok(BuildMode::None),
        _ => Err(format!("invalid mode `{s}`, expected greedy|selective|min-fill|no-elim|none")),
    }
}

// parse a node order: loop-greedy | min-cut.
pub fn parse_node_order(s: &str) -> Result<NodeOrder, String> {
    match s.to_lowercase().as_str() {
        "loop-greedy" | "loopgreedy" | "loop" => Ok(NodeOrder::LoopGreedy),
        "min-cut" | "mincut" | "min"          => Ok(NodeOrder::MinCut),
        "given"                               => Ok(NodeOrder::Given),
        _ => Err(format!("invalid node order `{s}`, expected loop-greedy|min-cut|given")),
    }
}

// parse a chunk strategy: frontier | boundary. (Manual is selected via `--cut`.)
pub fn parse_chunk_strategy(s: &str) -> Result<ChunkStrategy, String> {
    match s.to_lowercase().as_str() {
        "frontier" => Ok(ChunkStrategy::Frontier),
        "boundary" => Ok(ChunkStrategy::Boundary),
        _ => Err(format!("invalid chunk strategy `{s}`, expected frontier|boundary")),
    }
}

// parse one manual cut `e,e,e` — a comma list of edge labels (the flag is repeatable for more cuts).
pub fn parse_cut(s: &str) -> Result<Vec<Edge>, String> {
    s.split(',').map(|e| e.trim().parse::<Edge>().map_err(|err| format!("`{e}`: {err}"))).collect()
}