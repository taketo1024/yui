use std::marker::PhantomData;
use std::ops::RangeInclusive;
use std::str::FromStr;
use yui_core::TeX;
use yui_core::{EucRing, EucRingOps};
use yui_homology::{ToSeqString, ToTableString};
use yui_kh::khi::{KhIChain, KhIHomology, ssi_invariants_via_cone};
use yui_kh::tng::builder::{SymBuildConfig, BuildMode, NodeOrder, CutOption};
use yui_link::InvLink;
use crate::app::args::*;
use crate::app::utils::*;
use crate::app::err::*;

pub fn dispatch(args: &Args) -> Result<String, Box<dyn std::error::Error>> {
    dispatch_eucring!(App, boot, args)
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

    #[arg(short = 'a', long)]
    pub show_alpha: bool,

    #[arg(short = 's', long)]
    pub show_ssi: bool,

    #[arg(short = 'n', long)]
    pub no_simplify: bool,

    #[arg(long, value_parser = parse_h_range)]
    pub h_range: Option<RangeInclusive<isize>>,

    #[arg(long, value_parser = parse_build_mode, default_value = "greedy")]
    pub mode: BuildMode,

    // crossing order: min-cut (default; bounds cutwidth) or given (PD order, debug).
    #[arg(long, value_parser = parse_node_order, default_value = "min-cut")]
    pub node_order: NodeOrder,

    // skip the half-build/τ-mirror preprocess (which materializes the unbridged off-axis product).
    #[arg(long)]
    pub no_preprocess: bool,

    // temporary A/B: build KhI as the cobordism-level cone (ConeBuilder) instead of the matrix cone.
    #[arg(long)]
    pub cob_cone: bool,

    // cone only: direct symmetry-broken emission (Sano2026 Prop 4.6).
    #[arg(long)]
    pub cone_direct: bool,

    // cap the per-elimination fill cost; survivors defer to the matrix reduction.
    #[arg(long)]
    pub elim_max_cost: Option<usize>,

    // take the strongly-invertible Whitehead double of the (mirrored) knot: clasp sign pos|neg.
    #[arg(long, value_parser = parse_clasp)]
    pub whitehead: Option<bool>,

    // chunking: `auto(k)` (cutwidth) or manual τ-symmetric edge-cut(s) `e,e,e;e,e,e`.
    #[arg(long, value_parser = parse_cut)]
    pub cut: Option<CutOption>,

    #[arg(short, long, default_value = "unicode")]
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
        let (h, t) = parse_pair::<R>(&self.args.c_value)?;

        ensure!(self.args.c_type == CType::F2, "Only `-t F2` is supported.");

        if self.args.reduced { 
            ensure!(t.is_zero(), "`t` must be zero for reduced.");
        }
        if self.args.show_alpha { 
            ensure!(t.is_zero(), "`t` must be zero to have alpha.");
        }
        if self.args.show_ssi {
            ensure!(!h.is_zero() && !h.is_unit(), "`h` must be non-zero, non-invertible to compute ssi.");
            ensure!(t.is_zero(), "`t` must be zero to compute ss.");
        }

        let l = load_sinv_knot(&self.args.link, self.args.mirror)?;
        let l = if let Some(positive) = self.args.whitehead {
            l.whitehead_double(positive, 0)
        } else {
            l
        };

        let config = SymBuildConfig {
            h_range: self.args.h_range.clone(), // open ends are clamped inside the build
            mode: self.args.mode,
            node_order: self.args.node_order,
            preprocess: !self.args.no_preprocess,
            cut: self.args.cut.clone().unwrap_or_default(),
            cone_direct: self.args.cone_direct,
            elim_max_cost: self.args.elim_max_cost,
            ..Default::default()
        };

        // ssi-only: canon classes ride the trans-free reduction as vectors (no table output).
        let ssi_only = self.args.show_ssi && !(self.args.show_gens || self.args.show_alpha);
        if ssi_only && self.args.cob_cone && !self.args.no_simplify {
            let ssi = ssi_invariants_via_cone(&l, &h, self.args.reduced, config);
            self.out(&format!("ssi = ({}, {})", ssi.0, ssi.1));
            return Ok(self.flush());
        }

        let khi = if self.args.no_simplify {
            KhIHomology::new_no_simplify(&l, &h, &t, self.args.reduced)
        } else {
            if self.args.cob_cone {
                KhIHomology::from_cone(&l, &h, &t, self.args.reduced, config)
            } else {
                KhIHomology::new_with_config(&l, &h, &t, self.args.reduced, config)
            }
        };

        let bigraded = h.is_zero() && t.is_zero() || 
            ["H", "0,T"].contains(&self.args.c_value.as_str());

        let table = if bigraded { 
            khi.to_table_string()
        } else { 
            khi.to_seq_string()
        };
        self.out(&table);

        if self.args.show_gens { 
            self.show_gens(&khi);
        }

        if self.args.show_alpha { 
            let zs = khi.canon_cycles();
            self.show_alpha(&khi, zs);
        }

        if self.args.show_ssi { 
            let zs = khi.canon_cycles();
            self.show_ssi(&l, &h, &khi, zs)?;
        }

        Ok(self.flush())
    }

    fn show_gens(&mut self, khi: &KhIHomology<R>) { 
        for &i in khi.support() {
            let h = &khi[i];
            if h.is_zero() { continue }

            self.out(&format!("KhI[{i}]: {}", h));

            let r = h.rank() + h.tors().len();
            for i in 0..r { 
                let z = h.generator(i);
                self.out(&format!("  {i}: {z}"));
            }
            self.out("");
        }
    }

    fn show_alpha(&mut self, khi: &KhIHomology<R>, zs: &[KhIChain<R>]) {
        for (i, z) in zs.iter().enumerate() {
            let h = khi.h_deg_of_chain(z);
            let v = khi[h].vectorize_euc(z);
            self.out(&format!("a[{i}] in KhI[{h}]: {}", vec2str(&v)));
            self.out(&format!("  {z}\n"));
        }
    }

    fn show_ssi(&mut self, l: &InvLink, c: &R, khi: &KhIHomology<R>, zs: &[KhIChain<R>]) -> Result<(), Box<dyn std::error::Error>> { 
        assert!(!c.is_unit() && !c.is_unit());

        use yui_kh::util::calc::div_vec;

        let l = l.inner();
        let w = l.writhe();
        let r = l.seifert_circles().len() as i32;

        for (i, z) in zs.iter().enumerate() {
            let h = &khi[khi.h_deg_of_chain(z)];
            let v = h.vectorize(z).subvec(0..h.rank());
            let d = div_vec(&v, c);

            ensure!(d.is_some(), "invalid div for a = {z}.");

            let d = d.unwrap();
            let s = 2 * d + w - r + 1;

            self.out(&format!("ss[{i}] = {s} (d = {d}, w = {w}, r = {r})"));
        }

        Ok(())
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
            c_type: CType::F2,
            c_value: "0".to_string(),
            ..Default::default()
        };
        let res = dispatch(&args);
        assert!(res.is_ok());
    }

    #[test]
    fn test2() { 
        let args = Args { 
            link: "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]".to_string(),
            c_type: CType::F2,
            c_value: "1".to_string(),
            mirror: true,
            reduced: true,
            ..Default::default()
        };
        let res = dispatch(&args);
        assert!(res.is_ok());
    }

    #[test]
    fn test_poly_h() { 
        let args = Args {
            link: "3_1".to_string(),
            c_type: CType::F2,
            c_value: "H".to_string(),
            ..Default::default()
        };
        let res = dispatch(&args);
        assert!(res.is_ok());
    }
}