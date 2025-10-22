#![allow(non_snake_case, unused)]
use itertools::Itertools;
use num_traits::Zero;

use log::*;
use yui::poly::HPoly;
use yui::num::FF2;
use yui::{EucRing, EucRingOps};
use yui_homology::{ChainComplexTrait, DisplaySeq, DisplayTable, SummandTrait};
use yui_kh::kh::KhChainExt;
use yui_kh::khi::internal::v2::builder::SymTngBuilder;
use yui_link::{Crossing, InvLink, Link};
use yui_kh::khi::{ssi_invariants, KhIComplex, KhIHomology};

type R = FF2;
type P = HPoly<'H', R>;

fn init_logger() { 
    use yui::util::log::init_simple_logger;
    init_simple_logger(log::LevelFilter::Info).unwrap();
}

fn ssi_of(khi: &KhIHomology<P>) -> (isize, isize) {
    assert_eq!(khi[0].rank(), 2);
    assert_eq!(khi[1].rank(), 2);
    let s0 = (khi[0].gen(0).q_deg() + khi[0].gen(1).q_deg()) / 2;
    let s1 = (khi[1].gen(0).q_deg() + khi[1].gen(1).q_deg()) / 2;
    (s0, s1)
}

// run test by:
// cargo test -r -- --exact [NAME] --nocapture --include-ignored

#[test]
#[ignore]
fn k15n_103488() { 
    init_logger();

    let l = InvLink::sinv_knot_from_code(
        [[1,11,2,10],[2,20,3,19],[5,17,6,16],[6,25,7,26],[9,22,10,23],[12,30,13,29],[14,8,15,7],[15,27,16,26],[18,4,19,3],[20,11,21,12],[21,1,22,30],[23,4,24,5],[24,18,25,17],[27,8,28,9],[28,14,29,13]],
    );

    let khi = KhIHomology::new(&l, &P::variable(), &P::zero(), false);
    let (s0, s1) = ssi_of(&khi);

    khi.gen_grid().print_table("i", "j");
    
    println!("s0 = {s0}, s1 = {s1}");
    assert_eq!(s0, 0);
    assert_eq!(s1, 2);
}

#[test]
#[ignore]
fn knotJ() { 
    init_logger();
    
    let l = InvLink::sinv_knot_from_code([
        [1,27,2,26],[19,2,20,3],[3,13,4,12],[4,31,5,32],[30,5,31,6],
        [13,7,14,6],[8,27,9,28],[9,1,10,34],[10,18,11,17],[24,11,25,12],
        [14,21,15,22],[28,16,29,15],[33,16,34,17],[18,26,19,25],[20,8,21,7],
        [29,23,30,22],[23,33,24,32]
    ]);

    let khi = KhIHomology::new(&l, &P::variable(), &P::zero(), false);
    let (s0, s1) = ssi_of(&khi);

    khi.gen_grid().print_table("i", "j");
    
    println!("s0 = {s0}, s1 = {s1}");
    assert_eq!(s0, 0);
    assert_eq!(s1, 2);
}

#[test]
#[ignore]
fn k9_46_interlock() { 
    init_logger();
    
    let l = InvLink::sinv_knot_from_code([
        [2,50,3,49],[4,48,5,47],[7,24,8,25],[9,52,10,53],[12,60,13,59],
        [14,58,15,57],[16,36,17,35],[18,34,19,33],[20,32,21,31],[21,40,22,41],
        [23,38,24,39],[25,6,26,7],[26,46,27,45],[28,44,29,43],[30,42,31,41],
        [32,20,33,19],[34,18,35,17],[37,54,38,55],[39,22,40,23],[42,30,43,29],
        [44,28,45,27],[46,6,47,5],[48,4,49,3],[50,2,51,1],[51,10,52,11],
        [53,8,54,9],[55,36,56,37],[56,16,57,15],[58,14,59,13],[60,12,1,11]
    ]);

    let khi = KhIHomology::new(&l, &P::variable(), &P::zero(), false);
    let (s0, s1) = ssi_of(&khi);

    khi.gen_grid().print_table("i", "j");
    
    println!("s0 = {s0}, s1 = {s1}");
    assert_eq!(s0, 0);
    assert_eq!(s1, 4);
}

#[test]
#[ignore]
fn k15n_interlock() { 
    init_logger();
    
    let l = InvLink::sinv_knot_from_code([
        [3,13,4,12],[4,68,5,67],[6,66,7,65],[9,34,10,35],[11,74,12,75],
        [13,3,14,2],[15,70,16,71],[18,82,19,81],[20,80,21,79],[22,50,23,49],
        [24,48,25,47],[26,44,27,43],[27,56,28,57],[29,58,30,59],[31,41,32,40],
        [33,52,34,53],[35,8,36,9],[36,64,37,63],[38,62,39,61],[41,31,42,30],
        [42,60,43,59],[45,55,46,54],[46,26,47,25],[48,24,49,23],[51,76,52,77],
        [53,32,54,33],[55,45,56,44],[57,28,58,29],[60,40,61,39],[62,38,63,37],
        [64,8,65,7],[66,6,67,5],[68,2,69,1],[69,14,70,15],[71,16,72,17],
        [73,83,74,82],[75,10,76,11],[77,50,78,51],[78,22,79,21],[80,20,81,19],
        [83,73,84,72],[84,18,1,17]
    ]);

    let khi = KhIHomology::new(&l, &P::variable(), &P::zero(), false);
    let (s0, s1) = ssi_of(&khi);

    khi.gen_grid().print_table("i", "j");
    
    println!("s0 = {s0}, s1 = {s1}");
    assert_eq!(s0, 0);
    assert_eq!(s1, 4);
}

#[test]
#[ignore]
fn knotJ_interlock_truncated() { 
    init_logger();
    
    let l = InvLink::sinv_knot_from_code([
        [1,17,2,16],[6,23,7,24],[7,115,8,114],[8,103,9,104],[9,3,10,2],
        [10,96,11,95],[14,47,15,48],[15,91,16,90],[18,113,19,114],[19,25,20,24],
        [20,5,21,6],[21,101,22,100],[26,112,27,111],[28,33,29,34],[29,73,30,72],
        [31,107,32,106],[34,72,35,71],[36,63,37,64],[37,43,38,42],[38,55,39,56],
        [39,83,40,82],[44,77,45,78],[45,61,46,60],[48,13,49,14],[49,93,50,92],
        [50,88,51,87],[56,41,57,42],[57,65,58,64],[58,85,59,86],[59,53,60,52],
        [61,77,62,76],[66,83,67,84],[67,55,68,54],[68,43,69,44],[69,63,70,62],
        [70,36,71,35],[74,107,75,108],[75,31,76,30],[78,53,79,54],[79,85,80,84],
        [80,65,81,66],[81,41,82,40],[86,52,87,51],[88,93,89,94],[89,13,90,12],
        [91,47,92,46],[94,12,95,11],[96,3,97,4],[97,103,98,102],[98,115,99,116],
        [99,23,100,22],[104,17,105,18],[105,1,106,120],[108,73,109,74],[109,33,110,32],
        [110,28,111,27],[116,101,117,102],[117,5,118,4],[118,25,119,26],[119,113,120,112]
    ]);

    let n = l.link().signed_crossing_nums().1 as isize;
    let (h, t) = (P::variable(), P::zero());
    let mut b = SymTngBuilder::new(&l, &h, &t, false);

    b.set_h_range(-n ..= 2);
    b.process_partial(0..26);
    b.process_partial(0..22);
    b.finalize();
    
    let c = b.into_khi_complex().truncated(-n..=2);
    c.check_d_all();

    let khi = KhIHomology::new(&l, &P::variable(), &P::zero(), false);
    let (s0, s1) = ssi_of(&khi);

    khi.gen_grid().print_table("i", "j");
    
    println!("s0 = {s0}, s1 = {s1}");
    assert_eq!(s0, 0);
    assert_eq!(s1, 4);
}