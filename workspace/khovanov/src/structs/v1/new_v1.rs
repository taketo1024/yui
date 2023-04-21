use yui_core::{Ring, RingOps};
use yui_link::Link;

use super::canon_cycle::CanonCycles;
use super::cube::KhCube;

use crate::{KhComplex, KhChain};

impl<R> KhComplex<R>
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn new_v1(l: &Link, h: &R, t: &R, reduced: bool) -> Self { 
        let deg_shift = Self::deg_shift_for(l, reduced);
        let canon_cycles = if t.is_zero() && l.is_knot() {
            let ori = if reduced { vec![true] } else { vec![true, false] };
            ori.into_iter().map(|o| 
                KhChain::canon_cycle(l, &R::zero(), h, o)
            ).collect()
        } else { 
            vec![]
        };
        let cube = KhCube::new(l, h, t);
        let complex = cube.as_complex(deg_shift.0, reduced);

        KhComplex::_new(complex, canon_cycles, reduced, deg_shift)
    }        
}
