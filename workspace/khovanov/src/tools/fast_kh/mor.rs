use yui_core::{Ring, RingOps};
use yui_lin_comb::LinComb;

use super::cob::{Cob, Dot, Bottom, CobComp};
use super::tng::TngComp;

pub type Mor = LinComb<Cob, i32>; // Z-linear combination of cobordisms.

pub trait MorTrait: Sized {
    fn is_invertible(&self) -> bool;
    fn inv(&self) -> Option<Self>;
    fn map_cob<F>(self, f: F) -> Self where F: Fn(&mut Cob);
    fn connect(self, c: &Cob) -> Self;
    fn connect_comp(self, c: &CobComp) -> Self;
    fn cap_off(self, b: Bottom, c: &TngComp, dot: Dot) -> Self;
    fn eval<R>(&self, h: &R, t: &R) -> R
    where R: Ring + From<i32>, for<'x> &'x R: RingOps<R>;
}

impl MorTrait for Mor {
    fn is_invertible(&self) -> bool { 
        self.len() == 1 && 
        self.iter().next().map(|(c, a)| 
            c.is_invertible() && a.is_unit()
        ).unwrap_or(false)
    }

    fn inv(&self) -> Option<Self> { 
        if let Some((Some(cinv), Some(ainv))) = self.iter().next().map(|(c, a)| 
            (c.inv(), a.inv())
        ) { 
            let inv = Mor::from((cinv, ainv));
            Some(inv)
        } else { 
            None
        }
    }

    fn map_cob<F>(self, f: F) -> Self 
    where F: Fn(&mut Cob) {
        self.into_map(|mut cob, r| { 
            f(&mut cob);
            if cob.is_zero() { 
                (cob, 0)
            } else { 
                (cob, r)
            }
        })
    }

    fn connect(self, c: &Cob) -> Self {
        self.map_cob(|cob| cob.connect(c.clone()) )
    }

    fn connect_comp(self, c: &CobComp) -> Self {
        self.map_cob(|cob| cob.connect_comp(c.clone()) )
    }

    fn cap_off(self, b: Bottom, c: &TngComp, dot: Dot) -> Self {
        self.map_cob(|cob| cob.cap_off(b, c, dot) )
    }

    fn eval<R>(&self, h: &R, t: &R) -> R
    where R: Ring + From<i32>, for<'x> &'x R: RingOps<R> {
        self.iter().map(|(c, &a)| { 
            R::from(a) * c.eval(h, t)
        }).sum()
    }
}