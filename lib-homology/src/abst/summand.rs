use yui_core::{Ring, RingOps};

pub trait SummandTrait
where Self::R: Ring, for<'x> &'x Self::R: RingOps<Self::R> {
    type R;

    fn rank(&self) -> usize;
    fn tors(&self) -> &[Self::R];

    fn dim(&self) -> usize { // not a good name...
        self.rank() + self.tors().len()
    }

    fn is_zero(&self) -> bool { 
        self.rank() == 0 && self.is_free()
    }

    fn is_free(&self) -> bool { 
        self.tors().is_empty()
    }

    fn display(&self) -> String { 
        crate::rmod_str(self.rank(), self.tors())
    }
}