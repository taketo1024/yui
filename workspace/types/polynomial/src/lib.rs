mod traits;
mod structs;

pub use traits::*;
pub use structs::*;

// univar (ordinary, Laurent)
pub type Poly  <const X: char, R> = PolyBase<Univar<X, usize>, R>;          
pub type LPoly <const X: char, R> = PolyBase<Univar<X, isize>, R>;

// bivar (ordinary, Laurent)
pub type Poly2 <const X: char, const Y: char, R> = PolyBase<BiVar<X, Y, usize>, R>;
pub type LPoly2<const X: char, const Y: char, R> = PolyBase<BiVar<X, Y, isize>, R>;

// multivar (ordinary, Laurent)
pub type PolyN <const X: char, R> = PolyBase<MultiVar<X, usize>, R>;
pub type LPolyN<const X: char, R> = PolyBase<MultiVar<X, isize>, R>;