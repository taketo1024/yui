mod traits;
mod structs;

pub use traits::*;
pub use structs::*;

// univar (ordinary, Laurent)
pub type Poly  <const X: char, R> = PolyBase<Mono<X, usize>, R>;          
pub type LPoly <const X: char, R> = PolyBase<Mono<X, isize>, R>;

// bivar (ordinary, Laurent)
pub type Poly2 <const X: char, const Y: char, R> = PolyBase<Mono2<X, Y, usize>, R>;
pub type LPoly2<const X: char, const Y: char, R> = PolyBase<Mono2<X, Y, isize>, R>;

// multivar (ordinary, Laurent)
pub type PolyN <const X: char, R> = PolyBase<Mono<X, MultiDeg<usize>>, R>;
pub type LPolyN<const X: char, R> = PolyBase<Mono<X, MultiDeg<isize>>, R>;