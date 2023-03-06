mod traits;
pub use traits::*;

mod structs;
pub use structs::*;

// univar (ordinary, Laurent)
pub type Poly  <const X: char, R> = PolyBase<Mono<X, usize>, R>;          
pub type LPoly <const X: char, R> = PolyBase<Mono<X, isize>, R>;

// bivar (ordinary, Laurent)
pub type Poly2 <const X: char, const Y: char, R> = PolyBase<Mono2<X, Y, usize>, R>;
pub type LPoly2<const X: char, const Y: char, R> = PolyBase<Mono2<X, Y, isize>, R>;

// multivar (ordinary, Laurent)
pub type MPoly <const X: char, R> = PolyBase<Mono<X, MDegree<usize>>, R>;
pub type MLPoly<const X: char, R> = PolyBase<Mono<X, MDegree<isize>>, R>;