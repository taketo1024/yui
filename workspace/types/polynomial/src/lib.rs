mod traits;
pub use traits::*;

mod structs;
pub use structs::*;

pub type Poly  <const X: char, R> = PolyBase<Mono<X, usize>, R>;          // univar
pub type LPoly <const X: char, R> = PolyBase<Mono<X, isize>, R>;          // univar, Laurent
pub type MPoly <const X: char, R> = PolyBase<Mono<X, MDegree<usize>>, R>; // multivar
pub type MLPoly<const X: char, R> = PolyBase<Mono<X, MDegree<isize>>, R>; // multivar, Laurent