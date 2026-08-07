//! Polynomials and their monomials, in one, two, three or indexed-many
//! variables, ordinary or Laurent.

mod var;
mod var2;
mod var3;
mod mdeg;
mod mvar;
mod mono;
mod poly;

pub use mdeg::MultiDeg;
pub use var::Var;
pub use var2::Var2;
pub use var3::Var3;
pub use mvar::MultiVar;
pub use mono::*;
pub use poly::*;