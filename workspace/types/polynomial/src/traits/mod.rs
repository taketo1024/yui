mod deg;
pub use deg::PolyDeg;
pub (crate) use deg::{impl_polydeg_signed, impl_polydeg_unsigned};

mod gen;
pub use gen::PolyGen;