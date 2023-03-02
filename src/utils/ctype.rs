use clap::ValueEnum;

#[allow(non_camel_case_types)]
#[derive(Clone, ValueEnum, Debug, derive_more::Display)]
pub enum CType { 
    Z, Q, 
    F2, F3, 
    ZPoly, QPoly, 
    F2Poly, F3Poly, 
    ZMpoly, QMpoly, 
    F2Mpoly, F3Mpoly, 
    Gauss, Eisen,
    None
}