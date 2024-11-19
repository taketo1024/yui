#![cfg(feature = "tex")]
pub trait TeX { 
    fn tex_math_symbol() -> String;
    fn tex_string(&self) -> String;
}