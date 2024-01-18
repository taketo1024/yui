#![cfg(feature = "tex")]
pub trait TeX { 
    fn to_tex_string(&self) -> String;
}