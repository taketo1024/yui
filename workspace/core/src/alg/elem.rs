use std::fmt::{Debug, Display};

pub trait ElemBase: 
    Default + 
    PartialEq + 
    Eq + 
    Clone + 
    Send + 
    Sync + 
    Display + 
    Debug + 
    'static
{}

impl<T> ElemBase for T where T: 
    Default + 
    PartialEq + 
    Eq + 
    Clone + 
    Send + 
    Sync + 
    Display + 
    Debug + 
    'static
{}

pub trait Elem: ElemBase { 
    fn math_symbol() -> String;
}