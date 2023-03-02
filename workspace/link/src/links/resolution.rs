#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Resolution { 
    Res0, Res1
}

use std::fmt::Display;

use Resolution::{Res0, Res1};

impl Resolution { 
    pub fn as_u8(&self) -> u8 {
        match self { 
            Res0 => 0,
            Res1 => 1
        }
    }

    pub fn is_zero(&self) -> bool {
        self == &Res0
    }
}

impl From<u8> for Resolution {
    fn from(a: u8) -> Self {
        match a { 
            0 => Res0,
            1 => Res1,
            _ => panic!()
        }
    }
}

impl Display for Resolution {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self { 
            Res0 => f.write_str("0"),
            Res1 => f.write_str("1")
        }
    }
}