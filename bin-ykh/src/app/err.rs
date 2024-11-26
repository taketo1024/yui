#[derive(Debug, derive_more::Display)]
pub struct Error { 
    pub msg: String
}

impl std::error::Error for Error {}

macro_rules! err {
    ($($arg:tt)*) => {{
        use crate::app::err::*;
        let e = Error{ msg: format!($($arg)*) };
        Err( e.into() )
    }}
}

macro_rules! ensure {
    ($cond:expr, $($arg:tt)*) => {{
        if !$cond { 
            return err!($($arg)*);
        }
    }}
}

pub(crate) use {err, ensure};