#[derive(Debug, derive_more::Display)]
pub struct Error { 
    pub msg: String
}

impl std::error::Error for Error {}

macro_rules! err {
    ($($arg:tt)*) => {{
        let msg = format!($($arg)*);
        let e = crate::utils::Error{ msg };
        Err( e.into() )
    }}
}

pub(crate) use err;