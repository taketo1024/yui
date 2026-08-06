//! [`ParseErr`]: the `FromStr::Err` of every math type in this crate.

use std::error::Error;
use std::fmt::{Display, Formatter, Result as FmtResult};

/// Failure to parse a math object from its string form. The message is carried
/// privately so variants can be added later without breaking the public API.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseErr {
    msg: String,
}

impl ParseErr {
    pub fn new(msg: impl Into<String>) -> Self {
        Self { msg: msg.into() }
    }

    /// The standard form: `cannot parse "{input}" as {target}`.
    pub fn invalid(input: &str, target: &str) -> Self {
        Self::new(format!("cannot parse \"{input}\" as {target}"))
    }

    pub fn msg(&self) -> &str {
        &self.msg
    }
}

impl Display for ParseErr {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "{}", self.msg)
    }
}

impl Error for ParseErr {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invalid_names_the_input_and_the_target() {
        let e = ParseErr::invalid("1/0", "Q");
        assert_eq!(e.to_string(), "cannot parse \"1/0\" as Q");
        assert_eq!(e.msg(), e.to_string());
    }

    #[test]
    fn is_a_std_error() {
        fn takes_err<E: Error>(_: E) {}
        takes_err(ParseErr::new("boom"));

        // the point of the newtype: `?` into a `dyn Error` chain, which `String` cannot do.
        fn f() -> Result<(), Box<dyn Error>> {
            Err(ParseErr::new("boom"))?
        }
        assert_eq!(f().unwrap_err().to_string(), "boom");
    }
}
