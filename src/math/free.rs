use std::fmt::{Display, Debug};
use std::hash::Hash;
use super::traits::Symbol;

pub trait FreeGenerator: Clone + PartialEq + Eq + Hash + Display + Debug + Send + Sync + Symbol {}