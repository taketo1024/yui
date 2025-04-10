use std::hash::Hash;
use crate::Elem;

pub trait Gen: Elem + Hash + Ord {}