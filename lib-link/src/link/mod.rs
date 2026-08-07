//! Planar diagrams: the [`Link`] type, its nodes and components, PD-code
//! conversion, and the builder that assembles a diagram port by port.

mod link;
mod link_ops;
mod node;
mod path;
mod pd_code;
mod builder;
mod construct;

pub use link::{Link, Edge, State, StateRepr};
pub use pd_code::PDCodeX;
pub use node::{Node, NodeType, Slot};
pub use path::Path;
pub use builder::{LinkBuilder, LinkError, Port};
