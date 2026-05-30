use std::fmt::Display;

use smallvec::SmallVec;

use crate::Edge;

/// Inline capacity for the variable-length `Path::Arc`/`Path::Circ` variants.
/// With `Edge = u8`, the 16-byte inline buffer (= heap repr's ptr+cap) fits
/// 16 elements — same struct size as `Vec<Edge>`, but skips allocation for
/// short paths.
pub type PathEdges = SmallVec<[Edge; 16]>;

/// An *oriented* connected component of a tangle: either an arc (`Arc`) or a
/// closed loop (`Circ`). `Path` compares equal as oriented edge sequences.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Path {
    Arc(PathEdges),
    Circ(PathEdges),
}

impl Path {
    pub fn arc<I>(edges: I) -> Self
    where I: IntoIterator<Item = Edge> {
        let edges: PathEdges = edges.into_iter().collect();
        assert!(!edges.is_empty());
        Self::Arc(edges)
    }

    pub fn circ<I>(edges: I) -> Self
    where I: IntoIterator<Item = Edge> {
        let edges: PathEdges = edges.into_iter().collect();
        assert!(!edges.is_empty());
        Self::Circ(edges)
    }

    pub fn is_arc(&self) -> bool { matches!(self, Self::Arc(_)) }
    pub fn is_circle(&self) -> bool { matches!(self, Self::Circ(_)) }

    pub fn contains(&self, e: Edge) -> bool {
        match self {
            Self::Arc(es) | Self::Circ(es) => es.contains(&e),
        }
    }

    pub fn len(&self) -> usize {
        match self {
            Self::Arc(es) | Self::Circ(es) => es.len(),
        }
    }

    pub fn edges(&self) -> &[Edge] {
        match self {
            Self::Arc(es) | Self::Circ(es) => &es[..],
        }
    }

    pub fn min_edge(&self) -> Edge {
        match self {
            Self::Arc(es) | Self::Circ(es) => *es.iter().min().unwrap(),
        }
    }

    pub fn end_pts(&self) -> Option<(Edge, Edge)> {
        match self {
            Self::Arc(es) => Some((es[0], *es.last().unwrap())),
            Self::Circ(_) => None,
        }
    }

    pub fn into_seq(self) -> PathEdges {
        match self {
            Self::Arc(es) | Self::Circ(es) => es,
        }
    }
}

impl Display for Path {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let c = self.edges().iter().map(|e| e.to_string()).collect::<Vec<_>>().join("-");
        if self.is_circle() {
            write!(f, "⚪︎({c})")
        } else {
            write!(f, "[{c}]")
        }
    }
}
