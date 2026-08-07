//! [`Path`]: a connected component of a diagram, either an arc or a circle,
//! held as an oriented sequence of edges.

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arc_and_circ() {
        let a = Path::arc([1, 2, 3]);
        assert!(a.is_arc());
        assert!(!a.is_circle());
        assert_eq!(a.len(), 3);
        assert_eq!(a.edges(), &[1, 2, 3]);

        let c = Path::circ([1, 2, 3]);
        assert!(c.is_circle());
        assert!(!c.is_arc());
        assert_eq!(c.edges(), a.edges());

        // same edges, different kind — not equal.
        assert_ne!(a, c);
    }

    #[test]
    #[should_panic]
    fn arc_rejects_empty() {
        let _ = Path::arc([]);
    }

    #[test]
    #[should_panic]
    fn circ_rejects_empty() {
        let _ = Path::circ([]);
    }

    #[test]
    fn contains_and_min_edge() {
        let p = Path::arc([4, 2, 7]);
        assert!(p.contains(2));
        assert!(!p.contains(3));
        assert_eq!(p.min_edge(), 2);
    }

    #[test]
    fn end_pts_only_for_arcs() {
        assert_eq!(Path::arc([4, 2, 7]).end_pts(), Some((4, 7)));
        assert_eq!(Path::arc([5]).end_pts(), Some((5, 5)));
        assert_eq!(Path::circ([4, 2, 7]).end_pts(), None);
    }

    #[test]
    fn equality_is_oriented() {
        // `Path` compares as an oriented sequence: neither reversal nor rotation is equal.
        let p = Path::circ([1, 2, 3]);
        assert_ne!(p, Path::circ([3, 2, 1]));
        assert_ne!(p, Path::circ([2, 3, 1]));
    }

    #[test]
    fn into_seq_keeps_the_order() {
        assert_eq!(Path::circ([4, 2, 7]).into_seq().to_vec(), vec![4, 2, 7]);
    }

    #[test]
    fn display() {
        assert_eq!(Path::arc([1, 2, 3]).to_string(),  "[1-2-3]");
        assert_eq!(Path::circ([1, 2, 3]).to_string(), "⚪︎(1-2-3)");
    }
}
