use petgraph::algo::toposort;
use petgraph::graphmap::DiGraphMap;

/// Incremental builder for a directed graph over `usize` nodes.
/// Consume with [`into_sorted`](Self::into_sorted) for a topological ordering.
pub struct TopSort {
    g: DiGraphMap<usize, ()>,
}

impl TopSort {
    pub fn new() -> Self {
        Self { g: DiGraphMap::new() }
    }

    pub fn add_node(&mut self, i: usize) {
        self.g.add_node(i);
    }

    /// Add a directed edge `from → to`. Endpoints are auto-inserted if not present.
    pub fn add_edge(&mut self, from: usize, to: usize) {
        self.g.add_edge(from, to, ());
    }

    /// Consume self and return a topological ordering, or an error if the graph has a cycle.
    pub fn into_sorted(self) -> Result<Vec<usize>, String> {
        toposort(&self.g, None)
            .map_err(|c| format!("Input contains cycle at {}.", c.node_id()))
    }
}

impl Default for TopSort {
    fn default() -> Self {
        Self::new()
    }
}

impl<I> FromIterator<(usize, I)> for TopSort
where I: IntoIterator<Item = usize> {
    fn from_iter<T: IntoIterator<Item = (usize, I)>>(iter: T) -> Self {
        let mut ts = TopSort::new();
        for (i, succs) in iter {
            ts.add_node(i);
            for j in succs {
                ts.add_edge(i, j);
            }
        }
        ts
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn top_sort<Itr, I>(itr: Itr) -> Result<Vec<usize>, String>
    where Itr: IntoIterator<Item = (usize, I)>, I: IntoIterator<Item = usize> {
        TopSort::from_iter(itr).into_sorted()
    }

    #[test]
    fn valid() {
        // example from https://en.wikipedia.org/wiki/Topological_sorting
        let tree = vec![
            (5,  vec![11]),
            (7,  vec![11, 8]),
            (3,  vec![8, 10]),
            (11, vec![2, 9, 10]),
            (8,  vec![9]),
            (2,  vec![]),
            (9,  vec![]),
            (10, vec![])
        ];
        let res = top_sort(tree.clone());
        assert!(res.is_ok());

        let res = res.unwrap();
        assert!(tree.iter().all(|(i, list)| {
            list.iter().all(|j| {
                res.iter().position(|i1| i1 == i).unwrap()
                < res.iter().position(|j1| j1 == j).unwrap()
            })
        }));
    }

    #[test]
    fn cycle() {
        let tree = vec![
            (0, vec![1]),
            (1, vec![0]),
        ];
        assert!(top_sort(tree).is_err());
    }

    #[test]
    fn contains_cycle() {
        let tree = vec![
            (0, vec![1]),
            (1, vec![2]),
            (2, vec![1]),
        ];
        assert!(top_sort(tree).is_err());
    }

    #[test]
    fn incremental() {
        let mut ts = TopSort::new();
        ts.add_edge(0, 1);
        ts.add_edge(1, 2);
        ts.add_edge(0, 2);
        let sorted = ts.into_sorted().unwrap();
        assert_eq!(sorted, vec![0, 1, 2]);
    }

    #[test]
    fn incremental_cycle() {
        let mut ts = TopSort::new();
        ts.add_edge(0, 1);
        ts.add_edge(1, 0);
        assert!(ts.into_sorted().is_err());
    }
}
