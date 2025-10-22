use petgraph::Graph;
use crate::{Link, Path, CrossingType};

pub fn seifert_graph(link: &Link) -> Graph<Path, usize> { 
    type G = Graph<Path, usize>;
    let s0 = link.ori_pres_state();
    let l0 = link.resolve_all_crossings(&s0);
    let mut graph = Graph::new();

    for c in l0.components() {
        graph.add_node(c);
    }

    let find_node = |graph: &G, e| { 
        graph.node_indices().find(|&i| 
            graph[i].contains(e)
        )
    };

    for (i, x) in l0.data().iter().enumerate() {
        let (e1, e2) = if x.ctype() == CrossingType::V {
            (x.edge(0), x.edge(1))
        } else { 
            (x.edge(0), x.edge(2))
        };
        let n1 = find_node(&graph, e1).unwrap();
        let n2 = find_node(&graph, e2).unwrap();
        graph.add_edge(n1, n2, i);
    }

    graph
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn seif_graph() { 
        let l = Link::trefoil();
        let g = seifert_graph(&l);
        assert_eq!(g.node_count(), 2);
        assert_eq!(g.edge_count(), 3);
    }
}