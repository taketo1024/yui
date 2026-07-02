// Constructions on `Link` kept out of the core type.

use itertools::Itertools;

use crate::{Link, Edge, LinkBuilder};

impl Link {
    // Connected sum at the two base points (a PD-built link defaults to its minimal edge).
    pub fn conn_sum(&self, other: &Link) -> Link {
        let self_e = self.base_pt().expect("self needs a base point");
        let other_e = other.base_pt().expect("other needs a base point");
        self.conn_sum_at(other, self_e, other_e)
    }

    // Connected sum at edges `self_e`/`other_e`: import both links, then re-splice the two edges
    // tail-to-head (self out → other in, other out → self in), so the orientations are compatible.
    pub fn conn_sum_at(&self, other: &Link, self_e: Edge, other_e: Edge) -> Link {
        assert!(self.is_oriented(),  "conn_sum requires an oriented link (self)");
        assert!(other.is_oriented(), "conn_sum requires an oriented link (other)");

        let mut b = LinkBuilder::new();
        let v1 = b.add_link(self);
        let v2 = b.add_link(other);

        // the (tail, head) ends of each spliced edge, as builder ports.
        let port = |verts: &[_], (i, s): (usize, usize)| (verts[i], s);
        let (t1, h1) = self.edge_ports(self_e);
        let (t2, h2) = other.edge_ports(other_e);
        let (t1, h1) = (port(&v1, t1), port(&v1, h1));
        let (t2, h2) = (port(&v2, t2), port(&v2, h2));

        b.disconnect(h1);   // free self_e's two ports
        b.disconnect(h2);   // free other_e's two ports
        b.connect(t1, h2);  // self tail → other head
        b.connect(t2, h1);  // other tail → self head

        b.build().unwrap()
    }

    // The (tail, head) ports of an oriented edge: the tail is where the strand exits its node,
    // the head where it enters the next (cf. `NodeOri::in_ports`).
    fn edge_ports(&self, e: Edge) -> ((usize, usize), (usize, usize)) {
        let (x, y) = self.nodes().enumerate().flat_map(|(i, n)|
            (0..4).filter(move |&s| n.edge(s) == e).map(move |s| (i, s))
        ).collect_tuple().unwrap_or_else(||
            panic!("edge {e} must appear exactly twice")
        );

        let is_in = |(i, s): (usize, usize)| {
            let ports = self.node(i).ori().in_ports().expect("edge_ports requires an oriented link");
            ports.contains(&s)
        };
        debug_assert!(is_in(x) != is_in(y), "edge {e} must have one head and one tail");
        if is_in(x) {
            (y, x)
        } else {
            (x, y)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Braid;
    use crate::misc::jones_polynomial;

    #[test]
    fn conn_sum_is_jones_multiplicative() {
        // unreduced Jones: Ṽ(K1 # K2) · Ṽ(unknot) = Ṽ(K1) · Ṽ(K2)
        let k1 = Link::test_data("3_1");
        let k2 = Link::test_data("4_1");
        let cs = k1.conn_sum(&k2); // at the default base points (edge 1 each)
        assert_eq!(cs.n_comps(), 1);
        assert_eq!(cs.n_crossings(), k1.n_crossings() + k2.n_crossings());
        assert!(cs.is_oriented());

        let (vcs, vu) = (jones_polynomial(&cs), jones_polynomial(&Link::unknot()));
        let (v1, v2) = (jones_polynomial(&k1), jones_polynomial(&k2));
        assert_eq!(&vcs * &vu, &v1 * &v2);

        // the connected sum does not depend on the chosen edges.
        let cs2 = k1.conn_sum_at(&k2, 4, 6);
        assert_eq!(jones_polynomial(&cs2), jones_polynomial(&cs));
    }

    #[test]
    fn conn_sum_of_braid_closures() {
        // braid closures orient downward (`ori = Down`), exercising the non-PD arms of
        // `NodeOri::in_ports` in edge_ports.
        let k1 = Braid::from([1, 1, 1]).closure();      // 3_1
        let k2 = Braid::from([1, -2, 1, -2]).closure(); // 4_1
        let cs = k1.conn_sum(&k2);
        assert_eq!(cs.n_comps(), 1);
        assert!(cs.is_oriented());

        let (vcs, vu) = (jones_polynomial(&cs), jones_polynomial(&Link::unknot()));
        let (v1, v2) = (jones_polynomial(&k1), jones_polynomial(&k2));
        assert_eq!(&vcs * &vu, &v1 * &v2);
    }
}
