use super::*;
use super::Resolution::{Res0, Res1};
use super::CrossingType::{Xp, Xn, H, V};

fn a_crossing(ctype:CrossingType) -> Crossing {
    Crossing{
        ctype, 
        edges: [0,1,2,3]
    }
}

#[test]
fn test_crossing_is_resolved() {
    let c = a_crossing(Xp);
    assert!(!c.is_resolved());

    let c = a_crossing(Xn);
    assert!(!c.is_resolved());

    let c = a_crossing(H);
    assert!(c.is_resolved());

    let c = a_crossing(V);
    assert!(c.is_resolved());
}

#[test]
fn test_crossing_resolve() {
    let mut c = a_crossing(Xp);
    
    c.resolve(Res0);
    assert!(c.is_resolved());
    assert_eq!(c.ctype, V);

    let mut c = a_crossing(Xp);
    
    c.resolve(Res1);
    assert!(c.is_resolved());
    assert_eq!(c.ctype, H);

    let mut c = a_crossing(Xn);
    
    c.resolve(Res0);
    assert!(c.is_resolved());
    assert_eq!(c.ctype, H);

    let mut c = a_crossing(Xn);
    
    c.resolve(Res1);
    assert!(c.is_resolved());
    assert_eq!(c.ctype, V);
}

#[test]
fn test_crossing_mirror() {
    let mut c = a_crossing(Xp);
    c.mirror();
    assert_eq!(c.ctype, Xn);

    let mut c = a_crossing(Xn);
    c.mirror();
    assert_eq!(c.ctype, Xp);

    let mut c = a_crossing(H);
    c.mirror();
    assert_eq!(c.ctype, H);

    let mut c = a_crossing(V);
    c.mirror();
    assert_eq!(c.ctype, V);
}

#[test]
fn test_crossing_pass() {
    let c = a_crossing(Xp);
    assert_eq!(c.pass(0), 2);
    assert_eq!(c.pass(1), 3);
    assert_eq!(c.pass(2), 0);
    assert_eq!(c.pass(3), 1);

    let c = a_crossing(Xn);
    assert_eq!(c.pass(0), 2);
    assert_eq!(c.pass(1), 3);
    assert_eq!(c.pass(2), 0);
    assert_eq!(c.pass(3), 1);

    let c = a_crossing(V);
    assert_eq!(c.pass(0), 3);
    assert_eq!(c.pass(1), 2);
    assert_eq!(c.pass(2), 1);
    assert_eq!(c.pass(3), 0);

    let c = a_crossing(H);
    assert_eq!(c.pass(0), 1);
    assert_eq!(c.pass(1), 0);
    assert_eq!(c.pass(2), 3);
    assert_eq!(c.pass(3), 2);
}

#[test]
fn test_empty_link() { 
    let l = Link::empty();
    assert_eq!(l.data.len(), 0);
}

#[test]
fn test_link_from_pd_code() { 
    let pd_code = [[1,2,3,4]];
    let l = Link::from(pd_code);
    assert_eq!(l.data.len(), 1);
    assert_eq!(l.data[0].ctype, Xn);
}

#[test]
fn test_link_is_empty() {
    let l = Link::empty();
    assert!(l.is_empty());

    let pd_code = [[1,2,3,4]];
    let l = Link::from(pd_code);
    assert!(!l.is_empty());
}

#[test]
fn test_link_crossing_num() {
    let l = Link::empty();
    assert_eq!(l.crossing_num(), 0);

    let pd_code = [[1,2,3,4]];
    let l = Link::from(pd_code);
    assert_eq!(l.crossing_num(), 1);
}

#[test]
fn test_link_next() {
    let pd_code = [[0,0,1,1]];
    let l = Link::from(pd_code);

    assert_eq!(l.next(0, 0), Some((0, 1)));
    assert_eq!(l.next(0, 1), Some((0, 0)));
    assert_eq!(l.next(0, 2), Some((0, 3)));
    assert_eq!(l.next(0, 3), Some((0, 2)));

    let pd_code = [[0,3,1,4],[3,2,2,1]];
    let l = Link::from(pd_code);

    assert_eq!(l.next(0, 0), None);
    assert_eq!(l.next(0, 2), Some((1, 3)));
    assert_eq!(l.next(1, 1), Some((1, 2)));
    assert_eq!(l.next(1, 0), Some((0, 1)));
    assert_eq!(l.next(0, 3), None);
}

#[test]
fn test_link_traverse() {
    let pd_code = [[0,0,1,1]];
    let l = Link::from(pd_code);
    let mut queue = vec![];
    l.traverse_edges(0, 0, |i, j| queue.push((i, j)));
    
    assert_eq!(queue, [(0, 0), (0, 3), (0, 0)]); // loop

    let pd_code = [[0,3,1,4],[3,2,2,1]];
    let l = Link::from(pd_code);
    let mut queue = vec![];
    l.traverse_edges(0, 0, |i, j| queue.push((i, j)));

    assert_eq!(queue, [(0, 0), (1, 3), (1, 2), (0, 1), (0, 3)]); // no loop
}

#[test]
fn test_link_crossing_signs() {
    let pd_code = [[0,0,1,1]];
    let l = Link::from(pd_code);
    assert_eq!(l.crossing_signs(), vec![1]);

    let pd_code = [[0,1,1,0]];
    let l = Link::from(pd_code);
    assert_eq!(l.crossing_signs(), vec![-1]);
}

#[test]
fn test_link_writhe() {
    let pd_code = [[0,0,1,1]];
    let l = Link::from(pd_code);
    assert_eq!(l.writhe(), 1);

    let pd_code = [[0,1,1,0]];
    let l = Link::from(pd_code);
    assert_eq!(l.writhe(), -1);
}

#[test]
fn test_components() {
    let pd_code = [[0,0,1,1]];
    let l = Link::from(pd_code);
    let comps = l.components();
    assert_eq!(comps, vec![Component{ edges: vec![0, 1], closed: true }]);

    let pd_code = [[0,3,1,4],[3,2,2,1]];
    let l = Link::from(pd_code);
    let comps = l.components();
    assert_eq!(comps, vec![Component{ edges: vec![0,1,2,3,4], closed: false }]);
}

#[test]
fn test_2comp_unlink() {
    let pd_code = [[1,2,3,4], [3,2,1,4]];
    let l = Link::from(pd_code);
    assert_eq!(l.crossing_num(), 2);
    assert_eq!(l.crossing_signs(), [-1, 1]);
    assert_eq!(l.writhe(), 0);
}

#[test]
fn test_mirror() { 
    let pd_code = [[0,0,1,1]];
    let mut l = Link::from(pd_code);
    assert_eq!(l.data[0].ctype, Xn);

    l.mirror();

    assert_eq!(l.data[0].ctype, Xp);
}