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