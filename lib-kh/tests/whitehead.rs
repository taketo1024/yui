// Khovanov-homology verification of the Link constructions (braid closure, twist_knot, Whitehead
// doubles): the bigraded Kh is diagram-independent, so two diagrams of the same knot give identical
// homology (the mirror is tried where chirality is not the point).

use yui_link::{Link, InvLink, Braid, PDCodeX};
use yui_kh::kh::KhHomology;

fn kh(l: &Link) -> KhHomology<i32> {
    KhHomology::new(l, &0, &0, false)
}

fn same_knot(a: &Link, b: &Link) -> bool {
    let ka = kh(a);
    ka.is_identical(&kh(b)) || ka.is_identical(&kh(&b.mirror()))
}

// Authoritative symmetric PD of Wh⁺(4_1), 18 crossings.
fn wh_4_1() -> Vec<PDCodeX> {
    vec![
        [1, 9, 2, 8], [2, 11, 3, 12], [5, 25, 6, 24], [6, 31, 7, 32], [9, 29, 10, 28],
        [13, 33, 14, 32], [14, 23, 15, 24], [17, 13, 18, 12], [18, 7, 19, 8], [21, 5, 22, 4],
        [22, 15, 23, 16], [25, 21, 26, 20], [26, 35, 27, 36], [27, 11, 28, 10], [29, 1, 30, 36],
        [30, 19, 31, 20], [33, 17, 34, 16], [34, 3, 35, 4],
    ]
}

// Authoritative symmetric PD of Wh⁺(negative trefoil) — positive clasp, −6 twists, 20 crossings.
fn wh_3_1() -> Vec<PDCodeX> {
    vec![
        [3, 25, 4, 24], [4, 37, 5, 38], [7, 14, 8, 15], [9, 12, 10, 13], [11, 31, 12, 30],
        [13, 8, 14, 9], [17, 39, 18, 38], [18, 23, 19, 24], [21, 7, 22, 6], [22, 15, 23, 16],
        [25, 3, 26, 2], [26, 19, 27, 20], [27, 34, 28, 35], [29, 32, 30, 33], [31, 11, 32, 10],
        [33, 28, 34, 29], [35, 21, 36, 20], [36, 1, 37, 2], [39, 17, 40, 16], [40, 5, 1, 6],
    ]
}

#[test]
fn braid_closures_are_kh_faithful() {
    assert!(same_knot(&Braid::from([1, 1, 1]).closure(), &Link::test_data("3_1")));
    assert!(same_knot(&Braid::from([1, -2, 1, -2]).closure(), &Link::test_data("4_1")));
}

#[test]
fn twist_knot_is_correct() {
    // twist_knot(n) vs the independent PD-code reference, by Khovanov homology (mirror allowed).
    let table = ["3_1", "4_1", "5_2", "6_1", "7_2"];
    for (i, name) in table.iter().enumerate() {
        let n = i as i32 + 1;
        assert!(same_knot(&Link::twist_knot(n), &Link::test_data(name)), "twist_knot({n}) ≠ {name}");
        // negative is the mirror: twist_knot(-1-k) = mirror twist_knot(k).
        let neg = -1 - n;
        assert!(same_knot(&Link::twist_knot(neg), &Link::test_data(name)), "twist_knot({neg}) ≠ mirror {name}");
    }
    // n = 0 and n = -1 are the unknot (diagrams with crossings, so compare by homology).
    let unknot = Braid::from([1]).closure();
    assert!(same_knot(&Link::twist_knot(0), &unknot), "twist_knot(0) ≠ unknot");
    assert!(same_knot(&Link::twist_knot(-1), &unknot), "twist_knot(-1) ≠ unknot");
}

#[test]
fn whitehead_double_of_unknot_is_unknot() {
    // the untwisted Whitehead double of the unknot is the unknot, from any companion diagram.
    let unknot = Braid::from([1]).closure();
    for word in [vec![1, -2], vec![1, 2]] {
        let k = Braid::from_iter(word).closure();
        for positive in [true, false] {
            assert!(same_knot(&k.whitehead_double(positive, 0), &unknot), "D±(unknot) ≠ unknot");
        }
    }
}

#[test]
fn whitehead_double_matches_references() {
    // D⁺ of the companion must be Kh-identical (no mirror fallback: the chirality is the point)
    // to the authoritative reference. test_data("3_1") is the negative trefoil.
    let cases = [
        ("4_1", Link::test_data("4_1"), wh_4_1()),
        ("neg 3_1", Link::test_data("3_1"), wh_3_1()),
    ];
    for (name, companion, ref_pd) in cases {
        let dplus = kh(&companion.whitehead_double(true, 0));
        assert!(dplus.is_identical(&kh(&Link::from_pd_code(ref_pd))), "D⁺({name}) ≠ reference Wh⁺");
    }
}

#[test]
fn si_whitehead_double_recovers_references() {
    // The strongly-invertible Whitehead double must be a valid InvLink (its construction installs
    // the standard strong inversion) whose underlying knot is Kh-identical to the references.
    // InvLink::test_data("3_1") is the positive trefoil, so mirror it.
    let cases = [
        ("4_1", InvLink::test_data("4_1"), wh_4_1()),
        ("neg 3_1", InvLink::test_data("3_1").mirror(), wh_3_1()),
    ];
    for (name, comp, ref_pd) in cases {
        let si = comp.whitehead_double(true, 0);
        let reference = Link::from_pd_code(ref_pd);
        assert!(kh(si.inner()).is_identical(&kh(&reference)), "SI D⁺({name}) ≠ reference Wh⁺");
    }
}
