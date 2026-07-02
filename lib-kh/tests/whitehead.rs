// Khovanov-homology verification of the Link constructions (braid closure, twist_knot, Whitehead
// doubles): the bigraded Kh is diagram-independent, so two diagrams of the same knot give identical
// homology (the mirror is tried where chirality is not the point).

use yui_link::{Link, InvLink, Braid};
use yui_kh::kh::KhHomology;

mod common;
use common::{wh_4_1, wh_3_1};

fn kh(l: &Link) -> KhHomology<i32> {
    KhHomology::new(l, &0, &0, false)
}

fn same_knot(a: &Link, b: &Link) -> bool {
    let ka = kh(a);
    ka.is_identical(&kh(b)) || ka.is_identical(&kh(&b.mirror()))
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
