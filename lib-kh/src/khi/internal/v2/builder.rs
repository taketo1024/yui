use std::collections::HashSet;
use std::ops::RangeInclusive;
use ahash::AHashMap;
use cartesian::cartesian;
use delegate::delegate;
use itertools::Itertools;
use log::info;
use rayon::prelude::*;
use yui::bitseq::{Bit, BitSeq};
use yui::{KeyedUnionFind, Ring, RingOps};
use yui_homology::DisplaySeq;
use yui_link::{Node, Edge, InvLink};

use crate::kh::{KhComplex, KhChainGen, KhTensor};
use crate::khi::KhIComplex;
use crate::kh::internal::v2::builder::{BuildElem, TngComplexBuilder};
use crate::kh::internal::v2::cob::LcCobTrait;
use crate::kh::internal::v2::tng::TngComp;
use crate::kh::internal::v2::tng_complex::{TngComplex, TngKey};

pub struct SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> {
    inner: TngComplexBuilder<R>,
    x_map: AHashMap<Node, Node>,
    e_map: AHashMap<Edge, Edge>,
    key_map: AHashMap<TngKey, TngKey>,
    pub auto_deloop: bool,
    pub auto_elim: bool
}

impl<R> SymTngBuilder<R> 
where R: Ring, for<'x> &'x R: RingOps<R> { 
    pub fn build_kh_complex(l: &InvLink, h: &R, t: &R, reduced: bool) -> KhComplex<R> { 
        let mut b = Self::new(l, h, t, reduced);
        b.preprocess();
        b.process_all();
        b.finalize();
        b.into_kh_complex()
    }

    pub fn build_khi_complex(l: &InvLink, h: &R, t: &R, reduced: bool) -> KhIComplex<R> { 
        let mut b = Self::new(l, h, t, reduced);
        b.preprocess();
        b.process_all();
        b.finalize();
        b.into_khi_complex()
    }

    pub fn new(l: &InvLink, h: &R, t: &R, reduced: bool) -> SymTngBuilder<R> { 
        assert!(l.link().data().iter().all(|x| !x.is_resolved()));
        assert!(!reduced || l.base_pt().is_some());

        let base_pt = if reduced { l.base_pt() } else { None };
        let inner = TngComplexBuilder::new(l.link(), h, t, base_pt);
        let x_map = l.link().data().iter().map(|x| (x.clone(), l.inv_x(x).clone())).collect();
        let e_map = l.link().edges().iter().map(|&e| (e, l.inv_e(e))).collect();

        Self::new_impl(inner, x_map, e_map)
    }

    fn new_impl(mut inner: TngComplexBuilder<R>, x_map: AHashMap<Node, Node>, e_map: AHashMap<Edge, Edge>) -> Self { 
        inner.auto_deloop = false;
        inner.auto_elim = false;

        let key_map = AHashMap::from_iter([(TngKey::init(), TngKey::init())]);
        let auto_deloop = true;
        let auto_elim = true;
        
        SymTngBuilder { inner, x_map, e_map, key_map, auto_deloop, auto_elim }
    }

    delegate! { 
        to self.inner { 
            pub fn complex(&self) -> &TngComplex<R>;
            pub fn crossings(&self) -> impl Iterator<Item = &Node>;
            pub fn set_crossings<I>(&mut self, crossings: I) where I: IntoIterator<Item = Node>;
            pub fn elements(&self) -> impl Iterator<Item = &BuildElem<R>>;
            pub fn set_elements<I>(&mut self, elements: I) where I: IntoIterator<Item = BuildElem<R>>;
            pub fn set_h_range(&mut self, h_range: RangeInclusive<isize>);
            pub(crate) fn stat(&self) -> String;
        }
    }

    pub fn preprocess(&mut self) { 
        assert_eq!(self.complex().dim(), 0, "must start from init state.");

        let elements = self.inner.take_elements();
        let off_axis = self.off_axis_crossings(true).into_iter().cloned().collect_vec();

        info!("({}) preprocess off-axis: {}", self.stat(), off_axis.len());

        let (c, tc, key_map, elements) = self.build_from_half(off_axis.iter(), elements);

        self.inner.remove_crossings(off_axis.iter());
        self.inner.connect(c);

        let off_axis = off_axis.iter().map(|x| self.inv_x(x).clone()).collect_vec();
        self.inner.remove_crossings(off_axis.iter());
        self.inner.connect(tc);
        
        self.key_map = key_map;
        self.set_elements(elements);

        info!("({}) preprocess done.", self.stat());
    }

    fn off_axis_crossings(&self, take_half: bool) -> Vec<&Node> { 
        let off_axis = self.crossings().filter(|&x|
            self.inv_x(x) != x
        ).collect_vec();

        if !take_half { 
            return off_axis
        }

        let is_adj = |x: &Node, y: &Node| {
            x.edges().iter().filter(|&&e| 
                self.inv_e(e) != e // no axis-crossing edge
            ).any(|e| 
                y.edges().contains(e)
            )
        };

        let mut u = KeyedUnionFind::from_iter(off_axis.iter().cloned());

        for (i, x) in off_axis.iter().enumerate() { 
            for j in 0 .. i { 
                let y = &off_axis[j];
                if is_adj(x, y) { 
                    u.union(&x, &y);
                }
            }
        }

        u.group().into_iter().fold(vec![], |mut res, next| { 
            if let Some(x) = next.iter().next() { 
                let tx = self.inv_x(x);
                if !res.contains(&tx) { 
                    res.extend(next.into_iter());
                }
            }
            res
        })
    }

    fn build_from_half<'a, I>(&self, crossings: I, elements: Vec<BuildElem<R>>) -> (TngComplex<R>, TngComplex<R>, AHashMap<TngKey, TngKey>, Vec<BuildElem<R>>) 
    where I: IntoIterator<Item = &'a Node> { 
        let (h, t) = self.complex().ht();
        let mut b = TngComplexBuilder::init(h, t, (0, 0), None);

        b.set_crossings(crossings.into_iter().cloned());
        b.set_elements(elements);
        b.process_all();

        // take results
        let keys = b.complex().keys().cloned().collect_vec();
        let elements = b.take_elements();
        let c = b.into_tng_complex();
        let tc = c.convert_edges(|e| self.inv_e(e));

        // build keys
        let keys = cartesian!(keys.iter(), keys.iter()).collect_vec();
        let key_map = keys.into_par_iter().map(|(k1, k2)| { 
            let k  = k1 + k2;
            let tk = k2 + k1;
            (k, tk)
        }).collect::<ahash::HashMap<_, _>>().into();

        // build elements
        let elements = elements.into_par_iter().map(|mut e| { 
            e.modify(|k, c| {
                let kk = k + k;
                let cc = c.into_map(|mut c, r| { 
                    let r = &r * &r;
                    let tc = c.convert_edges(|e| self.inv_e(e));
                    c.connect(tc);
                    (c, r)
                });
                (kk, cc)
            });
            e
        }).collect::<Vec<_>>();

        (c, tc, key_map, elements)
    }

    pub fn process_all(&mut self) { 
        info!("({}) process {} crossings", self.stat(), self.crossings().count());
        
        while let Some(x) = self.inner.choose_next() { 
            let tx = self.inv_x(&x);
            if &x == tx { 
                self.append_on_axis(&x);
            } else { 
                self.append_off_axis(&x, &tx.clone());
            }
        }
    }

    fn append_on_axis(&mut self, x: &Node) { 
        info!("({}) append on-axis: {x}", self.stat());

        self.inner.append_prepare(&x);

        let c = self.complex().make_x(x);
        let key_map = if !x.is_resolved() { 
            [Bit::Bit0, Bit::Bit1].map(|b| { 
                let k = TngKey { state: BitSeq::from(b), label: KhTensor::empty() };
                (k, k)
            }).into_iter().collect()
        } else { 
            let k = TngKey::init();
            [(k, k)].into_iter().collect()
        };

        self.connect(c, key_map);
    }

    fn append_off_axis(&mut self, x: &Node, tx: &Node) { 
        assert_eq!(self.inv_x(x), tx);

        info!("({}) append off-axis: {x}, {tx}", self.stat());

        self.inner.append_prepare(&x);
        self.inner.append_prepare(&tx);

        let c = { 
            let mut c = self.complex().make_x(x);
            c.append(tx);
            c
        };
        let key_map = if !x.is_resolved() { 
            [
                ([0, 0], [0, 0]),
                ([1, 0], [0, 1]),
                ([0, 1], [1, 0]),
                ([1, 1], [1, 1])
            ].map(|(b0, b1)| { 
                let k = TngKey { state: BitSeq::from_iter(b0), label: KhTensor::empty() };
                let l = TngKey { state: BitSeq::from_iter(b1), label: KhTensor::empty() };
                (k, l)
            }).into_iter().collect()
        } else {
            let k = TngKey::init();
            [(k, k)].into_iter().collect()
        };

        self.connect(c, key_map);
    }

    pub fn deloop_all(&mut self, allow_based: bool) { 
        for i in self.complex().h_range() { 
            self.deloop_in(i, allow_based);
        }
    }

    fn deloop_in(&mut self, i: isize, allow_based: bool) {
        let mut keys = self.complex().keys_of(i).filter(|k| 
            self.complex().vertex(k).tng().contains_circle()
        ).cloned().collect::<HashSet<_>>();

        if keys.is_empty() { return }

        info!("({}) C[{i}]: {}, deloop: {}.", self.stat(), self.complex().rank(i), self.inner.count_loops_in(i, allow_based));
        let before = self.complex().rank(i);

        while let Some((k, r)) = self.inner.find_loop(keys.iter(), allow_based) { 
            keys.remove(&k);
            keys.remove(&self.inv_key(&k));

            let updated = self.deloop_equiv(&k, r);
            
            keys.extend(updated.into_iter().filter(|k| self.complex().contains_key(k)));
        }

        let after = self.complex().rank(i);
        info!("({}) -> C[{i}]: {} (+{}).", self.stat(), after, after - before);

        if self.auto_elim { 
            self.eliminate_in(i - 1);
            self.eliminate_in(i);
        }
    }

    pub fn deloop_equiv(&mut self, k: &TngKey, r: usize) -> Vec<TngKey> { 
        if self.is_sym_key(&k) { 
            let c = self.complex().vertex(&k).tng().comp(r);
            if self.is_sym_comp(c) {
                // symmetric loop on symmetric key
                self.deloop_on_axis_sym(&k, r)
            } else {
                // asymmetric loop on symmetric key
                self.deloop_on_axis_asym(&k, r)
            }
        } else { 
            // (symmetric or asymmetric) loop on asymmetric key
            self.deloop_off_axis(&k, r)
        }
    }

    fn deloop_on_axis_sym(&mut self, k: &TngKey, r: usize) -> Vec<TngKey> {
        let c = self.complex().vertex(k).tng().comp(r);

        assert!(self.is_sym_key(k));
        assert!(self.is_sym_comp(c));

        let updated = self.inner.deloop(k, r);

        self.remove_key_pair(k);

        for &k_new in updated.iter() { 
            self.add_key_pair(k_new, k_new);
        }

        updated
    }

    #[allow(non_snake_case)]
    fn deloop_on_axis_asym(&mut self, k: &TngKey, r: usize) -> Vec<TngKey> {
        let c = self.complex().vertex(k).tng().comp(r);

        assert!(self.is_sym_key(k));
        assert!(!self.is_sym_comp(c));
        assert!(!self.complex().contains_base_pt(c));

        //          ⚪︎1 | ⚪︎1
        //  ⚪︎1 | ⚪︎X  <-->  ⚪︎X | ⚪︎1
        //          ⚪︎X | ⚪︎X

        let tc = c.convert_edges(|e| self.inv_e(e));

        let ks = self.inner.deloop(k, r);

        let (k_X, k_1) = (ks[0], ks[1]);
        let (k_XX, k_X1) = { 
            let tr = self.complex().vertex(&k_X).tng().index_of(&tc).unwrap();
            let tks = self.inner.deloop(&k_X, tr);
            (tks[0], tks[1])
        };
        let (k_1X, k_11) = { 
            let tr = self.complex().vertex(&k_1).tng().index_of(&tc).unwrap();
            let tks = self.inner.deloop(&k_1, tr);
            (tks[0], tks[1])
        };

        self.remove_key_pair(k);

        self.add_key_pair(k_XX, k_XX);
        self.add_key_pair(k_X1, k_1X);
        self.add_key_pair(k_11, k_11);

        vec![k_XX, k_X1, k_1X, k_11]
    }

    #[allow(non_snake_case)]
    fn deloop_off_axis(&mut self, k: &TngKey, r: usize) -> Vec<TngKey> {
        let c = self.complex().vertex(k).tng().comp(r);

        assert!(!self.is_sym_key(k));

        //  ⚪︎1 | ..  <-->  .. | ⚪︎1
        //  ⚪︎X | ..  <-->  .. | ⚪︎X

        let tk = *self.inv_key(k);
        let tc = c.convert_edges(|e| self.inv_e(e));
        let tr = self.complex().vertex(&tk).tng().index_of(&tc).unwrap();

        let mut ks = self.inner.deloop(k, r);
        let mut tks = self.inner.deloop(&tk, tr);

        self.remove_key_pair(k);

        for (&k_new, &tk_new) in Iterator::zip(ks.iter(), tks.iter()) { 
            self.add_key_pair(k_new, tk_new);
        }

        ks.append(&mut tks);
        ks
    }

    pub fn eliminate_all(&mut self) { 
        for i in self.complex().h_range() { 
            self.eliminate_in(i)
        }
    }

    fn eliminate_in(&mut self, i: isize) { 
        let mut keys = self.complex().keys_of(i).filter(|k| 
            self.complex().keys_out_from(k).find(|l|
                self.is_equiv_inv_edge(k, l)
            ).is_some()
        ).cloned().collect::<HashSet<_>>();

        if keys.is_empty() { return }

        info!("({}) C[{i}]: {}, elim targets: {}.", self.stat(), self.complex().rank(i), keys.len());
        let before = self.complex().rank(i);

        while let Some((k, l, _)) = self.choose_pivot(keys.iter()) { 
            let (k, l) = (*k, *l);
            let tk = *self.inv_key(&k);

            self.eliminate_equiv(&k, &l);
            
            keys.remove(&k);
            keys.remove(&tk);
        }            

        let after = self.complex().rank(i);
        info!("({}) -> C[{i}]: {} (-{}).", self.stat(), after, before - after);
    }

    pub fn eliminate_equiv(&mut self, i: &TngKey, j: &TngKey) {
        assert_eq!(self.is_sym_key(i), self.is_sym_key(j));
        assert!(self.complex().has_edge(&i, &j));

        if self.is_sym_key(i) { 
            self.inner.eliminate(i, j);
        } else { 
            let ti = *self.inv_key(i);
            let tj = *self.inv_key(j);

            assert!(self.complex().has_edge(&ti, &tj));

            self.inner.eliminate(i, j);
            self.inner.eliminate(&ti, &tj);
        }

        self.remove_key_pair(i);
        self.remove_key_pair(j);
    }

    fn choose_pivot<'a, I>(&self, keys: I) -> Option<(&'a TngKey, &TngKey, usize)> 
    where I: IntoIterator<Item = &'a TngKey> { 
        keys.into_iter().filter_map(move |k|
            self.choose_pivot_col(k).map(move |(l, s)| (k, l, s))
        ).min_by_key(|(_, _, s)| *s)
    }

    fn choose_pivot_col(&self, k: &TngKey) -> Option<(&TngKey, usize)> { 
        self.complex().keys_out_from(k).filter_map(|l|
            self.is_equiv_inv_edge(k, l).then(|| {
                let s = self.inner.edge_weight(k, l);
                (l, s)
            })
        )
        .min_by_key(|(_, s)| *s)
    }

    fn is_equiv_inv_edge(&self, i: &TngKey, j: &TngKey) -> bool { 
        let f = self.complex().edge(i, j);
        f.is_invertible() && self.is_equiv_edge(i, j)
    }

    fn is_equiv_edge(&self, i: &TngKey, j: &TngKey) -> bool { 
        if self.is_sym_key(i) && self.is_sym_key(j) { 
            true
        } else if !self.is_sym_key(i) && !self.is_sym_key(j) { 
            //  i - - -> j 
            //    \   /   
            //      /     : not allowed
            //    /   \   
            // ti - - -> tj
            let ti = self.inv_key(i);
            let tj = self.inv_key(j);

            !self.complex().keys_into(j).contains(ti) && 
            !self.complex().keys_into(tj).contains(i)
        } else { 
            false
        }
    }

    pub fn finalize(&mut self) {
        if self.complex().is_completely_delooped() { 
            return
        }

        info!("({}) finalize", self.stat());

        self.deloop_all(false);
        self.deloop_all(true);

        assert!(self.complex().is_completely_delooped());
    }

    pub fn process_partial<I>(&mut self, indices: I)
    where I: IntoIterator<Item = usize> { 
        self.set_elements([]); // TODO

        // extract target crossings
        let indices = indices.into_iter().collect::<HashSet<_>>();
        let target = self.crossings().enumerate().filter(|(i, _)|
            indices.contains(&i)
        ).map(|(_, x)| 
            x.clone()
        ).collect_vec();
        
        self.inner.remove_crossings(target.iter());

        let (h, t) = self.complex().ht();
        let mut b = Self::new_impl(
            TngComplexBuilder::init(h, t, (0, 0), None),
            self.x_map.clone(),
            self.e_map.clone()
        );

        b.set_crossings(target);
        b.preprocess();
        b.process_all();

        let key_map = std::mem::take(&mut b.key_map);
        let c = b.into_tng_complex();

        info!("connect ({}) <- ({})", self.stat(), c.stat());

        self.connect(c, key_map);

        info!("connected ({})", self.stat());
    } 

    fn connect(&mut self, c: TngComplex<R>, key_map: AHashMap<TngKey, TngKey>) { 
        let (left, right) = self.inner.connect_init(c);
        let l_key_map = std::mem::take(&mut self.key_map);
        let r_key_map = key_map;

        let h_range = self.inner.current_h_range();

        for i in h_range.clone() { 
            self.inner.complex_mut().connect_vertices(&left, &right, i);
        }

        self.key_map = cartesian!(
            l_key_map.iter(),
            r_key_map.iter()
        ).map(|((k1, l1), (k2, l2))|
            (k1 + k2, l1 + l2)
        ).filter(|(k, l)|
            self.complex().contains_key(&k) && 
            self.complex().contains_key(&l)
        ).collect();

        for i in h_range { 
            self.inner.complex_mut().connect_edges(&left, &right, i);
            if self.auto_deloop { 
                self.deloop_in(i, false);
            }
        }

        // TODO merge elements
    }

    pub fn into_tng_complex(self) -> TngComplex<R> { 
        self.inner.into_tng_complex()
    }

    pub fn into_kh_complex(self) -> KhComplex<R> {
        self.inner.into_kh_complex()
    }

    pub fn into_khi_complex(mut self) -> KhIComplex<R> {
        assert!(self.complex().is_completely_delooped());
        
        let deg_shift = self.complex().deg_shift();
        let key_map = std::mem::take(&mut self.key_map);

        let map = move |x: &KhChainGen| -> KhChainGen { 
            let k = TngKey::from(x);
            let tk = key_map[&k];
            let tx = tk.as_gen(deg_shift);
            tx
        };

        let c = self.into_kh_complex();

        info!("build KhI complex...");
        let c = KhIComplex::from_kh_complex(c, map);
        info!("  done\n{}", c.display_seq("i"));

        c
    }

    fn inv_x(&self, x: &Node) -> &Node { 
        &self.x_map[x]
    }

    fn inv_e(&self, e: Edge) -> Edge { 
        self.e_map[&e]
    }

    fn inv_key(&self, k: &TngKey) -> &TngKey { 
        &self.key_map[&k]
    }

    fn add_key_pair(&mut self, k: TngKey, tk: TngKey) { 
        self.key_map.insert(k, tk);
        if k != tk { 
            self.key_map.insert(tk, k);
        }
    }

    fn remove_key_pair(&mut self, k: &TngKey) { 
        let tk = self.key_map.remove(k).unwrap();
        if k != &tk { 
            self.key_map.remove(&tk);
        }
    }

    fn is_sym_key(&self, k: &TngKey) -> bool { 
        self.inv_key(k) == k
    }

    fn is_sym_comp(&self, c: &TngComp) -> bool { 
        &c.convert_edges(|e| self.inv_e(e)) == c
    }

    #[allow(unused)]
    fn print_keys(&self) {
        let mut done = HashSet::new();
        for k in self.key_map.keys().sorted() { 
            if done.contains(&k) { continue }

            let tk = self.inv_key(k);
            if k == tk {
                println!("{}", self.complex().vertex(k));
            } else { 
                println!("{} ↔ {}", self.complex().vertex(k), self.complex().vertex(tk));
            }

            done.insert(k);
            done.insert(tk);
        }
        println!();
    }

    #[allow(unused)]
    fn validate_equiv(&self) {
        for k in self.complex().keys().sorted() { 
            assert!(self.key_map.contains_key(k), "no inv-key for {k}");
            let tk = self.inv_key(k);

            for l in self.complex().keys_out_from(k) { 
                assert!(self.key_map.contains_key(l), "no inv-key for {l}");
                let tl = self.inv_key(l);

                assert!(self.complex().has_edge(tk, tl));

                let f = self.complex().edge(k, l);
                let tf = self.complex().edge(tk, tl);

                assert_eq!(&f.convert_edges(|e| self.inv_e(e)), tf);
            }
        }
    }
}

#[cfg(test)]
#[allow(unused)]
mod tests { 
    use crate::khi::KhIHomology;

    use super::*;
    use num_traits::Zero;

    use yui::num::FF2;
    use yui::poly::HPoly;
    use yui_homology::{ChainComplexTrait, DisplaySeq, DisplayTable, SummandTrait};

    #[test]
    fn test_kh_3_1() { 
        let l = InvLink::load("3_1").unwrap();
        let (h, t) = (FF2::zero(), FF2::zero());

        let c = SymTngBuilder::build_kh_complex(&l, &h, &t, false);
        c.check_d_all();

        let h = c.inner().homology();

        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 2);
    }

    #[test]
    fn test_khi_3_1() { 
        let l = InvLink::load("3_1").unwrap();
        let (h, t) = (FF2::zero(), FF2::zero());

        let c = SymTngBuilder::build_khi_complex(&l, &h, &t, false);
        c.check_d_all();;

        let h = c.homology();

        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 4);
        assert_eq!(h[4].rank(), 2);
    }

    #[test]
    fn no_preprocess() { 
        let l = InvLink::load("3_1").unwrap();
        let (h, t) = (FF2::zero(), FF2::zero());

        let mut b = SymTngBuilder::new(&l, &h, &t, false);
        b.process_all();

        let c = b.into_khi_complex();
        c.check_d_all();;

        let h = c.inner().homology();

        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 2);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 4);
        assert_eq!(h[4].rank(), 2);
    }

    #[test]
    fn process_partial() { 
        let l = InvLink::sinv_knot_from_code([
            [6,9,7,10],[8,1,9,2],[14,7,1,8], // upper
            [3,13,4,12],[10,5,11,6],[11,3,12,2],[13,5,14,4], // lower
        ]); // 6_3

        let (h, t) = (FF2::zero(), FF2::zero());
        let mut b = SymTngBuilder::new(&l, &h, &t, false);

        b.auto_elim = false;
        b.auto_deloop = false;

        b.process_partial(0..3);
        b.process_partial(0..4);
        b.finalize();

        let c = b.into_khi_complex();
        c.check_d_all();

        assert!(c.canon_cycles().is_empty()); // TODO

        let h = c.homology();
        assert_eq!(h[-3].rank(), 2);
        assert_eq!(h[-2].rank(), 6);
        assert_eq!(h[-1].rank(), 8);
        assert_eq!(h[ 0].rank(), 10);
        assert_eq!(h[ 1].rank(), 10);
        assert_eq!(h[ 2].rank(), 8);
        assert_eq!(h[ 3].rank(), 6);
        assert_eq!(h[ 4].rank(), 2);
    }

    #[test]
    fn no_auto_deloop() { 
        let l = InvLink::load("3_1").unwrap();
        let (h, t) = (FF2::zero(), FF2::zero());

        let mut b = SymTngBuilder::new(&l, &h, &t, false);
        b.auto_deloop = false;
        b.preprocess();
        b.process_all();

        assert!(!b.complex().is_completely_delooped());

        b.finalize();

        assert!(b.complex().is_completely_delooped());

        let c = b.into_kh_complex();
        assert_eq!(c[0].rank(), 2);
        assert_eq!(c[1].rank(), 0);
        assert_eq!(c[2].rank(), 2);
        assert_eq!(c[3].rank(), 2);
        
        let h = c.homology();
        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 2);
    }

    #[test]
    fn no_auto_elim() { 
        let l = InvLink::load("3_1").unwrap();
        let (h, t) = (FF2::zero(), FF2::zero());

        let mut b = SymTngBuilder::new(&l, &h, &t, false);
        b.auto_elim = false;
        b.preprocess();
        b.process_all();

        assert!(b.complex().is_completely_delooped());

        let c = b.into_kh_complex();
        assert_eq!(c[0].rank(), 4);
        assert_eq!(c[1].rank(), 6);
        assert_eq!(c[2].rank(), 12);
        assert_eq!(c[3].rank(), 8);
        
        let h = c.homology();
        assert_eq!(h[0].rank(), 2);
        assert_eq!(h[1].rank(), 0);
        assert_eq!(h[2].rank(), 2);
        assert_eq!(h[3].rank(), 2);
    }

    #[test]
    fn h_range() { 
        let l = InvLink::load("6_3").unwrap();
        let (h, t) = (FF2::zero(), FF2::zero());
        let h_range = -2..=2;

        let mut b = SymTngBuilder::new(&l, &h, &t, false);
        b.set_h_range(h_range.clone());
        b.preprocess();
        b.process_all();
        b.finalize();

        let c = b.into_khi_complex().truncated(-1..=2);
        c.check_d_all();

        let h = c.homology().truncated(0..=1);
        assert_eq!(h[0].rank(), 10);
        assert_eq!(h[1].rank(), 10);
    }
}