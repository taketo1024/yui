use std::collections::HashMap;
use std::ops::RangeInclusive;

use itertools::Itertools;
use num_integer::Integer;
use serde::ser::SerializeSeq;
use yui::poly::{LPoly, LPoly2, LPoly3, Var3, Var2, Mono};

pub type QPoly = LPoly<'q', i32>;
pub type QAPoly = LPoly2<'q', 'a', i32>;
pub type TQAPoly = LPoly3<'t', 'q', 'a', i32>;

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct KRHomologyStr(HashMap<(isize, isize, isize), Option<usize>>);

impl From<HashMap<(isize, isize, isize), usize>> for KRHomologyStr {
    fn from(value: HashMap<(isize, isize, isize), usize>) -> Self {
        Self(value.into_iter().map(|(k, v)| (k, Some(v))).collect())
    }
}

impl From<HashMap<(isize, isize, isize), Option<usize>>> for KRHomologyStr {
    fn from(value: HashMap<(isize, isize, isize), Option<usize>>) -> Self {
        Self(value.into_iter().collect())
    }
}

impl FromIterator<((isize, isize, isize), usize)> for KRHomologyStr {
    fn from_iter<T: IntoIterator<Item = ((isize, isize, isize), usize)>>(iter: T) -> Self {
        Self::from(iter.into_iter().collect::<HashMap<_, _>>())
    }
}

impl FromIterator<((isize, isize, isize), Option<usize>)> for KRHomologyStr {
    fn from_iter<T: IntoIterator<Item = ((isize, isize, isize), Option<usize>)>>(iter: T) -> Self {
        Self::from(iter.into_iter().collect::<HashMap<_, _>>())
    }
}

impl KRHomologyStr { 
    pub fn indices(&self) -> impl Iterator<Item = &(isize, isize, isize)> { 
        self.0.keys()
    }

    pub fn iter(&self) -> impl Iterator<Item = (&(isize, isize, isize), &Option<usize>)> { 
        self.0.iter()
    }

    pub fn is_determined(&self) -> bool { 
        self.0.values().all(Option::is_some)
    }

    pub fn iter_determined(&self) -> impl Iterator<Item = (&(isize, isize, isize), &usize)> { 
        self.0.iter().filter_map(|(idx, r)| r.as_ref().map(|r| (idx, r)))
    }

    pub fn non_determined(&self) -> impl Iterator<Item = &(isize, isize, isize)> { 
        self.0.iter().filter_map(|(idx, r)| 
            if r.is_some() { None } else { Some(idx) }
        )
    }

    pub fn get(&self, idx: (isize, isize, isize)) -> Option<usize> {
        // Some(Some(r)) -> Some(r) // determined
        // Some(None)    -> None    // non-determined
        // None          -> Some(0) // determined
        if let Some(e) = self.0.get(&idx) { 
            e.to_owned()
        } else { 
            Some(0)
        }
    }

    pub fn set(&mut self, idx: (isize, isize, isize), value: usize) {
        if value > 0 { 
            self.0.insert(idx, Some(value));
        } else { 
            self.0.remove(&idx);
        }
    }

    pub fn mirror(&self) -> KRHomologyStr { 
        self.iter().map(|(idx, &v)| ((-idx.0, -idx.1, -idx.2), v)).collect()
    }
    
    pub fn total_rank(&self) -> usize { 
        self.iter_determined().map(|(_, r)| r).sum()
    }

    pub fn delta(&self) -> Vec<isize> { 
        self.indices().map(|(i, j, k)| i + j + k).unique().sorted().collect_vec()
    }

    pub fn is_thin(&self) -> bool { 
        self.delta().len() <= 1
    }

    pub fn poincare_poly(&self) -> TQAPoly { 
        self.iter_determined().map(|(idx, &r)| { 
            let &(i, j, k) = idx;
            let x = Var3::from(((k - j)/2, i, j));
            let r = r as i32;
            (x, r)
        }).collect()
    }
    
    pub fn homfly_poly(&self) -> QAPoly { 
        self.poincare_poly().into_iter().map(|(x, r)|{
            let (t, i, j) = x.deg();
            let x = Var2::from((i, j)); // q^i a^j
            let r = if t.is_even() { r } else { -r }; // t ↦ -1
            (x, r)
        }).collect::<QAPoly>()
    }
    
    pub fn qpoly_table(&self) -> String {
        let polys = self.iter_determined().into_group_map_by(|(idx, _)|
            (idx.1, idx.2) // (j, k)
        ).into_iter().map(|(jk, list)| { 
            let q = QPoly::mono;
            let elems = list.into_iter().map(|(idx, &a)| {
                let i = idx.0;
                let a = a as i32;
                (q(i), a) // a.q^i
            });
            let p = QPoly::from_iter(elems);
            (jk, p)
        }).collect::<HashMap<(isize, isize), QPoly>>();

        let j_range = range(polys.keys().map(|idx| idx.0)).step_by(2);
        let k_range = range(polys.keys().map(|idx| idx.1)).rev().step_by(2);
    
        yui::util::format::table("k\\j", k_range, j_range, |&k, &j| { 
            if let Some(p) = polys.get(&(j, k)) { 
                p.to_string()
            } else { 
                ".".to_string()
            }
        })
    }

    pub fn delta_table(&self) -> String {
        range(self.delta()).step_by(2).map(|d| { 
            let indices = self.indices().filter(|(i, j, k)| i + j + k == d).collect_vec();
            let i_range = range(indices.iter().map(|idx| idx.0)).step_by(2);
            let j_range = range(indices.iter().map(|idx| idx.1)).rev().step_by(2);
        
            format!("Δ: {d}\n") + 
            &yui::util::format::table("j\\i", j_range, i_range, |&j, &i| { 
                let k = d - (i + j);
                if let Some(p) = self.0.get(&(i, j, k)) { 
                    if let Some(r) = p { r.to_string() } else { "?".to_string() }
                } else { 
                    ".".to_string()
                }
            })
        }).join("\n")
    }
}

impl serde::Serialize for KRHomologyStr {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where S: serde::Serializer {
        let mut seq = serializer.serialize_seq(None)?;
        for ((i, j, k), r) in self.iter() { 
            let e = (i, j, k, r);
            seq.serialize_element(&e)?;
        }
        seq.end()
    }
}

impl<'de> serde::Deserialize<'de> for KRHomologyStr {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where D: serde::Deserializer<'de> {
        let seq: Vec<(isize, isize, isize, Option<usize>)> = Vec::deserialize(deserializer)?;
        let str = seq.into_iter().map(|(i, j, k, r)| ((i, j, k), r)).collect();
        Ok(str)
    }
}

fn range<Itr>(itr: Itr) -> RangeInclusive<isize>
where Itr: IntoIterator<Item = isize> {
    if let Some((l, r)) = itr.into_iter().fold(None, |res, i| { 
        if let Some((mut l, mut r)) = res { 
            if i < l { l = i }
            if r < i { r = i }
            Some((l, r))
        } else { 
            Some((i, i))
        }
    }) { 
        l..=r
    } else { 
        0..=0
    }
}

#[cfg(test)]
mod tests {
    use yui::hashmap;
    use super::*;
 
    #[test]
    fn non_determined() { 
        let s = KRHomologyStr::from(hashmap!{             
            (0,4,0)  => Some(1),
            (-2,6,0) => Some(1),
            (-4,4,4) => Some(1),
            (4,4,-4) => Some(1),
            (2,6,-4) => None
        });

        assert!(!s.is_determined());
        assert_eq!(s.total_rank(), 4);
        assert_eq!(s.non_determined().collect_vec(), vec![&(2,6,-4)]);
    }

    #[test]
    fn serialize_non_determined() { 
        let s = KRHomologyStr::from(hashmap!{             
            (0,4,0)  => Some(1),
            (-2,6,0) => Some(1),
            (-4,4,4) => Some(1),
            (4,4,-4) => Some(1),
            (2,6,-4) => None
        });

        let ser = serde_json::to_string(&s).unwrap();
        let des: KRHomologyStr = serde_json::from_str(&ser).unwrap();
        assert_eq!(s, des);
    }
}