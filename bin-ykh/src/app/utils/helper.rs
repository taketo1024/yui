#![allow(unused)]

use crate::app::err::*;
use std::str::FromStr;
use itertools::Itertools;
use num_traits::Zero;
use yui_core::{Ring, RingOps};
use yui_link::{Edge, InvLink, Link};
use yui_matrix::sparse::SpVec;

pub fn measure<F, Res>(proc: F) -> (Res, std::time::Duration) 
where F: FnOnce() -> Res { 
    let start = std::time::Instant::now();
    let res = proc();
    let time = start.elapsed();
    (res, time)
}

pub fn guard_panic<F, R>(f: F) -> Result<R, Box<dyn std::error::Error>>
where F: FnOnce() -> Result<R, Box<dyn std::error::Error>> + std::panic::UnwindSafe {
    std::panic::catch_unwind(|| {
        f()
    }).unwrap_or_else(|e| {
        let info = match e.downcast::<String>() {
            Ok(v) => *v,
            Err(e) => match e.downcast::<&str>() {
                Ok(v) => v.to_string(),
                _ => "Unknown Source of Error".to_owned()
            }
        };
        err!("panic: {info}")
    })
}

pub fn load_link(input: &String, mirror: bool) -> Result<Link, Box<dyn std::error::Error>> { 
    type PDCode = Vec<[Edge; 4]>;
    
    let l = { 
        if let Ok(pd_code) = serde_json::from_str::<PDCode>(input) { 
            Link::from_pd_code(pd_code)
        } else if let Ok(link) = Link::load(input) { 
            link
        } else { 
            return err!("invalid input link: '{}'", input);
        }
    };

    if mirror { 
        Ok(l.mirror())
    } else { 
        Ok(l)
    }
}

pub fn load_sinv_knot(input: &String, mirror: bool) -> Result<InvLink, Box<dyn std::error::Error>> { 
    type PDCode = Vec<[Edge; 4]>;
    
    let l = { 
        if let Ok(pd_code) = serde_json::from_str::<PDCode>(input) { 
            InvLink::sinv_knot_from_code(pd_code)
        } else if let Ok(link) = InvLink::load(input) { 
            link
        } else { 
            return err!("invalid input link: '{}'", input);
        }
    };

    if mirror { 
        Ok(l.mirror())
    } else { 
        Ok(l)
    }
}

pub fn parse_pair<R: FromStr + Zero>(s: &String) -> Result<(R, R), Box<dyn std::error::Error>> { 
    if let Ok(c) = R::from_str(s) { 
        return Ok((c, R::zero()))
    }

    let r = regex::Regex::new(r"^(.+),(.+)$").unwrap();
    if let Some(m) = r.captures(s) { 
        let (s1, s2) = (&m[1], &m[2]);
        if let (Ok(a), Ok(b)) = (R::from_str(s1), R::from_str(s2)) {
            return Ok((a, b))
        }
    }

    err!("cannot parse '{}' as {}.", s, std::any::type_name::<R>())
}

pub fn vec2str<R>(v: &SpVec<R>) -> String 
where R: Ring + ToString, for<'x> &'x R: RingOps<R> { 
    format!("({})", v.to_dense().iter().join(", "))
}

pub fn csv_writer(path: &String) -> Result<csv::Writer<std::fs::File>, Box<dyn std::error::Error>> { 
    use std::fs::OpenOptions;

    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .open(path)?;

    let wtr = csv::Writer::from_writer(file);

    Ok(wtr)
}