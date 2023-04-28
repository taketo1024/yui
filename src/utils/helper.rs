use super::err;
use std::str::FromStr;
use num_traits::Zero;
use yui_link::{Link, Edge};

const RESOURCE_DIR: &str = "resources";

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

pub fn load_json<Output>(path: &String) -> Result<Output, Box<dyn std::error::Error>>
where Output: serde::de::DeserializeOwned { 
    let json = std::fs::read_to_string(path)?;
    let data: Output = serde_json::from_str(&json)?;
    Ok(data)
}

pub fn load_link(link: &String, mirror: bool) -> Result<Link, Box<dyn std::error::Error>> { 
    type PDCode = Vec<[Edge; 4]>;
    
    let l = { 
        let path = format!("{}/links/{}.json", RESOURCE_DIR, link);
        
        if let Ok(pd_code) = serde_json::from_str::<PDCode>(&link) { 
            Link::from_pd_code(pd_code)
        } else if std::path::Path::new(&path).exists() { 
            Link::load(&path)?
        } else { 
            return err!("invalid input link: '{}'", link);
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
    if let Some(m) = r.captures(&s) { 
        let (s1, s2) = (&m[1], &m[2]);
        if let (Ok(a), Ok(b)) = (R::from_str(s1), R::from_str(s2)) {
            return Ok((a, b))
        }
    }

    err!("cannot parse '{}' as {}.", s, std::any::type_name::<R>())
}

pub fn write_csv(path: &String, records: Vec<&String>) -> Result<(), Box<dyn std::error::Error>> { 
    use std::fs::OpenOptions;

    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .append(true)
        .open(&path)?;

    let mut wtr = csv::Writer::from_writer(file);

    wtr.write_record(records)?;
    wtr.flush()?;

    Ok(())
}