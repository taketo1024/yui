use std::str::FromStr;
use indexmap::IndexMap;
use log::{info, error};
use yui_link::{Link, Edge};
use yui_core::{EucRing, EucRingOps};
use yui_khovanov::ss_invariant;
use crate::utils::*;

pub fn run(name: String, l_str: Option<String>, c_value: String, c_type: CType) -> Result<String, Box<dyn std::error::Error>> {
    info!("compute ss: {name}, c = {c_value}");

    let l = load_link(&name, &l_str)?;
    let s = guard_panic(|| 
        dispatch_eucring!(c_type, compute_ss, &l, &c_value)
    )?;

    info!("{name}: s = {s} (c = {c_value})");

    Ok(s.to_string())
}

pub fn run_batch(c_value: String, c_type: CType, data: String, output: Option<String>) -> Result<String, Box<dyn std::error::Error>> {
    let data: IndexMap<String, Vec<[Edge; 4]>> = load_json(&data)?;
    let mut all_res = String::from("");

    for (name, code) in data { 
        let l = Link::from(&code);
        let res = guard_panic(|| 
            dispatch_eucring!(c_type, compute_ss, &l, &c_value)
        );

        let s = if let Ok(s) = res { 
            s.to_string()
        } else { 
            if let Err(e) = res {
                error!("{}", e);
                eprintln!("{e}");
            }
            "!".to_string()            
        };
        
        if let Some(output) = &output { 
            write_csv(output, vec![&name, &s])?;
        }

        if !all_res.is_empty() { 
            all_res += "\n";
        }
        
        all_res += &format!("{name}: {s}");
    }
    
    Ok(all_res)
}

fn compute_ss<R>(l: &Link, c_value: &String) -> Result<i32, Box<dyn std::error::Error>>
where R: EucRing + FromStr, for<'x> &'x R: EucRingOps<R> { 
    let Ok(c) = R::from_str(c_value) else { 
        return err!("cannot parse c: '{}' as type: {}.", c_value, std::any::type_name::<R>())
    };

    let ss = ss_invariant(l, &c, true);
    Ok(ss)
}

#[cfg(test)]
mod tests { 
    use super::*;

    #[test]
    fn test1() { 
        let name = "3_1".to_string();
        let c_value = "2".to_string();
        let c_type = CType::Z;

        let res = run(name, None, c_value, c_type);
        assert!(res.is_ok());
    }

    #[test]
    fn test2() { 
        let name = "4_1".to_string();
        let c_value = "3".to_string();
        let c_type = CType::Z;

        let res = run(name, None, c_value, c_type);
        assert!(res.is_ok());
    }
}