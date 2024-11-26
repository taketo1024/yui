use super::app::{AppErr, err};
use log::info;

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

pub fn load<D>(name: &str) -> Result<D, Box<dyn std::error::Error>>
where for<'de> D: serde::Deserialize<'de> { 
    File::Result(name).read()
}

const RES_DIR: &str = "results";
const TMP_DIR: &str = "tmp";

pub enum File<'a> { 
    Result(&'a str),
    Tmp(&'a str)
}

impl<'a> File<'a> { 
    pub fn dir(&self) -> String { 
        let proj_dir = std::env!("CARGO_MANIFEST_DIR");
        let dir = match self {
            File::Result(_) => RES_DIR,
            File::Tmp(_)    => TMP_DIR,
        };
        format!("{proj_dir}/{dir}")
    }

    pub fn path(&self) -> String { 
        let dir = self.dir();
        let name = match self {
            File::Result(name) | 
            File::Tmp(name) => name,
        };
        let ext = "json";
        format!("{dir}/{name}.{ext}")
    }

    pub fn exists(&self) -> bool {
        std::path::Path::new(&self.path()).exists()
    }
    
    pub fn read<D>(&self) -> Result<D, Box<dyn std::error::Error>>
    where for<'de> D: serde::Deserialize<'de> { 
        let str = std::fs::read_to_string(&self.path())?;
        let data = serde_json::from_str::<D>(&str)?;
        Ok(data)
    }
    
    pub fn write<D>(&self, data: D) -> std::io::Result<()>
    where D: serde::Serialize {
        let dir = self.dir();
        let dir_path = std::path::Path::new(&dir);
        if !dir_path.exists() { 
            std::fs::create_dir(&dir_path)?;
        }

        if self.exists() { 
            info!("write (overwrite): {}", self.path());
        } else { 
            info!("write: {}", self.path());
        }    
        let json = serde_json::to_string(&data)?;
        std::fs::write(&self.path(), &json)?;
        Ok(())
    }

    pub fn delete(&self) -> std::io::Result<()> { 
        if self.exists() { 
            info!("delete: {}", self.path());
            std::fs::remove_file(&self.path())?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests { 
    use yui::hashmap;
    use yui_kr::KRHomologyStr;
    use super::*;

    #[test]
    fn data_exists() { 
        assert!( File::Result("3_1").exists());
        assert!(!File::Result("3_2").exists());
    }

    #[test]
    fn load() { 
        let res: Result<KRHomologyStr, _> = File::Result("3_1").read();
        
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), KRHomologyStr::from(hashmap!{ 
            (0,4,-2) => 1,
            (-2,2,2) => 1,
            (2,2,-2) => 1
        }));
    }
}