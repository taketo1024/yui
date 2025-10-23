use std::sync::Mutex;

pub struct SyncCounter { 
    count: Mutex<usize>
}

impl SyncCounter { 
    pub fn new() -> Self { 
        let count = Mutex::new(0);
        Self { count }
    }

    pub fn incr(&self) -> usize { 
        let mut c = self.count.lock().unwrap();
        *c += 1;
        *c
    }

    pub fn count(&self) -> usize { 
        *self.count.lock().unwrap()
    }
}