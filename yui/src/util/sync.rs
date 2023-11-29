use std::sync::Mutex;

pub fn sync_counter() -> impl Fn() -> usize { 
    let col_count = Mutex::new(0);
    move || {
        let mut c = col_count.lock().unwrap();
        *c += 1;
        *c
    }
}