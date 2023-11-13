use std::hash::Hash;
use yui::Elem;

pub trait OrdForDisplay {
    fn cmp_for_display(&self, other: &Self) -> std::cmp::Ordering;
}

impl<T> OrdForDisplay for T where T: Ord { 
    fn cmp_for_display(&self, other: &Self) -> std::cmp::Ordering {
        std::cmp::Ord::cmp(self, other)
    }
}

pub trait SortForDisplay { 
    fn sort_for_display(&mut self);
}

impl<T> SortForDisplay for Vec<T> where T: OrdForDisplay { 
    fn sort_for_display(&mut self) {
        self.sort_by(T::cmp_for_display)
    }
}

pub trait Gen: Elem + Hash + OrdForDisplay {}