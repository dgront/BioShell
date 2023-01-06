use std::ops::Range;

pub trait Energy<S> {
    fn energy(&self, system: &S) -> f64;
    fn energy_by_pos(&self, system: &S, pos: usize) -> f64;
    fn delta_energy_by_range(&self, old_system: &S, new_system: &S, pos: &Range<usize>) -> (f64, f64);
    /// Returns the name of this energy function
    /// The returned name may be used to identify this energy, e.g. to name a column in a score table
    fn name(&self) -> String;
}

pub trait System: Clone {
    fn size(&self) -> usize;
    fn copy_from(&mut self, i:usize, rhs: &Self);
}