use crate::Coordinates;

pub trait Energy {
    fn energy(&self, system: &Coordinates) -> f64;
    fn energy_by_pos(&self, system: &Coordinates, pos:usize) -> f64;
    fn delta_energy_by_pos(&self, old: &Coordinates, pos:usize, new: &Coordinates) -> f64;
}