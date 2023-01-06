use std::ops::Range;

use bioshell_sim::{System, Energy};
use bioshell_cartesians::{CartesianSystem};


// pub trait Energy {
//     fn energy(&self, system: &System) -> f64;
//     fn energy_by_pos(&self, system: &System, pos: usize) -> f64;
//     fn delta_energy_by_range(&self, old: &System, pos: &Range<usize>, new: &System) -> f64;
//     /// Returns the name of this energy function
//     /// The returned name may be used to identify this energy, e.g. to name a column in a score table
//     fn name(&self) -> String;
// }

// pub struct ZeroEnergy{}
//
// impl Energy for ZeroEnergy {
//     #[allow(unused)]
//     fn energy(&self, system: &System) -> f64 { 0.0 }
//     #[allow(unused)]
//     fn energy_by_pos(&self, system: &System, pos: usize) -> f64 { 0.0 }
//     #[allow(unused)]
//     fn delta_energy_by_range(&self, old: &System, moved: &Range<usize>, new: &System) -> f64 { 0.0 }
//     #[allow(unused)]
//     fn name(&self) -> String { String::from("ZeroEnergy") }
// }

/// Stores energy components and their weights
/// Allows one calculate the total energy methods for a given Coordinates
#[derive(Default)]
pub struct TotalEnergy<S> {
    pub weights: Vec<f64>,
    pub components: Vec<Box<dyn Energy<S>>>,
}


impl<S> TotalEnergy<S> {

    pub fn add_component(&mut self, en: Box<dyn  Energy<S>>, w: f64) -> &mut TotalEnergy<S> {
        self.weights.push(w);
        self.components.push(en);
        return self;
    }
}

impl<S> Energy<S> for TotalEnergy<S> {

    fn energy(&self, system: &S) -> f64 {
        let mut total : f64 = 0.0;
        let en_w = self.components.iter().zip(self.weights.iter());
        for (en, w) in en_w {
            total += w * en.energy(system);
        }
        return total;
    }

    fn energy_by_pos(&self, system: &S, pos: usize) -> f64 {
        let mut total : f64 = 0.0;
        let en_w = self.components.iter().zip(self.weights.iter());
        for (en, w) in en_w {
            total += w * en.energy_by_pos(system, pos);
        }
        return total;
    }

    fn delta_energy_by_range(&self, old: &S, new: &S, pos: &Range<usize>) -> (f64,f64) {
        let (mut total_old, mut total_new) = (0.0, 0.0);
        let en_w = self.components.iter().zip(self.weights.iter());
        for (en, w) in en_w {
            let (o,n) = en.delta_energy_by_range(&old, &new, &pos);
            total_new += w * n;
            total_old += w * o;
        }
        return (total_old, total_new);
    }

    fn name(&self) -> String { String::from("TotalEnergy") }
}
