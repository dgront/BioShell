use std::ops::Range;

use crate::System;

pub trait Energy {
    fn energy(&self, system: &System) -> f64;
    fn energy_by_pos(&self, system: &System, pos: usize) -> f64;
    fn delta_energy_by_range(&self, old: &System, pos: &Range<usize>, new: &System) -> f64;
    /// Returns the name of this energy function
    /// The returned name may be used to identify this energy, e.g. to name a column in a score table
    fn name(&self) -> String;
}

pub struct ZeroEnergy{}

impl Energy for ZeroEnergy {
    #[allow(unused)]
    fn energy(&self, system: &System) -> f64 { 0.0 }
    #[allow(unused)]
    fn energy_by_pos(&self, system: &System, pos: usize) -> f64 { 0.0 }
    #[allow(unused)]
    fn delta_energy_by_range(&self, old: &System, moved: &Range<usize>, new: &System) -> f64 { 0.0 }
    #[allow(unused)]
    fn name(&self) -> String { String::from("ZeroEnergy") }
}

/// Stores energy components and their weights
/// Allows one calculate the total energy methods for a given Coordinates
#[derive(Default)]
pub struct TotalEnergy {
    pub weights: Vec<f64>,
    pub components: Vec<Box<dyn Energy>>,
}


impl TotalEnergy {

    pub fn add_component(&mut self, en: Box<dyn  Energy>, w: f64) -> &mut TotalEnergy {
        self.weights.push(w);
        self.components.push(en);
        return self;
    }
}

impl Energy for TotalEnergy {

    fn energy(&self, system: &System) -> f64 {
        let mut total : f64 = 0.0;
        let en_w = self.components.iter().zip(self.weights.iter());
        for (en, w) in en_w {
            total += w * en.energy(system);
        }
        return total;
    }

    fn energy_by_pos(&self, system: &System, pos: usize) -> f64 {
        let mut total : f64 = 0.0;
        let en_w = self.components.iter().zip(self.weights.iter());
        for (en, w) in en_w {
            total += w * en.energy_by_pos(system, pos);
        }
        return total;
    }

    fn delta_energy_by_range(&self, old: &System, pos: &Range<usize>, new: &System) -> f64 {
        let mut total : f64 = 0.0;
        let en_w = self.components.iter().zip(self.weights.iter());
        for (en, w) in en_w {
            total += w * en.delta_energy_by_range(&old, &pos, &new);
        }
        return total;
    }

    fn name(&self) -> String { String::from("TotalEnergy") }
}
