use std::ops::Range;

use bioshell_numerical::Vec3;
use crate::Coordinates;

pub trait Energy {
    fn energy(&self, system: &Coordinates) -> f64;
    fn energy_by_pos(&self, system: &Coordinates, pos: usize) -> f64;
    fn delta_energy_by_range(&self, old: &Coordinates, pos: &Range<usize>, new: &Coordinates) -> f64;
}

pub struct ZeroEnergy{}

impl Energy for ZeroEnergy {
    fn energy(&self, system: &Coordinates) -> f64 { 0.0 }

    fn energy_by_pos(&self, system: &Coordinates, pos: usize) -> f64 { 0.0 }

    fn delta_energy_by_range(&self, old: &Coordinates, moved: &Range<usize>, new: &Coordinates) -> f64 { 0.0 }
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

    fn energy(&self, system: &Coordinates) -> f64 {
        let mut total : f64 = 0.0;
        let en_w = self.components.iter().zip(self.weights.iter());
        for (en, w) in en_w {
            total += w * en.energy(system);
        }
        return total;
    }

    fn energy_by_pos(&self, system: &Coordinates, pos: usize) -> f64 {
        let mut total : f64 = 0.0;
        let en_w = self.components.iter().zip(self.weights.iter());
        for (en, w) in en_w {
            total += w * en.energy_by_pos(system, pos);
        }
        return total;
    }

    fn delta_energy_by_range(&self, old: &Coordinates, pos: &Range<usize>, new: &Coordinates) -> f64 {
        let mut total : f64 = 0.0;
        let en_w = self.components.iter().zip(self.weights.iter());
        for (en, w) in en_w {
            total += w * en.delta_energy_by_range(&old, &pos, &new);
        }
        return total;
    }
}
