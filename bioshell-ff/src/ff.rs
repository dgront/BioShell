use std::ops::Range;

use bioshell_sim::{Energy};

/// Stores energy components and their weights
/// Allows one calculate the total energy methods for a given Coordinates
pub struct TotalEnergy<S> {
    pub weights: Vec<f64>,
    pub components: Vec<Box<dyn Energy<S>>>,
}


impl<S> TotalEnergy<S> {

    pub fn new() -> TotalEnergy<S> {
        TotalEnergy{weights: vec![], components: vec![] }
    }

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
