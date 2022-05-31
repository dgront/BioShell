use std::ops::Range;

use crate::ff::Energy;
use crate::Coordinates;

pub struct SimpleHarmonic {
    /// equilibrium distance
    d0: f32,
    /// force constant
    k: f32
}
impl SimpleHarmonic {
    pub fn new(d0: f32, k: f32) -> SimpleHarmonic {
        SimpleHarmonic { d0, k}
    }
}

impl Energy for SimpleHarmonic {

    fn energy(&self, system: &Coordinates) -> f64 {
        let mut en:f64 = 0.0;
        for i in 1..system.size() {
            if system[i].chain_id == system[i-1].chain_id {
                let mut v:f64 = (system.distance_square(i, i - 1).sqrt() - self.d0) as f64;
                en += v*v;
            }
        }

        return en * self.k as f64;
    }

    fn energy_by_pos(&self, system: &Coordinates, pos: usize) -> f64 {
        let mut en: f64 = 0.0;
        if pos > 0 && system[pos].chain_id == system[pos - 1].chain_id {
            en += (system.distance_square(pos, pos - 1).sqrt() - self.d0) as f64;
        }
        if pos < system.size() - 1 && system[pos].chain_id == system[pos + 1].chain_id {
            en += (system.distance_square(pos, pos + 1).sqrt() - self.d0) as f64;
        }

        return en * self.k as f64;
    }

    fn delta_energy_by_range(&self, old: &Coordinates, pos: &Range<usize>, new: &Coordinates) -> f64 {

        return 0.0;
    }
}