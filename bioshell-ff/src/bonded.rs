use std::ops::Range;

use bioshell_cartesians::CartesianSystem;
use bioshell_sim::{Energy, System};

pub struct SimpleHarmonic {
    /// equilibrium distance
    d0: f64,
    /// force constant
    k: f64,
}
impl SimpleHarmonic {
    pub fn new(d0: f64, k: f64) -> SimpleHarmonic {
        SimpleHarmonic { d0, k }
    }
}

macro_rules! spring_kernel {
    ($system:expr, $i:expr, $j:expr, $d0:expr, $en:expr, $OP:tt) => {
        let v:f64 = ($system.closest_distance_square($i, $j).sqrt() - $d0) as f64;
        $en $OP v*v;
    }
}

impl Energy<CartesianSystem> for SimpleHarmonic {
    fn energy(&self, system: &CartesianSystem) -> f64 {
        let mut en: f64 = 0.0;
        for i in 1..system.size() {
            if system.coordinates()[i].chain_id == system.coordinates()[i - 1].chain_id {
                spring_kernel!(system.coordinates(),i, i-1, self.d0, en, +=);
            }
        }

        return en * self.k as f64;
    }

    fn energy_by_pos(&self, system: &CartesianSystem, pos: usize) -> f64 {
        let mut en: f64 = 0.0;
        if pos > 0 && system.coordinates()[pos].chain_id == system.coordinates()[pos - 1].chain_id {
            spring_kernel!(system.coordinates(),pos, pos-1, self.d0, en, +=);
        }
        if pos < system.size() - 1
            && system.coordinates()[pos].chain_id == system.coordinates()[pos + 1].chain_id
        {
            spring_kernel!(system.coordinates(),pos, pos+1, self.d0, en, +=);
        }

        return en * self.k as f64;
    }

    fn name(&self) -> String {
        String::from("SimpleHarmonic")
    }

    fn energy_by_range(&self, system: &CartesianSystem, range: &Range<usize>) -> f64 {
        let mut total_en: f64 = 0.0;
        let start: usize = if range.start > 0 {
            range.start - 1
        } else {
            range.start
        };
        let end: usize = if range.end >= system.size() - 1 {
            system.size() - 1
        } else {
            range.end + 1
        };
        for ipos in start..end {
            if system.coordinates()[ipos].chain_id == system.coordinates()[ipos + 1].chain_id {
                spring_kernel!(system.coordinates(), ipos, ipos+1, self.d0, total_en, +=);
            }
        }

        return total_en * self.k;
    }
}
