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

macro_rules! spring_kernel {
    ($system:expr, $i:expr, $j:expr, $d0:expr, $en:expr, $OP:tt) => {
        let v:f64 = ($system.closest_distance_square($i, $j).sqrt() - $d0) as f64;
        $en $OP v*v;
    }
}

impl Energy for SimpleHarmonic {

    fn energy(&self, system: &Coordinates) -> f64 {
        let mut en:f64 = 0.0;
        for i in 1..system.size() {
            if system[i].chain_id == system[i-1].chain_id {
                spring_kernel!(system,i, i-1, self.d0, en, +=);
            }
        }

        return en * self.k as f64;
    }

    fn energy_by_pos(&self, system: &Coordinates, pos: usize) -> f64 {
        let mut en: f64 = 0.0;
        if pos > 0 && system[pos].chain_id == system[pos - 1].chain_id {
            spring_kernel!(system,pos, pos-1,self.d0, en, +=);
        }
        if pos < system.size() - 1 && system[pos].chain_id == system[pos + 1].chain_id {
            spring_kernel!(system,pos, pos-1,self.d0, en, +=);
        }

        return en * self.k as f64;
    }

    fn delta_energy_by_range(&self, old: &Coordinates, pos: &Range<usize>, new: &Coordinates) -> f64 {

        let mut en: f64 = 0.0;
        let start: usize = if pos.start > 0 { pos.start - 1 } else { pos.start };
        let end: usize = if pos.end >= old.size() - 1 { old.size() - 1 } else { pos.end + 1 };
        for ipos in start..end {
            if old[ipos].chain_id == old[ipos + 1].chain_id {
                spring_kernel!(new, ipos, ipos+1, self.d0, en, +=);
                spring_kernel!(old, ipos, ipos+1, self.d0, en, -=);
            }
        }

        return en * self.k as f64;
    }
}