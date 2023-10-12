use log::debug;
use crate::{NonBondedEnergyKernel, SurpassAlphaSystem};

pub struct CaContactEnergy {
    e_cont: f64,
    e_rep: f64,
    r_max: f64,
    r_min: f64,
    r_rep: f64,
    i_rep_2: f64,
    i_min_2: f64,
    i_max_2: f64,
}

impl CaContactEnergy {
    pub fn new(system: &SurpassAlphaSystem, e_rep: f64, e_cont: f64, r_rep: f64, r_min: f64, r_max: f64) -> CaContactEnergy {
        let mut i_rep_2 = system.real_to_int(r_rep) as f64;
        i_rep_2 *= i_rep_2;
        let mut i_min_2 = system.real_to_int(r_min) as f64;
        i_min_2 *= i_min_2;
        let mut i_max_2 = system.real_to_int(r_max) as f64;
        i_max_2 *= i_max_2;

        CaContactEnergy{ e_cont, e_rep, r_max, r_min, r_rep, i_rep_2, i_min_2, i_max_2 }
    }
}

impl NonBondedEnergyKernel for CaContactEnergy {

    #[inline(always)]
    fn energy_for_distance_squared(&self, i2: f64) -> f64 {
        if i2 >= self.i_max_2 { 0.0 }
        else {
            if i2 < self.i_rep_2 { self.e_rep}
            else {
                if i2 > self.i_min_2 { self.e_cont }
                else { 0.0 }
            }
        }
    }

    fn distance_cutoff(&self) -> f64 { self.r_max }
}
