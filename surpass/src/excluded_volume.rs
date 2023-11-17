use crate::{NonBondedEnergyKernel, SurpassAlphaSystem};


pub struct ExcludedVolume {
    e_penalty: f64,
    r_rep: f64,
    i_rep_2: f64,
}

impl ExcludedVolume {
    pub fn new(system: &SurpassAlphaSystem, r_rep: f64, e_rep: f64) -> ExcludedVolume {
        let mut i_rep_2 = system.real_to_int(r_rep) as f64;
        i_rep_2 *= i_rep_2;

        ExcludedVolume{ e_penalty: e_rep, r_rep, i_rep_2}
    }

    /// Returns the excluded volume distance that is currently used
    pub fn repulsion_cutoff(&self) -> f64 { self.r_rep }

    /// Returns the excluded volume penalty value that is currently used
    pub fn repulsion_energy(&self) -> f64 { self.e_penalty }

}

impl NonBondedEnergyKernel for ExcludedVolume {

    #[inline(always)]
    fn energy_for_residue_pair(&self, i2: f64) -> f64 {
        if i2 < self.i_rep_2 { self.e_penalty }
        else { 0.0 }
    }

    fn distance_cutoff(&self) -> f64 { self.r_rep }
}