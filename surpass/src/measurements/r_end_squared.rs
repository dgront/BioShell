use crate::measurements::SystemMeasurement;
use crate::SurpassAlphaSystem;

pub struct REndSquared {which_chain: usize}

impl REndSquared {
    pub fn new(which_chain: usize) -> REndSquared { REndSquared {which_chain} }
}

impl SystemMeasurement<f64> for REndSquared {
    fn measure(&self, system: &SurpassAlphaSystem) -> f64 {
        let chain_atoms = system.chain_atoms(self.which_chain);
        system.distance_squared(chain_atoms.start, chain_atoms.end - 1)
    }

    fn header(&self) -> String { String::from("r-end-squared\n") }
}