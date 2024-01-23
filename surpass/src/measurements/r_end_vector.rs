use bioshell_pdb::calc::Vec3;
use crate::measurements::SystemMeasurement;
use crate::SurpassAlphaSystem;

/// End-to-end vector measured for a single chain of a [`SurpassAlphaSystem`](SurpassAlphaSystem)
pub struct REndVector { which_chain: usize }

impl REndVector {

    /// Creates a new end-to-end vector measurement for the `which_chain` chain of a system
    pub fn new(which_chain: usize) -> REndVector { REndVector { which_chain } }
}

impl SystemMeasurement<Vec3> for REndVector {

    fn measure(&self, system: &SurpassAlphaSystem) -> Vec3 {
        let chain_atoms = system.chain_atoms(self.which_chain);
        let v_start = system.ca_to_vec3(chain_atoms.start);
        let mut v_end = system.ca_to_nearest_vec3( chain_atoms.end - 1, chain_atoms.start);
        v_end -= &v_start;
        return v_end;
    }

    fn header(&self) -> String { String::from("   r_end-X     r_end-Y    r_end-Z") }
}