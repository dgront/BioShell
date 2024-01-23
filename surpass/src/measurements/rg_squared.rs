use bioshell_pdb::calc::Vec3;
use crate::measurements::SystemMeasurement;
use crate::{calculate_cm, SurpassAlphaSystem};

pub struct RgSquared {which_chain: usize}

impl RgSquared {
    pub fn new(which_chain: usize) -> RgSquared { RgSquared {which_chain} }
}

impl SystemMeasurement<f64> for RgSquared {
    fn measure(&self, system: &SurpassAlphaSystem) -> f64 {
        let chain_atoms = system.chain_atoms(self.which_chain);
        let cm = calculate_cm(system, self.which_chain);
        let (mut s, mut n, mut cc) = (0.0, 0.0, 0.0);
        let mut atom = Vec3::from_float(0.0);
        for i_atom in chain_atoms.clone() {
            n += 1.0;
            system.set_ca_to_nearest_vec3(i_atom,chain_atoms.start, &mut atom);
            cc = atom.x - cm.x;
            s += cc * cc;
            cc = atom.y - cm.y;
            s += cc * cc;
            cc = atom.z - cm.z;
            s += cc * cc;
        }

        return s/n;
    }

    fn header(&self) -> String { String::from("gyration-radius-squared\n") }
}