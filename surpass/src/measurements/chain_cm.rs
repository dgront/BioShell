use bioshell_pdb::calc::Vec3;
use crate::measurements::SystemMeasurement;
use crate::{calculate_cm, SurpassAlphaSystem};

/// Center-of-Mass (CM) vector measured for a single chain of a [`SurpassAlphaSystem`](SurpassAlphaSystem)
pub struct ChainCM {which_chain: usize}

impl ChainCM {

    /// Creates a new Center-of-Mass measurement for the `which_chain` chain of a system
    pub fn new(which_chain: usize) -> ChainCM { ChainCM {which_chain} }
}

impl SystemMeasurement<Vec3> for ChainCM {

    fn measure(&self, system: &SurpassAlphaSystem) -> Vec3 { calculate_cm(system, self.which_chain) }

    fn header(&self) -> String { String::from("   cm-X     cm-Y    cm-Z") }
}