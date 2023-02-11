use std::any::Any;
use crate::{CartesianSystem, coordinates_to_pdb};
use bioshell_sim::{Energy, System, Observer};


pub struct PdbTrajectory {
    pub fname: String,
    pub if_append: bool,
    i_model: usize
}


impl PdbTrajectory {
    pub fn new(fname: String, if_append: bool) ->PdbTrajectory { PdbTrajectory{ fname, if_append, i_model: 0 } }
}

impl Observer for PdbTrajectory {
    type O = CartesianSystem;

    fn observe(&mut self, object: &Self::O) {
        coordinates_to_pdb(&object.coordinates(),self.i_model as i16,
                           &self.fname, self.if_append);
        self.i_model += 1;
    }

    fn flush(&mut self) {}

    fn name(&self) -> &str { "PdbTrajectory" }

    fn as_any(&self) -> &dyn Any { self }
}