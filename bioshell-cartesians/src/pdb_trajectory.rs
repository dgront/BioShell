use crate::{write_coordinates_to_pdb, CartesianSystem};
use bioshell_sim::{Observer};
use std::any::Any;

/// Observes conformations of a [`CartesianSystem`](CartesianSystem).
///
/// Each [`observe()`](observe()) call records atomic coordinates as a single frame in the PDB format.
pub struct PdbTrajectory {
    pub fname: String,
    pub if_append: bool,
    i_model: usize,
}

impl PdbTrajectory {
    /// Creates a new  [`PdbTrajectory`](PdbTrajectory) observer that will write each frame into a trajectory file.
    ///
    /// # Arguments
    /// * `fname` - name of the output file
    /// * `if_append` - if true, the new frames will be appended to an existing file (if found); otherwise
    ///     the existing file will be wiped off
    pub fn new(fname: String, if_append: bool) -> PdbTrajectory {
        PdbTrajectory {
            fname,
            if_append,
            i_model: 0,
        }
    }
}

impl Observer for PdbTrajectory {
    type S = CartesianSystem;

    fn observe(&mut self, object: &Self::S) {
        write_coordinates_to_pdb(
            &object.get_coordinates(),
            self.i_model as i16,
            &self.fname,
            self.if_append,
        );
        self.i_model += 1;
    }

    fn flush(&mut self) {}

    fn name(&self) -> &str {
        "PdbTrajectory"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}