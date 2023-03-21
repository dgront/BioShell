use crate::{coordinates_to_pdb, gyration_squared, r_end_squared, CartesianSystem};
use bioshell_core::utils::out_writer;
use bioshell_sim::{Energy, Observer, System};
use std::any::Any;
use std::io::Write;

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
        coordinates_to_pdb(
            &object.coordinates(),
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

pub struct GyrationSquared {
    pub out_fname: String,
    pub if_append: bool,
    i_model: usize,
}

impl GyrationSquared {
    pub fn new(fname: String, if_append: bool) -> GyrationSquared {
        GyrationSquared {
            out_fname: fname,
            if_append,
            i_model: 0,
        }
    }
}

impl Observer for GyrationSquared {
    type S = CartesianSystem;

    fn observe(&mut self, object: &Self::S) {
        let coords = object.coordinates();
        let mut out_writer = out_writer(&self.out_fname, self.if_append);
        out_writer
            .write(format!("{:.6} ", self.i_model).as_bytes())
            .ok();
        for ic in 0..coords.count_chains() {
            let rg = gyration_squared(&object.coordinates(), ic);
            out_writer.write(format!("{:>10.3} ", rg).as_bytes()).ok();
        }
        out_writer.write("\n".as_bytes()).ok();
        self.i_model += 1;
    }

    fn flush(&mut self) {}

    fn name(&self) -> &str {
        "GyrationSquared"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

pub struct REndSquared {
    pub out_fname: String,
    pub if_append: bool,
    i_model: usize,
}

impl REndSquared {
    pub fn new(fname: String, if_append: bool) -> REndSquared {
        REndSquared {
            out_fname: fname,
            if_append,
            i_model: 0,
        }
    }
}

impl Observer for REndSquared {
    type S = CartesianSystem;

    fn observe(&mut self, object: &Self::S) {
        let coords = object.coordinates();
        let mut out_writer = out_writer(&self.out_fname, self.if_append);
        out_writer
            .write(format!("{:.6} ", self.i_model).as_bytes())
            .ok();
        for ic in 0..coords.count_chains() {
            let rg = r_end_squared(&object.coordinates(), ic);
            out_writer.write(format!("{:>10.3} ", rg).as_bytes()).ok();
        }
        out_writer.write("\n".as_bytes()).ok();
        self.i_model += 1;
    }

    fn flush(&mut self) {}

    fn name(&self) -> &str {
        "REndSquared"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}
