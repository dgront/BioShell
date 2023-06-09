use crate::{compute_gyration_squared, CartesianSystem};
use bioshell_core::utils::out_writer;
use bioshell_sim::{ Observer};
use std::any::Any;
use std::io::Write;

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
        let coords = object.get_coordinates();
        let mut out_writer = out_writer(&self.out_fname, self.if_append);
        out_writer
            .write(format!("{:.6} ", self.i_model).as_bytes())
            .ok();
        for _ic in 0..coords.get_chains_count() {
            let rg = compute_gyration_squared(&object.get_coordinates());
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