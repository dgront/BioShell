use bioshell_cartesians::CartesianSystem;
use bioshell_core::utils::out_writer;
use bioshell_sim::Observer;

pub struct BondViolationObserver {
    pub out_fname: String,
    pub if_append: bool,
    i_model: usize,
}

impl BondViolationObserver {
    pub fn new(fname: String, if_append: bool) -> BondViolationObserver {
        BondViolationObserver {
            out_fname: fname,
            if_append,
            i_model: 0,
        }
    }
}

impl Observer for BondViolationObserver {
    type S = CartesianSystem;

    fn observe(&mut self, object: &Self::S) {
        let coords = object.get_coordinates();
        let mut out_writer = out_writer(&self.out_fname, self.if_append);
        out_writer
            .write(format!("{:.6} ", self.i_model).as_bytes())
            .ok();
        for _ic in 0..coords.get_chains_count() {
            // HERE !!!
            let result: f64 = 0.12345;
            out_writer.write(format!("{:>10.3} ", result).as_bytes()).ok();
        }
        out_writer.write("\n".as_bytes()).ok();
        self.i_model += 1;
    }

    fn flush(&mut self) {}

    fn name(&self) -> &str {
        "BondViolationObserver"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}