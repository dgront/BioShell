use bioshell_cartesians::CartesianSystem;
use bioshell_core::utils::out_writer;
use bioshell_sim::Observer;

pub struct BondViolationObserver {
    pub out_fname: String,
    pub if_append: bool,
    i_model: usize,
    longest_bond_length: f64,
}

impl BondViolationObserver {
    pub fn new(fname: String, if_append: bool) -> BondViolationObserver {
        BondViolationObserver {
            out_fname: fname,
            if_append,
            i_model: 0,
            longest_bond_length: 0.0,
        }
    }
}

impl Observer for BondViolationObserver {
    type S = CartesianSystem;

    fn observe(&mut self, object: &Self::S) {
        let coords = object.get_coordinates();
        let mut out_writer = out_writer(&self.out_fname, self.if_append);
        out_writer.write(format!("{:.6} ", self.i_model).as_bytes()).ok();

        let num_atoms = coords.get_size();
        let mut longest_bond_length = 0.0;

        for i in 0..num_atoms - 1 {
            let j = i + 1;
            let distance_squared = coords.get_closest_distance_square(i, j);
            let distance = distance_squared.sqrt();

            if distance > longest_bond_length {
                longest_bond_length = distance;
            }
        }

        out_writer
            .write(format!("{:>10.3} ", longest_bond_length).as_bytes())
            .ok();

        out_writer.write("\n".as_bytes()).ok();
        self.i_model += 1;
        self.longest_bond_length = longest_bond_length;
    }

    fn flush(&mut self) {}

    fn name(&self) -> &str {
        "BondViolationObserver"
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }
}
