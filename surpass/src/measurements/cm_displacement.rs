use std::collections::VecDeque;
use bioshell_io::out_writer;
use bioshell_pdb::calc::Vec3;
use crate::{calculate_cm, SurpassAlphaSystem};

pub struct CMDisplacement {
    box_length: f64,
    n_samples: usize,
    t_max: usize,
    displacements: Vec<f64>,
    cm_by_chain: Vec<VecDeque<Vec3>>,
    file_name: String
}

macro_rules! closest_dx_squared {
    ($coord_a: expr, $coord_b: expr, $box_len: expr) => {
        {
            let mut dx = ($coord_a - $coord_b).abs();
            if dx > $box_len/2.0 { dx = $box_len - dx}
            dx*dx
        }
    }
}
impl CMDisplacement {
    pub fn new(box_length: f64, n_chains: usize, t_max: usize, file_name: &str) -> CMDisplacement {
        CMDisplacement{
            box_length, cm_by_chain: vec![VecDeque::with_capacity(t_max); n_chains],
            n_samples: 0, t_max, displacements: vec![0.0; t_max],
            file_name: file_name.to_string(),
        }
    }

    pub fn observe(&mut self, system: &SurpassAlphaSystem) {
        for i_chain in 0..system.count_chains() {
            // ---------- compute the CM vector for the i-th chain
            let cm_vec = calculate_cm(system, i_chain);
            // ---------- if not all the CM vectors were collected: just push the next one to the front
            if self.cm_by_chain[i_chain].len() < self.t_max {
                self.cm_by_chain[i_chain].push_front(cm_vec);
                continue;
            }
            // ---------- if we have enough CM vectors collected - compute displacement between the new observations and all the CMs previously recorded
            for i_time in 0..self.t_max {
                let dx = closest_dx_squared!(cm_vec.x, self.cm_by_chain[i_chain][i_time].x, self.box_length);
                let dy = closest_dx_squared!(cm_vec.y, self.cm_by_chain[i_chain][i_time].y, self.box_length);
                let dz = closest_dx_squared!(cm_vec.z, self.cm_by_chain[i_chain][i_time].z, self.box_length);
                self.displacements[i_time] += dx + dy + dz;
            }
            self.cm_by_chain[i_chain].pop_back();
            self.cm_by_chain[i_chain].push_front(cm_vec);
            self.n_samples += 1;
        }
    }
}

impl Drop for CMDisplacement {
    fn drop(&mut self) {
        let mut stream = out_writer(&self.file_name, false);
        for i in 0..self.t_max {
            stream.write(format!("{:4} {:}\n", i+1, self.displacements[i] / self.n_samples as f64).as_bytes());
        }
    }
}