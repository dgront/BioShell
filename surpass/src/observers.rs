use std::collections::VecDeque;
use std::error::Error;
use std::io::Write;
use std::marker::PhantomData;
use bioshell_io::out_writer;
use bioshell_pdb::calc::Vec3;
use crate::{calculate_cm, SurpassAlphaSystem};

pub struct RecordMeasurements<T, M: SystemMeasurement<T>> {
    fname: String,
    measurements: Vec<M>,
    phantom: PhantomData<T>
}

impl<T: std::fmt::Display, M: SystemMeasurement<T>> RecordMeasurements<T, M> {
    pub fn new(fname: &str, measurements: Vec<M>) -> Result<RecordMeasurements<T, M>, Box<dyn Error>> {
        let mut stream = out_writer(&fname, false);
        stream.write(b"#")?;
        for o in &measurements {
            stream.write(o.header().as_bytes())?;
        }
        stream.write(b"\n");
        Ok(RecordMeasurements { fname: fname.to_string(), measurements, phantom: Default::default() })
    }

    pub fn observe(&self, system: &SurpassAlphaSystem) -> Result<(), Box<dyn Error>> {
        let mut stream = out_writer(&self.fname, true);
        for o in &self.measurements {
            let oi = o.measure(system);
            stream.write(format!("{:}", oi).as_bytes())?;
        }
        stream.write(b"\n");
        Ok(())
    }
}

/// Measures a property of a [`SurpassAlphaSystem`](SurpassAlphaSystem)
///
/// At each call of  the [`SystemMeasurement::measure()`](SystemMeasurement::measure()) method,
/// a property of the current conformation of a [`SurpassAlphaSystem`](SurpassAlphaSystem) is
/// evaluated and returned
pub trait SystemMeasurement<T> {
    fn measure(&self, system: &SurpassAlphaSystem) -> T;
    fn header(&self) -> String;
}

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

pub struct REndSquared {which_chain: usize}

impl REndSquared {
    pub fn new(which_chain: usize) -> REndSquared { REndSquared {which_chain} }
}

impl SystemMeasurement<f64> for REndSquared {
    fn measure(&self, system: &SurpassAlphaSystem) -> f64 {
        let chain_atoms = system.chain_residues(self.which_chain);
        system.distance_squared(chain_atoms.start, chain_atoms.end - 1)
    }

    fn header(&self) -> String { String::from("r-end-squared\n") }
}

pub struct RgSquared {which_chain: usize}

impl RgSquared {
    pub fn new(which_chain: usize) -> RgSquared { RgSquared {which_chain} }
}

impl SystemMeasurement<f64> for RgSquared {
    fn measure(&self, system: &SurpassAlphaSystem) -> f64 {
        let chain_atoms = system.chain_residues(self.which_chain);
        let cm = calculate_cm(system, self.which_chain);
        let (mut s, mut n, mut cc) = (0.0, 0.0, 0.0);
        let mut atom = Vec3::from_float(0.0);
        for i_atom in chain_atoms.clone() {
            n += 1.0;
            system.set_atom_to_nearest_vec3(i_atom,chain_atoms.start, &mut atom);
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
