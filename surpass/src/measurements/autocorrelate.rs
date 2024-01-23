use std::collections::VecDeque;
use bioshell_io::out_writer;
use bioshell_pdb::calc::Vec3;
use crate::measurements::SystemMeasurement;
use crate::SurpassAlphaSystem;

pub struct AutocorrelateVec3Measurements<M: SystemMeasurement<Vec3>> {
    measurement: M,
    i_chain: usize,
    correlations: Vec<f64>,
    observations: VecDeque<Vec3>,
    fname: String,
    n_samples: usize
}

impl<M: SystemMeasurement<Vec3>> AutocorrelateVec3Measurements<M> {

    pub fn new(measurement: M, t_max: usize, file_name: &str) -> AutocorrelateVec3Measurements<M> {
        return AutocorrelateVec3Measurements{
            measurement, i_chain: 0,
            correlations: vec![0.0; t_max],
            observations: Default::default(),
            fname: file_name.to_string(), n_samples: 0
        };
    }

    pub fn observe(&mut self, system: &SurpassAlphaSystem) {
        // --- the number of points to compute product with
        let t_max = self.correlations.len();
        // --- take the measurement of the current conformation
        let v = self.measurement.measure(system);
        if self.n_samples < self.correlations.len() {
            self.observations.push_front(v);
        } else {
            // --- compute products, update statistics
            for i_time in 0..t_max {
                self.correlations[i_time] += Vec3::dot(&self.observations[i_time], &v);
            }
            self.observations.pop_back();
            self.observations.push_front(v);
        }
        self.n_samples += 1;
    }

    pub fn observed_values(&self) -> Vec<f64> {
        self.correlations.iter().map(|vi| vi / self.n_samples as f64).collect()
    }

    pub fn n_observations(&self) -> usize { self.n_samples }

    pub fn write(&mut self) {
        let mut stream = out_writer(&self.fname, false);
        let data = self.observed_values();
        stream.write(format!("# n_samples: {}\n", self.n_samples).as_bytes()).expect("Can't write to a file");
        for i in 0..data.len() {
            stream.write(format!("{:4} {:}\n", i + 1, data[i]/data[0]).as_bytes()).expect("Can't write to a file");
        }
    }
}


impl<M: SystemMeasurement<Vec3>> Drop for AutocorrelateVec3Measurements<M> {
    fn drop(&mut self) { self.write(); }
}