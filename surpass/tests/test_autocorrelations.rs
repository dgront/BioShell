use rand::Rng;
use bioshell_pdb::calc::Vec3;
use surpass::{AutocorrelateVec3Measurements, SurpassAlphaSystem, SystemMeasurement};

struct Vec3Measurement { last_vec: Vec3 }

impl SystemMeasurement<Vec3> for Vec3Measurement {

    fn measure(&self, system: &SurpassAlphaSystem) -> Vec3 { return system.ca_to_vec3(0) }

    fn header(&self) -> String { unimplemented!() }
}

fn step(step_size: f64, system: &mut SurpassAlphaSystem) {
    let mut v = system.ca_to_vec3(0);
    let mut rng = rand::thread_rng();
    v.x += rng.gen::<f64>() * step_size - step_size/2.0;
    v.y += rng.gen::<f64>() * step_size - step_size/2.0;
    v.z += rng.gen::<f64>() * step_size - step_size/2.0;
    v.normalize();
    system.vec3_to_ca(0, &v);
}


#[test]
fn test_autocorrelations() {
    // --- the system is not used
    let mut s = SurpassAlphaSystem::by_length(&[5], 100.0);
    let m = Vec3Measurement{ last_vec: Default::default() };
    let t_max = 100;
    let n_samples = 100000;
    let mut autocorr = AutocorrelateVec3Measurements::new(m, t_max, "out");
    for i in 0..n_samples {
        autocorr.observe(&s);
        step(0.5, &mut s);
    }

    assert_eq!(autocorr.n_observations(), n_samples);
    let obs = autocorr.observed_values();
    assert_eq!(obs.len(), t_max);
    assert!(obs[0] > 0.9);
    assert!(obs[t_max - 1] < 0.5);
    // for i in 0..t_max {
    //     println!("{} {}", i+1, obs[i]);
    // }
}