
pub trait Distance {
    fn evaluate(&self, pi: &Vec<f64>, pj: &Vec<f64>) -> f64;
}

pub struct Euclidean {}

impl Distance for Euclidean  {
    fn evaluate(&self, pi: &Vec<f64>, pj: &Vec<f64>) -> f64 {
        let mut d: f64 = 0.0;
        for i in 0..pi.len() {
            let t = pi[i] - pj[i];
            d += t * t;
        }
        d.sqrt()
    }
}