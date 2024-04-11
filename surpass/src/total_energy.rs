use crate::{MoveProposal, SurpassAlphaSystem, SurpassEnergy};

pub struct TotalEnergy {
    components: Vec<Box<dyn SurpassEnergy>>,
    weights: Vec<f64>
}

impl TotalEnergy {
    pub fn new() -> Self { TotalEnergy{ components: vec![], weights: vec![] } }

    pub fn add_component(&mut self, energy_term: Box<dyn SurpassEnergy>, weight: f64) {
        self.components.push(energy_term);
        self.weights.push(weight);
    }
}

impl SurpassEnergy for TotalEnergy {
    fn evaluate(&self, conf: &SurpassAlphaSystem) -> f64 {
        let mut total = 0.0;
        for (w, en) in self.weights.iter().zip(&self.components) {
            total += w * en.evaluate(conf);
        }
        return total;
    }

    fn evaluate_delta(&self, conf: &SurpassAlphaSystem, move_prop: &MoveProposal) -> f64 {
        let mut total = 0.0;
        for (w, en) in self.weights.iter().zip(&self.components) {
            total += w * en.evaluate_delta(&conf, move_prop);
        }
        return total;
    }
}