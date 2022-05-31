use std::ops::Range;
use rand::Rng;

use bioshell_ff::{Coordinates, Energy, TotalEnergy, ZeroEnergy};


pub trait Sampler {
    fn run(&mut self, system: &mut Coordinates, n_steps:i32);
}

pub struct IsothermalMC {
    pub temperature: f64,
    pub energy: Box<dyn Energy>,
    movers: MoversSet,
}

impl IsothermalMC {
    pub fn new(temperature:f64) -> IsothermalMC {
        IsothermalMC{temperature, movers: Default::default(), energy: Box::new(ZeroEnergy{}) }
    }

    pub fn add_mover(&mut self, perturb_fn: Box<dyn Fn(&mut Coordinates,f32) -> Range<usize>>, move_range: f32) {
        self.movers.add_mover(perturb_fn,move_range);
    }

}

impl Sampler for IsothermalMC {
    fn run(&mut self, system: &mut Coordinates, n_steps: i32) {

        let mut rng = rand::thread_rng();
        let mut future_system = system.clone();
        for i_step in 0..n_steps {                            // --- iterate over steps
            for i_atom in 0..system.size() {                // --- iteration for moves within a single MC step
                // --- call a mover to modify a pose
                let moved_range = self.movers.make_move(&mut future_system);
                let delta_en: f64 = self.energy.delta_energy_by_range(system, &moved_range, &future_system);
                if delta_en <= 0.0 || rng.gen_range(0.0..1.0) <= (delta_en/self.temperature).exp() {
                    // --- update mover counts, copy future_pose on current_pose to make the move
                } else {
                    // --- update mover failures, copy current_pose on future_pose to clear the move
                }
            }
        }
    }
}


// ----------------- private content -----------------
pub struct AdaptiveMoverStats {
    n_succ:i32,
    n_failed:i32,
    move_range:f32,
    max_move_range:f32
}

impl AdaptiveMoverStats {
    pub fn new(max_move_range:f32) -> AdaptiveMoverStats {
        AdaptiveMoverStats { n_succ: 0, n_failed: 0,
            move_range: max_move_range/2.0, max_move_range }
    }

    pub fn add_success(&mut self) { self.n_succ+=1; }
    pub fn add_failure(&mut self) { self.n_failed+=1; }
    pub fn success_rate(&self) -> f64 {
        self.n_succ as f64 / (self.n_succ as f64 + self.n_failed as f64)
    }
}

#[derive(Default)]
pub struct MoversSet {
    pub stats: Vec<AdaptiveMoverStats>,
    pub proposals: Vec<Box<dyn Fn(&mut Coordinates,f32) -> Range<usize>>>,
}

impl MoversSet {
    pub fn add_mover(&mut self, perturb_fn: Box<dyn Fn(&mut Coordinates,f32) -> Range<usize>>, move_range: f32) {
        self.proposals.push(perturb_fn);
        self.stats.push(AdaptiveMoverStats::new(move_range));
    }

    pub fn make_move(&mut self, future: &mut Coordinates) -> Range<usize> {
        self.proposals[0](future, self.stats[0].move_range)
    }
}