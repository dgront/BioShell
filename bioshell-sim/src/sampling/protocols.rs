use std::ops::Range;
use rand::Rng;

use bioshell_ff::{Energy, System, ZeroEnergy};
use crate::sampling::movers::change_volume;


pub trait Sampler {
    fn energy(&self, system: &System) -> f64;
    fn run(&mut self, system: &mut System, n_steps:i32) -> f64;
}

pub enum Ensemle {
    NVT,
    NPT,
}

pub struct IsothermalMC {
    pub temperature: f64,
    pub pressure: f64,
    pub ensemble: Ensemle,
    pub energy: Box<dyn Energy>,
    movers: MoversSet,
}

impl IsothermalMC {
    /// Creates a new isothermal MC sampler
    /// The pressure value will be used only in the case for Ensemble::NPT simulation
    pub fn new(temperature:f64, ensemble:Ensemle, pressure: f64) -> IsothermalMC {
        IsothermalMC{temperature,ensemble, pressure,
            movers: Default::default(), energy: Box::new(ZeroEnergy{}) }
    }

    pub fn add_mover(&mut self, perturb_fn: Box<dyn Fn(&mut System,f32) -> Range<usize>>, move_range: f32) {
        self.movers.add_mover(perturb_fn,move_range);
    }
}

impl Sampler for IsothermalMC {

    fn energy(&self, system: &System) -> f64 { self.energy.energy(system) }

    fn run(&mut self, system: &mut System, n_steps: i32) -> f64 {

        let mut rng = rand::thread_rng();
        let mut future_system = system.clone();
        let mut n_succ: f64 = 0.0;
        for _i_step in 0..n_steps {                            // --- iterate over steps
            for _i_atom in 0..system.size() {                // --- iteration for moves within a single MC step
                // --- call a mover to modify a pose
                let moved_range = self.movers.make_move(&mut future_system);
                let delta_en: f64 = self.energy.delta_energy_by_range(system, &moved_range, &future_system);
                // --- hard rule to skip Metropolis criterion when energy grows to fast - it would be rejected anyways
                #[cfg(debug_assertions)]
                    if delta_en < 1000.0 {
                        // println!("en old:");
                        let en_old = self.energy.energy(&system);
                        // println!("en new:");
                        let en_new = self.energy.energy(&future_system);
                        let total_delta = en_new - en_old;
                        if f64::abs(total_delta-delta_en) > 0.01 {
                            let str = format!("Inconsistent energy! Total {en_old} -> {en_new} with delta = {total_delta}, local delta: {delta_en} after move {:?}",moved_range);
                            panic!("{}", str);
                        }
                    }

                if delta_en <= 0.0 || rng.gen_range(0.0..1.0) <= (-delta_en/self.temperature).exp() {
                    // --- update mover counts, copy future_pose on current_pose to make the move
                    for ipos in moved_range.start..moved_range.end + 1 {
                        system.copy(ipos, &future_system);
                    }
                    self.movers.accepted();
                    n_succ += 1.0;
                } else {
                    // --- update mover failures, copy current_pose on future_pose to clear the move
                    for ipos in moved_range.start..moved_range.end + 1 {
                        future_system.copy(ipos, &system);
                    }
                    self.movers.cancelled();
                }
            } // --- a single MC step done
            // ---------- Here we apply special moves for ensembles other than NVT
            match self.ensemble {
                // ---------- Volume change should be applie in the case of NPT ensemble
                Ensemle::NPT => {
                    const P_CONST : f64 = 1000.0 * 10e-30 / 1.380649e-23;    // converts [Pa/m^3] to [kPa/A^3], divided by Boltzmann constant
                    let moved_range = change_volume(&mut future_system, 0.1);
                    let delta_e = self.energy.energy(&future_system) - self.energy.energy(&system);
                    let pdv = self.pressure*(future_system.volume() - system.volume()) as f64 * P_CONST;
                    let ln_n = self.temperature * system.size() as f64 * ((future_system.volume() / system.volume())as f64).ln();
                    let delta = delta_e + pdv - ln_n;
                    if delta <= 0.0 || rng.gen_range(0.0..1.0) <= (-delta / self.temperature).exp() {
                        println!("changing volume from {} to {}", system.volume(), future_system.volume());
                        println!("deltas {} {} {}", delta_e, pdv, ln_n);
                        system.set_box_len(future_system.box_len());
                    } else {
                        future_system.set_box_len(system.box_len());
                    }
                }
                // ---------- NVT ensemble requires nothing extra
                Ensemle::NVT => {}
            }
        } // --- all inner MC cycles done

        self.movers.adapt_movers();

        n_succ / (n_steps * system.size() as i32) as f64
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
    pub fn success_rate(&self) -> f64 { self.n_succ as f64 / (self.n_succ as f64 + self.n_failed as f64) }

    pub fn adapt_range(&mut self) {
        let rate = self.success_rate();
        if rate < 0.35 { self.move_range *= 0.9 }
        if rate > 0.45 && self.move_range < self.max_move_range { self.move_range *= 1.1 }
    }
}

#[derive(Default)]
pub struct MoversSet {
    pub stats: Vec<AdaptiveMoverStats>,
    pub proposals: Vec<Box<dyn Fn(&mut System,f32) -> Range<usize>>>,
    recent_mover: usize
}

impl MoversSet {
    pub fn add_mover(&mut self, perturb_fn: Box<dyn Fn(&mut System,f32) -> Range<usize>>, move_range: f32) {
        self.proposals.push(perturb_fn);
        self.stats.push(AdaptiveMoverStats::new(move_range));
    }

    pub fn make_move(&mut self, future: &mut System) -> Range<usize> {
        let mut rng = rand::thread_rng();
        self.recent_mover = rng.gen_range(0..(self.stats.len()));
        self.proposals[self.recent_mover](future, self.stats[self.recent_mover].move_range)
    }
    pub fn accepted(&mut self) { self.stats[self.recent_mover].add_success(); }
    pub fn cancelled(&mut self) { self.stats[self.recent_mover].add_failure(); }

    pub fn adapt_movers(&mut self) {
        for i in 0..self.stats.len() {
            self.stats[i].adapt_range();
        }
    }
}