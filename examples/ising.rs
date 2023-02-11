use std::fmt;
use std::ops::Range;
use std::time::Instant;
use rand::{Rng, SeedableRng};
use rand::rngs::{SmallRng};

use clap::{Parser};

use bioshell_sim::{Energy, ObserversSet, System};
use bioshell_montecarlo::{Sampler, AcceptanceStatistics, IsothermalMC, Mover, AcceptanceCriterion};


#[derive(Parser, Debug)]
#[clap(name = "ising")]
#[clap(about = "Simple MC simulation of an Ising spin system", long_about = None)]
struct Args {
    /// loads a staring conformation from a text file
    #[clap(short, long, default_value = "", short='f')]
    infile: String,
    /// temperature of the simulation
    #[clap(short, long, default_value_t = 1.70)]
    temperature: f64,
    /// Inner cycles
    #[clap(short, long, default_value_t = 100)]
    inner: usize,
    /// outer cycles
    #[clap(short, long, default_value_t = 100)]
    outer: usize,
    /// prefix for output file names
    #[clap(long, default_value = "")]
    prefix: String,
}

const SIZE: usize = 200;

#[derive(Clone, Debug)]
struct IsingSystem { state: [[i8; SIZE]; SIZE] }

impl IsingSystem {

    pub fn new() -> IsingSystem { IsingSystem{state:  [[0; SIZE]; SIZE]} }

    pub fn randomise(&mut self) {
        let mut rng: SmallRng = SmallRng::from_entropy();
        for i in 0..SIZE {
            for j in 0..SIZE {
                let s: i8 = rng.gen();
                self.state[i][j] = if s < 0 { -1 } else { 1 };
            }
        }
    }
}

impl System for IsingSystem {
    fn size(&self) -> usize { SIZE*SIZE }

    fn copy_from(&mut self, pos: usize, rhs: &Self) {
        let i = pos / SIZE;
        let j = pos % SIZE;

        self.state[i][j] = rhs.state[i][j];
    }
}

impl fmt::Display for IsingSystem {

    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let mut out: String = String::new();
        for i in 0..SIZE {
            for j in 0..SIZE {
                match self.state[i][j] {
                    -1 => {out += "_"},
                    1 => {out += "#"}
                    _ => {out += " "}
                }
            }
            out += "\n";
        }

        write!(f, "{}", out)
    }
}

struct IsingEnergy;

impl IsingEnergy {

    pub fn energy_by_ij(s: &IsingSystem, i: usize, j: usize) -> i32 {
        let top = match i {
            0 => s.state[SIZE - 1][j],
            _ => s.state[i-1][j]
        };
        let bottom = match i + 1 {
            SIZE => s.state[0][j],
            _    => s.state[i+1][j]
        };
        let left = match j {
            0 => s.state[i][SIZE-1],
            _ => s.state[i][j-1]
        };
        let right = match j + 1 {
            SIZE => s.state[i][0],
            _    => s.state[i][j+1]
        };
        (-s.state[i][j] * (top + bottom + left + right)) as i32
    }
}

impl Energy<IsingSystem> for IsingEnergy {
    fn energy(&self, system: &IsingSystem) -> f64 {
        let mut total: f64 = 0.0;
        for i in 0..SIZE {
            for j in 0..SIZE {
                total += IsingEnergy::energy_by_ij(system, i, j) as f64;
            }
        }
        return total;
    }

    fn energy_by_pos(&self, system: &IsingSystem, pos: usize) -> f64 {
        let i = pos / SIZE;
        let j = pos % SIZE;
        return IsingEnergy::energy_by_ij(system, i, j) as f64;
    }

    fn energy_by_range(&self, system: &IsingSystem, range: &Range<usize>) -> f64 {
        let mut total: f64 = 0.0;
        for i in range.start..(range.end+1) {
                total += self.energy_by_pos(system, i);
        }
        return total;
    }

    fn name(&self) -> String { String::from("IsingEnergy") }
}

/// Flips a random spin
struct FlipMover {
    succ_rate: AcceptanceStatistics,
    rng: SmallRng,
}

impl FlipMover {
    pub fn new() -> FlipMover {
        FlipMover{ succ_rate: Default::default(), rng: SmallRng::from_entropy() }
    }
}

impl Mover<IsingSystem, IsingEnergy> for FlipMover {

    fn perturb(&mut self, system: &mut IsingSystem, energy: &IsingEnergy, acc: &mut dyn AcceptanceCriterion) -> Option<Range<usize>> {
        let i: usize = self.rng.gen_range(0..SIZE);
        let j: usize = self.rng.gen_range(0..SIZE);
        let delta_e: f64 = -2.0 * IsingEnergy::energy_by_ij(system, i, j) as f64;
        if acc.check(0.0, delta_e) {
            system.state[i][j] *= -1;
            self.succ_rate.n_succ += 1;
            return Some(i..j)
        }
        self.succ_rate.n_failed += 1;
        return None;
    }

    fn acceptance_statistics(&self) -> AcceptanceStatistics { self.succ_rate.clone() }

    fn max_range(&self) -> f64 { 1.0 }

    fn set_max_range(&mut self, _new_val: f64) {}
}

pub fn main() {

    let args = Args::parse();
    let temperature: f64 = args.temperature;    // --- Temperature of the isothermal simulation (in the energy units)
    let prefix = args.prefix;
    let tra_fname = format!("{}_tra.txt", &prefix);
    let final_fname = format!("{}_final.txt", &prefix);

    // ---------- Create system's coordinates
    let mut ising: IsingSystem = IsingSystem::new();
    ising.randomise();

    // ---------- Create the Ising energy function
    let en = IsingEnergy{};

    // ---------- Create a sampler and add a mover into it
    let mut simple_sampler: IsothermalMC<IsingSystem, IsingEnergy> = IsothermalMC::new(temperature);
    simple_sampler.add_mover(Box::new(FlipMover::new()));

    // ---------- Run the simulation!
    let mut observations: ObserversSet<IsingSystem> = ObserversSet::new();
    simple_sampler.run_simulation(args.inner, args.outer, &mut ising, &en, &mut observations);
}
