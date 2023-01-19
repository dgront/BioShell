use std::env;
use std::time::Instant;
use std::ops::Range;
use rand::Rng;

use log::{info};

use clap::{Parser};

use bioshell_cartesians::{Coordinates, CartesianSystem, coordinates_to_pdb, square_grid_atoms,
                          NbList, ArgonRules, pdb_to_coordinates};
use bioshell_sim::{Energy, System};

use bioshell_ff::nonbonded::{PairwiseNonbonded, PairwiseNonbondedEvaluator, SimpleContact};
use bioshell_montecarlo::{Sampler, AcceptanceStatistics, Mover, MCProtocol, MetropolisCriterion, AdaptiveMCProtocol, MoversSet, AcceptanceCriterion};

#[derive(Parser, Debug)]
#[clap(name = "disks")]
#[clap(version = "0.2")]
#[clap(about = "Simple MC simulation of disks in 2D", long_about = None)]
#[clap(after_help = "Example:\n\tdisks -d 0.4 -n 400 -i 1000 -o 1000 -t 0.42")]
struct Args {
    /// staring conformation in the PDB format
    #[clap(short, long, default_value = "", short='f')]
    infile: String,
    /// temperature of the simulation
    #[clap(short, long, default_value_t = 1.0, short='t')]
    temperature: f64,
    /// density of the simulated system
    #[clap(short, long, default_value_t = 0.4, short='d')]
    density: f64,
    /// Number of disks in a box to simulate
    #[clap(short, long, default_value_t = 400, short='n')]
    ndisks: usize,
    /// the number of inner MC cycles
    #[clap(short, long, default_value_t = 100)]
    inner: usize,
    /// the number of outer MC cycles
    #[clap(short, long, default_value_t = 100)]
    outer: usize,
    /// prefix for output file names
    #[clap(long, default_value = "")]
    prefix: String,
    /// don't run any simulation, just write energy value for each disk of the starting conformation
    #[clap(long)]
    rescore: bool,
    /// thickness of the buffer zone used by non-bonded list to hash interactions
    #[clap(long, default_value_t = 4.0)]
    buffer: f64,
}

fn box_width(disc_radius: f64, n_discs: usize, density: f64) -> f64 {

    let v: f64 = std::f64::consts::PI * disc_radius.powi(2);
    (n_discs as f64 * v/density).powf(0.5)
}

/// Moves a random disc
struct DiskMover {
    max_step: f64,
    succ_rate: AcceptanceStatistics
}

impl DiskMover {
    pub fn new(max_range: f64) -> DiskMover {
        DiskMover{ max_step: max_range, succ_rate: Default::default() }
    }
}

impl Mover<CartesianSystem, PairwiseNonbondedEvaluator<SimpleContact>, MetropolisCriterion> for DiskMover {

    fn perturb(&mut self, system: &mut CartesianSystem,
               energy: &PairwiseNonbondedEvaluator<SimpleContact>, acc: &mut MetropolisCriterion) -> Option<Range<usize>> {

        let mut rng = rand::thread_rng();
        // ---------- index of a disk we move
        let i_moved = rng.gen_range(0..system.size());
        // ---------- backup its coordinates and record current energy
        let old_x: f64 = system.coordinates()[i_moved].x;
        let old_y: f64 = system.coordinates()[i_moved].y;
        let old_en: f64 = energy.energy_by_pos(system, i_moved);

        // ---------- test the energy consistency - in debug build only
        #[cfg(debug_assertions)]
            let total_en_before: f64 = energy.energy(system);

        // ---------- actually make the move and update the NBL (mandatory!!)
        system.add(i_moved,rng.gen_range(-self.max_step..self.max_step),
                   rng.gen_range(-self.max_step..self.max_step), 0.0);
        system.update_nbl(i_moved);

        // ---------- calculate the new energy
        let new_en: f64 = energy.energy_by_pos(system, i_moved);

        // ---------- check the MC acceptance criterion
        if acc.check(old_en, new_en) {
            system.update_nbl(i_moved);

            // ---------- energy consistency test - part two
            #[cfg(debug_assertions)] {
                let total_en_after: f64 = energy.energy(system);
                let total_delta: f64 = total_en_after - total_en_before;
                let local_delta = new_en - old_en;
                if f64::abs(total_delta - local_delta) > 0.01 {

                    let str = format!("Inconsistent energy! Total {total_en_before} -> \
                        {total_en_after} with delta = {total_delta}, local delta: {local_delta} \
                        after moveing disk {}", i_moved);
                    panic!("{}", str);
                }
            }

            self.succ_rate.n_succ += 1;
            return Option::from(i_moved..i_moved);
        } else {
            // ---------- undo the move
            self.succ_rate.n_failed += 1;
            system.set(i_moved, old_x, old_y, 0.0);
            system.update_nbl(i_moved);
            return Option::None;
        }
    }

    fn acceptance_statistics(&self) -> AcceptanceStatistics { self.succ_rate.clone() }

    fn max_range(&self) -> f64 { return self.max_step; }

    fn set_max_range(&mut self, new_val: f64) { self.max_step = new_val; }
}

fn rescore(conformation: &CartesianSystem, energy: &PairwiseNonbondedEvaluator<SimpleContact>) {
    let mut total: f64 = 0.0;
    for pos in 0..conformation.size() {
        let en: f64 = energy.energy_by_pos(conformation, pos);
        total += en;
        println!("{} {}", pos, en);
    }
    println!("total, by-pos: {} {}", energy.energy(conformation), total/2.0);
}

pub fn main() {
    const R: f64 = 4.0;
    const W: f64 = 1.0;
    const MAX_MOVE_RANGE: f64 = 1.0;

    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    let args = Args::parse();
    // ---------- Temperature of the isothermal simulation (in the energy units)
    let temperature: f64 = args.temperature;
    let density: f64 = args.density;
    let prefix = args.prefix;
    let tra_fname = format!("{}_tra.pdb", &prefix);
    let final_fname = format!("{}_final.pdb", &prefix);
    let n_atoms: usize;
    let box_length: f64;
    let buffer_thickness = args.buffer.max(MAX_MOVE_RANGE);

    // ---------- Create system's coordinates
    let mut coords: Coordinates;
    if args.infile == "" {
        n_atoms = args.ndisks;
        box_length = box_width(R, n_atoms, density);
        coords = Coordinates::new(n_atoms);
        coords.set_box_len(box_length);
        square_grid_atoms(&mut coords);
    } else {
        let res = pdb_to_coordinates(&args.infile).ok();
        match res {
            Some(x) => { coords = x; },
            _ => panic!("Can't read from file >{}<", args.infile)
        };
        n_atoms = coords.size();
        box_length = box_width(R, n_atoms, density);
        coords.set_box_len(box_length);
    }
    coordinates_to_pdb(&coords, 1, tra_fname.as_str(), false);

    info!("Simulation settings:
        temperature: {:.3}
        density:     {}
        no. disks:   {}
        box length:  {:.3}
        NBL buffer:  {:.3}", temperature, density, n_atoms, box_length, buffer_thickness);

    // ---------- Create system's list of neighbors
    let nbl:NbList = NbList::new(R+W as f64,buffer_thickness,Box::new(ArgonRules{}));

    // ---------- Create the system
    let mut system: CartesianSystem = CartesianSystem::new(coords, nbl);

    // ---------- Create energy function
    let kernel = SimpleContact::new(R,R,R+W,100000.0,-1.0);
    let pairwise = PairwiseNonbondedEvaluator::new(R+W,kernel);

    // ---------- Just rescore, no simulation this time
    if args.rescore {
        rescore(&system, &pairwise);
        return;
    }

    // ---------- Create a sampler and add a mover into it
    let mut simple_sampler: MCProtocol<CartesianSystem, PairwiseNonbondedEvaluator<SimpleContact>, MetropolisCriterion> =
            MCProtocol::new(MetropolisCriterion::new(temperature));
    simple_sampler.add_mover(Box::new(DiskMover::new(MAX_MOVE_RANGE)));

    // ---------- Decorate the sampler into an adaptive MC protocol
    let mut sampler = AdaptiveMCProtocol::new(Box::new(simple_sampler));
    sampler.target_rate = 0.4;

    // ---------- Run the simulation!
    let start = Instant::now();
    let mut recent_acceptance = AcceptanceStatistics::default();
    for i in 0..args.outer {
        let stats = sampler.get_mover(0).acceptance_statistics();
        let max_move_range = sampler.get_mover(0).max_range();
        sampler.make_sweeps(args.inner,&mut system, &pairwise);
        coordinates_to_pdb(&system.coordinates(), i as i16, tra_fname.as_str(), true);
        let f_succ = stats.recent_success_rate(&recent_acceptance);
        recent_acceptance = stats;
        println!("{} {} {:.3} {:.3}  {:.2?}", i, pairwise.energy(&system)/system.size() as f64, f_succ,
                 max_move_range, start.elapsed());
    }
    coordinates_to_pdb(&system.coordinates(),1,final_fname.as_str(), false);
}
