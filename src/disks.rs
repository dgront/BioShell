use std::time::Instant;
use std::ops::Range;
use rand::Rng;

use bioshell_cartesians::{Coordinates, CartesianSystem, coordinates_to_pdb, square_grid_atoms,
                          NbList, ArgonRules};
use bioshell_sim::{Energy, System};

use bioshell_ff::nonbonded::{PairwiseNonbondedEvaluator, SimpleContact};
use bioshell_montecarlo::{Sampler, AcceptanceStatistics, Mover, MCProtocol, MetropolisCriterion,
                          AdaptiveMCProtocol, MoversSet};

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

impl Mover<CartesianSystem> for DiskMover {

    fn perturb(&mut self, system: &mut CartesianSystem) -> Range<usize> {
        let mut rng = rand::thread_rng();
        let i_moved = rng.gen_range(0..system.size());
        system.add(i_moved,rng.gen_range(-self.max_step..self.max_step),
                   rng.gen_range(-self.max_step..self.max_step), 0.0);

        i_moved..i_moved
    }

    fn acceptance_statistics(&self) -> AcceptanceStatistics { self.succ_rate.clone() }

    fn add_success(&mut self) { self.succ_rate.n_succ += 1; }

    fn add_failure(&mut self) { self.succ_rate.n_failed += 1; }

    fn max_range(&self) -> f64 { return self.max_step; }

    fn set_max_range(&mut self, new_val: f64) { self.max_step = new_val; }
}

pub fn main() {
    const R: f64 = 3.0;
    const W: f64 = 1.0;

    const N_SMALL_CYCLES: usize = 100;
    const N_LARGE_CYCLES: usize = 1000;

    let n_atoms: usize = 50 * 50;
    let density: f64 = 0.4;

    // ---------- Create system's coordinates
    let mut coords = Coordinates::new(n_atoms);
    coords.set_box_len(box_width(R,n_atoms, density));
    square_grid_atoms(&mut coords);
    coordinates_to_pdb(&coords, 1, "disks.pdb", false);

    // ---------- Create system's list of neighbors
    let nbl:NbList = NbList::new(R+W as f64,6.0,Box::new(ArgonRules{}));

    // ---------- Create the system
    let mut system: CartesianSystem = CartesianSystem::new(coords, nbl);

    // ---------- Create energy function
    let energy: Box<dyn Energy<CartesianSystem>> = Box::new(PairwiseNonbondedEvaluator::new(R+W,
        Box::new(SimpleContact::new(R,R,R+W,1000.0,-1.0)) ));

    // ---------- Create a sampler and add a mover into it
    let mut simple_sampler: MCProtocol<MetropolisCriterion,CartesianSystem> =
            MCProtocol::new(MetropolisCriterion::new(1.0));
    simple_sampler.add_mover(Box::new(DiskMover::new(3.0)));

    // ---------- Decorate the sampler into an adaptive MC protocol
    let mut sampler = AdaptiveMCProtocol::new(Box::new(simple_sampler));
    sampler.target_rate = 0.4;

    // ---------- Run the simulation!
    let start = Instant::now();
    let mut recent_acceptance = AcceptanceStatistics::default();
    for i in 0..N_LARGE_CYCLES {
        let stats = sampler.get_mover(0).acceptance_statistics();
        sampler.make_sweeps(N_SMALL_CYCLES,&mut system, &energy);
        coordinates_to_pdb(&system.coordinates(), i as i16, "disks.pdb", true);
        let f_succ = stats.recent_success_rate(&recent_acceptance);
        recent_acceptance = stats;
        println!("{} {} {}  {:.2?}", i, energy.energy(&system)/system.size() as f64, f_succ, start.elapsed());
    }
}
