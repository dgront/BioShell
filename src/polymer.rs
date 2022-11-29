use std::time::Instant;
use std::ops::Range;

use clap::{Parser};

use bioshell_numerical::Vec3;
use bioshell_core::structure::Coordinates;
use bioshell_core::calc::structure::{gyration_squared, r_end_squared};
use bioshell_core::io::pdb::{coordinates_to_pdb, pdb_to_coordinates};

use bioshell_ff::{Energy, TotalEnergy, System};
use bioshell_ff::bonded::SimpleHarmonic;
use bioshell_ff::nonbonded::{SimpleContact, PairwiseNonbondedEvaluator, NbList, PolymerRules};
use bioshell_sim::generators::random_chain;
use bioshell_sim::sampling::movers::{single_atom_move};
use bioshell_sim::sampling::protocols::{Ensemle, IsothermalMC, Sampler};

#[derive(Parser, Debug)]
#[clap(name = "polymer")]
#[clap(version = "0.2")]
#[clap(about = "Simple polymer model", long_about = None)]
/// simulation of a simple polymer model
/// say polymer -h to see options
struct Args {
    /// staring conformation in the PDB format
    #[clap(short, long, default_value = "", short='f')]
    infile: String,
    /// temperature of the simulation
    #[clap(short, long, default_value_t = 1.0, short='t')]
    temperature: f64,
    /// Number of beads of a polymer chain
    #[clap(short, long, default_value_t = 100, short='n')]
    nbeads: usize,
    /// Inner cycles
    #[clap(short, long, default_value_t = 100)]
    inner: i32,
    /// outer cycles
    #[clap(short, long, default_value_t = 100)]
    outer: i32,
    /// prefix for output file names
    #[clap(long, default_value = "")]
    prefix: String,
}

pub fn main() {
    const E_REP: f64 = 4.25;
    const E_FROM: f64 = 4.5;
    const E_TO: f64 = 6.0;
    const E_VAL: f64 = -1.0;
    const REP_VAL: f64 = 1000.0;
    const L: f64 = 900.0;

    let args = Args::parse();

    let temperature: f64 = args.temperature;
    // ---------- Create system's coordinates
    let mut n_beads: usize = 0;
    let mut coords = Coordinates::new(n_beads);

    if args.infile == "" {
        n_beads = args.nbeads;
        coords = Coordinates::new(n_beads);
        coords.set_box_len(L);
        let start: Vec3 = Vec3::new(L/2.0, L/2.0, L/2.0);
        random_chain(3.8, E_FROM as f64, &start, &mut coords);
    } else {
        let res = pdb_to_coordinates(&args.infile).ok();
        match res {
            Some(x) => {coords = x},
            _ => panic!("Can't read from file >{}<", args.infile)
        };
        coords.set_box_len(L);
    }

    // ---------- Create system's list of neighbors
    let nbl:NbList = NbList::new(E_TO,4.0,Box::new(PolymerRules{}));
    // ---------- Create the system
    let mut system: System = System::new(coords, nbl);

    let prefix = args.prefix;
    let tra_fname = format!("{}tra.pdb", &prefix);
    coordinates_to_pdb(&system.coordinates(), 0, tra_fname.as_str(), true);

    // ---------- Contact energy
    let contacts = PairwiseNonbondedEvaluator::new(E_TO as f64,
            Box::new(SimpleContact::new(E_REP,E_FROM,E_TO,REP_VAL,E_VAL)));
    // ---------- Harmonic energy (i.e. springs between beads)
    let harmonic = SimpleHarmonic::new(3.8,5.0);

    let mut total = TotalEnergy::default();
    total.add_component(Box::new(harmonic), 1.0);
    total.add_component(Box::new(contacts), 1.0);
    println!("# starting energy: {}", total.energy(&system));

    let mut sampler = IsothermalMC::new(temperature, Ensemle::NVT, 1.0);
    sampler.energy = Box::new(total);    // --- The total has been moved to a box within the sampler
    let m: Box<dyn Fn(&mut System,f64) -> Range<usize>> = Box::new(single_atom_move);
    sampler.add_mover(m,3.0);
    // let m: Box<dyn Fn(&mut Coordinates,f64) -> Range<usize>> = Box::new(perturb_chain_fragment);
    // sampler.add_mover(m,5.0);

    let start = Instant::now();
    println!("#cycle   energy  f_acc   r_end_sq     rg_sq   time");
    for i in 0..args.outer {
        let f_succ = sampler.run(&mut system, args.inner);
        coordinates_to_pdb(&system.coordinates(), (i+1) as i16, tra_fname.as_str(), true);
        println!("{:6} {:9.3} {:5.3} {:>10.3} {:>10.3}  {:.2?}", i, sampler.energy(&system), f_succ,
                 r_end_squared(&system.coordinates(), 0), gyration_squared(&&system.coordinates(), 0), start.elapsed());
    }

    coordinates_to_pdb(&system.coordinates(), 1i16, format!("{}final.pdb", &prefix).as_str(), false);
}
