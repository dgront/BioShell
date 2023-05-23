use std::any::Any;
use clap::builder::TypedValueParser;
use std::env;
use std::f32::consts::PI;
use std::time::Instant;

use clap::Parser;
use log::info;

use bioshell_cartesians::movers::{CrankshaftMove, TerminalMove};//import single-atom-mover
use bioshell_cartesians::observers::{GyrationSquared, PdbTrajectory, REndSquared};
use bioshell_cartesians::{
    box_width, coordinates_to_pdb, pdb_to_coordinates, CartesianSystem, Coordinates, NbList,
    PolymerRules, RandomChain,
};
use bioshell_core::utils::out_writer;

use bioshell_sim::{Energy, Observer, ObserversSet, System};

use bioshell_ff::nonbonded::{PairwiseNonbondedEvaluator, SimpleContact};
use bioshell_montecarlo::{AdaptiveMCProtocol, IsothermalMC, Sampler, StepwiseBuilder};

use bioshell_ff::bonded::SimpleHarmonic;
use bioshell_ff::TotalEnergy;

#[derive(Parser, Debug)]
#[clap(name = "polymer")]
#[clap(version = "0.2")]
#[clap(about = "Simple polymer model", long_about = None)]
/// simulation of a simple polymer model
/// say polymer -h to see options
struct Args {
    /// staring conformation in the PDB format
    #[clap(short, long, default_value = "", short = 'f')]
    infile: String,
    /// temperature of the simulation
    #[clap(short, long, default_value_t = 1.0, short = 't')]
    temperature: f64,
    /// Number of beads of a polymer chain
    #[clap(short, long, default_value_t = 100, short = 'n')]
    n_beads: usize,
    /// by-volume density of the simulated system
    #[clap(short, long, default_value_t = 0.4, short = 'd')]
    density: f64,
    /// Inner cycles
    #[clap(short, long, default_value_t = 100)]
    inner: usize,
    /// outer cycles
    #[clap(short, long, default_value_t = 100)]
    outer: usize,
    /// prefix for output file names
    #[clap(long, default_value = "")]
    prefix: String,
    /// thickness of the buffer zone used by non-bonded list to hash non-bonded interactions
    #[clap(long, default_value_t = 4.0)]
    buffer: f64,
}

pub fn main() {
    // --- This represents the repulsion distance between the beads in the simulation.
    const E_REP: f64 = 4.25;
    // It is the distance range where a contact energy is awarded (to model the interactions between particles).
    const E_FROM: f64 = 4.5;
    const E_TO: f64 = 6.0;
    /*
    The program  assigned an energy value of E_VAL = -1.0, when a distance between beads falls in the
    range between E_FROM and E_TO. When the distance is shorted than E_REP, the repulsion energy
    REP_VAL = 1000.0 is applied
    */
    const E_VAL: f64 = -1.0;
    const REP_VAL: f64 = 1000.0;
    const MAX_MOVE_RANGE: f64 = (15.0 * 180.0 / PI) as f64;    // --- max move range is 15 degrees covnerted to  radians

    if env::var("RUST_LOG").is_err()
    {
        env::set_var("RUST_LOG", "info")
    }
    env_logger::init();

    let args = Args::parse();
    // ---------- Temperature of the isothermal simulation (in the energy units)
    let temperature: f64 = args.temperature; // --- Temperature of the isothermal simulation (in the energy units)
    let density: f64 = args.density; // --- density of the system
    let box_length: f64;
    let prefix = args.prefix; // --- prefix for the output files
    let tra_f_name = format!("{}tra.pdb", &prefix);

    // ---------- Non-bonded List buffer radius in Angstroms
    let buffer_thickness = args.buffer.max(5.0);

    // ---------- Create system's coordinates
    let n_beads: usize;//declare a variable to keep track of number of beads
    let mut coords: Coordinates;//declare a variable for the coordinates of polymer beads

    if args.infile == "" {//if 'infile' is empty,
        n_beads = args.n_beads;//set the number of beads
        coords = Coordinates::new(n_beads);//create `n_beads` number of coordinates
        box_length = box_width(E_REP, n_beads, density);//calculate 3D box length
        coords.set_box_len(box_length);//set box length to the calculated length
    }
    else
    {
        let res = pdb_to_coordinates(&args.infile).ok();//converts Result type to Option-type
        match res
        {   //Some() is used to extract the value from an Option when it is not None.
            Some(x) => {
                coords = x;
            }
            _ => panic!("Can't read from file >{}<", args.infile),
        };
        n_beads = coords.size();//get number of atoms
        box_length = box_width(E_REP, n_beads, density);//calculate box length
        coords.set_box_len(box_length);//set coords to box length
    }

    info!(
        "Simulation settings:
        temperature: {:.3}
        no. beads:   {}
        box length:  {:.3}
        NBL buffer:  {:.3}",
        temperature,
        n_beads,
        box_length,
        buffer_thickness
    );

    // ---------- Create a non-bonded list of neighbors
    let nbl: NbList = NbList::new(E_TO, buffer_thickness, Box::new(PolymerRules {}));
    
	// ---------- Create the system
    let mut system: CartesianSystem = CartesianSystem::new(coords, nbl);
    
	// ---------- Contact energy
    let contact_kernel = SimpleContact::new(E_REP, E_FROM, E_TO, REP_VAL, E_VAL);
    let contacts: PairwiseNonbondedEvaluator<SimpleContact> =
        PairwiseNonbondedEvaluator::new(E_TO as f64, contact_kernel);
    
	// ---------- Harmonic energy (i.e. springs between beads)
    let harmonic = SimpleHarmonic::new(3.8, 1.0);
    
	// ---------- Total energy contains the contacts energy and bond springs
    let mut total = TotalEnergy::new();
    total.add_component(Box::new(harmonic), 1.0);
    total.add_component(Box::new(contacts), 1.0);
    if system.size() == 0 {
        RandomChain::default().build(&mut system, &total);
    }

    // ---------- Show the starting point
    coordinates_to_pdb(&system.coordinates(), 0, tra_f_name.as_str(), false);
    println!("# starting energy: {}", total.energy(&system));

    // ---------- Create a sampler and add a mover into it
    let mut simple_sampler: IsothermalMC<CartesianSystem, TotalEnergy<CartesianSystem>> =
        IsothermalMC::new(temperature);
    simple_sampler.add_mover(Box::new(CrankshaftMove::new(MAX_MOVE_RANGE)));
    simple_sampler.add_mover(Box::new(TerminalMove::new(MAX_MOVE_RANGE)));

    // ---------- Decorate the sampler into an adaptive MC protocol
    let mut sampler = AdaptiveMCProtocol::new(Box::new(simple_sampler));
    sampler.target_rate = 0.4;

    let mut observations: ObserversSet<CartesianSystem> = ObserversSet::new();
    observations.add_observer(Box::new(PdbTrajectory::new(tra_f_name, true)), 1);
    observations.add_observer(
        Box::new(GyrationSquared::new("rg.dat".to_string(), true)),
        1,
    );
    observations.add_observer(Box::new(REndSquared::new("r2.dat".to_string(), true)), 1);
    observations.add_observer(Box::new(BondViolationObserver::new("bond_violations.dat".to_string(), true)), 1);
    sampler.run_simulation(
        args.inner,
        args.outer,
        &mut system,
        &total,
        &mut observations,
    );

    coordinates_to_pdb(
        &system.coordinates(),
        1i16,
        format!("{}final.pdb", &prefix).as_str(),
        false,
    );
}

pub struct BondViolationObserver {
    pub out_fname: String,
    pub if_append: bool,
    i_model: usize,
}

impl BondViolationObserver {
    pub fn new(fname: String, if_append: bool) -> BondViolationObserver {
        BondViolationObserver {
            out_fname: fname,
            if_append,
            i_model: 0,
        }
    }
}

impl Observer for BondViolationObserver {
    type S = CartesianSystem;

    fn observe(&mut self, object: &Self::S) {
        let coords = object.coordinates();
        let mut out_writer = out_writer(&self.out_fname, self.if_append);
        out_writer
            .write(format!("{:.6} ", self.i_model).as_bytes())
            .ok();
        for ic in 0..coords.count_chains() {
            // HERE !!!
            let result: f64 = 0.12345;
            out_writer.write(format!("{:>10.3} ", result).as_bytes()).ok();
        }
        out_writer.write("\n".as_bytes()).ok();
        self.i_model += 1;
    }

    fn flush(&mut self) {}

    fn name(&self) -> &str {
        "BondViolationObserver"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}