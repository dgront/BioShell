use rand::Rng;
use std::env;
use std::time::Instant;
use std::io::stdout;
use std::io::{BufWriter,Write};
use std::path::Path;
use std::fs::{File};
use std::ops::Range;

use bioshell_numerical::Vec3;
use bioshell_ff::{Coordinates, Energy, TotalEnergy, ZeroEnergy};
use bioshell_ff::bonded::SimpleHarmonic;
use bioshell_ff::nonbonded::SimpleContact;
use bioshell_sim::generators::random_chain;
use bioshell_sim::sampling::movers::single_atom_move;
use bioshell_sim::sampling::protocols::{IsothermalMC, Sampler};

pub fn main() {
    let args: Vec<String> = env::args().collect();

    let n_res :usize = 100;
    let mut coords = Coordinates::new(n_res);
    random_chain(3.8,4.5,&mut coords);
    chain_to_pdb(&coords,"1.pdb");

    let harmonic = SimpleHarmonic::new(3.8,1.0);
    let contacts = SimpleContact::new(3.5,4.5,6.0,100.0,-1.0);

    let mut total = TotalEnergy::default();
    total.add_component(Box::new(harmonic),1.0).add_component(Box::new(contacts),1.0);
    println!("{}", total.energy(&coords));

    let mut sampler = IsothermalMC::new(1.0);
    sampler.energy = Box::new(total);    // --- The total has been moved to a box within the sampler
    let m: Box<dyn Fn(&mut Coordinates,f32) -> Range<usize>> = Box::new(single_atom_move);
    sampler.add_mover(m,1.0);

    sampler.run(&mut coords, 10);
}

pub fn chain_to_pdb(chain: &Coordinates, out_fname: &str) {

    let mut out_writer = match out_fname {
        "" => Box::new(stdout()) as Box<dyn Write>,
        _ => {
            let path = Path::new(out_fname);
            Box::new(File::create(&path).unwrap()) as Box<dyn Write>
        }
    };

    out_writer.write(b"MODEL    0\n").ok();
    for i in 0..chain.size() {
        out_writer.write(format!("ATOM   {:4}{}  ALA A{:4}    {:8.3}{:8.3}{:8.3}  1.00 99.88           C\n",
                                 i+1, " CA ", i+1, chain[i].x, chain[i].y, chain[i].z).as_bytes()).ok();
    }
    out_writer.write(b"ENDMDL\n").ok();
}