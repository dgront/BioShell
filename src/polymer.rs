use rand::Rng;
use std::env;
use std::time::Instant;
use std::io::stdout;
use std::io::{BufWriter,Write};
use std::path::Path;
use std::fs::{File};

use bioshell_numerical::Vec3;
use bioshell_ff::{Coordinates, Energy};
use bioshell_ff::bonded::SimpleHarmonic;
use bioshell_ff::nonbonded::SimpleContact;
use bioshell_sim::generators::random_chain;


pub fn main() {
    let args: Vec<String> = env::args().collect();

    let n_res :usize = 100;
    let mut coords = Coordinates::new(n_res);
    random_chain(3.8,4.5,&mut coords);
    chain_to_pdb(&coords,"1.pdb");

    let harmonic = SimpleHarmonic::new(3.8,1.0);
    let contacts = SimpleContact::new(3.5,4.5,6.0,100.0,-1.0);
    let en_h = harmonic.energy(&coords);
    let en_c = contacts.energy(&coords);
    println!("{} {}", en_h, en_c);
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