use std::env;
use clap::error::ContextKind::Usage;
use bioshell_pdb::{load_pdb_file, PdbAtom};
use bioshell_pdb::calc::distance;
use bioshell_pdb::pdb_atom_filters::{AlwaysPass, IsBackbone, IsCA, PdbAtomPredicate};

fn atom_descr(a: &PdbAtom) -> String {
    format!("{} {} {:4} {}",a.chain_id, a.res_name, a.res_seq, a.name)
}

const USAGE: &str = "Simple program to compute a matrix of inter-atomic distances for a given PDB file

    distance_map [--ca|--bb] input.pdb
";
fn main() {
    let args: Vec<String> = env::args().collect();
    let mut selector: Box<dyn PdbAtomPredicate> = Box::new(AlwaysPass);
    let mut fname = "";
    match args.len() {
        1 => { panic!("No input PDB file given") }
        2 => { fname = &args[1] }
        3 => {
            match args[1].as_str() {
                "--ca" => { selector = Box::new(IsCA); }
                "--bb" => { selector = Box::new(IsBackbone); }
                _ => { panic!("Unknown PDB atom predicate; run `distance_map -h` for help")}
            }
            fname = &args[2];
        }
        _ => { panic!("Too many arguments; run `distance_map -h` for help") }
    };
    if fname == "--help" || fname == "-h" {
        println!("{:?}", Usage);
        return;
    }
    let mut strctr = load_pdb_file(&fname).unwrap();

    for ai in strctr.atoms().iter().filter(|&a|selector.check(a)) {
        for aj in strctr.atoms().iter().filter(|&a|selector.check(a)) {
            if ai.res_seq == aj.res_seq && ai.i_code==aj.i_code { break }
            println!("{} {} : {}", atom_descr(ai), atom_descr(aj), distance(ai, aj));
        }
    }
}
