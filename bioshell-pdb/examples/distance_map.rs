use std::env;
use bioshell_pdb::{load_pdb_file, PdbAtom};
use bioshell_pdb::calc::distance;
use bioshell_pdb::pdb_atom_filters::{IsBackbone, IsCA, PdbAtomPredicate};

fn atom_descr(a: &PdbAtom) -> String {
    format!("{} {} {:4} {}",a.chain_id, a.res_name, a.res_seq, a.name)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut strctr = load_pdb_file(&args[1]).unwrap();
    strctr.drop_ligands();

    // let selector = IsCA;
    let selector = IsBackbone;
    for ai in strctr.atoms().iter().filter(|&a|selector.check(a)) {
        for aj in strctr.atoms().iter().filter(|&a|selector.check(a)) {
            println!("{} {} : {}", atom_descr(ai), atom_descr(aj), distance(ai, aj));
        }
    }
}
