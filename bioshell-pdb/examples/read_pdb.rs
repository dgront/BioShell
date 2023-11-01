use std::env;
use bioshell_pdb::{load_pdb_file};

fn main(){
    let args: Vec<String> = env::args().collect();
    let strctr = load_pdb_file(&args[1]).unwrap();
    println!("Chains: {:?}",strctr.chain_ids());

    for chain_id in strctr.chain_ids() {
        println!("{}", strctr.sequence(&chain_id));
    }
}
