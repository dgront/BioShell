use std::env;
use bioshell_pdb::{load_pdb};

fn main(){
    let args: Vec<String> = env::args().collect();
    let mut strctr = load_pdb(&args[1]).unwrap();
    println!("Chains: {:?}",strctr.chain_ids());

    for chain_id in strctr.chain_ids() {
        println!("{}", strctr.sequence(&chain_id));
    }
}
