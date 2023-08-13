use std::env;
use bioshell_pdb::{Structure, load_pdb};

fn main(){
    let args: Vec<String> = env::args().collect();
    let mut strctr = load_pdb(&args[1]).unwrap();
    println!("Chains: {:?}",strctr.chain_ids());

    for chain_id in strctr.chain_ids() {
        let res_ids = strctr.chain_residue_ids(&chain_id);
        println!("\nChain: {} has {} residues:", chain_id, res_ids.len());
        for rid in  res_ids {
            print!("{} ",&rid);
        }
    }
}
