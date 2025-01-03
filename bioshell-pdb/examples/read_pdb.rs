use std::env;
use bioshell_pdb::{Deposit};

fn main(){
    unsafe {
        if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    }
    env_logger::init();

    let args: Vec<String> = env::args().collect();
    let deposit = Deposit::from_file(&args[1]).unwrap();
    let strctr = deposit.structure().unwrap();
    println!("Chains: {:?}",strctr.chain_ids());

    for chain_id in strctr.chain_ids() {
        println!("{}", strctr.sequence(&chain_id));
    }

    println!("Are the atoms sorted: {}", if strctr.is_sorted() {"YES"} else {"NO"});
}
