use std::env;
use bioshell_pdb::{Deposit, EntityType, PolymerEntityType};

const USAGE: &str = "Simple program to select the longest chain from each entity found in a given .cif file

    best_entity_chain input.cif
";

fn main() {
    unsafe {
        if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    }
    env_logger::init();

    let args: Vec<String> = env::args().collect();
    let mut fname: &str = "-h";
    match args.len() {
        1 => { eprintln!("No input CIF file given! Use distance_map -h for help") }
        2 => { fname = &args[1] }
        _ => { eprintln!("Too many input arguments! Use distance_map -h for help")}
    }
    if fname == "--help" || fname == "-h" {
        eprintln!("{:?}", USAGE);
        return;
    }

    let deposit = Deposit::from_file(&fname).unwrap();
    let strctr = deposit.structure().unwrap();
    for (_e_name, entity) in deposit.entities() {
        if let EntityType::Polymer(poly_type) = entity.entity_type() {
            if poly_type!=PolymerEntityType::PolypeptideL { continue; }
            let mut seq = vec![];
            for chain_id in entity.chain_ids() {
                seq.push((strctr.sequence(&chain_id), chain_id.clone()));
            }
            if seq.len() > 1 { seq.sort_by(|a, b| b.0.len().cmp(&a.0.len())); }

            println!(">{}\n{}", &seq[0].0.description(),&seq[0].0.to_string(0));
            println!(">{} secondary\n{}", &seq[0].0.description(),strctr.secondary(&seq[0].1).to_string());
        }
    }
}
