use std::collections::HashMap;
use std::env;
use bioshell_pdb::{Deposit, PDBError};
use bioshell_seq::chemical::{ResidueTypeManager, ResidueTypeProperties};

fn main() -> Result<(), PDBError> {

    unsafe {
        if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    }
    env_logger::init();

    let args: Vec<String> = env::args().collect();
    let deposit = Deposit::from_file(&args[1]).unwrap();
    let date = deposit.dep_date.clone().unwrap();
    let strctr = deposit.structure()?;

    let mut d_aa_counts: HashMap<String, usize> = HashMap::default();
    strctr.residues().iter()
        .map(|res_id| strctr.residue_type(res_id).unwrap())
        .filter(|rtype| rtype.chem_compound_type.is_D_peptide_linking())
        .fold(&mut d_aa_counts, |acc, rtype| {
            *acc.entry(rtype.code3()).or_insert(0) += 1;
            acc
        });
    let all_aa_cnt = strctr.residues().iter()
        .map(|res_id| strctr.residue_type(res_id).unwrap())
        .filter(|rtype| rtype.chem_compound_type().is_peptide_linking()).count()  as f64;
    let d_aa_cnt  = d_aa_counts.values().sum::<usize>() as f64;
    if d_aa_counts.is_empty() {
        println!("{}: no D-amino acids found in the structure.", strctr.id_code);
    } else {
        print!("{} ({date}) {:5.1}%: ", strctr.id_code, (d_aa_cnt / all_aa_cnt) * 100.0);
        for (aa, count) in &d_aa_counts {
            print!("{}:{} ", aa, count);
        }
        println!();
    }

    Ok(())
}
