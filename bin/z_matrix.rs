use std::env;
use log::{info};
use bioshell_chem::{Molecule, molecule_from_cif};
use bioshell_chem::icoords::KinematicAtomChain;
use bioshell_io::open_file;
use bioshell_pdb::read_cif_monomers;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    let args: Vec<String> = env::args().collect();

    let strctr = read_cif_monomers(open_file(&args[1])?, None)?;
    info!("{}", format!("Number of monomers found: {}", strctr.count_residues()));
    let pos = strctr.atoms().iter().map(|a| a.pos).collect::<Vec<_>>();

    let molcle = molecule_from_cif(open_file(&args[1])?)?;
    let kt = KinematicAtomChain::from_molecule(&molcle, 0, 1, 2)?;
    assert_eq!(kt.len(), pos.len(), "Number of atoms in the molecule and the structure should be the same");

    Ok(())
}
