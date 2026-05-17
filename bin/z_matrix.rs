use std::env;
use log::{info};
use bioshell_chem::{Molecule, molecule_from_cif};
use bioshell_chem::ChemErrors::InvalidAtomIndex;
use bioshell_chem::icoords::KinematicAtomChain;
use bioshell_core::io::open_file;
use bioshell_pdb::read_cif_monomers;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    let args: Vec<String> = env::args().collect();

    // ---------- Prepare the molecule and the kinematic atom chain
    let molcle = molecule_from_cif(open_file(&args[1])?)?;
    let kac = KinematicAtomChain::from_molecule(&molcle, 0, 1, 2)?;

    // ---------- Load the structure and get the source Cartesian coordinates
    let strctr = read_cif_monomers(open_file(&args[1])?, None)?;
    info!("{}", format!("Number of monomers found: {}", strctr.count_residues()));
    let cartesian = strctr.atoms().iter().map(|a| a.pos).collect::<Vec<_>>();

    assert_eq!(kac.len(), cartesian.len(), "Number of atoms in the molecule and the structure should be the same");

    // ---------- Compute the internal coordinates and print the z-matrix
    let icoords = kac.get_icoords(&cartesian);

    for i in 0..icoords.len() {
        let ipos = &icoords[i];
        let a = molcle.get_atom(i).ok_or_else(||InvalidAtomIndex(i))?;
        println!("{:3} {:2} {:3} {:3} {:3} : {:5.3} {:6.2} {:7.2}",
                 i, a.element(), kac[i].i, kac[i].j, kac[i].k, ipos.d, ipos.alpha.to_degrees(), ipos.phi.to_degrees());
    }
    Ok(())
}
