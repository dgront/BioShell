use std::env;
use log::{info};
use bioshell_chem::{ChemErrors, Molecule, molecule_from_cif};
use bioshell_chem::ChemErrors::InvalidAtomIndex;
use bioshell_chem::icoords::{InternalCoordinate, KinematicAtomChain};
use bioshell_core::io::open_file;
use bioshell_pdb::read_cif_monomers;

use std::io::Write;

fn write_z_matrix<W: Write>(mut out: W, molecule: &Molecule, kac: &KinematicAtomChain, icoords: &[InternalCoordinate], ) -> Result<(), ChemErrors> {

    if icoords.len() != kac.len() {
        return Err(ChemErrors::IncorrectNumberOfAtoms(icoords.len(), kac.len()));
    }

    for i in 0..icoords.len() {
        let ipos = &icoords[i];
        let ka = &kac[i];
        let atom = molecule.get_atom(ka.atom)
            .ok_or(InvalidAtomIndex(ka.atom))?;

        writeln!(
            out,
            "{:3} {:2} {:3} {:3} {:3} : {:5.3} {:6.2} {:7.2}",
            ka.atom,
            atom.element(),
            ka.i,
            ka.j,
            ka.k,
            ipos.d,
            ipos.alpha.to_degrees(),
            ipos.phi.to_degrees()
        )?;
    }

    Ok(())
}

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

    write_z_matrix(std::io::stdout(), &molcle, &kac, &icoords?)?;

    Ok(())
}
