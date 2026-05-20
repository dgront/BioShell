use std::env;
use bioshell_chem::{ChemErrors, load_molecule, Molecule, molecule_from_cif};
use bioshell_chem::ChemErrors::InvalidAtomIndex;
use bioshell_chem::icoords::{InternalCoordinate, KinematicAtomTree};

use std::io::Write;
use clap::builder::TypedValueParser;

fn write_z_matrix<W: Write>(mut out: W, molecule: &Molecule, kac: &KinematicAtomTree, icoords: &[InternalCoordinate], ) -> Result<(), ChemErrors> {

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
            ka.a,
            ka.b,
            ka.c,
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
    let molcle = load_molecule(&args[1])?;
    let kac = KinematicAtomTree::from_molecule(&molcle, 0, 1, 2)?;

    // ---------- Extract the source Cartesian coordinates
    let cartesian = molcle.atoms().map(|a|a.pos().clone()).collect::<Vec<_>>();

    // ---------- Compute the internal coordinates and print the z-matrix
    let icoords = kac.get_icoords(&cartesian);

    write_z_matrix(std::io::stdout(), &molcle, &kac, &icoords?)?;

    Ok(())
}
