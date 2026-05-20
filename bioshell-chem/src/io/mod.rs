mod sdf;

use std::path::Path;
pub use sdf::*;

mod cif;
pub use cif::*;

mod mol2;
pub use mol2::*;

mod itp;
pub use itp::*;

use crate::{ChemErrors, Molecule};
use std::str::FromStr;
use bioshell_core::io::open_file;

/// Loads a molecule from a file based on its extension.
///
/// This function checks the extension of the provided file name
/// and calls the appropriate reader:
/// - for `.sdf`, `.sdf.gz`, `mol` or `mol.gz` : [`molecule_from_sdf()`]
/// - for `.cif` or `.cif.gz` : [`molecule_from_cif()`]
/// - for `.mol2` or `.mol2.gz` : [`molecule_from_mol2()`]
/// - for `.itp` or `.itp.gz` : [`molecule_from_itp()`]
///
/// If the extension is not recognized, the function results in an error.
///
/// # Example
/// ```
/// use bioshell_chem::{load_molecule, ChemErrors};
/// # fn main() -> Result<(), ChemErrors> {
/// let toluene_cif = load_molecule("./tests/test_files/MBN.cif")?;
/// let toluene_sdf = load_molecule("./tests/test_files/toluene.sdf")?;
/// assert!(toluene_cif.is_isomorphic_to(&toluene_sdf));
/// # Ok(())
/// # }
/// ```
pub fn load_molecule(fname: &str) -> Result<Molecule, ChemErrors> {

    let path = Path::new(fname);

    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .ok_or_else(|| ChemErrors::UnknownFileFormat(fname.to_string()))?
        .to_ascii_lowercase();

    let reader = open_file(fname)?;

    match ext.as_str() {
        "sdf" | "mol" | "mol.gz" | "sdf.gz" => molecule_from_sdf(reader),
        "cif" | "cif.gz" => molecule_from_cif(reader),
        "mol2" | "mol2.gz" => molecule_from_mol2(reader),
        "itp" | "itp.gz" => molecule_from_itp(reader),
        _ => Err(ChemErrors::UnknownFileFormat(fname.to_string())),
    }
}

pub(super) fn parse<T:FromStr>(token: &str, token_name: &str) -> Result<T, ChemErrors> {
    token.trim().parse::<T>()
        .map_err(|_| ChemErrors::NumericParsingError(token_name.to_string(), token.to_string()))
}

pub(super) fn parse_substr<T:FromStr>(token: &str, from:usize, to:usize, token_name: &str) -> Result<T, ChemErrors> {
    let substr = token.get(from..to)
        .ok_or_else(|| ChemErrors::LineTooShort(token_name.to_string(), to+1))?;
    parse(substr, token_name)
}
