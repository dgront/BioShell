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

pub fn load_molecule(fname: &str) -> Result<Molecule, ChemErrors> {

    let path = Path::new(fname);

    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .ok_or_else(|| ChemErrors::UnknownFileFormat(fname.to_string()))?
        .to_ascii_lowercase();

    let reader = open_file(fname)?;

    match ext.as_str() {
        "sdf" | "mol"| "sdf.gz" => molecule_from_sdf(reader),
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
