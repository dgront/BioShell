use thiserror::Error;
use crate::ResidueId;

#[derive(Debug, Error)]
pub enum ParseError {
    #[error("Invalid PDB file format")]
    /// Invalid format of a PDB line
    InvalidFormat,
    #[error("I/O error: {0}")]
    /// I/O error occured while reading a  PDB file
    Io(#[from] std::io::Error),
    #[error("Residue not found: {res_id}")]
    /// Residue corresponding to a given `res_id` could not be located
    NoSuchResidue {res_id: ResidueId},
    #[error("Atom not found: {atom_name} in the residue {res_id}")]
    /// Atom named `atom_name` could not be located in a residue `res_id`
    NoSuchAtom {atom_name: String, res_id: ResidueId},
    #[error("Residue type not registered: {res_type}")]
    /// Unknown 3-letter residue code: `res_type`
    UnknownResidueType {res_type: String},
}