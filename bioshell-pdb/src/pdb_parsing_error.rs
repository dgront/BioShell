use thiserror::Error;
use crate::ResidueId;

/// Errors that may appear while using the bioshell-pdb crate
#[derive(Debug, Error)]
pub enum PDBError {
    #[error("Invalid PDB file format")]
    /// Invalid format of a PDB line
    InvalidFormat,
    #[error("I/O error: {0}")]
    /// I/O error occurred while reading a  PDB file
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
    #[error("TER record has incorrect format: {ter_line}")]
    /// A TER record line has incorrect format: `ter_line`
    IncorrectlyFormattedTER {ter_line: String},
    #[error("Can't find a file with input parameters: {fname} - check your BIOSHEL_DB_PATH system variable")]
    /// Missing parameters' file, that should be found in a BioShell's database
    MissingBioShellFile {fname: String},
}