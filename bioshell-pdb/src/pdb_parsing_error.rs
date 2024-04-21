use thiserror::Error;
use crate::ResidueId;

/// Errors that may appear while using the bioshell-pdb crate
#[derive(Debug, Error)]
pub enum PDBError {
    #[error("Invalid PDB file format: {broken_pdb_line}")]
    /// Invalid format of a PDB line
    InvalidFormat {broken_pdb_line: String},

    #[error("CIF input contains more than one data block")]
    /// CIF input contains more than one data block
    ExtraDataBlock,

    #[error("I/O error")]
    /// I/O error occurred while reading a PDB or a CIF file
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

    #[error("Error while creating an InternalAtomDefinition struct")]
    /// Can't create `InternalAtomDefinition`
    InternalAtomDefinitionError {error: String},

    #[error("Can't find atom named: {atom_name}. Has it been defined in the residue {residue_index}?")]
    /// Can't find an atom defined for a given residue
    DefinedAtomNotFound {residue_index: usize, atom_name: String},
    
    #[error("Can't find a residue {residue_index}?")]
    /// Can't find a residue for a given index
    ResidueNotDefined {residue_index: usize},
}