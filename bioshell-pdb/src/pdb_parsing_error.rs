use thiserror::Error;
use bioshell_cif::CifError;
use crate::ResidueId;

/// Errors that may appear while using the bioshell-pdb crate
#[derive(Debug, Error)]
pub enum PDBError {
    #[error("Invalid PDB file format: {broken_pdb_line}")]
    /// Invalid format of a PDB line
    InvalidFormat {broken_pdb_line: String},

    #[error("General I/O error occurred while reading a PDB or mmCIF file")]
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
    
    #[error("Can't find a residue {residue_index}")]
    /// Can't find a residue for a given index
    ResidueNotDefined {residue_index: usize},

    #[error("The number of atoms in model {model_index} is different than the number of atoms in the first model")]
    /// A requested model has different number of atoms that the first model
    WrongAtomsNumberInModel {model_index: usize},

    #[error("Unknown chemical compound type {compound_type} used to define monomer {compound_id}")]
    /// Can't parse a chemical component type
    IncorrectCompoundTypeName {compound_id: String, compound_type: String},

    #[error("Cif parser returned an error")]
    /// A Cif parser returned an error which is stored inside this enum variant
    CifParsingError(#[from] CifError),

    #[error("The following string {data_value} can't be parsed into an enum {enum_name} variant")]
    /// A given string can't be parsed into an enum variant
    CantParseEnumVariant {data_value: String, enum_name: String},
}