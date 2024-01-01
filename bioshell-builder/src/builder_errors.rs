use thiserror::Error;

/// Errors that may appear while using the bioshell-builder crate
#[derive(Debug, Error)]
pub enum BuilderError {

    #[error("Can't find a file with input parameters: {fname} - check your BIOSHEL_DB_PATH system variable")]
    /// Missing parameters' file, that should be found in a BioShell's database
    MissingDataFile {fname: String},

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