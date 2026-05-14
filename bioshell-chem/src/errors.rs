use thiserror::Error;
use bioshell_cif::CifError;

/// Error type for molecule graph operations.
#[derive(Debug, Error)]
pub enum ChemErrors {

    /// Invalid atomic number
    #[error("Invalid atomic number: {0}")]
    InvalidAtomicNumber(u8),

    /// Can't parse the given string as an element symbol.
    #[error("Unknown element symbol: {0}")]
    UnknownElement(String),

    /// Can't find an atom with the given index in the molecule.
    #[error("atom index out of bounds: {0}")]
    AtomIndexOutOfBounds(usize),

    /// Can't add an atom with the given index to the molecule, because it already exists.
    #[error("atom index {0} already exists in this molecule")]
    DuplicateAtomId(usize),

    /// Can't find a bond between the given atom indices.
    #[error("bond not found between atom indices {0} and {1}")]
    BondNotFound(usize, usize),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Cif error: {0}")]
    CifError(#[from] CifError),

    #[error("Error while parsing CIF value {0} for field {1}")]
    CifParsingError(String, String),

    #[error("Received unexpected atom name: {0}")]
    UnknownAtomName(String),

    #[error("missing required CIF field or column: {0}")]
    MissingCifField(String),
}