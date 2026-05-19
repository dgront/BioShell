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
    DuplicatedAtomIndex(usize),

    /// Can't find any atom with the given index.
    #[error("Can't find any atom by the: index {0}")]
    InvalidAtomIndex(usize),

    /// THe number of atoms in the molecule is different from the expected one.
    #[error("Expected {0} atoms but found: {1}")]
    IncorrectNumberOfAtoms(usize, usize),

    /// Invalid reference atoms; the three reference atoms must be different, but the given indexes are not.
    #[error("Invalid reference atoms: {0}, {1}, {2}. All the three indexes must be different.")]
    InvalidReferenceAtoms(usize, usize, usize),

    /// Can't find a bond between the given atom indices.
    #[error("bond not found between atom indices {0} and {1}")]
    BondNotFound(usize, usize),

    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Can't recognize file format based on the extension: {0}")]
    UnknownFileFormat(String),

    #[error("Cif error: {0}")]
    CifError(#[from] CifError),

    #[error("Error while parsing CIF value {0} for field {1}")]
    NumericParsingError(String, String),

    #[error("Received unexpected atom name: {0}")]
    UnknownAtomName(String),

    #[error("missing required CIF field or column: {0}")]
    MissingCifField(String),

    /// Input line is too short
    #[error("Line too short, expected {1} characters: {0}")]
    LineTooShort(String, usize),

    /// Incorrect SDF format
    #[error("Incorrect SDF line: {0}")]
    IncorrectSdfFormat(String),

    /// Incorrect Gromacs topology ITP format
    #[error("Incorrect Gromacs topology ITP format: {0}")]
    IncorrectItpFormat(String),

    #[error("invalid MOL2 atom line: {0}")]
    InvalidMol2AtomLine(String),

    #[error("invalid MOL2 bond line: {0}")]
    InvalidMol2BondLine(String),

    #[error("invalid MOL2 atom id: {0}")]
    InvalidMol2AtomId(String),

    #[error("incorrect SMILES string: {smiles}; {reason}")]
    IncorrectSmilesString {
        smiles: String,
        reason: String,
    },
}