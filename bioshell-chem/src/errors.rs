use thiserror::Error;

/// Error type for molecule graph operations.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum ChemErrors {
    #[error("atom index out of bounds: {0}")]
    AtomIndexOutOfBounds(usize),

    #[error("atom index {0} already exists in this molecule")]
    DuplicateAtomId(usize),

    #[error("bond not found between atom indices {0} and {1}")]
    BondNotFound(usize, usize),
}