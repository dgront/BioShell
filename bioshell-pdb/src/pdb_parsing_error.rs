use thiserror::Error;
use crate::ResidueId;

#[derive(Debug, Error)]
pub enum ParseError {
    #[error("Invalid PDB file format")]
    InvalidFormat,
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("Residue not found: {res_id}")]
    NoSuchResidue {res_id: ResidueId},
    #[error("Residue type not registered: {res_type}")]
    UnknownResidueType {res_type: String},
}