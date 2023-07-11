use thiserror::Error;

#[derive(Debug, Error)]
pub enum PdbParseError {
    #[error("Invalid PDB file format")]
    InvalidFormat,
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}