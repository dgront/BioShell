use thiserror::Error;


/// Errors that may appear while using the bioshell-clustering crate
#[derive(Debug, Error)]
pub enum ClusteringError {
    #[error("Data parsing error: {reason}; incorrect data: {data}")]
    /// A given file is neither a PDB nor an mmCIF file
    InvalidDataFormat { reason: String, data: String },

    #[error("General I/O error occurred while reading an input file")]
    /// I/O error occurred while reading an input file
    Io(#[from] std::io::Error),

    #[error("Error returned by the CSV parser")]
    /// Error returned by the CSV parser
    CsvError(#[from] csv::Error),
}