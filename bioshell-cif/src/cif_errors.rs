use thiserror::Error;

/// Errors that may appear while using the bioshell-cif crate
#[derive(Debug, Error)]
pub enum CifError {

    #[error("General I/O error occurred while reading a CIF file")]
    /// I/O error occurred while reading a PDB or a CIF file
    Io(#[from] std::io::Error),

    #[error("Data key of a loop block not found in CIF input: {item_key}")]
    /// A loop block miss a mandatory column identified by the key: `item_key`
    MissingCifLoopKey {item_key: String},

    #[error("CIF input contains more than one data block")]
    /// CIF input contains more than one data block
    ExtraDataBlock,

    #[error("The following mandatory data item can't be found in a CIF input: {item_key}")]
    /// CIF input contains more than one data block
    MissingCifDataKey {item_key: String},
}