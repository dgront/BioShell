use thiserror::Error;

/// Errors that may appear while running a SURPASS simulation
#[derive(Debug, Error)]
pub enum SurpassError {
    
    #[error("Invalid secondary structure type:  {ss_type}")]
    /// Invalid character used to denote a secondary structure type
    InvalidSecondaryStructureCode {ss_type: char},

    #[error("Invalid secondary structure character {ss_type} found in the input secondary structure string:  {ss_string}")]
    /// Invalid character found in an input secondary structure string
    InvalidSecondaryStructureString {ss_type: char, ss_string: String},
}