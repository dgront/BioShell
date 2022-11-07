//! Module that defines chemical entities important for bioinformatics and molecular modeling
//!
//! # Residue typing
//!
//! Residue typing provides a unique integer ID for each chemical group or molecule that can be found
//! in a PDB or FASTA file. While the latter assume only 20 _standard_ amino acids and 5 bases,
//! PDB files may contain a very wide variety of chemical components. In order to map them efficiently,
//! the standard monomers are listed in ``StandardResidueType`` enum while all other can be defined
//! as a ``ResidueType`` struct.
//!
mod residue_types;

pub use residue_types::{ResidueType, ResidueTypeProperties, ResidueTypeManager, MonomerType, StandardResidueType};