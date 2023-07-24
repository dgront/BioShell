//! Module that defines chemical entities important for bioinformatics and molecular modeling
//!
//! # Residue typing
//!
//! Residue typing describes each chemical group or molecule that can be found
//! in a PDB or a FASTA file. While the latter files assume only 20 _standard_ amino acids and 5 nucleotides,
//! PDB files may contain a very wide variety of chemical components. Many of these components are
//! derived from _standard_ amino acids or bases by a chemical modification. Residue typing engine
//! facilitates conversion from a modified residue type into it's original (_parent_) counterpart.
//!
//! Since the set of _standard_ monomer is closed, all of them have been listed in [`StandardResidueType`](StandardResidueType) enum;
//! on the other hand [`ResidueType`](ResidueType) structs may be created dynamically.
//!
//! Finally, a [`ResidueTypeManager`](ResidueTypeManager) can provide a unique integer ID for each [`ResidueType`](ResidueType) instance
//!
mod residue_types;

pub use residue_types::{ResidueType, ResidueTypeProperties, ResidueTypeManager, MonomerType, StandardResidueType};