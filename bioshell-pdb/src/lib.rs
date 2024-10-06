//! Efficient and clean library for processing biomacromolecular structures.
//!
//! # Loading PDB files
//! Structures can be loaded from a PDB file with [`load_pdb_file()`](crate::load_pdb_file())
//! or directly from a PDB-formatted buffer with [`load_pdb_reader()`](crate::load_pdb_reader()).
//! When successful, each of these method returns a [`Structure`](crate::Structure) object that
//! holds all the atoms parsed from the input.
//!
//! # Selecting chains, residues and atoms
//! The `bioshell-pdb` crate provides access to the vector of atoms for a given [`Structure`](crate::Structure)
//! which may be processed as any Rust [`Iterator`](std::iter::Iterator) method, e.g. filtered with
//! [`filter()`](std::iter::Iterator::filter()) or mapped with [`map()`](std::iter::Iterator::map()).
//! The [`pdb_atom_filters`](crate::pdb_atom_filters) module provides several predicates for such applications.
//! See also the documentation for [`Structure`](crate::Structure) struct for more examples.
//!
//! # Structural calculation
//! The [`calc`](crate::calc) module provides functions to calculate structural properties,
//! such as distances, planar or dihedral angles.

#![allow(clippy::needless_return)]
mod structure;
mod pdb_header;
mod secondary_structure;
mod pdb_title;
mod pdb_parsing_error;
mod pdb_helix;
mod pdb_sheet;
mod pdb_atom;
mod assertions;
pub mod pdb_atom_filters;
pub mod residue_filters;
pub mod calc;

mod residue_id;
mod load_pdb;
mod exp_data;
mod remarks;
mod unit_cell;
mod load_cif;
mod entity;
pub(crate) mod crate_utils;

pub use structure::Structure;
pub use secondary_structure::{SecondaryStructureTypes,SecondaryStructure};
pub use load_pdb::{load_pdb_file, load_pdb_reader, is_pdb_file, find_pdb_file_name};
pub use assertions::*;
pub use pdb_parsing_error::PDBError;
pub use pdb_atom::{PdbAtom, same_residue_atoms};
pub use residue_id::{ResidueId, residue_id_from_ter_record};
pub use pdb_helix::PdbHelix;
pub use pdb_sheet::PdbSheet;
pub use exp_data::{ExperimentalMethod};
pub use unit_cell::{UnitCell};
pub use load_cif::{load_cif_file, load_cif_reader, is_cif_file, find_cif_file_name};
pub use entity::{EntitySource, EntityType, Entity};

mod macros {
    #[macro_export]
    macro_rules! value_or_missing_key_pdb_error {
        ($cif_data_block:expr, $data_key:expr, $casted_type:ty) => {
            $cif_data_block.get_item::<$casted_type>($data_key)
                .ok_or( MissingCifDataKey{ data_name: $data_key.to_string()}).map_err(|e|CifParsingError(e))?
        };
    }
}


