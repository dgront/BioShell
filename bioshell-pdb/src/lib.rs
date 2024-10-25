//! Efficient and clean library for processing biomacromolecular structures.
//!
//! # Loading CIF and PDB deposits
//! Biomacromolecular deposits (either in mmCIF or PDB file format) can be directly loaded into a [`Deposit`](Deposit) struct:
//!
//!```no_run
//! use bioshell_pdb::Deposit;
//! let deposit = Deposit::from_file("2gb1.cif");
//!```
//! the file format is determined automatically by its header. A [`Deposit`](Deposit) struct
//! can be also directly loaded from a data buffer; in such a case user must specify the file format
//! by calling the appropriate method:
//!```
//! # use std::io::BufReader;
//! # use bioshell_pdb::{Deposit, PDBError};
//! # fn main() -> Result<(), PDBError> {
//! # let pdb_data = include_str!("../tests/test_files/2gb1.pdb");
//! let reader = BufReader::new(pdb_data.as_bytes());
//! let deposit_2gb1 = Deposit::from_pdb_reader(reader)?;
//! # let cif_data = include_str!("../tests/test_files/2fdo.cif");
//! let reader = BufReader::new(cif_data.as_bytes());
//! let deposit_2fdo = Deposit::from_cif_reader(reader)?;
//! # Ok(())
//! # }
//! ```
//! Once successfully loaded, it provides access to the information stored in a file, such as
//! [`resolution`](Deposit::title), [`UnitCell`](UnitCell) or  [`resolution`](Deposit::resolution).
//!```
//! # use std::io::BufReader;
//! # use bioshell_pdb::{Deposit, PDBError};
//! # fn main() -> Result<(), PDBError> {
//! # let cif_data = include_str!("../tests/test_files/2fdo.cif");
//! # let reader = BufReader::new(cif_data.as_bytes());
//! let deposit_2fdo = Deposit::from_cif_reader(reader)?;
//! println!("Title: {}", deposit_2fdo.title.unwrap());
//! # Ok(())
//! # }
//! ```
//!
//!
//! # Structures
//!
//! A [`Deposit`](Deposit) struct can obviously provide also the atomic coordinates, which are
//! stored in a [`Structure`](Structure) struct:
//! ```
//! # use std::io::BufReader;
//! # use bioshell_pdb::{Deposit, PDBError};
//! # fn main() -> Result<(), PDBError> {
//! # let cif_data = include_str!("../tests/test_files/2gb1.cif");
//! # let reader = BufReader::new(cif_data.as_bytes());
//! # let deposit = Deposit::from_cif_reader(reader)?;
//! let strctr = deposit.structure();
//! # Ok(())
//! # }
//! ```
//!
//! # Entities
//!
//! When loaded from an mmCIF file, a [`Deposit`](Deposit) struct provides also detailed information
//! about [`Entities`](Entity) found in that deposit:
//! ```
//! # use std::io::BufReader;
//! # use bioshell_pdb::{Deposit, PDBError};
//! # fn main() -> Result<(), PDBError> {
//! # let cif_data = include_str!("../tests/test_files/2fdo.cif");
//! # let reader = BufReader::new(cif_data.as_bytes());
//! let deposit = Deposit::from_cif_reader(reader)?;
//! println!("Number of entities: {}", deposit.count_entities());
//! # assert_eq!(deposit.count_entities(), 2);
//! let first_entity = deposit.entity("1");
//! println!("First entity in chains: {:?}", first_entity.chain_ids());
//! # assert_eq!(first_entity.chain_ids(), &vec!["B", "A"]);
//! # Ok(())
//! # }
//! ```
//!
//!
//! # Selecting chains, residues and atoms
//! The `bioshell-pdb` crate provides access to the vector of atoms for a given [`Structure`](Structure)
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

pub mod monomers;
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
mod deposit;
mod ligands;

pub use structure::Structure;
pub use deposit::*;
pub use ligands::*;
pub use secondary_structure::{SecondaryStructureTypes,SecondaryStructure};
pub use load_pdb::{is_pdb_file, find_pdb_file_name};
pub use pdb_parsing_error::PDBError;
pub use pdb_atom::{PdbAtom, same_residue_atoms};
pub use residue_id::{ResidueId};
pub use pdb_helix::PdbHelix;
pub use pdb_sheet::PdbSheet;
pub use exp_data::{ExperimentalMethod};
pub use unit_cell::{UnitCell};
pub use load_cif::{is_cif_file, find_cif_file_name};
pub use entity::{EntitySource, EntityType, Entity, PolymerEntityType};

mod macros {
    #[macro_export]
    macro_rules! value_or_missing_key_pdb_error {
        ($cif_data_block:expr, $data_key:expr, $casted_type:ty) => {
            $cif_data_block.get_item::<$casted_type>($data_key)
                .ok_or( MissingCifDataKey{ data_name: $data_key.to_string()}).map_err(|e|CifParsingError(e))?
        };
    }
}


