//! Efficient and clean library for processing biomacromolecular structures.
//!
//! # Loading PDB files
//! Structures can be loaded from a PDB file with [`load_pdb_file()`](bioshell-pdb::load_pdb_file())
//! or directly from a PDB-formatted buffer with [`load_pdb_reader()`](bioshell-pdb::load_pdb_reader())
//!
mod structure;
mod pdb_compound;
mod pdb_header;
mod pdb_sequence_of_residue;
mod pdb_source;
mod pdb_title;
pub mod pdb_parsing_error;
mod pdb_helix;
mod pdb_sheet;
pub mod pdb_helix_line_parser;
pub mod pdb_sheet_line_parser;
mod pdb_atom;

pub mod pdb_atom_filters;
mod residue_id;
pub mod calc;

pub use structure::{Structure, load_pdb_file, load_pdb_reader};

pub use pdb_atom::PdbAtom;
pub use residue_id::{ResidueId, residue_id_from_ter_record};
pub use pdb_compound::PdbCompound;
pub use pdb_header::PdbHeader;
pub use pdb_source::PdbSource;
pub use pdb_title::PdbTitle;
pub use pdb_helix::PdbHelix;
pub use pdb_sheet::PdbSheet;


