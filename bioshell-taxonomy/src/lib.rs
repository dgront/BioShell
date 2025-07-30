//! A fast and lightweight library for working with the NCBI Taxonomy database in Rust.
//!
//! This crate provides efficient tools for loading and querying taxonomic data
//! from the [NCBI taxonomy dump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/).
//! It supports scientific names, common names, and synonyms, allowing fast resolution
//! of taxonomic IDs (taxids), lineage queries, and rank-based lookups.
//!
//! ## Features
//! - Assign taxonomic rank using the compact `Rank` enum.
//! - Map species names (scientific or common) to taxonomic IDs.
//! - Traverse lineage from any node to the root.
//! - List all names (including synonyms) for a given taxid.
//! - Search nodes by rank (e.g., all kingdoms).
//!
//! ## Command line application
//! Basic functionality of the library has been provided as a command line application ``taxonomy``
//! that can be found in ``target/release/`` folder once the project is built. Manual for this app
//! can be found in [this documentation](crate::documentation#taxonomy_app)
//!
//! ## Examples
//!
//! ```rust
//! use bioshell_taxonomy::{Taxonomy, Rank};
//! # use std::error::Error;
//! # fn main() -> Result<(), Box<dyn Error>> {
//! // Load from the official NCBI file (downloaded separately)
//! let taxonomy = Taxonomy::load_from_tar_gz("./tests/test_files/test_taxdump.tar.gz")?;
//!
//! // Lookup by species name (case-sensitive)
//! let taxid = taxonomy.taxid("Homo sapiens").unwrap();
//! assert_eq!(taxid, 9606);
//!
//! // Traverse lineage from taxid
//! let lineage = taxonomy.lineage(taxid);
//! for node in &lineage {
//!     println!("{} ({:?})", node.name, node.rank);
//! }
//! // Is a human a vertebrate?
//! assert!(lineage.iter().any(|n| n.name=="Vertebrata"));
//!
//! // Find the node corresponding to a given taxid
//! let species_node = taxonomy.node(9606).unwrap();
//! assert_eq!(species_node.name, "Homo sapiens");
//!
//! // List all known names for the taxid, like "Homo sapiens", "human"
//! for name in taxonomy.names(9606) {
//!     println!("Known name: {}", name);
//! }
//!
//! // List all Kingdoms
//! for kingdom in taxonomy.nodes().filter(|n| n.rank==Rank::Kingdom) {
//!     println!("Kingdom: {}", kingdom.name);
//! }
//! let kingdom_cnt = taxonomy.nodes().filter(|n| n.rank==Rank::Kingdom).count();
//! # assert_eq!(kingdom_cnt, 2);
//! # Ok(())
//! # }
//! ```
//!
//! ## Minimum Supported Rust Version
//! Rust 1.65 or later
//!
//! ## License
//! Licensed under Apache-2.0
pub mod documentation;

mod taxonomy;
pub use taxonomy::*;

mod rank;
pub use rank::*;

mod taxonomy_matcher;
pub use taxonomy_matcher::*;
