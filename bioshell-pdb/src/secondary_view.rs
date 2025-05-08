use crate::{HasCartesians, PDBError, ResidueId, SecondaryStructure, SecondaryStructureTypes, Structure};
use crate::calc::{SubstructureAxis, Vec3};

/// Iterates over [`SecondarySegment`](SecondarySegment)s of a given chain.
///
/// ```
/// # use bioshell_pdb::{Deposit, PDBError, ResidueId, SecondaryStructureTypes};
/// # fn main() -> Result<(), PDBError> {
/// use bioshell_pdb::{assert_delta, SecondaryView};
/// let deposit = Deposit::from_file("./tests/test_files/2gb1.cif")?;
/// let mut strctr = deposit.structure().unwrap();
///
/// let secondary_view = SecondaryView::new(&strctr, "A");
/// assert_eq!(secondary_view.strands().count(), 4);
///
/// let helix = secondary_view.helices().next().unwrap();
/// let axis = helix.ca_axis()?;
/// assert_delta!(axis.length(), 21.9, 0.1);
/// # Ok(())
/// # }
/// ```
pub struct SecondaryView<'a> {
    structure: &'a Structure,
    annotation: SecondaryStructure,
    residues: Vec<&'a ResidueId>, // Cached for slicing later
}

impl<'a> SecondaryView<'a> {
    pub fn new(structure: &'a Structure, chain_id : &str) -> Self {
        let residues: Vec<&'a ResidueId> = structure.residues().into_iter().collect();
        let secondary = structure.secondary(chain_id);
        Self { structure, annotation: secondary, residues }
    }

    /// Returns an iterator over secondary structure segments
    pub fn segments(&'a self) -> Vec<SecondarySegment<'a>> {
        self.annotation
            .ranges()
            .iter()
            .map(move |range| {
                let slice = &self.residues[range.first..=range.last];
                SecondarySegment {
                    kind: range.kind,
                    residues: slice,
                    structure: self.structure,
                }
            })
            .collect()
    }

    /// Provides an iterator over helices of any type
    ///
    pub fn helices(&'a self) -> impl Iterator<Item = SecondarySegment<'a>> {
        self.segments().into_iter().filter(|seg| seg.kind.hec_code() == b'H')
    }

    /// Provides an iterator over all beta strands
    pub fn strands(&'a self) -> impl Iterator<Item = SecondarySegment<'a>> {
        self.segments()
            .into_iter()
            .filter(|seg| matches!(seg.kind, SecondaryStructureTypes::Strand(_)))
    }
}

/// Lists all the residues that belong to a secondary structure segment.
///
/// [`SecondarySegment`](SecondarySegment)s are provided by the [`SecondaryView`](SecondaryView) struct.
#[derive(Clone)]
pub struct SecondarySegment<'a> {
    /// The type of the secondary structure
    pub kind: SecondaryStructureTypes,
    /// Borrowed slice of residue references
    residues: &'a [&'a ResidueId],
    /// The reference to the structure is required to access the coordinates
    structure: &'a Structure,
}

impl<'a> SecondarySegment<'a> {

    /// Returns an immutable slice of the residues in this segment
    pub fn residues(&self) -> &[&'a ResidueId] {
        self.residues
    }

    pub fn ca_axis(&self) -> Result<SubstructureAxis, PDBError> {
        let mut ca: Vec<&Vec3> = vec![];
        for res in self.residues {
            let ca_i = self.structure.atom(res, " CA ")?;
            ca.push(ca_i.position());
        }
        Ok(SubstructureAxis::from_3d_points(&ca))
    }
}

use std::fmt;

impl<'a> fmt::Display for SecondarySegment<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let (Some(first), Some(last)) = (self.residues.first(), self.residues.last()) {
            write!(f, "{}: {} {}", self.kind, first, last)
        } else {
            write!(f, "{}: (empty segment)", self.kind)
        }
    }
}