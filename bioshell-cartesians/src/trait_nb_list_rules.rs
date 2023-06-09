use crate::{Coordinates};

/// Rules that define which atoms and which atom pairs will be excluded from hashing
///
/// This trait defines exclusion rules for a neighbor list structs, such as  [`NbList`](NbList).
/// By default, a neighbor list provides all spatial neighbors of a given atom for efficient
/// evaluation of pairwise interactions.  An object derived from this [`NbListRules`] trait asks
/// a neighbor list to omit some of them.
///
/// For example, atoms that are directly connected with a covalent bond are typically excluded
/// from non-bonded energy evaluation. ``NbListRules::if_pair_excluded(i, j)`` should return
/// ``true`` in such cases. To exclude a given atom `i` from any energy evaluation,
/// ``NbListRules::if_atom_excluded(i)`` can be used.
pub trait NbListRules {
    /// Says if an atom is excluded from any interactions
    fn if_atom_excluded(&self, coordinates: &Coordinates, i_atom: usize) -> bool;

    /// Says if a given pair of atoms is excluded from any interactions
    fn if_pair_excluded(&self, coordinates: &Coordinates, i_atom: usize, j_atom: usize) -> bool;

    /// Each NbListRules must provide a way to clone its boxed instance
    fn box_clone(&self) -> Box<dyn NbListRules>;
}