use crate::calc::Vec3;

/// A trait for objects that provide access to 3D Cartesian coordinates as a `Vec3`.
///
/// This allows generalizing algorithms over types like `Vec3`, `PdbAtom`, etc.
pub trait HasCartesians {
    /// Returns a reference to the Cartesian coordinate represented by a `Vec3`.
    fn position(&self) -> &Vec3;
}
