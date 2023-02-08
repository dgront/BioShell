use std::ops::Range;

/// Defines the way how a system's energy is evaluated in a BioShell's simulation.
///
pub trait Energy<S> {
    /// Evaluates the total energy of a given system, taking into account all its components
    fn energy(&self, system: &S) -> f64;

    /// Evaluates the energy of a single component (atom, spin, residue) of a given system
    fn energy_by_pos(&self, system: &S, pos: usize) -> f64;

    /// Evaluates the energy of a contiguous range of system's components.
    fn energy_by_range(&self, system: &S, range: &Range<usize>) -> f64;

    /// Returns the name of this energy function.
    /// The returned name may be used to identify this energy, e.g. to name a column in a score table
    fn name(&self) -> String;
}

/// Defines the basic properties of a simulated system.
///
/// BioShell's approach to molecular simulations assumes each modelled system comprises a number
/// of interaction centers: spins, atoms, etc.
pub trait System: Clone {
    /// Returns the current size of the modelled system
    fn size(&self) -> usize;

    /// Changes the current size of the modelled system
    fn set_size(&mut self, new_size: usize);

    /// Returns the maximum size of the modelled system
    fn capacity(&self) -> usize;

    /// Sets i-th component of this system by copying its DOFs from a given `rhs` object
    fn copy_from(&mut self, i:usize, rhs: &Self);
}
