/// A lightweight representation of a chemical atom.
///
/// The struct stores an internal atom index, atomic number, and formal charge.
/// It does not store coordinates, atom names or residue information.
///
/// # Examples
///
/// ```
/// # use bioshell_chem::Atom;
/// let carbon = Atom::neutral(0, 6);
/// let charged_oxygen = Atom::charged(1, 8, -1);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Atom {
    /// Internal atom index.
    ///
    /// This index can be used to preserve the original atom numbering
    /// from an input file or to identify the atom inside a molecule.
    index: usize,
    atomic_number: u8,
    formal_charge: i8,
}


impl Atom {
    /// Creates a new atom with an explicit formal charge.
    ///
    /// # Examples
    ///
    /// ```
    /// # use bioshell_chem::Atom;
    /// let carbon = Atom::charged(0, 6, 0);
    /// let charged_oxygen = Atom::charged(1, 8, -1);
    /// ```
    pub fn charged(index: usize, atomic_number: u8, formal_charge: i8) -> Self {
        Self { index, atomic_number, formal_charge}
    }

    /// Creates a neutral atom.
    ///
    /// This is a convenience constructor equivalent to calling
    /// `Atom::new(index, atomic_number, 0)`.
    ///
    /// # Examples
    ///
    /// ```
    /// # use bioshell_chem::Atom;
    /// let carbon = Atom::neutral(0, 6);
    /// ```
    pub fn neutral(index: usize, atomic_number: u8) -> Self {
        Self { index, atomic_number, formal_charge: 0}
    }

    /// Returns the internal atom index.
    pub fn index(&self) -> usize { self.index }

    /// Returns the atomic number.
    pub fn atomic_number(&self) -> u8 { self.atomic_number }

    /// Returns the formal charge.
    pub fn formal_charge(&self) -> i8 {
        self.formal_charge
    }

    /// Sets the formal charge.
    pub fn set_formal_charge(&mut self, formal_charge: i8) {
        self.formal_charge = formal_charge;
    }
}