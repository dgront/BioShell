use crate::Element;

/// A lightweight representation of a chemical atom.
///
/// The struct stores an internal atom index, atomic number, and formal charge.
/// It does not store coordinates, atom names or residue information.
///
/// # Examples
///
/// ```
/// # use bioshell_chem::{Atom, Element};
/// let carbon = Atom::neutral(0, Element::C);
/// let charged_oxygen = Atom::charged(1, Element::O, -1);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Atom {
    /// Internal atom index.
    ///
    /// This index can be used to preserve the original atom numbering
    /// from an input file or to identify the atom inside a molecule.
    index: usize,
    element: Element,
    formal_charge: i8,
}


impl Atom {
    /// Creates a new atom with an explicit formal charge.
    ///
    /// # Examples
    ///
    /// ```
    /// # use bioshell_chem::{Atom, Element};
    /// let carbon = Atom::charged(0, Element::C, 0);
    /// let charged_oxygen = Atom::charged(1, Element::O, -1);
    /// ```
    pub fn charged(index: usize, element: Element, formal_charge: i8) -> Self {
        Self { index, element, formal_charge}
    }

    /// Creates a neutral atom.
    ///
    /// # Examples
    /// ```
    /// use bioshell_chem::{Atom, Element};
    /// let carbon = Atom::neutral(0, Element::C);
    /// ```
    pub fn neutral(index: usize, element: Element) -> Self {
        Self { index, element, formal_charge: 0}
    }

    /// Returns `true` if this atom is a hydrogen.
    ///
    /// Note that we consider an atom to be *heavy* if it's not a hydrogen.
    pub fn if_hydrogen(self) -> bool {
        self.element == Element::H
    }

    /// Returns the internal atom index.
    pub fn index(&self) -> usize { self.index }

    /// Returns the element type of this atom.
    pub fn element(&self) -> Element { self.element }

    /// Returns the atomic number.
    pub fn atomic_number(&self) -> u8 { self.element.atomic_number() }

    /// Returns the formal charge.
    pub fn formal_charge(&self) -> i8 {
        self.formal_charge
    }

    /// Sets the formal charge.
    pub fn set_formal_charge(&mut self, formal_charge: i8) {
        self.formal_charge = formal_charge;
    }
}