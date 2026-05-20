use bioshell_core::Vec3;
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
#[derive(Debug, Clone, Copy)]
pub struct Atom {
    /// Internal atom index.
    index: usize,
    /// coordinates of this atom
    pos: Vec3,
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
        Self { index, pos: Default::default(), element, formal_charge}
    }

    /// Creates a neutral atom.
    ///
    /// # Examples
    /// ```
    /// use bioshell_chem::{Atom, Element};
    /// let carbon = Atom::neutral(0, Element::C);
    /// ```
    pub fn neutral(index: usize, element: Element) -> Self {
        Self { index, pos: Default::default(), element, formal_charge: 0}
    }

    pub fn x(&self) -> f64 { self.pos.x }
    pub fn y(&self) -> f64 { self.pos.y }
    pub fn z(&self) -> f64 { self.pos.z }

    /// Provides read-only access to the coordinates of this atom
    pub fn pos(&self) -> &Vec3 { &self.pos }

    /// Provides mutable access to the coordinates of this atom
    pub fn pos_mut(&mut self) -> &mut Vec3 { &mut self.pos }

    /// Set the new vector of coordinates for this atom
    pub fn set_pos(mut self, new_pos: &Vec3)  {
        self.pos.set(new_pos);
    }

    /// Set new coordinates for this atom
    pub fn set_pos3(&mut self, x:f64, y:f64, z:f64)  {
        self.pos.x = x;
        self.pos.y = y;
        self.pos.z = z;
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
    pub fn charge(&self) -> i8 {
        self.formal_charge
    }

    /// Sets the formal charge.
    pub fn set_charge(&mut self, formal_charge: i8) {
        self.formal_charge = formal_charge;
    }
}