use std::fmt;
use std::str::FromStr;

/// Defines the types of chemical bonds.
///
/// A bond type, internally represented as a `u8`, can be converted to and from a string
/// as well as to a MOL2 code.
///
/// # Example
/// ```
/// use bioshell_chem::BondType;
/// let b = BondType::from_code("DOUB");
/// assert_eq!(b, BondType::Double);
/// assert_eq!(b.mol2_code(), "2");
/// assert_eq!(b.name(), "double");
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum BondType {
    /// A dummy bond, most likely not a chemical bond.
    Dummy = 0,
    /// A single bond.
    Single = 1,
    /// A double bond.
    Double = 2,
    /// A triple bond.
    Triple = 3,
    /// An aromatic bond.
    Aromatic = 4,
    /// A disulfide bond.
    Disulfide = 5,
    /// An amide bond.
    Amide = 6,
    /// A hydrogen bond.
    Hydrogen = 7,
    /// Any other chemical bond.
    Unknown = 8,
}

impl BondType {
    /// Returns a string denoting a given bond type, as in MOL2 format.
    pub fn mol2_code(self) -> &'static str {
        match self {
            BondType::Single => "1",
            BondType::Double => "2",
            BondType::Triple => "3",
            BondType::Aromatic => "ar",
            BondType::Hydrogen => "hy",
            BondType::Dummy => "du",
            BondType::Disulfide => "1",
            BondType::Amide => "am",
            BondType::Unknown => "un",
        }
    }

    /// Returns a human-readable bond name.
    pub fn name(self) -> &'static str {
        match self {
            BondType::Dummy => "dummy",
            BondType::Single => "single",
            BondType::Double => "double",
            BondType::Triple => "triple",
            BondType::Aromatic => "aromatic",
            BondType::Disulfide => "disulfide",
            BondType::Amide => "amide",
            BondType::Hydrogen => "hydrogen",
            BondType::Unknown => "unknown",
        }
    }

    /// Returns a bond type for a given string code.
    pub fn from_code(code: &str) -> Self {
        match code.trim().to_ascii_uppercase().as_str() {
            "1" | "SING" | "SINGLE" => BondType::Single,
            "2" | "DOUB" | "DOUBLE" => BondType::Double,
            "3" | "TRIP" | "TRIPLE" => BondType::Triple,
            "AR" | "AROM" | "AROMATIC" => BondType::Aromatic,
            "DU" | "DUMMY" => BondType::Dummy,
            "HY" | "HYDROGEN" => BondType::Hydrogen,
            "AM" | "AMIDE" => BondType::Amide,
            "SS" | "DISULFIDE" => BondType::Disulfide,
            _ => BondType::Unknown,
        }
    }

    /// Returns the bond order as a floating point value.
    ///
    /// Single and double bonds have orders of 1.0 and 2.0, as expected;
    /// aromatic and amide bonds have an order of 1.5.
    pub fn order(self) -> f32 {
        match self {
            BondType::Dummy => 0.0,
            BondType::Single => 1.0,
            BondType::Double => 2.0,
            BondType::Triple => 3.0,
            BondType::Aromatic => 1.5,
            BondType::Disulfide => 1.0,
            BondType::Amide => 1.5,
            BondType::Hydrogen => 0.0,
            BondType::Unknown => 0.0,
        }
    }
}

impl FromStr for BondType {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(BondType::from_code(s))
    }
}

impl fmt::Display for BondType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}