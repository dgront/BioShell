use crate::SecondaryStructureTypes::{Coil, LeftAlphaHelix, LeftGammaHelix, LeftOmegaHelix, Polyproline, RibbonHelix, Right3_10Helix, RightAlphaHelix, RightGammaHelix, RightOmegaHelix, RightPiHelix, Strand};

/// Defines possible secondary structure types.
///
/// This list is based mostly on the PDB classification and contains all the helical types defined by the PDB standard,
/// extended by a single [`Strand`](SecondaryStructureTypes::Strand) type. Finally, the
/// [`SecondaryStructureTypes`](SecondaryStructureTypes) enum provides also the
/// [`Coil`](Coil) type, which is not a secondary structure type per se,
/// nevertheless is a necessary symbol to provide a secondary structure string for a polypeptide chain.
#[repr(i8)]
#[derive(PartialEq, Debug)]
pub enum SecondaryStructureTypes {
    RightAlphaHelix = 1,
    RightOmegaHelix = 2,
    RightGammaHelix = 3,
    RightPiHelix = 4,
    Right3_10Helix = 5,
    LeftAlphaHelix = 6,
    LeftOmegaHelix = 7,
    LeftGammaHelix = 8,
    RibbonHelix = 9,
    Polyproline = 10,
    Strand = 11,
    Coil = 12,
}


impl SecondaryStructureTypes {
    /// Returns a single byte in the 3-letter `HEC` code that describes this secondary structure element type
    ///
    /// All helical types are converted to 'H' character, while 'E' and 'C' are used to represent a strand
    /// and a coil residue, respectively.
    pub fn hec_code(&self) -> u8 {
        match *self {
            RightAlphaHelix => b'H',
            RightOmegaHelix => b'H',
            RightGammaHelix => b'H',
            RightPiHelix => b'H',
            Right3_10Helix => b'H',
            LeftAlphaHelix => b'H',
            LeftOmegaHelix => b'H',
            LeftGammaHelix => b'H',
            RibbonHelix => b'H',
            Polyproline => b'H',
            Strand => b'E',
            Coil => b'C',
        }
    }

    pub fn from_pdb_class(class: usize) -> SecondaryStructureTypes {
        match class {
            1 => RightAlphaHelix,
            2 => RightOmegaHelix,
            3 => RightPiHelix,
            4 => RightGammaHelix,
            5 => Right3_10Helix,
            6 => LeftAlphaHelix,
            7 => LeftOmegaHelix,
            8 => LeftGammaHelix,
            9 => RibbonHelix,
            10 => Polyproline,
            11 => Strand,
            _ => Coil
        }
    }
}

#[derive(Default, Clone, Debug)]
/// Secondary structure for a polypeptide chain
///
pub struct SecondaryStructure {
    /// a secondary structure is represented as a vector of SecondaryStructureTypes elements
    sec_str: Vec<u8>,
}

impl SecondaryStructure {
    /// Creates a [`SecondaryStructure`](SecondaryStructure) object for a chain of `len` positions.
    ///
    /// All the positions are assigned to  `Coil`
    pub fn new(len: usize) -> SecondaryStructure { SecondaryStructure{ sec_str: vec![12; len] } }

    /// Creates a new [`SecondaryStructure`](SecondaryStructure) object from a vector of [`SecondaryStructureTypes`](SecondaryStructureTypes) codes
    ///
    /// The input vector should comprise byte values corresponding to [`SecondaryStructureTypes`](SecondaryStructureTypes) order numbers
    /// as they are defined by that enum
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{SecondaryStructure, SecondaryStructureTypes};
    /// let codes = vec![12, 12, 1, 1, 1, 1, 1, 1, 12, 12];
    /// let sec_str = SecondaryStructure::from_attrs(codes);
    /// assert_eq!(sec_str.sec_str(0), SecondaryStructureTypes::Coil);
    /// assert_eq!(sec_str.sec_str(2), SecondaryStructureTypes::RightAlphaHelix);
    /// assert_eq!(sec_str.to_string(), String::from("CCHHHHHHCC"));
    /// ```
    pub fn from_attrs(sec_str: Vec<u8>)  -> SecondaryStructure { SecondaryStructure{ sec_str} }

    /// Creates a new [`SecondaryStructure`](SecondaryStructure) object by consuming a given string in the 3-letter `HEC` code
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{SecondaryStructure, SecondaryStructureTypes};
    /// let sec_str = SecondaryStructure::from_hec_string("CCHHHHHHCC");
    /// assert_eq!(sec_str.sec_str(0), SecondaryStructureTypes::Coil);
    /// assert_eq!(sec_str.sec_str(2), SecondaryStructureTypes::RightAlphaHelix);
    /// assert_eq!(sec_str.to_string(), String::from("CCHHHHHHCC"));
    /// ```
    pub fn from_hec_string(sec_str: &str) -> SecondaryStructure {
        let data = sec_str.bytes().map(|e|
            match e {
                b'H' => 1,
                b'E' => 11,
                _ => 12
            }).collect();
        return SecondaryStructure::from_attrs(data);
    }

    /// Provide secondary structure type at a given position of the polypeptide chain
    pub fn sec_str(&self, pos: usize) -> SecondaryStructureTypes {
        SecondaryStructureTypes::from_pdb_class(self.sec_str[pos] as usize)
    }

    /// Return the length of this secondary structure
    pub fn len(&self) -> usize { self.sec_str.len() }

    /// Returns a `HEC` code at a given position in this [`SecondaryStructure`](SecondaryStructure)
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{SecondaryStructure, SecondaryStructureTypes};
    /// let sec_str = SecondaryStructure::from_hec_string("CCHHHHHHCC");
    /// assert_eq!(sec_str.hec_code(0), b'C');
    /// assert_eq!(sec_str.hec_code(2), b'H');
    /// ```
    pub fn hec_code(&self, pos:usize) -> u8 {
        SecondaryStructure::pdb_code_to_hec_code(&self.sec_str[pos])
    }

    /// Creates a string representing this secondary structure.
    ///
    /// The returned string contains only the `H`, `E` and `C` characters
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{SecondaryStructure, SecondaryStructureTypes};
    /// let sec_str = SecondaryStructure::from_hec_string("CCHHHHHHCC");
    /// assert_eq!(sec_str.to_string(), String::from("CCHHHHHHCC"));
    /// ```
    pub fn to_string(&self) -> String {
        String::from_utf8(self.sec_str.iter().map(|e| SecondaryStructure::pdb_code_to_hec_code(e)).collect()).unwrap()
    }

    fn pdb_code_to_hec_code(pdb_code: &u8) -> u8 {
        let sec_str_type = SecondaryStructureTypes::from_pdb_class(*pdb_code  as usize);
        SecondaryStructureTypes::hec_code(&sec_str_type)
    }
}