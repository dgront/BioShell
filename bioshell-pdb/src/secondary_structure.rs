use crate::SecondaryStructureTypes::{Coil, LeftAlphaHelix, LeftGammaHelix, LeftOmegaHelix, Polyproline, RibbonHelix, Right3_10Helix, RightAlphaHelix, RightGammaHelix, RightOmegaHelix, RightPiHelix, Strand};

/// Defines possible secondary structure types.
///
/// This list is based mostly on the PDB classification and contains all the helical types defined by the PDB standard,
/// extended by a single [`Strand`](SecondaryStructureTypes::Strand) type. Finally, the
/// [`SecondaryStructureTypes`](SecondaryStructureTypes) enum provides also the
/// [`Coil`](Coil) type, which is not a secondary structure type per se,
/// nevertheless is a necessary symbol to provide a secondary structure string for a polypeptide chain.
#[repr(u8)]
#[derive(PartialEq, Debug, Clone)]
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
    sec_str: Vec<SecondaryStructureTypes>,
}

impl SecondaryStructure {
    /// Creates a [`SecondaryStructure`](SecondaryStructure) object for a chain of `len` positions.
    ///
    /// All the positions are assigned to  `Coil`
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{SecondaryStructure, SecondaryStructureTypes};
    /// let sec_str = SecondaryStructure::new(10);
    /// assert_eq!(sec_str.len(), 10);
    /// assert_eq!(sec_str.hec_code(5), b'C');
    /// ```
    pub fn new(len: usize) -> SecondaryStructure {
        SecondaryStructure{ sec_str: vec![Coil; len] }
    }

    /// Creates a new [`SecondaryStructure`](SecondaryStructure) object from a vector of [`SecondaryStructureTypes`](SecondaryStructureTypes) codes
    ///
    /// The input vector should contain a respective [`SecondaryStructureTypes`](SecondaryStructureTypes)
    /// enum variant for each residue in the sequence for which the  [`SecondaryStructure`] struct is created.
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{SecondaryStructure, SecondaryStructureTypes};
    /// use bioshell_pdb::SecondaryStructureTypes::{Coil, RightAlphaHelix};
    /// let codes = vec![Coil, Coil, RightAlphaHelix, RightAlphaHelix, RightAlphaHelix, Coil, Coil];
    /// let sec_str = SecondaryStructure::from_attrs(codes);
    /// assert_eq!(sec_str.sec_str(0), Coil);
    /// assert_eq!(sec_str.sec_str(2), RightAlphaHelix);
    /// assert_eq!(sec_str.to_string(), String::from("CCHHHCC"));
    /// ```
    pub fn from_attrs(sec_str: Vec<SecondaryStructureTypes>)  -> SecondaryStructure { SecondaryStructure{ sec_str} }

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
                b'H' => RightAlphaHelix,
                b'E' => Strand,
                _ => Coil
            }).collect();
        return SecondaryStructure::from_attrs(data);
    }

    /// Provide secondary structure type at a given position of the polypeptide chain
    pub fn sec_str(&self, pos: usize) -> SecondaryStructureTypes { self.sec_str[pos].clone() }

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
    pub fn hec_code(&self, pos:usize) -> u8 { self.sec_str[pos].hec_code() }

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
        String::from_utf8(self.sec_str.iter().map(|e| e.hec_code()).collect()).unwrap()
    }

    // fn pdb_code_to_hec_code(pdb_code: &u8) -> u8 {
    //     let sec_str_type = SecondaryStructureTypes::from_pdb_class(*pdb_code  as usize);
    //     SecondaryStructureTypes::hec_code(&sec_str_type)
    // }
}