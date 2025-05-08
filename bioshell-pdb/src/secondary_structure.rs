use crate::SecondaryStructureTypes::{Coil, LeftAlphaHelix, LeftGammaHelix, LeftOmegaHelix, Polyproline, RibbonHelix, Right3_10Helix, RightAlphaHelix, RightGammaHelix, RightOmegaHelix, RightPiHelix, Strand};

/// Defines possible secondary structure types.
///
/// Each enum variant (except for [`Coil`](Coil)) stores the order number (i.e. the index) of a repsective helix or strand.
/// Since the index is stored as `u8` type, the values will run from `0` again when there is more than 255
/// secondary structure elements in a [`Structure`](crate::Structure).
///
/// The enum variants are based mostly on the PDB classification and contains all the helical types defined by the PDB standard,
/// extended by a single [`Strand`](Strand) type. Finally, the
/// [`SecondaryStructureTypes`](SecondaryStructureTypes) enum provides also the
/// [`Coil`](Coil) type, which is not a secondary structure type per se,
/// nevertheless is a necessary symbol to provide a secondary structure string for a polypeptide chain.
#[repr(u8)]
#[derive(PartialEq, Debug, Clone, Copy, Eq)]
pub enum SecondaryStructureTypes {
    RightAlphaHelix(u8),
    RightOmegaHelix(u8),
    RightGammaHelix(u8),
    RightPiHelix(u8),
    Right3_10Helix(u8),
    LeftAlphaHelix(u8),
    LeftOmegaHelix(u8),
    LeftGammaHelix(u8),
    RibbonHelix(u8),
    Polyproline(u8),
    Strand(u8),
    Coil = 12,
}

impl SecondaryStructureTypes {
    /// Returns a single byte in the 3-letter `HEC` code that describes this secondary structure element type
    ///
    /// All helical types are converted to 'H' character, while 'E' and 'C' are used to represent a strand
    /// and a coil residue, respectively.
    pub fn hec_code(&self) -> u8 {
        match *self {
            RightAlphaHelix(_) => b'H',
            RightOmegaHelix(_) => b'H',
            RightGammaHelix(_) => b'H',
            RightPiHelix(_) => b'H',
            Right3_10Helix(_) => b'H',
            LeftAlphaHelix(_) => b'H',
            LeftOmegaHelix(_) => b'H',
            LeftGammaHelix(_) => b'H',
            RibbonHelix(_) => b'H',
            Polyproline(_) => b'H',
            Strand(_) => b'E',
            Coil => b'C',
        }
    }

    /// Assigns the proper [`SecondaryStructureTypes`](SecondaryStructureTypes) enum variant from a PDB class
    ///
    /// # Arguments
    ///  - `class` - PDB class of the secondary structure element: one of the helix classes, a strand or coil; see
    ///  [PDB documentation](https://www.wwpdb.org/documentation/file-format-content/format33/sect5.html#HELIX)
    ///  - `sse_index` - the index (e.g. the order number from 0) of the secondary structure element in the structure
    pub fn from_pdb_class(class: usize, sse_index: u8) -> SecondaryStructureTypes {
        match class {
            1 => RightAlphaHelix(sse_index),
            2 => RightOmegaHelix(sse_index),
            3 => RightPiHelix(sse_index),
            4 => RightGammaHelix(sse_index),
            5 => Right3_10Helix(sse_index),
            6 => LeftAlphaHelix(sse_index),
            7 => LeftOmegaHelix(sse_index),
            8 => LeftGammaHelix(sse_index),
            9 => RibbonHelix(sse_index),
            10 => Polyproline(sse_index),
            11 => Strand(sse_index),
            _ => Coil
        }
    }
}

use std::fmt;
use std::mem::discriminant;

impl fmt::Display for SecondaryStructureTypes {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            RightAlphaHelix(_) => "Right-handed alpha-helix",
            RightOmegaHelix(_) => "Right-handed omega-helix",
            RightGammaHelix(_) => "Right-handed gamma-helix",
            RightPiHelix(_) => "Right-handed pi-helix",
            Right3_10Helix(_) => "Right-handed 3-10 helix",
            LeftAlphaHelix(_) => "Left-handed alpha-helix",
            LeftOmegaHelix(_) => "Left-handed omega-helix",
            LeftGammaHelix(_) => "Left-handed gamma-helix",
            RibbonHelix(_) => "Ribbon helix",
            Polyproline(_) => "Polyproline helix",
            Strand(_) => "Strand",
            Coil => "Coil",
        };
        write!(f, "{}", name)
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
    /// let codes = vec![Coil, Coil, RightAlphaHelix(0), RightAlphaHelix(0), RightAlphaHelix(0), Coil, Coil];
    /// let sec_str = SecondaryStructure::from_attrs(codes);
    /// assert_eq!(sec_str.sec_str(0), Coil);
    /// assert!(matches!(sec_str.sec_str(2), RightAlphaHelix(0)));
    /// assert_eq!(sec_str.to_string(), String::from("CCHHHCC"));
    /// ```
    pub fn from_attrs(sec_str: Vec<SecondaryStructureTypes>)  -> SecondaryStructure { SecondaryStructure{ sec_str} }

    /// Creates a new [`SecondaryStructure`](SecondaryStructure) object by consuming a given string in the 3-letter `HEC` code
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{SecondaryStructure, SecondaryStructureTypes, SecondaryStructureTypes::*};
    /// let result = SecondaryStructure::from_hec_string("CHHHCHHHCEEC");
    /// assert_eq!(result.hec_code(1), b'H');
    /// assert!(matches!(result.sec_str(2), RightAlphaHelix(0)));
    /// # let expected = vec![Coil, RightAlphaHelix(0), RightAlphaHelix(0), RightAlphaHelix(0), Coil,
    /// #     RightAlphaHelix(1), RightAlphaHelix(1), RightAlphaHelix(1), Coil, Strand(2), Strand(2), Coil];
    /// # for (i, sec) in expected.iter().enumerate() {
    /// #     assert_eq!(result.sec_str(i), *sec);
    /// # }
    /// ```
    pub fn from_hec_string(sec_str: &str) -> SecondaryStructure {
        let mut data = Vec::with_capacity(sec_str.len());
        let mut current_char: Option<u8> = None;
        let mut segment_index: u8 = 0;

        for &byte in sec_str.as_bytes() {
            match byte {
                b'H' | b'E' => {
                    if Some(byte) != current_char {
                        current_char = Some(byte);
                        segment_index += 1;
                    }

                    let ss = match byte {
                        b'H' => RightAlphaHelix(segment_index - 1),
                        b'E' => Strand(segment_index - 1),
                        _ => unreachable!(),
                    };

                    data.push(ss);
                }
                _ => {
                    current_char = None;
                    data.push(Coil);
                }
            }
        }

        SecondaryStructure::from_attrs(data)
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

    /// Breaks a secondary structure into secondary structure elements (SSEs).
    ///
    /// Each SSE is defined as a block of consecutive residues with the same secondary structure type. Therefore,
    /// when an alpha-helix is directly followed by a pi-helix, they are considered as two different SSEs,
    /// even though they both are marked as `'H'` in the `HEC` code.
    ///
    /// # Examples
    ///
    /// Loops are recognized as SSEs by this method, so the string `"CCHHHHHHCCHHHHC"` will be broken into 5 ranges:
    /// ```
    /// use bioshell_pdb::{SecondaryStructure, SecondaryRange, SecondaryStructureTypes};
    /// let sec_str = SecondaryStructure::from_hec_string("CCHHHHHHCCHHHHC");
    /// let blocks = sec_str.ranges();
    /// assert_eq!(blocks.len(), 5);
    /// assert_eq!(blocks[1], SecondaryRange::new(SecondaryStructureTypes::RightAlphaHelix(0), 2, 7));
    /// assert_eq!(blocks[4].kind, SecondaryStructureTypes::Coil);
    /// ```
    /// Segments of the same type but different index are split into multiple ranges:
    /// ```
    /// use bioshell_pdb::{SecondaryStructure, SecondaryRange, SecondaryStructureTypes, SecondaryStructureTypes::*};
    /// let sec_seq = vec![Coil, RightAlphaHelix(0), RightAlphaHelix(1), Coil];
    /// let sec_str = SecondaryStructure::from_attrs(sec_seq);
    /// let blocks = sec_str.ranges();
    /// assert_eq!(blocks.len(), 4);
    /// assert_eq!(&sec_str.to_string(), "CHHC");
    /// ```
    pub fn ranges(&self) -> Vec<SecondaryRange> {
        let mut result = Vec::new();

        if self.sec_str.is_empty() { return result; }

        let mut start = 0;
        let mut current_val = self.sec_str[0].clone();

        for (i, &val) in self.sec_str.iter().enumerate().skip(1) {
            if !Self::same_segment(&current_val, &val) {
                result.push(SecondaryRange::new(current_val, start, i - 1));
                current_val = val;
                start = i;
            }
        }

        // Push the final group
        result.push(SecondaryRange::new(current_val, start, self.sec_str.len() - 1));

        result
    }
    fn same_segment(a: &SecondaryStructureTypes, b: &SecondaryStructureTypes) -> bool {

        if discriminant(a) != discriminant(b) {
            return false;
        }
        match (a, b) {
            (Coil, Coil) => true,
            _ => a == b, // same variant type, so compare payloads
        }
    }
}

/// A range of residues spanning a single Secondary Structure Element: a helix, strand, etc
///
/// This struct provides only the boundaries of the segment in a respective proteins chain, not the residues themselves.
/// A vector of such ranges can be obtained by calling the [`SecondaryStructure::ranges()`](SecondaryStructure::ranges) method.
///
/// To access the [`Residues`](crate::ResidueId)s of a segment, use the [`SecondarySegment`](crate::SecondarySegment),
/// provided by [`SecondaryView`](crate::SecondaryView) struct.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SecondaryRange {
    /// The type of the secondary structure
    pub kind: SecondaryStructureTypes,
    /// The first residue of the segment
    pub first: usize,
    /// The last residue of the segment (inclusive)
    pub last: usize
}

impl SecondaryRange {
    pub fn new(kind: SecondaryStructureTypes, first: usize, last: usize) -> SecondaryRange {
        SecondaryRange { kind, first, last }
    }
}

