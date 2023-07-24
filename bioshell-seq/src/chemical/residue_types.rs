use std::collections::HashMap;
use std::str::SplitWhitespace;

/// Defines the types of monomers - residue types that are biomolecular building blocks.
///
/// The set of possible types is a shortened version of the full list of possibilities found in PDB files.
/// Each of these types has a single-letter code assigned (such as ``'P'`` for proteins and ``'S'`` for sacharides);
/// these codes are used by [`try_from()`](ResidueType::try_from) method to create a [`ResidueType`](ResidueType)
#[derive(Debug, Ord, PartialOrd, Eq, PartialEq, Clone)]
#[repr(i8)]
pub enum MonomerType {
    /// standard amino acid, single-character code: ``'P'``
    PeptideLinking = 0,
    /// standard deoxy-nucleotide, single-character code: ``'D'``
    DNALinking,
    /// standard nucleotide, single-character code: ``'R'``
    RNALinking,
    /// a sugar residue, single-character code: ``'S'``
    Sacharide,
    /// most likely a ligand, single-character code: ``'N'``
    NonPolymer,
    /// anything else, single-character code: ``'O'``
    OTHER
}

impl TryFrom<char> for MonomerType {
    type Error = String;

    /// Returns a ``MonomerType`` for its one-letter code
    ///
    /// # Example
    /// ```rust
    /// use bioshell_seq::chemical::MonomerType;
    /// assert_eq!(MonomerType::try_from('P').unwrap(), MonomerType::PeptideLinking);
    /// ```
    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            'P' => Ok(MonomerType::PeptideLinking),
            'D' => Ok(MonomerType::DNALinking),
            'R' => Ok(MonomerType::RNALinking),
            'S' => Ok(MonomerType::Sacharide),
            'N' => Ok(MonomerType::NonPolymer),
            'O' => Ok(MonomerType::OTHER),
            _ => {Err(format!("Can't find a MonomerType for the one-letter code: {}!", value))}
        }
    }
}

/// Provides basic information for a given residue type found in a biomolecule
pub trait ResidueTypeProperties {

    /// single-letter code of this monomer
    fn code1(&self) -> char;

    /// integer ID  of this monomer
    fn id(&self) -> u16;

    /// three-letter code  of this monomer
    fn code3(&self) -> String;

    /// chemical type of this monomer
    fn chem_compound_type(&self) -> MonomerType;
}

/// Defines a residue type, i.e. a small molecule or its fragment that can be found in biomolecular data
///
/// Currently there is nearly 39 thousand residue types (_aka_ monomers) found in PDB files, all listed in
/// this [Components-pub.cif](http://ligand-expo.rcsb.org/dictionaries/Components-pub.cif) file.
///
/// The purpose of the ``ResidueType`` is to provide conversion between these non-standard residue types
/// and their _canonical_ variants they originate from. Note that for majority of cases such a conversion
/// is not possible and the ``UNK`` standard residue (``X``) is used as the _canonical_ variant.
///
/// # Example
/// ```rust
/// use bioshell_seq::chemical::{MonomerType, ResidueType, ResidueTypeManager, StandardResidueType};
/// let aln = ResidueType::try_from(String::from("ALN A P")).unwrap();
/// assert_eq!(aln.code3, String::from("ALN"));
/// assert_eq!(aln.parent_type, StandardResidueType::ALA);
/// assert_eq!(aln.chem_compound_type, MonomerType::PeptideLinking);
/// ```
#[derive(Debug, Clone)]
pub struct ResidueType {
    /// three-letter code of this ``ResidueType``, such as ``"ALA"`` for alanine
    pub code3: String,
    /// _"standard"_ residue that this residue type is a variant of
    pub parent_type: StandardResidueType,
    /// chemical type of this residue type says whether a monomer can form a protein or nucleic acid
    pub chem_compound_type: MonomerType,
}

impl TryFrom<String> for ResidueType {
    type Error = &'static str;

    /// Parses a string into a ``ResidueType`` instance.
    ///
    /// Expected string format: ``"ALN A P"``, where ``"ALN"`` is a three-letter code (
    /// in this example its _NAPHTHALEN-2-YL-3-ALANINE_, represented by ``"ALN"``),
    /// ``'A'`` is the one-letter code of its parent standard residue type (which is alanine)
    /// and ``'P'`` encodes is its monomer type (``MonomerType::PeptideLinking`` in this case)
    ///
    /// # Examples
    /// ```rust
    /// use bioshell_seq::chemical::ResidueType;
    /// let aln = ResidueType::try_from(String::from("ALN A P")).unwrap();
    /// ```
    fn try_from(value: String) -> Result<Self, Self::Error> {
        let mut tokens: SplitWhitespace = value.split_whitespace();
        let code3 = String::from(tokens.next().unwrap());
        let code1: char = tokens.next().unwrap().chars().next().unwrap();
        let type1: char = tokens.next().unwrap().chars().next().unwrap();
        let parent_type = StandardResidueType::try_from(code1).unwrap();
        let chem_compound_type = MonomerType::try_from(type1).unwrap();
        return Ok(ResidueType { code3, parent_type, chem_compound_type });
    }
}

/// Provides unique integer ID for each registered ResidueType
///
/// The standard residue types (amino aicds, bases, etc) are registered by default. User can add
/// additional, non-standard residue types into a manager.
///
/// # Examples
/// ```rust
/// use bioshell_seq::chemical::{ResidueType, ResidueTypeManager};
///
/// let mut mgr = ResidueTypeManager::new();
/// // --- This should pass, as all standard residue types (including alanine) are preloaded by a constructor
/// let ala = mgr.by_code3(&String::from("ALA"));
/// assert!(ala.is_some());
/// // --- There are 29 standard residue types
/// assert_eq!(mgr.count(), 29);
/// // --- ALN hasn't been inserted yet
/// assert!(mgr.by_code3(&String::from("ALN")).is_none());
/// let aln = ResidueType::try_from(String::from("ALN A P")).unwrap();
/// mgr.register_residue_type(aln);
/// assert!(mgr.by_code3(&String::from("ALN")).is_some());
/// // --- ALN residue type has been registered at index the 29
/// assert_eq!(mgr.index(&String::from("ALN")).unwrap(), 29);
/// ```
pub struct ResidueTypeManager {
    registered_types: Vec<ResidueType>,
    by_code_3: HashMap<String, usize>
}

impl ResidueTypeManager {
    /// Creates a new residue type manager
    ///
    /// ``ResidueType`` objects corresponding to standard StandardResidueType enum values are
    /// automatically created and registered in this manager
    pub fn new() ->  ResidueTypeManager{
        let mut out = ResidueTypeManager{registered_types: vec![], by_code_3: HashMap::new()};

        for srt in StandardResidueType::TYPES {
            let rt = ResidueType{code3: srt.code3(), parent_type: srt,
                chem_compound_type: srt.chem_compound_type()};
            out.by_code_3.insert(srt.code3(), out.registered_types.len());
            out.registered_types.push(rt);
        }
        return out;
    }

    /// Counts the residue types registered in this manager
    pub fn count(&self) -> usize { self.registered_types.len() }

    /// Register a new residue type in this manager.
    ///
    /// If the monomer (identified by its three-letter code) already exists in this manager,
    /// it will not be inserted, i.e. this method does not replace an old monomer with a new one.
    pub fn register_residue_type(&mut self, res_type: ResidueType) -> bool {
        let key = &res_type.code3;
        if self.by_code_3.contains_key(key) { return false; }
        self.by_code_3.insert(key.clone(), self.registered_types.len());
        self.registered_types.push(res_type);
        return true;
    }

    /// Provides a monomer for a given three-letter code.
    ///
    /// Returns an option that contains the monomer or ``None`` if it hasn't been registered
    ///
    /// # Examples
    /// ```rust
    /// use bioshell_seq::chemical::{ResidueType, ResidueTypeManager};
    /// let aln = ResidueType::try_from(String::from("ALN A P")).unwrap();
    /// let mut mgr = ResidueTypeManager::new();
    /// // This should pass, as all standard residue types are preloaded by a constructor
    /// let ala = mgr.by_code3(&String::from("ALA"));
    /// assert!(ala.is_some());
    /// // --- ALN hasn't been inserted yet
    /// assert!(mgr.by_code3(&String::from("ALN")).is_none());
    /// mgr.register_residue_type(aln);
    /// assert!(mgr.by_code3(&String::from("ALN")).is_some());
    /// ```
    pub fn by_code3(&self, code3: &String) -> Option<&ResidueType> {
        if self.by_code_3.contains_key(code3) {
            return Some(&self.registered_types[self.by_code_3[code3]]);
        } else {
            return None;
        }
    }

    /// Returns an integer index identifying a residue type.
    ///
    /// Returns ``None`` if the residue type hasn't been registered by this manager.
    pub fn index(&self, code3: &String) -> Option<usize> {
        if self.by_code_3.contains_key(code3) {
            return Some(self.by_code_3[code3]);
        } else {
            return None;
        }
    }
}


macro_rules! define_res_types {

    ($($name:ident $id:literal $one_letter_code:literal $three_letters_code:literal $res_type:ident),*) => {

        /// Lists all the 20 standard amino acids, 5 standard DNA/RNA building blocks and a few _special_ types
        #[derive(Debug, Ord, PartialOrd, Eq, PartialEq, Clone, Copy)]
        #[repr(i16)]
        #[allow(non_camel_case_types)]
        pub enum StandardResidueType {
            $(
                $name = $id,
            )*
        }

        impl ResidueTypeProperties for StandardResidueType {

            fn code1(&self) -> char {
                match *self {
                    $(
                        Self::$name => $one_letter_code
                    ),*
                }
            }

            fn id(&self) -> u16 {
                match *self {
                    $(
                        Self::$name => (Self::$name) as u16,
                    )*
                }
            }

            fn code3(&self) -> String {
                match *self {
                    $(
                        Self::$name => String::from($three_letters_code),
                    )*
                }
            }

            fn chem_compound_type(&self) -> MonomerType {
                match *self {
                    $(
                        Self::$name => MonomerType::$res_type,
                    )*
                }
            }
        }

        impl TryFrom<char> for StandardResidueType {
            type Error = String;

            /// Returns a [`StandardResidueType`](StandardResidueType) enum for a given one-letter code
            fn try_from(value: char) -> Result<Self, Self::Error> {
                match value {
                    $(
                        $one_letter_code => Ok(StandardResidueType::$name),
                    )*
                    _ => {Err(format!("Can't find an amino acid for the one-letter code: {}!", value))}
                }
            }
        }

        impl StandardResidueType {
            /// contains all standard residue types, gathered in an array to iterate over them easily
            ///
            /// # Examples
            /// ``` rust
            /// use bioshell_seq::chemical::{StandardResidueType, ResidueTypeProperties, MonomerType};
            ///
            /// // ---------- Iterate over standard amino acid enum types
            /// let mut n_aa: i8 = 0;
            /// for srt in StandardResidueType::TYPES {
            ///     if srt.chem_compound_type() == MonomerType::PeptideLinking { n_aa += 1; }
            /// }
            /// assert_eq!(n_aa, 21);       // 20 standard amino acids + UNK
            /// ```
            pub const TYPES: [Self; 29] = [
            $(
                Self::$name,
            )*
            ];
        }
    };
}

define_res_types! {
    ALA 0 'A' "ALA" PeptideLinking,
    ARG 1 'R' "ARG" PeptideLinking,
    ASN 2 'N' "ASN" PeptideLinking,
    ASP 3 'D' "ASP" PeptideLinking,
    CYS 4 'C' "CYS" PeptideLinking,
    GLN 5 'Q' "GLN" PeptideLinking,
    GLU 6 'E' "GLU" PeptideLinking,
    GLY 7 'G' "GLY" PeptideLinking,
    HIS 8 'H' "HIS" PeptideLinking,
    ILE 9 'I' "ILE" PeptideLinking,
    LEU 10 'L' "LEU" PeptideLinking,
    LYS 11 'K' "LYS" PeptideLinking,
    MET 12 'M' "MET" PeptideLinking,
    PHE 13 'F' "PHE" PeptideLinking,
    PRO 14 'P' "PRO" PeptideLinking,
    SER 15 'S' "SER" PeptideLinking,
    THR 16 'T' "THR" PeptideLinking,
    TRP 17 'W' "TRP" PeptideLinking,
    TYR 18 'Y' "TYR" PeptideLinking,
    VAL 19 'V' "VAL" PeptideLinking,
    UNK 20 'X' "UNK" PeptideLinking,
    A 21 'a' "A" DNALinking,
    C 22 'c' "C" DNALinking,
    G 23 'g' "G" DNALinking,
    T 24 't' "T" DNALinking,
    U 25 'u' "U" RNALinking,
    GAP 26 '-' "GAP" OTHER,
    GPE 27 '_' "GPE" OTHER,
    UNL 28 'Z' "UNL" OTHER
}
