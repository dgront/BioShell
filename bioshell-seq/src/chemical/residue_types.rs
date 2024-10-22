use std::collections::HashMap;
use std::sync::{Mutex, MutexGuard};
use clap::__macro_refs::once_cell::sync::Lazy;

/// Defines the types of monomers - residue types that are biomolecular building blocks.
///
/// The set of possible types has been taken from [*this list*](https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp.type.html):
/// This list will most likely grow over time and it is not recommended to exhaustively match against it.
#[derive(Debug, Ord, PartialOrd, Eq, PartialEq, Clone)]
#[repr(i8)]
pub enum MonomerType {
    // Peptide-related
    PeptideLinking = 0,
    LPeptideLinking,
    LPeptideCOOH,
    LPeptideNH3,
    DBetaPeptide,
    CGammaLinking,
    DGammaPeptide,
    CDeltaLinking,
    DPeptideCOOH,
    DPeptideNH3,
    DPeptideLinking,
    LBetaPeptideCGammaLinking,
    LGammaPeptideCDeltaLinking,
    PeptideLike,

    // RNA-related
    RNALinking,
    LRNALinking,
    RNAOH3PrimeTerminus,
    RNAOH5PrimeTerminus,

    // DNA-related
    DNALinking,
    LDNALinking,
    DNAOH3PrimeTerminus,
    DNAOH5PrimeTerminus,

    // Saccharides
    Saccharide,
    LSaccharide,
    LSaccharideAlphaLinking,
    LSaccharideBetaLinking,
    DSaccharide,
    DSaccharideAlphaLinking,
    DSaccharideBetaLinking,

    // Other
    NonPolymer,
    Other,
}

impl MonomerType {
    /// Returns `true` if the monomer can form a peptide bond
    ///
    /// # Examples
    /// ``` rust
    /// use bioshell_seq::chemical::MonomerType::{DSaccharide, LPeptideLinking, PeptideLinking};
    /// assert!(PeptideLinking.is_peptide_linking());
    /// assert!(LPeptideLinking.is_peptide_linking());
    /// assert!(! DSaccharide.is_peptide_linking());
    /// ```
    pub fn is_peptide_linking(&self) -> bool {
        matches!(self, MonomerType::PeptideLinking | MonomerType::LPeptideLinking | MonomerType::LPeptideCOOH | MonomerType::LPeptideNH3 | MonomerType::DPeptideCOOH | MonomerType::DPeptideNH3 | MonomerType::DPeptideLinking | MonomerType::PeptideLike)
    }
    /// Returns `true` if the monomer can form a nucleic acid
    ///
    /// # Examples
    /// ``` rust
    /// use bioshell_seq::chemical::MonomerType::{LRNALinking, RNALinking, DNALinking, DSaccharide};
    /// assert!(LRNALinking.is_peptide_linking());
    /// assert!(RNALinking.is_peptide_linking());
    /// assert!(DNALinking.is_peptide_linking());
    /// assert!(! DSaccharide.is_peptide_linking());
    /// ```
    pub fn is_nucleic_linking(&self) -> bool {
        matches!(self, MonomerType::RNALinking | MonomerType::LRNALinking
            | MonomerType::DNALinking | MonomerType::LDNALinking
            | MonomerType::DNAOH3PrimeTerminus | MonomerType::DNAOH5PrimeTerminus
            | MonomerType::RNAOH3PrimeTerminus | MonomerType::RNAOH5PrimeTerminus)
    }
}

impl TryFrom<&str> for MonomerType {
    type Error = String;

    /// Returns a [`MonomerType`](MonomerType) enum for a given string
    ///
    /// The input string is trimmed, whitespace and dashes are removed, and the string is converted to uppercase
    /// to prevent errors.
    /// # Examples
    /// ``` rust
    /// use bioshell_seq::chemical::MonomerType;
    /// assert!(MonomerType::try_from("L-peptide linking").is_ok());
    /// assert!(MonomerType::try_from("'L-peptide linking'").is_ok());
    /// assert!(MonomerType::try_from("'L-saccharide, alpha linking'").is_ok());
    /// ```
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        // Remove whitespace and capitalize the input string
        let cleaned_value: String = value.trim().to_ascii_uppercase().chars()
            .filter(|c| !c.is_whitespace() && *c != ',' && *c != '-' && *c != '\'').collect();

        // Match the cleaned string to the enum variants
        if let Some(monomer) = MONOMER_MAP.get(cleaned_value.as_str()) {
            return Ok(monomer.clone())
        } else {
            for (name, monomer) in MONOMER_MAP.iter() {
                if cleaned_value.contains(name) {
                    return Ok(monomer.clone())
                }
            }
        }
        Err(format!("Unknown chemical: {}", value))
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
/// let aln = ResidueType::from_attrs("ALN", StandardResidueType::ALA, MonomerType::PeptideLinking);
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

impl ResidueType {
    pub fn from_attrs(code3: &str, parent_type: StandardResidueType, chem_compound_type: MonomerType) -> Self {
        ResidueType { code3: code3.to_string(), parent_type, chem_compound_type, }
    }
}

impl ResidueTypeProperties for ResidueType {
    fn code1(&self) -> char { self.parent_type.code1() }
    fn id(&self) -> u16 { self.parent_type.id() }
    fn code3(&self) -> String { self.code3.clone() }
    fn chem_compound_type(&self) -> MonomerType { self.chem_compound_type.clone() }
}

/// Provides unique integer ID for each registered ResidueType
///
/// The standard residue types (amino acids, bases, etc) are registered by default. User can add
/// additional, non-standard residue types into a manager.
///
/// # Examples
/// ```rust
/// use bioshell_seq::chemical::{MonomerType, ResidueType, ResidueTypeManager, StandardResidueType};
///
/// let mut mgr = ResidueTypeManager::get();
/// // --- This should pass, as all standard residue types (including alanine) are preloaded by a constructor
/// let ala = mgr.by_code3(&String::from("ALA"));
/// assert!(ala.is_some());
/// // --- There are 33 standard residue types
/// assert_eq!(mgr.count(), 33);
/// // --- ALN hasn't been inserted yet
/// assert!(mgr.by_code3(&String::from("ALN")).is_none());
/// let aln = ResidueType::from_attrs("ALN", StandardResidueType::ALA, MonomerType::PeptideLinking);
/// mgr.register_residue_type(aln);
/// assert!(mgr.by_code3(&String::from("ALN")).is_some());
/// // --- ALN residue type has been registered at index the 33
/// assert_eq!(mgr.index(&String::from("ALN")).unwrap(), 33);
/// ```
pub struct ResidueTypeManager {
    registered_types: Vec<ResidueType>,
    by_code_3: HashMap<String, usize>
}

impl ResidueTypeManager {
    /// Creates a new residue type manager.
    ///
    /// ResidueTypeManager implements the singleton pattern: there is only one instance of this struct
    /// allowed, created as a static object. Simply put: the [`new()`] method is private; you should
    /// use [`KnownResidueTypes`] object instead of creating a new one.
    ///
    /// ``ResidueType`` objects corresponding to standard StandardResidueType enum values are
    /// automatically created and registered in this manager
    pub(crate) fn new() ->  ResidueTypeManager {
        let mut out = ResidueTypeManager{registered_types: vec![], by_code_3: HashMap::new()};

        for srt in StandardResidueType::TYPES {
            let rt = ResidueType{code3: srt.code3(), parent_type: srt,
                chem_compound_type: srt.chem_compound_type()};
            out.by_code_3.insert(srt.code3(), out.registered_types.len());
            out.registered_types.push(rt);
        }
        return out;
    }

    /// Provides a singleton instance of a [`ResidueTypeManager`] object
    ///
    /// # Example
    /// ```
    /// use bioshell_seq::chemical::ResidueTypeManager;
    /// let rt_mgr = ResidueTypeManager::get();
    /// assert!(rt_mgr.by_code3("ALA").is_some());
    /// ```
    pub fn get() -> MutexGuard<'static, ResidueTypeManager> {
        KNOWN_RESIDUE_TYPES.lock().unwrap()
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
    /// use bioshell_seq::chemical::{ResidueTypeManager};
    /// let mut mgr = ResidueTypeManager::get();
    /// // This should pass, as all standard residue types are preloaded by a constructor
    /// let ala = mgr.by_code3("ALA");
    /// assert!(ala.is_some());
    /// ```
    pub fn by_code3(&self, code3: &str) -> Option<&ResidueType> {
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

/// Residue type manager singleton provides information about registered residues on demand.
///
/// **Note**: This object may be globally accessed by calling [`ResidueTypeManager::get()`] method.
/// ```
/// use bioshell_seq::chemical::{KNOWN_RESIDUE_TYPES, ResidueTypeProperties};
/// let res_manager = KNOWN_RESIDUE_TYPES.lock().unwrap();
/// // ---------- try unknown residue, which should return `None`
/// let unknown_res = res_manager.by_code3("XYZ");
/// assert!(unknown_res.is_none());
/// // ---------- ALA should exist
/// let ala_res = res_manager.by_code3("ALA");
/// assert!(ala_res.is_some());
/// assert_eq!(ala_res.unwrap().parent_type.code1(),'A');
/// ```
pub static KNOWN_RESIDUE_TYPES: Lazy<Mutex<ResidueTypeManager>>
    = Lazy::new(|| Mutex::new(ResidueTypeManager::new()));

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
            pub const TYPES: [Self; 33] = [
            $(
                Self::$name,
            )*
            ];
        }
    };
}

impl StandardResidueType {
    /// Iterates over the 20 standard amino acid types.
    ///
    /// The iteration does not include the `UNK` residue type.
    ///
    /// # Examples
    /// ``` rust
    /// use bioshell_seq::chemical::StandardResidueType;
    /// assert_eq!(StandardResidueType::amino_acids().count(), 20);
    /// ```
    pub fn amino_acids() -> impl Iterator<Item = &'static StandardResidueType> {
        StandardResidueType::TYPES.iter()
            .filter(|srt| srt.chem_compound_type() == MonomerType::PeptideLinking && srt.code1() != 'X')
    }
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
    A 21 'a' "A" RNALinking,
    C 22 'c' "C" RNALinking,
    G 23 'g' "G" RNALinking,
    T 24 't' "T" RNALinking,
    U 25 'u' "U" RNALinking,
    DA 26 'a' "DA" DNALinking,
    DC 27 'c' "DC" DNALinking,
    DG 28 'g' "DG" DNALinking,
    DT 29 't' "DT" DNALinking,
    GAP 30 '-' "GAP" Other,
    GPE 31 '_' "GPE" Other,
    UNL 32 'Z' "UNL" Other
}


// ----------- the following map is used by TryFrom<&str> for MonomerType
static MONOMER_MAP: Lazy<HashMap<&'static str, MonomerType>> = Lazy::new(|| {
    let mut map = HashMap::new();
    map.insert("DBETAPEPTIDE", MonomerType::DBetaPeptide);
    map.insert("CGAMMALINKING", MonomerType::CGammaLinking);
    map.insert("DGAMMAPEPTIDE", MonomerType::DGammaPeptide);
    map.insert("CDELTALINKING", MonomerType::CDeltaLinking);
    map.insert("DPEPTIDECOOH", MonomerType::DPeptideCOOH);
    map.insert("DPEPTIDENH3", MonomerType::DPeptideNH3);
    map.insert("DPEPTIDELINKING", MonomerType::DPeptideLinking);
    map.insert("DSACCHARIDE", MonomerType::DSaccharide);
    map.insert("DSACCHARIDEALPHALINKING", MonomerType::DSaccharideAlphaLinking);
    map.insert("DSACCHARIDEBETALINKING", MonomerType::DSaccharideBetaLinking);
    map.insert("DNAOH3PRIMETERMINUS", MonomerType::DNAOH3PrimeTerminus);
    map.insert("DNAOH5PRIMETERMINUS", MonomerType::DNAOH5PrimeTerminus);
    map.insert("DNALINKING", MonomerType::DNALinking);
    map.insert("LDNALINKING", MonomerType::LDNALinking);
    map.insert("LRNALINKING", MonomerType::LRNALinking);
    map.insert("LBETAPEPTIDECGAMMALINKING", MonomerType::LBetaPeptideCGammaLinking);
    map.insert("LGAMMAPEPTIDECDELTALINKING", MonomerType::LGammaPeptideCDeltaLinking);
    map.insert("LPEPTIDECOOH", MonomerType::LPeptideCOOH);
    map.insert("LPEPTIDENH3", MonomerType::LPeptideNH3);
    map.insert("LPEPTIDELINKING", MonomerType::LPeptideLinking);
    map.insert("LSACCHARIDE", MonomerType::LSaccharide);
    map.insert("LSACCHARIDEALPHALINKING", MonomerType::LSaccharideAlphaLinking);
    map.insert("LSACCHARIDEBETALINKING", MonomerType::LSaccharideBetaLinking);
    map.insert("RNAOH3PRIMETERMINUS", MonomerType::RNAOH3PrimeTerminus);
    map.insert("RNAOH5PRIMETERMINUS", MonomerType::RNAOH5PrimeTerminus);
    map.insert("RNALINKING", MonomerType::RNALinking);
    map.insert("NONPOLYMER", MonomerType::NonPolymer);
    map.insert("OTHER", MonomerType::Other);
    map.insert("PEPTIDELINKING", MonomerType::PeptideLinking);
    map.insert("PEPTIDELIKE", MonomerType::PeptideLike);
    map.insert("SACCHARIDE", MonomerType::Saccharide);
    map
});