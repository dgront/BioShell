use std::str::SplitWhitespace;

#[derive(Debug, Ord, PartialOrd, Eq, PartialEq, Clone)]
#[repr(i8)]
pub enum MonomerType {
    /// standard amino acid, single-character code: 'P'
    PeptideLinking = 0,
    /// standard deoxy-nucleotide, single-character code: 'D'
    DNALinking,
    /// standard nucleotide, single-character code: 'R'
    RNALinking,
    /// a sugar residue, single-character code: 'S'
    Sacharide,
    /// most likely a ligand, single-character code: 'N'
    NonPolymer,
    /// anything else, single-character code: 'O'
    OTHER
}

impl TryFrom<char> for MonomerType {
    type Error = &'static str;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            'P' => Ok(MonomerType::PeptideLinking),
            'D' => Ok(MonomerType::DNALinking),
            'R' => Ok(MonomerType::RNALinking),
            'S' => Ok(MonomerType::Sacharide),
            'N' => Ok(MonomerType::NonPolymer),
            'O' => Ok(MonomerType::OTHER),
            _ => {Err(&*format!("Can't find a MonomerType for the one-letter code: {}!", value))}
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
pub struct ResidueType {
    /// three-letter code such as ALA for alanine
    pub code3: String,
    /// "standard" residue that this residue type is a variant of
    pub parent_type: StandardResidueType,
    /// chemical type of this residue type
    pub chem_compound_type: MonomerType,
}


pub struct ResidueTypeManager {
    registered_types: Vec<ResidueType>
}

impl TryFrom<String> for ResidueType {
    type Error = &'static str;

    /// Parses a string into a ResidueType instance.
    /// Expected string format: ``"ALN A P"``
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

impl ResidueTypeManager {
    pub fn new() ->  ResidueTypeManager{
        let mut out = ResidueTypeManager{registered_types:vec![]};

        for rt in StandardResidueType::TYPES {
            out.registered_types.push(ResidueType{code3:rt.code3(), parent_type:rt,
                chem_compound_type:rt.chem_compound_type()});
        }
        return out;
    }

    pub fn count(&self) -> usize { self.registered_types.len() }

    pub fn register_residue_type(&mut self, res_type: ResidueType) {
        self.registered_types.push(res_type);
    }
}


macro_rules! define_res_types {

    ($($name:ident $id:literal $one_letter_code:literal $three_letters_code:literal $res_type:ident),*) => {
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
            type Error = &'static str;

            /// Returns a [`StandardResidueType`](StandardResidueType) enum for a given one-letter code
            fn try_from(value: char) -> Result<Self, Self::Error> {
                match value {
                    $(
                        $one_letter_code => Ok(StandardResidueType::$name),
                    )*
                    _ => {Err(&*format!("Can't find an amino acid for the one-letter code: {}!", value))}
                }
            }
        }

        impl StandardResidueType {
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
