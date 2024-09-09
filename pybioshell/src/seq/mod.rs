use pyo3::prelude::*;
use pyo3::types::PyString;

// Defines the types of monomers - residue types that are biomolecular building blocks.
#[derive(Debug, Ord, PartialOrd, Eq, PartialEq, Clone, Copy)]
#[repr(i8)]
#[pyclass]
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

#[pymethods]
impl MonomerType {
    #[new]
    fn new(value: i8) -> Self {
        match value {
            0 => MonomerType::PeptideLinking,
            1 => MonomerType::LPeptideLinking,
            2 => MonomerType::LPeptideCOOH,
            3 => MonomerType::LPeptideNH3,
            4 => MonomerType::DBetaPeptide,
            5 => MonomerType::CGammaLinking,
            6 => MonomerType::DGammaPeptide,
            7 => MonomerType::CDeltaLinking,
            8 => MonomerType::DPeptideCOOH,
            9 => MonomerType::DPeptideNH3,
            10 => MonomerType::DPeptideLinking,
            11 => MonomerType::LBetaPeptideCGammaLinking,
            12 => MonomerType::LGammaPeptideCDeltaLinking,
            13 => MonomerType::PeptideLike,
            14 => MonomerType::RNALinking,
            15 => MonomerType::LRNALinking,
            16 => MonomerType::RNAOH3PrimeTerminus,
            17 => MonomerType::RNAOH5PrimeTerminus,
            18 => MonomerType::DNALinking,
            19 => MonomerType::LDNALinking,
            20 => MonomerType::DNAOH3PrimeTerminus,
            21 => MonomerType::DNAOH5PrimeTerminus,
            22 => MonomerType::Saccharide,
            23 => MonomerType::LSaccharide,
            24 => MonomerType::LSaccharideAlphaLinking,
            25 => MonomerType::LSaccharideBetaLinking,
            26 => MonomerType::DSaccharide,
            27 => MonomerType::DSaccharideAlphaLinking,
            28 => MonomerType::DSaccharideBetaLinking,
            29 => MonomerType::NonPolymer,
            30 => MonomerType::Other,
            _ => panic!("Invalid MonomerType value"),
        }
    }
}

/// Returns the integer value of the MonomerType
#[pyfunction]
pub fn get_monomer_type_value(monomer_type: MonomerType) -> i8 {
    monomer_type as i8
}


