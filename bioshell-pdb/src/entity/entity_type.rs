use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;
use once_cell::sync::Lazy;
use crate::PDBError;
use crate::PDBError::CantParseEnumVariant;

/// Represents the different types of entities in CIF data.
///
/// For more information, see the [definition of `_entity.type`](https://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Items/_entity.type.html).
///
/// # Example
/// ```
/// use bioshell_pdb::EntityType;
/// let enity_type: EntityType = "non-polymer".parse().unwrap();
/// assert_eq!(enity_type, EntityType::NonPolymer);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EntityType {
    /// The entity is a polymer.
    Polymer(PolymerEntityType),
    /// The entity is a non-polymer (e.g., a small molecule or ligand).
    NonPolymer,
    /// The entity is water.
    Water,
    /// The entity is branched, typically a polysaccharide.
    Branched,
}

impl fmt::Display for EntityType {

    /// Returns a string representation of the [`EntityType`].
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let description = match self {
            EntityType::Polymer(_) => "Polymer",
            EntityType::NonPolymer => "Non-polymer",
            EntityType::Water => "Water",
            EntityType::Branched => "Branched",
        };
        write!(f, "{}", description)
    }
}

impl FromStr for EntityType {
    type Err = PDBError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "polymer" => Ok(EntityType::Polymer(PolymerEntityType::Other)),
            "non-polymer" => Ok(EntityType::NonPolymer),
            "water" => Ok(EntityType::Water),
            "branched" => Ok(EntityType::Branched),
            _ => Err(CantParseEnumVariant{ data_value: s.to_string(), enum_name: "EntityType".to_string() }),
        }
    }
}


/// Enum representing different types of polymers for _entity_poly.type in an mmCIF file.
///
/// # Example
/// ```
/// use bioshell_pdb::PolymerEntityType;
/// let src_method: PolymerEntityType = "Polypeptide(L)".parse().unwrap();
/// assert_eq!(src_method, PolymerEntityType::PolypeptideL);
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash, Copy)]
pub enum PolymerEntityType {
    /// L polypeptide
    PolypeptideL,
    /// D polypeptide
    PolypeptideD,
    /// DNA polymer.
    DNA,
    /// RNA polymer.
    RNA,
    /// L polysaccharide polymer.
    PolysaccharideL,
    /// D polysaccharide polymer.
    PolysaccharideD,
    /// Represents a peptide nucleic acid.
    PeptideNucleicAcid,
    /// Represents a hybrid polymer type.
    Other,
}

impl FromStr for PolymerEntityType {
    type Err = PDBError;
    /// Parse a string and return the corresponding [`PolymerEntityType`] enum variant.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let tmp = s.trim_matches('\'').to_lowercase();
        // --- First try to find an exact match
        if let Some(&entity_type) = ENTITY_TYPE_MAP.get(tmp.as_str()) {
            return Ok(entity_type);
        }
        // --- If no exact match is found, try to find a substring match
        for (key, &entity_type) in ENTITY_TYPE_MAP.iter() {
            if s.contains(key) { return Ok(entity_type); }
        }
        // --- If neither exact nor substring match is found, return an error
        Err(CantParseEnumVariant {
            data_value: s.to_string(),
            enum_name: "PolymerEntityType".to_string(),
        })
    }
}

/// Returns a string representation of the [`PolymerEntityType`].
impl fmt::Display for PolymerEntityType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let description = match self {
            PolymerEntityType::PolypeptideL => "L polypeptide",
            PolymerEntityType::PolypeptideD => "D polypeptide",
            PolymerEntityType::DNA => "DNA polymer",
            PolymerEntityType::RNA => "RNA polymer",
            PolymerEntityType::PolysaccharideL => "L polysaccharide polymer",
            PolymerEntityType::PolysaccharideD => "D polysaccharide polymer",
            PolymerEntityType::PeptideNucleicAcid => "Peptide nucleic acid",
            PolymerEntityType::Other => "Other hybrid polymer type",
        };
        write!(f, "{}", description)
    }
}


// --- entity type map used to implement FromStr
static ENTITY_TYPE_MAP: Lazy<HashMap<&'static str, PolymerEntityType>> = Lazy::new(|| {
    let mut m = HashMap::new();
    m.insert("polypeptide(l)", PolymerEntityType::PolypeptideL);
    m.insert("polypeptide(d)", PolymerEntityType::PolypeptideD);
    m.insert("dna", PolymerEntityType::DNA);
    m.insert("rna", PolymerEntityType::RNA);
    m.insert("polysaccharide(l)", PolymerEntityType::PolysaccharideL);
    m.insert("polysaccharide(d)", PolymerEntityType::PolysaccharideD);
    m.insert("polydeoxyribonucleotide", PolymerEntityType::DNA);
    m.insert("polyribonucleotide", PolymerEntityType::RNA);
    m.insert("peptide nucleic acid", PolymerEntityType::PeptideNucleicAcid);
    m.insert("other", PolymerEntityType::Other);
    m
});