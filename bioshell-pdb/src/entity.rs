use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;
use bioshell_cif::{CifData, CifTable};
use crate::{PDBError};
use crate::PDBError::{CantParseEnumVariant, CifParsingError, InconsistentEntity};
use bioshell_cif::CifError::{ItemParsingError};
use bioshell_seq::chemical::{ResidueType, ResidueTypeManager};

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


/// Represents the different sources an entity in CIF data may be obtained from.
///
/// For more information, see the [definition of `_entity.src_method`](https://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Items/_entity.src_method.html).
///
/// # Example
/// ```
/// use bioshell_pdb::EntitySource;
/// let src_method: EntitySource = "nat".parse().unwrap();
/// assert_eq!(src_method, EntitySource::Natural);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EntitySource {
    /// entity isolated from a natural source
    Natural,
    /// entity isolated from a genetically manipulated source.
    GeneticallyManipulated,
    /// entity obtained synthetically
    Synthesized,
}


impl FromStr for EntitySource {
    type Err = PDBError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "nat" => Ok(EntitySource::Natural),
            "man" => Ok(EntitySource::GeneticallyManipulated),
            "syn" => Ok(EntitySource::Synthesized),
            _ => Err(CantParseEnumVariant{ data_value: s.to_string(), enum_name: "EntitySource".to_string() }),
        }
    }
}


/// Represents an entity in an mmCIF file.
#[derive(Debug, Clone)]
pub struct Entity {
    /// Unique identifier for the entity.
    id: String,
    /// Description of the entity.
    description: String,
    /// Type of the entity.
    entity_type: EntityType,
    /// Source method of the entity.
    src_method: EntitySource,
    /// Molecular mass of the entity in Daltons.
    formula_weight: f64,
    /// Which chains contain that entity?
    chain_ids: Vec<String>,
    /// Monomers (or ligands) comprising this entity.
    monomer_sequence: Vec<ResidueType>,
}


impl Entity {
    /// Constructs a new `Entity` instance.
    ///
    /// # Arguments
    ///
    /// * `id` - A string slice that holds the ID of the entity.
    /// * `description` - A string slice that holds a description of the entity.
    /// * `entity_type` - The type of the entity, represented by the `EntityType` enum.
    /// * `src_method` - The source method of the entity,
    /// * `formula_weight` - An optional string slice that holds the sequence information of the entity.
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{Entity, EntitySource, EntityType, PolymerEntityType};
    /// let entity = Entity::from_strings("1", "HIV protease", "polymer", "man", 10916.0).unwrap();
    /// assert_eq!(entity.id(), "1");
    /// assert_eq!(entity.description(), "HIV protease");
    /// assert_eq!(entity.entity_type(), EntityType::Polymer(PolymerEntityType::Other));
    /// assert_eq!(entity.src_method(), EntitySource::GeneticallyManipulated);
    /// assert_eq!(entity.formula_weight(), 10916.0);
    /// ```
    pub fn from_strings(id: &str, description: &str, entity_type: &str,
                        src_method: &str, formula_weight: f64, ) -> Result<Entity, PDBError> {
        Ok(Entity {
            id: id.to_string(), description: description.to_string(),
            entity_type:EntityType::from_str(entity_type)?,
            src_method: EntitySource::from_str(src_method)?, formula_weight,
            chain_ids: vec![],
            monomer_sequence: vec![],
        })
    }

    /// Returns the ID of the entity.
    pub fn id(&self) -> &str { &self.id }

    /// Provides identifiers of the chains that contain this entity.
    ///
    /// A single entity, such as a protein molecule, can be present in multiple chains. The sequence
    /// of this entity should be the same in all chains, in practice however it often differs due to
    /// experimental conditions.
    ///
    /// # Example
    /// ```
    /// use bioshell_pdb::{PDBError, load_cif_reader};
    /// # fn main() -> Result<(), PDBError> {
    /// let cif_data = include_str!("../tests/test_files/2fdo.cif");
    /// let strctr = load_cif_reader(cif_data.as_bytes())?;
    /// let entity = strctr.entity("1");
    /// assert_eq!(entity.chain_ids(), &vec!["B", "A"]);
    /// assert_eq!(entity.monomer_sequence().len(), 94);
    /// # Ok(())
    /// # }
    /// ```
    ///
    pub fn chain_ids(&self) -> &Vec<String> { &self.chain_ids }

    /// Provides monomers (or ligands) comprising this entity.
    ///
    /// A single entity, such as a protein molecule, can be present in multiple chains.
    pub fn monomer_sequence(&self) -> &Vec<ResidueType> { &self.monomer_sequence }

    /// Returns the description of the entity.
    pub fn description(&self) -> &str { &self.description }

    /// Returns the type of the entity.
    ///
    /// An entity can be a [`Polymer`](EntityType::Polymer), a [`NonPolymer`](EntityType::NonPolymer), or [`Water`](EntityType::Water).
    pub fn entity_type(&self) -> EntityType { self.entity_type }

    /// Returns the source method of the entity.
    ///
    /// An entity could be obtained as [`Natural`](EntitySource::Natural),  [`GeneticallyManipulated`](EntitySource::GeneticallyManipulated), or [`Synthesized`](EntitySource::Synthesized).
    pub fn src_method(&self) -> EntitySource { self.src_method }

    /// Returns the molecular mass of the entity in Daltons.
    pub fn formula_weight(&self) -> f64 { self.formula_weight }

    /// Creates a list of `Entity` instances from the given `CifData`.
    ///
    /// # Example
    /// ```
    /// use std::io::BufReader;
    /// use bioshell_cif::read_cif_buffer;
    /// use bioshell_pdb::Entity;
    /// let input_cif = r#"data_1O56
    /// loop_
    /// _entity.id
    /// _entity.type
    /// _entity.src_method
    /// _entity.pdbx_description
    /// _entity.formula_weight
    /// 1 polymer syn "DNA (5'-CD(*AP*AP*AP*)-3')" 894.663
    /// 2 water   nat water                        18.015
    /// "#;
    /// let data_block = &read_cif_buffer(&mut BufReader::new(input_cif.as_bytes())).unwrap()[0];
    /// let entities = Entity::from_cif_data(data_block).unwrap();
    /// ```
    pub fn from_cif_data(cif_data: &CifData) -> Result<HashMap<String, Entity>, PDBError> {

        let mut entity_map = HashMap::default();
        // ---------- extract entity table from cif data into a hashmap
        let entity_table = CifTable::new(cif_data, "_entity.",
            ["id", "pdbx_description", "type", "src_method", "formula_weight",])?;
        for tokens in entity_table.iter() {
            if let Ok(mass) = tokens[4].parse::<f64>() {
                entity_map.insert(tokens[0].to_string(),Entity::from_strings(tokens[0], tokens[1], tokens[2], tokens[3], mass)?);
            } else {
                return Err(CifParsingError(ItemParsingError{
                    item: tokens[4].to_string(),
                    type_name: "f64".to_string(),
                    details: "_entity.formula_weight can't be parsed".to_string(),
                }));
            }
        }

        // ---------- extract also the vector of polymer entities
        let polymers = PolymerEntity::from_cif_data(cif_data)?;

        // ---------- copy data from polymer entities to the corresponding entity
        for polymer in &polymers {
            let id = &polymer.entity_id;
            if let Some(entity) = entity_map.get_mut(id) {
                entity.chain_ids = polymer.chain_ids.clone();
                entity.monomer_sequence = polymer.monomer_sequence.clone();
                entity.entity_type = EntityType::Polymer(polymer.poly_type);
            } else {
                return Err(InconsistentEntity {
                    entity_id: id.clone(),
                    details: "PolymerEntity found, but Entity missing".to_string(),
                });
            }
        }

        // ---------- Now extract nonpolymer entities, which are not mandatory
        if let Ok(nonpoly_table) = CifTable::new(cif_data, "_pdbx_nonpoly_scheme",
                [".entity_id", ".mon_id", ".pdb_strand_id",]) {
            // ---------- to find the residue type, we need the residue type manager
            let mgr = ResidueTypeManager::get();
            for [entity_id, monomer_id, chain_id] in nonpoly_table.iter() {
                let chain_id_str = chain_id.to_string();
                if let Some(entity) = entity_map.get_mut(entity_id) {
                    // --- add a chain for this entity
                    if !entity.chain_ids.contains(&chain_id_str) {
                        entity.chain_ids.push(chain_id_str);
                    }
                    // --- residue type for that non-polymer molecule
                    let m = mgr.by_code3(monomer_id).ok_or(PDBError::UnknownResidueType{res_type: monomer_id.to_string()})?;
                    entity.monomer_sequence.push(m.clone());
                    // --- if it's water, set the non-polymer type accordingly
                    if entity.entity_type == EntityType::NonPolymer || monomer_id == "HOH" {
                        entity.entity_type = EntityType::Water;
                    }
                } else {
                    return Err(InconsistentEntity {
                        entity_id: entity_id.to_string(),
                        details: "NonPolymerEntity found, but Entity missing".to_string(),
                    });
                }
            }
        }

        Ok(entity_map)
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
    /// Parse a string and return the corresponding `EntityPolyType` enum variant.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let tmp = s.trim_matches('\'').to_lowercase();
        match tmp.as_str() {
            "polypeptide(l)" => Ok(PolymerEntityType::PolypeptideL),
            "polypeptide(d)" => Ok(PolymerEntityType::PolypeptideD),
            "dna" => Ok(PolymerEntityType::DNA),
            "rna" => Ok(PolymerEntityType::RNA),
            "polysaccharide(l)" => Ok(PolymerEntityType::PolysaccharideL),
            "polysaccharide(d)" => Ok(PolymerEntityType::PolysaccharideD),
            "polydeoxyribonucleotide" => Ok(PolymerEntityType::DNA),
            "polyribonucleotide" => Ok(PolymerEntityType::RNA),
            "other" => Ok(PolymerEntityType::Other),
            _ => Err(CantParseEnumVariant{ data_value: s.to_string(), enum_name: "EntityPolyType".to_string() }),
        }
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

pub struct PolymerEntity {
    entity_id: String,
    poly_type: PolymerEntityType,
    chain_ids: Vec<String>,
    monomer_sequence: Vec<ResidueType>,
}

impl PolymerEntity {
    pub fn from_str(entity_id: &str, poly_type: &str, chain_ids: Vec<String>, monomer_sequence: Vec<String>) -> Result<Self, PDBError> {
        let poly_type = PolymerEntityType::from_str(poly_type)?;
        let mgr = ResidueTypeManager::get();
        let mut monomers = vec![];
        for code in monomer_sequence {
            let m = mgr.by_code3(&code).ok_or(PDBError::UnknownResidueType{res_type: code.clone()})?;
            monomers.push(m.clone());
        }
        Ok(PolymerEntity {
            entity_id: entity_id.to_string(),
            poly_type,
            chain_ids,
            monomer_sequence: monomers,
        })
    }

    /// Load a list of polymer entities from a mmCIF file.
    ///
    /// ```
    /// use std::io::BufReader;
    /// use bioshell_cif::read_cif_buffer;
    /// use bioshell_pdb::{PolymerEntity,PDBError};
    /// # fn main() -> Result<(), PDBError> {
    /// let cif_data = include_str!("../tests/test_files/2gb1.cif");
    /// let reader = BufReader::new(cif_data.as_bytes());
    /// let cif_data = read_cif_buffer(reader).unwrap();
    /// let entities = PolymerEntity::from_cif_data(&cif_data[0])?;
    /// // --- 2GB1 protein has only one chain of 56 residues
    /// assert_eq!(entities.len(), 1);
    /// assert_eq!(entities[0].monomer_sequence().len(), 56);
    /// # Ok(())
    /// # }
    /// ```
    pub fn from_cif_data(cif_data: &CifData) -> Result<Vec<PolymerEntity>, PDBError> {
        // ------------- HasMap to store data extracted from cif_data
        let mut entities: HashMap<String, (String, Vec<String>, Vec<String>)> = HashMap::new();
        // ------------- Extract the list of all polymer entities
        if let Ok(entity_table)  = CifTable::new(cif_data, "_entity_poly.", ["entity_id", "type", "pdbx_strand_id"]) {
            for [entity_id, poly_type, chain_ids] in entity_table.iter() {
                // --- chain IDs are separated by commas like "A,B,C", so we need to split them
                let chain_ids: Vec<String> = chain_ids.split(',').map(|s| s.to_string()).collect();
                entities.insert(entity_id.to_string(), (poly_type.to_string(), chain_ids, Vec::new()));
            }
        } else {
            return Ok(vec![]);        // --- no polymer entities found!
        }
        // ------------- Extract the sequence for each of the polymer entities
        let entity_table = CifTable::new(cif_data, "_entity_poly_seq.",["entity_id", "mon_id"])?;
        for [entity_id, mon_id] in entity_table.iter() {
            let entity_id = entity_id.to_string();
            let mon_id = mon_id.to_string();
            if let Some((_, _, seq)) = entities.get_mut(&entity_id) {
                seq.push(mon_id);
            }
        }
        // ------------- Convert the data into PolymerEntity objects and return them
        let mut poly_entities: Vec<PolymerEntity> = vec![];
        for (entity_id, (poly_type, chain_ids, monomer_sequence)) in entities {
            poly_entities.push(PolymerEntity::from_str(&entity_id, &poly_type, chain_ids, monomer_sequence).unwrap());
        }

        Ok(poly_entities)
    }

    pub fn entity_id(&self) -> &str { &self.entity_id }

    pub fn monomer_sequence(&self) -> &Vec<ResidueType> { &self.monomer_sequence }

    pub fn entity_type(&self) -> &PolymerEntityType { &self.poly_type }

    pub fn chain_ids(&self) -> &Vec<String> { &self.chain_ids }
}