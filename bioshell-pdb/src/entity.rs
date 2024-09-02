use std::str::FromStr;
use bioshell_cif::{cif_columns_by_name, CifData, CifLoop, CifError};
use crate::{PDBError, value_or_missing_key_pdb_error};
use crate::PDBError::{CantParseEnumVariant, CifParsingError};
use bioshell_cif::CifError::{ItemParsingError, MissingCifDataKey};

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
    Polymer,
    /// The entity is a non-polymer (e.g., a small molecule or ligand).
    NonPolymer,
    /// The entity is water.
    Water,
}


impl FromStr for EntityType {
    type Err = PDBError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "polymer" => Ok(EntityType::Polymer),
            "non-polymer" => Ok(EntityType::NonPolymer),
            "water" => Ok(EntityType::Water),
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
    /// use bioshell_pdb::{Entity, EntitySource, EntityType};
    /// let entity = Entity::from_strings("1", "HIV protease", "polymer", "man", 10916.0);
    /// assert_eq!(entity.id(), "1");
    /// assert_eq!(entity.description(), "HIV protease");
    /// assert_eq!(entity.entity_type(), EntityType::Polymer);
    /// assert_eq!(entity.src_method(), EntitySource::GeneticallyManipulated);
    /// assert_eq!(entity.formula_weight(), 10916.0);
    /// ```
    pub fn from_strings(id: &str, description: &str, entity_type: &str,
                        src_method: &str, formula_weight: f64, ) -> Result<Entity, PDBError> {
        Ok(Entity {
            id: id.to_string(), description: description.to_string(),
            entity_type:EntityType::from_str(entity_type)?,
            src_method: EntitySource::from_str(src_method)?, formula_weight,
        })
    }

    /// Returns the ID of the entity.
    pub fn id(&self) -> &str { &self.id }
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
    /// let data_block = &read_cif_buffer(&mut BufReader::new(input_cif.as_bytes()))[0];
    /// let entities = Entity::from_cif_data(data_block).unwrap();
    /// ```
    pub fn from_cif_data(cif_data: &CifData) -> Result<Vec<Entity>, PDBError> {
        let mut entities = Vec::new();

        if let Some(entity_loop) = cif_data.first_loop("_entity.id") {
            cif_columns_by_name!(EntityData, "_entity.id",
                "_entity.pdbx_description", "_entity.type",
                "_entity.src_method", "_entity.formula_weight",
            );
            let extractor = EntityData::new(entity_loop)?;
            let mut tokens = [""; 5];
            for row in entity_loop.rows() {
                extractor.data_items(&row, &mut tokens);
                if let Ok(mass) = tokens[4].parse::<f64>() {
                    entities.push(Entity::from_strings(tokens[0], tokens[1], tokens[2], tokens[3], mass)?);
                } else {
                    return Err(CifParsingError(ItemParsingError{
                        item: tokens[4].to_string(),
                        type_name: "f64".to_string(),
                        details: "_entity.formula_weight can't be parsed".to_string(),
                    }));
                }
            }
            Ok(entities)
        } else {
            let id = value_or_missing_key_pdb_error!(cif_data, "_entity.id", String);
            let desc = value_or_missing_key_pdb_error!(cif_data, "_entity.pdbx_description", String);
            let entity_type = value_or_missing_key_pdb_error!(cif_data, "_entity.type", String);
            let src_method = value_or_missing_key_pdb_error!(cif_data, "_entity.src_method", String);
            let weight = value_or_missing_key_pdb_error!(cif_data, "_entity.formula_weight", f64);
            entities.push(Entity::from_strings(&id, &desc, &entity_type, &src_method, weight)?);
            Ok(entities)
        }
    }

}