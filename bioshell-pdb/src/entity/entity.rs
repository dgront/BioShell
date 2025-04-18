use std::collections::HashMap;
use std::str::FromStr;
use bioshell_cif::{CifData, CifTable};
use crate::{EntityType, PDBError, PolymerEntityType};
use crate::PDBError::{CantParseEnumVariant, InconsistentEntity, NoSuchChain, UnknownResidueType};
use bioshell_seq::chemical::{ResidueType, ResidueTypeManager};


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
///
/// # What is an entity?
///
/// An entity represents a chemically distinct object. For example, a hemoglobin deposit (PDB code: 4esa)
/// deposit contains six chemical entities:
/// ```
/// # use bioshell_pdb::PDBError;
/// # fn main() -> Result<(), PDBError> {
/// # use std::io::BufReader;
/// # use bioshell_pdb::Deposit;
/// # let pdb_data = include_str!("../../tests/test_files/4esa.cif");
/// let reader = BufReader::new(pdb_data.as_bytes());
/// let deposit_4esa = Deposit::from_cif_reader(reader)?;
/// assert_eq!(deposit_4esa.count_entities(), 6);
/// # Ok(())
/// # }
/// ```
/// Two of them are polymer entities: alpha and beta hemoglobin chains (entities named `"1"` and `"2"`), respectively.
/// The former one appears in chains `A` and `C` while latter - in chains `B` and `D`.
/// ```
/// # use bioshell_pdb::PDBError;
/// # fn main() -> Result<(), PDBError> {
/// # use std::io::BufReader;
/// # use bioshell_pdb::{Deposit, EntityType};
/// # use bioshell_pdb::PDBError::NoSuchEntity;
/// # use bioshell_pdb::PolymerEntityType::PolypeptideL;
/// # let pdb_data = include_str!("../../tests/test_files/4esa.cif");
/// # let reader = BufReader::new(pdb_data.as_bytes());
/// # let deposit_4esa = Deposit::from_cif_reader(reader)?;
/// let entity_1 = deposit_4esa.entity("1").ok_or_else(||NoSuchEntity{ entity_id: "1".to_string() })?;
/// assert_eq!(entity_1.entity_type(), EntityType::Polymer(PolypeptideL));
/// assert_eq!(entity_1.chain_ids(), &["A", "C"]);
///
/// let entity_2 = deposit_4esa.entity("2");
/// assert_eq!(entity_2.entity_type(), EntityType::Polymer(PolypeptideL));
/// assert_eq!(entity_2.chain_ids(), &["B", "D"]);
/// # Ok(())
/// # }
/// ```
/// There are also four non-polymer entities in the `4esa` deposit,
/// e.g. entity `"4"` comprises 408 water molecules while entity `"6"` four heme prostetic groups.
/// Each of these four HEM groups is placed in a different chain.
///
/// # Entities, chains and their sequences
///
/// The chains `"B"` and `"D"` of the `4esa` deposit contain a protein entity `"2"`. The entity sequence
/// comprises 146 amino acid residues:
/// ```txt
/// VEWTDQERATISSIFGSLDYDDIGPKALSRCLIVYPWTQRHFGSFGNLYNAEAIIGNQKVA
/// AHGIKVLHGLDRAVKNMDNIKEIYAELSILHSEKLHVDPDNFKLLADCLTIVVAAKMGSGF
/// NPGTQATFQKFLAVVVSALGKQYH
/// ```
/// The sequence can be accessed by the [`entity_monomers()`](Entity::entity_monomers()) method
/// and is provided as a vector of [ResidueType] objects.
/// The sequence observed for chain `"B"` may be accessed by calling [`chain_monomers()`](Entity::chain_monomers())
/// The sequence is shorter by two amino acid residues, since the C-terminal tyrosine (`Y`) and histidine (`H`)
/// hasn't been resolved experimentally. The missing residues are denoted in the chain sequence by the special [`GAP`](GAP)
/// residue type.
/// ```
/// # use bioshell_pdb::PDBError;
/// # fn main() -> Result<(), PDBError> {
/// # use std::io::BufReader;
/// # use bioshell_pdb::{Deposit, EntityType};
/// # use bioshell_pdb::PDBError::NoSuchEntity;
/// # use bioshell_pdb::PolymerEntityType::PolypeptideL;
/// # use bioshell_seq::chemical::StandardResidueType::GAP;
/// # let pdb_data = include_str!("../../tests/test_files/4esa.cif");
/// # let reader = BufReader::new(pdb_data.as_bytes());
/// # let deposit_4esa = Deposit::from_cif_reader(reader).ok_or(NoSuchEntity{ entity_id: "1".to_string() })?;;
/// let entity_2 = deposit_4esa.entity("2");
/// assert_eq!(entity_2.entity_monomers().len(), 146);
/// assert_eq!(entity_2.chain_monomers("B")?.iter().filter(|m| m.parent_type==GAP).count(), 2);
/// # Ok(())
/// # }
/// ```
///
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
    /// Monomers (or ligands) comprising each chain of this entity.
    chain_sequences: HashMap<String, Vec<ResidueType>>,
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
            chain_sequences: Default::default(),
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
    /// use bioshell_pdb::{PDBError, Deposit};
    /// # fn main() -> Result<(), PDBError> {
    /// use bioshell_pdb::Deposit;
    /// use bioshell_pdb::PDBError::NoSuchEntity;
    /// let cif_data = include_str!("../../tests/test_files/2fdo.cif");
    /// let deposit = Deposit::from_cif_reader(cif_data.as_bytes())?;
    /// let entity = deposit.entity("1").ok_or_else(|| NoSuchEntity{ entity_id: "1".to_string() })?;;
    /// assert_eq!(entity.chain_ids(), &vec!["B", "A"]);
    /// assert_eq!(entity.entity_monomers().len(), 94);
    /// # Ok(())
    /// # }
    /// ```
    ///
    pub fn chain_ids(&self) -> &Vec<String> { &self.chain_ids }

    /// Provides monomers (or ligands) comprising this entity.
    ///
    /// A single entity, such as a protein molecule, can be present in multiple chains.
    pub fn entity_monomers(&self) -> &Vec<ResidueType> { &self.monomer_sequence }

    /// Provides monomers (or ligands) comprising this entity as observed in a particular chain.
    ///
    /// The returned vector may contain [`GAP`](GAP) monomers when a residue present in the respective
    /// entity can't be observed in a given chain
    pub fn chain_monomers(&self, chain_id: &str) -> Result<&Vec<ResidueType>, PDBError> {

        if let Some(chain_seq) = self.chain_sequences.get(chain_id) {
            return Ok(chain_seq);
        } else {
            Err(NoSuchChain { chain_id: chain_id.to_string() })
        }
    }

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
            let mass = tokens[4].parse::<f64>().unwrap_or(0.0);
            entity_map.insert(tokens[0].to_string(),Entity::from_strings(tokens[0], tokens[1], tokens[2], tokens[3], mass)?);
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

        // --- load sequences for each chain as they are defined in entity definitions, including gaps
        let chains_in_entities = load_chain_residue_types(cif_data)?;
        for (entity_id, chains) in chains_in_entities {
            let entity = entity_map.get_mut(&entity_id).unwrap();
            entity.chain_sequences = chains;
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

/// Loads residue types for each polymer entity
///
/// This function loads polymer entities ONLY!
fn load_chain_residue_types(cif_data: &CifData) -> Result<HashMap<String,HashMap<String, Vec<ResidueType>>>, PDBError> {

    let res_mgr = ResidueTypeManager::get();
    let gap = res_mgr.by_code3("GAP").unwrap();
    let mut out = HashMap::default();
    if let Ok(seq_table) = CifTable::new(cif_data, "_pdbx_poly_seq_scheme",
                                             [".entity_id", ".pdb_mon_id", ".pdb_strand_id",]) {
        for [id, monomer, chain] in seq_table.iter() {
            if !out.contains_key(id) { out.insert(id.to_string(), HashMap::default()); }
            let map = out.get_mut(id).unwrap();
            if !map.contains_key(chain) { map.insert(chain.to_string(), vec![]); }
            let vec = map.get_mut(chain).unwrap();
            if monomer == "?" || monomer == "." {
                vec.push(gap.clone());
            } else {
                if let Some(res_type) = res_mgr.by_code3(monomer) { vec.push(res_type.clone()); }
                else {
                    return Err(UnknownResidueType { res_type: monomer.to_string() })
                }
            }
        }
    }
    return Ok(out);
}



struct PolymerEntity {
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

    pub fn from_cif_data(cif_data: &CifData) -> Result<Vec<PolymerEntity>, PDBError> {
        // ------------- HasMap to store data extracted from cif_data
        let mut entities: HashMap<String, (String, Vec<String>, Vec<String>)> = HashMap::new();
        let mut entity_order: Vec<String> = vec![];
        // ------------- Extract the list of all polymer entities
        if let Ok(entity_table)  = CifTable::new(cif_data, "_entity_poly.", ["entity_id", "type", "pdbx_strand_id"]) {
            for [entity_id, poly_type, chain_ids] in entity_table.iter() {
                // --- chain IDs are separated by commas like "A,B,C", so we need to split them
                // --- the first character is a semicolon, so we need to skip it
                let chains_fixed = chain_ids.trim_matches(';').trim();
                let chain_ids: Vec<String> =  chains_fixed.split(',').map(|s| s.trim().to_string()).collect();
                entities.insert(entity_id.to_string(), (poly_type.to_string(), chain_ids, Vec::new()));
                // --- record the order in which the entities were found
                entity_order.push(entity_id.to_string());
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

        // --- Sort entities based on the order in entity_order
        // Create a HashMap for quick index lookups
        let order_map: HashMap<_, _> = entity_order.iter().enumerate().map(|(i, id)| (id, i)).collect();
        // Sort entities based on the order in entity_order
        poly_entities.sort_by_key(|e| order_map.get(&e.entity_id).copied().unwrap_or(usize::MAX));

        Ok(poly_entities)
    }
}

#[cfg(test)]
mod test_poly_entity {
    use std::io::BufReader;
    use bioshell_cif::read_cif_buffer;
    use crate::{PDBError};
    use super::PolymerEntity;
    #[test]
    fn polymer_6qyd_entities() -> Result<(), PDBError> {

        let reader = BufReader::new(entities_6qyd.as_bytes());
        let cif_data = read_cif_buffer(reader).unwrap();
        let entities = PolymerEntity::from_cif_data(&cif_data[0]).unwrap();
        assert_eq!(entities.len(), 2);
        assert_eq!(entities[0].entity_id, "1");
        assert_eq!(entities[0].chain_ids.len(), 235);
        assert_eq!(entities[0].chain_ids[0], "1A");
        assert_eq!(entities[0].chain_ids[234], "9E");
        assert_eq!(entities[1].entity_id, "2");
        assert_eq!(entities[1].chain_ids.len(), 165);
        Ok(())
    }

    const entities_6qyd: &str = "data_6QYD
#
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.nstd_linkage
_entity_poly.nstd_monomer
_entity_poly.pdbx_seq_one_letter_code
_entity_poly.pdbx_seq_one_letter_code_can
_entity_poly.pdbx_strand_id
_entity_poly.pdbx_target_identifier
1 'polypeptide(L)' no no
;MRIT
;
;MRIT
;
;1A,1B,1C,1D,1E,1F,2A,2B,2C,2D,2E,3A,3B,3C,3D,3E,3F,4A,4B,4C,4D,4E,4F,5A,5B,5C,5D,5E,5F,6A,6B,6C,6D,6E,6F,7A,7B,7C,7D,7E,8A,8B,8C,8D,8E,8F,9A,1G,1H,1I,1J,1K,1L,2F,2G,2H,2I,2J,3G,3H,3I,3J,3K,3L,4G,4H,4I,4J,4K,4L,5G,5H,5I,5J,5K,5L,6G,6H,6I,6J,6K,6L,7F,7G,7H,7I,7J,8G,8H,8I,8J,8K,8L,9B,1M,1N,1O,1P,1Q,1R,2K,2L,2M,2N,2O,3M,3N,3O,3P,3Q,3R,4M,4N,4O,4P,4Q,4R,5M,5N,5O,5P,5Q,5R,6M,6N,6O,6P,6Q,6R,7K,7L,7M,7N,7O,8M,8N,8O,8P,8Q,8R,9C,1S,1T,1U,1V,1W,1X,2P,2Q,2R,2S,2T,3S,3T,3U,3V,3W,3X,4S,4T,4U,4V,4W,4X,5S,5T,5U,5V,5W,5X,6S,6T,6U,6V,6W,6X,7P,7Q,7R,7S,7T,8S,8T,8U,8V,8W,8X,9D,1Y,1Z,1a,1b,1c,1d,2U,2V,2W,2X,2Y,3Y,3Z,3a,3b,3c,3d,4Y,4Z,4a,4b,4c,4d,5Y,5Z,5a,5b,5c,5d,6Y,6Z,6a,6b,6c,6d,7U,7V,7W,7X,7Y,8Y,8Z,8a,8b,8c,8d,9E
;
?
2 'polypeptide(L)' no no
;MMVS
;
;MMVS
;
;1e,1f,1g,1h,2Z,2a,2b,2c,2d,3e,3f,3g,3h,4e,4f,5e,5f,6e,6f,6g,6h,7Z,7a,7b,7c,7d,8e,8f,8g,8h,8i,8j,9F,1i,1j,1k,1l,2e,2f,2g,2h,2i,3i,3j,3k,3l,4g,4h,5g,5h,6i,6j,6k,6l,7e,7f,7g,7h,7i,8k,8l,8m,8n,8o,8p,9G,1m,1n,1o,1p,2j,2k,2l,2m,2n,3m,3n,3o,3p,4i,4j,5i,5j,6m,6n,6o,6p,7j,7k,7l,7m,7n,8q,8r,8s,8t,8u,8v,9H,1q,1r,1s,1t,2o,2p,2q,2r,2s,3q,3r,3s,3t,4k,4l,5k,5l,6q,6r,6s,6t,7o,7p,7q,7r,7s,8w,8x,8y,8z,9K,9L,9I,1u,1v,1w,1x,2t,2u,2v,2w,2x,3u,3v,3w,3x,4m,4n,5m,5n,6u,6v,6w,6x,7t,7u,7v,7w,7x,9M,9N,9O,9P,9Q,9R,9J
;
?
#
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
1 1   MET n
1 2   ARG n
1 3   ILE n
1 4   THR n
2 1   MET n
2 2   MET n
2 3   VAL n
2 4   SER n
";
}