use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;
use once_cell::sync::Lazy;
use bioshell_cif::{CifData, CifTable};
use crate::{PDBError};
use crate::PDBError::{CantParseEnumVariant, InconsistentEntity, NoSuchChain, UnknownResidueType};
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
/// # let pdb_data = include_str!("../tests/test_files/4esa.cif");
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
/// # use bioshell_pdb::PolymerEntityType::PolypeptideL;
/// # let pdb_data = include_str!("../tests/test_files/4esa.cif");
/// # let reader = BufReader::new(pdb_data.as_bytes());
/// # let deposit_4esa = Deposit::from_cif_reader(reader)?;
/// let entity_1 = deposit_4esa.entity("1");
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
/// comprises 144 amino acid residues:
/// ```txt
/// VEWTDQERATISSIFGSLDYDDIGPKALSRCLIVYPWTQRHFGSFGNLYNAEAIIGNQKVA
/// AHGIKVLHGLDRAVKNMDNIKEIYAELSILHSEKLHVDPDNFKLLADCLTIVVAAKMGSGF
/// NPGTQATFQKFLAVVVSALGKQY
/// ```
/// The sequence can be accessed by the [`entity_monomers()`](Entity::entity_monomers()) method
/// and is provided as a vector of [ResidueType] objects.
/// The sequence observed for chain `"D"` may be accessed by calling [`chain_monomers()`](Entity::chain_monomers())
/// and is indeed identical with the respective entity definition. The sequence
/// of chain `"B"` however is one amino acid residue shorter since the N-terminal tyrosine (`Y`) hasn't been
/// resolved experimentally. The missing residue is denoted in the chain sequence by the special [`GAP`](GAP)
/// residue type.
/// ```
/// # use bioshell_pdb::PDBError;
/// # fn main() -> Result<(), PDBError> {
/// # use std::io::BufReader;
/// # use bioshell_pdb::{Deposit, EntityType};
/// # use bioshell_pdb::PolymerEntityType::PolypeptideL;
/// # use bioshell_seq::chemical::StandardResidueType::GAP;
/// # let pdb_data = include_str!("../tests/test_files/4esa.cif");
/// # let reader = BufReader::new(pdb_data.as_bytes());
/// # let deposit_4esa = Deposit::from_cif_reader(reader)?;
/// let entity_2 = deposit_4esa.entity("2");
/// assert_eq!(entity_2.entity_monomers().len(), 146);
/// assert_eq!(entity_2.chain_monomers("B")?.iter().filter(|m| m.parent_type==GAP).count(), 1);
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
    /// let cif_data = include_str!("../tests/test_files/2fdo.cif");
    /// let deposit = Deposit::from_cif_reader(cif_data.as_bytes())?;
    /// let entity = deposit.entity("1");
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
            let mut entity = entity_map.get_mut(&entity_id).unwrap();
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

fn load_chain_residue_types(cif_data: &CifData) -> Result<HashMap<String,HashMap<String, Vec<ResidueType>>>, PDBError> {

    let res_mgr = ResidueTypeManager::get();
    let gap = res_mgr.by_code3("GAP").unwrap();
    let mut out = HashMap::default();
    if let Ok(seq_table) = CifTable::new(cif_data, "_pdbx_poly_seq_scheme",
                                             [".entity_id", ".pdb_mon_id", ".pdb_strand_id",]) {
        for [id, monomer, chain] in seq_table.iter() {
            if !out.contains_key(id) { out.insert(id.to_string(), HashMap::default()); }
            let mut map = out.get_mut(id).unwrap();
            if !map.contains_key(chain) { map.insert(chain.to_string(), vec![]); }
            let mut vec = map.get_mut(chain).unwrap();
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