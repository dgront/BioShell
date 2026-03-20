use thiserror::Error;
use bioshell_cif::CifError;
use crate::ResidueId;

/// Errors that may appear while using the bioshell-pdb crate
#[derive(Debug, Error)]
pub enum PDBError {
    #[error("Invalid file format: {file_name}")]
    /// A given file is neither a PDB nor an mmCIF file
    InvalidFileFormat {file_name: String},

    #[error("Invalid PDB file format: {broken_pdb_line}")]
    /// Invalid format of a PDB line
    InvalidPdbLineFormat {broken_pdb_line: String},

    #[error("General I/O error occurred while reading a PDB or mmCIF file")]
    /// I/O error occurred while reading a PDB or a CIF file
    Io(#[from] std::io::Error),

    #[error("Residue not found: {res_id}")]
    /// Residue corresponding to a given `res_id` could not be located
    NoSuchResidue {res_id: ResidueId},

    #[error("Chain not found: {chain_id}")]
    /// Chain corresponding to a given `chain_id` could not be located
    NoSuchChain {chain_id: String},

    #[error("Entity not found: {entity_id}")]
    /// Can't find entity for the given `entity_id` string
    NoSuchEntity {entity_id: String},

    #[error("Residue is a terminal one: {res_id}, proceeding or following residue is a terminal can't be located")]
    /// Residue following or proceeding the given `res_id` could not be located
    ResidueAtTerminus {res_id: ResidueId},

    #[error("Atom not found: {atom_name} in the residue {res_id}")]
    /// Atom named `atom_name` could not be located in a residue `res_id`
    NoSuchAtom {atom_name: String, res_id: ResidueId},

    #[error("Residue type not registered: {res_type}")]
    /// Unknown 3-letter residue code: `res_type`
    UnknownResidueType {res_type: String},

    #[error("Can't find a file with input parameters: {fname} - check your BIOSHEL_DB_PATH system variable")]
    /// Missing parameters' file, that should be found in a BioShell's database
    MissingBioShellFile {fname: String},

    #[error("Error while creating an InternalAtomDefinition struct")]
    /// Can't create `InternalAtomDefinition`
    InternalAtomDefinitionError {error: String},

    #[error("Can't find atom named: {atom_name}. Has it been defined in the residue {residue_index}?")]
    /// Can't find an atom defined for a given residue
    DefinedAtomNotFound {residue_index: usize, atom_name: String},
    
    #[error("Can't find a residue {residue_index}")]
    /// Can't find a residue for a given index
    ResidueNotDefined {residue_index: usize},

    #[error("The number of atoms in model {model_index} is different than the number of atoms in the first model")]
    /// A requested model has different number of atoms that the first model
    WrongAtomsNumberInModel {model_index: usize},

    #[error("Unknown chemical compound type {compound_type} used to define monomer {compound_id}")]
    /// Can't parse a chemical component type
    IncorrectCompoundTypeName {compound_id: String, compound_type: String},

    #[error("Cif parser returned an error")]
    /// A Cif parser returned an error which is stored inside this enum variant
    CifParsingError(#[from] CifError),

    #[error("The following string {data_value} can't be parsed into an enum {enum_name} variant")]
    /// A given string can't be parsed into an enum variant
    CantParseEnumVariant {data_value: String, enum_name: String},

    #[error("The following entity {entity_id} is inconsistent with other data entries: {details}")]
    /// A given entity is inconsistent with other data entries
    InconsistentEntity {entity_id: String, details: String},

    #[error("Error while parsing a residue id string: {residue_id}")]
    /// Error while parsing a residue id string
    ResidueIdParsingError {residue_id: String},

    #[error("No PDB / CIF data loaded")]
    /// No PDB / CIF data loaded
    NoStructureDataLoaded,

    #[error("Can't download the {pdb_id} deposit for from RCSB: {reason}")]
    /// Can't download mmCIF file from RCSB website a deposit
    CantDownladFromRCSB{pdb_id: String, reason: String},

    #[error("structure cannot be converted to PDB-compatible form: {reason}")]
    PdbConversionNotPossible {
        reason: PdbConversionImpossibleReason,
    },
}

/// Provides detailed information why a [`Structure`] can't be converted to
#[derive(Debug, Error, Clone)]
pub enum PdbConversionImpossibleReason {
    #[error("too many chains for PDB: {chains} (max {max}); chain IDs must be 1 character")]
    TooManyChains { chains: usize, max: usize },

    #[error("too many atoms for PDB atom serial field: {atoms} (max {max})")]
    TooManyAtoms { atoms: usize, max: usize },

    #[error("atom name too long for PDB (max {max} chars): '{atom_name}'")]
    AtomNameTooLong { atom_name: String, max: usize },

    #[error("residue name too long for PDB (max {max} chars): '{res_name}'")]
    ResidueNameTooLong { res_name: String, max: usize },

    #[error("element symbol too long for PDB (max {max} chars): '{element}'")]
    ElementTooLong { element: String, max: usize },

    #[error("charge field too long for PDB (max {max} chars): '{charge}'")]
    ChargeTooLong { charge: String, max: usize },

    #[error("residue numbering overflow after renumbering; exceeded PDB limit {max} within a chain")]
    ResidueNumberOverflow { max: i32 },
}