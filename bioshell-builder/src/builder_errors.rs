use thiserror::Error;

/// Errors that may appear while using the bioshell-builder crate
#[derive(Debug, Error)]
pub enum BuilderError {


    #[error("Error while creating an InternalAtomDefinition struct")]
    /// Can't parse `InternalAtomDefinition` data; make sure the line has proper columns
    InternalAtomDefinitionError {error: String},
    #[error("Can't find atom named: {atom_name}. Has it been defined in the residue {residue_index}?")]
    /// Can't find an atom defined for a given residue
    DefinedAtomNotFound {residue_index: usize, atom_name: String},
    #[error("Can't find a residue {residue_index}?")]
    /// Can't find a residue for a given index
    ResidueNotDefined {residue_index: usize},
    #[error("Can't find dihedral angle named: {dihedral_name}. Has it been defined in the residue {residue_index}?")]
    /// Can't find a named dihedral angle.
    ///
    /// Make sure its name is properly spelled and it has actually been defined in the respective InternalAtomDefinition entry
    DihedralAngleNotFound {residue_index: usize, dihedral_name: String},

}