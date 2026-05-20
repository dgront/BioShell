use std::collections::HashMap;
use std::io::BufRead;

use bioshell_cif::{read_cif_buffer};
use bioshell_cif::CifTable;
use crate::parse;
use crate::{Atom, BondType, ChemErrors, Element, Molecule};

/// Reads a molecule topology from a  `.cif` file.
///
/// Example file for ethanol molecule can be found on (the RCSB website)[https://files.rcsb.org/ligands/view/EOH.cif]
///
/// Such a file can be loaded as:
/// ```
/// use bioshell_chem::{ChemErrors};
/// # fn main() -> Result<(), ChemErrors> {
/// use bioshell_core::io::open_file;
/// use bioshell_chem::molecule_from_cif;
/// let reader = open_file("./tests/test_files/EOH.cif")?;
/// let mol = molecule_from_cif(reader)?;
/// assert_eq!(mol.count_atoms(), 9);
/// assert_eq!(mol.count_bonds(), 8);
/// # Ok(())
/// # }
/// ```
pub fn molecule_from_cif<R: BufRead>(reader: R) -> Result<Molecule, ChemErrors> {

    let data_blocks = read_cif_buffer(reader)?;
    let mol_data = &data_blocks[0];     // --- read the first data block, if there are multiple ones, we ignore the rest

    let mol_name: String = mol_data.get_item("_chem_comp.name")
        .ok_or_else(|| ChemErrors::MissingCifField("_chem_comp.name".to_string()))?;
    let mut molecule = Molecule::new(&mol_name);
    molecule.code = mol_data.get_item("_chem_comp.id");

    let atoms_table = CifTable::new(&mol_data, "_chem_comp_atom",
        ["pdbx_ordinal", "atom_id", "model_Cartn_x", "model_Cartn_y", "model_Cartn_z", "type_symbol","charge"])?;
    let mut name_to_idx: HashMap<String, usize> = HashMap::new();
    for [pdbx_ordinal, atom_id, x, y, z, atom_type, charge] in atoms_table.iter() {
        let vx = parse(x, "model_Cartn_x")?;
        let vy = parse(y, "model_Cartn_x")?;
        let vz = parse(z, "model_Cartn_x")?;
        let e: Element = atom_type.parse()?;
        let idx: usize = pdbx_ordinal.parse()
            .map_err(|_| ChemErrors::NumericParsingError("_chem_comp_atom.pdbx_ordinal".to_string(), pdbx_ordinal.to_string()))?;
        let q: i8 = charge.parse()
            .map_err(|_| ChemErrors::NumericParsingError("_chem_comp_atom.charge".to_string(), charge.to_string()))?;
        let mut a = Atom::charged(idx-1, e, q);
        a.set_pos3(vx, vy, vz);
        molecule.add_atom(a)?;
        name_to_idx.insert(atom_id.to_string(), idx-1);
    }

    let bonds_table = CifTable::new(&mol_data, "_chem_comp_bond",
                                    ["atom_id_1", "atom_id_2", "value_order", "pdbx_aromatic_flag"])?;
    for [id_1, id_2, order, is_aromatic] in bonds_table.iter() {
        let idx_1 = name_to_idx.get(id_1)
            .ok_or_else(|| ChemErrors::UnknownAtomName(id_1.to_string()))?;
        let idx_2 = name_to_idx.get(id_2)
            .ok_or_else(|| ChemErrors::UnknownAtomName(id_2.to_string()))?;
        let bond_type = if is_aromatic == "Y" {
            BondType::Aromatic
        } else {
            BondType::from_code(order)
        };
        molecule.bind_atoms(*idx_1, *idx_2, bond_type)?;
    }
    return Ok(molecule);
}

