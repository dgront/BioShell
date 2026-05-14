use std::collections::HashMap;
use std::io::BufRead;
use bioshell_cif::{read_cif_buffer};
use bioshell_cif::CifTable;
use crate::{Atom, BondType, ChemErrors, Element, Molecule};

pub fn molecule_from_cif<R: BufRead>(reader: R) -> Result<Molecule, ChemErrors> {

    let data_blocks = read_cif_buffer(reader)?;
    let mol_data = &data_blocks[0];     // --- read the first data block, if there are multiple ones, we ignore the rest

    let mol_name: String = mol_data.get_item("_chem_comp.name")
        .ok_or_else(|| ChemErrors::MissingCifField("_chem_comp.name".to_string()))?;
    let mut molecule = Molecule::new(&mol_name);

    let atoms_table = CifTable::new(&mol_data, "_chem_comp_atom",
                    ["pdbx_ordinal", "atom_id", "type_symbol","charge"])?;
    let mut name_to_idx: HashMap<String, usize> = HashMap::new();
    for [pdbx_ordinal, atom_id, atom_type, charge] in atoms_table.iter() {
        let e: Element = atom_type.parse()?;
        let idx: usize = pdbx_ordinal.parse()
            .map_err(|_| ChemErrors::CifParsingError("_chem_comp_atom.pdbx_ordinal".to_string(), pdbx_ordinal.to_string()))?;
        let q: i8 = charge.parse()
            .map_err(|_| ChemErrors::CifParsingError("_chem_comp_atom.charge".to_string(), charge.to_string()))?;
        molecule.add_atom(Atom::charged(idx-1, e.atomic_number(), q))?;
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