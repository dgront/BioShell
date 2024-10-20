#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::io::BufReader;
    use bioshell_cif::read_cif_buffer;
    use bioshell_pdb::{Entity, EntityType, load_cif_reader, PDBError, PolymerEntity, PolymerEntityType};
    use bioshell_pdb::PolymerEntityType::{PolypeptideL, DNA};

    #[allow(non_upper_case_globals)]
    const cif_2gb1:  &str = include_str!("./test_files/2gb1.cif");
    #[allow(non_upper_case_globals)]
    const cif_5edw:  &str = include_str!("./test_files/5edw.cif");

    #[test]
    fn test_entity_parsing() {
        let test_data = [(cif_2gb1, 1), (cif_5edw, 6)];
        for (data, num_entities) in test_data.iter() {
            let reader = BufReader::new(data.as_bytes());
            let cif_data = read_cif_buffer(reader).unwrap();
            assert_eq!(cif_data.len(), 1);
            let entities = Entity::from_cif_data(&cif_data[0]).unwrap();
            assert_eq!(entities.len(), *num_entities);
        }
    }

    #[test]
    fn check_5edw_entities() {
        let reader = BufReader::new(cif_5edw.as_bytes());
        let cif_data = read_cif_buffer(reader).unwrap();
        let entities = Entity::from_cif_data(&cif_data[0]).unwrap();

        let mut keys: Vec<&str> = entities.keys().map(|k| k.as_str()).collect();
        keys.sort();
        let expected_keys = vec!["1", "2", "3", "4", "5", "6"];
        assert_eq!(keys, expected_keys);

        let mut counts: HashMap<EntityType, usize> = [
            (EntityType::Polymer(PolypeptideL), 0),
            (EntityType::Polymer(DNA), 0),
            (EntityType::NonPolymer, 0),
            (EntityType::Water, 0)].iter().cloned().collect();
        for (_,entity) in entities.iter() {
            counts.entry(entity.entity_type()).and_modify(|e| *e += 1);
        }
        assert_eq!(counts[&EntityType::Polymer(PolypeptideL)], 1);
        assert_eq!(counts[&EntityType::Polymer(DNA)], 2);
        assert_eq!(counts[&EntityType::NonPolymer], 2);
        assert_eq!(counts[&EntityType::Water], 1);
    }

    #[test]
    fn load_polymer_entities() -> Result<(), PDBError> {
        let reader = BufReader::new(cif_5edw.as_bytes());
        let cif_data = read_cif_buffer(reader).unwrap();
        let entities = PolymerEntity::from_cif_data(&cif_data[0])?;
        assert_eq!(entities.len(), 3);
        let result = vec![(341, PolypeptideL), (19, DNA), (13, DNA)];
        for entity in entities {
            let id = entity.entity_id().parse::<usize>().unwrap();
            assert_eq!(result[id-1].0, entity.monomer_sequence().len());
            assert_eq!(&result[id-1].1, entity.entity_type());
        }
        Ok(())
    }

    #[test]
    fn entities_5edw_structure() -> Result<(), PDBError> {
        let cif_data = include_str!("../tests/test_files/5edw.cif");
        let strctr = load_cif_reader(cif_data.as_bytes())?;

        assert_eq!(strctr.entities().count(), 6);

        let mut protein_entities = 0;
        let mut dna_entities = 0;
        for (id, entity) in strctr.entities() {
            if entity.entity_type() == EntityType::Polymer(PolypeptideL) {
                protein_entities += 1;
                assert_eq!(entity.chain_ids(), &["A"]);
                assert_eq!(entity.monomer_sequence().len(), 341);
            }
            if entity.entity_type() == EntityType::Polymer(DNA) { dna_entities += 1; }
        }
        assert_eq!(protein_entities, 1);
        assert_eq!(dna_entities, 2);
        Ok(())
    }

    #[test]
    fn entities_2fdo_structure() -> Result<(), PDBError> {
        let cif_data = include_str!("../tests/test_files/2fdo.cif");
        let strctr = load_cif_reader(cif_data.as_bytes())?;

        assert_eq!(strctr.entities().count(), 2);
        let entity = strctr.entity("1");
        assert_eq!(entity.entity_type(), EntityType::Polymer(PolypeptideL));

        // --- full chains, including water molecules
        assert_eq!(strctr.chain_residue_ids("B").len(), 108);
        assert_eq!(strctr.chain_residue_ids("A").len(), 102);

        // --- amino acids only
        let n_aa = strctr.chain_residue_ids("B").iter()
            .filter(|ri| strctr.residue_type(ri).unwrap().chem_compound_type.is_peptide_linking())
            .count();
        assert_eq!(n_aa, 94);
        let n_aa = strctr.chain_residue_ids("A").iter()
            .filter(|ri| strctr.residue_type(ri).unwrap().chem_compound_type.is_peptide_linking())
            .count();
        assert_eq!(n_aa, 93);
        Ok(())
    }
}


