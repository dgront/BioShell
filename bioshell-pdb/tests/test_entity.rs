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
        let expected_keys = vec!["1", "2", "3", "4", "5", "6"];
        for (i, entity) in entities.iter().enumerate() {
            assert_eq!(entity.id(), expected_keys[i]);
        }
        let mut counts: HashMap<EntityType, usize> = [
            (EntityType::Polymer, 0),
            (EntityType::NonPolymer, 0),
            (EntityType::Water, 0)].iter().cloned().collect();
        for entity in entities.iter() {
            counts.entry(entity.entity_type()).and_modify(|e| *e += 1);
        }
        assert_eq!(counts[&EntityType::Polymer], 3);
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

        let mut polymer_entities = 0;
        for (id, entity) in strctr.entities() {
            if entity.entity_type() == EntityType::Polymer { polymer_entities += 1; }
        }
        assert_eq!(polymer_entities, 3);
        Ok(())
    }
}


