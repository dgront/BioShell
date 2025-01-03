#[cfg(test)]
mod tests {
    use std::io::BufReader;
    use bioshell_cif::{is_cif_file, read_cif_buffer};
    use bioshell_pdb::{Deposit, UnitCell};

    #[allow(non_upper_case_globals)]
    const cif_2gb1:  &str = include_str!("./test_files/2gb1.cif");

    #[allow(non_upper_case_globals)]
    const cif_2fdo:  &str = include_str!("./test_files/2fdo.cif");

    #[test]
    fn load_2gb1_from_cif() {
        let reader = BufReader::new(cif_2gb1.as_bytes());
        let strctr = Deposit::from_cif_reader(reader).unwrap().structure().unwrap();
        assert_eq!(strctr.count_atoms(), 855);
        assert_eq!(strctr.count_chains(), 1);
        assert_eq!(strctr.count_models(), 1);
        assert_eq!(strctr.count_residues(), 56);
    }

    #[test]
    fn load_2fdo_from_cif() {
        let reader = BufReader::new(cif_2fdo.as_bytes());
        let deposit = Deposit::from_cif_reader(reader).unwrap();
        let strctr = deposit.structure().unwrap();
        assert_eq!(strctr.count_atoms(), 1456);
        assert_eq!(strctr.count_chains(), 2);
        assert_eq!(strctr.count_models(), 1);
        assert_eq!(strctr.count_residues(), 214);
        assert!(deposit.methods.len() > 0);
        assert!(deposit.resolution.is_some());
        assert!(deposit.r_factor.is_some());
        assert!(deposit.r_free.is_some());
    }

    #[test]
    fn unit_cell_from_cif() {
        let missing_data = "data_cryst_cell
            _cell.length_b                         86.70
            _cell.length_c                         46.27
            _cell.angle_alpha                      90.00
            _cell.angle_beta                       90.00
            _cell.angle_gamma                      90.00
            _symmetry.space_group_name_H-M         'C 1 21 1'
            _cell.Z_PDB                            1
        ";
        let data_block = &read_cif_buffer(&mut BufReader::new(missing_data.as_bytes())).unwrap();
        let data_block = &data_block[0];
        let uc = UnitCell::from_cif_data(data_block);
        assert!(uc.is_err());
    }

    #[test]
    fn test_if_cif() {
        let try_2gb1 = is_cif_file("./tests/test_files/2gb1.cif");
        assert!(try_2gb1.is_ok());
        assert!(try_2gb1.unwrap());
        let try_2gb1 = is_cif_file("./tests/test_files/2gb1.pdb");
        assert!(try_2gb1.is_ok());
        assert!(!try_2gb1.unwrap());
        let try_2gb1 = is_cif_file("./tests/test_files/2gb1.blabla");
        assert!(try_2gb1.is_err());
    }
}