#[cfg(test)]
mod tests {
    use std::io::BufReader;
    use bioshell_cif::read_cif_buffer;
    use bioshell_pdb::{load_cif_reader, UnitCell};

    #[allow(non_upper_case_globals)]
    const cif_2gb1:  &str = include_str!("./test_files/2gb1.cif");
    #[test]
    fn load_2gb1_from_cif() {
        let reader = BufReader::new(cif_2gb1.as_bytes());
        let strctr = load_cif_reader(reader).unwrap();
        assert_eq!(strctr.count_atoms(), 855);
        assert_eq!(strctr.count_residues(), 56);
        assert_eq!(strctr.count_chains(), 1);
        assert_eq!(strctr.count_models(), 1);
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
        let data_block = read_cif_buffer(&mut BufReader::new(missing_data.as_bytes())).unwrap();
        let data_block= &data_block[0];
        let uc = UnitCell::from_cif_data(data_block);
        assert!(uc.is_err());
    }
}