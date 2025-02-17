

#[cfg(test)]
mod tests {
    use std::io::BufReader;
    use bioshell_cif::*;
    #[test]
    fn read_cif_file() {
        let mut reader = BufReader::new(ALA_CIF.as_bytes());
        let data_blocks = read_cif_buffer(&mut reader);
        assert!(data_blocks.is_ok());
        let data_blocks = data_blocks.unwrap();
        assert_eq!(data_blocks.len(), 1);
        assert_eq!(data_blocks[0].loop_blocks().count(), 5);
        // println!("{}",data_blocks[0]);

        let mut reader = BufReader::new(BB_.as_bytes());
        let data_blocks = read_cif_buffer(&mut reader);
        assert!(data_blocks.is_ok());
        let data_blocks = data_blocks.unwrap();
        assert_eq!(data_blocks.len(), 1);
        assert_eq!(data_blocks[0].loop_blocks().count(), 1);
        let a_loop = data_blocks[0].loop_blocks().next();
        assert!(a_loop.is_some());
        assert_eq!(a_loop.unwrap().column_names().count(), 9);
        // println!("{}",data_blocks[0]);
    }

    #[test]
    fn test_data_formatter() {
        let cif_block = "data_ALA
_chem_comp.id                                    ALA
_chem_comp.name                                  ALANINE
_chem_comp.type                                  'L-PEPTIDE LINKING'
_chem_comp.pdbx_type                             ATOMP
";
        let expected = "data_ALA
_chem_comp.id        ALA
_chem_comp.name      ALANINE
_chem_comp.pdbx_type ATOMP
_chem_comp.type      'L-PEPTIDE LINKING'
";
        let mut reader = BufReader::new(cif_block.as_bytes());
        let data_blocks = read_cif_buffer(&mut reader);
        assert!(data_blocks.is_ok());
        let data_blocks = data_blocks.unwrap();
        assert_eq!(data_blocks.len(), 1);
        assert_eq!(data_blocks[0].name(), "ALA");
        let out = format!("{}", data_blocks[0]);
        let mut lines: Vec<&str> = out.lines().collect::<Vec<&str>>()[1..].to_vec();
        lines.sort();
        let out = lines.iter().fold(String::from("data_ALA"), |a, b| a + b + "\n");
        assert_eq!(out, expected);
    }

    #[test]
    fn test_loop_formatter() {
        let cif_block = "data_some_name
    loop_
    _first_column
    _second_column
    'value A' 1
    'value B' 2
    'value C' 2

    ";
        let expected = "loop_
_first_column
_second_column
'value A' 1
'value B' 2
'value C' 2

";
        let mut reader = BufReader::new(cif_block.as_bytes());
        let data_blocks = read_cif_buffer(&mut reader);
        assert!(data_blocks.is_ok());
        let data_blocks = data_blocks.unwrap();
        assert_eq!(data_blocks.len(), 1);
        assert_eq!(data_blocks[0].name(), "some_name");
        let a_loop = data_blocks[0].loop_blocks().next().unwrap();
        let out = format!("{}", a_loop);

        assert_eq!(out, expected);
    }

    #[test]
    fn build_cif_loop() {
        let mut data_loop = CifLoop::new(&["_symmetry_equiv_pos_site_id"]);
        assert_eq!(data_loop.count_columns(), 1);
        let _ = data_loop.add_column("_symmetry_equiv_pos_as_xyz");
        assert_eq!(data_loop.count_columns(), 2);
        let _ = data_loop.add_data_row(vec!["1", "x,y,z"].iter().map(|&s| s.to_string()).collect());
        let _ = data_loop.add_data_row(vec!["2", "x,y,z"].iter().map(|&s| s.to_string()).collect());
        assert_eq!(data_loop.count_rows(), 2);
        *data_loop.entry_mut(0, "_symmetry_equiv_pos_as_xyz").unwrap() = "-x,-y,-z".to_string();
        println!("{}", data_loop);
        let txt = format!("{}", data_loop);
        let expected = "loop_
_symmetry_equiv_pos_site_id
_symmetry_equiv_pos_as_xyz
1 -x,-y,-z
2 x,y,z

";
        assert_eq!(txt, expected);
    }

    #[test]
    fn test_cif_table() -> Result<(), CifError> {
        let cif_str = "data_loop
loop_
_symmetry_equiv_pos.site_id
_symmetry_equiv_pos.as_xyz
1 -x,-y,-z
2 x,y,z
";
        let mut reader = BufReader::new(cif_str.as_bytes());
        let data_block = &read_cif_buffer(&mut reader)?[0];
        let cif_table = CifTable::new(data_block, "_symmetry_equiv_pos",["site_id", "as_xyz"])?;
        assert_eq!(cif_table.iter().count(), 2);
        Ok(())
    }

    #[test]
    fn parse_edge_cases() {

        // ---------- multiline data value
        let mut reader = BufReader::new(edge_multiline_value.as_bytes());
        let data_blocks = read_cif_buffer(&mut reader);
        assert!(data_blocks.is_ok());
        let data_blocks = data_blocks.unwrap();
        assert_eq!(data_blocks.len(), 1);
        assert_eq!(data_blocks[0].data_items().len(), 1);

        // ---------- data value in the next line after data key
        let mut reader = BufReader::new(edge_data_in_next_line.as_bytes());
        let data_blocks = read_cif_buffer(&mut reader).unwrap();
        assert_eq!(data_blocks.len(), 1);

        // ---------- data name as a single token in a line just after a loop
        let mut reader = BufReader::new(edge_data_name_after_loop.as_bytes());
        let data_blocks = read_cif_buffer(&mut reader).unwrap();
        assert_eq!(data_blocks.len(), 1);
    }

    #[allow(non_upper_case_globals)]
    static edge_multiline_value: &'static str = r#"data_ALA
_chem_comp.name
;
multiline
;
    #"#;

    #[allow(non_upper_case_globals)]
    static edge_data_in_next_line: &'static str = r#"data_1crr
_struct.entry_id                  1CRR
_struct.title
'THE SOLUTION STRUCTURE AND DYNAMICS OF RAS P21.'
_struct.pdbx_model_details        ?
"#;

    #[allow(non_upper_case_globals)]
    static edge_data_name_after_loop: &'static str = r#"data_2jnw
_struct.entry_id                  2JNW
loop_
_atom_type.symbol
C
H
N
_pdbx_nmr_refine.details
;XPLOR-NIH was used.
;
"#;

    #[allow(non_upper_case_globals)]
    static ALA_CIF: &'static str = r#"data_ALA
#
_chem_comp.id                                    ALA
_chem_comp.name                                  ALANINE
_chem_comp.type                                  "L-PEPTIDE LINKING"
_chem_comp.pdbx_type                             ATOMP
_chem_comp.formula                               "C3 H7 N O2"
_chem_comp.mon_nstd_parent_comp_id               ?
_chem_comp.pdbx_synonyms                         ?
_chem_comp.pdbx_formal_charge                    0
_chem_comp.pdbx_initial_date                     1999-07-08
_chem_comp.pdbx_modified_date                    2011-06-04
_chem_comp.pdbx_ambiguous_flag                   N
_chem_comp.pdbx_release_status                   REL
_chem_comp.pdbx_replaced_by                      ?
_chem_comp.pdbx_replaces                         ?
_chem_comp.formula_weight                        89.093
_chem_comp.one_letter_code                       A
_chem_comp.three_letter_code                     ALA
_chem_comp.pdbx_model_coordinates_details        ?
_chem_comp.pdbx_model_coordinates_missing_flag   N
_chem_comp.pdbx_ideal_coordinates_details        ?
_chem_comp.pdbx_ideal_coordinates_missing_flag   N
_chem_comp.pdbx_model_coordinates_db_code        ?
_chem_comp.pdbx_subcomponent_list                ?
_chem_comp.pdbx_processing_site                  RCSB
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_align
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
_chem_comp_atom.model_Cartn_x
_chem_comp_atom.model_Cartn_y
_chem_comp_atom.model_Cartn_z
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
_chem_comp_atom.pdbx_component_atom_id
_chem_comp_atom.pdbx_component_comp_id
_chem_comp_atom.pdbx_ordinal
ALA N   N   N 0 1 N N N 2.281  26.213 12.804 -0.966 0.493  1.500  N   ALA 1
ALA CA  CA  C 0 1 N N S 1.169  26.942 13.411 0.257  0.418  0.692  CA  ALA 2
ALA C   C   C 0 1 N N N 1.539  28.344 13.874 -0.094 0.017  -0.716 C   ALA 3
ALA O   O   O 0 1 N N N 2.709  28.647 14.114 -1.056 -0.682 -0.923 O   ALA 4
ALA CB  CB  C 0 1 N N N 0.601  26.143 14.574 1.204  -0.620 1.296  CB  ALA 5
ALA OXT OXT O 0 1 N Y N 0.523  29.194 13.997 0.661  0.439  -1.742 OXT ALA 6
ALA H   H   H 0 1 N N N 2.033  25.273 12.493 -1.383 -0.425 1.482  H   ALA 7
ALA H2  HN2 H 0 1 N Y N 3.080  26.184 13.436 -0.676 0.661  2.452  H2  ALA 8
ALA HA  HA  H 0 1 N N N 0.399  27.067 12.613 0.746  1.392  0.682  HA  ALA 9
ALA HB1 1HB H 0 1 N N N -0.247 26.699 15.037 1.459  -0.330 2.316  HB1 ALA 10
ALA HB2 2HB H 0 1 N N N 0.308  25.110 14.270 0.715  -1.594 1.307  HB2 ALA 11
ALA HB3 3HB H 0 1 N N N 1.384  25.876 15.321 2.113  -0.676 0.697  HB3 ALA 12
ALA HXT HXT H 0 1 N Y N 0.753  30.069 14.286 0.435  0.182  -2.647 HXT ALA 13
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
_chem_comp_bond.pdbx_ordinal
ALA N   CA  SING N N 1
ALA N   H   SING N N 2
ALA N   H2  SING N N 3
ALA CA  C   SING N N 4
ALA CA  CB  SING N N 5
ALA CA  HA  SING N N 6
ALA C   O   DOUB N N 7
ALA C   OXT SING N N 8
ALA CB  HB1 SING N N 9
ALA CB  HB2 SING N N 10
ALA CB  HB3 SING N N 11
ALA OXT HXT SING N N 12
#
loop_
_pdbx_chem_comp_descriptor.comp_id
_pdbx_chem_comp_descriptor.type
_pdbx_chem_comp_descriptor.program
_pdbx_chem_comp_descriptor.program_version
_pdbx_chem_comp_descriptor.descriptor
ALA SMILES           ACDLabs              10.04 "O=C(O)C(N)C"
ALA SMILES_CANONICAL CACTVS               3.341 "C[C@H](N)C(O)=O"
ALA SMILES           CACTVS               3.341 "C[CH](N)C(O)=O"
ALA SMILES_CANONICAL "OpenEye OEToolkits" 1.5.0 "C[C@@H](C(=O)O)N"
ALA SMILES           "OpenEye OEToolkits" 1.5.0 "CC(C(=O)O)N"
ALA InChI            InChI                1.03  "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1"
ALA InChIKey         InChI                1.03  QNAYBMKLOCPYGJ-REOHCLBHSA-N
#
loop_
_pdbx_chem_comp_identifier.comp_id
_pdbx_chem_comp_identifier.type
_pdbx_chem_comp_identifier.program
_pdbx_chem_comp_identifier.program_version
_pdbx_chem_comp_identifier.identifier
ALA "SYSTEMATIC NAME" ACDLabs              10.04 L-alanine
ALA "SYSTEMATIC NAME" "OpenEye OEToolkits" 1.5.0 "(2S)-2-aminopropanoic acid"
#
loop_
_pdbx_chem_comp_audit.comp_id
_pdbx_chem_comp_audit.action_type
_pdbx_chem_comp_audit.date
_pdbx_chem_comp_audit.processing_site
ALA "Create component"  1999-07-08 RCSB
ALA "Modify descriptor" 2011-06-04 RCSB
#"#;

    const BB_: &str = "data_bb_
loop_
_res_name
_atom_a_name
_atom_b_name
_atom_c_name
_atom_d_name
_c_d_bond_length
_b_c_d_planar_angle
_a_b_c_d_dihedral_angle
_dihedral_angle_name
'bb ' ' N  ' ' CA ' ' C  ' ' N  ' 1.328685 114.0  180.0 psi
'bb ' ' CA ' ' C  ' ' N  ' ' CA ' 1.458001 123.0  180.0 omega
'bb ' ' C  ' ' N  ' ' CA ' ' C  ' 1.523258 110.0 -180.0 phi
'bb ' ' N  ' ' CA ' ' C  ' ' O  ' 1.231015 121.0  180.0 -
#";
}