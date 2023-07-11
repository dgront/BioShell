use bioshell_pdb::Pdb;

fn main()-> std::io::Result<()> {
    let pdb_file_path = "bioshell-pdb/tests/test_files/16pk.pdb";
    let pdb = Pdb::from_file(pdb_file_path).unwrap();

    let _atoms_list_length = pdb.get_atoms_list().len();//3167

    let first_atom = &pdb.get_atoms_list()[727];
    //"ATOM    728 HH11AARG A  98       7.165  44.290  18.691  0.50  0.00           H"
    let _atom_serial_no = first_atom.get_atom_serial_no();// Some(0)
    let _atom_symbol = first_atom.get_atom_symbol();// "N"
    let _alt_loc_indicator = first_atom.get_alt_loc_indicator(); //" "
    let _residue_name = first_atom.get_residue_name();// "VAL");
    let _chain_id = first_atom.get_chain_name(); //"A"
    let _residue_no = first_atom.get_residue_no();// Some(1));
    let _insertion_code = first_atom.get_insertion_code();// " ");
    let _x = first_atom.get_coordinate().x;// 11.54
    let _y = first_atom.get_coordinate().y;// 11.88
    let _z = first_atom.get_coordinate().z;// 7.95
    let _occupancy = first_atom.get_occupancy();// Some(1.0));
    let _temperature_factor = first_atom.get_temperature_factor();// Some(0.0));
    let _atom_symbol = first_atom.get_atom_symbol(); //"N"
    let _charge_of_the_atom = first_atom.get_charge_of_the_atom();// " "
    let _protein_name = first_atom.get_protein_name(); //"4HHB"

    Ok(())
}