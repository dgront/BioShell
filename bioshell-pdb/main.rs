use bioshell_pdb::Pdb;

fn main()-> std::io::Result<()> {
    let pdb_file_path = "bioshell-pdb/tests/test_files/16pk.pdb";
    let pdb = Pdb::from_file(pdb_file_path).unwrap();

    println!("Total atoms : {:?}", pdb.get_atoms_list().len());//3167

    let first_atom = &pdb.get_atoms_list()[727];
    //"ATOM    728 HH11AARG A  98       7.165  44.290  18.691  0.50  0.00           H"
    println!("atom_serial_no : {:?}",&first_atom.atom_serial_no);// Some(0)
    println!("atom_symbol : {:?}",&first_atom.atom_symbol);// "N"
    println!("alt_loc_indicator : {:?}",&first_atom.alt_loc_indicator); //" "
    println!("residue_name : {:?}",&first_atom.residue_name);// "VAL");
    println!("chain_name : {:?}",&first_atom.chain_name); //"A"
    println!("residue_no : {:?}",&first_atom.residue_no);// Some(1));
    println!("insertion_code : {:?}",&first_atom.insertion_code);// " ");
    println!("x : {:?}",first_atom.get_coordinate().x);// 11.54
    println!("y : {:?}",first_atom.get_coordinate().y);// 11.88
    println!("z : {:?}",first_atom.get_coordinate().z);// 7.95
    println!("occupancy : {:?}",&first_atom.occupancy);// Some(1.0));
    println!("temperature_factor : {:?}",&first_atom.temperature_factor);// Some(0.0));
    println!("charge_of_the_atom : {:?}",&first_atom.charge_of_the_atom);// " "

    //let _protein_name = first_atom.get_protein_name(); //"4HHB"

    Ok(())
}