use std::env;
use bioshell_pdb::{Pdb, PdbAtom};
use std::io::{self, Read};

fn main(){
    let atom_line = "HETATM 3169  PB  BIS A 499     -10.950  43.372   1.324  1.00 18.00           P  ";
    let pdbAtom = PdbAtom::parse(atom_line);

    println!("1. atom_serial_no_1: {:?}", pdbAtom.atom_serial_no_1);//Some(1)
    println!("2. atom_symbol_2{:?}", pdbAtom.atom_symbol_2);//N
    println!("3. atom_position_3{:?}", pdbAtom.atom_position_3);// ""
    println!("4. atom_no_in_the_branch_4{:?}", pdbAtom.atom_no_in_the_branch_4);// Some(0)
    println!("5. {:?}", pdbAtom.connected_to_atom_no_5);// Some(0)
    println!("6. {:?}", pdbAtom.alt_loc_indicator_6);//""
    println!("7. {:?}", pdbAtom.residue_name_7);//VAL
    println!("8. {:?}", pdbAtom.chain_name_8);//A
    println!("9. {:?}", pdbAtom.residue_no_9);//Some(1)
    println!("10. {:?}", pdbAtom.insertion_code_10);//""
    println!("11. {:?}", pdbAtom.coordinate_11);//[19.323 29.727 42.781] >0,0,0<
    println!("12. {:?}", pdbAtom.occupancy_12);//Some(1.0)
    println!("13. {:?}", pdbAtom.temperature_factor_13);//Some(49.05)
    println!("14. {:?}", pdbAtom.segment_identifier_14);//""
    println!("15. {:?}", pdbAtom.segment_identifier_symbol_15);//"N"
    println!("16. {:?}", pdbAtom.charge_of_the_atom_16);//""
    println!("17. {:?}", pdbAtom.is_hetero_atom);//false
    println!("18. {:?}", pdbAtom.protein_name);//""
}

/*
fn main(){
    // let atom_line="ATOM      1  N   VAL A   1      19.323  29.727  42.781  1.00 49.05           N  ";
    let atom_line="ATOM      1  N   VAL A   1      19.323  29.727  42.781  1.00 49.05           N  ";
    let pdbAtom = PdbAtom::parse(atom_line);

    println!("1. atom_serial_no_1: {:?}", pdbAtom.atom_serial_no_1);//Some(1)
    println!("2. atom_symbol_2{:?}", pdbAtom.atom_symbol_2);//N
    println!("3. atom_position_3{:?}", pdbAtom.atom_position_3);// ""
    println!("4. atom_no_in_the_branch_4{:?}", pdbAtom.atom_no_in_the_branch_4);// Some(0)
    println!("5. {:?}", pdbAtom.connected_to_atom_no_5);// Some(0)
    println!("6. {:?}", pdbAtom.alt_loc_indicator_6);//""
    println!("7. {:?}", pdbAtom.residue_name_7);//VAL
    println!("8. {:?}", pdbAtom.chain_name_8);//A
    println!("9. {:?}", pdbAtom.residue_no_9);//Some(1)
    println!("10. {:?}", pdbAtom.insertion_code_10);//""
    println!("11. {:?}", pdbAtom.coordinate_11);//[19.323 29.727 42.781] >0,0,0<
    println!("12. {:?}", pdbAtom.occupancy_12);//Some(1.0)
    println!("13. {:?}", pdbAtom.temperature_factor_13);//Some(49.05)
    println!("14. {:?}", pdbAtom.segment_identifier_14);//""
    println!("15. {:?}", pdbAtom.segment_identifier_symbol_15);//"N"
    println!("16. {:?}", pdbAtom.charge_of_the_atom_16);//""
    println!("17. {:?}", pdbAtom.is_hetero_atom);//false
    println!("18. {:?}", pdbAtom.protein_name);//""
}
*/


/*
fn main()-> std::io::Result<()> {
    if let Ok(current_dir) = env::current_dir() {
        println!("Current working directory: {}", current_dir.display());
    }
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
    let mut buffer = [0; 1];
    println!("Press any key to continue...");
    io::stdin().read_exact(&mut buffer).unwrap();

    Ok(())

    
}
*/
