use bioshell_numerical::Vec3;
use bioshell_pdb::pdb_atom::PdbAtom;


fn main() {
    let mut atom = PdbAtom::new();
    atom.atom_serial_no = Some(2831);
    atom.atom_symbol = "O".to_string();//O
    atom.atom_position = "E".to_string();//E
    atom.atom_no_in_the_branch = Some(1);//1
    atom.residue_name = "GLN".to_string();//
    atom.chain_name = "A".to_string();
    atom.residue_no = Some(294);
    atom.coordinate = Vec3::new(-27.117, 12.343, 28.479);
    atom.occupancy = Some(1.0);
    atom.temperature_factor = Some(9.58);
    atom.charge_of_the_atom = "O".to_string();

    //////////////////////////2831,O,E,1,0,,GLN,A,294,,-27.117,12.343,28.479,1,9.58,,,O
    //assert_eq!(csv_string, "2831,O,E,1,#,#,GLN,A,294,#,-27.117,12.343,28.479,1,9.58,#,O");
    let csv_string = atom.to_csv_string();
}