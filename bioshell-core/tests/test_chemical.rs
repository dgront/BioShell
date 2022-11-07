use bioshell_core::chemical::{ResidueType, ResidueTypeManager, MonomerType, StandardResidueType, ResidueTypeProperties};

#[test]
fn test_residue_types() {

    // ---------- Create a few standard amino acid types
    let cys: StandardResidueType = StandardResidueType::CYS;
    let trp: StandardResidueType = StandardResidueType::TRP;

    // ---------- Check their properties
    assert_eq!(cys, StandardResidueType::CYS);
    assert_eq!(cys.code1(),'C');
    assert_eq!(cys.code3(),"CYS");
    assert_eq!(trp.code3(),"TRP");
    assert_eq!(cys.chem_compound_type(),MonomerType::PeptideLinking);
    assert_eq!(trp.chem_compound_type(),cys.chem_compound_type());

    // ---------- Check conversion 'C' -> StandardResidueType::CYS
    assert_eq!(StandardResidueType::try_from('C').unwrap(), cys);
    // ---------- This conversion should fail!
    let result = StandardResidueType::try_from('B');
    assert!(result.is_err());

    let aln = ResidueType{parent_type: StandardResidueType::ALA, code3:String::from("ALN"), chem_compound_type: MonomerType::PeptideLinking};
    let mut mgr = ResidueTypeManager::new();
    assert_eq!(mgr.count(), 29);
    mgr.register_residue_type(aln);
}
