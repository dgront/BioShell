use bioshell_seq::chemical::{ResidueType, ResidueTypeManager, MonomerType, StandardResidueType, ResidueTypeProperties};

#[test]
fn test_residue_types() {

    // ---------- Iterate over standard amino acid enum types
    let mut n_aa: i8 = 0;
    for srt in StandardResidueType::TYPES {
        if srt.chem_compound_type() == MonomerType::PeptideLinking { n_aa += 1; }
    }
    assert_eq!(n_aa, 21);       // 20 standard amino acids + UNK

    // ---------- Fetch a few standard amino acid types
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

    // ---------- Create a new residue type
    let alm = ResidueType::from_attrs("ALM", StandardResidueType::ALA, MonomerType::PeptideLinking);
    // ---------- Create a new RT manager, it should preload all the 29 standard residue types
    let mut mgr = ResidueTypeManager::get();
    assert_eq!(mgr.count(), 33);

    // ---------- ALA should be already in the manager
    let ala = mgr.by_code3(&String::from("ALA"));
    assert!(ala.is_some());

    // ---------- ... but ALN hasn't been registered yet
    assert!(mgr.by_code3(&String::from("ALN")).is_none());
    mgr.register_residue_type(alm);
    let aln = ResidueType::from_attrs("ALN", StandardResidueType::ALA, MonomerType::PeptideLinking);
    mgr.register_residue_type(aln);
    assert_eq!(mgr.count(), 35);

    // ---------- code1() for ALN
    assert_eq!(mgr.by_code3("ALN").unwrap().parent_type.code1(), 'A');
}
