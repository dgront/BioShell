use bioshell_seq::scoring::{SubstitutionMatrix, SubstitutionMatrixList};

#[test]
fn test_predefined_matrix() {
    let blosum80 = SubstitutionMatrix::load(SubstitutionMatrixList::BLOSUM80);
    assert_eq!(blosum80.score_by_aa('X' as u8, 'X' as u8), -1);
    assert_eq!(blosum80.score_by_aa('C' as u8, 'C' as u8), 13);
    assert_eq!(blosum80.score_by_aa('W' as u8, 'W' as u8), 16);
    assert_eq!(blosum80.score_by_aa('A' as u8, 'W' as u8), -5);
}