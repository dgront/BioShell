#[cfg(test)]
mod tests_bond_types {
    use bioshell_chem::BondType;

    #[test]
    fn test_formatting() {
        assert_eq!(BondType::from_code("1"), BondType::Single);
        assert_eq!(BondType::from_code("2"), BondType::Double);
        assert_eq!(BondType::from_code("3"), BondType::Triple);
        assert_eq!(BondType::from_code("ar"), BondType::Aromatic);
        assert_eq!(BondType::from_code("hy"), BondType::Hydrogen);
        assert_eq!(BondType::from_code("du"), BondType::Dummy);

        assert_eq!(BondType::from_code("DOUB"), BondType::Double);
        assert_eq!(BondType::from_code("TRIP"), BondType::Triple);
        assert_eq!(BondType::from_code("SING"), BondType::Single);

        assert_eq!(BondType::from_code("unknown-code"), BondType::Unknown);

        assert_eq!(BondType::Double.mol2_code(), "2");
        assert_eq!(BondType::Aromatic.mol2_code(), "ar");
        assert_eq!(BondType::Dummy.mol2_code(), "du");

        assert_eq!(BondType::Double.name(), "double");
        assert_eq!(BondType::Aromatic.to_string(), "aromatic");
    }

    #[test]
    fn test_bond_order() {
        assert_eq!(BondType::Single.order(), 1.0);
        assert_eq!(BondType::Double.order(), 2.0);
        assert_eq!(BondType::Triple.order(), 3.0);
        assert_eq!(BondType::Aromatic.order(), 1.5);
    }
}