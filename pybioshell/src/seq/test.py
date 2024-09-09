import residue_types

def test_monomer_type():
    # Accède à un monomer type depuis le module
    monomer1 = residue_types.MonomerType.PeptideLinking
    monomer2 = residue_types.MonomerType.DBetaPeptide
    monomer3 = residue_types.MonomerType.CDeltaLinking
    monomer4 = residue_types.MonomerType.DPeptideCOOH
    monomer5 = residue_types.MonomerType.DNALinking

    
    # Obtient la valeur associée à ce monomer type
    value1 = residue_types.get_monomer_type_value(monomer1)
    value2 = residue_types.get_monomer_type_value(monomer2)
    value3 = residue_types.get_monomer_type_value(monomer3)
    value4 = residue_types.get_monomer_type_value(monomer4)
    value5 = residue_types.get_monomer_type_value(monomer5)
    
    # Affiche le monomer type et sa valeur
    print(f"MonomerType: {monomer1}, Value: {value1}")
    print(f"MonomerType: {monomer2}, Value: {value2}")
    print(f"MonomerType: {monomer3}, Value: {value3}")
    print(f"MonomerType: {monomer4}, Value: {value4}")
    print(f"MonomerType: {monomer5}, Value: {value5}")

if __name__ == "__main__":
    test_monomer_type()