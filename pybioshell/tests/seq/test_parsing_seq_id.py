from bioshell.seq import parse_sequence_id

test_data = """sp|P0A3D1|PETF_ECOLI 2Fe-2S ferredoxin OS=Escherichia coli OX=562
CYP7598A1_Helminthora_divaricata red algae, Nemaliales) GCA_978020535.1 OZ391815.1_000176.1 30% to Chondrus crispus XP_005715170.1
CYP716A41v2_Bupleurum_chinense (Apiales) PZ062902.1, 98% to CYP716A41_Bupleurum_chinense
""".splitlines()

for desc_line in test_data:
    ids = parse_sequence_id(desc_line)
    print(ids)