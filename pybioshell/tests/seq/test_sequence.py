from bioshell.seq import Sequence, parse_sequence_id

seq = Sequence("1clf:A", "AYKIADSCVSCGACASECPVNAISQGDSIFVIDADTCIDCGNCANVCPVGAPVQE")
assert len(seq) == 55

seq = Sequence("sp|P0A3D1|PETF_ECOLI 2Fe-2S ferredoxin [Escherichia coli] OX=562",
        "MVFVKCGIPVDYVCPGKEVLHQCGHCPDRGEESAMKGVVKIANTDETTVAGELWVCEYTNDRIGEKLAVKEYGEPDVILLTRQGQCGGRVLLTVRGQKAEKEKENVTVISNYPEGPD")
print(seq.full_id(False))
assert seq.first_id(False) == "sp|P0A3D1"
assert seq.full_id(False) == "sp|P0A3D1|PETF_ECOLI|taxid=562 [organism=Escherichia coli]"

expected_id_types = ["SwissProt", "UniProtEntry", "TaxId", "Organism"]
ids_list = parse_sequence_id(seq.description)
for id in ids_list:
    print(id.kind, id.value)
    assert id.kind in expected_id_types
    if id.kind == "Organism":
        assert id.value == "Escherichia coli"