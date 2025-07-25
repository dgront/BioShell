from bioshell.seq import Sequence

seq = Sequence("1clf:A", "AYKIADSCVSCGACASECPVNAISQGDSIFVIDADTCIDCGNCANVCPVGAPVQE")
print(seq)
assert len(seq) == 55

seq = Sequence("sp|P0A3D1|PETF_ECOLI 2Fe-2S ferredoxin OS=Escherichia coli OX=562",
        "MVFVKCGIPVDYVCPGKEVLHQCGHCPDRGEESAMKGVVKIANTDETTVAGELWVCEYTNDRIGEKLAVKEYGEPDVILLTRQGQCGGRVLLTVRGQKAEKEKENVTVISNYPEGPD")
assert seq.species == "Escherichia coli"
assert seq.id == "sp|P0A3D1|PETF_ECOLI"