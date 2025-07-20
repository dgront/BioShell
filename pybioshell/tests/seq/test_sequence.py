from bioshell.seq import Sequence

seq = Sequence("> 1clf:A", "AYKIADSCVSCGACASECPVNAISQGDSIFVIDADTCIDCGNCANVCPVGAPVQE")
print(seq)
assert len(seq) == 55
