from bioshell.seq import Sequence, FastaIterator

n_seq = 0
for seq in FastaIterator("test_files/fdx_2fe_2s.fasta"):
    n_seq += 1
    print(seq)
    assert isinstance(seq, Sequence)
assert n_seq == 32