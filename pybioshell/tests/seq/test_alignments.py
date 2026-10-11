from bioshell.seq import Sequence, SubstitutionMatrix, global_alignment

fdx1 = Sequence("> 2506930501 Pseudomonas syringae CC1557 (CC1557)",
            "MISLLVDGEAVEPWQGDTFMVSLLTRVSALRQSEFDPGRRAGFCLMGACQDCWAH")

fdx2 = Sequence("> 2632118236 Burkholderia pseudomallei TSV202",
        "MSEPIRFYHRRAIREVSGANVTRTVLQYLREDAHCTGTKEGCAEGDCGACTVVVGELNDAGGVEF")

blosum62 = SubstitutionMatrix("BLOSUM62")
ali = global_alignment(fdx1, fdx2, blosum62, -10, -1)
aligned_query, aligned_template = ali.aligned_sequences(fdx1, fdx2)
assert aligned_query.seq == "MISLL--VDGEAVEPWQGDTFMVSLLTRVSALRQSEFDPGRRAGFCLMGACQDCWA------------H"
assert aligned_template.seq == "MSEPIRFYHRRAIREVSGANVTRTVL---QYLREDAHCTGTKEG-CAEGDCGACTVVVGELNDAGGVEF"

print(f"calculated alignment:\n{aligned_query.seq}\n{aligned_template.seq}")
