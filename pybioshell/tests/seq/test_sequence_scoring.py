from bioshell import seq

# Built-in matrix
blosum = seq.SubstitutionMatrix("BLOSUM62")

# Custom matrix loaded from a file
custom = seq.SubstitutionMatrix.from_file("test_files/BLOSUM50")

print(blosum.name)  # BLOSUM62
print(custom.name)  # BLOSUM50

print(custom.score("A", "W"))