import sys, re, logging
from collections import defaultdict

from bioshell.seq import Sequence, FastaIterator

def split_fasta_by_species(in_fname: str):
    # Dictionary to accumulate sequences per filename
    sequences_by_file = defaultdict(list)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    filename_sanitizer = re.compile(r'[\\/*?:"<>| ]')

    # Collect sequences in memory, grouped by sanitized filename
    cnt = 0
    for seq in FastaIterator(in_fname):
        fasta_entry = f">{seq.description}\n{seq.seq}"
        sequences_by_file[seq.species].append(fasta_entry)
        cnt += 1
        if cnt % 10000 == 0:
            logging.info(f"{cnt} sequences loaded")

    # Write each group of sequences to its corresponding file
    for fname, entries in sequences_by_file.items():
        fname = filename_sanitizer.sub("_", fname)
        with open(fname + ".fasta", "w") as f:
            f.write("\n".join(entries) + "\n")

if __name__ == "__main__":
    split_fasta_by_species(sys.argv[1])