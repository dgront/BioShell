An iterator over sequences in a FASTA file.

Yields `Sequence` objects one by one; reading only a fragment
of the input file at a time. The sequences are parsed with one of the following ``parsing_mode`` strategies:

- ``"raw"`` - as is, without any quality control
- ``"protein"`` - characters other than denoting amino acids will be removed
- ``"protein*"`` - as in  ``"protein"`` but stop codon denoted as ``'*'`` is allowed in a sequence
- ``"nucleic"`` - characters other than denoting nucleic acids will be removed

Example
--------

.. code-block:: python

    from bioshell.seq import FastaIterator

    n_seq = 0
    for seq in FastaIterator("test_files/fdx_2fe_2s.fasta", "protein"):
        print(seq.id)
        n_seq += 1
    print(f"Processed {n_seq} sequences")
