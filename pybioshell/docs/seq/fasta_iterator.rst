An iterator over sequences in a FASTA file.

Yields `Sequence` objects one by one; reading only a fragment
of the input file at a time.

Example
--------

.. code-block:: python

    from bioshell.seq import FastaIterator

    n_seq = 0
    for seq in FastaIterator("test_files/fdx_2fe_2s.fasta"):
        print(seq.id)
        n_seq += 1
    print(f"Processed {n_seq} sequences")
