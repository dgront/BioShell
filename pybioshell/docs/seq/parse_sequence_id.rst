Extract sequence identifiers from a free-text description string.

This function scans the input for standard database identifiers such as those from
PDB, SwissProt/TrEMBL, UniRef, RefSeq, GenBank/EMBL, Ensembl, and NCBI GI.
It returns all matches found, each represented as a [`SeqID`] object

Example
--------

.. code-block:: python

    from bioshell.seq import parse_sequence_id

    descr = "sp|P0A3D1|PETF_ECOLI 2Fe-2S ferredoxin [Escherichia coli] OX=562"
    id_list = parse_sequence_id(descr)
    assert id_list[0].kind == "SwissProt"
    assert id_list[-1].kind == "TaxId"