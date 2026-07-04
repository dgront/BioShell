Create a new PDB-id from a given string.

Example
--------

.. code-block:: python

    from pybioshell.seq import seq_id

    pdb_id = seq_id.pdb("1a2b")
    print(pdb_id)
    assert pdb_id.value == "1a2b"