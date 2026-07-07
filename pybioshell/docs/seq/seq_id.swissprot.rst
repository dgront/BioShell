Create a new SwissProt id from a given string.

`SwissProt` id actually is the `UniProtKB` id with the `sp` prefix to indicate that it is a SwissProt entry - the curated part of the UniProt database.

Example
--------

.. code-block:: python

    from pybioshell.seq import seq_id

    swp_id = seq_id.swissprot("P12345")
    assert str(swp_id) == "sp|P12345"
    assert swp_id.value == "P12345"
    assert swp_id.kind == "SwissProt"