Create a new SwissProt id from a given string.

Example
--------

.. code-block:: python

    from pybioshell.seq import seq_id

    swp_id = seq_id.swissprot("P12345")
    print(swp_id)
    assert swp_id.value == "P12345"
    assert swp_id.kind == "SwissProt"