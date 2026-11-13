Instantiate the scan for running the NBSOS
==========================================

Instantiation
-------------

First, instantiate the scan class. This class helps running the various elements of the NBSOS, e.g., data retrieval, suitability analysis, clustering, and benefit assessments. To instantiate the scan class, provide the previously defined scan context and the list of selected NbS families, as exemplified below.

.. code-block:: python

    # instantiate scan
    my_scan = Scan(
        scan_context=scan_context,
        nbs_list=nbs_list
    )


Assign the storage handler
--------------------------

Before running the scan, assign a storage handler to the scan instance. The storage handler needs to be instantiated prior to assignment to the scan instance.

.. code-block:: python

    # pass storage handler to scan. This will also intialize the handler
    my_scan.set_storage_options(storage_handler=my_storage_handler)
