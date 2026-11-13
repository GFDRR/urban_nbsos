Import the NBSOS
================

To start, import relevant classes from the NBSOS package. Typically, these include:

* The :py:class:`nbsos.urban_v1.ScanContext` class, that implements variables for defining the spatial context of a scan.
* The `NbsTypes` enumeration, to define the NBS families of interest to be considered in the scan.
* The `Benefits` enumeration, that lists the benefits available for (quantitative) assessment in the NBSOS.
* The :py:class:`nbsos.urban_v1.ScanBenefit` class, that is used to define benefits of interest to be considered in the scan.
* A storage/data handler, to handle transfer of data. Here, the :py:class:`nbsos.urban_v1.DirectDownloadHandler` class is used.
* The :py:class:`nbsos.urban_v1.Scan` class, which is the main class in the package used to run a scan.

.. code-block::  python

    # import relevant classes
    from nbsos.urban_v1 import NbsTypes, Benefits, ScanBenefit
    from nbsos.urban_v1 import ScanContext
    from nbsos.urban_v1 import DirectDownloadHandler
    from nbsos.urban_v1 import Scan
