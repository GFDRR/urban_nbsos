NBS families
============

The NBS families to be considered in the scan must be provided to the scan as a list of NbsTypes, as exemplified below. The :py:class:`nbsos.urban_v1.NbsTypes` enumeration defines all types of NBS families supported by the NBSOS.

.. code-block:: python

    # define NBS families to consider in the scan
    nbs_list = [
        NbsTypes.OpenGreenSpace,
        NbsTypes.GreenCorridors,
        NbsTypes.UrbanForest,
        NbsTypes.InlandWetlands,
        NbsTypes.BioretentionAreas
    ]
