Run the suitability model
=========================

To run the suitability model, use the :py:meth:`nbsos.urban_v1.Scan.run_suitability` method of the scan instance, as exemplified below. This will create suitability masks for each of the NBS families included in the scan.

.. code-block:: python

    # run the suitability model
    my_scan.run_suitability()


.. admonition:: Technical note

    Within the data repository, the suitability masks are stored in the `suitability` folder, following the naming convention `<nbs_family>_suitability.tif`, where `<nbs_family>` corresponds to the specific NBS family for which the suitability mask was generated (see also :py:class:`nbsos.urban_v1.NbsTypes`).


.. caution::

    The suitability analysis considers various environmental and climatic factors to determine suitable areas for each of the NBS families. This includes factors such as slope, land cover, or soil type. From the factors considered, particularly annual average precipitation may play a significant role in determining suitability. In water-scarce contexts, lack of precipitation or sufficient water supply may hinder NBS implementation; in such cases, the NBSOS may fail to identify suitable areas for certain NBS families, such as urban forests. 

    To force the NBSOS to identify suitable areas regardless of annual average precipitation, a parameter *ignore_annual_precipitation* is provided for the run_suitability method. You can set this parameter to `True` to bypass the consideration of annual average precipitation in the suitability analysis. By default, this argument is set to `False`. Note however that in the delivery of results, the potential need for irrigation should be taken into account.
    
