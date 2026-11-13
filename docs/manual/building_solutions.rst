Building solutions
==================

The building solutions NBS familiy is run differently to other NBS families. It focuses on the integration of nature-based solutions within the built environment, thus referring to the potential of green roofs,  green walls, etc. to reduce runoff for pluvial flood benefits and/or reduce urban heat stress. 

This NBS family has no specific suitability rules. It builds on the flood and/or heat priority areas, and is thus run separately from other NBS families through the :py:meth:`nbsos.urban_v1.Scan.run_building_solutions` method.


.. code-block:: python

    # run building solutions NBS family at the end of the scan, once priority areas are identified
    my_scan.run_building_solutions()


The :py:meth:`nbsos.urban_v1.Scan.run_building_solutions` determines a building priority score for flood and/or heat, depending on the defined benefits of interest. Other potential benefits of building solutions are currently not considered in the scoring. 

The building priority indicates where building solutions could be most effective to help mitigate pluvial flooding and/or heat stress and is determined at the neighborhood (hexagonal tiling) level by combining the following criteria:

- Areas with a high density of buildings;
- Areas with limited creation opportunities to implement other NBS families;
- Areas with high priority for applying NBS to reduce flooding or heat.

This building solutions priority score is the higher ...

- the higher the share (density) of buildings;
- the lower the share of protection opportunities;
- the lower the share of NBS creation opportunities;
- the higher the flood or heat priority.  

By default, building priority is communicated along a qualitative scale, denoting very low, low, medium, high, and very high levels of priority, as follows. A five-breaks Jenks natural breaks classification schema is used to categorize building priority. 

.. list-table:: Classified building priority score
   :header-rows: 1
   :widths: 20 40

   * - Class Value
     - Description

   * - 1
     - Very low priority
   * - 2
     - Low priority
   * - 3
     - Medium priority
   * - 4
     - High priority
   * - 5
     - Very high priority



.. admonition:: Technical note

    Binned rasters of building priority scores are stored in the `output/building_solutions` folder within the data repository, organized by hazard type, following the naming convention `building_solutions_<hazard type>_binned.tif`. A shapefile is additionally written to this folder to support the inspection of the score criteria at the neighbourhood level.  
