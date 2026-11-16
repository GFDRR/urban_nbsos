Identify priority areas
=======================

To identify exposure of population to hazards, and determine priority areas for NBS implementation to reduce hazard exposure based on demand as a function of exposure, the :py:meth:`nbsos.urban_v1.Scan.run_priority` method is used. This method determines expsoure and priority areas based on the hexagonal tiling of the study area and  requires that creation benefit models were run previously.


.. list-table:: Exposure and priority area estimation
   :header-rows: 1
   :widths: 20 40 40

   * - Benefit
     - Exposure
     - Priority area

   * - Pluvial flooding
     - Population exposed to flooding, as a function of population and annual flood probability estimated at the neighbourhood level, i.e., the hexagonal tiling.
     - Using results about exposure to flood hazard distribution and a simple hydrological model to generate a series of watersheds, at the neighbourhood level, to account for downstream effects of runoff reduction.
   * - Heat stress
     - Population potentially exposed to heat stress, determined as a function of population and the estimated annual average air temperature derived from applying a geographically-weighted and climate space-weighted regression model by Hooker, Duveiller, and Cescatti (2018).
     - Higher priority corresponds to higher exposure (in-situ).
   * - Health and recreation
     - Population without access to min. 1.0 ha of green space within 300 m.
     - Higher priority corresponds to higher exposure (in-situ).
   * - Social cohesion
     - Population without access to min. 0.5 ha of green space within 300 m.
     - Higher priority corresponds to higher exposure (in-situ).
   * - Erosion
     - Population in erosion-prone areas.
     - Higher priority corresponds to higher exposure (in-situ).


For the flood hazard, areas of high exposure to flooding and priority areas to address flood exposure may not spatially coincide. This is because it is often more effective to implement NBS for flood mitigation at upstream locations of flood-prone areas rather than directly at the locations experiencing flooding. For the other benefits, areas of high exposure are tpically also areas where NBS should be implemented to address the demand for reducing hazard exposure. This is indicated in the above table through the term 'in-situ'. 

By default, exposure and priority are both communicated along a qualitative scale, denoting very low, low, medium, high, and very high levels of exposure or priority, as follows. A five-breaks Jenks natural breaks classification schema is used to categorize exposure and priority. 

.. list-table:: Classified exposure and priority
   :header-rows: 1
   :widths: 20 40

   * - Class Value
     - Description

   * - 1
     - Very low exposure or priority
   * - 2
     - Low exposure or priority
   * - 3
     - Medium exposure or priority
   * - 4
     - High exposure or priority
   * - 5
     - Very high exposure or priority



.. code-block:: python

    # identify priority areas
    my_scan.run_priority(
        benefit_parameters=benefits
    )

.. admonition:: Technical note

    Binned rasters of prioritty scores are stored in the `output/priority` folder within the data repository, organized by hazard type, following the naming convention `<hazard type>_priority_binned.tif`. For hazards where exposure and priority areas are not considered in-situ, binned exposure is stored in addition, following the naming convention `<hazard type>_exposure_binned.tif`.
