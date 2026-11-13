Multi-Criteria Analysis
=======================

The NBSOS uses a multi-criteria analysis (MCA) approach to identify optimal NBS that provide the overall highest benefits where spatial suitability of multiple NBS overlaps. The MCA is run through the :py:meth:`nbsos.urban_v1.Scan.run_multicriteria_analysis` method. This method:

.. _mca:

- Estimates the relative performance of NBS families, thus classifying the raw creation benefits and writing binned creation benefit rasters (cf. :ref:`here <benefits-classes>`);
- Applies a weighted-summation approach to identify optimal NBS based on the combined performance across multiple benefits, taking into account specific weighting scenarios;
- Identifies the locations of optimal NBS providing the highest overall benefits (by default: the 20% highest benefits).

The :py:meth:`nbsos.urban_v1.Scan.run_multicriteria_analysis` method allows users to run the MCA with different weighting scenarios, including:

- Equal weights, i.e., all benefits of interest receive the same weight. The equal weights scenario is run by setting ``assess_equal_weights_scenario=True``;
- Benefit-heavy weights, i.e., where one benefit of interest is given a higher weight compared to others, and where the remaining co-benefits receive an equal weight. The weight for that focal benefit is specified using the ``heavy_focal_weight`` parameter. Benefit-heavy scenarios are run by setting ``assess_benefit_heavy_scenario=True``. In so-doing, each benefit of interest is considered as the focal benefit in turn;
- Benefit-only weights, i.e., where only one benefit of interest is considered, thus receiving a weight of 100%, while other benefits are ignored. The benefit-only scenario is run by setting ``assess_benefit_only_scenario=True``.


.. code-block:: python

   # run multicriteria analysis
   my_scan.run_multicriteria_analysis(
       benefit_parameters=benefits,
       heavy_focal_weight=0.8,
       assess_equal_weights_scenario=True,
       assess_benefit_heavy_scenario=True,
       assess_benefit_only_scenario=True,
       custom_scenarios=None,
   )

In addition to these pre-defined scenarios, users can also define custom weighting scenarios by providing a dictionary to the ``custom_scenarios`` parameter. Each key in the dictionary represents the scenario name, and the corresponding value is another dictionary specifying the weights for each benefit of interest, as exemplified below.


.. code-block:: python

   custom_scenarios = {
       'flood-heat-heavy' : {
           Benefits.PluvialFlooding : 0.4,
           Benefits.UrbanHeat : 0.4,
           Benefits.Recreation : 0.1,
           Benefits.Biodiversity : 0.1,
       }
   }

   my_scan.run_multicriteria_analysis(
       benefit_parameters=benefits,
       heavy_focal_weight=0.8,
       assess_equal_weights_scenario=True,
       assess_benefit_heavy_scenario=True,
       assess_benefit_only_scenario=True,
       custom_scenarios=custom_scenarios,
   )


By default, the MCA identifies the optimal NBS providing the 20% highest benefits, overall, i.e., values above the 80th percentile. The  :py:meth:`nbsos.urban_v1.Scan.run_multicriteria_analysis` method allows users to customize this threshold, or use various thresholds simultaneously, by specifying desired percentiles as a list, e.g., ``percentiles=[70, 80, 90]`` provided as argument to the ``highest_benefit_percentiles`` parameter, as follows. In the example, the top 20%, 10%, 5% and 2.5% of benefits will be identified.

.. code-block:: python

   # Run MCA with custom highest benefits percentiles
   my_scan.run_multicriteria_analysis(
      benefit_parameters=benefits,
      heavy_focal_weight=0.85,
      assess_equal_weights_scenario=True,
      assess_benefit_heavy_scenario=True,
      assess_benefit_only_scenario=True,
      custom_scenarios= None,
      highest_benefits_percentiles=[80, 90, 95, 97.5]
   )

In the optimal NBS raster layers, pixel values denote associated NBS families, as follows.

.. list-table:: NBS Families in Optimal NBS Raster Layers
   :header-rows: 1
   :widths: 15 85

   * - Pixel Value
     - NBS Family

   * - 10
     - Bioretention Areas
   * - 11
     - Green Corridors
   * - 12
     - Inland Wetlands
   * - 13
     - Open Green Spaces
   * - 14
     - River/Stream Renaturation
   * - 15
     - Terraces and Slopes
   * - 16
     - Urban Farming
   * - 17
     - Urban Forests



.. admonition:: Technical note

    The relative performances of NBS families, i.e., binned creation benefits are stored in the `outputs/nbs_compare` folder within the data repository, organized by benefit type, following the naming convention: `<nbs family>_<benefit type>_binned.tif`. Optimal NBS, including optimal NBS providing the highest benefits, are stored in the `outputs/mca` folder, organized by assessed scenarios. For each scenario, a subfolder is created. Within that folder, the optimum MCA solutions in the form of optimal NBS and optimal NBS providing highest benefits are stored following the naming convention `<benefit type>_<scenario type>_optimum_nbs.tif` and `<benefit type>_<scenario type>_optimum_nbs_highest_<percentile>_benefits.tif`.