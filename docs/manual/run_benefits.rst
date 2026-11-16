Estimate creation and protection benefits
=========================================

Based on identified NBS protection and creation opportunities, the NBSOS allows for the assessment of the associated benefits. These benefits are categorized into creation benefits and protection benefits. 

- Creation benefits are derived from the implementation of NBS and are assessed for each NBS family included in the analysis based on the determined NBS creation masks;
- Protection benefits are derived from existing urban green spaces and are assessed based on the NBS protection mask.


Creation benefits
-----------------

To assess creation benefits, the :py:meth:`nbsos.urban_v1.Scan.run_creation_benefits` method is used. It iterates over all benefits of interest, and assesses creation benefits for each NBS family of interest.

.. code-block:: python

   # run creation benefit models
   my_scan.run_creation_benefits(
      benefit_parameters=benefits,
   )

.. _benefits-classes:

Creation benefits are communicated along a qualitative scale, indicating low, medium, and high potential to produce desired benefits based on the assessed NBS creation opportunities, as follows.

.. list-table:: Classified NBS creation benefits
   :header-rows: 1
   :widths: 30 70

   * - Class Value
     - Estimated Creation Benefit

   * - 1
     - Low creation benefit
   * - 2
     - Medium creation benefit
   * - 3
     - High creation benefit

By default, bins are based on percentiles, where low benefits correspond to the 25% lowest values; medium benefits correspond to the 50% middle values; and high benefits correspond to the 25% highest values.

.. admonition:: Technical note

   Running the :py:meth:`nbsos.urban_v1.Scan.run_creation_benefits` method generates creation benefits for each NBS family of interest. Corresponding raw benefit data is stored across several folders within the data repository, depending on the type of benefit assessed. Various ancillary data may additionally be written by the benefit models. 
   
   This data is consolidated at a later step through the :py:meth:`nbsos.urban_v1.Scan.run_multicriteria_analysis` method, which produces the final binned creation benefit rasters (cf. :ref:`here <mca>`). The binned creation benefits are stored in the `outputs/nbs_compare` folder within the data repository, organized by benefit type, following the naming convention: `<nbs family>_<benefit type>_binned.tif`.

Protection benefits
-------------------

To assess protection benefits, the :py:meth:`nbsos.urban_v1.Scan.run_protection_benefits` method of the scan instance is used. This method iterates over benefits of interest, and assesses the protection benefits associated with existing urban green spaces based on the NBS protection mask for the following benefits:

- Pluvial flooding
- Heat stress
- Health and recreation
- Social cohesion
- Biodiversity
- Erosion

.. code-block:: python

   # run protection benefit models
   my_scan.run_protection_benefits(
      benefit_parameters=benefits
   )


Protection benefits are communicated along a qualitative scale, indicating low, medium, and high NBS protection benefits based on the assessed NBS protection opportunities, as follows. 

.. list-table:: Classified NBS protection benefits
   :header-rows: 1
   :widths: 30 70

   * - Class Value
     - Estimated Protection Benefit

   * - 1
     - Low protection benefit
   * - 2
     - Medium protection benefit
   * - 3
     - High protection benefit


By default, bins are based on percentiles, where low benefits correspond to the 25% lowest values; medium benefits correspond to the 50% middle values; and high benefits correspond to the 25% highest values.


.. admonition:: Technical note

   Binned protection benefits are stored in the `outputs/protection` folder within the data repository, organized by benefit type. The following naming convention applies: `nbs_protection_<benefit type>_benefits_binned.tif`.