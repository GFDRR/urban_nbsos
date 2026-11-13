Benefits
========

All quantitative benefits supported by the NBSOS are defined in the :py:class:`nbsos.urban_v1.Benefits` enumeration. However, in a given scan, only a subset of these benefits may be of interest, for example, as discussed in the inception. Therefore, the benefits to be assessed need to be defined per scan. 


Benefits of interest are defined as a list of :py:class:`nbsos.urban_v1.ScanBenefit` instances. The :py:class:`nbsos.urban_v1.ScanBenefit` class allows defining the benefit to be assessed, as well as specifying parameters for its assessment if default parameters shall be overridden for a given scan. To do so, the class provides a parameter `params` as a means to specify parameters for the assessment of each benefit. This parameter takes a dictionary of parameter names and values, as exemplified below.

.. note::
    
    To assess pluvial flood benefits (:py:class:`nbsos.urban_v1.Benefits.PluvialFlooding`), the return periods to include in the assessment need to be specified for each scan. This is done using the `return_periods` parameter within the `params` dictionary of the :py:class:`nbsos.urban_v1.ScanBenefit` instance, passing a list of return periods of interest as argument, as exemplified below.
    

.. code-block:: python

    # define benefits of interest
    benefits = [
        ScanBenefit(Benefits.PluvialFlooding, params={'return_periods' : [20, 50, 100, 200, 500, 1500]}),
        ScanBenefit(Benefits.UrbanHeat, params={'temperature_resampling_method' : Resampling.cubic}),
        ScanBenefit(Benefits.Recreation),
        ScanBenefit(Benefits.Biodiversity),
        ScanBenefit(Benefits.AvoidedErosion),
    ]

.. admonition:: Parameters for benefit assessment

    Sensible default parameters are specified for all benefit modules by default. The parameters that may be provided as key-value pairs in params for each benefit are documented here:

    * :py:meth:`nbsos.urban_v1.Scan.get_pluvial_flood_creation_benefits`
    * :py:meth:`nbsos.urban_v1.Scan.get_heat_creation_benefits`
    * :py:meth:`nbsos.urban_v1.Scan.get_recreation_creation_benefits`
    * :py:meth:`nbsos.urban_v1.Scan.get_flood_demand`
    * :py:meth:`nbsos.urban_v1.Scan.get_biodiversity_creation_benefits`
    * :py:meth:`nbsos.urban_v1.Scan.get_erosion_creation_benefits`
    * :py:meth:`nbsos.urban_v1.Scan.get_job_creation_creation_benefits`



