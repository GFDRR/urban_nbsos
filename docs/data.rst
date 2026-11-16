Data requirements
=================

Global datasets
---------------

The NBSOS uses a range of globally available geospatial and Earth Observation datasets, as summarized in :ref:`table-data`. This data is globally available free of cost, and is retrieved from Google Earth Engine (GEE).


.. _table-data:

.. list-table:: NBSOS Data Requirements
   :header-rows: 1
   :widths: 30 30 15 35

   * - Dataset
     - Reference
     - Resolution
     - Used in the NBSOS for

   * - CHIRPS Precipitation Pentad/Daily (Version 2.0)
     - Funk et al., 2015
     - 5566 m
     - Annual average precipitation; Days of heavy precipitation

   * - FABDEM
     - Hawker et al., 2022
     - 30 m
     - Digital elevation model; Slope estimation

   * - MERIT Hydro (HAND)
     - Yamazaki et al., 2019
     - 93 m
     - Height above nearest drainage estimation

   * - iSDAsoil USDA
     - Hengl et al., 2021
     - 30 m
     - Soil properties; Soil permeability; Hydrologic soil class (within Africa)

   * - OpenLandMap Soil Texture Class
     - Hengl, 2018
     - 250 m
     - Soil properties; Soil permeability; Hydrologic soil group (outside Africa)

   * - SoilGrids v2.0
     - Poggio et al., 2021
     - 250 m
     - Soil properties

   * - Harmonized Sentinel-2 MSI Level-2A
     - ESA/Copernicus
     - 10–60 m (band-specific)
     - Bare soil frequency

   * - JRC Monthly Water History v1.4
     - Pekel et al., 2016
     - 30 m
     - Water frequency

   * - JRC Global Surface Water v1.4
     - Pekel et al., 2016
     - 30 m
     - Surface water characterization

   * - ESA WorldCover v100
     - ESA / Zanaga et al., 2021
     - 10 m
     - Land cover

   * - GLAD LCLUC, LAND
     - Potapov et al., 2022
     - 30 m
     - Land cover

   * - Google Buildings V3 Polygons
     - Google Research / Sirko et al., 2021
     - Vector / 10 m rasterized
     - Building footprints

   * - World Settlement Footprint
     - DLR / Marconcini et al., 2020
     - 10 m
     - World Settlement Footprint 2019

   * - World Settlement Footprint Evolution
     - DLR / Marconcini et al., 2021
     - 30 m
     - World Settlement Footprint evolution

   * - WRI Tree Cover
     - WRI Datalab
     - 10 m
     - Identification of tree cover share

   * - HRSL Population
     - Facebook Connectivity Lab and CIESIN, 2016
     - 30 m
     - Population distribution

   * - MODIS Terra LST Emissivity Daily
     - NASA LP DAAC / USGS EROS
     - 1000 m
     - Estimation of annual average air temperature

   * - GWR Regression
     - Hooker et al., 2018
     - 5566 m
     - Statistical downscaling/regression

   * - Global Rainfall Erosivity 2017
     - EU JRC / Panagos et al., 2017
     - 1000 m
     - Rainfall erosivity estimation

   * - Soil erodibility estimation
     - Salhi et al. 2023
     - 250 m
     - Soil erodibility estimation

   * - NDVI
     - Various satellite sources
     - Variable
     - Vegetation productivity / greenness

   * - Road features
     - OpenStreetMap
     - Vector / 10 m rasterized
     - Road network and identification of major roads

   * - Water features
     - OpenStreetMap
     - Vector / 10 m rasterized
     - Preparation of a water body mask



Flood hazard data
-----------------

The assessment of (pluvial or fluvial) flood risk and NbS benefits for (pluvial) flood risk mitigation requires corresponding flood hazard data, i.e. flood extent and associated return period of the mapped flood events. This data needs to be obtained for the area of interest.
