Populate the data repository
============================

Running the :py:meth:`nbsos.urban_v1.Scan.run_preparation` method populates the local NBSOS data repository with relevant data, i.e., data processed in and retrieved from GEE, OSM, etc. Depending on the type of data handler used, this method in principle conducts the following steps:

- Copy the area of interest shapefile into the local data repository, and prepare other required ancillary data such as a hexagonal tiling layer as units of analysis;
- Prepare indicators in Google Earth Engine and download these indicators to the local data repository;
- Download data and indicators from other data sources, e.g., OpenStreetMap; 
- Process flood hazard data and store it into the local data repository;
- Determine additional indicators at the client-side. 

.. caution::
   
   Processing data in Google Earth Engine for download, particularly for larger areas of interest, may take a considerable amount of time.

.. code-block:: python

    # prepare indicators and download to local device
    my_scan.run_preparation()


Earth Observation (EO) Indicators
---------------------------------

Retrieved Earth Observation (EO) indicators are located in the `eo_indicators` folder within the data repository, and include the data listed in :ref:`eo-indicators`.

.. _eo-indicators:

.. list-table:: Earth Observation (EO) Indicators
   :header-rows: 1
   :widths: 35 65

   * - EO indicator filename
     - Description

   * - ``Annual_prec.tif``
     - Annual average precipitation

   * - ``Bs_freq.tif``
     - Bare soil frequency

   * - ``Building_mask.tif``
     - Building footprints, generated from the integration of Google Buildings and World Settlement Footprint 2019

   * - ``Building_min10_mask.tif``
     - Building mask, eroded by 1 pixel (10 m), used for suitability assessment

   * - ``Dem.tif``
     - Digital Elevation Model, scan-native resolution (10 m)

   * - ``Dem_30.tif``
     - Digital Elevation Model, 30 m resolution

   * - ``Google_buildings.tif``
     - Building footprints, sourced from Google (Sirko et al., 2021)

   * - ``Hand.tif``
     - Height above nearest drainage

   * - ``Heavy_prec.tif``
     - Days with heavy precipitation (≥ 10 mm)

   * - ``Hrsl_pop.tif``
     - HRSL population layer

   * - ``Hydrosoil.tif``
     - Hydrosoil group

   * - ``Land_cover.tif``
     - Land cover

   * - ``Lcc.tif``
     - Land capability class

   * - ``Lppi.tif``
     - Productivity performance indicator

   * - ``LST day/night tif files``
     - Monthly average daytime (nighttime) land surface temperature

   * - ``Major_roads_10m_mask.tif``
     - Rasterized major roads

   * - ``Road_mask.tif``
     - Rasterized roads (OSM road features)

   * - ``Rusle_k.tif``
     - Soil erodibility (RUSLE K factor)

   * - ``Rusle_ls.tif``
     - Slope-length/gradient (RUSLE LS factor)

   * - ``Rusle_r.tif``
     - Rainfall erosivity (RUSLE R factor)

   * - ``Slope.tif``
     - Slope (degrees)

   * - ``Soil_perm.tif``
     - Soil permeability class

   * - ``Soil_properties.tif``
     - Soil properties (multiband image with bands for texture, sand content, pH, silt content, clay content, bulk density, soil carbon, soil depth, and stone content)

   * - ``Water_50m_mask.tif``
     - Areas within 50 m of identified water bodies

   * - ``Water_mask.tif``
     - Water bodies mask, derived from the integration of land cover and rasterized OSM water features

   * - ``Wri_trees.tif``
     - Tree cover ratio, sourced from WRI

   * - ``Wsf_evo.tif``
     - World Settlement Footprint evolution (1980–2015)

   * - ``Wsf2019.tif``
     - World Settlement Footprint 2019


Ancillary data
--------------

The run method will also copy the AOI shapefile into the data repository, and create a hexagonal tiling of the AOI for subsequent analysis steps. For example, the hexagons are the basis for assessing NBS flood, heat, or biodiversity benefits.

.. admonition:: Technical note

    The AOI and hexagonal tiling layer for a scan are stored in the `ancillary` folder within the data repository. The AOI file follows the naming convention `<city name>_aoi.shp`, and the hexagonal tiling layer follows the naming convention `<city name>_hex.shp`. A rasterized version of that layer is stored in the same folder, following the naming convention `<city name>_aoi_hex.tif`.

Flood hazard data
-----------------

Flood hazard data, if relevant, is consolidated into the data repository within the `fathom` folder.
    