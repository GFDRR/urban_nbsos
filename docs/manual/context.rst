
Configure the scan
==================

Define the scan context
-----------------------

The :py:class:`nbsos.urban_v1.ScanContext` class is used to specify environmental variables relevant for the scan. This includes the local path to the NBSOS data repository, and information to structure the data repository (country name and name of the city of the scan, and path to the prepared AOI shapefile). Other relevant parameters for running the scan include whether the scan location is in Africa, which determines which datasets are used for the scan in question.

.. code-block:: python

    # define the context for the scan
    scan_context: ScanContext = ScanContext()
    scan_context.local_root = 'path/to/nbsos/data'
    scan_context.country = '<country name>'
    scan_context.city = '<city name>'
    scan_context.location_in_africa = <True or False>
    scan_context.local_area_of_interest_file = '/path/to/aoi_shapefile'


Define flood hazard data location
---------------------------------

If flood benefits are to be assessed, flood hazard data must be provided for the scan. Typically, lower (e.g., 1 in 20 years, 1 in 50 years), intermediate (e.g., 1 in 100 years, 1 in 200 years), and higher (e.g., 1 in 500 years, 1 in 1500 years) return periods are included in a scan, as available. 

The location of these data must be specified, and a mapping of filenames and associated return periods be provided. This mapping ensures that the scan can correctly associate each flood hazard file with its corresponding return period. This mapping is defined as a dictionary where the keys are the filenames of the flood hazard data files located at the specified path, and the values are the corresponding return periods. 


.. code-block:: python

    # define the path to flood hazard data
    pluvial_flood_hazard_data_path = '/path/to/flood/hazard/data'

    # map files to return periods
    pluvial_flood_hazard_data_mapping = {
        'filename.tif' : 20,
        'filename.tif' : 50,
        'filename.tif' : 100,
        'filename.tif' : 200,
        'filename.tif' : 500,
        'filename.tif' : 1500
    }
