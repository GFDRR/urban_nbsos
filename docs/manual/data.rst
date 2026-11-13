
Data handlers
=============

Data handlers (or storage handlers) support the processing and transfer/download of data that is used in the NBSOS. 

The data handler is responsible for

* preparing and managing the processing of geospatial and Earth Observation data in Google Earth Engine (GEE);
* transferring data between GEE and the local device (the local data repository);
* retrieving data from the OpenStreetMap (OSM);
* copying flood hazard data into the local data repository;
* processing downloaded data locally to determine a range of additional indicators.



Available data handlers
-----------------------


At this time, the NBSOS supports two types of data handler modules: 

* A *GcsHandler* that supports the transfer of data between Google Earth Engine (GEE) and the local device via Google Cloud Storage (GCS); and
* A *DirectDownloadHandler* that supports the direct download of data from GEE to the local device.

.. warning::
    The prepared path to the flood hazard data folder and the mapping of specific files to their associated return periods must be provided to the either storage handler at the time of instantiation. 

Using the GcsHandler
^^^^^^^^^^^^^^^^^^^^

The GcsHandler transfers data between GEE and the local device via GCS. I.e., data for the AOI is first exported from GEE into a user-specified GCS bucket, in which data is organized by country and city name, similar to the local data repository. Once all indicators are exported to GCS, the data is then downloaded from GCS to the local device. 

The use of the GcsHandler (:py:class:`nbsos.urban_v1.GcsHandler`) is demonstrated in the following example. Note that the GcsHandler requires a path to a credentials file and the name of the GCS bucket to be specified after instantiation.

.. code-block:: python

    # instantiate GCS storage handler
    my_storage_handler = GcsHandler(
        scan_context=scan_context,
        flood_hazard_data_folder=pluvial_flood_hazard_data_path,
        file_to_probability_mapping=pluvial_flood_hazard_data_mapping
    )
    # add required information, depending on storage handler chosen
    storage_handler.credentials = 'path/to/credentials/file'
    storage_handler.bucket_name = 'gcs-bucket-name'


Using the DirectDownloadHandler
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The DirectDownloadHandler transfers data directly from GEE to the local device, bypassing the need for an intermediate GCS bucket, as well as bypassing an intermediate export to Google Drive. The use of the GcsHandler (:py:class:`nbsos.urban_v1.DirectDownloadHandler`) is demonstrated in the following example. Note that the DirectDownloadHandler requires the name of the GEE project to be specified after instantiation. This project is used for processing the data in GEE. 



.. code-block:: python
        
    # instantiate storage handler for direct download
    my_storage_handler = DirectDownloadHandler(
        scan_context=scan_context,
        flood_hazard_data_folder=pluvial_flood_hazard_data_path,
        file_to_probability_mapping=pluvial_flood_hazard_data_mapping
    )
    # set name of project used to process data in GEE
    my_storage_handler.project_name = 'project-name'


.. warning::
    The DirectDownloadHandler implements an inofficial export mechanism for transferring data from GEE to the local device. Users should be aware that this method may not be officially supported and could be subject to changes. This method limits the size of data that can be exported in a single operation, and thus affects the potential extent of the AOI. Other limitations like quota restrictions and processing time in GEE may additionally apply.



Assigning a storage handler to the scan instance
------------------------------------------------

Assign the instantiated storage handler to the scan instance as exemplified below.

.. code-block:: python

    # pass storage handler to scan. This will also intialize the handler
    my_scan.set_storage_options(storage_handler=my_storage_handler)

