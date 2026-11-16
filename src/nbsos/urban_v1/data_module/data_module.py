from enum import Enum
from typing import Protocol, List, Dict
from abc import abstractmethod
import logging
from tqdm import tqdm
from tqdm.auto import tqdm as tqauto
import time
import ee
import geopandas as gpd
import rioxarray
import os
import requests

from .. import google_auth
from .. import get_local_indicators
from .. import get_external_data
from .. import get_indicators
from .. import file_structuring
from ..scan_context import ScanContext
from .. import gis_tools

from ..heat_module.heat_module import HeatModule
from ..biodiversity_module.biodiversity_module import BiodiversityModule
from ..erosion_module import ErosionModule
from ..jobs_module import JobCreationModule


logger = logging.getLogger('nbsos.scan')


class DataModule(Protocol):
    
    local_root: str = None
    country: str = None
    city: str = None
    location_in_africa: bool = None
    area_of_interest: str = None
    scan_context: ScanContext = None
    flood_hazard_data_folder: str = None
    file_to_probability_mapping: Dict[str, int] = None

    def __init__(self, 
                 scan_context: ScanContext, 
                 flood_hazard_data_folder: str, 
                 file_to_probability_mapping: Dict[str, int],
                 overpass_url: str = "https://overpass-api.de/api",
                 overpass_timeout: int = 180
                 ):
        
        logger.info('Initializing storage handler')
        self.local_root = scan_context.local_root
        self.country = scan_context.country
        self.city = scan_context.city
        self.area = scan_context.local_area_of_interest_file
        self.location_in_africa = scan_context.location_in_africa

        self.flood_hazard_data_folder = flood_hazard_data_folder 
        self.file_to_probability_mapping = file_to_probability_mapping
        self.overpass_url = overpass_url
        self.overpass_timeout = overpass_timeout

    
    def make_local_folders(self):
        logger.info(f"Creating local folder structure in {self.local_root}")
        file_structuring.make_local_folders(self.local_root, self.country, self.city)
        file_structuring.move_area(self.local_root, self.country, self.city, self.area)

    def get_temperature_regression_coefficients(self, local_folder: str):

        # create directory
        os.makedirs(local_folder, exist_ok=True)

        # download data from paper DOI
        files = {
            "Coeffs_GWR_b0.nc" : "https://figshare.com/ndownloader/files/12983822",
            "Coeffs_GWR_b1.nc" : "https://figshare.com/ndownloader/files/12984128",
            "Coeffs_GWR_b2.nc" : "https://figshare.com/ndownloader/files/12984044",
            "Coeffs_CSWR_b0.nc" : "https://figshare.com/ndownloader/files/12983702",
            "Coeffs_CSWR_b1.nc" : "https://figshare.com/ndownloader/files/12983705",
            "Coeffs_CSWR_b2.nc" : "https://figshare.com/ndownloader/files/12984131"
        }

        for filename, url in tqauto(files.items()):
            response = requests.get(url) 
            with open(os.path.join(local_folder, filename), 'wb') as f:
                f.write(response.content)
        

    def get_external_data(self):
        
        # get external data from OSM
        logger.info("Download OSM data")
        get_external_data.get_osm(self.area, self.local_root, self.country, self.city, self.overpass_url, self.overpass_timeout)
        # make tessellation for AOI
        logger.info("Getting H3 hexagons for area of interest")
        get_external_data.get_hexagons(self.local_root, self.country, self.city)
        

    def get_client_indicators(self):
        
        logger.info("Rasterizing Hexagon tiles")
        
        # rasterize hex layer
        gis_tools.rasterize(
            os.path.join(self.local_root, self.country, self.city, 'ancillary', f'{self.city}_hex.shp'), 
            os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'land_cover.tif'), 
            os.path.join(self.local_root, self.country, self.city, 'ancillary', f'{self.city}_aoi_hex.tif'),
            'neigh_id'
        )

        # extract flood hazard data
        logger.info("Extracting local flood hazard data")
        # extract flood hazard data only if a folder has been provided. If this is not the case, skip this step.
        if self.flood_hazard_data_folder is not None:
            get_external_data.extract_flood_hazard_data(self.flood_hazard_data_folder, self.country, self.city, self.local_root, self.file_to_probability_mapping)
        else:
            logger.info("No flood hazard data folder provided. Skipping flood hazard data extraction.")

        logger.info('Calculate additional indicators using client-side Python')
        locals = get_local_indicators.local_indicators(self.local_root, self.country, self.city)
        locals.get_local_indicators()
       
    @abstractmethod
    def initalize_storage(self): raise NotImplementedError

    @abstractmethod
    def gee_to_storage(self, 
        download_lppi_only: bool = False, 
        include_basic_indicators: bool = True, 
        include_extended_indicators: bool = True, 
        include_rusle: bool = True, 
        include_LST: bool = True): raise NotImplementedError

    @abstractmethod
    def storage_to_local(self): raise NotImplementedError

    @abstractmethod
    def client_compute_indicators(self): raise NotImplementedError

    @abstractmethod
    def archive_to_storage(self): raise NotImplementedError





class GcsHandler(DataModule):
    
    _bucket_name = None
    _credentials = None
    _bucket = None


    def __init__(self, 
                 scan_context: ScanContext, 
                 flood_hazard_data_folder: str = None, 
                 file_to_probability_mapping: Dict[str, int] = None,
                 overpass_url: str = "https://overpass-api.de/api",
                 overpass_timeout: int = 180
                 ):
        
        super().__init__(scan_context, flood_hazard_data_folder, file_to_probability_mapping, overpass_url, overpass_timeout)
        
    @property
    def bucket_name(self):
        return self._bucket_name
    @bucket_name.setter
    def bucket_name(self, value):
        self._bucket_name = value

    @property
    def credentials(self):
        return self._credentials    
    @credentials.setter
    def credentials(self, value):
        self._credentials = value

    @property
    def bucket(self):
        return self._bucket
    @bucket.setter
    def bucket(self, value):
        self._bucket = value

    def initalize_storage(self):
        """Initializes the storage handler for Google Cloud Storage (GCS). The method sets up the local folder structure,
           and starts the Google Cloud session.
        """

        # initialize folder structure upon starting a scan
        super().make_local_folders()

        self.start_google_session()
        logger.info(f'GCS Data Module initialized using bucket {self.bucket_name}')
        

    def gee_to_storage(self, 
                       download_lppi_only: bool = False, 
                       include_basic_indicators: bool = True, 
                        include_extended_indicators: bool = True, 
                        include_rusle: bool = True, 
                        include_LST: bool = True):
        logger.info("Running GEE to Storage")
        self.gee_to_gcs(download_lppi_only=download_lppi_only,
                        include_basic_indicators=include_basic_indicators,
                        include_extended_indicators=include_extended_indicators,
                        include_rusle=include_rusle,
                        include_LST=include_LST)


    def storage_to_local(self):
        logger.info("Running Storage to Local")
        self.get_local_indicators()
        super().get_external_data()

    def client_compute_indicators(self):
        super().get_client_indicators()

    def archive_to_storage(self):
        self.upload_to_gcs()


    def start_google_session(self):
        logger.info("Connecting to Google Cloud Storage and EE")
        S = google_auth.Authenticate(self._credentials)
        S.start_google_cloud_session(self.bucket_name)
        S.initialize_ee()
        self._bucket = S.bucket

    def gee_to_gcs(self, 
                download_lppi_only: bool = False, 
                include_basic_indicators: bool = True, 
                include_extended_indicators: bool = True, 
                include_rusle: bool = True, 
                include_LST: bool = True) -> List[any]:

        crs = gis_tools.get_area_crs(self.area)
        crs_str = f'EPSG:{crs}'

        tasks = self.get_gee_indicators(crs_str, 
                                        download_lppi_only=download_lppi_only,
                                        include_basic_indicators=include_basic_indicators,
                                        include_extended_indicators=include_extended_indicators,
                                        include_rusle=include_rusle,
                                        include_LST=include_LST)
        
        total_task_count = len(tasks)
        logger.info(f'Total number of tasks: {total_task_count}')

        known_tasks = []

        with tqauto(total=total_task_count, position=0, leave=True) as pbar:

            while True:

                # check remaining tasks status
                for t in tasks:
                    current_task_id = t.id
                    status = t.status()

                    if status['state'] == 'COMPLETED':
                        if current_task_id not in known_tasks:
                            known_tasks.append(current_task_id)
                            pbar.update(1)
                            tqauto.write(f'{status["description"]} export is {status["state"]}')

                    elif status['state'] == 'FAILED':
                        if current_task_id not in known_tasks:
                            known_tasks.append(current_task_id)
                            pbar.update(1)
                            tqauto.write(f'{status["description"]} export is {status["state"]}')
                    
                    else:
                        continue

                if len(known_tasks) == total_task_count:
                    break

                time.sleep(10)

        logger.info("All tasks completed.")
        


    def get_gee_indicators(self,             
            crs, 
            download_lppi_only: bool = False,
            include_basic_indicators: bool = True, 
            include_extended_indicators: bool = True, 
            include_rusle: bool = True, 
            include_LST: bool = True
        ) -> List[any]:

        logger.info('Getting all necessary eo-indicators from GEE')
        i = get_indicators.Indicators(self.country, self.city, self.local_root, self.bucket, africa=self.location_in_africa)
        
        if not download_lppi_only:
            images = i.download_all(
                include_basic_indicators=include_basic_indicators,
                include_extended_indicators=include_extended_indicators,
                include_rusle=include_rusle,
                include_LST=include_LST
            )
        elif download_lppi_only:
            images = i.download_lppi()

        logger.info(f'Exporting all GEE eo-indicators to {self.bucket_name} in GCS')
        logger.info(f'This could take a while. Check the status of downloads.')
        exporting = get_indicators.ExportIndicators(self.bucket_name, self.country, self.city, i.aoi)
        tasks = exporting.export_all(images, crs)
        return tasks

    def get_local_indicators(self):
        logger.info('Download local copy of GEE eo-indicators')
        gcs_loc = get_local_indicators.GcsToLocal(self.bucket, self.country, self.city, self.local_root)
        gcs_loc.get_eo_indicators()        

    def upload_to_gcs(self):
        logger.info(f'upload data to GCS')
        file_structuring.upload_folder_tif(self.bucket, self.local_root, self.country, self.city, 'eo_indicators')
        file_structuring.upload_folder_tif(self.bucket, self.local_root, self.country, self.city, 'suitability')
        file_structuring.upload_folder_tif(self.bucket, self.local_root, self.country, self.city, 'cn')
        file_structuring.upload_folder(self.bucket, self.local_root, self.country, self.city, 'rec_benefit')
        file_structuring.upload_folder(self.bucket, self.local_root, self.country, self.city, 'osm')
        file_structuring.upload_folder(self.bucket, self.local_root, self.country, self.city, 'ancillary')

        # archive benefit modules
        file_structuring.upload_folder_tif(self.bucket, self.local_root, self.country, self.city, 'flood_benefit')
        file_structuring.upload_folder_from_path(self.bucket, self.local_root, self.country, self.city, ErosionModule.get_module_pathname())
        file_structuring.upload_folder_from_path(self.bucket, self.local_root, self.country, self.city, HeatModule.get_module_pathname())        
        file_structuring.upload_folder_from_path(self.bucket, self.local_root, self.country, self.city, BiodiversityModule.get_module_pathname())
        file_structuring.upload_folder_from_path(self.bucket, self.local_root, self.country, self.city, JobCreationModule.get_module_pathname())

        file_structuring.upload_folder_from_path(self.bucket, self.local_root, self.country, self.city, 'output', exclude_folders=['nbs_compare', 'mca'])

class DirectDownloadHandler(DataModule):
    
    project: str = None

    def __init__(self, 
                 scan_context: ScanContext, 
                 flood_hazard_data_folder: str = None, 
                 file_to_probability_mapping: Dict[str, int] = None,
                 overpass_url: str = "https://overpass-api.de/api",
                 overpass_timeout: int = 180
                 ):
        
        super().__init__(scan_context, flood_hazard_data_folder, file_to_probability_mapping, overpass_url, overpass_timeout)

    # implement project name as property to be more consistent with GcsDataModule
    @property
    def project_name(self) -> str:
        return self.project
    @project_name.setter
    def project_name(self, value: str):
        self.project = value

    def initalize_storage(self):
        """Initializes the storage handler for direct download from Google Earth Engine. The method sets up the local folder structure,
           and connects to the Google Earth Engine project.
        """

        # prepare the environment
        # initialize folder structure upon starting a scan
        super().make_local_folders()

        self.start_google_session()
        logger.info(f'Storage handler initialized using project name {self.project_name}')

    def client_compute_indicators(self):
        super().get_client_indicators()

    def start_google_session(self):
        logger.info("Connecting to  GEE")

        gg  = google_auth.AuthenticateIndividual(self.project)
        gg.initialize_ee_indiviudal(self.project)
        pass

    def gee_to_storage(self, 
            overwrite_files: bool = False, 
            download_lppi_only: bool = False,
            include_basic_indicators: bool = True, 
            include_extended_indicators: bool = True, 
            include_rusle: bool = True, 
            include_LST: bool = True
        ):
        logger.info("Running Storage to Local")        
        i = get_indicators.Indicators(self.country, self.city, self.local_root, None,
                                      africa=self.location_in_africa, project=self.project_name)
        
        if not download_lppi_only:        
            images = i.download_all(
                include_basic_indicators=include_basic_indicators,
                include_extended_indicators=include_extended_indicators,
                include_rusle=include_rusle,
                include_LST=include_LST
            )
        elif download_lppi_only:
            images = i.download_lppi()
        
        crs = gis_tools.get_area_crs(self.area)
        logger.info(f'Exporting all GEE eo-indicators to {self.local_root}')

        for image in tqauto(images):
            i.export_ee(image, crs, over_write=overwrite_files)

    def storage_to_local(self):
        super().get_external_data()
    
    def archive_to_storage(self):
        logger.info('There is currently no archiving implemented for this storage handler.')