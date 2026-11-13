from . import loggerconfig

import os
import logging
import shutil
from pathlib import Path
import geopandas as gpd

# setup the module logger logger
logger = logging.getLogger('nbsos.scan')



def make_local_folders(local_root, country, city):
    """
    Create lcoal folder structure for NBSOS input and output data
    :param local_root: root folder for nbsos projects
    :param country: name of the country
    :param city: name of the city
    :return:
    """

    try:
        
        folders = [
            'eo_indicators',
            'fathom',
            'osm',
            'logs',
            'flood_benefit',
            'rec_benefit',
            'output',
            'maps',
            'cn',
            'ancillary',
            'suitability',
            'cluster_maps'
        ]
        # create all subfolders
        for folder in folders:
            folder = os.path.join(local_root, country, city, folder)
            os.makedirs(folder, exist_ok=True)
            logger.info(f"{folder} created locally")

        # create subfolders for mca and demand in output folder
        os.makedirs(os.path.join(local_root, country, city, 'output', 'mca'), exist_ok=True)
        os.makedirs(os.path.join(local_root, country, city, 'output', 'priority'), exist_ok=True)
        os.makedirs(os.path.join(local_root, country, city, 'output', 'protection'), exist_ok=True)

        # create a folder for masks to be consumed
        os.makedirs(os.path.join(local_root, country, city, 'masks'), exist_ok=True)


        # create subfolders in mca for various weighting schemas
        mca_folders = [
            'equal',
            'benefit_heavy',
            'benefit_only'
        ]
        for folder in mca_folders:
            folder = os.path.join(local_root, country, city, 'output', 'mca', folder)            
            os.makedirs(folder, exist_ok=True)
            # in addition, create a subfolder that will contain all weighted matrices of the nbs so that we have less chaos in the output folder
            folder2 = os.path.join(local_root, country, city, 'output', 'mca', folder, 'nbs_weighted')  
            os.makedirs(folder2, exist_ok=True)
            logger.info(f"{folder} created locally")


    except Exception as e:
        logger.error('Could not create local folder structure')
        raise e

def move_area(local_root, country, city, area):
    """
    move the aoi shapefile to the correct location for further processing
    local_root = local root for storing all nbsos projects
    country = country
    city =  city,
    area = locally stored copy of the aoi shapefile
    """
    area_folder = os.path.dirname(os.path.abspath(area))    
    stem_name = Path(area).stem

    # get all files of the shapefile dataset 
    # remove file extension and check against complete filename, otherwise, similar AOI files that contain name will also be copied, overwriting the intended AOI file 
    area_files = [f for f in os.listdir(area_folder) if Path(f).stem.lower() == (stem_name.lower())]    
    area_files = [f for f in area_files if not f.endswith('lock')]
    
    
    ancillary = os.path.join(local_root, country, city, 'ancillary')
    for area_file in area_files:

        # rename stem of area_file to assert city_aoi format 
        # this is needed in case aoi file had different format
        current_path = Path(area_file)
        new_filename = f'{city}_aoi{current_path.suffix}'

        area_path = os.path.join(area_folder, area_file)
        new_path = os.path.join(ancillary, new_filename)
        logger.info(f'copy {area_path} to {new_path}:')
        try:
            shutil.copy(area_path, new_path)
        except Exception as e:
            logger.error(f'could not copy {area_path} to {new_path}: {e}')

    # assert single-part polygons for the AOI shapefile
    aoi_shapefile = os.path.join(ancillary, f'{city}_aoi.shp')
    aoi_gdf = gpd.read_file(aoi_shapefile)
    print("Exploding AOI")
    aoi_gdf = aoi_gdf.explode(index_parts=False)    
    aoi_gdf.to_file(aoi_shapefile)


def upload_folder_tif(bucket, local_root, country, city, folder):
    loc_folder = os.path.join(local_root, country, city, folder)
    suitability_files = [f for f in os.listdir(loc_folder) if f.endswith('.tif')]
    locals = [os.path.join(loc_folder, f) for f in suitability_files]
    clouds = [f'nbsos/{country}/{city}/{folder}/{f}' for f in suitability_files]
    for local, cloud in zip(locals, clouds):
        logger.info(f'Exporting {local} to GCS: {cloud}')
        try:
            blob = bucket.blob(cloud)
            blob.upload_from_filename(local)
        except Exception as e:
            logger.error(f'Exporting {local} to GCS: {cloud} failed: {e}', exc_info=True)

def upload_folder(bucket, local_root, country, city, folder):
    loc_folder = os.path.join(local_root, country, city, folder)
    folder_files = os.listdir(loc_folder)
    locals = [os.path.join(loc_folder, f) for f in folder_files]
    clouds = [f'nbsos/{country}/{city}/{folder}/{f}' for f in folder_files]
    for local, cloud in zip(locals, clouds):
        logger.info(f'Exporting {local} to GCS: {cloud}')
        try:
            blob = bucket.blob(cloud)
            blob.upload_from_filename(local)
        except Exception as e:
            logger.error(f'Exporting {local} to GCS: {cloud} failed: {e}', exc_info=True)

def upload_folder_from_path(bucket, local_root, country, city, folder, exclude_folders=None, exclude_files_with_ending=['.aux.xml']):
    
    # construct local path, and iterate over files and directories in path
    loc_folder = os.path.join(local_root, country, city, folder)
    for root, dirs, files in os.walk(loc_folder):
        
        # root is the absolute path to the current folder.
        # see if we should exclude this folder from upload, if so, skip to next iteration
        current_folder = os.path.basename(root)
        print(f'Processing folder: {current_folder}')

        # exclude folders is a list of folder names. Exclude this folder
        # and any of its subfolders from upload
        # test if current folder is a subfolder of any excluded folder
        # get all parent folders of current folder
        parent_folders = Path(root).parents
        # get all parent folders until current city folder name is reached
        parents = []
        for f in parent_folders:
            if f.name == city:
                break
            parents.append(f.name)
        
        parents = list(reversed(parents))  # reverse to get correct order from city to current folder
        if exclude_folders is not None and (any(f in exclude_folders for f in parents) or current_folder in exclude_folders):
            print(f'Excluding folder {current_folder} from upload because it is a (sub-)folder of an excluded folder')
            continue

        # upload all files into corresponding cloud path, maintaining folder structure in cloud
        # exlude files with certain endings, e.g. .aux.xml files that are created when processing geotiffs in GEE, but are not needed for further processing and take up a lot of space
        filtered_files = [f for f in files if not any(f.endswith(ending) for ending in exclude_files_with_ending)]
        local_files = [os.path.join(root, file) for file in filtered_files]
        
        cloud_path = f'nbsos/{country}/{city}/'
        for parent in parents:
            cloud_path += f'{parent}/'
        cloud_path += f'{current_folder}'

        
        cloud_files = [f'{cloud_path}/{f}' for f in filtered_files]

        print(root)

        for local, cloud in zip(local_files, cloud_files):
            logger.info(f'Exporting {local} to GCS: {cloud}')
            try:
                #pass
                blob = bucket.blob(cloud)
                blob.upload_from_filename(local)
            except Exception as e:
                logger.error(f'Exporting {local} to GCS: {cloud} failed: {e}', exc_info=True)