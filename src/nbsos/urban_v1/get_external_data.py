from . import loggerconfig
from . import gis_tools

import osmnx
import osmnx.settings
from shapely.geometry import MultiPolygon, shape, Polygon
import fiona
import warnings
import os
import logging
import h3pandas
import geopandas as gpd
from shapely.ops import unary_union
import rioxarray as rds
from typing import Dict
import sys

logger = logging.getLogger('nbsos.scan')

def get_osm_part(Multi, tags, name, output_folder):

    # get osm gdf
    logger.info(f'Downloading {name} OSM data')
    
    try:
        osm = osmnx.features_from_polygon(Multi, tags)
    
    except Exception as e:
        logger.error(f'Could not download OSM {name} OSM data')
        raise Exception(f'Could not download OSM {name} OSM data: {e}')
    
    # download successful. Proceed with processing
    df_colnames = list(osm.columns)
    
    selected_cols = []
    for col in ['name', 'geometry', 'highway']:
        if col in df_colnames:
            selected_cols.append(col)
    osm = osm[selected_cols]


    # save output
    for geomtype in osm.geom_type.unique(): # save shapefile for each feature type
        file_name = f'osm_{name}_{geomtype}.shp'
        full_path = os.path.join(output_folder, file_name)
        osm[osm.geom_type == geomtype].to_file(full_path)



def get_osm(area: str, 
            root_folder: str, 
            country: str, 
            city: str, 
            overpass_url: str = "https://overpass-api.de/api",
            overpass_timeout: int = 180,
            ):
    
    """Download OSM data.

    :param area: Path to the shapefile of the area of interest (AOI)
    :type area: str
    :param root_folder: Root folder for local storage
    :type root_folder: str
    :param country: Country name
    :type country: str
    :param city: City name
    :type city: str
    :param overpass_url: URL of the Overpass API endpoint, defaults to "https://overpass-api.de/api"
    :type overpass_url: str, optional
    :param overpass_timeout: Timeout for Overpass API requests in seconds, defaults to 180
    :type overpass_timeout: int, optional

    """

    warnings.filterwarnings(action='ignore', category=DeprecationWarning)
    warnings.filterwarnings(action='ignore', category=FutureWarning)

    output_folder = os.path.join(root_folder, country, city, 'osm')
    tags_water = {
        'waterway': True,
        'water': True,
        'natural': ['water', 'wetland', 'bay']
    }

    tags_roads = {
        'highway': True
    }

    # tags_buildings = {
    #     'building': True
    # }

    args = [{'name': 'water',
             'tags': tags_water},
            {'name': 'roads',
             'tags': tags_roads},
            # {'name': 'buildings',
            #  'tags': tags_buildings}
            ]


    df = gpd.read_file(area)
    df = df.to_crs(4326)
    # explore to assert single-part polygons
    df = df.explode(index_parts=False)
    df.to_file(os.path.join(output_folder, 'temp.shp'))

    overpass_url = overpass_url.strip()
    print(f"Configuring overpass endpoint as {overpass_url} with timeout {overpass_timeout} seconds")

    # OSMnx v2 uses overpass_url; older versions used overpass_endpoint.
    osmnx.settings.overpass_url = overpass_url
    if hasattr(osmnx.settings, 'overpass_endpoint'):
        osmnx.settings.overpass_endpoint = overpass_url
    osmnx.settings.requests_timeout = overpass_timeout
    osmnx.settings.use_cache = True
    osmnx.settings.log_console = True

    Multi = None
    with fiona.open(os.path.join(output_folder, 'temp.shp')) as ar:
        # convert to single mult-part polygon, but this fails if it is already a multi-polygon        
        Multi = MultiPolygon([shape(pol['geometry']) for pol in ar])

    for arg in args:
        get_osm_part(Multi, arg['tags'], arg['name'], output_folder)


def extract_flood_hazard_data(flood_hazard_data_folder, country, city, root_folder, file_to_probability_mapping: Dict[str, int] = None):
    
    # create the temp path if it doesn't exist
    temp_folder = os.path.join(root_folder, country, city, 'temp')
    os.makedirs(temp_folder, exist_ok=True)  
    
    dest_folder = os.path.join(root_folder, country, city, 'fathom')
    src_folder = flood_hazard_data_folder
    names = [file for file in os.listdir(src_folder) if file.endswith(".tif")]
    mask_file = os.path.join(root_folder, country, city, 'eo_indicators', 'dem.tif')

    ds_template = rds.open_rasterio(mask_file)

    for name in names:
        logger.info(f'Extracting {name} to {dest_folder}')
        try:
            old_file = os.path.join(src_folder, name)

            # use new_name for destination file. Map new_name to probability, if mapping is provided
            new_name = name
            
            if file_to_probability_mapping is not None:
                # get the probability from the mapping in a failsafe way
                prob = file_to_probability_mapping.get(name)
                if prob is not None:
                    # create new name with probability
                    new_name = f'P_1in{prob}.tif'
                else:
                    # if we have provided a mapping schema, but the file is not in the mapping, exit with error
                    logger.error(f'File {name} not in mapping schema. Please check the mapping schema.')
                    sys.exit(1)

            new_file = os.path.join(dest_folder, new_name)

            # match template file using rioxarray: this seems to handle different CRS better?
            ds_src = rds.open_rasterio(old_file)
            ds_matched = ds_src.rio.reproject_match(ds_template, align_on_exact=True)
            ds_matched.rio.to_raster(new_file, driver='GTiff', dtype='float32')
        
        except Exception as e:
            logger.warning(f'could not extract {name} to {dest_folder}: {e}')


def get_hexagons(local_root, country, city, resolution=8):

    aoi = os.path.join(local_root, country, city, 'ancillary', f'{city}_aoi.shp')
    output = os.path.join(local_root, country, city, 'ancillary', f'{city}_hex.shp')

    aoi_gdf = gpd.read_file(aoi)

    ### add a buffer of 1km to ensure all area is covered with hexagons
    # temporarily transform to projected mercator
    aoi_gdf_buf = aoi_gdf.to_crs(3395)
    aoi_gdf_buf['geometry'] = aoi_gdf_buf.geometry.buffer(1000)

    # transform back to 4326
    aoi_gdf_buf = aoi_gdf_buf.to_crs(4326)

    # create hexagons for aoi
    hexagons = aoi_gdf_buf.h3.polyfill_resample(resolution)
    # common projection needed for clipping: assert same crs as input layer
    hexagons = hexagons.to_crs(aoi_gdf.crs)

    # clip hexagons with aoi
    clipped_hexagons = gpd.clip(hexagons, aoi_gdf)
    # add neigh_id index for later processing
    idx = range(1, len(clipped_hexagons) + 1)
    clipped_hexagons['neigh_id'] = idx
    clipped_hexagons = clipped_hexagons[['neigh_id', 'geometry']]
    # save output file
    clipped_hexagons.to_file(output)

