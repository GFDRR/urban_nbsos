from .gis_tools import export_raster
from . import loggerconfig


import rasterio as rio
import distancerasters as dr
from scipy.ndimage import binary_dilation, binary_erosion
import numpy as np
import os
import geopandas as gpd
import logging
import rioxarray

logger = logging.getLogger('nbsos.scan')


class GcsToLocal:

    def __init__(self, bucket, country, city, local_root):
        self.bucket = bucket
        self.country = country
        self.city = city
        self.local_root = local_root

    def get_eo_indicators(self):
        prefix = f'nbsos/{self.country}/{self.city}/eo_indicators/'
        local_folder = os.path.join(self.local_root, self.country, self.city, 'eo_indicators')
        os.makedirs(local_folder, exist_ok=True)
        self.get_gcs_to_local(prefix, local_folder)

    def get_gcs_to_local(self, prefix: str, local_folder: str):
        blobs = self.bucket.list_blobs(prefix=prefix)
        for blob in blobs:
            filename = blob.name.split('/')[-1]
            logger.info(f'Downloading local copy of {filename} to {local_folder}')
            try:
                blob.download_to_filename(os.path.join(local_folder, filename))
            except Exception as e:
                logger.error(f'Could not download local copy of {filename} to {local_folder}: {e}')

class local_indicators:

    def __init__(self, local_root, country, city):
        self.local_root = local_root
        # self.bucket = bucket
        self.country = country
        self.city = city
        self.land_cover_path = os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'land_cover.tif')
        # get template info from land cover raster
        with rio.open(self.land_cover_path) as src:
            mask = src.read(1)
            self.affine = src.transform
            self.shape = mask.shape
            self.crs = src.crs


    def get_water_mask(self):
        """
        Get 10m mask for water pixels using OSM and ESA world cover
        :return: 10m water mask
        """

        # get water mask
        with rio.open(self.land_cover_path) as src:
            lc = src.read(1)
            lc_water = np.where(lc == 6, 1, 0)

        water_lines = os.path.join(self.local_root, self.country, self.city, 'osm/osm_water_LineString.shp')
        water_polygons = os.path.join(self.local_root, self.country, self.city, 'osm/osm_water_Polygon.shp')
        
        gdf_lines = gpd.read_file(water_lines).to_crs(self.crs)
        gdf_polygons = gpd.read_file(water_polygons).to_crs(self.crs)

        lines, _ = dr.rasterize(gdf_lines, affine=self.affine, shape=self.shape)
        polygons, _ = dr.rasterize(gdf_polygons, affine=self.affine, shape=self.shape)
        
        all_water = lc_water + lines + polygons


        return all_water


    def get_building_mask(self):

        google_buildings = os.path.join(self.local_root, self.country, self.city, 'eo_indicators/google_buildings.tif')
        wsf_buildings = os.path.join(self.local_root, self.country, self.city, 'eo_indicators/wsf2019.tif')
        #osm_buildings = os.path.join(self.local_root, self.country, self.city, 'osm/osm_buildings_Polygon.shp')

        #osm, _ = dr.rasterize(osm_buildings, affine=self.affine, shape=self.shape)

        with rio.open(google_buildings) as src:
            google = src.read(1)
        with rio.open(wsf_buildings) as src:
            wsf = src.read(1)

        all_buildings = google + wsf #+ osm # USE osm if there is no google buildings

        return all_buildings

    def get_major_road_mask(self):
        osm_roads = os.path.join(self.local_root, self.country, self.city, 'osm/osm_roads_LineString.shp')
        major_roads = ['motorway', 'primary', 'primary_link', 'road',
                       'secondary', 'secondary_link', 'tertiary',
                       'tertiary_link', 'trunk', 'trunk_link']

        roads_gdf = gpd.read_file(osm_roads).to_crs(self.crs)
        major_roads_gdf = roads_gdf[roads_gdf['highway'].isin(major_roads)]
        lines, _ = dr.rasterize(major_roads_gdf, affine=self.affine, shape=self.shape)
        return lines

    def get_road_mask(self):
        osm_roads = os.path.join(self.local_root, self.country, self.city, 'osm/osm_roads_LineString.shp')
        roads_gdf = gpd.read_file(osm_roads)
        roads_gdf = roads_gdf.to_crs(self.crs)
        lines, _ = dr.rasterize(roads_gdf, affine=self.affine, shape=self.shape)
        return lines

    def dilate(self, image, iterations):
        dilation = binary_dilation(image, iterations=iterations).astype('uint8')
        return dilation

    def erode(self, image, iterations):
        erosion = binary_erosion(image, iterations=iterations).astype('uint8')
        return erosion

    def water(self):
        logger.info('creating water mask using OSM and ESA wrold cover')
        try:
            water = self.get_water_mask()
            out = os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'water_mask.tif')
            export_raster(water, self.land_cover_path, out)
        except Exception as e:
            logger.error(f'Could not create and export water mask: {e}')

    def water_50(self):
        logger.info('creating water mask with 50m buffer using OSM and ESA world cover')
        try:
            water = self.get_water_mask()
            water_50 = self.dilate(water, 5)
            out = os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'water_50m_mask.tif')
            export_raster(water_50, self.land_cover_path, out)
        except Exception as e:
            logger.error(f'Could not create and export 50m buffer water mask: {e}')

    def buildings(self):
        logger.info('creating building mask using OSM and Google Buildings V3')
        try:
            buildings = self.get_building_mask()
            out = os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'building_mask.tif')
            export_raster(buildings, self.land_cover_path, out)
        except Exception as e:
            logger.error(f'Could not create and export building mask using OSM and Google Buildings V3: {e}')

    def buildings_min10(self):
        logger.info('creating building mask using OSM and Google Buildings V3 with 10m erosion for green corridors')
        try:
            buildings = self.get_building_mask()
            buildings_min10 = self.erode(buildings, 1)
            out = os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'building_min10_mask.tif')
            export_raster(buildings_min10, self.land_cover_path, out)
        except Exception as e:
            logger.error(f'Could not create building mask using OSM and Google Buildings V3 with 10m erosion for'
                         f'green corridors: {e}')

    def roads(self):
        logger.info('creating roads mask using OSM')
        try:
            roads = self.get_road_mask()
            out = os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'road_mask.tif')
            export_raster(roads, self.land_cover_path, out)
        except Exception as e:
            logger.error(f'Could not create roads mask: {e}')

    def roads_10(self):
        logger.info('creating major roads mask with 10m buffer for green corridors')
        try:
            roads = self.get_major_road_mask()
            out = os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'major_roads_10m_mask.tif')
            roads_10 = self.dilate(roads, 1)
            export_raster(roads_10, self.land_cover_path, out)
        except Exception as e:
            logger.error(f'Could not create and export major roads mask with 10m buffer for green corridors: {e}')

    def get_local_indicators(self):
        self.water()
        self.water_50()
        self.buildings()
        self.buildings_min10()
        self.roads_10()
        self.roads()




