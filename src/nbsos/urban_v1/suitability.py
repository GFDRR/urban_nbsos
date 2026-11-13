from .gis_tools import export_raster
from . import loggerconfig

#!/usr/bin/env python
import os
import rasterio as rio
import numpy as np
from skimage.measure import label, regionprops_table
from skimage.util import map_array
import geopandas as gpd
from rasterio import features
import logging
import h3pandas
import rioxarray
from rasterio.enums import Resampling

logger = logging.getLogger('nbsos.scan')

class Suitability:

    def __init__(self, country, city, local_root, ignore_annual_precipitation: bool = False):
        self.country = country
        self.city = city
        self.local_root = local_root
        self.ignore_annual_precipitation = ignore_annual_precipitation
        indicator_path = os.path.join(local_root, country, city, 'eo_indicators')

        # conformation message regarding considered suitability criteria
        logger.info(f"Annual precipitation suitability criterion is {'ignored' if self.ignore_annual_precipitation else 'considered'}")


        wsf_mask = os.path.join(indicator_path, 'wsf2019.tif')
        with rio.open(wsf_mask) as src:
            self.wsf = src.read(1)

        # Building mask with Google buildings and OSM
        building_mask_p = os.path.join(indicator_path, 'building_mask.tif')
        with rio.open(building_mask_p) as src:
            b_mask = src.read(1)
            self.building_mask = np.where((b_mask == 1) | (b_mask == 2), 1, 0)

        building_min_10_mask_p = os.path.join(indicator_path, 'building_min10_mask.tif')
        with rio.open(building_min_10_mask_p) as src:
            self.building_min10_mask = src.read(1)

        # annual precipitation from CHIPRS
        annual_precipitation_p = os.path.join(indicator_path, 'annual_prec.tif')
        with rio.open(annual_precipitation_p) as src:
            self.annual_precipitation = src.read(1)

        # DEM from FABDEM
        dem_p = os.path.join(indicator_path, 'dem.tif')
        with rio.open(dem_p) as src:
            self.dem = src.read(1)

        # Slope in degress from FABDEM
        slope_deg_p = os.path.join(indicator_path, 'slope.tif')
        with rio.open(slope_deg_p) as src:
            self.slope_deg = src.read(1)

        # Soil properties, multiband raster
        soil_properties = os.path.join(indicator_path, 'soil_properties.tif')
        with rio.open(soil_properties) as src:
            self.soil_texture = src.read(get_id_from_description(src, 'texture'))
            self.soil_sand = src.read(get_id_from_description(src, 'sand'))
            self.soil_ph = src.read(get_id_from_description(src, 'ph'))
            self.soil_silt = src.read(get_id_from_description(src, 'silt'))
            self.soil_clay = src.read(get_id_from_description(src, 'clay'))
            self.soil_bulk_density = src.read(get_id_from_description(src, 'bulk_density'))
            self.soil_carbon = src.read(get_id_from_description(src, 'soil_carbon'))
            self.soil_depth = src.read(get_id_from_description(src, 'depth'))
            self.soil_stone_content = src.read(get_id_from_description(src, 'stone_content'))

        try:
            # Mask of water with 5 pixel buffer (50m)
            water_50m_mask_p = os.path.join(indicator_path, 'water_50m_mask.tif')
            with rio.open(water_50m_mask_p) as src:
                self.water_50m_mask = src.read(1)

        except Exception as e:
            self.water_50m_mask = np.zeros(self.dem.shape)

        try:
            # Mask of water
            water_mask_p = os.path.join(indicator_path, 'water_mask.tif')
            with rio.open(water_mask_p) as src:
                self.water_mask = src.read(1)
        except Exception as e:
            self.water_mask = np.zeros(self.dem.shape)


        # road mask from OSM
        road_mask_p = os.path.join(indicator_path, 'road_mask.tif')
        with rio.open(road_mask_p) as src:
            self.road_mask = src.read(1)

        # major road mask from OSM with 10m buffer
        road_mask_10_p = os.path.join(indicator_path, 'major_roads_10m_mask.tif')
        with rio.open(road_mask_10_p) as src:
            self.major_roads_10m_mask = src.read(1)

        # Land cover from ESA world Cover
        self.land_cover_p = os.path.join(indicator_path, 'land_cover.tif')
        with rio.open(self.land_cover_p) as src:
            self.land_cover = src.read(1)

        # Land cover from ESA world Cover
        hydrosoil_p = os.path.join(indicator_path, 'hydrosoil.tif')
        with rio.open(hydrosoil_p) as src:
            self.hydrosoil = src.read(1)

        lcc_p = os.path.join(indicator_path, 'lcc.tif')
        with rio.open(lcc_p) as src:
            self.lcc = src.read(1)

        heavy_prec_p = os.path.join(indicator_path, 'heavy_prec.tif')
        with rio.open(heavy_prec_p) as src:
            self.heavy_prec = src.read(1)

        hand_p = os.path.join(indicator_path, 'hand.tif')
        with rio.open(hand_p) as src:
            self.hand = src.read(1)

    def river_renaturation(self):
        logger.info('Identify suitability for river renaturation')
        try:
            riv_nat = np.where(
                (self.building_mask < 1) &
                (self.ignore_annual_precipitation or self.annual_precipitation > 450) &
                #(self.dem < 1500) &
                (self.slope_deg <= 35) &
                ((self.soil_depth >= 75) | (self.soil_depth == 0)) &
                ((self.soil_ph >= 4) | (self.soil_ph == 0)) &
                (self.soil_bulk_density <= 1.6) &
                (self.water_50m_mask == 1) &
                (self.water_mask < 1) &
                (self.road_mask != 1),
                1, 0
            )
            if np.nanmax(riv_nat) != 1:
                logger.warning('no suitability found for river renaturation')
            riv_nat_segmented = self.linear_segmentation(riv_nat)
        except Exception as e:
            logger.error(f'Could not perform suitability identification for river renatration: {e}')

        return riv_nat_segmented

    def green_corridors(self):
        logger.info('Identify suitability for green corridors')
        try:
            green_cor = np.where(
                (self.major_roads_10m_mask == 1) &
                (self.ignore_annual_precipitation or self.annual_precipitation > 450) &
                #(self.dem < 1500) &
                (self.slope_deg <= 35) &
                ((self.soil_depth >= 75) | (self.soil_depth == 0)) &
                (((self.soil_ph >= 4) & (self.soil_ph < 8)) | (self.soil_ph == 0)) &
                (self.soil_bulk_density <= 1.6) &
                (self.water_mask < 1) &
                (self.building_min10_mask != 1),
                1, 0
            )
            if np.nanmax(green_cor) != 1:
                logger.warning('no suitability found for green corridors')
            green_cor_segmented = self.linear_segmentation(green_cor)
        except Exception as e:
            logger.error(f'Could not perform suitability identification for green corridors: {e}')

        return green_cor_segmented

    def urban_forest(self, min_size=25, max_size=1000000):
        logger.info('Identify suitability for urban forest')
        try:
            urb_for = np.where(
                (self.building_mask < 1) &
                (self.ignore_annual_precipitation or self.annual_precipitation > 450) &
                #(self.dem < 1500) &
                (self.slope_deg <= 35) &
                ((self.soil_depth >= 75) | (self.soil_depth == 0)) &
                (((self.soil_ph >= 4) & (self.soil_ph < 8)) | (self.soil_ph == 0)) &
                (self.soil_bulk_density <= 1.6) &
                (self.water_mask < 1) &
                (self.land_cover != 3) &
                (self.land_cover != 7) &
                (self.road_mask != 1),
                1, 0
            )
            if np.nanmax(urb_for) != 1:
                logger.warning('no suitability found for green urban forests')
            urb_for_f = filter_image(urb_for, min_size, max_size)

            return urb_for_f
        except Exception as e:
            logger.error(f'Could not perform suitability identification for urban forest: {e}')

    def bioretention_area(self, min_size=10, max_size=1000):
        logger.info('Identify suitability for bioretention area')
        try:
            bio_ret = np.where(
                (self.building_mask < 1) &
                (self.ignore_annual_precipitation or self.annual_precipitation > 450) &
                #(self.dem < 1500) &
                ((self.soil_depth >= 75) | (self.soil_depth == 0)) & # added by Alida's comment date: 20230918
                (self.slope_deg < 5) &
                (self.hydrosoil <= 3) &
                (self.water_mask < 1) &
                (self.road_mask != 1),
                1, 0
            )
            if np.nanmax(bio_ret) != 1:
                logger.warning('no suitability found for bioretention area')
            bio_ret_f = filter_image(bio_ret, min_size, max_size=max_size)

            return bio_ret_f
        except Exception as e:
            logger.error(f'Could not perform suitability identification for bioretention area: {e}')

    def inland_wetland(self, min_size=10, max_size=1000):
        logger.info('Identify suitability for inland wetland')
        try:
            inl_wet = np.where(
                (self.building_mask < 1) &
                (self.ignore_annual_precipitation or self.annual_precipitation > 450) &
                #(self.dem < 1500) &
                ((self.soil_depth >= 75) | (self.soil_depth == 0)) & # added by Alida's comment date: 20230918
                (self.slope_deg < 5) &
                (self.hydrosoil > 3) &
                (self.water_mask < 1) &
                (self.hand <= 20) &
                (self.road_mask != 1),
                1, 0
            )
            if np.nanmax(inl_wet) != 1:
                logger.warning('no suitability found for inland wetland')
            inl_wet_f = filter_image(inl_wet, min_size, max_size=max_size)

            return inl_wet_f
        except Exception as e:
            logger.error(f'Could not perform suitability identification for bioretention area: {e}')


    def urban_farming(self, min_size=10, max_size=1000):
        logger.info('Identify suitability for urban farming')
        try:
            urb_far = np.where(
                (self.building_mask < 1) &
                (self.ignore_annual_precipitation or self.annual_precipitation > 450) &
                (self.dem < 1500) &
                (self.slope_deg <= 20) &
                (self.soil_texture != 1) &
                (self.soil_texture != 12) &
                (self.water_mask < 1) &
                (self.lcc <= 7) &
                (self.road_mask != 1),
                1, 0
            )
            if np.nanmax(urb_far) != 1:
                logger.warning('no suitability found for urban farming')
            urb_far_f = filter_image(urb_far, min_size, max_size=max_size)

            return urb_far_f
        except Exception as e:
            logger.error(f'Could not perform suitability identification for urban farming: {e}')

    def open_green_space(self, min_size=10, max_size=1000):
        logger.info('Identify suitability for open green space')
        try:
            open_gs = np.where(
                (self.building_mask < 1) &
                (self.ignore_annual_precipitation or self.annual_precipitation > 450) &
                #(self.heavy_prec >= 2) & IAS commented out heavy precipitation for open green space on 231002
                #(self.dem < 1500) &
                (self.slope_deg <= 35) &
                ((self.soil_depth >= 75) | (self.soil_depth == 0)) &
                (((self.soil_ph >= 4) & (self.soil_ph < 8)) | (self.soil_ph == 0)) &
                (self.soil_bulk_density <= 1.6) &
                (self.water_mask < 1) &
                (self.road_mask != 1),
                1, 0
            )
            if np.nanmax(open_gs) != 1:
                logger.warning('no suitability found for open green space')
            open_gs_f = filter_image(open_gs, min_size, max_size=max_size)

            return open_gs_f
        except Exception as e:
            logger.error(f'Could not perform suitability indentification for open green space: {e}')

    
    def terraces_and_slopes(self):
        logger.info('Identify suitability for terraces and slopes')
        try:

            terr_slopes = np.where(
               (self.road_mask == 0) &
               (self.major_roads_10m_mask == 0) &
               ((self.slope_deg > 10) & (self.slope_deg < 45)) &
               (self.wsf == 0),
               1, 0
            )
            if np.nanmax(terr_slopes) != 1:
                logger.warning('no suitability found for terraces and slopes')

            return terr_slopes
        
        except Exception as e:
            logger.error(f'Could not perform suitability indentification for terraces and slopes: {e}')


    def linear_segmentation(self, image):
        hexagons_path = os.path.join(self.local_root, self.country, self.city, 'ancillary', f'{self.city}_hex.shp')
        hexagons = gpd.read_file(hexagons_path)
        hex = os.path.join(self.local_root, self.country, self.city, 'ancillary', f'{self.city}_aoi_hex.tif')
        with rio.open(self.land_cover_p) as src:
            meta = src.meta.copy()
            meta.update(compress='lzw')
        
        with rio.open(hex) as src:
            hex_ras = src.read(1)

        segmented = np.where(image == 1, hex_ras, 0)

        return segmented

    def get_images(self, nbs_list):
        images = []
        if 'open_green_space' in nbs_list:
            d = {'name': 'open_green_space',
                 'image': self.open_green_space()}
            images.append(d)
        if 'urban_forest' in nbs_list:
            d = {'name': 'urban_forest',
                 'image': self.urban_forest()}
            images.append(d)
        if 'urban_farming' in nbs_list:
            d = {'name': 'urban_farming',
                 'image': self.urban_farming()}
            images.append(d)
        if 'bioretention_area' in nbs_list:
            d = {'name': 'bioretention_area',
                 'image': self.bioretention_area()}
            images.append(d)
        if 'inland_wetland' in nbs_list:
            d = {'name': 'inland_wetland',
                 'image': self.inland_wetland()}
            images.append(d)
        if 'river_renaturation' in nbs_list:
            d = {'name': 'river_renaturation',
                 'image': self.river_renaturation()}
            images.append(d)
        if 'green_corridors' in nbs_list:
            d = {'name': 'green_corridors',
                 'image': self.green_corridors()}
            images.append(d)        
        if 'terraces_and_slopes' in nbs_list:
            d = {'name': 'terraces_and_slopes',
                 'image': self.terraces_and_slopes()}
            images.append(d)


        return images

    def export_local(self, nbs_list):
        images = self.get_images(nbs_list)
        for image in images:
            logger.info(f'Exporting suitability for {image["name"]}')
            try:
                out = os.path.join(self.local_root, self.country, self.city, 'suitability',
                                   f'{image["name"]}_suitability.tif')
                export_raster(image['image'], self.land_cover_p, out)
            except Exception as e:
                logger.error(f'Could not export {image["name"]}: {e}')



def get_id_from_description(src, description):
    descr = src.descriptions
    ids = [d+1 for d in range(len(descr))]
    d_dict = {key: val for key, val in zip(descr, ids)}
    value = d_dict[description]

    return value


def filter_image(image, min_size, max_size=None):
    lab = label(image)
    table = regionprops_table(
        lab,
        properties=('label', 'area')
    )
    condition = (table['area'] >= min_size)
    if max_size:
        condition = ((table['area'] >= min_size) & (table['area'] <= max_size))

    # zero out labels not meeting condition
    input_labels = table['label']
    output_labels = input_labels * condition

    new_image = map_array(
        lab, input_labels, output_labels
    )

    if np.nanmax(new_image) < 1:
        logger.warning('no suitability found after filtering')

    # image = np.where(new_image > 0, 1, 0)
    return new_image
