from . import gis_tools

import numpy as np
import rasterio as rio
import os
import rioxarray as riox
from rasterio.enums import Resampling

class HeatProtection:

    def __init__(self, local_root, country, city, protection_mask_type: str, use_wri_tree_layer: bool = True):
        
        self.local_root = local_root
        self.country = country
        self.city = city
        self.hrsl_pop = os.path.join(local_root, country, city, 'eo_indicators', 'hrsl_pop.tif')
        self.land_cover = os.path.join(local_root, country, city, 'eo_indicators', 'land_cover.tif')
        self.airT = os.path.join(local_root, country, city, 'Heat', 'outputs', 'airT_10m.tif')
        self.airT_res = os.path.join(local_root, country, city, 'Heat', 'outputs', 'airT_10m_resampled.tif')
        self.heat_protection = os.path.join(local_root, country, city, 'output', 'protection', 'nbs_protection_heat_benefits.tif')
        self.protection_mask_type = protection_mask_type
        self.use_wri_tree_layer = use_wri_tree_layer

        if self.protection_mask_type == "suitability":
            self.protection = os.path.join(local_root, country, city, 'suitability', 'nbs_protection.tif')
        elif self.protection_mask_type == "lppi":
            self.protection = os.path.join(local_root, country, city, 'suitability', 'nbs_protection_clustering.tif')

        ds_air_temp_orig = riox.open_rasterio(self.airT)
        ds_air_temp_reproj = ds_air_temp_orig.rio.reproject_match(riox.open_rasterio(self.hrsl_pop), resampling=Resampling.nearest)       
        gis_tools.export_raster_dataset(ds_air_temp_reproj, self.airT_res)

    def check_wri(self):
        """
        Checks for the presence of tree coverage data in a  WRI (World Resources Institute)
        tree data GeoTIFF file and evaluates its contents.

        This method attempts to load tree coverage data. It calculates
        the sum of all non-NaN values within the data and checks whether tree coverage data is
        present. If valid data exists, the method returns the tree coverage array; otherwise, it
        returns False.

        Raises:
            RasterioIOError: If the GeoTIFF file cannot be opened or does not exist.

        Attributes:
            wri (str): The file path to the WRI tree coverage GeoTIFF file.
            land_cover (str): Path denoting the land cover tree data for a specific country and
            city.

        Returns:
            numpy.ndarray or None: Returns a multi-dimensional array representation of the tree
            coverage data if valid data is present, otherwise returns False.
        """
        self.wri = os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'wri_trees.tif') # self.land_cover = 
        try:
            with rio.open(self.wri) as src:
                wri = src.read(1)
            
            wri_sum = np.nansum(wri)
            if wri_sum > 0:
                return wri
            else:
                return None
        
        except rio.errors.RasterioIOError:
            return None

    def calc_heat_protection(self):
        """
        Calculate the heat protection benefits of a given area based on protection raster,
        population raster, air temperature raster, and land cover raster.
        The output raster reflects the heat protection value derived from these parameters.

        Attributes:
            protection (str): Filepath to the protection raster, indicating areas currently
                under greenspace.
            hrsl_pop (str): Filepath to the population raster.
            airT_res (str): Filepath to the air temperature raster, generally defining
                creation opportunities.
            land_cover (str): Filepath to the land cover raster, indicating land classifications.
            heat_protection (str): Output filepath where the calculated heat protection
                raster will be saved.

        Raises:
            Exception: Will be raised if the raster files provided fail to open or if errors
                occur during calculations.

        Notes:
            This method utilizes rasterio for importing raster datasets and NumPy for spatial
            operations. Additionally, helper functions from gis_tools are used for focal
            statistics calculation and raster export.
        """
        with rio.open(self.protection) as src:
            # protection raster, i.e. where there is current greenspace
            protection = src.read(1)
        
        if self.protection_mask_type == "suitability":
            protection = np.where(protection > 0, 1.0, 0)  # make sure that we have the correct value for protection
        elif self.protection_mask_type == "lppi":
            protection = np.where(protection == 1, 1.0, 0)  # make sure that we have the correct value for protection
        
        with rio.open(self.hrsl_pop) as src:
            # population raster
            hrsl_pop = src.read(1)
        with rio.open(self.airT_res) as src:
            airt = src.read(1)
        with rio.open(self.land_cover) as src:
            land_cover = src.read(1)

        
        focal_pop = gis_tools.focal_stat(hrsl_pop, 9, in_file='np')
        
        protection_benefit = airt * focal_pop * protection
        
        wri = self.check_wri()
        if wri is not None and self.use_wri_tree_layer:
            protection_benefit = protection_benefit * (wri + 10.0)
        
        else:
            # double heat benefits at land cover equals trees
            lc3 = np.where(land_cover == 3, 2.0, 1.0)
            protection_benefit = protection_benefit * lc3

            
        gis_tools.export_raster(protection_benefit, self.hrsl_pop, self.heat_protection)
