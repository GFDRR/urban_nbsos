from . import gis_tools
from .scan_context import ScanContext


import numpy as np
import pandas as pd
import rasterio as rio
import os
import random
import jenkspy
from typing import List, Dict, Tuple
import rasterio
from tqdm.auto import tqdm as tqauto
from tqdm import tqdm
import geopandas as gpd
import logging



logger = logging.getLogger('nbsos.scan')


class CreationCluster:

    def __init__(self, local_root, country, city):
        self.local_root = local_root
        self.country = country
        self.city = city
        self.lppi_p = os.path.join(local_root, country, city, 'eo_indicators', 'lppi.tif')
        self.bs_freq_p = os.path.join(local_root, country, city, 'eo_indicators', 'bs_freq.tif')
        self.nbs_status = None

        self.nan_mask: np.ndarray = None

        ## Create a random sample of 10000 pixels for Jenks Natural Breaks data clustering
        with rio.open(self.lppi_p) as src:
            lppi = src.read(1)                        
            self.nan_mask = np.isnan(lppi)
            self.rst_shape = lppi.shape

            # correct for changes in lppi layer introduced recently
            lppi_flat = lppi[~self.nan_mask].flatten()
            
        self.sample_pixels = random.sample(range(len(lppi_flat)), 10000)


    def get_jenks_clusters(self, path, inverse: bool = False):
        with rio.open(path) as src:
            image = src.read(1)  
        shape = image.shape          

        # flatten image for sampling
        image_flat = image[~self.nan_mask].flatten()
        
        # sample values
        image_sample = image_flat[self.sample_pixels]
        # determine breaks
        image_breaks = jenkspy.jenks_breaks(image_sample, n_classes=6)

        vals = [1, 2, 3, 4, 5]
        if inverse:
            vals = [5, 4, 3, 2, 1]

        out_image = np.zeros(shape, dtype=np.float32)
        out_image = np.where( (image < image_breaks[2]), vals[0], out_image)
        out_image = np.where( (image >= image_breaks[2]) & (image < image_breaks[3]), vals[1], out_image)
        out_image = np.where( (image >= image_breaks[3]) & (image < image_breaks[4]), vals[2], out_image)
        out_image = np.where( (image >= image_breaks[4]) & (image < image_breaks[5]), vals[3], out_image)
        out_image = np.where( (image >= image_breaks[5]), vals[4], out_image)

        return out_image

    def get_jenks_image(self, path, sample_n=10000):

        with rio.open(path) as src:
            image = src.read(1)            

        # flatten image for sampling
        image_flat = image[~self.nan_mask].flatten()
        
        # sample values
        image_sample = image_flat[self.sample_pixels]
        # determine breaks
        image_breaks = jenkspy.jenks_breaks(image_sample, n_classes=5)
        # fill with cluster value where 1 = most degraded land and 5 = least degraded land
        image_clusters = np.where(
            image < image_breaks[1], 1, np.where(
                image < image_breaks[2], 2 , np.where(
                    image < image_breaks[3], 3, np.where(
                        image < image_breaks[4], 4, 5
                    ))))

        return image_clusters

    def get_suitability_mask(self):
        s_folder = os.path.join(self.local_root, self.country, self.city, 'suitability')
        suitability_masks = [os.path.join(s_folder, f) for f in os.listdir(s_folder) if f.endswith('suitability.tif')]
        masks = []
        for raster in suitability_masks:
            with rio.open(raster) as src:
                masks.append(src.read(1))
        stack = np.where(np.nansum(masks, 0) > 0, 1, 0)
        return stack

    def combine_clusters(self, threshold=7):
        lppi_clusters = self.get_jenks_image(self.lppi_p)
        bs_freq_clusters = self.get_jenks_image(self.bs_freq_p)
        jenks_clusters = lppi_clusters + (6 - bs_freq_clusters)
        clusters_protection_creation = np.where(jenks_clusters <= threshold, 2, 1)
        return clusters_protection_creation

    def make_clusters_for_exploration(self) -> np.ndarray:
        
        template = os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'lppi.tif')
        with rasterio.open(template) as src:
            meta = src.meta.copy()
           
        meta.update({
            'dtype': np.float32,
            'nodata': np.nan
        })

        lppi_clusters = self.get_jenks_clusters(self.lppi_p)
        bs_freq_clusters = self.get_jenks_clusters(self.bs_freq_p, inverse=True)
        result = lppi_clusters + bs_freq_clusters
        result[self.nan_mask] = np.nan

        mask_path = os.path.join(self.local_root, self.country, self.city, 'cluster_maps', 'clusters.tif')
        gis_tools.export_raster(result, template, mask_path, meta)

        return result
    
   
    def get_protection_mask(self, threshold=7):
        nbs_suitability_f = os.path.join(self.local_root, self.country, self.city,
                                       'suitability')
        nbs_suitability = [os.path.join(nbs_suitability_f,f) for f in os.listdir(nbs_suitability_f) if
                             f.endswith('suitability.tif')][0]

        bs_protection_p = os.path.join(self.local_root, self.country, self.city,
                                      'suitability', f"nbs_protection.tif")

        clusters = self.combine_clusters(threshold=threshold)
        nbs = self.get_suitability_mask()
        nbs_status = np.where(nbs == 1, clusters, 0)
        self.nbs_status = nbs_status

        nbs_protection = np.where(nbs_status == 1, 1, 0)
        gis_tools.export_raster(nbs_protection, nbs_suitability, bs_protection_p)

        
        bs_clusters = os.path.join(self.local_root, self.country, self.city,
                                      'suitability', f"nbs_protection_clustering.tif")
        gis_tools.export_raster(clusters, nbs_suitability, bs_clusters)


    def nbs_creation(self, nbs):
        nbs_suitability = os.path.join(self.local_root, self.country, self.city,
                                       'suitability', f"{nbs}_suitability.tif")
        nbs_creation_p = os.path.join(self.local_root, self.country, self.city,
                                       'suitability', f"{nbs}_creation.tif")
        nbs_creation_perc_p = os.path.join(self.local_root, self.country, self.city,
                                      'suitability', f"{nbs}_creation_percentage.tif")

        with rio.open(nbs_suitability) as src:
            suitability = src.read(1).astype(np.int32)

        unique, counts = np.unique(suitability, return_counts=True)
        df = pd.DataFrame({'id':unique, 'total':counts})

        nbs_creation = np.where((self.nbs_status == 2) & (suitability > 0), suitability, 0).astype(np.int32)
        nbs_creation_layer = np.where((self.nbs_status == 2) & (suitability > 0), 1, 0)

        if np.nanmax(nbs_creation_layer) == 0:
            logger.info(f'No creation for {nbs}.')

        gis_tools.export_raster(nbs_creation_layer, nbs_suitability, nbs_creation_p)

    
    def apply_creation_masks(self, nbs_list: List[str], mask_layer_name: str) -> None:
        
        """Apply masks to creation opportunities"""
        # get mask from shapefile
        mask_path = os.path.join(self.local_root, self.country, self.city, 'masks', f"{mask_layer_name}")

        # template is landcover 
        templ = self.lppi_p
        meta = gis_tools.template(templ)
        

        # import shapefile and match projection in layer crs
        mask_gdf = gpd.read_file(mask_path)
        mask_gdf = mask_gdf.to_crs(meta['crs'])

        # rasterize shapefile
        rasterized_mask_path = os.path.join(self.local_root, self.country, self.city, 'masks', f"{mask_layer_name.replace('.shp', '.tif')}")
        gis_tools.rasterize_geom(mask_gdf, templ, rasterized_mask_path, burn_value=1, fill_value=0)
        
        with rio.open(rasterized_mask_path) as r:
            mask_mtx = r.read(1)
            

        for nbs in nbs_list:

            logger.info(f"Applying mask {mask_path} to {nbs} creation opportunities")            
            nbs_creation_path = os.path.join(self.local_root, self.country, self.city, 'suitability', f"{nbs}_creation.tif")

            with rio.open(nbs_creation_path) as r:
                creation_mtx = r.read(1)

            nbs_creation_masked = np.where(mask_mtx == 1, 0, creation_mtx)
            gis_tools.export_raster(nbs_creation_masked, templ, nbs_creation_path)




def upload_folder(bucket, local_root, country, city, folder):
    loc_folder = os.path.join(local_root, country, city, folder)
    suitability_files = [f for f in os.listdir(loc_folder) if f.endswith('.tif')]
    locals = [os.path.join(loc_folder, f) for f in suitability_files]
    clouds = [f'nbsos/{country}/{city}/{folder}/{f}' for f in suitability_files]
    for local, cloud in zip(locals, clouds):
        blob = bucket.blob(cloud)
        blob.upload_from_filename(local)

