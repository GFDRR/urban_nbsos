from . import gis_tools
import faulthandler
faulthandler.enable()
import rasterio as rio
import numpy as np

from pysheds.grid import Grid
from skimage.measure import label, regionprops, regionprops_table
from skimage.util import map_array
import matplotlib
import geopandas as gpd
from rasterio import features
import concurrent.futures
import multiprocessing as mp
import gc

from tqdm import tqdm
from tqdm.auto import tqdm as tqauto

import os
matplotlib.use('TkAgg')

import faulthandler
faulthandler.enable()

class Shedder:

    def __init__(self, dem_path):
        """
        Calculate catchment areas for pour point using pysheds
        :param dem_path:
            Path to DEM
        """

        nodata_value = np.array([np.nan], dtype=np.float32)[0]
        self.grid = Grid.from_raster(dem_path, nodata = nodata_value)
        self.dem = self.grid.read_raster(dem_path, nodata=self.grid.nodata)
        self.dirmap = (64, 128, 1, 2, 4, 8, 16, 32)

    def handle_dem(self):
        """
        Prepare DEM for calculations (fill depressions and pits etc.)
        :return:
        """
        pit_filled_dem = self.grid.fill_pits(self.dem)
        flooded_dem = self.grid.fill_depressions(pit_filled_dem)
        self.inflated_dem = self.grid.resolve_flats(flooded_dem)

    def get_acc(self):
        """
        Get the accumulation grid
        :return:
        """
        self.fdir = self.grid.flowdir(self.inflated_dem, dirmap=self.dirmap)
        self.acc = self.grid.accumulation(self.fdir, dirmap=self.dirmap)
        return self.acc

    def get_catchment(self, x, y):
        catch = self.grid.catchment(x=x, y=y, fdir=self.fdir, dirmap=self.dirmap, xytype='index')
        return catch




class Storage:

    def __init__(self, local_root, country, city, nbs):
        self.root = os.path.join(local_root, country, city)
        self.city = city
        self.nbs_mask = os.path.join(self.root, 'suitability', f'{nbs}_creation.tif')
        self.cn_baseline = os.path.join(self.root, 'cn', f'cn_baseline.tif')
        self.dem_path = os.path.join(self.root, 'eo_indicators', f'dem.tif')
        with rio.open(self.nbs_mask) as src:
            self.nbs = src.read(1)
        with rio.open(self.cn_baseline) as src:
            self.cn = src.read(1)
        with rio.open(self.dem_path) as src:
            self.dem = src.read(1)
        self.shedder = Shedder(self.dem_path)
        self.labels = []
        self.invest = InvestFlood()

    def get_acc(self):
        self.shedder.handle_dem()
        acc = self.shedder.get_acc()
        self.acc = np.array(acc)
        return self.acc


    def label_image(self, min_size):
        l = label(self.nbs)
        table = regionprops_table(
            l,
            properties=('label', 'area')
        )
        condition = (table['area'] >= min_size)

        # zero out labels not meeting condition
        input_labels = table['label']
        output_labels = input_labels * condition

        new_image = map_array(
            l, input_labels, output_labels
        )
        return new_image

    def label_hex(self):
        """
        This method takes two arrays of the hexagon id  and labels,
        converts them to string arrays,
        concatenates them with an underscore in between,
        finds unique combinations, maps original arrays to corresponding unique indexes,
        and returns the mapped array.
        """
        labels_large = self.label_image(100)
        labels_small = self.label_image(10)
        hex = self.rasterize_hexagons()
        hex_labels = np.where(labels_large > 0, hex, 0)

        hex_str = hex_labels.astype(str)
        labels_str = labels_small.astype(str)
        result_array = np.char.add(np.char.add(hex_str, '_'), labels_str)
        unique_strings, unique_indexes = np.unique(result_array, return_inverse=True)
        mapped_array = unique_indexes.reshape(result_array.shape)

        return mapped_array

    def get_hexagons(self, resolution=9):
        aoi = os.path.join(self.root, 'ancillary', f'{self.city}_aoi.shp')
        output = os.path.join(self.root, 'ancillary', f'{self.city}_hex_{resolution}.shp')

        aoi_gdf = gpd.read_file(aoi)

        ### add a buffer of 1km to ensure all area is covered with hexagons
        # temporarily transform to projected mercator
        aoi_gdf_buf = aoi_gdf.to_crs(3395)
        aoi_gdf_buf['geometry'] = aoi_gdf_buf.geometry.buffer(1000)
        # transform back to 4326
        aoi_gdf_buf = aoi_gdf_buf.to_crs(4326)

        # create hexagons for aoi
        hexagons = aoi_gdf_buf.h3.polyfill_resample(resolution)

        # clip hexagons with aoi
        clipped_hexagons = gpd.clip(hexagons, aoi_gdf)
        # add neigh_id index for later processing
        idx = range(1, len(clipped_hexagons) + 1)
        clipped_hexagons['neigh_id'] = idx
        clipped_hexagons = clipped_hexagons[['neigh_id', 'geometry']]
        # save output file
        clipped_hexagons.to_file(output)
        return output

    def rasterize_hexagons(self):
        hexagons_path = self.get_hexagons()
        land_cover_p = os.path.join(self.root, 'eo_indicators', f'land_cover.tif')
        hexagons = gpd.read_file(hexagons_path)
        hex = os.path.join(self.root, 'ancillary', f'{self.city}_aoi_hex.tif')
        with rio.open(land_cover_p) as src:
            meta = src.meta.copy()
            meta.update(compress='lzw',
                        dtype='uint32')

        
        with rio.open(hex) as src:
            hex_ras = src.read(1)

        return hex_ras

    def get_label_names(self):
        self.nbs_image = self.label_hex()
        self.labels = np.unique(self.nbs_image)[1:]
        self.acc = self.get_acc()


    def get_all_runoff(self, divider=0.8, max_rel=0.03, max_a=100):
        
        # result raster initialized with zero
        sum_raster = np.zeros(self.nbs_image.shape)

        # iterate over neighbourhoods
        with tqauto(total=len(self.labels), desc="Calculating runoff") as p:

            def get_single_area(label1, max_a, max_rel, divider):
                try:
                    # get nbs cluster
                    acc = self.acc
                    label_mask = self.nbs_image == label1
                    if not np.any(label_mask):
                        tqauto.write(f"Skipping label {label1}: no cluster pixels")
                        return None

                    flat_indices = np.flatnonzero(label_mask)
                    label_acc = acc.ravel()[flat_indices]
                    if label_acc.size == 0 or label_acc.max() == 0:
                        tqauto.write(f"Skipping label {label1}: no valid accumulation pixel")
                        return None

                    max_idx_in_label = np.argmax(label_acc)
                    idx = np.unravel_index(flat_indices[max_idx_in_label], acc.shape)
                    # skip pour points on the raster boundary (pysheds segfaults on edge cells)
                    if idx[0] == 0 or idx[1] == 0 or idx[0] == acc.shape[0] - 1 or idx[1] == acc.shape[1] - 1:
                        tqauto.write(f"Skipping label {label1}: pour point on raster boundary {idx}")
                        return None

                    # calculate catchment and cn
                    catch = self.shedder.get_catchment(idx[1], idx[0])
                    catch2 = np.array(catch)

                    valid_cn_mask = (catch2 > 0) & (self.cn > 0)
                    if not np.any(valid_cn_mask):
                        tqauto.write(f"Skipping label {label1}: no valid CN values in catchment")
                        return None

                    mean_cn = self.cn[valid_cn_mask].mean()

                    # calculate maximum catchment area
                    nbs_area = np.count_nonzero(label_mask)
                    max_area = min(nbs_area / max_rel, max_a)
                    area = min(np.count_nonzero(valid_cn_mask), max_area)

                    # recalculate CN
                    c = 0.0138 * mean_cn - 0.4586
                    new_c = c * divider
                    new_cn = 69.908 * new_c + 34.457

                    # Calculate runoff difference with baseline
                    old_runoff = self.invest._calc_runoff(mean_cn)
                    new_runoff = self.invest._calc_runoff(new_cn)
                    dif_runoff = old_runoff - new_runoff
                    tot_runoff = dif_runoff * area
                    runoff = tot_runoff / nbs_area
                    return (label1, runoff)
                
                except Exception as e:
                    tqauto.write(f"Error processing label {label1}: {e}")
                    return None


            with concurrent.futures.ThreadPoolExecutor(max_workers=mp.cpu_count()) as executor:
                futures = [executor.submit(get_single_area, label, max_a, max_rel, divider) for label in self.labels]

                # Iterate over futures as they complete
                try:
                    for future in concurrent.futures.as_completed(futures):
                        result = future.result() 
                        if result is not None:
                            label1, runoff = result
                            sum_raster[self.nbs_image == label1] += runoff

                        # Ensure each completed future can be reclaimed immediately.
                        del result
                        p.update(1)

                except Exception as e:
                    print(f"Generated an exception: {e}")
                    return None

        # Explicitly release task artifacts before returning large output array.
        del futures
        gc.collect()

        return sum_raster

    def calc_benefits(self, sum_raster, demand_ras):
        with rio.open(demand_ras) as src:
            demand = src.read(1)

        if sum_raster.shape[0] < demand.shape[0]:
            sum_raster = np.vstack((sum_raster, np.zeros((demand.shape[0] - sum_raster.shape[0], demand.shape[1]))))
        elif sum_raster.shape[0] > demand.shape[0]:
            demand = np.vstack((demand, np.zeros((sum_raster.shape[0] - demand.shape[0], sum_raster.shape[1]))))

        if sum_raster.shape[1] < demand.shape[1]:
            sum_raster = np.hstack((sum_raster, np.zeros((demand.shape[0], demand.shape[1] - sum_raster.shape[1]))))
        elif sum_raster.shape[1] > demand.shape[1]:
            demand = np.hstack((demand, np.zeros((sum_raster.shape[0], sum_raster.shape[1] - demand.shape[1]))))

        benefit = sum_raster * demand
        return benefit




class InvestFlood:
    """
    Calculates runoff, retention, runoff volume and retention volume. Input needed
    is storm depth in mm and a raster with CN-numbers. All equations are based on
    Urban Invest flood mitigation module, found at
    https://invest-userguide.readthedocs.io/en/latest/urban_flood_mitigation.html

    """

    def __init__(self):
        pass

    def _calc_smax(self, cn, p_1=25400, p_2=254):
        c = p_1 / cn
        smax = c - p_2
        return smax

    def _calc_runoff(self,cn, p=200, lam=0.2):
        """
        Conditional function from InVest to calculate runoff as a function of
        Curve Number.
        https://invest-userguide.readthedocs.io/en/latest/urban_flood_mitigation.html
        If storm depth (p) > Smax * Lambda runoff is 0
        :param p: Storm depth in mm
        :param lam: Lambda parameter, set default to 0.2 by InVest
        :param smax: np array with Smax calculated by calc_smax
        :return: np array with runoff in mm
        """
        smax = self._calc_smax(cn)
        new_smax = smax * lam
        a = (p - new_smax) ** 2
        b = p + (1 - lam) * smax

        q = a / b
        if (p > new_smax) & (q > 0):
            q = q
        else:
            q = 0

        return q

