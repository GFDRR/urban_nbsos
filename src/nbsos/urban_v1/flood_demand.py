from . import gis_tools

import numpy as np
import rasterstats
from pysheds.grid import Grid
from tqdm import tqdm
from tqdm.auto import tqdm as tqauto
import rasterio as rio
import os
import geopandas as gpd
import pandas as pd
import exactextract as ee
import rioxarray
from typing import Dict, List
from rasterio.warp import Resampling

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

    def filter_acc(self, min=1000):
        """
        Filter acculmulation grid with given miniuum
        :param min: minimum of accumulation (unit is cells)
        :return:
        """
        acc = np.asarray(self.acc)
        acc_min = np.where(acc < min, 0, acc)
        return acc_min

    def sum_catchments(self, acc_points):
        """
        A little tricky method to sum the values of all catchments
        :param acc_points: points with accumulation higher than threshold
        :return:
        """
        bbox = self.grid.bbox
        
        # get all points in input file that are within the bbox (if there aer points without it things go wrong)
        ## internal pysheds coords = [a for a in acc_points if (bbox[0] < a['x'] < bbox[2]) & (bbox[1] < a['y'] < bbox[3])]
        
        min_x = min(bbox[0], bbox[2])
        max_x = max(bbox[0], bbox[2])
        min_y = min(bbox[1], bbox[3])
        max_y = max(bbox[1], bbox[3])

        coords = [a for a in acc_points if (min_x < a['x'] < max_x) & (min_y < a['y'] < max_y)]
        neigh_ids = list((set([a['neigh_id'] for a in coords])))

        #start the process with the first neighborhood to have a canvas
        first_neigh = [record for record in coords if record['neigh_id'] == neigh_ids[0]]
        first_neigh_record = first_neigh[0]

        col, row = self.grid.nearest_cell(first_neigh_record['x'], first_neigh_record['y'])
        catch1 = self.grid.catchment(x=col,#first_neigh_record['x'],
                                     y=row,#first_neigh_record['y'],
                                     fdir=self.fdir, dirmap=self.dirmap,
                                     xytype='index') # ,nodata_out=np.int64(0)
        
        # add other points in first neighborhood
        for rec in first_neigh[1:]:
            
            col_rec, row_rec = self.grid.nearest_cell(rec['x'], rec['y'])
            catch2 = self.grid.catchment(x=col_rec,#rec['x'],
                                         y=row_rec,#rec['y'],
                                         fdir=self.fdir, dirmap=self.dirmap,
                                         xytype='index') # coordinate # ,nodata_out=np.int64(0)
            catch1 = catch1 | catch2
        service_area1 = catch1.astype(int) * first_neigh_record['risk']

        # add other neighborhoods
        for neigh_id in tqauto(neigh_ids[1:]):
            try:
                records = [record for record in coords if record['neigh_id'] == neigh_id]
                first_neigh_record = records[0]
                col_f, row_f = self.grid.nearest_cell(first_neigh_record['x'], first_neigh_record['y'])
                catch1 = self.grid.catchment(x=col_f, #first_neigh_record['x'],
                                             y=row_f, #first_neigh_record['y'],
                                             fdir=self.fdir, dirmap=self.dirmap,
                                             xytype='index') # coordinate # ,nodata_out=np.int64(0)
                for rec in records[1:]:
                    col_rec, row_rec = self.grid.nearest_cell(rec['x'], rec['y'])
                    catch2 = self.grid.catchment(x=col_rec,#rec['x'],
                                                 y=row_rec,#rec['y'],
                                                 fdir=self.fdir, dirmap=self.dirmap,
                                                 xytype='index') # ,nodata_out=np.int64(0)
                    catch1 = catch1 | catch2
                service_area2 = catch1.astype(int) * first_neigh_record['risk']
                service_area1 += service_area2

            except Exception as e:
                print(e)
                continue

        return service_area1


class AnnualFlood:
    def __init__(self, local_root, country, city, flood_return_periods: List[int], **kwargs):

        # construct flood paths from return periods and the naming schema
        self.flood_paths = {}
        for return_period in flood_return_periods:
            self.flood_paths[f'P_1in{return_period}.tif'] = return_period
   
        self.root = os.path.join(local_root, country, city)
        self.neighborhoods_path = os.path.join(self.root, 'ancillary', f'{city}_hex.shp')
        self.neighborhoods = os.path.join(self.root, 'ancillary', f'{city}_hex_new.shp')
        
        neighborhoods = gpd.read_file(self.neighborhoods_path)
        neighborhoods = neighborhoods.to_crs(4326)
        neighborhoods['neigh_id'] = neighborhoods.index + 1
        neighborhoods.to_file(self.neighborhoods)
        
        self.flood_folder = os.path.join(self.root, 'fathom')
        self.temp_path = os.path.join(local_root, country, city, 'temp')
        os.makedirs(self.temp_path, exist_ok=True)  # create the temp path if it doesn't exist

    def get_annual_flood_prop(self):
        """
        Get annual flood probability from fathom return period rasters
        NB: here we only look at area flooded or not flooded and calculate the annual probability
        for each pixel to be flooded with he inverted return period
        :return:
        """
        f = f'{list(self.flood_paths.items())[0][0]}'# .tif'
        path = os.path.join(self.flood_folder, f)
        with rio.open(path) as src:
            self.annual_flood = np.zeros(src.read(1).shape)
        
        for key, val in tqauto(self.flood_paths.items()):
            f = f'{key}' # .tif'
            path = os.path.join(self.flood_folder, f)
            with rio.open(path) as src:
                fl_raw = src.read(1) # give new raster the value of annual flood probability

            fl_raw1 = np.where(fl_raw == -32767, 0, fl_raw)
            fl_annual = np.where(fl_raw1 > 0, 1 / val, 0)
            self.annual_flood += fl_annual
        self.annual_flood_path = os.path.join(self.flood_folder, 'annual_flood.tif')

        # output is exported for use with rasterstats zonalstats
        gis_tools.export_raster(self.annual_flood, path, self.annual_flood_path)

    def get_population(self):
        """
        Get sum of population per neighborhood
        :param in_shp: shapefile with neighborhoods
        :return:
        """
        
        pop_ras = os.path.join(self.root, 'eo_indicators', 'hrsl_pop.tif')
        
        tmp = gpd.read_file(self.neighborhoods)
        tmp_ras = rioxarray.open_rasterio(pop_ras)
        tmp_ras = tmp_ras.rio.reproject(tmp.crs)

        pop_shp = ee.exact_extract(tmp_ras, tmp, 'sum', include_geom=True, include_cols=['neigh_id'])

        return pop_shp

    def get_flood(self):
        """
        Get mean annual flood probability per neighborhood
        :param in_shp: shapefile with neighborhoods
        :return:
        """

        tmp = gpd.read_file(self.neighborhoods)
        tmp_ras = rioxarray.open_rasterio(self.annual_flood_path)
        tmp_ras = tmp_ras.rio.reproject(tmp.crs)

        flood_shp = ee.exact_extract(tmp_ras, tmp, 'mean', include_geom=True, include_cols=['neigh_id'])

        return flood_shp

    def get_hazard_exposure(self, out_path, key_column='neigh_id'):
        self.get_annual_flood_prop()
        pop = self.get_population()
        flood = self.get_flood()
        
        pop_gdf = gpd.GeoDataFrame.from_features(pop)        
        flood_gdf = gpd.GeoDataFrame.from_features(flood)

        # merge population and flood hazard
        cols_to_use = pop_gdf.columns.difference(flood_gdf.columns)
        c_list = cols_to_use.to_list()
        c_list.append(key_column)
        mer = flood_gdf.merge(pop_gdf[c_list], on=key_column)

        # calcualte risk by exposure * hazard
        mer['risk'] = mer['sum'] * mer['mean']

        # set crs to wgs_84
        mer = mer.set_crs('epsg:4326')

        # export to shapefile
        mer.to_file(out_path)


class SpatialRisk:

    def __init__(self, local_root, country, city, flood_return_periods: List[int] = None):
        """
        """
        self.temp_path = os.path.join(local_root, country, city, 'temp')
        os.makedirs(self.temp_path, exist_ok=True)  # create the temp path if it doesn't exist
        
        self.dem_path = os.path.join(local_root, country, city, 'eo_indicators', 'dem_30.tif')
        self.shedder = Shedder(self.dem_path)

        self.neighborhoods = os.path.join(local_root, country, city, 'ancillary', f'{city}_hex.shp')
        self.neighborhoods_new = os.path.join(local_root, country, city, 'ancillary', f'{city}_hex_new.shp')


        self.out_mask = os.path.join(local_root, country, city, 'eo_indicators', 'dem.tif')
        self.acc_out = os.path.join(self.temp_path, 'acc_min_500.tif')
        self.risk = AnnualFlood(local_root, country, city, flood_return_periods)

    def prepare_shedder(self):
        self.shedder.handle_dem()
        self.shedder.get_acc()

    def calc_risk(self):
        self.risk.get_hazard_exposure(self.neighborhoods)

    def get_acc(self, min_reach=100):
        filtered_acc = self.shedder.filter_acc(min=min_reach)
        gis_tools.export_raster(filtered_acc, self.dem_path, self.acc_out)

    def get_points_dict(self, min_risk=1):

        out_fn = os.path.join(self.temp_path, 'id_neigh.tif')
        out_fn2 = os.path.join(self.temp_path, 'risk_neigh.tif')

        gis_tools.rasterize(self.neighborhoods, self.dem_path, out_fn2, 'risk')
        gis_tools.rasterize(self.neighborhoods, self.dem_path, out_fn, 'neigh_id')

        with rio.open(out_fn) as src:
            neigh_id = src.read(1)

        with rio.open(out_fn2) as src:
            pop = src.read(1)

        with rio.open(self.acc_out) as src:
            acc = src.read(1)

        acc_binom = np.where(acc > 1, 1, 0)
        pop_acc = pop * acc_binom
        id_acc = neigh_id * acc_binom

        gis_tools.export_raster(id_acc, self.dem_path, os.path.join(self.temp_path, 'id_acc_binom.tif'))
        gis_tools.export_raster(pop_acc, self.dem_path, os.path.join(self.temp_path, 'pop_acc_binom.tif'))

        points_id = gis_tools.get_points(os.path.join(self.temp_path, 'id_acc_binom.tif'), 'neigh_id')
        points_pop = gis_tools.get_points(os.path.join(self.temp_path, 'pop_acc_binom.tif'), 'risk')
        acc_val = gis_tools.get_points(self.acc_out, 'acc')

        pop_df = pd.DataFrame(points_pop)
        id_df = pd.DataFrame(points_id)
        acc = pd.DataFrame(acc_val)

        new_df = pd.merge(pop_df, id_df, how='left', on=['x', 'y'])
        new_df = pd.merge(new_df, acc, how='left', on=['x', 'y'])

        new_df['new_acc'] = new_df['acc'] / 10
        new_df['new_acc_int'] = new_df['new_acc'].astype(int)

        ac_unique = new_df.drop_duplicates(subset=['neigh_id', 'new_acc_int'])
        points_filtered = ac_unique[ac_unique['risk'] > min_risk]

        df_dict = points_filtered.to_dict('records')

        return df_dict

    def get_spatial_risk(self, df_dict):
        risk = self.shedder.sum_catchments(df_dict)
        return risk

    def get_risk(self, risk_out, min_reach=500):
        
        risk_out_temp = os.path.join(self.temp_path, 'temp_risk.tif')
        self.prepare_shedder()
        self.calc_risk()
        self.get_acc(min_reach=min_reach)
        df_dict = self.get_points_dict()
        

        risk = self.get_spatial_risk(df_dict)
        gis_tools.export_raster(risk, self.dem_path, risk_out_temp)
        gis_tools.reproj_match(risk_out_temp, self.out_mask, risk_out)

    @staticmethod
    def get_flood_protection_benefits(local_root, country_name, city_name, protection_mask_type: str) -> None:
        """
        Get flood protection benefits
        """
       
        # import protection opportunities
        if protection_mask_type == 'suitability':
            protection_opportunity_filename = os.path.join(local_root, country_name, city_name, 'suitability', 'nbs_protection.tif')
        elif protection_mask_type == 'lppi':
            protection_opportunity_filename = os.path.join(local_root, country_name, city_name, 'suitability', 'nbs_protection_clustering.tif')
            
        ds_protection = rioxarray.open_rasterio(protection_opportunity_filename)
        mtx_prot = ds_protection.values[0]
        # reclasssify as needed
        if protection_mask_type == "suitability":
            mtx_prot = np.where(mtx_prot > 0, 1.0, 0)  # make sure that we have the correct value for protection
        elif protection_mask_type == "lppi":
            mtx_prot = np.where(mtx_prot == 1, 1.0, 0)  # make sure that we have the correct value for protection
        ds_protection.values[0] = mtx_prot

        # import flood risk file
        flood_risk_file = os.path.join(local_root, country_name, city_name, 'output', 'flood_risk.tif')
        ds_flood_risk = rioxarray.open_rasterio(flood_risk_file)

        # match projection to protection datasource
        ds_flood_risk = ds_flood_risk.rio.reproject_match(ds_protection, resampling=Resampling.nearest, align_on_exact=True)

        out_path = os.path.join(local_root, country_name, city_name, 'output', 'protection', 'nbs_protection_pluvial_flood_benefits.tif')
        ds_protection = ds_flood_risk * ds_protection
        ds_protection.rio.to_raster(out_path, driver='GTiff', dtype='float32')
    