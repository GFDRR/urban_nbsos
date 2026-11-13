from ..scan_context import ScanContext

from .heat_logging import heat_module_logger
from .heat_exports import RegressionPredictors, NbsTreeCoverage

from ..gis_tools import make_gaussian_kernel, make_exact_zonal_statistics, export_raster_dataset, clip, rasterize_df, get_area_crs, reproj_match, rasterize
from ..nbsos_exports import NbsTypes, LandCoverClass

import exactextract as ee
import numpy as np
import rasterio
import rioxarray
import xarray as xr
from rasterio.enums import Resampling
from scipy.ndimage import convolve
import geopandas as gpd
from tqdm import tqdm
from tqdm.auto import tqdm as tqauto
import os
import logging
import warnings
from typing import Dict, List, Tuple
from builtins import open as fopen
import pandas as pd
from skimage import util

import statsmodels.api as sm


# modules for spatial regression
from spreg import ML_Error, GMM_Error
from libpysal.io import open
from libpysal import weights

warnings.filterwarnings("ignore")


class HeatModule():

    
    # data directory name. This is used for writing exports
    module_pathname: str = "Heat"
    nbs_list: List[NbsTypes] = None

    # variables to retrieve data. Those need to be passed to the module at a minimum
    local_root: str = None
    country_name: str = None
    city_name: str = None

    # AOI, hexagons, other ancillary data
    area_of_interest_filename: str = None
    hexagon_filename: str = None
    population_filename: str = None
    baseline_land_cover_filename: str = None
    wri_tree_layer_filename: str = None
    

    # determined and imported datasets
    ds_landcover: xr.Dataset = None
    ds_tree_layer: xr.Dataset = None
    tree_layer_is_valid: bool = None
    ds_mean_air_temperature: xr.Dataset = None
    ds_weighted_population: xr.Dataset = None
    hexagon_layer: gpd.GeoDataFrame = None

    recommended_crs: str = None #"EPSG:3857"

    # parameters for population Gauss weights 
    kernel_radius: float = None
    kernel_sigma: float = None
    # parameters for regression
    y_name: str = None
    x_names: List[str] = None

    regression_weights: None
    regression_predictors: List[RegressionPredictors] = None
    regression_betas: List[float] = None

    # support use of predefined air temperature dataset, if specified
    air_temperature_path: str = None

    def __init__(self, local_root: str, country: str, city: str):
        
        # we assign the pathname to the module from the static method, so that it is consistent across all instances
        # and can also be accessed from other classes/modules from the scan.
        self.module_pathname = HeatModule.get_module_pathname()

        self.local_root = local_root
        self.country_name = country
        self.city_name = city

        required_dirs = [
            'outputs',
            os.path.join('outputs', 'heat_benefits')            
        ]
        for f in required_dirs:
            self.make_folder(f)

        # get started with prepring input data
        heat_module_logger.info(f"Initializing heat scan for {self.city_name}, {self.country_name}")

        # map filenames and assert file existence
        self.map_nbsos_files()

        # determine recommended crs
        self.recommended_crs = get_area_crs(self.area_of_interest_filename)

    
    def map_nbsos_files(self) -> None:
        
        # the AOI file should be local_root/country/city/ancillary/city_AOI.shp
        self.area_of_interest_filename = os.path.join(self.local_root, self.country_name, self.city_name, 'ancillary', f'{self.city_name}_AOI.shp')
        heat_module_logger.info(f"Mapping AOI to file {self.area_of_interest_filename}")
        if not os.path.exists(self.area_of_interest_filename):
            error_msg = f"{self.area_of_interest_filename} not found."
            heat_module_logger.error(error_msg)
            raise FileExistsError(error_msg)
        
        # hexagons should be local_root/country/city/ancillary/city_hex.shp
        self.hexagon_filename = os.path.join(self.local_root, self.country_name, self.city_name, 'ancillary', f'{self.city_name}_hex.shp')
        heat_module_logger.info(f"Mapping hexagons to file {self.hexagon_filename}")
        if not os.path.exists(self.hexagon_filename):
            error_msg = f"{self.hexagon_filename} not found."
            heat_module_logger.error(error_msg)
            raise FileExistsError(error_msg)
        
        # population should be local_root/country/city/eo_indicators/hrsl_pop.tif
        self.population_filename = os.path.join(self.local_root, self.country_name, self.city_name, 'eo_indicators', 'hrsl_pop.tif')
        heat_module_logger.info(f"Mapping population to file {self.hexagon_filename}")
        if not os.path.exists(self.population_filename):
            error_msg = f"{self.population_filename} not found."
            heat_module_logger.error(error_msg)
            raise FileExistsError(error_msg)
        
        # land cover should be local_root/country/city/eo_indicators/land_cover.tif
        self.baseline_land_cover_filename = os.path.join(self.local_root, self.country_name, self.city_name, 'eo_indicators', 'land_cover.tif')
        heat_module_logger.info(f"Mapping baseline land cover to file {self.baseline_land_cover_filename}")
        if not os.path.exists(self.baseline_land_cover_filename):
            error_msg = f"{self.baseline_land_cover_filename} not found."
            heat_module_logger.error(error_msg)
            raise FileExistsError(error_msg)
        
        # tree layer should be local_root/country/city/eo_indicators/wri_trees.tif
        self.wri_tree_layer_filename = os.path.join(self.local_root, self.country_name, self.city_name, 'eo_indicators', f'wri_trees.tif')
        heat_module_logger.info(f"Mapping WRI tree layer to file {self.wri_tree_layer_filename}")
        if not os.path.exists(self.wri_tree_layer_filename):
            error_msg = f"{self.wri_tree_layer_filename} not found."
            heat_module_logger.error(error_msg)
            raise FileExistsError(error_msg)

        # dem 
        self.dem_filename = os.path.join(self.local_root, self.country_name, self.city_name, 'eo_indicators', f'dem.tif')
        heat_module_logger.info(f"Mapping DEM to file {self.dem_filename}")
        if not os.path.exists(self.dem_filename):
            error_msg = f"{self.dem_filename} not found."
            heat_module_logger.error(error_msg)
            raise FileExistsError(error_msg)


    def make_folder(self, folder_name):
        if folder_name is not None:
            target_folder = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, folder_name)
        else:
            target_folder = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname)
        os.makedirs(target_folder, exist_ok=True)

    

    def estimate_weighted_population(self, kernel_radius: float, kernel_sigma: float):

        heat_module_logger.info("Estimating weighted population")

        # import original pop file and project to crs and resample
        ds_population = rioxarray.open_rasterio(self.population_filename)
        
        # mask population outside the AOI 
        ds_aoi_hex = rioxarray.open_rasterio(os.path.join(self.local_root, self.country_name, self.city_name, 'ancillary', f'{self.city_name}_aoi_hex.tif'))
        aoi_mask = np.isin(ds_aoi_hex.data[0], [0], invert=False)
        ds_population.data[0][aoi_mask] = np.nan


        # resample: avoid inconsistent x/y resolution, therefore, adopt determine new resolution along x, and use same scale factor for y
        scale_factor = 3
        population_source_resolution = ds_population.rio.resolution()
        new_resolution = (population_source_resolution[0] * scale_factor)

        population_downsampled = ds_population.rio.reproject(
            dst_crs=ds_population.rio.crs,
            resolution=new_resolution,
            resampling=Resampling.average,
        )
        
        # reproject to recommended crs
        population_downsampled = population_downsampled.rio.reproject(f'EPSG:{self.recommended_crs}')

        # make kernel
        downscaled_resolution = population_downsampled.rio.resolution()
        population_kernel = make_gaussian_kernel(d=kernel_radius, sigma=kernel_sigma, rs=downscaled_resolution)
        population_kernel = population_kernel/np.max(population_kernel) # see R script

        # conduct convolution
        xds_downscaled_mtx = np.nan_to_num(population_downsampled.to_numpy())
        xds_downscaled_mtx = convolve(xds_downscaled_mtx[0], weights=population_kernel, mode='constant', cval=0)
        # replace old data in dataset with modified numpy array
        population_downsampled.data = xds_downscaled_mtx.reshape(1, xds_downscaled_mtx.shape[0], xds_downscaled_mtx.shape[1])

        # reproject to original resolution and export
        self.ds_weighted_population = population_downsampled.rio.reproject_match(ds_population, resampling=Resampling.nearest)
        out_path = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, "weighted_population.tif")
        export_raster_dataset(self.ds_weighted_population, out_path)


    def estimate_mean_air_temperature(self, regression_coefficients_path: str, resampling_method: Resampling = Resampling.cubic):

        heat_module_logger.info("Estimating mean air temperature")
        heat_module_logger.info("Resampling method: " + str(resampling_method))

        self.ds_mean_air_temperature = HeatModule.generate_air_temperature_raster(
            local_root=self.local_root,
            country_name=self.country_name,
            city_name=self.city_name,
            regression_coefficients_path=regression_coefficients_path,
            resampling_method=resampling_method,
            out_path=None
        )    

               

    @staticmethod
    def generate_air_temperature_raster(local_root: str,
                                        country_name: str,
                                        city_name: str, 
                                        regression_coefficients_path: str = None, 
                                        resampling_method: Resampling = Resampling.cubic,
                                        out_path: str = None):

        regression_coefficients_path = os.path.join(local_root, "data", HeatModule.get_module_pathname().lower()) if regression_coefficients_path is None else regression_coefficients_path
        
                    
        # Weighting Scheme for combining Climate Space Weighted Regressions and Geographically Weighted Regression
        gwr_weights = [0.6002, 0.5691, 0.5655, 0.4320, 0.3490, 0.2861, 0.3407, 0.3780, 0.3080, 0.3332, 0.3790, 0.4904]
        cswr_weights = [0.3997, 0.4309, 0.4345, 0.5680, 0.6509, 0.7139, 0.6593, 0.6219, 0.6920, 0.6668, 0.6209, 0.5096]
        months = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

        # load the coefficient files
        beta0_gwr = rioxarray.open_rasterio(os.path.join(regression_coefficients_path, "Coeffs_GWR_b0.nc"))
        beta1_gwr = rioxarray.open_rasterio(os.path.join(regression_coefficients_path, "Coeffs_GWR_b1.nc"))
        beta2_gwr = rioxarray.open_rasterio(os.path.join(regression_coefficients_path, "Coeffs_GWR_b2.nc"))
        beta0_cswr = rioxarray.open_rasterio(os.path.join(regression_coefficients_path, "Coeffs_CSWR_b0.nc"))
        beta1_cswr = rioxarray.open_rasterio(os.path.join(regression_coefficients_path, "Coeffs_CSWR_b1.nc"))
        beta2_cswr = rioxarray.open_rasterio(os.path.join(regression_coefficients_path, "Coeffs_CSWR_b2.nc"))

        # make a raster stack, to later average computed monthly tair; 
        # iterate over months, and compute monthly tair
                
        air_temperatures = None
        for idx, month in tqauto(enumerate(months), total=len(months)):
            
            month_no = idx + 1
            month_name = month

            # for now, make filename using apparent naming convention; revise as discussed later
            filename_day = os.path.join(local_root, country_name, city_name, "eo_indicators", f"LST_{city_name}_day_{month_name}_{month_no:02d}.tif")
            filename_night = os.path.join(local_root, country_name, city_name, "eo_indicators", f"LST_{city_name}_night_{month_name}_{month_no:02d}.tif")
            
            # import LST and force WGS84 projection
            ds_day = rioxarray.open_rasterio(filename_day).rio.reproject("EPSG:4326")
            ds_night = rioxarray.open_rasterio(filename_night).rio.reproject("EPSG:4326")

            # match coefficients with LST data
            b0 = beta0_gwr[[idx]].rio.clip_box(*ds_day.rio.bounds()).rio.reproject_match(ds_day, resampling=Resampling.nearest)
            b1 = beta1_gwr[[idx]].rio.clip_box(*ds_day.rio.bounds()).rio.reproject_match(ds_day, resampling=Resampling.nearest)
            b2 = beta2_gwr[[idx]].rio.clip_box(*ds_day.rio.bounds()).rio.reproject_match(ds_day, resampling=Resampling.nearest)
            gwr = (b0 + (b1 * ds_night) + (b2 * ds_day)) - 273.15

            b0 = beta0_cswr[[idx]].rio.clip_box(*ds_day.rio.bounds()).rio.reproject_match(ds_day, resampling=Resampling.nearest)
            b1 = beta1_cswr[[idx]].rio.clip_box(*ds_day.rio.bounds()).rio.reproject_match(ds_day, resampling=Resampling.nearest)
            b2 = beta2_cswr[[idx]].rio.clip_box(*ds_day.rio.bounds()).rio.reproject_match(ds_day, resampling=Resampling.nearest) 
            cswr = (b0 + (b1 * ds_night) + (b2 * ds_day)) - 273.15

            current_tair = gwr_weights[idx] * gwr + cswr_weights[idx] * cswr
            if air_temperatures is None:
                air_temperatures = current_tair
            
            else:
                air_temperatures = xr.concat([air_temperatures, current_tair], dim="month")

        ds_mean_air_temperature = air_temperatures.mean(dim="month", skipna=True, keep_attrs=True).transpose('band', 'y', 'x')

        # resample annual mean airT raster to 10 m using population raster, clip to AOI
        lppi_filename = os.path.join(local_root, country_name, city_name, 'eo_indicators', 'lppi.tif')
        ds_lppi = rioxarray.open_rasterio(lppi_filename)

        ds_mean_air_temperature = ds_mean_air_temperature.rio.reproject_match(ds_lppi, resampling=resampling_method)

        # import aoi to clipping/masking
        area_of_interest_filename = os.path.join(local_root, country_name, city_name, 'ancillary', f'{city_name}_AOI.shp')
        shape_aoi = gpd.read_file(area_of_interest_filename)
        shape_aoi = shape_aoi.to_crs(ds_mean_air_temperature.rio.crs)

        # add all_touched=True to ensure all pixels touched by the AOI geometry are included; otherwise, shape may be incorrect
        # this changes the shape of the resulting raster, use drop=False to retain original shape
        ds_mean_air_temperature = ds_mean_air_temperature.rio.clip(shape_aoi['geometry'], all_touched=True, drop=False)
                
        # replace infinity values with nan values
        ds_mean_air_temperature.values[0][np.isinf(ds_mean_air_temperature.values[0])] = np.nan
        # get current nodata value from ds_mean_air_tempature
        current_nodata = ds_mean_air_temperature.rio.nodata
        ds_mean_air_temperature.rio.write_nodata(np.nan, inplace=True)
        ds_mean_air_temperature = ds_mean_air_temperature.where(ds_mean_air_temperature != current_nodata, np.nan)


        # also write result to disk
        ds_out = ds_lppi.copy()
        ds_out.values[0] = ds_mean_air_temperature.values[0]
        
        if out_path is None:
            out_path = os.path.join(local_root, country_name, city_name, HeatModule.module_pathname, 'outputs', 'airT_10m.tif')
        export_raster_dataset(ds_out, out_path)

        return ds_mean_air_temperature

    @staticmethod 
    def interpolate_temperature_raster(
        local_root: str, 
        country_name: str,
        city_name: str,
        temperature_in_path: str = None,        
        temperature_out_path: str = None,
        upscale_factor: float = 0.01,
        interpolation_distance: float = 1000,
        regression_coefficients_path: str = None, 
        resampling_method: Resampling = Resampling.nearest,
    ):
        """Generate and/or interpolate air temperature data. If a temperature_in_path is provided, the raster is read, upscaled, and interpolated over the specified interpolation_distance. If no temperature_in_path is provided, the air temperature raster will be generated.

        Args:
            local_root (str): Local root directory of the scan's data repository.
            country_name (str): Country name (path in local root)
            city_name (str): City name (path in country folder within local root)
            temperature_in_path (str, optional): Path to the input temperature raster. If None, the raster will be generated from LST data. Defaults to None.
            temperature_out_path (str, optional): Path to save the interpolated temperature raster. If None, the default path will automatically be used. Defaults to None.
            upscale_factor (float, optional): Factor by which to upscale the raster. Defaults to 0.01.
            interpolation_distance (float, optional): Distance over which to interpolate missing values. Defaults to 1000.
            regression_coefficients_path (str, optional): Path to the regression coefficients file. Defaults to None.
            resampling_method (Resampling, optional): Resampling method to use. Defaults to Resampling.nearest.
        """

        # if temperature_in_path is None, we will generate the air temperature raster first
        if temperature_in_path is None:
            HeatModule.generate_air_temperature_raster(
                local_root=local_root,
                country_name=country_name,
                city_name=city_name,
                regression_coefficients_path=regression_coefficients_path,
                resampling_method=resampling_method,
                out_path=None
            )

            default_in_path = os.path.join(local_root, country_name, city_name, HeatModule.module_pathname, 'outputs', 'airT_10m.tif')
            ds_tair = rioxarray.open_rasterio(default_in_path)

        else:
            ds_tair = rioxarray.open_rasterio(temperature_in_path)


        hex_path = os.path.join(local_root, country_name, city_name, 'ancillary', f'{city_name}_aoi_hex.tif')
        ds_hex = rioxarray.open_rasterio(hex_path)
        nan_mask = np.isnan(ds_hex.values[0])

        new_width = int(ds_tair.rio.width * upscale_factor)
        new_height = int(ds_tair.rio.height * upscale_factor)
        xds_upsampled = ds_tair.rio.reproject(
            ds_tair.rio.crs,
            shape=(new_height, new_width),
            resampling=Resampling.cubic,
        )

        xds_upsampled = xds_upsampled.rio.pad_box(
            minx=xds_upsampled.rio.bounds()[0] - interpolation_distance,  
            miny=xds_upsampled.rio.bounds()[1] - interpolation_distance,  
            maxx=xds_upsampled.rio.bounds()[2] + interpolation_distance,  
            maxy=xds_upsampled.rio.bounds()[3] + interpolation_distance,  
            constant_values=float("nan")    # fill value for new cells
        )

        ds_tmp = xds_upsampled.rio.interpolate_na(method='nearest')
        ds_tmp = ds_tmp.rio.reproject_match(ds_hex, resampling=Resampling.cubic)
        ds_tmp.values[0][nan_mask] = np.nan

        # close tair dataset
        ds_tair.close()

        # export result to out_path
        if temperature_out_path is None:
            temperature_out_path = os.path.join(local_root, country_name, city_name, HeatModule.module_pathname, 'outputs', 'airT_10m.tif')
        ds_tmp.rio.to_raster(temperature_out_path)


    def make_hexagons_for_baseline(self, use_wri_tree_layer: bool = True):
        # make zonal stats for hexagons using exaxt extract method
        
        # import hexagons
        self.hexagon_layer = gpd.read_file(self.hexagon_filename)
        self.hexagon_layer = self.hexagon_layer.to_crs(self.ds_weighted_population.rio.crs)

        # use exact extract method implemented in gistools
        # extract population sum
        self.hexagon_layer = make_exact_zonal_statistics(self.hexagon_layer, self.ds_weighted_population, {'pop_sum': 'sum'})  
        # extract mean air temperature
        self.hexagon_layer = make_exact_zonal_statistics(self.hexagon_layer, self.ds_mean_air_temperature, {'airT_mean' : 'mean'})
        # extract mean elevation if elevation covariate is included
        if self.elevation_covariate:
            self.hexagon_layer = make_exact_zonal_statistics(self.hexagon_layer, self.ds_dem, {'dem_mean' : 'mean'})


        # compute priority
        self.hexagon_layer['priority'] = self.hexagon_layer['airT_mean'] * self.hexagon_layer['pop_sum']
        # write out priority layer as tif raster
        rasterize_df(self.hexagon_layer, self.baseline_land_cover_filename, os.path.join(self.local_root, self.country_name, self.city_name, 'output', 'heat_demand.tif'), 'priority')

        # import landcover and reproject layer
        self.ds_landcover = rioxarray.open_rasterio(self.baseline_land_cover_filename).rio.reproject(f'EPSG:{self.recommended_crs}', resampling=Resampling.nearest)	
        
        # open tree layer and reproject to match landcover layer
        self.ds_tree_layer = rioxarray.open_rasterio(self.wri_tree_layer_filename).rio.reproject_match(self.ds_landcover, resampling=Resampling.nearest)

        # reproject other relevant layers to common, projected crs to account for area in sq. meters
        self.hexagon_layer = self.hexagon_layer.to_crs(self.recommended_crs)

        # determine validity of tree layer
        # assess whether WRI layer contains any values > 0 other than nan or 0
        # use tree layer only if specified by user
        if use_wri_tree_layer:
            heat_module_logger.info("Assessing validity of WRI Tree layer for use in regression model")
            wri_data = self.ds_tree_layer.values[0]
            tree_layer_relevant_values = wri_data[(np.isin(wri_data, [np.nan], invert=True)) & (wri_data > 0)]
            self.tree_layer_is_valid = any(tree_layer_relevant_values)
        else:
            heat_module_logger.info("WRI Tree layer will not be used in regression model as specified by user")
            self.tree_layer_is_valid = False
        
        
        if self.tree_layer_is_valid:
            heat_module_logger.info("WRI Tree layer is used in assessment")
        else:
            heat_module_logger.warning("WRI Tree layer is not used in assessment")

        # determine share of landuse cover
        landcover_shares = ee.exact_extract(self.ds_landcover, self.hexagon_layer, ["unique", "frac"], include_cols=['id'], output='pandas')
        landcover_shares = landcover_shares.explode(['unique', 'frac'])
        landcover_shares = landcover_shares.pivot(columns=['unique'], values='frac')
        landcover_shares = landcover_shares.fillna(0.0)
        # test a potential fix 
        landcover_shares = landcover_shares.astype(float)

        landcover_shares.columns = [f"class_{int(col)}" for col in landcover_shares.columns]
        landcover_shares['cover_sum'] = landcover_shares.sum(axis=1)

        # overwrite tree cover (class 3) by WRI layer if applicable
        if self.tree_layer_is_valid:
            heat_module_logger.info("Extracting tree coverage values from WRI Tree layer, overwriting class_3")            
            wri_means = ee.exact_extract(self.ds_tree_layer, self.hexagon_layer, "mean", include_cols=['id'], output='pandas')
            landcover_shares['class_3'] = (wri_means['mean'] / 100)

        # make predictors for baseline
        predictor_set = self.make_predictor_set(landcover_shares, nbs_type=None)

        # merge predictor set into hexagon layer
        # at this stage, merging by index should be ok as the order of hexagons is not affected
        self.hexagon_layer = self.hexagon_layer.merge(predictor_set, left_index=True, right_index=True, how='left')

        # add specific id for spatial regression
        self.hexagon_layer['reg_idx'] = self.hexagon_layer.index
        self.hexagon_layer['reg_id'] = self.hexagon_layer['neigh_id']
        
        # filter hexagons already at this stage for inf and nan values
        self.hexagon_layer = self.hexagon_layer.replace(-np.inf, np.nan)
        self.hexagon_layer = self.hexagon_layer.dropna()

        # reset index of hexagon layer following nan filtering
        self.hexagon_layer = self.hexagon_layer.reset_index(drop=True)

        # write to disk as well
        self.hexagon_layer.to_file(os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'outputs', f'hexagons_{self.city_name}.shp'))

    def make_predictor_set(self, landcover_shares: pd.DataFrame, nbs_type: NbsTypes = None) -> pd.DataFrame:

        # build predictor data as dictionary
        data_dict = {}
        
        if RegressionPredictors.BuiltCoverage in self.regression_predictors:
            data_dict[RegressionPredictors.BuiltCoverage.value] = landcover_shares[f'class_{LandCoverClass.Urban.value}']

        if RegressionPredictors.CropCoverage in self.regression_predictors:
            data_dict[RegressionPredictors.CropCoverage.value] = landcover_shares[f'class_{LandCoverClass.Cropland.value}']

        if RegressionPredictors.TreeCoverage in self.regression_predictors:
            data_dict[RegressionPredictors.TreeCoverage.value] = landcover_shares[f'class_{LandCoverClass.Forest.value}']  

        if RegressionPredictors.ShrubCoverage in self.regression_predictors:
            data_dict[RegressionPredictors.ShrubCoverage.value] = landcover_shares[f'class_{LandCoverClass.Shrubs.value}']

        if RegressionPredictors.GrassCoverage in self.regression_predictors:
            data_dict[RegressionPredictors.GrassCoverage.value] = landcover_shares[f'class_{LandCoverClass.Grassland.value}']

        if RegressionPredictors.WaterCoverage in self.regression_predictors:
            data_dict[RegressionPredictors.WaterCoverage.value] = landcover_shares[f'class_{LandCoverClass.Waterbodies.value}']

        if RegressionPredictors.BareCoverage in self.regression_predictors:
            data_dict[RegressionPredictors.BareCoverage.value] = landcover_shares[f'class_{LandCoverClass.Bare.value}']

        if nbs_type is not None:
            # check if the NbS type is represented as column in the df. This might not be the case if suitability is entirely 0 for the NbS in question. 
            if f'class_{nbs_type.get_land_cover_class_value()}' in landcover_shares.columns:
                nbs_specific_tree_cover_ratio = NbsTreeCoverage[nbs_type.name].value
                nbs_tree_cover_share = landcover_shares[f'class_{nbs_type.get_land_cover_class_value()}'] * nbs_specific_tree_cover_ratio
                # add nbs-determined share to overall share (forest/class-3)
                tree_coverage_col = RegressionPredictors.TreeCoverage.value
                data_dict[tree_coverage_col] = data_dict[tree_coverage_col] + nbs_tree_cover_share
        


        # create DataFrame from dictionary with hexagon_layer index
        dest_df = pd.DataFrame(data_dict, index=self.hexagon_layer.index)

        return dest_df

    def make_regression(self):
        
        # make regression using Spatial Durbin Error Model on the baseline hexagon layer
        heat_module_logger.info("Generating spatial weights")
        

        # define names of (in)dependent variables in data frame
        self.y_name = 'airT_mean'
        self.x_names = [p.value for p in self.regression_predictors]  # by default, built/urban and tree cover (either from landcover, or from WRI layer)
        if self.elevation_covariate:
            self.x_names.append('dem_mean')

        w_name = 'queen_weights'
        ds_name = 'hexagon_layer'

        # create spatial weights
        self.regression_weights = weights.Queen.from_dataframe(self.hexagon_layer, use_index=False, ids='reg_id') 
        self.regression_weights.transform = 'r' # row transform weights

        self.hexagon_layer.to_file(os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'outputs', 'reg_data.shp'))    

        y = self.hexagon_layer[self.y_name]
        y = y.astype(float)
        
        x = self.hexagon_layer[self.x_names]
        x = sm.add_constant(x)        
        x = x.astype(float)
        
        # this is the spatial durbin error model producing results in line with R implementation
        heat_module_logger.info("Computing Spatial Durbin Error regression model")        
        gns1 = ML_Error(y, x, w=self.regression_weights, slx_lags = 1, add_wy = False, name_w=w_name, name_ds=ds_name, method='full')

        # print model summary and redirect to file as well for later inspection
        # this is mainly meant as a means to document what has been done
        print(gns1.summary)
        with fopen(os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'outputs', 'regression_results.txt'), 'w') as f:
            f.write(str(gns1.summary))

        # get relevant coefficients to determine model validity
        self.regression_betas = gns1.betas

        # iterate over predictors included in the model
        # when we iterate, betas should be 
        # [0] = constant
        # [idx+1] ... [n+1] = predictors direct effects  
        # [idx+1+n] ... [idx+1+n] = predictors indirecteffects

        predictor_count = len(self.regression_predictors)
        # if elevetation is included as covariate, this will be the last predictor, so we need to add an additional offset for indirect effect
        if self.elevation_covariate:
            predictor_count += 1
        
        for idx, predictor in enumerate(self.regression_predictors):
            
            idx_direct_effect = idx + 1 # offset due to 0 being constant
            idx_indirect_effect = idx_direct_effect + predictor_count # offset due to direct effects
            
            # print out the determined marginal effects again for confirmation purposes
            heat_module_logger.info(f"Validation of {predictor.value}: direct={self.regression_betas[idx_direct_effect]}, indirect={self.regression_betas[idx_indirect_effect]}")

            # check validity of built coverage
            if predictor == RegressionPredictors.BuiltCoverage:
                if self.regression_betas[idx_direct_effect] < 0:
                    heat_module_logger.warning("Direct marginal effect for built coverage is negative. Setting to 0.001")
                    self.regression_betas[idx_direct_effect] = 0.001
                else:
                    heat_module_logger.info("Direct marginal effect for built coverage is positive in regression model - No correction applied.")
                
                if self.regression_betas[idx_indirect_effect] < 0:
                    heat_module_logger.warning("Indirect marginal effect for built coverage is negative. Setting to 0.0001")
                    self.regression_betas[idx_indirect_effect] = 0.0001
                else:
                    heat_module_logger.info("Indirect marginal effect for built coverage is positive in regression model - No correction applied.")
            
            # check validity of tree coverage
            if predictor == RegressionPredictors.TreeCoverage:
                if self.regression_betas[idx_direct_effect] > 0:
                    heat_module_logger.warning("Direct marginal effect for tree coverage is positive. Setting to -0.001")
                    self.regression_betas[idx_direct_effect] = -0.001
                else:
                    heat_module_logger.info("Direct marginal effect for tree coverage is negative in regression model - No correction applied.")
                
                if self.regression_betas[idx_indirect_effect] > 0:
                    heat_module_logger.warning("Indirect marginal effect for tree coverage is positive. Setting to -0.0001")
                    self.regression_betas[idx_indirect_effect] = -0.0001
                else:
                    heat_module_logger.info("Indirect marginal effect for tree coverage is negative in regression model - No correction applied.")

            # also document the coefficients used for predictions
            with fopen(os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'outputs', 'corrected_regression_betas.txt'), 'w') as f:
                f.write(str(self.regression_betas))




    def asess_scenarios(self, nbs_list: List[NbsTypes]):

        # create lagged variables on baseline conditions: subsequently, we can overwrite non-lagged variables with scenario-specific, updated values 
        lagged_predictors = (
            self.hexagon_layer[self.x_names]
            .apply(
                # make spatially lagged variables
                lambda y: weights.spatial_lag.lag_spatial(self.regression_weights, y)
            )
            .rename(
                # Rename the spatial lag, adding w_ to the original name
                columns=lambda c: "w_" + c
            )
        )

        # this makes data for baseline to apply regression model for estimation of current
        nbs_prediction_data = self.hexagon_layer[self.x_names]

        # merge using key, to assert correct mapping
        # insert key to lagged predictors; this dataset corresponds to baseline and should not be affected by NBS modelling
        lagged_predictors['reg_id_r'] = self.hexagon_layer['reg_id']
        lagged_predictors.to_excel(os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'outputs', f'lagged_predictors.xlsx'))
        
        nbs_prediction_data['reg_id_r'] = self.hexagon_layer['reg_id']
        nbs_prediction_data = nbs_prediction_data.merge(lagged_predictors, on='reg_id_r', how='left')

        # remove merge col
        nbs_prediction_data = nbs_prediction_data.drop(columns=['reg_id_r'])
        nbs_prediction_data.to_excel(os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'outputs', f'lagged_predictors_merged.xlsx'), index=False)

        self.hexagon_layer['y_current'] = self.regression_betas[0] + nbs_prediction_data @ self.regression_betas[1:-1]
        self.hexagon_layer['e_current'] = self.hexagon_layer['airT_mean'] - self.hexagon_layer['y_current']

        # iterate over NBS types and apply changes to the hexagon layer
        for current_nbs in tqauto(nbs_list, total=len(nbs_list)):

            nbs_y_fieldname = f'y_{current_nbs.get_name()}'

            # load corresponding creation layer
            scenario_filename = os.path.join(self.local_root,self.country_name,self.city_name, 'suitability', f"{current_nbs.get_name()}_creation.tif")
            ds_scenario = rioxarray.open_rasterio(scenario_filename).rio.reproject_match(self.ds_landcover, resampling=Resampling.nearest)
            
            # assert that scenario creation layer has opportunities for creation.
            if np.nansum(ds_scenario.data[0]) > 0:

                # burn nbs creation into landcover layer
                scenario_landcover = self.ds_landcover.copy()                            
                scenario_landcover.data[0][ds_scenario.data[0] == 1] = current_nbs.get_land_cover_class_value()

                # export for review
                out_path = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'outputs', f"landcover_{current_nbs.get_name()}.tif")
                export_raster_dataset(scenario_landcover, out_path=out_path)

                # use exact_extract to determine the share of the NBS in the hexagons
                scenario_landcover_shares = ee.exact_extract(scenario_landcover, self.hexagon_layer, ["unique", "frac"], include_cols=['id'], output='pandas')
                scenario_landcover_shares = scenario_landcover_shares.explode(['unique', 'frac'])
                scenario_landcover_shares = scenario_landcover_shares.pivot(columns=['unique'], values='frac')
                scenario_landcover_shares = scenario_landcover_shares.fillna(0.0)
                scenario_landcover_shares.columns = [f"class_{int(col)}" for col in scenario_landcover_shares.columns]

                non_lagged_predictors = self.make_predictor_set(scenario_landcover_shares, nbs_type=current_nbs) # decide whether we want to store those to the hexagon layer or not
                non_lagged_predictors['reg_id_r'] = self.hexagon_layer['reg_id']
                if self.elevation_covariate:
                    non_lagged_predictors['dem_mean'] = self.hexagon_layer['dem_mean']

                nbs_prediction_data = non_lagged_predictors.merge(lagged_predictors, on='reg_id_r', how='left')
                nbs_prediction_data = nbs_prediction_data.drop(columns=['reg_id_r'])

                # sanity checking if somewhere, tree cover is now > 100 (this should not happen)
                max_values = nbs_prediction_data[[pred.value for pred in self.regression_predictors]].max()
                if max(max_values) > 1.00:
                    heat_module_logger.error(f"Coverage of some predictor classes exceeds 100% in some hexagons: {max_values}")

                # apply regression model to estimate air temperature for current nbs
                self.hexagon_layer[nbs_y_fieldname] = self.regression_betas[0] + nbs_prediction_data @ self.regression_betas[1:-1]
                # determine impact in terms of difference to baseline: nbs - baseline
                self.hexagon_layer[f'{current_nbs.get_name()}_temp_diff'] = self.hexagon_layer[nbs_y_fieldname] - self.hexagon_layer['y_current']

            else:

                heat_module_logger.warning(f"Scenario {current_nbs.get_name()} has no creation opportunities in the AOI.")
                self.hexagon_layer[f'{current_nbs.get_name()}_temp_diff'] = 0


        
    def determine_benefits(self, nbs_list: List[NbsTypes]):
        
        # get hexagon layer
        ds_hex = rioxarray.open_rasterio(os.path.join(self.local_root, self.country_name, self.city_name, 'ancillary', f'{self.city_name}_aoi_hex.tif')) 
        
        # determine impact for each nbs
        heat_module_logger.info("Determining and exporting heat benefits")
        for current_nbs in tqauto(nbs_list, total=len(nbs_list)):
            
            # impact is equal to -1 * temperature * population * t_diff
            nbs_benefit_fieldname = f'{current_nbs.get_name()}_heat_benefits'
            self.hexagon_layer[nbs_benefit_fieldname] = -1 * self.hexagon_layer[self.y_name] * self.hexagon_layer['pop_sum'] * self.hexagon_layer[f'{current_nbs.get_name()}_temp_diff']

            # write hex
            out_layer = self.hexagon_layer.copy()
            out_layer.to_crs(ds_hex.rio.crs, inplace=True)
            rasterize_df(out_layer, self.baseline_land_cover_filename, os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'outputs', 'heat_benefits', f'{current_nbs.get_name()}_heat_benefits.tif'), nbs_benefit_fieldname, fill_value=0, nodata_value=np.nan)
                        
        # write hexagon layer as a reference
        heat_module_logger.info("Exporting hexagon layer to geopackage")
        self.hexagon_layer.to_file(os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'outputs', f'predictions_{self.city_name}.gpkg'), driver="GPKG", layer='predictions')


    def run_module(self,  
            nbs_list: List[NbsTypes],          
            kernel_radius: float,
            kernel_sigma: float,
            temperature_resampling_method: Resampling,
            regression_coefficients_path: str,
            predictors: List[RegressionPredictors] = [RegressionPredictors.BuiltCoverage, RegressionPredictors.TreeCoverage],
            elevation_covariate: bool = False,
            use_wri_tree_layer: bool = True,
            air_temperature_path: str = None,
            skip_air_temperature_computation: bool = False,
        ):

        # store internal parameters
        self.regression_predictors = predictors

        # store elevation covariate parameter for later use in regression and scenario assessment
        self.elevation_covariate = elevation_covariate

        # estimate weighted population
        self.estimate_weighted_population(
            kernel_radius=kernel_radius,
            kernel_sigma=kernel_sigma
        )

        # estimate mean air temperature, if no dataset specified by user
        self.air_temperature_path = air_temperature_path

        if air_temperature_path is None and not skip_air_temperature_computation:
            heat_module_logger.info("Generating air temperature raster from LST data")
            self.estimate_mean_air_temperature(
                regression_coefficients_path=regression_coefficients_path,
                resampling_method=temperature_resampling_method
            )
        else:
            
            # import existing datasets, as specified
            if skip_air_temperature_computation:
                heat_module_logger.info("Skipping air temperature computation as specified by user - importing from default location instead")
                default_path = os.path.join(self.local_root, self.country_name, self.city_name, HeatModule.module_pathname, 'outputs', 'airT_10m.tif')
                self.ds_mean_air_temperature = rioxarray.open_rasterio(default_path)
            elif air_temperature_path is not None:
                heat_module_logger.info(f"Importing air temperature raster from user-specified path: {air_temperature_path}")
                self.ds_mean_air_temperature = rioxarray.open_rasterio(air_temperature_path)

        if self.elevation_covariate:
            self.ds_dem = rioxarray.open_rasterio(self.dem_filename).rio.reproject_match(self.ds_mean_air_temperature, resampling=Resampling.bilinear)

        # make hexagons for baseline and subsequent analysis
        self.make_hexagons_for_baseline(use_wri_tree_layer=use_wri_tree_layer)

        # make regression
        self.make_regression()

        # assess scenarios: predict air temperature for different scenarios
        self.asess_scenarios(nbs_list=nbs_list)	

        # determine benefits
        self.determine_benefits(nbs_list=nbs_list)

    @staticmethod
    def get_module_pathname() -> str:
        return "Heat"

    @staticmethod
    def run(
        local_root: str,
        country: str,
        city: str,
        nbs_list: List[NbsTypes],
        kernel_radius: float,
        kernel_sigma: float,
        temperature_resampling_method: Resampling,
        regression_coefficients_path: str,
        predictors: List[RegressionPredictors] = [RegressionPredictors.BuiltCoverage, RegressionPredictors.TreeCoverage],
        elevation_covariate: bool = False,
        use_wri_tree_layer: bool = True,
        air_temperature_path: str = None,
        skip_air_temperature_computation: bool = False,
        ):

        module_instance = HeatModule(
            local_root=local_root, 
            country=country, 
            city=city
        )

        module_instance.run_module(
            nbs_list=nbs_list,
            kernel_radius=kernel_radius,
            kernel_sigma=kernel_sigma,
            temperature_resampling_method=temperature_resampling_method,
            regression_coefficients_path=regression_coefficients_path,
            predictors=predictors,
            elevation_covariate=elevation_covariate,
            use_wri_tree_layer=use_wri_tree_layer,
            air_temperature_path=air_temperature_path,
            skip_air_temperature_computation=skip_air_temperature_computation,
        )

# expose static methods at module level
generate_air_temperature_raster = HeatModule.generate_air_temperature_raster
interpolate_air_temperature_raster = HeatModule.interpolate_temperature_raster