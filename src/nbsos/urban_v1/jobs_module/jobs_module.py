from .jobs_exports import Jobs_Nbs_Weights

from ..nbsos_exports import NbsTypes, LandCoverClass
from ..nbsos_exports import LandCoverClass
from ..scan_context import ScanContext
from ..gis_tools import export_raster_dataset, export_raster, make_exact_zonal_statistics, rasterize_df, make_gaussian_kernel, get_area_crs

from scipy.ndimage import convolve
import rioxarray as riox
import os
import numpy as np
import rasterio
from typing import List
import rioxarray
from rasterio.warp import calculate_default_transform, reproject, Resampling
from tqdm.auto import tqdm as tqauto
import geopandas as gpd
import exactextract as ee

class JobCreationModule:

    nbs_list: List[NbsTypes] = None

    def __init__(self, scan_context: ScanContext):
        
        self.local_root = scan_context.local_root
        self.country_name = scan_context.country
        self.city_name = scan_context.city

        target_folder = os.path.join(self.local_root, self.country_name, self.city_name, self.get_module_pathname())
        os.makedirs(target_folder, exist_ok=True) 

        # determine recommended crs
        self.area_of_interest_filename = os.path.join(self.local_root, self.country_name, self.city_name, 'ancillary', f'{self.city_name}_aoi.shp')
        self.recommended_crs = get_area_crs(self.area_of_interest_filename)



    def run_module(
            self, 
            nbs_list: List[NbsTypes],
            kernel_radius: float,
            kernel_sigma: float,
        ):
        
        self.nbs_list = nbs_list
        self.kernel_radius = kernel_radius
        self.kernel_sigma = kernel_sigma
        
        print("Running Jobs Creation module for NbS types: ", [nbs.get_name() for nbs in self.nbs_list])

        # make weighted population
        self.estimate_weighted_population(kernel_radius=1000, kernel_sigma=1000/3)

        # run job creation benefits for each nbs
        if self.nbs_list is not None:
            for nbs in tqauto(self.nbs_list):
                tqauto.write(f'Estimating job creation benefits for NbS: {nbs.get_name()}')
                self.estimate_job_creation_benefits_for_nbs(nbs)





    def estimate_weighted_population(self, kernel_radius: float, kernel_sigma: float):
        """Gaussian kernel of weighted population"""
        
        hrsl_path = os.path.join(self.local_root, self.country_name, self.city_name, "eo_indicators", "hrsl_pop.tif")
        self.ds_population = rioxarray.open_rasterio(hrsl_path, masked=True)
        
        # mask population outside the AOI 
        ds_aoi_hex = rioxarray.open_rasterio(os.path.join(self.local_root, self.country_name, self.city_name, 'ancillary', f'{self.city_name}_aoi_hex.tif'))
        self.aoi_mask = np.isin(ds_aoi_hex.data[0], [0], invert=False)
        self.ds_population.values[0][self.aoi_mask] = np.nan

        # aggregate: by factor of 3
        scale_factor = 3
        population_source_resolution = self.ds_population.rio.resolution()
        new_resolution = (population_source_resolution[0] * scale_factor)
        
        population_downsampled = self.ds_population.rio.reproject(
            dst_crs=self.ds_population.rio.crs,
            resolution=new_resolution,
            resampling=Resampling.average,
        ).rio.reproject(f'EPSG:{self.recommended_crs}')

        # create gaussian kernel on resampled raster
        # scale by max value
        downscaled_resolution = population_downsampled.rio.resolution()
        population_kernel = make_gaussian_kernel(d=kernel_radius, sigma=kernel_sigma, rs=downscaled_resolution)
        population_kernel = population_kernel/np.nanmax(population_kernel) 
        
        # conduct convolution
        xds_downscaled_mtx = np.nan_to_num(population_downsampled.to_numpy())
        xds_downscaled_mtx = convolve(xds_downscaled_mtx[0], weights=population_kernel, mode='constant', cval=0)
        # replace old data in dataset with modified numpy array
        population_downsampled.data = xds_downscaled_mtx.reshape(1, xds_downscaled_mtx.shape[0], xds_downscaled_mtx.shape[1])

        # reproject to original resolution and export
        self.ds_weighted_population = population_downsampled.rio.reproject_match(self.ds_population, resampling=Resampling.nearest)
        out_path = os.path.join(self.local_root, self.country_name, self.city_name, self.get_module_pathname(), "weighted_population.tif")
        export_raster_dataset(self.ds_weighted_population, out_path)


    def estimate_job_creation_benefits_for_nbs(self, nbs: NbsTypes):

        # import suitability for current nbs
        suitability_path = os.path.join(self.local_root,self.country_name,self.city_name, 'suitability', f"{nbs.get_name()}_creation.tif")
        suitability_mtx = rioxarray.open_rasterio(suitability_path).values[0]
        
        # get weighting factor
        nbs_weights_factor = Jobs_Nbs_Weights[nbs.name].value

        # jobs benefits as suitability * weighted population * nbs weight factor
        jobs_benefits_mtx = suitability_mtx * self.ds_weighted_population.values[0] * nbs_weights_factor
        ds_jobs_benefits = self.ds_population.copy()
        ds_jobs_benefits.values[0] = jobs_benefits_mtx

        # export jobs benefits raster
        out_path = os.path.join(self.local_root, self.country_name, self.city_name, self.get_module_pathname(), f"{nbs.get_name()}_jobs_benefits.tif")
        export_raster_dataset(ds_jobs_benefits, out_path)



    @staticmethod
    def get_module_pathname() -> str:
        return "JobsCreation"
    
    @staticmethod
    def run(
            scan_context: ScanContext,
            nbs_list: List[NbsTypes],
            kernel_radius: float = 1000, 
            kernel_sigma: float = 1000/3):

        module_instance = JobCreationModule(
            scan_context=scan_context
        )

        module_instance.run_module(
            nbs_list=nbs_list,
            kernel_radius=kernel_radius,
            kernel_sigma=kernel_sigma,           
        )