import rioxarray as riox
import os
import numpy as np
import rasterio
from typing import List
import rioxarray
from rasterio.warp import calculate_default_transform, reproject, Resampling
from ..nbsos_exports import NbsTypes, LandCoverClass
from ..nbsos_exports import LandCoverClass

from .erosion_exports import Nbs_Erosion_Weights, Lulc_Erosion_Weights

from ..gis_tools import export_raster, make_exact_zonal_statistics, rasterize_df, export_raster_dataset, make_gaussian_kernel, get_area_crs
from scipy.ndimage import convolve
from tqdm.auto import tqdm as tqauto
import geopandas as gpd
import exactextract as ee

class ErosionModule:

    nbs_list: List[NbsTypes] = None

    def __init__(self, local_root: str, country_name: str, city: str, kernel_radius: float = 1000, kernel_sigma: float = 1000/3, use_population_weighting: bool = True):
        

        self.module_pathname = ErosionModule.get_module_pathname()
        
        self.local_root = local_root
        self.country_name = country_name
        self.city_name = city
        self.use_population_weighting = use_population_weighting

        # parameters for population Gauss weights 
        self.kernel_radius = kernel_radius
        self.kernel_sigma = kernel_sigma

        # set some default filenames
        self.population_filename = os.path.join(self.local_root, self.country_name, self.city_name, 'eo_indicators', 'hrsl_pop.tif')
        self.hexagon_filename = os.path.join(self.local_root, self.country_name, self.city_name, 'ancillary', f'{self.city_name}_hex.shp')
        self.area_of_interest_filename = os.path.join(self.local_root, self.country_name, self.city_name, 'ancillary', f'{self.city_name}_AOI.shp')

        self.recommended_crs = get_area_crs(self.area_of_interest_filename)


        required_dirs = [
            'output',
            os.path.join('output', 'erosion_benefits')            
        ]
        for f in required_dirs:
            self.make_folder(f)





    def make_baseline_weight_maps(self):
                
        land_cover_mtx = self.ds_land_cover.values[0]
        self.baseline_weight_map = np.full_like(land_cover_mtx, fill_value=1, dtype=np.float32)

        for lulc in LandCoverClass:
            
            # current name and lu class
            current_land_cover_name = lulc.name
            current_land_cover_class = lulc.value
            
            # Get corresponding weight and make baseline weight map
            current_weight = Lulc_Erosion_Weights[current_land_cover_name].value
            self.baseline_weight_map[land_cover_mtx == current_land_cover_class] = current_weight

        # write out priority layer as tif raster
        out_path = os.path.join(self.local_root, self.country_name, self.city_name, self.get_module_pathname(), 'output', 'baseline_weights.tif')
        export_raster(self.baseline_weight_map, self.rusle_k_path, out_path)

    def estimate_weighted_population(self, kernel_radius: float, kernel_sigma: float):

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
        out_path = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, "pop_weight.tif")
        export_raster_dataset(self.ds_weighted_population, out_path)
        self.pop_weights = rioxarray.open_rasterio(out_path).values[0]





    def run_module(self, nbs_list: List[NbsTypes]):

        self.nbs_list = nbs_list

        # import land cover, which we use as a reference as well later
        land_cover_path = os.path.join(self.local_root, self.country_name, self.city_name, 'eo_indicators', 'land_cover.tif')
        self.ds_land_cover = rioxarray.open_rasterio(land_cover_path)

        self.rusle_k_path = os.path.join(self.local_root, self.country_name, self.city_name, 'eo_indicators', 'rusle_k.tif')
        with rasterio.open(self.rusle_k_path) as src:
            profile = src.meta
            rusle_k_mtx = src.read(1)

        rusle_r_path = os.path.join(self.local_root, self.country_name, self.city_name, 'eo_indicators', 'rusle_r.tif')
        self.rusle_r_mtx = riox.open_rasterio(rusle_r_path).values[0]
        rusle_ls_path = os.path.join(self.local_root, self.country_name, self.city_name, 'eo_indicators', 'rusle_ls.tif')
        rusle_ls_mtx = riox.open_rasterio(rusle_ls_path).values[0]

        # compute rkls raster
        self.rusle_k_ls_mtx = rusle_k_mtx * rusle_ls_mtx
        out_path = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'output', 'k_ls.tif')
        export_raster(self.rusle_k_ls_mtx, self.rusle_k_path, out_path)
        self.rkls_mtx = self.rusle_r_mtx * self.rusle_k_ls_mtx # this is potential maximum erosion, as it is not considering vegetation effects

        # make population
        self.estimate_weighted_population(kernel_radius=self.kernel_radius, kernel_sigma=self.kernel_sigma)

        # make c map for baseline
        self.make_baseline_weight_maps()

        # avoided erosion for baseline and assumed erosion rate
        # mask rusle r
        self.rusle_r_mtx[np.isnan(self.pop_weights)] = np.nan
        self.baseline_erosion_rate = np.full_like(self.rusle_r_mtx, fill_value=np.nan, dtype=np.float32)
        self.baseline_erosion_rate = np.divide(self.rusle_r_mtx, self.baseline_weight_map, out=self.baseline_erosion_rate, where=self.baseline_weight_map!=0) 
        self.baseline_erosion_rate = self.baseline_erosion_rate * self.rusle_k_ls_mtx # this is actual erosion rate, considering vegetation benefits

        # determine avoided erosion as difference of rkls and r/v * k ls 
        avoided_erosion_baseline = self.rkls_mtx - self.baseline_erosion_rate
        
        out_profile = profile.copy()
        out_profile['nodata'] = np.nan        

        # export avoided erosion for baseline 
        out_path = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'output', 'baseline_avoided_erosion.tif')
        export_raster(avoided_erosion_baseline, self.rusle_k_path, out_path, out_profile)
        
        # export baseline actual erosion rate 
        out_path = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'output', 'baseline_erosion_rate.tif')
        export_raster(self.baseline_erosion_rate, self.rusle_k_path, out_path, out_profile)
        
        # erosion priority is erosion rate baseline * population
        baseline_priority = self.baseline_erosion_rate
        if self.use_population_weighting:
            baseline_priority = baseline_priority * self.pop_weights
            
        priority_out_path = os.path.join(self.local_root, self.country_name, self.city_name, 'output', 'priority', 'erosion_demand.tif')
        export_raster(baseline_priority, self.rusle_k_path, priority_out_path, out_profile)

        if self.nbs_list is not None:
            for nbs in tqauto(self.nbs_list):
                tqauto.write(f'Processing erosion for NbS: {nbs.get_name()}')
                self.run_erosion_model_for_nbs(nbs, out_profile)
    

    def run_erosion_model_for_nbs(self, nbs: NbsTypes, out_profile):

        # create a copy of the baseline weight map
        nbs_weights_map = np.zeros_like(self.baseline_weight_map, dtype=np.float32)

        # get the corresponding NBS erosion weights
        nbs_weights_factor = Nbs_Erosion_Weights[nbs.name].value

        # import suitability for current nbs
        suitability_path = os.path.join(self.local_root,self.country_name,self.city_name, 'suitability', f"{nbs.get_name()}_creation.tif")
        suitability_mtx = riox.open_rasterio(suitability_path).values[0]
        
        # update the weight map for the areas where the land cover matches the NbS land cover class
        nbs_weights_map[suitability_mtx == 1] = nbs_weights_factor
        out_path = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'output', 'erosion_benefits', f'{nbs.get_name()}_weights.tif')
        export_raster(nbs_weights_map, self.rusle_k_path, out_path, out_profile)

        # compute erosion benefits
        nbs_erosion_benefits_mtx = self.baseline_erosion_rate * nbs_weights_map

        if self.use_population_weighting:
            nbs_erosion_benefits_mtx = nbs_erosion_benefits_mtx * self.pop_weights

        # export rasters
        out_path = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'output', 'erosion_benefits', f'diff_{nbs.get_name()}_erosion.tif')
        export_raster(nbs_erosion_benefits_mtx, self.rusle_k_path, out_path, out_profile)
        
    def make_folder(self, folder_name):
        if folder_name is not None:
            target_folder = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, folder_name)
        else:
            target_folder = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname)
        os.makedirs(target_folder, exist_ok=True) 

    @staticmethod
    def get_module_pathname() -> str:
        return "Erosion"
    
    @staticmethod
    def run(
        local_root: str,
        country_name: str,
        city: str,
        nbs_list: List[NbsTypes],
        use_population_weighting: bool,
        kernel_radius: float = 1000,
        kernel_sigma: float = 1000/3,
        ):

        module_instance = ErosionModule(
            local_root=local_root, 
            country_name=country_name, 
            city=city,
            use_population_weighting=use_population_weighting,
            kernel_radius=kernel_radius,
            kernel_sigma=kernel_sigma
        )

        module_instance.run_module(
            nbs_list=nbs_list           
        )

    @staticmethod
    def get_erosion_protection_benefits(local_root, country_name, city_name, protection_mask_type: str, use_population_weighting: bool = True) -> None:
        """
        Get the benefits of erosion protection. This is based on avoided erosion for baseline, i.e., valuing the benefits of existing greenery in 
        the avoidance of soil loss.
        """
        
        # import protection opportunities        
        if protection_mask_type == 'suitability':
            protection_opportunity_filename = os.path.join(local_root, country_name, city_name, 'suitability', 'nbs_protection.tif')
        elif protection_mask_type == 'lppi':
            protection_opportunity_filename = os.path.join(local_root, country_name, city_name, 'suitability', 'nbs_protection_clustering.tif')
            
        ds_protection_opportunities = rioxarray.open_rasterio(protection_opportunity_filename)
        mtx_prot = ds_protection_opportunities.values[0]
        
        # reclasssify as needed
        if protection_mask_type == "suitability":
            mtx_prot = np.where(mtx_prot > 0, 1.0, 0)  # make sure that we have the correct value for protection
        elif protection_mask_type == "lppi":
            mtx_prot = np.where(mtx_prot == 1, 1.0, 0)  # make sure that we have the correct value for protection
        ds_protection_opportunities.values[0] = mtx_prot
        
        baseline_protection_indicator = os.path.join(local_root, country_name, city_name, ErosionModule.get_module_pathname(), 'output', 'baseline_erosion_rate.tif') # this is actual erosion for baseline in line with the manuscript
        ds_erosion_protection_baseline = rioxarray.open_rasterio(baseline_protection_indicator)

        # match projections to protection datasource
        ds_erosion_protection_baseline = ds_erosion_protection_baseline.rio.reproject_match(ds_protection_opportunities, resampling=Resampling.nearest, align_on_exact=True)

        # compute protection benefits
        ds_erosion_rate_protection_benefits = ds_erosion_protection_baseline * ds_protection_opportunities
        
        if use_population_weighting:
            print("Applying population weighting")
            pop_path = os.path.join(local_root, country_name, city_name, ErosionModule.get_module_pathname(), "pop_weight.tif") 
            with rasterio.open(pop_path) as src:
                mtx_pop = src.read(1)

            ds_erosion_rate_protection_benefits.values[0] = ds_erosion_rate_protection_benefits.values[0] * mtx_pop

        # export rasters
        out_path_erosion_rate_protection_benefits = os.path.join(local_root, country_name, city_name, 'output', 'protection', 'nbs_protection_erosion_rate_benefits.tif')
        ds_erosion_rate_protection_benefits.rio.to_raster(out_path_erosion_rate_protection_benefits, driver='GTiff', dtype='float32')