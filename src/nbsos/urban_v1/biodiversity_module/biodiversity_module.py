##########################################################################################################
# Nature-Based Solutions Opportunity Scan: Biodiversity Module
# Author: Sebastian Scheuer 
##########################################################################################################

# module-specific imports
from .biodiversity_exports import LandCover, OpportunityType, Metrics, LandCoverClassVegetationProperties, HabitatThreat, PhysiognomicLandscapeType, Metrics
from .biodiversity_results import MetricsContainer, RasterContainer, CreationOpportunityImpactsContainer, RasterSet, MetricsSet, CreationOpportunityImpactsSet, ResultCollection
from .biodiversity_assessment import AssessmentContainer
from .biodiversity_spatial import Spatial, PatchContiguity
from .biodiversity_indicators import BiodiversityIndicators
from .biodiversity_logging import biodiversity_logger

from ..gis_tools import export_raster

# general imports
from contextlib import contextmanager 
from tqdm import tqdm
from tqdm.auto import tqdm as tqauto
from typing import List, Tuple, Dict
from enum import Enum
import rasterio
from rasterio.io import MemoryFile
import rasterio.features
import rasterio.mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from pathlib import Path
import rioxarray

import numpy as np
import numpy.ma as ma

from scipy import ndimage

from skimage import util

import pandas as pd
import geopandas as gpd

import functools

import os
import logging
import warnings
warnings.filterwarnings("ignore")

def log_function(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as exc:
            print(args, kwargs, repr(exc))

    return wrapper

class BiodiversityModule:

    # data directory name. This is used for writing exports
    module_pathname: str = "Biodiversity"

    # variables to retrieve data. Those need to be passed to the module at a minimum
    local_root: str = None
    country_name: str = None
    city_name: str = None

    # internal handling of data and model specifications
    # parameters needed in the assessment
    recommended_epsg: int = None
    resolution: Tuple[float,float] = None
    nodata_values: List[int] = None
    connectivity_threshold: float = 100.0 # in meters

    # AOI, hexagons, roads
    area_of_interest_filename: str = None
    hexagon_filename: str = None
    mtx_hexagons: np.ndarray = None    
    roads_filename: str = None
    geodf_roads: gpd.GeoDataFrame = None

    # land cover (also used as reference layer)
    baseline_habitat_types: List[LandCover] = None
    landscape_types_to_assess: List[PhysiognomicLandscapeType] = None
    baseline_land_cover_filename: str = None
    land_cover_dataset_reader = None
    projected_raster_metadata = None
    mtx_clipped_land_cover: np.ndarray = None
    
    src_dataset = None
    src_raster_metadata = None 
    src_raster_transform = None
    src_raster_bounds = None
    src_raster_crs = None
        
    # habitat quality variables
    habitat_threats: List[HabitatThreat] = None
    habitat_sensitivity = None

    # overall results will be kept track in sort of dict.
    assessment_results: List[ResultCollection] = []
    nbs_scenarios_for_assessment = None
    optimization_criterion: Metrics = None

    df_impact_comparison: Dict[PhysiognomicLandscapeType, pd.DataFrame] = None

    def __init__(self, local_root: str, country: str, city: str):

        self.module_pathname = BiodiversityModule.get_module_pathname()

        # store paths and location names to retrieve data from corresponding folder hierarchy
        self.local_root = local_root
        self.country_name = country
        self.city_name = city

        # results of impact assessment will be stored in a dict
        self.df_impact_comparison = {}
        
        # prepare logging: create module pathname
        self.make_folder(None)

        required_dirs = [
            'input_data',
            'outputs'
        ]
        for f in required_dirs:
            self.make_folder(f)


        # log_path = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname)
        # logging.basicConfig(filename=f'{log_path}/{self.city_name}_biodiversity_scan.log', format='%(asctime)s - %(levelname)s - %(message)s', datefmt="%Y-%m-%d %H:%M:%S", filemode='w', encoding='utf-8', level=logging.INFO)
        # logging_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
        # logging_stream_handler = logging.StreamHandler()
        # logging_stream_handler.setFormatter(logging_formatter)        
        # biodiversity_logger.addHandler(logging_stream_handler)

        # # get started with prepring input data
        # biodiversity_logger.info(f"Initializing biodiversity scan for {self.city_name}, {self.country_name}")
    
    def make_folder(self, folder_name):
        if folder_name is not None:
            target_folder = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, folder_name)
        else:
            target_folder = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname)
        os.makedirs(target_folder, exist_ok=True)


    def set_environment(self, 
                        baseline_habitat_types: List[LandCover] = None, 
                        landscape_types: List[PhysiognomicLandscapeType] = [PhysiognomicLandscapeType.Any], 
                        sensitivity = None, 
                        threats = None, 
                        nodata_values: List[int] = [0], 
                        connectivity_threshold: float = 100,
                        road_types = ['motorway', 'primary', 'construction', 'secondary', 'tertiary']) -> None:
        
        # use default set of habitat types (forest, shrubs, grassland, wetlands) if not explicitly provided
        if baseline_habitat_types is None:
            baseline_habitat_types = [LandCover.FOREST, LandCover.SHRUBS, LandCover.GRASSLAND, LandCover.WETLANDS]
        self.baseline_habitat_types = baseline_habitat_types

        # physiognomic landscape types to assess. by default, assess "Any", i.e., do not differentiate potential habitats and therefore, assess overall connectivity
        self.landscape_types_to_assess = landscape_types  
        self.road_types = road_types

        # use default threats and sensitivities, if not explicitly provided
        if sensitivity is None:
            sensitivity = [
                LandCover.URBAN,
                LandCover.CROPLAND,
                LandCover.FOREST,
                LandCover.SHRUBS,
                LandCover.GRASSLAND,
                LandCover.WATERBODIES,
                LandCover.WETLANDS,
                LandCover.BARE
            ]
        self.habitat_sensitivity = sensitivity
        
        if threats is None:
            threats = [
                HabitatThreat.CROPLAND,
                HabitatThreat.URBAN,
                HabitatThreat.PRIMARY_ROAD,
                HabitatThreat.SECONDARY_ROAD,
                HabitatThreat.TERTIARY_ROAD
            ]
        self.habitat_threats = threats
        
        # data pre-processing
        self.nodata_values = nodata_values

        # default connectivity threshold
        self.connectivity_threshold = connectivity_threshold
        
        # the following method intializes all the filenames needed
        self.map_nbsos_files()

        self.recommended_epsg = Spatial.get_area_crs(self.area_of_interest_filename)
        biodiversity_logger.info(f"Recommended EPSG for project area is {self.recommended_epsg}")
        
        # reproject input data sources into module input data folder
        output_folder = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'input_data')
        
        # for now, reproject input files; later on, probably reproject in-memory 
        self.area_of_interest_filename = self.project_shapefile_datasource(self.area_of_interest_filename)
        self.hexagon_filename = self.project_shapefile_datasource(self.hexagon_filename)
        self.roads_filename = self.project_shapefile_datasource(self.roads_filename)

        # store reference to the original raster metadata 
        with rasterio.open(self.baseline_land_cover_filename) as src_file:
            self.src_raster_metadata = src_file.meta.copy()
        
        # now continue projecting the land cover to recommended EPSG/CRS 
        # this will also update the internal filename to point to the projected raster
        self.baseline_land_cover_filename = Spatial.project_raster_datasource(self.baseline_land_cover_filename, self.recommended_epsg, output_folder)

        # finally, do some further initial processing to prepare the data needed for the module
        shape_of_aoi = gpd.read_file(self.area_of_interest_filename)
        
        # get land cover dataset reader and get metadata for raster
        # reader as well as raster_metadata refer to the projected raster
        self.land_cover_dataset_reader = rasterio.open(self.baseline_land_cover_filename)
        self.projected_raster_metadata = self.land_cover_dataset_reader.meta.copy()
        self.projected_raster_metadata.update({'dtype' : np.float32}) # this is basically for all exports

        # get resolution
        self.resolution = self.land_cover_dataset_reader.res
        biodiversity_logger.info(f"Resolution for current land cover layer is {self.resolution}")

        # clip to aoi
        self.mtx_clipped_land_cover, clip_affine_transform = Spatial.clip_to_area_of_interest(self.land_cover_dataset_reader, shape_of_aoi['geometry'])

        # make a rasterized version of the hexagon file using neigh_id as raster value
        shape_of_hex = gpd.read_file(self.hexagon_filename)
        self.mtx_hexagons = Spatial.vector_to_raster(shape_geodf=shape_of_hex, out_shape=self.mtx_clipped_land_cover.shape, out_transform=self.land_cover_dataset_reader.transform, fill_attribute='neigh_id')
        
        self.write_raster_file(os.path.join('input_data', f'{self.city_name}_hexagons.tif'), self.mtx_hexagons) # more for debug 
        
        
        biodiversity_logger.info("Tessellation layer successfully processed.")

        # preprocess road barriers for subsequent integration as barriers into the dataset
        if self.roads_filename is not None:
        
            # open shapefile
            roads = gpd.read_file(self.roads_filename)
            self.geodf_roads = BiodiversityModule.preprocess_road_barriers(roads, road_types=self.road_types)
            
            # mostly for debugging purposes, export preprocessed shape
            # later on, remove this code again
            self.geodf_roads.to_file(f"{self.local_root}/{self.country_name}/{self.city_name}/biodiversity/input_data/processed_roads.shp", driver='ESRI Shapefile')                        
            biodiversity_logger.info("Road barriers successfully processed.")
        


    

    def project_shapefile_datasource(self, in_file: str, to_crs = None) -> str:
        """Reprojects a shapefile to the recommended EPSG and returns the path of the reprojected file. 

        :param in_file: Shapefile to reproject.
        :type in_file: str
        :return: Path to the reprojected shapefile.
        :rtype: str
        """

        dst_crs = self.recommended_epsg if to_crs is None else to_crs
                
        infile_path = Path(in_file)
        output_folder = f"{self.local_root}/{self.country_name}/{self.city_name}/biodiversity/input_data"
        outfile_path = f'{output_folder}/{infile_path.name}'

        indata: gpd.GeoDataFrame = gpd.read_file(in_file)
        biodiversity_logger.info(f"Reprojecting {in_file} from EPSG:{indata.crs.to_epsg()} to EPSG:{dst_crs}")
        outdata = indata.to_crs(dst_crs)
        outdata.to_file(outfile_path)

        return outfile_path




    def map_nbsos_files(self) -> None:
        
        # the AOI file should be local_root/country/city/ancillary/city_AOI.shp
        self.area_of_interest_filename = os.path.join(self.local_root, self.country_name, self.city_name, 'ancillary', f'{self.city_name}_AOI.shp')
        biodiversity_logger.info(f"Mapping AOI to file {self.area_of_interest_filename}")
        if not os.path.exists(self.area_of_interest_filename):
            error_msg = f"{self.area_of_interest_filename} not found."
            biodiversity_logger.error(error_msg)
            raise FileExistsError(error_msg)
        
        # hexagons should be local_root/country/city/ancillary/city_hex.shp
        self.hexagon_filename = os.path.join(self.local_root, self.country_name, self.city_name, 'ancillary', f'{self.city_name}_hex.shp')
        biodiversity_logger.info(f"Mapping hexagons to file {self.hexagon_filename}")
        if not os.path.exists(self.hexagon_filename):
            error_msg = f"{self.hexagon_filename} not found."
            biodiversity_logger.error(error_msg)
            raise FileExistsError(error_msg)
        
        # land cover should be local_root/country/city/eo_indicators/land_cover.tif
        self.baseline_land_cover_filename = os.path.join(self.local_root, self.country_name, self.city_name, 'eo_indicators', 'land_cover.tif')
        biodiversity_logger.info(f"Mapping baseline land cover to file {self.baseline_land_cover_filename}")
        if not os.path.exists(self.baseline_land_cover_filename):
            error_msg = f"{self.baseline_land_cover_filename} not found."
            biodiversity_logger.error(error_msg)
            raise FileExistsError(error_msg)
        
        # roads should be in local_root/country/city/osm/osm_roads_LineString.shp        
        self.roads_filename = os.path.join(self.local_root, self.country_name, self.city_name, 'osm', 'osm_roads_LineString.shp')
        biodiversity_logger.info(f"Mapping OSM roads to file {self.roads_filename}")
        if not os.path.exists(self.roads_filename):
            error_msg = f"{self.roads_filename} not found."
            biodiversity_logger.error(error_msg)
            raise FileExistsError(error_msg)



    @staticmethod
    def preprocess_road_barriers(
        geodf_roads: gpd.GeoDataFrame,
        road_types = ['motorway', 'primary', 'construction', 'secondary', 'tertiary']) -> gpd.GeoDataFrame:
        
        # according to Evariste's method:
        # Motorway were considered to have 30 m diameters​
        # Primary and construction road categories to have 15 m diameter​
        # Secondary and tertiary road category to have 10 m diameter​
        # Other road categories were not considered as fragmenting geometries
        
        # confirm road types used
        print(f"Preprocessing road barriers for road types {road_types}")

        # consequently, filter down the roads to types needed
        geodf_roads = geodf_roads[geodf_roads['highway'].isin(road_types)].copy()
        
        # define diameter values
        diameter_values = {
            'motorway' : 30,
            'primary' : 15,
            'construction' : 15,
            'secondary' : 10,
            'tertiary' : 10    
        }

        # map diameters to road types and buffer geometry using half the diameter width
        geodf_roads["diameter"] = geodf_roads["highway"].map(diameter_values)
        geodf_roads['geometry'] = geodf_roads.geometry.buffer(distance=geodf_roads["diameter"]/2)
        return geodf_roads


    # Note that all outputs should be directed to local_root/country/city/outputs/ and here, likely, a folder biodiversity
    # def write_scenario_summaries(self):

    #     out_df: pl.DataFrame = None
    #     scenario: AssessmentContainer
    #     for scenario in self.assessment_results:
    #         tmp_dict = scenario.get_summary()            
    #         tmp_df = pl.DataFrame(tmp_dict)

    #         out_df = tmp_df if out_df is None else pl.concat([out_df, tmp_df]) 

    #     out_df.write_excel(os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'outputs', 'summary_of_scenarios.xlsx'))


    def write_excel_file(self, filename: str, df: pd.DataFrame) -> None:
        
        out_path = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, filename)
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        # write dataframe to excel
        df.to_excel(out_path, index=False)


    def write_raster_file(self, filename: str, mtx: np.ndarray, raster_metadata = None, match_nbsos_reference_layer: bool = False) -> None:
        
        # reference metadata for the module
        raster_metadata = self.projected_raster_metadata if raster_metadata is None else raster_metadata 

        # make filepath and assert existence of target directory
        out_path = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, filename)
        os.makedirs(os.path.dirname(out_path), exist_ok=True)

        if match_nbsos_reference_layer:

            ref_source = self.land_cover_dataset_reader        
            template_meta = self.src_raster_metadata.copy()
            template_meta.update({
                'dtype' : raster_metadata['dtype'],
            })

            out_ref = rasterio.open(out_path, 'w', **template_meta)
            rasterio.warp.reproject(
                source=mtx,
                destination=rasterio.band(out_ref, 1),
                src_transform=ref_source.transform,
                src_crs=ref_source.crs,
                dst_transform=out_ref.transform,
                dst_crs=out_ref.crs,
                resampling=Resampling.nearest
            )
            out_ref.close()        

        else:
            # write the raster out
            with rasterio.open(out_path, "w", **raster_metadata) as writer:
                writer.write(mtx, 1)
        
    @staticmethod
    def get_biodiversity_protection_benefits(local_root, country_name, city_name, protection_mask_type: str) -> None:
        """
        Get the benefits of biodiversity protection
        """
        
        scenario_name = 'baseline'
        
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
            
        # import habitat quality file
        habitat_quality_filename = os.path.join(local_root, country_name, city_name, BiodiversityModule.get_module_pathname(), 'outputs', 'habitat_quality', scenario_name, f'habitat_quality_{scenario_name}.tif')
        ds_quality = rioxarray.open_rasterio(habitat_quality_filename)

        # match projection to protection datasource
        ds_quality = ds_quality.rio.reproject_match(ds_protection, resampling=Resampling.nearest, align_on_exact=True)

        out_path = os.path.join(local_root, country_name, city_name, 'output', 'protection', 'nbs_protection_biodiversity_benefits.tif')
        ds_protection = ds_quality * ds_protection
        ds_protection.rio.to_raster(out_path, driver='GTiff', dtype='float32')
    
        #mtx_protection = np.multiply(mtx_protection, assessment_result.habitat_quality.habitat_quality_matrix)
        # self.write_raster_file(os.path.join('outputs', 'habitat_quality', assessment_result.name, f'intersected_habitat_quality_{assessment_result.name}.tif'), mtx_protection.astype(np.float32), raster_metadata=meta, match_nbsos_reference_layer=True)
        
        # # additionally assess average habitat quality
        # patch_count, mtx_patches = Spatial._compute_patches(mtx_protection, PatchContiguity.Rook) 
        # patch_metrics: pd.DataFrame = BiodiversityIndicators.mean_of_patch(mtx_patches, mtx_protection) # 'label', 'intensity_mean'

        # mtx_avg_quality = np.zeros(mtx_protection.shape, np.float32) 
        # util.map_array(mtx_patches.astype(np.int32), np.array(patch_metrics['label'].values.tolist()), np.array(patch_metrics['intensity_mean'].values.tolist()), out=mtx_avg_quality)
        # self.write_raster_file(os.path.join('outputs', 'habitat_quality', assessment_result.name, f'averaged_habitat_quality_{assessment_result.name}.tif'), mtx_avg_quality.astype(np.float32), raster_metadata=meta, match_nbsos_reference_layer=True)


    def write_habitat_quality_rasters(self):

        meta = self.projected_raster_metadata.copy()
        meta.update({'dtype' : np.float32})

        assessment_result: ResultCollection
        for assessment_result in tqauto(self.assessment_results, desc="Exporting habitat quality"):
            if assessment_result.contains_habitat_quality():
                # write habitat quality rasters as WGS84
                self.write_raster_file(os.path.join('outputs', 'habitat_quality', assessment_result.name, f'habitat_suitability_{assessment_result.name}.tif'), assessment_result.habitat_quality.habitat_suitability_matrix.astype(np.float32), raster_metadata=meta, match_nbsos_reference_layer=True)
                self.write_raster_file(os.path.join('outputs', 'habitat_quality', assessment_result.name, f'sum_of_degradation_{assessment_result.name}.tif'), assessment_result.habitat_quality.degradation_matrix.astype(np.float32), raster_metadata=meta, match_nbsos_reference_layer=True)
                self.write_raster_file(os.path.join('outputs', 'habitat_quality', assessment_result.name, f'habitat_quality_{assessment_result.name}.tif'), assessment_result.habitat_quality.habitat_quality_matrix.astype(np.float32), raster_metadata=meta, match_nbsos_reference_layer=True)
                
                


    def write_connectivity_rasters(self):

        # export results
        # these rasters are explicitly stored only for baseline
        # baseline is the first result stored in the list of results
        baseline_result: ResultCollection = self.assessment_results[0]
        
        # iterate over the landscape types assessed under basleine                
        for landscape_type, raster_results in baseline_result.rasters.get_results().items():
                
                # write results into appropriate locations
                indicator_name = self.optimization_criterion.name
                self.write_raster_file(os.path.join('outputs', 'connectivity', landscape_type.name, f'baseline_{indicator_name}_in_{landscape_type.name}.tif'), raster_results.scaled_indicator_matrix.astype(np.float32), match_nbsos_reference_layer=True)

                # write some further outputs for connected components, for debugging and inspection purposes
                if self.optimization_criterion is Metrics.ScaledConnectedEffectiveMeshSize or self.optimization_criterion is Metrics.ConnectedEffectiveMeshSize:
                    
                    self.write_raster_file(os.path.join('outputs', 'connectivity', landscape_type.name, f'baseline_component_patch_count_in_{landscape_type.name}.tif'), raster_results.component_patch_count.astype(np.float32), match_nbsos_reference_layer=True)
                    self.write_raster_file(os.path.join('outputs', 'connectivity', landscape_type.name, f'baseline_component_label_in_{landscape_type.name}.tif'), raster_results.component_labels.astype(np.float32), match_nbsos_reference_layer=True)


    def determine_impacts(self) -> None:
        """Iterate over physiognomic landscape types to be assessed, and write out impacts as estimated
        """        
        biodiversity_logger.info("Determine creation opportunity impacts")

        for landscape_type in tqauto(self.landscape_types_to_assess, desc="Determine impacts"):
            self.determine_impact(landscape_type)

    def determine_impact(self, landscape_type: PhysiognomicLandscapeType, write_spreadsheet: bool = True) -> None:
        
        # compare differences to baseline based on hexagons
        # get baseline: this is always the first result set.
        baseline_results: ResultCollection = self.assessment_results[0]

        df_impact_comparison = pd.DataFrame({'hex_label' : [int(f) for f in np.unique(self.mtx_hexagons[self.mtx_hexagons > 0])]})
        df_impact_comparison['baseline'] = baseline_results.metrics.get_result_for_type(landscape_type).metrics[self.optimization_criterion]

        meta = self.projected_raster_metadata.copy()
        meta.update({
            'dtype' : np.float32,
        })

        # iterate over subsequent scenarios and compare them to baseline
        for i in range(1,len(self.assessment_results)):
            
            # scenario_result is actually a Tuple of four items
            scenario_result: ResultCollection = self.assessment_results[i]

            nbs_impacts: CreationOpportunityImpactsContainer = scenario_result.impacts.get_result_for_type(landscape_type)
            
            hex_vals = [float(f) for f in nbs_impacts.get_results().keys()]            
            ind_vals = [f.metrics[self.optimization_criterion] for f in nbs_impacts.get_results().values()]

            col_name = scenario_result.name 
            scenario_df = pd.DataFrame({'hex_label' : hex_vals, col_name : ind_vals})
            df_impact_comparison = pd.merge(df_impact_comparison, scenario_df, how='left', on='hex_label')
            df_impact_comparison[f"d_{col_name}"] = df_impact_comparison[col_name] - df_impact_comparison['baseline']         
            df_impact_comparison.fillna(0, inplace=True)
            
            # now that we have the actual differences we can write these as outputs using skimage
            mtx_scenario_impact = np.zeros(self.mtx_hexagons.shape, np.float32)
            util.map_array(self.mtx_hexagons.astype(np.int32), np.array(df_impact_comparison['hex_label'].values.tolist()), np.array(df_impact_comparison[f"d_{col_name}"].values.tolist()), out=mtx_scenario_impact)            
            
            # write out differences to baseline as dataframe
            df_impact_comparison.to_excel(os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'outputs', 'connectivity', landscape_type.name, f'diffs_{scenario_result.name}_in_{landscape_type.name}.xlsx'), index=False)

            # write out differences to baseline as raster
            self.write_raster_file(os.path.join('outputs', 'connectivity', landscape_type.name, f'diffs_{scenario_result.name}_in_{landscape_type.name}.tif'), mtx_scenario_impact, raster_metadata=meta, match_nbsos_reference_layer=True)
        


        self.df_impact_comparison[landscape_type] = df_impact_comparison
        if write_spreadsheet:
            self.write_excel_file(os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'outputs', 'connectivity', landscape_type.name, f'diffs_in_{landscape_type.name}.xlsx'), df_impact_comparison)
            


    def determine_optimums(self):
        """Identify max. of optimization criterion at the hexagon level, and identify associated scenario/creation opportunity
        """
        biodiversity_logger.info(f"Determine creation opportunity type by optimization criterion {self.optimization_criterion.value}")
        for landscape_type in tqauto(self.landscape_types_to_assess, desc="Determine optimums"):
            self.determine_optimum(landscape_type)
        
    def determine_optimums_from_spreadsheet(self, nbs_scenarios: List[OpportunityType], landscape_type: PhysiognomicLandscapeType, write_spreadsheet: bool = True) -> None:

        biodiversity_logger.info(f"Determine creation opportunity type by optimization criterion from external dataframe")
        self.nbs_scenarios_for_assessment = nbs_scenarios
        
        # import files according to scenarios
        df_impacts = None
        for scenario in self.nbs_scenarios_for_assessment:
            in_filename = os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'outputs', 'connectivity', landscape_type.name, f"diffs_{scenario.name}_in_{landscape_type.name}.xlsx")
            tmp_df = pd.read_excel(in_filename)
            if df_impacts is None:
                df_impacts = tmp_df
            else:
                tmp_df = tmp_df.drop(columns=['baseline'])
                df_impacts = pd.merge(left=df_impacts, right=tmp_df, on="hex_label")

        #df_impacts = pd.read_excel(filename)
        # save import to internal storage
        self.df_impact_comparison[landscape_type] = df_impacts   

        # call optimization method
        self.determine_optimum(landscape_type=landscape_type, write_spreadsheet=write_spreadsheet)

    def determine_optimum(self, landscape_type: PhysiognomicLandscapeType, write_spreadsheet: bool = True):

        # determine max benefit based on df
        # mappings references cols of interest and is used for value replacement later on
        scenario_raster_value_map = {
            'baseline' : 1
        }
        if self.nbs_scenarios_for_assessment is not None:
            for nbs_type in self.nbs_scenarios_for_assessment:
                scenario_raster_value_map[nbs_type.name] = nbs_type.get_land_cover_class_value()

        # determine dataframe of interest, according to landscape type
        df_of_interest = self.df_impact_comparison[landscape_type]
        # copy relevant columns/vars from dataframe of interest
        columns_to_evaluate = [f for f in scenario_raster_value_map.keys()]
        df_impacts = df_of_interest[['hex_label'] + columns_to_evaluate].copy()

        # identify optimums in terms of max(optimization criterion value) and corresponding scenario/creation opportunity
        df_impacts['max_criterion_value'] = df_impacts[columns_to_evaluate].max(axis=1)
        df_impacts['optimum_nbs'] = df_impacts[columns_to_evaluate].idxmax(axis=1)
        df_impacts['optimum_nbs_class'] = df_impacts['optimum_nbs'].map(scenario_raster_value_map)

        # write optimum optimization criterion value to disk
        mtx_max_criterion_value = np.zeros(self.mtx_hexagons.shape, np.float32)
        meta = self.projected_raster_metadata.copy()
        meta.update({
            'dtype' : np.float32
        })

        util.map_array(self.mtx_hexagons.astype(np.int32), np.array(df_impacts['hex_label'].values.tolist()), np.array(df_impacts['max_criterion_value'].values.tolist()), out=mtx_max_criterion_value)
        self.write_raster_file(os.path.join('outputs', 'optimum', landscape_type.name, f'max_criterion_value_in_{landscape_type.name}.tif'), mtx_max_criterion_value, raster_metadata=meta, match_nbsos_reference_layer=True)

        # write optimum nbs result to disk
        mtx_optimum_solution = np.zeros(self.mtx_hexagons.shape, np.int32)
        meta = self.projected_raster_metadata.copy()
        meta.update({
            'dtype' : np.int32,
        })

        util.map_array(self.mtx_hexagons.astype(np.int32), np.array(df_impacts['hex_label'].values.tolist()), np.array(df_impacts['optimum_nbs_class'].values.tolist()), out=mtx_optimum_solution)            
        self.write_raster_file(os.path.join('outputs', 'optimum', landscape_type.name, f'optimum_nbs_in_{landscape_type.name}.tif'), mtx_optimum_solution, raster_metadata=meta, match_nbsos_reference_layer=True)
        
        if write_spreadsheet:
            self.write_excel_file(os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'outputs', 'optimum', landscape_type.name, f"max_criterion_value_in_{landscape_type.name}.xlsx"), df_impacts)

    def get_preprocessed_roads(self) -> gpd.GeoDataFrame:
        return self.geodf_roads


    def assess_scenario(self, 
                        nbs_scenario: Tuple[OpportunityType, np.ndarray] = None, 
                        sensitivity = None, 
                        threats = None, 
                        assess_habitat_quality: bool = True, 
                        connectivity_threshold: float = 100,
                        min_pixels: int = None) -> ResultCollection:
        
        # keep results and scenario-specific data in encapsulated form
        # assert single-responsibility principle for methods by implementing them within the proper scope
        current_scenario = AssessmentContainer(
            mtx_land_cover_baseline=self.mtx_clipped_land_cover, 
            mtx_hexagons=self.mtx_hexagons,
            nbs_scenario=nbs_scenario, 
            habitat_types=self.baseline_habitat_types, 
            landscape_types=self.landscape_types_to_assess,
            road_barriers=self.get_preprocessed_roads(), 
            sensitivity=sensitivity,
            threats=threats,
            resolution=self.resolution, 
            transform=self.land_cover_dataset_reader.transform,
            require_optimization_criterion=self.optimization_criterion,
            assess_habitat_quality=assess_habitat_quality, 
            connectivity_threshold=connectivity_threshold,
            min_pixels=min_pixels
        )      

        # results is a Tuple with various results collected from assessment  
        current_results: ResultCollection = current_scenario.assess()    

        # debug outputs can go here
        return current_results



    def run_module(self, nbs_scenarios: List[OpportunityType] = None, optimization_criterion: Metrics = None, assess_habitat_quality: bool = True, connectivity_threshold: float = 100, min_pixels: int = None) -> None:
        
        # store optimization criterion
        self.optimization_criterion = optimization_criterion

        # assess baseline
        # for baseline, also take care of habitat quality
        # baseline results is a tuple of 4
        baseline_results: ResultCollection = self.assess_scenario(
            sensitivity=self.habitat_sensitivity, 
            threats=self.habitat_threats, 
            assess_habitat_quality=assess_habitat_quality,
            connectivity_threshold=connectivity_threshold,
            min_pixels=min_pixels
        )

        

        # to show some progress, test-print a summary of the metrics
        baseline_results.metrics.print_summary()
        self.assessment_results.append(baseline_results)
        
        # take care of exports etc.
        if assess_habitat_quality:
            self.write_habitat_quality_rasters()
        
        # iterate over creation opportunities, if any
        if nbs_scenarios is not None:

            self.nbs_scenarios_for_assessment = nbs_scenarios

            scenario: OpportunityType
            for scenario in nbs_scenarios:
                
                # we have to project this layer into the input_data, and then import the projected raster for subsequent analysis
                # suitability layers for nbs should be local_root/country/city/suitability/{intervention type}_creation.tif 
                scenario_filename = os.path.join(self.local_root,self.country_name,self.city_name, 'suitability', f"{scenario.name}.tif")
                
                # assert file reprojection to match template, which is the projected land cover raster
                reprojected_file = Spatial.match_reference_raster(scenario_filename, self.projected_raster_metadata.copy(), os.path.join(self.local_root, self.country_name, self.city_name, self.module_pathname, 'input_data'))

                # import scenario to consume in the scenario
                mtx_scenario = rasterio.open(reprojected_file).read(1)
                scenario_result: ResultCollection = self.assess_scenario(nbs_scenario=(scenario, mtx_scenario), connectivity_threshold=connectivity_threshold, min_pixels=min_pixels)
                self.assessment_results.append(scenario_result)

        

        # write connectivity raster for landscape types assessed under baseline
        # note again that no raster results are stored for the creation opportunities
        self.write_connectivity_rasters()

        if self.nbs_scenarios_for_assessment is not None:
            # once all scenarios are assessed, continue to assess impacts and write results
            self.determine_impacts()
            # identify optima
            self.determine_optimums()

        else:
            biodiversity_logger.warning("Assessing of impacts skipped - no scenarios.")

        # we seem to be done at this point.
        biodiversity_logger.info(f"Assessment of potential biodiversity benefits for {self.city_name} in {self.country_name} completed. Have a nice day!")

    @staticmethod
    def get_module_pathname() -> str:
        return "Biodiversity"

    @staticmethod
    def run(
        local_root: str, 
        country: str, 
        city: str, 
        opportunities: List[OpportunityType] = None,
        optimization_criterion: Metrics = Metrics.ScaledConnectedEffectiveMeshSize,
        land_covers_considered_habitat: List[LandCover] = None, 
        landscape_types_to_assess: List[PhysiognomicLandscapeType] = [PhysiognomicLandscapeType.Any], 
        assess_habitat_quality: bool = True,
        habitat_relative_sensitivity = None, 
        habitat_threats = None,
        nodata_values: List[int] = [0],
        connectivity_threshold: float = 100,
        road_types = ['motorway', 'primary', 'construction', 'secondary', 'tertiary'],
        min_pixels: int = None
        ):
       
        print(f"Running biodiversity module for {city} in {country} with optimization criterion {optimization_criterion.value} and minimum pixel threshold {min_pixels} and connectivity threshold {connectivity_threshold}")
        
        # instantiate the biodiversity module
        module_instance = BiodiversityModule(
            local_root=local_root,
            country=country,
            city=city
        )

        # set data and environment. If arguments not specified by user, default values will be passed or generated, if None is passed. 
        module_instance.set_environment(
            baseline_habitat_types=land_covers_considered_habitat,
            landscape_types=landscape_types_to_assess,
            sensitivity=habitat_relative_sensitivity,
            threats=habitat_threats,
            nodata_values=nodata_values,
            connectivity_threshold=connectivity_threshold,
            road_types=road_types
        )

        # run the biodiversity module
        module_instance.run_module(
            nbs_scenarios=opportunities,
            optimization_criterion=optimization_criterion,
            assess_habitat_quality=assess_habitat_quality,
            connectivity_threshold=connectivity_threshold,
            min_pixels=min_pixels
        )
