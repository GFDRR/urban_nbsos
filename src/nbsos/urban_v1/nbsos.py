from .scan_context import ScanContext
from .data_module.data_module import DataModule

from .suitability import Suitability
from .clustering import CreationCluster
from .priority import get_priority_scores, get_nbs_mca, write_binned_raster

from .scs_cn import CN
from .flood_demand import SpatialRisk
from .pluvial_flood import Runoff
from .water_storage import Storage
from .recreation import get_300_demand, get_300_protection
from .heat_module.heat_module import HeatModule
from .heat_module.heat_exports import RegressionPredictors
from .heat_protection import HeatProtection
from .biodiversity_module.biodiversity_module import BiodiversityModule
from .biodiversity_module.biodiversity_exports import PhysiognomicLandscapeType, OpportunityType, Metrics, LandCover, HabitatThreat
from .erosion_module import ErosionModule
from .jobs_module.jobs_module import JobCreationModule

from .loggerconfig import config_logger, save_local_logfile, save_logfile_to_cloud
from .nbsos_exports import NbsTypes, Benefits, ScanBenefit

from rasterio.enums import Resampling

import mapclassify
import geopandas as gpd
import exactextract as ee
from .gis_tools import rasterize_df
from .gis_tools import make_exact_zonal_statistics
from .gis_tools import export_raster
import rioxarray
import rasterio
from typing import List, Dict, Tuple
from tqdm.auto import tqdm as tqauto
import numpy as np
import logging
import os
import sys

from importlib.metadata import version


# ******************************************************************************************************************** #
#                                                  RUN NBSOS                                                           #
# ******************************************************************************************************************** #

class Scan:


    def __init__(self, nbs_list: List[NbsTypes], scan_context: ScanContext):
        """Initialize the NBSOS scan.

        :param nbs_list: List of NbS families of interest for the scan.
        :type nbs_list: List[NbsTypes]
        :param scan_context: Scan context containing environment parameters for the scan.
        :type scan_context: ScanContext
        """
        self.local_root = scan_context.local_root
        self.country = scan_context.country
        self.city = scan_context.city
        self.africa = scan_context.location_in_africa
        self.area = scan_context.local_area_of_interest_file
        self.scan_context = scan_context

        # initialize logging
        logger = logging.getLogger('nbsos.scan')
        logging_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt="%Y-%m-%d %H:%M:%S")
        logging_stream_handler = logging.StreamHandler()
        logging_stream_handler.setFormatter(logging_formatter)        
        
        self.logger = config_logger(logger, mode='w')
        self.logger.addHandler(logging_stream_handler)  # add the main logger to the console output        
        self.logger.info(f"Starting NBSOS for {self.city}, {self.country}")

        self.nbs_list = nbs_list

        # confirmation message to indicate the package version being used
        print(f"Urban NBSOS version: {version('nbsos.urban_v1')}")
        



    def suitability(self, ignore_annual_precipitation: bool = False):
        self.logger.info('Calculate suitability for nbs list')
        su = Suitability(self.country, self.city, self.local_root, ignore_annual_precipitation=ignore_annual_precipitation)
        
        if self.nbs_list is not None:
            nbs_of_interest = [t.get_name() for t in self.nbs_list]
        else:
            nbs_of_interest = None

        su.export_local(nbs_of_interest)


    def explore_clusters(self) -> Tuple[CreationCluster, np.ndarray]:
        self.logger.info(f'Clustering LPPI and bare soil frequency for exploration')
        c = CreationCluster(self.local_root, self.country, self.city)
        mask = c.make_clusters_for_exploration()
        return c, mask        

    def clustering(self, threshold: int = 5) -> CreationCluster:
        self.logger.info(f'Perform clustering based on LCCI and bare soil frequency with threshold {threshold}')
        c = CreationCluster(self.local_root, self.country, self.city)
        # write protection layers
        c.get_protection_mask(threshold=threshold)
        # write creation layers
        for nbs in self.nbs_list:
            c.nbs_creation(nbs.get_name())

        return c
        
        
        
    def write_combined_creation_mask(self):

        # write combined creation opportunity layer
        # open hex raster to get metadata and mask AOI by zero values in the hex raster
        template_raster = os.path.join(self.local_root, self.country, self.city, 'ancillary', f'{self.city}_aoi_hex.tif')
        with rasterio.open(template_raster) as rst:
            out_meta = rst.meta.copy()
            out_shape = rst.shape
            
        out_meta.update({
            'nodata' : np.nan,
            'dtype' : np.float32
        })

        # get creation masks for each nbs
        # determine combined areas for masking purposes
        self.logger.info('Write combined creation opportunity layer')
        mtx_combined_suitability = np.zeros(out_shape, dtype=np.float32)
        for nbs in self.nbs_list:
            src_nbs_creation = rasterio.open(os.path.join(self.local_root, self.country, self.city, 'suitability', f'{nbs.get_name()}_creation.tif'))
            mtx_nbs_creation = src_nbs_creation.read(1)
            mtx_combined_suitability = np.add(mtx_combined_suitability, mtx_nbs_creation)
        
        mtx_combined_suitability[mtx_combined_suitability == 0] = np.nan        
        out_raster = os.path.join(self.local_root, self.country, self.city, 'output', 'mca', f'suitabilities_overlaid.tif')
        export_raster(mtx_combined_suitability, template_raster, out_raster, out_meta)
        
        mtx_combined_suitability = np.where(mtx_combined_suitability > 0, 1.0, 0.0)
        mtx_combined_suitability[mtx_combined_suitability == 0] = np.nan
        out_raster = os.path.join(self.local_root, self.country, self.city, 'output', 'mca', f'suitabilities.tif')
        export_raster(mtx_combined_suitability, template_raster, out_raster, out_meta)


    def get_cns(self):
        self.logger.info(f'Calculating curve numbers')
        # TODO: Change default path ; deploy as heat data
        cns = CN(self.local_root, self.country, self.city)
        cns.get_cns([t.get_name() for t in self.nbs_list])

    def get_flood_demand(self, return_periods: List[int] = None):
        self.logger.info(f'Calculating flood mititgation demand')

        # set default values for flood paths
        if return_periods is None:
            return_periods = [5, 10, 20, 50, 100, 200, 500, 1000]
           
        risk_out = os.path.join(self.local_root, self.country, self.city, 'output', 'flood_risk.tif')
        flood_demand_calc = SpatialRisk(self.local_root, self.country, self.city, return_periods)
        flood_demand_calc.get_risk(risk_out)

    def get_runoff(self):
        self.logger.info(f'Calculating flood runoff reduction potential')
        runoff = Runoff(self.local_root, self.country, self.city)
        runoff.get_runoff()

    def get_storage(self):
        self.logger.info(f'Calculating flood storage potential')
        demand_ras = os.path.join(self.local_root, self.country, self.city, 'output', 'flood_risk.tif')
        nbss = [nbs.get_name() for nbs in self.nbs_list if (nbs == NbsTypes.BioretentionAreas) | (nbs == NbsTypes.OpenGreenSpace) | (nbs == NbsTypes.InlandWetlands)]
        for nbs in nbss:
            output = os.path.join(self.local_root, self.country, self.city, 'flood_benefit', f'flood_benefit_{nbs}_storage.tif')
            s = Storage(self.local_root, self.country, self.city, nbs)
            s.get_label_names()
            if nbs == 'open_green_space':
                new_runoff = s.get_all_runoff(divider=0.8, max_rel=0.03, max_a=100)  ### bioretention: 0.6, 0.05, ogs: 0.8, 0.03, max_a: 100
                new_benefits = s.calc_benefits(new_runoff, demand_ras)
            else:
                new_runoff = s.get_all_runoff(divider=0.6, max_rel=0.05, max_a=100)  ### bioretention: 0.6, 0.05, ogs: 0.8, 0.03, max_a: 100 # changed from 400
                new_benefits = s.calc_benefits(new_runoff, demand_ras)

            export_raster(new_benefits, demand_ras, output)


    def get_recreation_creation_benefits(self, 
                       requested_benefit: Benefits,
                       min_size=None, 
                       bioretention_area=1, 
                       open_green_space=3, 
                       green_corridors=2,
                       urban_farming=3, 
                       urban_forest=2, 
                       river_renaturation=1,
                       river_floodplain=1, 
                       inland_wetland=1,
                       terraces_and_slopes=1
                       ):
        """Estimate recreation, health and social cohesion benefits of NbS.

        :param min_size: Minimum size of NbS/greenspaces to be considered for benefits. By default, a minimum size of 1ha is considered, hence, defaults to None (in units of pixels). Value is chosen automatically depending on requested type of benefit, but can be overridden may any value passed that is not none.
        :type min_size: int, optional
        :param bioretention_area: Importance of bioretention areas for providing benefits, defaults to 1.
        :type bioretention_area: int, optional
        :param open_green_space: Importance of open green spaces for providing benefits, defaults to 3.
        :type open_green_space: int, optional
        :param green_corridors: Importance of green corridors for providing benefits, defaults to 2.
        :type green_corridors: int, optional
        :param urban_farming: Importance of urban farming for providing benefits, defaults to 3.
        :type urban_farming: int, optional
        :param urban_forest: Importance of urban forest for providing benefits, defaults to 2.
        :type urban_forest: int, optional
        :param river_renaturation: Importance of river renaturation for providing benefits, defaults to 1.
        :type river_renaturation: int, optional
        :param river_floodplain: Importance of river floodplain for providing benefits, defaults to 1.
        :type river_floodplain: int, optional
        :param inland_wetland: Importance of inland wetlands for providing benefits, defaults to 1.
        :type inland_wetland: int, optional
        :param terraces_and_slopes: Importance of terraces and slopes for providing benefits, defaults to 1.
        :type terraces_and_slopes: int, optional
        """
        
        
        mask = os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'dem.tif')
        loop_list = [
            {'nbs': 'bioretention_area',
             'weight': bioretention_area},
            {'nbs': 'open_green_space',
             'weight': open_green_space},
            {'nbs': 'green_corridors',
             'weight': green_corridors},
            {'nbs': 'urban_farming',
             'weight': urban_farming},
            {'nbs': 'urban_forest',
             'weight': urban_forest},
            {'nbs': 'river_renaturation',
             'weight': river_renaturation},
            {'nbs': 'inland_wetland',
             'weight': inland_wetland},
            {'nbs': 'river_floodplain',
             'weight': river_floodplain},
            {'nbs': 'terraces_and_slopes',
             'weight': terraces_and_slopes}
        ]

        requested_nbs_types = [t.get_name() for t in self.nbs_list]
        nbss = [nbs for nbs in loop_list if nbs['nbs'] in requested_nbs_types]

        # DIFFERENTIATE BETWEEN SOCIAL COHESION AND RECREATION                        
        if requested_benefit is Benefits.Recreation:
            self.logger.info(f'Assessing recreation benefits')
            file_prefix = "recreation"
            if min_size is None:
                min_size = 100        
        elif requested_benefit is Benefits.SocialCohesion:
            self.logger.info(f'Assessing social cohesion benefits')
            file_prefix = "social_cohesion"
            if min_size is None:
                min_size = 50

        for nbs in tqauto(nbss):
            nbs_output = os.path.join(self.local_root, self.country, self.city, 'rec_benefit', f'{file_prefix}_benefit_{nbs["nbs"]}.tif')
            nbs_potential, demand = get_300_demand(self.local_root, self.country,
                                                              self.city, nbs['nbs'], nbs['weight'],
                                                              min_size=min_size)
            
            # write recreation benefits for current nbs
            export_raster(nbs_potential, mask, nbs_output)

            # make exact zones for benefits
            shapes = os.path.join(self.local_root, self.country, self.city, 'ancillary', f'{self.city}_hex.shp')
            gdf = gpd.read_file(shapes)

            ds_nbs_output = rioxarray.open_rasterio(nbs_output)
            gdf = gdf.to_crs(ds_nbs_output.rio.crs)

            tmp_stats = ee.exact_extract(ds_nbs_output, gdf, 'sum', output='pandas')
            gdf['tmp'] = tmp_stats['sum']    
            
            rasterize_df(
                gdf,
                os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'land_cover.tif'),
                os.path.join(self.local_root, self.country, self.city, 'rec_benefit', f'{file_prefix}_benefit_{nbs["nbs"]}_hex.tif'),
                'tmp'
            )
        
        # write recreation demand
        demand_output = os.path.join(self.local_root, self.country, self.city, 'output', f'{file_prefix}_demand.tif')
        export_raster(demand, mask, demand_output)

        # make tiled version of recreation demand
        gdf = gpd.read_file(shapes)
        ds_demand = rioxarray.open_rasterio(demand_output)
        gdf = gdf.to_crs(ds_demand.rio.crs)
        tmp_stats = ee.exact_extract(ds_demand, gdf, 'sum', output='pandas')
        gdf['tmp'] = tmp_stats['sum']    
        
        rasterize_df(
            gdf,
            os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'land_cover.tif'),
            os.path.join(self.local_root, self.country, self.city, 'output', f'{file_prefix}_demand_hex.tif'),
            'tmp'
        )



    def get_biodiversity_creation_benefits(
            self, 
            assess_habitat_quality: bool = True,
            optimization_criterion: Metrics = Metrics.ScaledConnectedEffectiveMeshSize,
            land_covers_considered_habitat: List[LandCover] = None, 
            landscape_types_to_assess: List[PhysiognomicLandscapeType] = [PhysiognomicLandscapeType.Any], 
            habitat_relative_sensitivity: List[LandCover] = None, 
            habitat_threats: List[HabitatThreat] = None,
            nodata_values: List[int] = [0],
            connectivity_threshold: float = 100,
            road_types = ['motorway', 'primary', 'construction', 'secondary', 'tertiary'],
            min_pixels: int = None):   
        
        """Assess habitat quality and connectivity as proxy indicators for supporting biodiversity. 
        Impacts of NbS on connectivity are assessed by comparing specified optimization criterion metrics with baseline.

        :param assess_habitat_quality: Assess habitat quality and degradation for baseline, defaults to True
        :type assess_habitat_quality: bool, optional
        :param optimization_criterion: Metric used for assessing NbS impacts, defaults to Metrics.ScaledEffectiveMeshSize
        :type optimization_criterion: Metrics, optional
        :param land_covers_considered_habitat: List of land cover types considered as habitats. If None, use default list, defaults to None
        :type land_covers_considered_habitat: List[LandCover], optional
        :param landscape_types_to_assess: Physiognomic landscape types to assess, defaults to [PhysiognomicLandscapeType.Any]
        :type landscape_types_to_assess: List[PhysiognomicLandscapeType], optional
        :param habitat_relative_sensitivity: List of Land Cover to consider in habitat quality assessment. If None, use default classes, defaults to None
        :type habitat_relative_sensitivity: List[LandCover], optional
        :param habitat_threats: List of HabitatThreats to consider in habitat quality assessment. If None, use default classes, defaults to None
        :type habitat_threats: List[HabitatThreat], optional
        :param nodata_values: List of values in the land cover raster to consider as nodata, defaults to [0]
        :type nodata_values: List[int], optional
        :param connectivity_threshold: Distance threshold for structural connectivity assessment, defaults to 100m.
        :type connectivity_threshold: float, optional
        :param road_types: List of road types to consider as barriers for connectivity assessment, defaults to ['motorway', 'primary', 'construction', 'secondary', 'tertiary']
        :type road_types: List[str], optional
        :param min_pixels: Minimum number of pixels for a habitat patch to be considered in the connectivity assessment. If None, no minimum is applied, defaults to None.
        :type min_pixels: int, optional
        """
                    
        

        self.logger.info("Starting assessment of biodiversity benefits.")

        # map nbs in scan list to opportunity types
        # this will be changed going forward to use enumerations
        opportunity_types = []
        if self.nbs_list is not None:
            for type in self.nbs_list:
                # match NBSOS NbsTypes to BiodiversityModule OpportunityType
                # this is done explicitly not by Enum name matching, as string values for 
                # NbsTypes might change 
                if type is NbsTypes.BioretentionAreas:
                    opportunity_types.append(OpportunityType.bioretention_area_creation)
                elif type is NbsTypes.GreenCorridors:
                    opportunity_types.append(OpportunityType.green_corridors_creation)
                elif type is NbsTypes.InlandWetlands:
                    opportunity_types.append(OpportunityType.inland_wetland_creation)
                elif type is NbsTypes.OpenGreenSpace:
                    opportunity_types.append(OpportunityType.open_green_space_creation)
                elif type is NbsTypes.RiverRenaturation:
                    opportunity_types.append(OpportunityType.river_renaturation_creation)
                elif type is NbsTypes.TerracesAndSlopes:
                    opportunity_types.append(OpportunityType.terraces_and_slopes_creation)
                elif type is NbsTypes.UrbanFarming:
                    opportunity_types.append(OpportunityType.urban_farming_creation)
                elif type is NbsTypes.UrbanForest:
                    opportunity_types.append(OpportunityType.urban_forest_creation)
                elif type is NbsTypes.RiverFloodplain:
                    opportunity_types.append(OpportunityType.river_floodplain_creation)
                #elif type is NbsTypes.BuildingSolutions:
                #    opportunity_types.append(OpportunityType.building_solutions_creation)

        # run the module
        BiodiversityModule.run(
            local_root=self.local_root,
            country=self.country,
            city=self.city,
            assess_habitat_quality=assess_habitat_quality,
            opportunities=opportunity_types,
            optimization_criterion=optimization_criterion,
            land_covers_considered_habitat=land_covers_considered_habitat,
            landscape_types_to_assess=landscape_types_to_assess,
            habitat_relative_sensitivity=habitat_relative_sensitivity,
            habitat_threats=habitat_threats,
            nodata_values=nodata_values,
            connectivity_threshold=connectivity_threshold,
            road_types=road_types,
            min_pixels=min_pixels
        )

    

    def get_job_creation_creation_benefits(
            self,
            nbs_list: List[NbsTypes] = None,
            kernel_radius: float = 1000,
            kernel_sigma: float = 1000/3):
        """Assess jobs creation benefits of NBS families
        
        :param nbs_list: List of NBS types to assess for job creation benefits. If None, all NBS types in the scan context will be assessed, defaults to None.
        :type nbs_list: List[NbsTypes], optional
        :param kernel_radius: Radius of the Gaussian kernel for population weighting, defaults to 1000m.
        :type kernel_radius: float, optional
        :param kernel_sigma: Standard deviation of the Gaussian kernel for population weighting, defaults to 1000/3.
        :type kernel_sigma: float, optional        
        """

        nbs_list = self.nbs_list if nbs_list is None else nbs_list

        self.logger.info("Calculating job creation benefits")
        JobCreationModule.run(
            scan_context=self.scan_context,
            nbs_list=nbs_list,
            kernel_radius=kernel_radius,
            kernel_sigma=kernel_sigma,
        )        

    def get_pluvial_flood_creation_benefits(
            self, 
            return_periods: List[int] = None,
        ):
        """
        Determine flood creation benefits.

        :param return_periods: List of return periods to consider for flood risk assessment, e.g. [10, 20, 50, 100, 200, 500, 1500] (in years).
        :type return_periods: List[int], optional
        """
        self.logger.info("Calculating pluvial flood creation benefits")
        # get curve numbers
        self.get_cns()
        # get demand for flood mitigation
        self.get_flood_demand(return_periods=return_periods)
        # get runoff reduction
        self.get_runoff()
        # get storage potential
        self.get_storage()


    def get_heat_creation_benefits(
            self, 
            kernel_radius: float = 1000, 
            kernel_sigma: float = 1000/3, 
            temperature_resampling_method: Resampling = Resampling.cubic,
            regression_coefficients_path: str = None, 
            nbs_list: List[ScanBenefit] = None, 
            predictors: List[RegressionPredictors] = [RegressionPredictors.BuiltCoverage, RegressionPredictors.TreeCoverage],
            elevation_covariate: bool = False,
            use_wri_tree_layer: bool = True,
            air_temperature_path: str = None,
            skip_air_temperature_computation: bool = False            
        ):
        """
        Assess urban heat mitigation benefits. By default, regression predictors include built-up coverage and tree coverage. Additional predictors
        may be defined as a list of RegressionPredictors. Geographically weighted and space-weighted regression coefficients that are used to estimate 
        air temperature are needed to run this benefit module. This data is retrieved on-demand.

        :param kernel_radius: Radius of gaussian kernel for determining population weight, defaults to 1000m.
        :type kernel_radius: float, optional
        :param kernel_sigma: Standard deviation of the Gaussian kernel for determining population weight, defaults to 1000/3.
        :type kernel_sigma: float, optional
        :param regression_coefficients_path: Path to geography- and space-weighted climate regression coefficients data, defaults to None. A default path will be used if None is provided, and data is retreived automatically on-demand.
        :type regression_coefficients_path: str, optional
        :param nbs_list: List of NbS to be assessed for heat mitigation benefits, defaults to None. If None, all NbS in the scan context will be used.
        :type nbs_list: List[ScanBenefit], optional
        :param predictors: List of predictors used in the regression, defaults to [RegressionPredictors.BuiltCoverage, RegressionPredictors.TreeCoverage].
        :type predictors: List[RegressionPredictors], optional,
        :param elevation_covariate: Whether to include elevation as covariate in the regression, defaults to False.
        :type elevation_covariate: bool, optional
        :param use_wri_tree_layer: Whether to use the WRI tree cover layer as predictor for tree coverage, defaults to True. If False, tree coverage will be determined based on ESA WorldCover land cover data.
        :type use_wri_tree_layer: bool, optional
        :param air_temperature_path: Path to air temperature data to use in place of internal computation, defaults to None. If None, internal model will be used and air temperature raster be written at the default location.
        :type air_temperature_path: str, optional
        :param skip_air_temperature_computation: Whether to skip the internal computation of air temperature, defaults to False. If True, air temperature will not be computed internally, but the default path for the air temperature raster will be used.
        :type skip_air_temperature_computation: bool, optional
        """

        self.logger.info("Starting assessment of heat benefits.")
        nbs_to_assess = self.nbs_list if nbs_list is None else nbs_list

        # point to a default data folder if no path is explicitly set by user
        regression_coefficients_path = os.path.join(self.local_root, "data", HeatModule.get_module_pathname().lower()) if regression_coefficients_path is None else regression_coefficients_path
        
        # check for availability of required files in this path, and download if not present
        coefficient_data_is_complete = True
        coefficients_data = [
            os.path.join(regression_coefficients_path, "Coeffs_GWR_b0.nc"),
            os.path.join(regression_coefficients_path, "Coeffs_GWR_b1.nc"),
            os.path.join(regression_coefficients_path, "Coeffs_GWR_b2.nc"),
            os.path.join(regression_coefficients_path, "Coeffs_CSWR_b0.nc"),
            os.path.join(regression_coefficients_path, "Coeffs_CSWR_b1.nc"),
            os.path.join(regression_coefficients_path, "Coeffs_CSWR_b2.nc")
        ]
        for f in coefficients_data:
            if not os.path.exists(f):
                coefficient_data_is_complete = False
                break

        if not os.path.exists(regression_coefficients_path) or not coefficient_data_is_complete:
            self.logger.warning(f"Regression coefficients not found in {regression_coefficients_path}. Downloading from GCS.")
            # download the coefficients
            self.storage_handler.get_temperature_regression_coefficients(regression_coefficients_path)
            
        HeatModule.run(
            local_root=self.local_root,
            country=self.country,
            city=self.city,
            kernel_radius=kernel_radius,
            kernel_sigma=kernel_sigma,
            temperature_resampling_method=temperature_resampling_method,
            regression_coefficients_path=regression_coefficients_path,
            nbs_list=nbs_to_assess,
            predictors=predictors,
            elevation_covariate=elevation_covariate,
            use_wri_tree_layer=use_wri_tree_layer,
            air_temperature_path=air_temperature_path,
            skip_air_temperature_computation=skip_air_temperature_computation
        )

    def get_erosion_creation_benefits(self,  nbs_list: List[ScanBenefit] = None, use_population_weighting: bool = True, kernel_radius: float = 1000, kernel_sigma: float = 1000/3):
        """Determine erosion based on the RUSLE equation. Here, avoided eorsion is assessed, i.e., R * K * LS * C. The P factor is not considered in this rapid assessment.        

        
        :param nbs_list: List of NBS to assess, defaults to None. If none, use all NBS used in the scan.
        :type nbs_list: List[ScanBenefit], optional
        :param use_population_weighting: Whether to use population weighting for erosion benefits, defaults to True. If True, population weighting is applied based on sum of population within neighbourhoods.
        :type use_population_weighting: bool, optional
        :param kernel_radius: Radius of the Gaussian kernel for population weighting, defaults to 1000m. Only used if use_population_weighting is True.
        :type kernel_radius: float, optional
        :param kernel_sigma: Standard deviation of the Gaussian kernel for population weighting, defaults to 1000/3. Only used if use_population_weighting is True.
        :type kernel_sigma: float, optional

        """
        self.logger.info("Calculating avoided erosion benefits")
        nbs_to_assess = self.nbs_list if nbs_list is None else nbs_list

        
        ErosionModule.run(
            local_root=self.local_root,
            country_name=self.country,
            city=self.city,
            nbs_list=nbs_to_assess,
            use_population_weighting=use_population_weighting,
            kernel_radius=kernel_radius,
            kernel_sigma=kernel_sigma
        )


    def get_heat_protection_benefits(self, protection_mask_type: str, use_wri_tree_layer: bool = True):
        self.logger.info("Calculating heat protection benefits")
        HeatProtection(self.local_root, self.country, self.city, protection_mask_type, use_wri_tree_layer).calc_heat_protection()
        
        in_file = os.path.join(self.local_root, self.country, self.city, 'output', 'protection', 'nbs_protection_heat_benefits.tif')
        write_binned_raster(
            in_file,
            os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'lppi.tif'),
            os.path.join(self.local_root, self.country, self.city, 'output', 'protection', 'nbs_protection_heat_benefits_binned.tif'),
            binning_method='percentiles',
            percentiles=[0, 25, 75, 100]
        )

    def get_recreation_protection_benefits(self, protection_mask_type: str, min_size=100):
        self.logger.info("Calculating recreation protection benefits")
        out_path = os.path.join(self.local_root, self.country, self.city, 'output', 'protection', 'nbs_protection_health_benefits.tif')
        get_300_protection(self.local_root, self.country, self.city, out_path, protection_mask_type, min_size=min_size)

        write_binned_raster(
            out_path,
            os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'lppi.tif'),
            os.path.join(self.local_root, self.country, self.city, 'output', 'protection', 'nbs_protection_health_benefits_binned.tif'),
            binning_method='percentiles',
            percentiles=[0, 25, 75, 100]
        )

    def get_social_cohesion_protection_benefits(self, protection_mask_type: str, min_size=50):
        self.logger.info("Calculating social cohesion protection benefits")
        out_path = os.path.join(self.local_root, self.country, self.city, 'output', 'protection', 'nbs_protection_cohesion_benefits.tif')
        get_300_protection(self.local_root, self.country, self.city, out_path, protection_mask_type, min_size=min_size) 

        write_binned_raster(
            out_path,
            os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'lppi.tif'),
            os.path.join(self.local_root, self.country, self.city, 'output', 'protection', 'nbs_protection_cohesion_benefits_binned.tif'),
            binning_method='percentiles',
            percentiles=[0, 25, 75, 100]
        )

    def get_biodiversity_protection_benefits(self, protection_mask_type: str):
        self.logger.info("Calculating biodiversity protection benefits")
        BiodiversityModule.get_biodiversity_protection_benefits(self.local_root, self.country, self.city, protection_mask_type)

        in_file = os.path.join(self.local_root, self.country, self.city, 'output', 'protection', 'nbs_protection_biodiversity_benefits.tif')
        write_binned_raster(
            in_file,
            os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'lppi.tif'),
            os.path.join(self.local_root, self.country, self.city, 'output', 'protection', 'nbs_protection_biodiversity_benefits_binned.tif'),
            binning_method='percentiles',
            percentiles=[0, 25, 75, 100]
        )

    def get_erosion_protection_benefits(self, protection_mask_type: str):
        self.logger.info("Calculating erosion protection benefits")
        ErosionModule.get_erosion_protection_benefits(self.local_root, self.country, self.city, protection_mask_type)

        in_file = os.path.join(self.local_root, self.country, self.city, 'output', 'protection', 'nbs_protection_erosion_rate_benefits.tif')
        write_binned_raster(
            in_file,
            os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'lppi.tif'),
            os.path.join(self.local_root, self.country, self.city, 'output', 'protection', 'nbs_protection_erosion_rate_benefits_binned.tif'),
            binning_method='percentiles',
            percentiles=[0, 25, 75, 100]
        )

    def get_pluvial_flood_protection_benefits(self, protection_mask_type: str):
        self.logger.info("Calculating pluvial flood protection benefits")
        SpatialRisk.get_flood_protection_benefits(self.local_root, self.country, self.city, protection_mask_type)

        in_file = os.path.join(self.local_root, self.country, self.city, 'output', 'protection', 'nbs_protection_pluvial_flood_benefits.tif')
        write_binned_raster(
            in_file,
            os.path.join(self.local_root, self.country, self.city, 'eo_indicators', 'lppi.tif'),
            os.path.join(self.local_root, self.country, self.city, 'output', 'protection', 'nbs_protection_pluvial_flood_benefits_binned.tif'),
            binning_method='percentiles',
            percentiles=[0, 25, 75, 100]
        )
        
         


    def get_priorities(
            self, 
            benefits_list: List[ScanBenefit],
            binning_method: str = 'jenks', 
            percentiles: List[int] = [0, 20, 40, 60, 80, 100]):
        
        get_priority_scores(
            self.local_root, 
            self.country, 
            self.city,             
            benefits_list, 
            binning_method=binning_method,
            percentiles=percentiles
        )

        

    def get_mca(self, 
                benefits_list: List[ScanBenefit], 
                heavy_focal_weight: float = 0.8, 
                percentiles: List[int] = [0, 25, 75, 100],
                assess_equal_weights_scenario: bool = True,
                assess_benefit_heavy_scenario: bool = True,
                assess_benefit_only_scenario: bool = True,
                custom_scenarios: Dict[str, Dict[Benefits, float]] = None,
                highest_benefits_percentiles: List[int] = [80],
                focus_area_paths: List[str] = None,
            ):
        """Run multicriteria analysis (MCA) for the selected benefits. The MCA will conduct equal weighting, benefit-heavy weighting (if more than one benefit in model) 
        and benefit-only weighting. Output are rasters with determined scores as well as categorized rasters with bins based on specified percentiles. Default percentiles are 
        the 0th, 25th, 75th and 100th percentile, assuming indicating low, medium and high benefits. 

        :param benefits_list: List of ScanBenefit objects defining benefits to be assessed.
        :type benefits_list: List[ScanBenefit]
        :param heavy_focal_weight: _description_, defaults to 0.8
        :type heavy_focal_weight: float, optional
        :param percentiles: _description_, defaults to [0, 25, 75, 100]
        :type percentiles: List[int], optional
        :param plot: Generate priority area plots, defaults to False.
        :type plot: bool, optional
        :param assess_equal_weights_scenario: Should the scenario with equal weights for all benefits be assessed, defaults to True.
        :type assess_equal_weights_scenario: bool, optional
        :param assess_benefit_heavy_scenario: Should the scenario with heavy weight for one benefit and light weights for other benefits be assessed, defaults to True.
        :type assess_benefit_heavy_scenario: bool, optional
        :param assess_benefit_only_scenario: Should the scenario with only one benefit considered and all other benefits not considered be assessed, defaults to True.
        :type assess_benefit_only_scenario: bool, optional
        :param custom_scenarios: List of custom scenarios to be assessed, where each scenario is defined as a dictionary with benefits as keys and their corresponding weights as values. If provided, these scenarios will be assessed in addition to the default scenarios (equal weighting, benefit-heavy and benefit-only), defaults to None.
        :type custom_scenarios: Dict[str, Dict[Benefits, float]], optional
        :param highest_benefits_percentiles: List of percentiles to be used for identifying highest benefit areas, defaults to [80]. If multiple percentiles are provided, the highest benefit areas will be categorized according to these percentiles.
        :type highest_benefits_percentiles: List[int], optional
        :param focus_area_paths: List of paths to shapefiles defining the focus areas, defaults to None.
        :type focus_area_paths: List[str], optional
        """
        get_nbs_mca(
            self.scan_context, 
            self.nbs_list, 
            benefits_list, 
            heavy_focal_weight=heavy_focal_weight,
            percentiles_for_binning=percentiles,
            assess_equal_weights_scenario=assess_equal_weights_scenario,
            assess_benefit_heavy_scenario=assess_benefit_heavy_scenario,
            assess_benefit_only_scenario=assess_benefit_only_scenario,
            custom_scenarios=custom_scenarios,
            highest_benefits_percentiles=highest_benefits_percentiles,
            focus_area_paths=focus_area_paths,
        )

    def save_logs(self, save_to_cloud=True):
        
        logfile_name = save_local_logfile(self.local_root, self.country, self.city) 
        if save_to_cloud:
            save_logfile_to_cloud(logfile_name, self.country, self.city, self.bucket)


    
    def set_storage_options(self, storage_handler: DataModule) -> None:
        """Assign storage handler to the NBSOS scan and initialize the handler. 

        :param storage_handler: _description_
        :type storage_handler: DataModule
        """
        self.storage_handler = storage_handler
        self.storage_handler.initalize_storage()


    def run_preparation(self, run_gee: bool = True, run_download: bool = True, run_client_computations: bool = True, download_lppi_only: bool = False, include_basic_indicators: bool = True, include_extended_indicators: bool = True, include_rusle: bool = True, include_LST: bool = True):
        """Populate local NBSOS data repository. Running GEE will process indicators and depending on the storage handler, transfer them to Gcs or the local device.
        Running download will download additional external data to the local device, and for the GcsHandler, also download the processed indicators.
        Running client-side computations will extract flood hazard data from the local device into the NBSOS data repository, and will determine various additional data using client-side Python.  

        :param run_gee: Process indicators in Google Earth Engine and pull processed indicators into the chosen storage provider or to the local device, defaults to True.
        :type run_gee: bool, optional
        :param run_download: Download indicators and other external data to the local device, defaults to True.
        :type run_download: bool, optional
        :param run_client_computations: Run client-side computation of indicators, defaults to True.
        :type run_client_computations: bool, optional
        :param download_lppi_only: Only download LPPI indicators, defaults to False.
        :type download_lppi_only: bool, optional
        :param include_basic_indicators: Include basic indicators in the GEE processing, defaults to True.
        :type include_basic_indicators: bool, optional
        :param include_extended_indicators: Include extended indicators in the GEE processing, defaults to True.
        :type include_extended_indicators: bool, optional
        :param include_rusle: Include RUSLE indicators in the GEE processing, defaults to True.
        :type include_rusle: bool, optional
        :param include_LST: Include Land Surface Temperature (LST) indicators in the GEE processing, defaults to True.
        :type include_LST: bool, optional
        """
        self.logger.info("Running scan preparation")
        
        if run_gee:
            self.logger.info("Processing indicators in GEE")
            self.storage_handler.gee_to_storage(download_lppi_only=download_lppi_only,
                                                include_basic_indicators=include_basic_indicators,
                                                include_extended_indicators=include_extended_indicators,
                                                include_rusle=include_rusle,
                                                include_LST=include_LST)
        if run_download:
            self.logger.info("Downloading external data")
            self.storage_handler.storage_to_local()
        if run_client_computations:
            self.logger.info("Running client-side computations")
            self.storage_handler.client_compute_indicators()

        self.logger.info("Scan preparation completed")



    def run_suitability(self, ignore_annual_precipitation: bool = False):
        """This will run the suitability model for the selected NbS types. For river floodplain creation,
        a flood hazard raster is required to delineate the floodplain. The respective flood hazard raster is 
        chosen by specifying its corresponding return period. 

        :param ignore_annual_precipitation: Ignore annual precipitation suitability criterion in the suitability assessment, defaults to False.
        :type ignore_annual_precipitation: bool, optional
        """

        self.logger.info("Running suitability model")
        self.suitability(ignore_annual_precipitation=ignore_annual_precipitation)
    
    def update_creation_mask(self):
        # update mask
        self.write_combined_creation_mask()

    def run_clustering(self, clustering_threshold: int, exclusion_masks : List[str] = None, update_creation_mask: bool = True):
        """
        Run clustering of bare soil frequency and LPPI to determine creation and protection opportunities.
        
        :param clustering_threshold: Threshold value to be used for clustering.
        :type clustering_threshold: int
        :param exclusion_masks: List of paths to shapefiles to be used as exclusion masks in clustering, defaults to None. If provided, creation opportunities within exclusion masks will be removed.
        :type exclusion_masks: List[str], optional
        :param update_creation_mask: Should the creation mask of all NBS be updated. Defaults to True
        :type update_creation_mask: bool, optional
        """
        
        self.logger.info(f"Running bare-soil/LPPI clustering with threshold {clustering_threshold}")
        c = self.clustering(threshold=clustering_threshold)

        # clustering completed. If any exclusion masks are provided, apply them to the creation masks of all nbs types.
        if exclusion_masks is not None:
            nbss = [nbs.get_name() for nbs in self.nbs_list]
            for mask_name in exclusion_masks:
                self.logger.info(f"Applying exclusion mask {mask_name} to creation opportunities")
                c.apply_creation_masks(nbss, mask_name)
        
        # update mask
        if update_creation_mask:
            self.update_creation_mask()


    def run_creation_benefits(self, 
            benefit_parameters: List[ScanBenefit]):
        """
        Run creation benefit models for the selected benefits and NbS families, using the specified clustering 
        threshold for identifying protection and creation opportunities.

        :param benefit_parameters: List of ScanBenefit objects defining benefits to be assessed.
        :type benefit_parameters: List[ScanBenefit]

        """

        self.logger.info(f"Running creation benefits")
        
        if self.nbs_list is None:
            self.logger.info("No NbS selected - creation benefits cannot be assessed.")
        else:
            self.logger.info(f"Selected NbS: {', '.join(map(str, [t.get_name() for t in self.nbs_list]))}")
        
        if benefit_parameters is None:
            self.logger.info("No benefits selected - stopping assessment at this point.")
        else:
            self.logger.info(f"Selected Benefits: {', '.join(map(str, [b.benefit.get_name() for b in benefit_parameters]))}")

            # iterate over benefits
            for benefit in benefit_parameters:

                if benefit.benefit is Benefits.Biodiversity:
                    self.get_biodiversity_creation_benefits(**benefit.parameters)
                elif benefit.benefit is Benefits.UrbanHeat:
                    self.get_heat_creation_benefits(**benefit.parameters)
                elif benefit.benefit is Benefits.Recreation:
                    self.get_recreation_creation_benefits(benefit.benefit, **benefit.parameters)
                elif benefit.benefit is Benefits.SocialCohesion:
                    self.get_recreation_creation_benefits(benefit.benefit, **benefit.parameters)
                elif benefit.benefit is Benefits.PluvialFlooding:
                    self.get_pluvial_flood_creation_benefits(**benefit.parameters)
                elif benefit.benefit is Benefits.AvoidedErosion:
                    self.get_erosion_creation_benefits(**benefit.parameters)
                elif benefit.benefit is Benefits.JobsCreation:
                    self.get_job_creation_creation_benefits(**benefit.parameters)
           
    

    def run_protection_benefits(self, benefit_parameters: List[ScanBenefit], protection_mask_type: str = "lppi"):
        """Run protection benefit models for the selected benefits and NbS families.

        :param benefit_parameters: List of ScanBenefit objects defining benefits to be assessed.
        :type benefit_parameters: List[ScanBenefit]
        :param protection_mask_type: Type of mask to be used for protection benefits assessment, either "suitability" or "lppi", defaults to "lppi". If "suitability" is chosen, the suitability-derived mask will be used, otherwise, the mask of protection opportunities resulting from the LPPI/bare soil frequency clustering.
        :type protection_mask_type: str, optional
        """
        self.logger.info("Running protection benefits")
        
        if self.nbs_list is None:
            self.logger.info("No NbS selected - creation benefits cannot be assessed.")
        else:
            self.logger.info(f"Selected NbS: {', '.join(map(str, [t.get_name() for t in self.nbs_list]))}")
                
        if benefit_parameters is None:
            self.logger.info("No benefits selected - stopping assessment at this point.")
        else:

            self.logger.info(f"Selected Benefits: {', '.join(map(str, [b.benefit.get_name() for b in benefit_parameters]))}")
            for benefit in benefit_parameters:
                if benefit.benefit is Benefits.UrbanHeat:      
                    use_wri_tree_layer = benefit.parameters.get('use_wri_tree_layer', True)              
                    self.get_heat_protection_benefits(protection_mask_type=protection_mask_type, use_wri_tree_layer=use_wri_tree_layer)
                elif benefit.benefit is Benefits.Recreation:
                    self.get_recreation_protection_benefits(protection_mask_type, **benefit.parameters)
                elif benefit.benefit is Benefits.SocialCohesion:
                    self.get_social_cohesion_protection_benefits(protection_mask_type, **benefit.parameters)    
                elif benefit.benefit is Benefits.Biodiversity:
                    self.get_biodiversity_protection_benefits(protection_mask_type)
                elif benefit.benefit is Benefits.PluvialFlooding:
                    self.get_pluvial_flood_protection_benefits(protection_mask_type)
                elif benefit.benefit is Benefits.AvoidedErosion:
                    self.get_erosion_protection_benefits(protection_mask_type)


    def run_priority(self, benefit_parameters: List[ScanBenefit], binning_method: str = 'jenks', percentiles: List[int] = [0, 20, 40, 60, 80, 100]):
        """Identify priority areas for the selected benefits. 
        This method will write outputs in both shapefile and raster format, including categorized rasters with bins according to the 
        specified percentiles.

        :param benefit_parameters: List of ScanBenefit objects defining benefits to be assessed.
        :type benefit_parameters: List[ScanBenefit]
        :param binning_method: Method used for binning (categorizing) priority scores, either 'jenks' or 'percentiles', defaults to "jenks".
        :param binning_method: str, optional
        :param percentiles: List of percentiles to be used for binning, defaults to [0, 20, 40, 60, 80, 100].
        :type percentiles: List[int], optional
        :param plot: Generate priority area plots, defaults to False.
        :param plot: bool, optional 
        """

        # get priorities
        self.get_priorities(benefit_parameters, binning_method=binning_method, percentiles=percentiles)

     
    def run_multicriteria_analysis(
            self, 
            benefit_parameters: List[ScanBenefit], 
            heavy_focal_weight: float = 0.8,
            percentiles: List[int] = [0, 25, 75, 100],
            assess_equal_weights_scenario: bool = True,
            assess_benefit_heavy_scenario: bool = True,
            assess_benefit_only_scenario: bool = True,
            custom_scenarios: Dict[str, Dict[Benefits, float]] = None,
            highest_benefits_percentiles: List[int] = [80],
            focus_area_paths: List[str] = None,
        ):
        
        """
        Run a multicriteria analysis for selected benefits. This method will determine highest benefits and the optimal NbS for a scenario of equal 
        importance of all benefits (equal weighting), and for scenarios to maximize a specific NbS benefit (benefit-only scenarios). If more than one 
        benefit is being assessed, in addition, benefit-heavy scenarios will be run, where a higher weight (focal weight) is given to a focal benefit compared to 
        the remaining benefits, that will be equally weighted. 

        The method will write outputs as rasters, including assessed scores and categorized rasters with bins according to the specified percentiles. Default bins
        are the 0th, 25th, 75th and 100th percentiles, which are assumed to represent low, medium and high benefits.
         
        :param benefit_parameters: List of ScanBenefit objects defining benefits to be assessed. For benefit-only and benefit-heavy scenarios, each benefit in the list will be treated as a focal benefit.
        :type benefit_parameters: List[ScanBenefit]
        :param heavy_focal_weight: Focal weight of focal benefit in benefit-heavy scenarios, defaults to 0.8.
        :type heavy_focal_weight: float, optional
        :param percentiles: List of percentiles used for binning benefits, defaults to List[0, 25, 75, 100].
        :type percentiles: List[int], optional
        :param assess_equal_weights_scenario: Should the scenario with equal weights for all benefits be assessed, defaults to True.
        :type assess_equal_weights_scenario: bool, optional
        :param assess_benefit_heavy_scenario: Should the scenario with heavy weight for one benefit and light weights for other benefits be assessed, defaults to True.
        :type assess_benefit_heavy_scenario: bool, optional
        :param assess_benefit_only_scenario: Should the scenario with only one benefit considered and all other benefits not considered be assessed, defaults to True.
        :type assess_benefit_only_scenario: bool, optional
        :param custom_scenarios: Dictionary defining custom scenarios with specific weights for each benefit, defaults to None.
        :type custom_scenarios: Dict[str, Dict[Benefits, float]], optional
        :param highest_benefits_percentiles: List of percentiles to determine highest benefits, defaults to [80].
        :type highest_benefits_percentiles: List[int], optional
        :param focus_area_paths: List of paths to shapefiles defining the focus areas, defaults to None.
        :type focus_area_paths: List[str], optional
        """

        # run mca that will create most of the key outputs
        self.get_mca(
            benefit_parameters, 
            heavy_focal_weight=heavy_focal_weight, 
            percentiles=percentiles,
            assess_equal_weights_scenario=assess_equal_weights_scenario,
            assess_benefit_heavy_scenario=assess_benefit_heavy_scenario,
            assess_benefit_only_scenario=assess_benefit_only_scenario,
            custom_scenarios=custom_scenarios,
            highest_benefits_percentiles=highest_benefits_percentiles,
            focus_area_paths=focus_area_paths)


    def run_building_solutions(
            self, 
            benefit_parameters: List[ScanBenefit], 
            jenks_classes: int = 5, 
            protection_mask_type: str = "lppi", 
            building_source: str = "building_mask.tif",
            exclude_non_built_areas: bool = False):
        """
        Run building solutions based on identified heat/flood priorities as a 
        function of building density and (lack of) protection and creation opportunities. 

        :param benefit_parameters: List of ScanBenefit objects defining benefits to be assessed. This will be used to evaluate if flood and/or heat benefits are of interest.
        :type benefit_parameters: List[ScanBenefit]
        :param jenks_classes: Number of classes to be used for categorizing building solution scores based on Jenks natural breaks, defaults to 5.
        :type jenks_classes: int, optional
        :param protection_mask_type: Type of mask to be used for protection benefits assessment, either "suitability" or "lppi", defaults to "lppi". If "suitability" is chosen, the suitability-derived mask will be used, otherwise, the mask of protection opportunities resulting from the LPPI/bare soil frequency clustering.
        :type protection_mask_type: str, optional
        :param exclude_non_built_areas: Should areas without any building coverage be excluded from the analysis, defaults to False. If True, only areas with some building coverage will be included in priorities, with non-built hexagons assigned nodata value.
        :type exclude_non_built_areas: bool, optional
        """
                
        self.logger.info("Running building solutions")
        
        # make base path to access required data
        base_path = os.path.join(
            self.scan_context.local_root, 
            self.scan_context.country, 
            self.scan_context.city
        )

        
        # create output path and assert the output path exists
        out_path = os.path.join(
            base_path,
            'output', 
            'building_solutions'
        )
        os.makedirs(out_path, exist_ok=True)

        # operate on the priorites shapefile written in priorities assessment
        # this saves us from rasterizing and or merging needs a bit downstream
        hexagons_path = os.path.join(
            base_path, 
            'output', 
            'priority', 
            f'{self.scan_context.city}_priorities.shp'
        )
        self.logger.info(f"Reading hexagons from {hexagons_path}")
        gdf_hexagons = gpd.read_file(hexagons_path)

        # determine building density 
        buildings_path = os.path.join(
            base_path, 
            'eo_indicators', 
            building_source
        )
        self.logger.info(f"Reading building data from {buildings_path}")

        ds_buildings = rioxarray.open_rasterio(buildings_path, masked=True)

        mtx_buildings = ds_buildings.values[0] 
        mtx_buildings = np.where(mtx_buildings > 0, 1, 0)  # make sure that we have the correct value for building footprints

        ds_buildings.values[0] = mtx_buildings  


        buildings_resolution = ds_buildings.rio.resolution()
        pix_area = abs(buildings_resolution[0]) * abs(buildings_resolution[1])

                
        # assert common projection in non-wgs84 (eo indicator layers should always be UTM)
        gdf_hexagons = gdf_hexagons.to_crs(ds_buildings.rio.crs)
        gdf_hexagons["hex_area"] = gdf_hexagons.geometry.area

        # determine building footprint cover share
        gdf_hexagons = make_exact_zonal_statistics(gdf_hexagons, ds_buildings, {'bldg_sum': 'sum'})
        gdf_hexagons['bldg_sum'] = gdf_hexagons['bldg_sum'] * pix_area  # convert from number of pixels to area in m2
        gdf_hexagons['bldg_dens'] = gdf_hexagons['bldg_sum'] / gdf_hexagons['hex_area']

        # determine share of protection and creation opportunities in each hexagon
        if protection_mask_type == "suitability":
            prot_file_name = 'nbs_protection.tif'
        elif protection_mask_type == "lppi":
            prot_file_name = 'nbs_protection_clustering.tif'


        nbs_protection_path = os.path.join(
            base_path, 
            'suitability', 
            prot_file_name
        )
        ds_protection = rioxarray.open_rasterio(nbs_protection_path, masked=True)
        protection_pixel_area = abs(ds_protection.rio.resolution()[0]) * abs(ds_protection.rio.resolution()[1])

        protection = ds_protection.values[0]  # extract array from xarray dataset, assuming single band raster with protection opportunities indicated by value of 1 and non-protection opportunities by value of 0        
        if protection_mask_type == "suitability":
            protection = np.where(protection == 1, 1, 0) # make sure that we have the correct value for protection
        elif protection_mask_type == "lppi":
            protection = np.where(protection == 1, 1, 0)  # make sure that we have the correct value for protection
        ds_protection.values[0] = protection  # update values in xarray dataset

        nbs_creation_path = os.path.join(
            base_path, 
            'output', 
            'mca', 
            'suitabilities.tif'
        )        
        ds_creation = rioxarray.open_rasterio(nbs_creation_path, masked=True)
        creation_pixel_area = abs(ds_creation.rio.resolution()[0]) * abs(ds_creation.rio.resolution()[1])

        # extract share of creation and protecton opportunities for hexagons
        gdf_hexagons = make_exact_zonal_statistics(gdf_hexagons, ds_protection, {'prot_sum' : 'sum'})
        gdf_hexagons['prot_sum'] = gdf_hexagons['prot_sum'] * protection_pixel_area  # convert from number of pixels to area in m2
        
        gdf_hexagons = make_exact_zonal_statistics(gdf_hexagons, ds_creation, {'creat_sum' :'sum'})
        gdf_hexagons['creat_sum'] = gdf_hexagons['creat_sum'] * creation_pixel_area  # convert from number of pixels to area in m2

        # sum opportunities and create the inverse
        gdf_hexagons['opp_sum'] = gdf_hexagons['prot_sum'] + gdf_hexagons['creat_sum']
        gdf_hexagons['opp_dens'] = gdf_hexagons['opp_sum'] / gdf_hexagons['hex_area']
        # adjust to rounding errors
        # replace all values in opp_dens that are above 1 (due to rounding errors) with 1
        gdf_hexagons.loc[gdf_hexagons['opp_dens'] > 1, 'opp_dens'] = 1

        gdf_hexagons['inv_opp'] = 1 - (gdf_hexagons['opp_dens'])

        # determine building opportunity score as product of building density and inverse of opportunity density
        gdf_hexagons['bldg_score'] = gdf_hexagons['bldg_dens'] * gdf_hexagons['inv_opp']

        # determine if heat and/or flooding are 
        # in the benefit parameters that are being assessed
        # these fields should already be contained in the priorities shapefile
        fields_to_run = []
        for scan_benefit in benefit_parameters:
            if scan_benefit.benefit is Benefits.UrbanHeat:
                fields_to_run.append('heat')
            if scan_benefit.benefit is Benefits.PluvialFlooding:
                fields_to_run.append('flood')
                
        # multiply with heat and/or flood priority
        for field in fields_to_run:

            # construct fieldnames
            out_score_fieldname = f'{field}_b_s'
            out_priority_fieldname = f'{field}_b_p'

            # multiply building score with priority score
            gdf_hexagons[out_score_fieldname] = gdf_hexagons[field] * gdf_hexagons['bldg_score']

            # replace NaN with 0 in output so that we can do the multiplication, these are areas where there is no heat priority, so it makes sense to have 0 here
            gdf_hexagons[out_score_fieldname] = gdf_hexagons[out_score_fieldname].fillna(0)
                        
            classifier = mapclassify.NaturalBreaks(gdf_hexagons[out_score_fieldname], k=jenks_classes)
            gdf_hexagons[out_priority_fieldname] = classifier.yb + 1 # add 1 to have the classes start at 1 instead of 0

            # mask non-built areas if specified
            # use a bit of tolerance to account for rounding errors
            if exclude_non_built_areas:
                gdf_hexagons.loc[gdf_hexagons['bldg_dens'] < 0.001, out_priority_fieldname] = np.nan


            # rasterize outputs, using buildings raster as template
            rasterize_df(
                polys=gdf_hexagons,
                templ=buildings_path,
                out_ras=os.path.join(out_path, f'building_solutions_{field}_binned.tif'),
                attr=out_priority_fieldname,
            )



        hex_out_path = os.path.join(
            out_path, 
            f'buildings.shp'
        )
        gdf_hexagons.to_file(hex_out_path)

    def archive_data(self):
        """
        Archive data to the cloud storage provider. This will upload all data in the local NBSOS data repository to the cloud storage provider.
        Note: This is only implemented for the GcsStorageHandler.
        
        """
        self.logger.info('Archiving data')
        self.storage_handler.archive_to_storage()

    def close(self, archive_log_to_cloud: bool = False):

        self.logger.info("Closing NBSOS scan")
        # now that we are done with the scan, take care of the logs
        self.save_logs(save_to_cloud=archive_log_to_cloud)

