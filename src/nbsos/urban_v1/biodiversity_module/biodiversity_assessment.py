
from .biodiversity_exports import VegetationLayer, PhysiognomicLandscapeType, LandCover, OpportunityType, Metrics, LandCoverClassVegetationProperties
from .biodiversity_results import MetricsContainer, RasterContainer, RasterSet, MetricsSet, CreationOpportunityImpactsContainer, HabitatQualityContainer, ResultCollection, CreationOpportunityImpactsSet
from .biodiversity_indicators import BiodiversityIndicators
from .biodiversity_spatial import Spatial
from .biodiversity_logging import biodiversity_logger

from typing import List, Tuple, Dict
import numpy as np
import geopandas as gpd
import rasterio
import rasterio.features
import rasterio.mask
import concurrent.futures
from enum import Enum
from skimage.measure import label, regionprops
import concurrent.futures
from tqdm import tqdm
from tqdm.auto import tqdm as tqauto
import multiprocessing as mp
import functools
import sys

def log_function(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as exc:
            print(args, kwargs, repr(exc))

    return wrapper


def _hexagon_worker(hex_label, land_cover_matrix, mtx_intervention, mtx_hexagons,
                    mtx_rasterized_roads, creation_opportunity_type, scenario_habitat_values,
                    landscape_types_to_assess, resolution, min_pixels,
                    require_optimization_criterion, connectivity_threshold):
    current_land_cover = land_cover_matrix.copy()
    current_land_cover[((mtx_intervention == 1) & (mtx_hexagons == hex_label))] = creation_opportunity_type.get_land_cover_class_value()
    current_land_cover, current_masked_roads = AssessmentContainer._burn_roads_into_land_cover(current_land_cover, mtx_rasterized_roads.copy(), scenario_habitat_values)
    for landscape_type in landscape_types_to_assess:
        _, current_stats = AssessmentContainer._assess_landscape_properties(current_land_cover, current_masked_roads, landscape_type, scenario_habitat_values, resolution, require_optimization_criterion, connectivity_threshold, min_pixels)
    return landscape_type, hex_label, current_stats


class AssessmentContainer:
    """This class encapsulates methods and data for the assessment of land-use scenario, including baseline.
    """
    resolution: Tuple[float,float] = None
    creation_opportunity_type = None
    scenario_name: str = None
    scenario_habitat_values: List[LandCoverClassVegetationProperties] = None
    landscape_types_to_assess: List[PhysiognomicLandscapeType] = None
    mtx_land_cover: np.ndarray = None
    mtx_intervention: np.ndarray = None
    
    require_optimization_criterion: Metrics = None
    connectivity_threshold: float = 100

    # barriers
    gdf_road_barriers: gpd.GeoDataFrame = None
    geo_transform = None
    
    # habitat quality
    assess_habitat_quality: bool = None    
    threats = None
    sensitivity = None 

    def __init__(
            self, mtx_land_cover_baseline: np.ndarray, 
            mtx_hexagons: np.ndarray, 
            nbs_scenario: Tuple[OpportunityType, np.ndarray], 
            habitat_types: List[LandCover], 
            landscape_types: List[PhysiognomicLandscapeType], 
            road_barriers: gpd.GeoDataFrame, 
            sensitivity, 
            threats, 
            resolution: Tuple[float,float], 
            transform, 
            require_optimization_criterion: Metrics = Metrics.ScaledEffectiveMeshSize,
            assess_habitat_quality: bool = True,
            connectivity_threshold: float = 100,
            min_pixels: int = None
            ) -> None:
        
        # as we modify land cover per scenario, make a copy of baseline
        # this copy is subsequently modified to integrate the nbs intervention pixels
        self.mtx_land_cover = mtx_land_cover_baseline.copy()
        self.mtx_hexagons = mtx_hexagons
        self.mtx_rasterized_roads = None
        self.require_optimization_criterion = require_optimization_criterion
        self.connectivity_threshold = connectivity_threshold
        self.sensitivity = sensitivity
        self.threats = threats        
        self.assess_habitat_quality = assess_habitat_quality
        self.resolution = resolution
        self.geo_transform = transform        

        self.min_pixels = min_pixels

        self.results_connectivity = {}



        # note that if nbs_scenario is None, we consider the scenario to be baseline
        self.name = 'baseline' if nbs_scenario is None else str(nbs_scenario[0].name)
        
        # store physiognomic landscape types to assess
        self.landscape_types_to_assess = landscape_types

        # store baseline habitat values
        self.scenario_habitat_values = [f.value.get_vegetation_properties() for f in habitat_types]
               
        # nbs_scenario is a Tuple[ScenarioNbsType, np.ndarray]
        if nbs_scenario is not None:                   
            self.creation_opportunity_type = nbs_scenario[0]   
            self.mtx_intervention = nbs_scenario[1]                       
            # append  scenario nbs class value to list of values            
            self.scenario_habitat_values.append(nbs_scenario[0].value)

        # take care of road barriers        
        if road_barriers is not None:
            # store original geodataframe for later re-use to compute threat layers
            self.gdf_road_barriers = road_barriers
            # make a matrix of rasterized road shapes and burn roads into scenario
            shapes = [geom for geom in self.gdf_road_barriers.geometry]
            # make raster of roads to be re-used across scenarios
            self.mtx_rasterized_roads = rasterio.features.rasterize(shapes=shapes, fill=0, default_value=1, all_touched=False, out_shape=self.land_cover_matrix.shape, dtype=np.float32, transform=self.geo_transform)    
            

    @staticmethod
    def _burn_roads_into_land_cover(mtx_land_cover: np.ndarray, mtx_roads:np.ndarray, habitat_types: List[LandCoverClassVegetationProperties]) -> np.ndarray:
        
        # mask pixels to non-habitat (infrastructure=9)            
        # however, burn roads only at pixel locations where roads do not intersect three-storey vegetation            
        # iterate over habitat types in scenario and burn over roads when considered not-fragmenting (canopy layer present)
        habitat_type: LandCoverClassVegetationProperties
        for habitat_type in habitat_types:
            if habitat_type.get_vegetation_layer_count() == 3:
                # "un-burn" roads across land cover where this class is present
                mtx_roads[mtx_land_cover == habitat_type.get_land_cover_class()] = 0
        
        # burn remaining roads as barrier into land cover matrix
        mtx_land_cover[((mtx_roads == 1) & (mtx_land_cover > 0))] = 9 
        return mtx_land_cover, mtx_roads

            
      
    def _assess_habitat_quality(self):
        # return results as container
        habitat_quality_results = HabitatQualityContainer()
        # assess habitat quality
        habitat_quality_results.habitat_suitability_matrix, habitat_quality_results.degradation_matrix, habitat_quality_results.habitat_quality_matrix = BiodiversityIndicators.habitat_quality(self.land_cover_matrix, self.gdf_road_barriers, self.sensitivity, self.threats, self.resolution, self.geo_transform)                
        return habitat_quality_results        

    @staticmethod
    def _assess_landscape_properties(
                                     mtx_land_cover: np.ndarray, 
                                     masked_roads: np.ndarray, 
                                     landscape_type: PhysiognomicLandscapeType, 
                                     habitat_values: List[LandCoverClassVegetationProperties],
                                     resolution: Tuple[float, float],
                                     require_optimization_criterion: Metrics,
                                     connectivity_threshold: float,
                                     min_pixels: int = None) -> Tuple[RasterContainer, MetricsContainer]:

        current_layers = RasterContainer(landscape_type, habitat_values)
        current_stats = MetricsContainer()

        current_stats.patch_count, current_stats.total_habitat_area, current_layers.patch_matrix, current_stats.patch_dataframe, mtx_patch_sizes = Spatial.determine_patches(mtx_land_cover, current_layers.habitats_in_layer, resolution=resolution, min_pixels=min_pixels)

        # non-connected components assessment
        indicator_method = None
        if require_optimization_criterion is Metrics.Coherence or require_optimization_criterion is Metrics.ScaledCoherence or require_optimization_criterion is Metrics.EffectiveMeshSize or require_optimization_criterion is Metrics.ScaledEffectiveMeshSize:

            if require_optimization_criterion is Metrics.Coherence or require_optimization_criterion is Metrics.ScaledCoherence:
                indicator_method = BiodiversityIndicators.coherence
            elif require_optimization_criterion is Metrics.EffectiveMeshSize or require_optimization_criterion is Metrics.ScaledEffectiveMeshSize:
                indicator_method = BiodiversityIndicators.effective_mesh_size
                            
            # compute coherence or effective mesh size at patch level
            current_layers.indicator_matrix = indicator_method(mtx_patch_sizes, current_stats.total_habitat_area)             
            # determine statistics according to reference paper in tabular form from table data on current scenario
            current_stats.patch_dataframe, current_stats.metrics = BiodiversityIndicators.determine_metrics(current_stats.patch_dataframe, current_stats.total_habitat_area)

        # adding connected components 
        elif require_optimization_criterion is Metrics.ConnectedEffectiveMeshSize or require_optimization_criterion is Metrics.ScaledConnectedEffectiveMeshSize:
            
            # first, determine connected components
            dilation_iterations = int((connectivity_threshold / resolution[0]))
            mtx_dilated_edges, component_count, current_layers.component_labels = Spatial.edge_dilation(current_layers.patch_matrix, masked_roads, dilation_iterations)
            current_layers.component_labels, current_layers.indicator_matrix, current_layers.component_patch_count, current_stats.patch_dataframe = BiodiversityIndicators.connected_effective_mesh_size(current_layers.patch_matrix, current_stats.patch_dataframe, mtx_patch_sizes, current_layers.component_labels, current_stats.total_habitat_area)
            current_stats.patch_dataframe, current_stats.component_dataframe, current_stats.metrics = BiodiversityIndicators.determine_metrics_for_connected_components(current_stats.patch_dataframe, current_stats.total_habitat_area)


        # determine the scaled indicator. actual scaling is only done if required, i.e., if the optimization criterion is scaled
        if require_optimization_criterion is Metrics.ScaledCoherence or require_optimization_criterion is Metrics.ScaledEffectiveMeshSize or require_optimization_criterion is Metrics.ScaledConnectedEffectiveMeshSize:
            # scale by vegetation structure
            scale_factor_matrix = BiodiversityIndicators.get_vegetation_layer_weight_mask(mtx_land_cover, habitat_values)
            # scale by vegetation structure
            current_layers.scaled_indicator_matrix = BiodiversityIndicators.vegetation_layer_scaling(current_layers.indicator_matrix, scale_factor_matrix)
            
            if require_optimization_criterion is Metrics.ScaledCoherence or require_optimization_criterion is Metrics.ScaledEffectiveMeshSize:
                # determine mean of scaled indicator per patch as substitute global statistics value 
                region_stats = BiodiversityIndicators.mean_of_patch(current_layers.patch_matrix, current_layers.scaled_indicator_matrix)
                current_stats.metrics[require_optimization_criterion] = region_stats['intensity_mean'].sum()
            
            elif require_optimization_criterion is Metrics.ScaledConnectedEffectiveMeshSize:
                # determine mean of scaled component sizes. If the complete component is a single land-use, component size is equal to unscaled.
                region_stats = BiodiversityIndicators.mean_of_patch(current_layers.component_labels, current_layers.scaled_indicator_matrix)
                current_stats.metrics[require_optimization_criterion] = region_stats['intensity_mean'].sum()

        elif require_optimization_criterion is Metrics.Coherence or require_optimization_criterion is Metrics.EffectiveMeshSize or require_optimization_criterion is Metrics.ConnectedEffectiveMeshSize:
            # since the module will write out the scaled indicator matrix, and use that also as a reference, assert that this matrix will be available
            # this is essentially then the same matrix that is stored.
            current_layers.scaled_indicator_matrix = current_layers.indicator_matrix

        return current_layers, current_stats



    @log_function
    def assess(self) -> ResultCollection:
        
        # Essentially, assess shall return the relevant results to be stored in the module,
        # so that after completion of the assessment, this container can be "thrown away"

        result_container = ResultCollection(self.name)
        
        # the following are basically only returned for baseline
        raster_results: RasterSet = RasterSet() if self.creation_opportunity_type is None else None
        metrics_results: MetricsSet = MetricsSet() if self.creation_opportunity_type is None else None
        habitat_quality_results: HabitatQualityContainer = self._assess_habitat_quality() if (self.assess_habitat_quality and self.sensitivity is not None and self.threats is not None) else None        
        # whereas the following is only returned for nbs intervention scenarios, as we do not want to store thousands of individual scenarios
        nbs_impacts: CreationOpportunityImpactsSet = CreationOpportunityImpactsSet(self.landscape_types_to_assess) if self.creation_opportunity_type is not None else None

        # this iterates over physiognomic landscape types, if more specified than any. 
        # Any would be overall connectivity including all habitats, comparable to prior methodology, and is the default. 
        if self.creation_opportunity_type is None:
            biodiversity_logger.info("Assessing baseline")
            # this is an assessment according to baseline.
            # therefore, we do not need to burn hexagons into the baseline land cover layer
            # simply burn roads
            mtx_land_cover, masked_roads = self._burn_roads_into_land_cover(self.land_cover_matrix.copy(), self.mtx_rasterized_roads.copy(), self.scenario_habitat_values)

            for landscape_type in self.landscape_types_to_assess:
                biodiversity_logger.info(f"Assessing landscape of type {landscape_type.name}")

                current_layers, current_stats = self._assess_landscape_properties(mtx_land_cover, masked_roads, landscape_type, self.scenario_habitat_values, self.resolution, self.require_optimization_criterion, self.connectivity_threshold, self.min_pixels)
                biodiversity_logger.info(f"{self.name} with {current_stats.patch_count} habitat patches and a total habitat area of {current_stats.total_habitat_area} m² in layer {landscape_type.name}") 

                raster_results.insert_result(landscape_type, current_layers)
                metrics_results.insert_result(landscape_type, current_stats)
            
        else:

            # this is an assessment of a scenario; here, we would need to iterate over the hexagons, integrating per hexagon the proposed measures
            # in this case, we will not write connectivity results, but throw away the resultcontainer each time
            # we are only interested in getting the metrics that will include the indicator for optimization
            # and we write this into a dict
           
            # get unique values of hexagons, excluding background
            # unique_hexagon_labels = np.unique(self.mtx_hexagons[self.mtx_hexagons > 0]) # deprecate after testing
            # however, get only those hexagons that are of interest for the current scenario, i.e., that intersect at least 1 cell of nbs creation presence
            intersecting_grid_cells = (np.isin(self.mtx_intervention, [1], invert=False) & (self.mtx_hexagons > 0))
            unique_hexagon_labels = np.unique(self.mtx_hexagons[intersecting_grid_cells])

            # iterate over hexagon labels 
            # determine later if we can really parallelize this. the fear is that too many rasters will be written into memory for parallelization
                    
            worker = functools.partial(
                _hexagon_worker,
                land_cover_matrix=self.land_cover_matrix,
                mtx_intervention=self.mtx_intervention,
                mtx_hexagons=self.mtx_hexagons,
                mtx_rasterized_roads=self.mtx_rasterized_roads,
                creation_opportunity_type=self.creation_opportunity_type,
                scenario_habitat_values=self.scenario_habitat_values,
                landscape_types_to_assess=self.landscape_types_to_assess,
                resolution=self.resolution,
                min_pixels=self.min_pixels,
                require_optimization_criterion=self.require_optimization_criterion,
                connectivity_threshold=self.connectivity_threshold
            )

            with tqauto(total=len(unique_hexagon_labels), desc="Iterating hexagons") as p:
                with concurrent.futures.ThreadPoolExecutor(max_workers=mp.cpu_count()) as executor:
                    futures = {executor.submit(worker, label): label for label in unique_hexagon_labels}
                    for future in concurrent.futures.as_completed(futures):
                        try:
                            result = future.result()
                            nbs_impacts.insert_into_results(result[0], result[1], result[2])
                        except Exception as e:
                            print(f"Generated an exception for hexagon {futures[future]}: {e}")
                        finally:
                            p.update(1)


            


        biodiversity_logger.info(f"Assessment of scenario {self.name} completed.")   
                
        result_container.rasters = raster_results
        result_container.metrics = metrics_results
        result_container.habitat_quality = habitat_quality_results 
        result_container.impacts = nbs_impacts
        return result_container

    

    @property
    def name(self) -> str:
        return self.scenario_name
    @name.setter
    def name(self, value) -> None:
        self.scenario_name = value
    
    @property
    def land_cover_matrix(self) -> np.ndarray:
        return self.mtx_land_cover
    @land_cover_matrix.setter
    def land_cover_matrix(self, value) -> None:
        self.mtx_land_cover = value
    
    @property
    def intervention_matrix(self) -> np.ndarray:
        return self.mtx_intervention
    @intervention_matrix.setter
    def intervention_matrix(self, value) -> None:
        self.mtx_intervention = value

