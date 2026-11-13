from .biodiversity_exports import VegetationLayer, PhysiognomicLandscapeType, LandCover, OpportunityType, Metrics, LandCoverClassVegetationProperties
from .biodiversity_indicators import BiodiversityIndicators
from .biodiversity_spatial import Spatial
from .biodiversity_logging import biodiversity_logger

from typing import List, Tuple, Dict
import numpy as np
import geopandas as gpd
import pandas as pd

# The following classes are implemented to make handling of the various outputs created across the biodiversity module easier.
# Generally, there are containers to store the following outputs: rasters, metrics, habitat quality, and creation opportunities
# For rasters, metrics, and impacts, there is typically a container per physiognomic landscape type assessed. An exception to this is the container for habitat quality.  
# Most containers are in turn encapsulated in sets. There are sets for rasters, metrics, and impacts. 
# Finally, all of these different results are in turn handled by a collection, that encapsulates all of these sets in a convenient way.


class MetricsContainer:
    """This class serves as container to store metrics determined for a specific scenario, i.e., baseline, or a creation opportunity. 
    """
    total_patch_count : int = None
    df_of_patches: pd.DataFrame = None
    df_of_components: pd.DataFrame = None
    dict_of_metrics = None
    habitat_area: float = None

    @property
    def metrics(self) -> Dict[any, float]:
        return self.dict_of_metrics
    @metrics.setter
    def metrics(self, value) -> None:
        self.dict_of_metrics = value

    @property
    def patch_dataframe(self) -> pd.DataFrame:
        return self.df_of_patches
    @patch_dataframe.setter
    def patch_dataframe(self, value) -> None:
        self.df_of_patches = value  

    @property 
    def component_dataframe(self) -> pd.DataFrame:
        return self.df_of_components
    @component_dataframe.setter 
    def component_dataframe(self, value) -> None:
        self.df_of_components = value

    @property
    def patch_count(self) -> int:
        return self.total_patch_count
    @patch_count.setter
    def patch_count(self, value) -> None:
        self.total_patch_count = value

    @property 
    def total_habitat_area(self) -> float:
        return self.habitat_area
    @total_habitat_area.setter
    def total_habitat_area(self, value) -> None:
        self.habitat_area = value




class RasterContainer:
    """This class serves as a container to store spatially explicit results determined for a specific scenario, i.e., baseline, or a creation opportunity. 
    """

    habitat_types: List[LandCoverClassVegetationProperties] = None
    mtx_patches: np.ndarray = None       
    mtx_indicator: np.ndarray = None
    mtx_indicator_scaled: np.ndarray = None
    mtx_connected_component_labels: np.ndarray = None
    mtx_component_patch_count: np.ndarray = None

    def __init__(self, landscape_type: PhysiognomicLandscapeType, habitat_types: List[LandCoverClassVegetationProperties]):
        self.habitats_in_layer = BiodiversityIndicators.get_habitats_of_physiognomic_type(landscape_type, habitat_types)

    @property
    def habitats_in_layer(self) -> List[LandCoverClassVegetationProperties]:
        return self.habitat_types
    @habitats_in_layer.setter
    def habitats_in_layer(self, value) -> None:
        self.habitat_types = value

    @property
    def scaled_indicator_matrix(self) -> np.ndarray:
        return self.mtx_indicator_scaled
    @scaled_indicator_matrix.setter
    def scaled_indicator_matrix(self, value) -> None:
        self.mtx_indicator_scaled = value

    @property
    def indicator_matrix(self) -> np.ndarray:
        return self.mtx_indicator
    @indicator_matrix.setter
    def indicator_matrix(self, value) -> None:
        self.mtx_indicator = value

    @property
    def patch_matrix(self) -> np.ndarray:
        return self.mtx_patches
    @patch_matrix.setter
    def patch_matrix(self, value) -> None:
        self.mtx_patches = value

    @property
    def component_labels(self) -> np.ndarray:
        return self.mtx_connected_component_labels
    @component_labels.setter
    def component_labels(self, value) -> None:
        self.mtx_connected_component_labels = value

    @property
    def component_patch_count(self) -> np.ndarray:
        return self.mtx_component_patch_count
    @component_patch_count.setter
    def component_patch_count(self, value) -> None:
        self.mtx_component_patch_count = value

class CreationOpportunityImpactsContainer:
    """This class serves as a container to store the impacts of creation opportunity compared to baseline, for each hexagon of the AOI.
    """

    result_items: Dict[int, MetricsContainer] = None

    def __init__(self):
        self.result_items = {}

    def insert_result(self, hexagon_label: int, result: MetricsContainer) -> None:
        self.result_items[hexagon_label] = result

    def get_results(self) -> Dict[int, MetricsContainer]:
        return self.result_items
    
    def get_metric_value_for_hexagon(self, hexagon_label, metric: Metrics) -> float:
        return self.result_items[hexagon_label].metrics.get(metric, None)



class RasterSet:
    """This class serves as a manager to handle raster containers, as assessed for each physiognomic landscape type.
    """
    result_items: Dict[PhysiognomicLandscapeType, RasterContainer] = None

    def __init__(self):
        self.result_items = {}

    def insert_result(self, landscape_type, results):
        self.result_items[landscape_type] = results
    
    def get_results(self) -> Dict[PhysiognomicLandscapeType, RasterContainer]:
        return self.result_items
    
    def get_result_for_type(self, landscape_type: PhysiognomicLandscapeType) -> RasterContainer:
        return self.result_items.get(landscape_type, None)
    
class MetricsSet:
    """This class serves as a manager to handle metrics containers, as assessed for each physiognomic landscape type.
    """
    result_items: Dict[PhysiognomicLandscapeType, MetricsContainer] = None

    def __init__(self):
        self.result_items = {}

    def insert_result(self, landscape_type, results):
        self.result_items[landscape_type] = results
    
    def get_results(self) -> Dict[PhysiognomicLandscapeType, MetricsContainer]:
        return self.result_items

    def get_result_for_type(self, landscape_type: PhysiognomicLandscapeType) -> MetricsContainer:
        return self.result_items.get(landscape_type, None)
    
    def print_summary(self):        
        for landscape_type, res in self.get_results().items():                  
           print(f"Summary for {landscape_type.name}:")
           for metric, value in res.metrics.items():
               print(f"{metric.value} = {value}")


class CreationOpportunityImpactsSet:
    """This class serves as a manager to handle creation opportunity impacts containers, as assessed for each pyhsiognomic landscape type. 
    """
    result_items: Dict[PhysiognomicLandscapeType, CreationOpportunityImpactsContainer] = None

    def __init__(self, landscape_types_to_assess: List[PhysiognomicLandscapeType]):
        self.result_items = {}
        for landscape_type in landscape_types_to_assess:
            self.result_items[landscape_type] = CreationOpportunityImpactsContainer()

    def insert_result(self, landscape_type: PhysiognomicLandscapeType, impact_results: CreationOpportunityImpactsContainer) -> None:
        self.result_items[landscape_type] = impact_results

    def insert_into_results(self, landscape_type: PhysiognomicLandscapeType, hexagon_label: int, result: MetricsContainer) -> None:
        self.result_items[landscape_type].insert_result(hexagon_label, result)

    def get_results(self) -> Dict[PhysiognomicLandscapeType, CreationOpportunityImpactsContainer]:
        return self.result_items

    def get_result_for_type(self, landscape_type: PhysiognomicLandscapeType) -> CreationOpportunityImpactsContainer:
        return self.result_items.get(landscape_type, None)

class HabitatQualityContainer:
    
    mtx_suitability: np.ndarray = None
    mtx_degradation: np.ndarray = None
    mtx_quality: np.ndarray = None

    @property
    def habitat_suitability_matrix(self) -> np.ndarray:
        return self.mtx_suitability
    @habitat_suitability_matrix.setter
    def habitat_suitability_matrix(self, value) -> None:
        self.mtx_suitability = value

    @property
    def degradation_matrix(self) -> np.ndarray:
        return self.mtx_degradation
    @degradation_matrix.setter
    def degradation_matrix(self, value) -> None:
        self.mtx_degradation = value

    @property
    def habitat_quality_matrix(self) -> np.ndarray:
        return self.mtx_quality
    @habitat_quality_matrix.setter
    def habitat_quality_matrix(self, value) -> None:
        self.mtx_quality = value


class ResultCollection:
    """This class serves as a container for all the results returned. 
    """

    scenario_name: str = None

    set_of_rasters: RasterSet = None
    set_of_metrics: MetricsSet = None
    set_habitat_quality: HabitatQualityContainer = None
    set_of_impacts: CreationOpportunityImpactsSet = None

    def __init__(self, scenario_name: str):
        self.scenario_name = scenario_name

    # getter and setter methods for retrieving data in container
    @property 
    def rasters(self) -> RasterSet:
        return self.set_of_rasters
    @rasters.setter
    def rasters(self, value) -> None:
        self.set_of_rasters = value

    @property
    def metrics(self) -> MetricsSet:
        return self.set_of_metrics
    @metrics.setter
    def metrics(self, value) -> None:
        self.set_of_metrics = value

    @property 
    def habitat_quality(self) -> HabitatQualityContainer:
        return self.set_habitat_quality
    @habitat_quality.setter
    def habitat_quality(self, value) -> None:
        self.set_habitat_quality = value

    @property
    def impacts(self) -> CreationOpportunityImpactsSet:
        return self.set_of_impacts
    @impacts.setter
    def impacts(self, value) -> None:
        self.set_of_impacts = value

    @property
    def name(self) -> str:
        return self.scenario_name
    
    def contains_habitat_quality(self) -> bool:
        return True if self.set_habitat_quality is not None else False
        
 

    

