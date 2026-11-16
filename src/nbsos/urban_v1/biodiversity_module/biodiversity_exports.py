from typing import Dict, List
from enum import Enum

class ThreatSourceType:
    Land_Cover = 0
    Road_Barriers = 1

class ThreatDecayFunction:
    Linear = 0
    Exponential = 1

class ThreatProperties:

    source_type = None
    threat_id = None
    maximum_distance = None
    relative_weight = None
    decay_function = None

    def __init__(self, source_type, id, max_distance, weight, decay):
        self.source_type = source_type
        self.threat_id = id
        self.maximum_distance = max_distance
        self.relative_weight = weight
        self.decay_function = decay

    def get_decay_function(self) -> ThreatDecayFunction:
        return self.decay_function
    def get_maximum_distance(self) -> float:
        return self.maximum_distance
    def get_relative_weight(self) -> float:
        return self.relative_weight
    def get_value_for_selection(self) -> any:
        return self.threat_id
    def get_source_type(self) -> ThreatSourceType:
        return self.source_type

class VegetationLayer(Enum):
    Lower = 0
    Middle = 1
    Upper = 2


class HabitatThreat(Enum):
    CROPLAND = ThreatProperties(ThreatSourceType.Land_Cover, 2, 4.97, 0.65, ThreatDecayFunction.Linear)
    URBAN = ThreatProperties(ThreatSourceType.Land_Cover, 1, 5.45, 1, ThreatDecayFunction.Exponential)
    PRIMARY_ROAD = ThreatProperties(ThreatSourceType.Road_Barriers, 'primary', 2.37, 0.98, ThreatDecayFunction.Exponential)
    SECONDARY_ROAD = ThreatProperties(ThreatSourceType.Road_Barriers, 'secondary', 1.325, 0.69, ThreatDecayFunction.Linear)
    TERTIARY_ROAD = ThreatProperties(ThreatSourceType.Road_Barriers, 'tertiary', 0.573, 0.44, ThreatDecayFunction.Exponential)

class LandCoverClassVegetationProperties:
    
    has_lower_layer  = None
    has_middle_layer = None
    has_canopy_layer = None
    land_cover_class = None

    def __init__(self, class_value, has_lower_layer: bool, has_middle_layer: bool, has_upper_layer: bool):
        self.land_cover_class = class_value
        self.has_lower_layer = has_lower_layer
        self.has_middle_layer = has_middle_layer
        self.has_canopy_layer = has_upper_layer

    def get_land_cover_class(self) -> int:
        return self.land_cover_class

    def get_vegetation_layer_count(self) -> int:
        return int(self.has_lower_layer) + int(self.has_middle_layer) + int(self.has_canopy_layer)
    
    def has_vegetation_layer(self, layer: VegetationLayer) -> bool:
        if layer is VegetationLayer.Lower:
            return self.has_layer_lower()
        elif layer is VegetationLayer.Middle:
            return self.has_layer_middle()
        elif layer is VegetationLayer.Upper:
            return self.has_layer_upper()

    def has_layer_lower(self) -> bool:
        return self.has_lower_layer
    def has_layer_middle(self) -> bool:
        return self.has_middle_layer
    def has_layer_upper(self) -> bool:
        return self.has_canopy_layer


class PhysiognomicLandscapeType(Enum):
    Any = 0
    OpenGrasslandLike = 1 # lower layer OR # lower layer + upper layer
    ShrublandLike = 2 # lower layer + middle layer
    ForestLike    = 3 # lower layer + middle layer + upper layer


    def get_landscape_type(habitat_type: LandCoverClassVegetationProperties):
        if habitat_type.has_layer_lower() and not habitat_type.has_layer_middle() and not habitat_type.has_layer_upper():
            return PhysiognomicLandscapeType.OpenGrasslandLike
        elif habitat_type.has_layer_lower() and habitat_type.has_layer_middle() and not habitat_type.has_layer_upper():
            return PhysiognomicLandscapeType.ShrublandLike
        elif habitat_type.has_layer_lower() and not habitat_type.has_layer_middle() and habitat_type.has_layer_upper():
            return PhysiognomicLandscapeType.OpenGrasslandLike
        elif habitat_type.has_layer_lower() and habitat_type.has_layer_middle() and habitat_type.has_layer_upper():
            return PhysiognomicLandscapeType.ForestLike
         

class LandCoverClassProperties:

    habitat_suitability: float = None
    habitat_sensitivities: Dict[HabitatThreat, float] = None
    vegetation_properties: LandCoverClassVegetationProperties = None

    def __init__(self, Land_cover_class: int, habitat_suitability: float, has_understory_layer: bool, has_midstory_layer: bool, has_canopy_layer: bool, dict_of_sensitivities):
        self.vegetation_properties = LandCoverClassVegetationProperties(Land_cover_class, has_understory_layer, has_midstory_layer, has_canopy_layer)
        self.habitat_suitability = habitat_suitability        
        self.habitat_sensitivities = dict_of_sensitivities
        
    def get_vegetation_properties(self) -> LandCoverClassVegetationProperties:
        return self.vegetation_properties
    
    def get_land_cover_class_value(self) -> int:
        return self.vegetation_properties.get_land_cover_class()

    def get_sensitivity_to_threat(self, threat_type: HabitatThreat) -> float:
        if not threat_type in self.habitat_sensitivities.keys():
            return 0
        return self.habitat_sensitivities[threat_type]
    
    def get_habitat_suitability(self):
        return self.habitat_suitability
    
    
class LandCover(Enum):

    URBAN = LandCoverClassProperties(1, 0.0306, False, False, False, {HabitatThreat.CROPLAND : 0, HabitatThreat.URBAN : 0, HabitatThreat.PRIMARY_ROAD : 0.165, HabitatThreat.SECONDARY_ROAD : 0.17, HabitatThreat.TERTIARY_ROAD : 0.182}) 
    CROPLAND = LandCoverClassProperties(2, 0.525, False, False, False, {HabitatThreat.CROPLAND : 0, HabitatThreat.URBAN : 0.47, HabitatThreat.PRIMARY_ROAD : 0.389, HabitatThreat.SECONDARY_ROAD : 0.55, HabitatThreat.TERTIARY_ROAD : 0.295}) 
    FOREST = LandCoverClassProperties(3, 0.938, True, True, True, {HabitatThreat.CROPLAND : 0.59, HabitatThreat.URBAN : 0.7966, HabitatThreat.PRIMARY_ROAD : 0.715, HabitatThreat.SECONDARY_ROAD : 0.608, HabitatThreat.TERTIARY_ROAD : 0.459}) 
    SHRUBS = LandCoverClassProperties(4, 0.712, True, True, False, {HabitatThreat.CROPLAND : 0.675, HabitatThreat.URBAN : 0.749, HabitatThreat.PRIMARY_ROAD : 0.593, HabitatThreat.SECONDARY_ROAD : 0.636, HabitatThreat.TERTIARY_ROAD : 0.416})
    GRASSLAND = LandCoverClassProperties(5, 0.596, True, False, False, {HabitatThreat.CROPLAND : 0.458, HabitatThreat.URBAN : 0.597, HabitatThreat.PRIMARY_ROAD : 0.567, HabitatThreat.SECONDARY_ROAD : 0.514, HabitatThreat.TERTIARY_ROAD : 0.357})
    WATERBODIES = LandCoverClassProperties(6, 0.945, False, False, False, {HabitatThreat.CROPLAND : 0.686, HabitatThreat.URBAN : 0.84, HabitatThreat.PRIMARY_ROAD : 0.6, HabitatThreat.SECONDARY_ROAD : 0.593, HabitatThreat.TERTIARY_ROAD : 0.446})
    WETLANDS = LandCoverClassProperties(7, 0.784, True, True, False, {HabitatThreat.CROPLAND : 0.595, HabitatThreat.URBAN : 0.915, HabitatThreat.PRIMARY_ROAD : 0.582, HabitatThreat.SECONDARY_ROAD : 0.573, HabitatThreat.TERTIARY_ROAD : 0.496})
    BARE = LandCoverClassProperties(8, 0.06, False, False, False, {HabitatThreat.CROPLAND : 0.008, HabitatThreat.URBAN : 0.083, HabitatThreat.PRIMARY_ROAD : 0.04, HabitatThreat.SECONDARY_ROAD : 0, HabitatThreat.TERTIARY_ROAD : 0})

    def get_land_cover_class_value(self) -> int:
        return self.value.get_land_cover_class_value()
    def get_vegetation_layer_count(self) -> int:
        return self.value.get_vegetation_properties().get_vegetation_layer_count()

class OpportunityType(Enum):
    
    bioretention_area_creation = LandCoverClassVegetationProperties(10, True, True, False)
    green_corridors_creation = LandCoverClassVegetationProperties(11, True, True, True) 
    inland_wetland_creation = LandCoverClassVegetationProperties(12, True, True, False) 
    open_green_space_creation = LandCoverClassVegetationProperties(13, True, False, True) 
    river_renaturation_creation = LandCoverClassVegetationProperties(14, True, True, True) 
    terraces_and_slopes_creation = LandCoverClassVegetationProperties(15, True, False, False) 
    urban_farming_creation = LandCoverClassVegetationProperties(16, True, False, False) 
    urban_forest_creation = LandCoverClassVegetationProperties(17, True, True, True) 
    river_floodplain_creation = LandCoverClassVegetationProperties(18, True, True, True) 
    building_solutions_creation = LandCoverClassVegetationProperties(19, False, False, False)
    
    def get_land_cover_class_value(self):
        return self.value.get_land_cover_class()    
    def get_vegetation_layer_count(self):
        return self.value.get_vegetation_layer_count()
    
    
    

class Metrics(Enum):
    
    EffectiveMeshSize = "effective mesh size"
    ScaledEffectiveMeshSize = "scaled effective mesh size"
    Coherence = "coherence"
    ScaledCoherence = "scaled coherence"
    LandscapeDivision = "landscape division"
    SplittingIndex = "splitting index"
    PatchCount = "patch count"
    ConnectedEffectiveMeshSize = "connected effective mesh size"
    ScaledConnectedEffectiveMeshSize = "scaled connected effective mesh size"
    ComponentCount = "component count"
    