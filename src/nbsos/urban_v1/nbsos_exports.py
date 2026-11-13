from enum import Enum
from typing import Dict


class NbsTypes(Enum):
    """These are nbs types with names (value[0]) aligned to suitability str used, and proposing common landcover class values to be used across scan benefit modules
    """
    BioretentionAreas = 'bioretention_area', 10
    GreenCorridors = 'green_corridors', 11
    InlandWetlands = 'inland_wetland', 12
    OpenGreenSpace = 'open_green_space', 13
    RiverRenaturation = 'river_renaturation', 14
    TerracesAndSlopes = 'terraces_and_slopes', 15
    UrbanFarming = 'urban_farming', 16
    UrbanForest = 'urban_forest', 17

    def get_name(self):
        return self.value[0]

    def get_land_cover_class_value(self):
        return self.value[1]



class Benefits(Enum):
    """Benefits supported by the NBSOS.
    """
    PluvialFlooding = 'pluvial flooding'
    Recreation = 'recreation'
    UrbanHeat = 'urban heat mitigation'
    Biodiversity = 'biodiversity'
    SocialCohesion = 'social cohesion'
    AvoidedErosion = 'avoided erosion'
    JobsCreation = 'jobs creation'

    def get_name(self):
        return self.value

class ScanBenefit():
    """Define benefits to be assessed in a scan, and parameters for the assessment, as needed.
    """

    benefit_type: Benefits
    params: Dict[str,any]

    def __init__(self, benefit_type: Benefits, params: Dict[str,any] = {}):
        """Instantiate the ScanBenefit class to define a benefit to be assessed in a scan.

        :param benefit_type: The type of benefit to be assessed.
        :type benefit_type: Benefits
        :param params: Parameters for assessing the benefit, defaults to {}.
        :type params: Dict[str,any], optional
        """
        self.benefit_type = benefit_type
        self.params = params

    @property
    def benefit(self) -> Benefits:
        """Type of benefit to be assessed.

        :return: Benefit type.
        :rtype: Benefits
        """
        return self.benefit_type
    
    @property
    def parameters(self) -> Dict[str,any]:
        """Parameters for assessing the benefit.

        :return: User-defined parameters for the benefit assessment, overriding any default parameters. 
        :rtype: Dict[str,any]
        """
        return self.params

class LandCoverClass(Enum):

    Urban = 1
    Cropland = 2 
    Forest = 3
    Shrubs = 4
    Grassland = 5 
    Waterbodies = 6
    Wetlands = 7
    Bare = 8
