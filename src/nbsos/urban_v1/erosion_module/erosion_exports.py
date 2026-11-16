from enum import Enum

class Lulc_Erosion_Weights(Enum):
    Urban = 10
    Cropland = 5.5 
    Forest = 10
    Shrubs = 7
    Grassland = 8
    Waterbodies = 1
    Wetlands = 7
    Bare = 4 # mean of sparse vegetation and bare soil as in irrigated

class Nbs_Erosion_Weights(Enum):
    BioretentionAreas = 1
    GreenCorridors = 1
    InlandWetlands = 1
    OpenGreenSpace = 2
    RiverRenaturation = 2
    TerracesAndSlopes = 3
    UrbanFarming = 2
    UrbanForest = 3
    RiverFloodplain = 1 
    BuildingSolutions = 1 