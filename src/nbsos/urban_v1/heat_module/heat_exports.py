from enum import Enum

class RegressionPredictors(Enum):
    BuiltCoverage = 'built_pct'
    CropCoverage = 'crops_pct'
    TreeCoverage = 'tree_pct'
    ShrubCoverage = 'shrub_pct'
    GrassCoverage = 'grass_pct'
    WaterCoverage = 'water_pct'
    BareCoverage = 'bare_pct'

class NbsTreeCoverage(Enum):
    BioretentionAreas = 0.05
    GreenCorridors = 1.00
    InlandWetlands = 0.50
    OpenGreenSpace = 0.30
    RiverRenaturation = 1.00
    TerracesAndSlopes = 0.20
    UrbanFarming = 0.05
    UrbanForest = 1.00
    RiverFloodplain = 0
    BuildingSolutions = 0