# Running the NBSOS

# nbsos-specific imports
# this assumes that urban_v1 is installed as a local package (pip install -e .) in src dir
# requirements.txt in the urban_v1 folder should list all required packages.
from nbsos.urban_v1 import Scan, ScanContext, DirectDownloadHandler, GcsHandler, NbsTypes, Benefits, ScanBenefit

# other imports
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

# General approach
# * Setup a scan context and other environment parameters;
# * Instantiate a new scan;
# * Retrieve and prepare flood hazard data (if required);
# * Download required data from Google Earth Engine, OSM, etc.;
# * Run suitability model
# * Run cluster exploration
# * Run clustering to create creation opportunity rasters and protection mask(s);
# * Run creation benefit models;
# * Run protection benefit models;
# * Run priority area identification;
# * Run multicriteria analysis;
# * If building solutions are relevant for the case, run building solutions method.


# define data and other parameters of scan, that will be used to instantiate the scan context
local_root = "path"
country_name = "name"
city_name = "name"
aoi_path = "path to wgs84 shape"
location_is_in_africa = False  # set as required

# set scan context
scan_context: ScanContext = ScanContext()
scan_context.local_root = local_root
scan_context.country = country_name
scan_context.city = city_name
scan_context.local_area_of_interest_file = aoi_path
scan_context.location_in_africa = location_is_in_africa


# To assess pluvial flood benefits, the scan requires flood extent/flood depth data 
# as flood hazard data (provided from Fathom). Defended pluvial flood data is used for 
# assessing flood benefits, but data may also be pre-processed, e.g., to combine fluvial 
# and pluvial flood data, as required for a case. Retrieved/pre-processed data will be 
# ingested from a scan from a defined location. Moreover, the return periods for the data 
# needs to be specified, using a Dict[str,int] (example below shows the naming convention 
# for Fathom, but this can obviously vary):

flood_data_path = "path to data"
flood_mapping = {
	"P_1in5.tif": 5,
	"P_1in10.tif": 10,
	"P_1in20.tif": 20,
	"P_1in50.tif": 50,
	"P_1in100.tif": 100,
	"P_1in200.tif": 200,
	"P_1in500.tif": 500,
	"P_1in1000.tif": 1000,
}

# define path to hazard extent used for floodplain delineation - this is needed for 
# the floodplain nbs family.
floodplain_reference_path = "path to tif used for floodplain delineation"

# define NbS of interest
# comment NBS families not relevant for a given case
nbs_list = [
	NbsTypes.BioretentionAreas,
	NbsTypes.UrbanForest,
	NbsTypes.GreenCorridors,
	NbsTypes.InlandWetlands,
	NbsTypes.OpenGreenSpace,
	NbsTypes.RiverFloodplain,
	NbsTypes.RiverRenaturation,
	NbsTypes.TerracesAndSlopes,
]

# define benefits
# Note that for the PluvialFlooding benefit, return_periods always need to be specified.
# for parameters to be set in params dicts, refer to the following methods in nbsos.py:
# * get_biodiversity_creation_benefits
# * get_heat_creation_benefits
# * get_recreation_creation_benefits (for recreation and social cohesion)
# * get_pluvial_flood_creation_benefits
# * get_job_creation_creation_benefits
# * get_erosion_creation_benefits

# define additional parameters as needed
benefits = [
	ScanBenefit(Benefits.UrbanHeat),
	ScanBenefit(Benefits.Recreation),
	ScanBenefit(Benefits.SocialCohesion),
	ScanBenefit(Benefits.PluvialFlooding, params={"return_periods": list(flood_mapping.values())}),
	ScanBenefit(Benefits.Biodiversity),
	ScanBenefit(Benefits.AvoidedErosion),
]

# instantiate scan
my_scan = Scan(
	scan_context=scan_context,
	nbs_list=nbs_list,
)

# instantiate storage handler; here, using the GcsHandler
my_storage_handler = GcsHandler(
	scan_context=scan_context,
	flood_hazard_data_folder=flood_data_path,
	file_to_probability_mapping=flood_mapping,
)
my_storage_handler.credentials = "path to credentials json"
my_storage_handler.bucket_name = "urban-nbs-scan-bucket"

# pass storage handler to scan. This will also intialize the handler
my_scan.set_storage_options(storage_handler=my_storage_handler)

# run analytics
# prepare indicators and download to local device
my_scan.run_preparation(
	run_gee=True,
	run_download=True,
	run_client_computations=True,
)

# run suitability model
my_scan.run_suitability(
	floodplain_path=floodplain_reference_path,
)

# run clustering for visual inspection
my_scan.explore_clusters()

# the cluster threshold derived from visual inspection to be inserted here, for example 8.
cluster_threshold = 8

# define masks for creation opportunities, as needed, 
# as a list of str referring to shapefiles in the masks folder ['shapefile1_name.shp', '...']
# set to None if no masks are needed (None is the default value)
masks = None

my_scan.run_clustering(
	clustering_threshold=cluster_threshold,
	exclusion_masks=masks,
)

# run creation benefit models
# see method in the nbsos.py
my_scan.run_creation_benefits(
	benefit_parameters=benefits,
)

# run protection benefit models
# by default, use lppi protection mask type
my_scan.run_protection_benefits(
	benefit_parameters=benefits,
)

# create priority areas
my_scan.run_priority(
	benefit_parameters=benefits,
)

# run multicriteria analysis, including comparion of relative performances of nbs per benefit, etc.
# adapt focal weight for benefit heavy-scenarios, as needed
focal_weight = 0.8

# specify whether to run scenarios (equal, heavy, only)
assess_equal_weights_scenario = True
assess_benefit_heavy_scenarios = True
assess_benefit_only_scenarios = True

# define custom scenarios, as needed.
# This is a dict of dicts (key refers to scenario name, and value to a dict of benefit-weight pairs).
# If not needed, set to None.
custom_scenarios = {
	"balanced": {
		Benefits.PluvialFlooding: 0.3,
		Benefits.AvoidedErosion: 0.3,
		Benefits.UrbanHeat: 0.1,
		Benefits.Biodiversity: 0.1,
		Benefits.Recreation: 0.1,
		Benefits.JobsCreation: 0.1,
	}
}

# in addition, specify percentile thresholds to define the highest X% benefits to be exported
# E.g., [80, 90, 95] for highest 20%, 10%, and 5% benefits, respectively.
# By default, 80th percentile, i.e., highest 20% benefits.
highest_benefits_percentiles = [80]

# run multicriteria analysis
my_scan.run_multicriteria_analysis(
	benefit_parameters=benefits,
	heavy_focal_weight=focal_weight,
	assess_equal_weights_scenario=assess_equal_weights_scenario,
	assess_benefit_heavy_scenario=assess_benefit_heavy_scenarios,
	assess_benefit_only_scenario=assess_benefit_only_scenarios,
	custom_scenarios=custom_scenarios,  # or set to None
	highest_benefits_percentiles=highest_benefits_percentiles,
)


# run building solutions, after clustering and priority area identification, if needed
my_scan.run_building_solutions(
	benefit_parameters=benefits,
	exclude_non_built_areas=False,
)

# close scan, save logs
my_scan.archive_data()
my_scan.close()

