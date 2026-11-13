from .soilwatch import s2_masks, S2Fcover, RUSLE_factors
from . import loggerconfig

import geemap
import os
import logging
import ee
import datetime
import math
import xarray as xr
import rioxarray
from tqdm.auto import tqdm as tqauto
import time

# logger = logging.getLogger(__name__)
# logger = loggerconfig.config_logger(logger)
logger = logging.getLogger('nbsos.scan')


class Indicators:
    CHIRPS_PENTAD = "UCSB-CHG/CHIRPS/PENTAD"
    CHIRPS_DAILY = "UCSB-CHG/CHIRPS/DAILY"
    DEM = "projects/sat-io/open-datasets/FABDEM"
    HAND = "MERIT/Hydro/v1_0_1"

    ISDASOIL_TEXTURECLASS = "ISDASOIL/Africa/v1/texture_class"
    ISDASOIL_SANDCONTENT = "ISDASOIL/Africa/v1/sand_content"
    ISDASOIL_PH = "ISDASOIL/Africa/v1/ph"
    ISDASOIL_SILTCONTENT = "ISDASOIL/Africa/v1/silt_content"
    ISDASOIL_CLAYCONTENT = "ISDASOIL/Africa/v1/clay_content"
    ISDASOIL_BULKDENSITY = "ISDASOIL/Africa/v1/bulk_density"
    ISDASOIL_CARBONORGANIC = "ISDASOIL/Africa/v1/carbon_organic"
    ISDASOIL_BEDROCKDEPTH = "ISDASOIL/Africa/v1/bedrock_depth"
    ISDASOIL_STONECONTENT = "ISDASOIL/Africa/v1/stone_content"

    OPENLANDMAP_TEXTURE = "OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02"
    SOILGRIDS_SANDCONTENT = "projects/soilgrids-isric/sand_mean"
    SOILGRIDS_PH = 'projects/soilgrids-isric/phh2o_mean'
    SOILGRIDS_SILTCONTENT = "projects/soilgrids-isric/silt_mean"
    SOILGRIDS_CLAY = "projects/soilgrids-isric/clay_mean"
    SOILGRIDS_BULKDENSITY = "projects/soilgrids-isric/bdod_mean"
    SOILGRIDS_CARBON = "projects/soilgrids-isric/soc_mean"
    SOILGRIDS_STONECONTENT = 'projects/soilgrids-isric/cfvo_mean'
    SOILGRIDS_DEPTH = "users/wb-nbs/soil/soilgrids_dtb"

    S2_COLLECTION = 'COPERNICUS/S2_SR_HARMONIZED'
    JRC_HIST_WATER = "JRC/GSW1_4/MonthlyHistory"
    ESA_WORLDCOVER = "ESA/WorldCover/v100"
    GLAD_LAND = "projects/glad/landBuffer4"
    GLAD_LULC = 'projects/glad/GLCLU2020/LCLUC_2020'
    GOOGLE_BUILDING = "GOOGLE/Research/open-buildings/v3/polygons"
    WSF = "projects/sat-io/open-datasets/WSF/WSF_2019"
    WSF_EVO = "projects/sat-io/open-datasets/WSF/WSF_EVO"
    WRI_TREES = 'projects/wri-datalab/TML'
    HRSL_POP = "projects/sat-io/open-datasets/hrsl/hrslpop"
    JRC_SURFACE_WATER = "JRC/GSW1_2/GlobalSurfaceWater"
    MODIS_LST = "MODIS/061/MOD11A1"


    def __init__(self, country, city, root_folder, bucket, africa=True, project=None):
        self.root = os.path.join(root_folder, country, city)
        self.bucket = bucket
        self.country = country
        self.city = city
        self.aoi = self.get_aoi()
        # to clip all data layers with
        self.clip_aoi = self.aoi.geometry().buffer(2000)
        self.africa = africa
        self.project = project

        # times
        today = datetime.date.today()
        # get the last full year
        self.end_year = today.year - 1

    def get_aoi(self):
        """
        Import the AOI shapefile into GEE. Shapefile path needs to be exactly:
        "/{root_folder}/{country}/{city}/[ancillary}/{city}_aoi.shp"
        :return: aoi GEE featurecollection object with AOI
        """
        aoi_path = os.path.join(self.root, 'ancillary', f'{self.city}_aoi.shp')
        logger.info(f"Importing {self.city} AOI into Earth Engine from {aoi_path}")
        try:
            aoi_tmp = geemap.shp_to_ee(aoi_path)

            geoms = aoi_tmp.geometry().geometries()
            aoi = ee.FeatureCollection(
                geoms.map(lambda g: ee.Feature(ee.Geometry(g)))
            )

        except Exception as e:
            logger.error(f"import {aoi_path} failed: {e}", exc_info=True)
            raise e

        return aoi

    def get_factor_K_Benavidez(self):

        try:
            # Salhi, A., El Hasnaoui, Y., Cutillas, P. P., & Heggy, E., 2023, Scientific Reports, 13, 13158.
            buffered_aoi = self.aoi.geometry().buffer(5000)

            soilTextureImage = ee.Image("OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02").clip(buffered_aoi)

            soil_raw = soilTextureImage.select('b0').rename('soil')
            K = soil_raw.expression( 
                "(b('soil') == 12) ? 0.0053" + #  Clay
                ": (b('soil') == 11) ? 0.0170" + # Silty clay
                ": (b('soil') == 10) ? 0.045" +  #  Sandy clay
                ": (b('soil') == 9) ? 0.050" +  #  Clay loam
                ": (b('soil') == 8) ? 0.0499" + #  Silty clay loam
                ": (b('soil') == 7) ? 0.0394" + #  Sandy clay loam
                ": (b('soil') == 6) ? 0.0264" + #  Loam
                ": (b('soil') == 5) ? 0.0423" + #  Silt loam
                ": (b('soil') == 4) ? 0.0394" + #  Silt
                ": (b('soil') == 3) ? 0.036" +  #  Sandy loam
                ": (b('soil') == 2) ? 0.0341" + #  Loamy sand
                ": (b('soil') == 1) ? 0.0288" + #  Sand
                ": 0").rename('K')
            
            return K.clip(buffered_aoi)
        
        except Exception as e:
            print(f"Failed to compute factor K: {e}")
            raise e

    def get_annual_precipitation(self, smoothing_years=3):
        """
        Get the annual precipitation rate for AOI using CHIRSP data
        Current repository is
        https://developers.google.com/earth-engine/datasets/catalog/UCSB-CHG_CHIRPS_DAILY
        Annual precipitation is mean annual precipitation of x number of years
        :param smoothing_years: set the nu,ber of years to use for mean annual precipitation
                default: 3
        :return: Annual precipitation GEE image object
        """
        # location of the data in GEE

        start_year = self.end_year - smoothing_years + 1

        # get list of years
        years = [*range(start_year, self.end_year + 1, 1)]

        # calculate annual precipitaiont for each year
        logger.info(f'Importing CHIRPS precipitation data from {self.CHIRPS_PENTAD}')
        try:
            precipitation = [ee.ImageCollection(self.CHIRPS_PENTAD).filterDate(
                f'{year}-01-01', f'{year}-12-31').filterBounds(
                self.clip_aoi).select('precipitation').sum() for year in years]

            # get average annual precipitation
            annual_prec = ee.ImageCollection.fromImages(precipitation).mean()
            clipped = annual_prec.clip(self.clip_aoi)
            return clipped
        except Exception as e:
            logger.error(f'Importing precipitation data from {self.CHIRPS_PENTAD} failed: {e}', exc_info=True)

    def get_heavy_prec(self, smoothing_years=3):
        """
        Get the annual number of heavy precipitation days (=> 10mm) for AOI using CHIRPS data
        Current repository is
        https://developers.google.com/earth-engine/datasets/catalog/UCSB-CHG_CHIRPS_DAILY
        Annual  is mean annual heavy precipitation days of x number of years
        :param smoothing_years: set the nu,ber of years to use for mean annual precipitation
                default: 3
        :return: Mean annual days of heavy precipitation GEE image object
        """


        start_year = self.end_year - smoothing_years + 1
        start_date = f'{start_year}-01-01'
        end_date = f'{self.end_year}-12-31'

        logger.info(f'Importing CHIRPS precipitation data for heavy rainfall days from {self.CHIRPS_DAILY}')
        try:
            prec_daily = ee.ImageCollection(self.CHIRPS_DAILY).filterDate(start_date, end_date).select(
                'precipitation').filterBounds( # filter for aoi
                self.clip_aoi).map( # count days over 10mm precipitation
                    lambda image: image.gt(10).selfMask()).reduce( # count the number of times > 10
                        ee.Reducer.count(), 4).divide(3) # divide by 3 for 3 years

            clipped = prec_daily.clip(self.clip_aoi)
            return clipped

        except Exception as e:
            logger.error(f'Importing precipitation data from {self.CHIRPS_DAILY} failed: {e}', exc_info=True)


    def get_lst(self):
        """
            This function retrieves Land Surface Temperature (LST) data for both day and night from
            MODIS satellite imagery by Earth Engine (EE) for a specified area of interest (AOI) and
            year range. The function computes monthly mean LST values, converts them to Kelvin,
            assigns descriptive names, and outputs them in a structured format.

            Returns:
                A list of dictionaries, each containing the monthly processed LST image, its
                descriptive name, and the associated scale.

            Parameters:
                self (object): Instance of the class containing configurations for MODIS LST
                               collection, AOI, relevant end year, and city name.

            Raises:
                None
        """
        end = f'{self.end_year}-12-31'
        start = f'{self.end_year - 4}-01-01'
        
        #add a 5km buffer aorund lst to be able to calculate later
        buf_geom = self.aoi.geometry().buffer(5000)
        
        modisLST_day = ee.ImageCollection(
            self.MODIS_LST).filterDate(
            start, end).select(
            'LST_Day_1km').filterBounds(
            buf_geom)

        modisLST_night = ee.ImageCollection(
            self.MODIS_LST).filterDate(
            start, end).select(
            'LST_Night_1km').filterBounds(
            buf_geom)


        months = range(1, 13)
        def calc_monthly_lst_day(ic, month, timeOfDay):
            monthLST = ic.filter(ee.Filter.calendarRange(month, month, 'month')).mean()
            monthLST = monthLST.clip(buf_geom)
            monthLST = monthLST.multiply(0.02)
            name_month = datetime.date(1900, month, 1).strftime('%B')
            month_nr = ee.Number(month).format('%02d').getInfo()
            name = f'LST_{self.city}_{timeOfDay}_{name_month}_{month_nr}'
            lst_dict = {'image': monthLST.toFloat(), 'name': name, 'scale': 1000,
                        'geometry': buf_geom}
            return lst_dict

        days = [calc_monthly_lst_day(modisLST_day, month, 'day') for month in months]
        nights = [calc_monthly_lst_day(modisLST_night, month, 'night') for month in months]
        all_lst = days + nights
        return all_lst

    def get_slope_rad(self):
        """
        Calculates the slope of a terrain in radians using a digital elevation model (DEM). Retrieves the DEM data and computes the
        slope using Earth Engine's terrain functions.

        Returns:
            ee.image.Image
            An Earth Engine Image object representing the slope data in radians.

        Raises:
            Exception
            If the slope calculation fails, the method logs the error and raises an exception.
        """
        logger.info(f'Calculating slope (radians) data...')
        try:
            dem = self.get_dem()
            slope_rad = ee.Terrain.slope(dem).multiply(math.pi / 180).clip(self.clip_aoi)

            return slope_rad

        except Exception as e:
            logger.error(f'Calculating slope (radians) data failed: {e}', exc_info=True)

    def get_aspect(self):
        """
        Calculates the aspect of a terrain in degrees using a digital elevation model (DEM). Retrieves the DEM data and computes the
        aspect using Earth Engine's terrain functions.

        Returns:
            ee.image.Image
            An Earth Engine Image object representing the aspect data in degrees.

        Raises:
            Exception
            If the aspect calculation fails, the method logs the error and raises an exception.
        """
        logger.info(f'Calculating aspect data...')
        try:
            dem = self.get_dem()
            aspect_deg = ee.Terrain.aspect(dem).clip(self.clip_aoi)

            return aspect_deg

        except Exception as e:
            logger.error(f'Calculating aspect data failed: {e}', exc_info=True)

    def get_dem(self):
        """
        Retrieves and processes Digital Elevation Model (DEM) data.

        Summary:
        This function attempts to import Digital Elevation Model (DEM) data from a specified
        source (FABDEM in this case) using the Earth Engine API. The function filters the DEM
        data to the area of interest (AOI), determines the projection of the DEM, mosaics the
        DEM, and sets the default projection for further use.

        Raises an error with a logged message if the import fails for any reason.

        Returns:
            ee.Image: A processed Earth Engine DEM image object.

        Raises:
            Exception: If any error occurs during the DEM data import or processing.
        """

        logger.info(f'Importing FABDEM DEM data from {self.DEM}')
        try:
            dem = ee.ImageCollection(self.DEM)
            dem_proj = dem.first().projection()
            dem = dem.mosaic().setDefaultProjection(dem_proj).clip(self.clip_aoi)

            return dem
        except Exception as e:
            logger.error(f'Importing DEM data from {self.DEM} failed: {e}', exc_info=True)

    def get_slope(self):
        """
        Calculates the slope of a terrain in degrees using a digital elevation model (DEM). Retrieves the DEM data and computes the
        slope using Earth Engine's terrain functions.

        Returns:
            ee.image.Image
            An Earth Engine Image object representing the slope data in degrees.

        Raises:
            Exception
            If the slope calculation fails, the method logs the error and raises an exception.
        """
        logger.info(f'Calculating slope data...')
        try:
            dem = self.get_dem()
            slope_deg = ee.Terrain.slope(dem).clip(self.clip_aoi)

            return slope_deg

        except Exception as e:
            logger.error(f'Calculating slope data failed: {e}', exc_info=True)

    def get_hand(self):
        """
        Retrieve the HAND (Height Above Nearest Drainage) data for the specified area
        of interest (AOI). This method imports HAND data using the Earth Engine API,
        clips it to the AOI's geometry, and returns the resulting image object. If there
        is an issue during the import, a warning is logged.
        Hydrologically adjusted elevations, also known as "hand" (height above the nearest drainage).
        The elevations are adjusted to satisfy the condition "downstream is not higher than its upstream" while
        minimizing the required modifications from the original DEM. The elevation above EGM96 geoid is
        represented in meters, and the vertical increment is set to 10cm. For detailed method,
        see [Yamazaki et al., 2012, WRR].

        Returns:
            ee.Image: The clipped HAND data image.

        Raises:
            Exception: If the HAND data import or processing fails.

        """
        logger.info(f'Importing HAND data from {self.HAND}')
        try:
            hand = ee.Image(self.HAND).select('hnd').clip(self.clip_aoi)
            return hand
        except Exception as e:
            logger.warning(f'Importing HAND data failed: {e}', exc_info=True)


    def get_soil_properties_africa(self):
        """
                Get the soil properties fro Africa using ISDASOIL. This method should only be used when
                AOI is in Africa.
                :return: GEE imageCollection of soil properties
        """
        logger.info("Assuming case study is in AFRICA")
        logger.info(f'Importing soil properites data from {self.ISDASOIL_TEXTURECLASS}')
        try:
            soil_properties = ee.Image(self.ISDASOIL_TEXTURECLASS).select(0).rename('texture').toFloat().addBands(
                ee.Image(self.ISDASOIL_SANDCONTENT).select(0).rename('sand')).toFloat().addBands(
                ee.Image(self.ISDASOIL_PH).select(0).divide(10).rename('ph')).toFloat().addBands(
                ee.Image(self.ISDASOIL_SILTCONTENT).select(0).rename('silt')).toFloat().addBands(
                ee.Image(self.ISDASOIL_CLAYCONTENT).select(0).rename('clay')).toFloat().addBands(
                ee.Image(self.ISDASOIL_BULKDENSITY).select(0).divide(100).rename('bulk_density')).toFloat().addBands(
                ee.Image(self.ISDASOIL_CARBONORGANIC).select(0).divide(10).exp().subtract(1).rename('soil_carbon')).toFloat().addBands(
                ee.Image(self.ISDASOIL_BEDROCKDEPTH).select(0).rename('depth')).toFloat().addBands(
                ee.Image(self.ISDASOIL_STONECONTENT).select(0).divide(10).exp().subtract(1).toFloat().rename('stone_content')).clip(
                self.clip_aoi)

            return soil_properties

        except Exception as e:
            logger.error(f'Importing Africa soil properties data failed: {e}', exc_info=True)


    def get_soil_permeability_africa(self):
        """
        Fetches and processes soil permeability data for the specified area of interest (AOI) in Africa
        based on predefined texture classes. The function retrieves the soil texture image, clips it to
        the AOI, and computes the permeability class using a combination of logical and arithmetic operations
        on the texture bands.

        Returns:
            ee.Image
            An Earth Engine image where each pixel represents a permeability class derived from
            the texture data. Permeability is categorized using specific texture class mappings.

        Raises:
        E   xception
            If an error occurs during the retrieval or processing of soil permeability data,
            it logs an error message and the exception details.
        """
        logger.info(f'Importing soil permeability data from {self.ISDASOIL_TEXTURECLASS}')
        try:
            texture = ee.Image(self.ISDASOIL_TEXTURECLASS).select(0).rename('texture').clip(self.clip_aoi)

            permeability_class = texture.eq(1).Or(texture.eq(2)).multiply(6).add(
                texture.eq(3).Or(texture.eq(5)).multiply(5)).add(
                texture.eq(4).Or(texture.eq(6)).multiply(4)).add(
                texture.eq(7).Or(texture.eq(8)).Or(texture.eq(10)).multiply(3)).add(
                texture.eq(9).Or(texture.eq(11)).multiply(2)).add(
                texture.eq(12))

            return permeability_class

        except Exception as e:
            logger.error(f"Calculating soil permeability failed: {e}", exc_info=True)

    def get_hydrosoil_group_africa(self):
        """
               Get the hydrosoil group (based on texture) Africa using ISDASOIL. This method
               should only be used when  AOI is in Africa.
               :return: GEE imageCollection of hydrosoil classes
        """
        logger.info(f'Importing and calculating Hydrosoil group with data from {self.ISDASOIL_TEXTURECLASS}')
        try:
            texture = ee.Image(self.ISDASOIL_TEXTURECLASS).select(0).rename('texture').clip(self.clip_aoi)

            hydrosoil_group = texture.eq(12).Or(texture.eq(11)).Or(texture.eq(9)).multiply(1).add(
                texture.eq(7).Or(texture.eq(8)).Or(texture.eq(10)).multiply(2)).add(
                texture.eq(6).multiply(3)).add(
                texture.eq(4).Or(texture.eq(5)).Or(texture.eq(3)).Or(texture.eq(1)).Or(texture.eq(2)).multiply(4))

            return hydrosoil_group

        except Exception as e:
            logger.warning(f"Calculating hydrosoil group failed: {e}", exc_info=True)

    def get_soil_properties(self):
        """
                Get the soil properties for non Africa using OpenLandMap. This method should only be used when
                AOI is not in Africa.
                :return: GEE imageCollection of soil properties
        """
        logger.info("Assuming case study is NOT in AFRICA")
        logger.info(f'Importing soil properties data from {self.OPENLANDMAP_TEXTURE} and {self.SOILGRIDS_SANDCONTENT}')
        try:
            soil_properties = ee.Image(self.OPENLANDMAP_TEXTURE).select('b0').toFloat().addBands(
                ee.Image(self.OPENLANDMAP_TEXTURE).select('b10')).toFloat().addBands(
                ee.Image(self.OPENLANDMAP_TEXTURE).select('b30')).reduce(
                ee.Reducer.mode()).rename('texture').toFloat().addBands(
                ee.Image(self.SOILGRIDS_SANDCONTENT).select('sand_0-5cm_mean').add(
                    ee.Image(self.SOILGRIDS_SANDCONTENT).select('sand_5-15cm_mean')).add(
                    ee.Image(self.SOILGRIDS_SANDCONTENT).select('sand_15-30cm_mean')
                ).divide(3).divide(10).rename('sand')).toFloat().addBands(
                ee.Image(self.SOILGRIDS_PH).select('phh2o_0-5cm_mean').add(
                    ee.Image(self.SOILGRIDS_PH).select('phh2o_5-15cm_mean')).add(
                    ee.Image(self.SOILGRIDS_PH).select('phh2o_15-30cm_mean')
                ).divide(3).divide(10).rename('ph')).toFloat().addBands(
                ee.Image(self.SOILGRIDS_SILTCONTENT).select('silt_0-5cm_mean').add(
                    ee.Image(self.SOILGRIDS_SILTCONTENT).select('silt_5-15cm_mean')).add(
                    ee.Image(self.SOILGRIDS_SILTCONTENT).select('silt_15-30cm_mean')
                ).divide(3).divide(10).rename('silt')).toFloat().addBands(
                ee.Image(self.SOILGRIDS_CLAY).select('clay_0-5cm_mean').add(
                    ee.Image(self.SOILGRIDS_CLAY).select('clay_5-15cm_mean')).add(
                    ee.Image(self.SOILGRIDS_CLAY).select('clay_15-30cm_mean')
                ).divide(3).divide(10).rename('clay')).toFloat().addBands(
                ee.Image(self.SOILGRIDS_BULKDENSITY).select('bdod_0-5cm_mean').add(
                    ee.Image(self.SOILGRIDS_BULKDENSITY).select('bdod_5-15cm_mean')).add(
                    ee.Image(self.SOILGRIDS_BULKDENSITY).select('bdod_15-30cm_mean')
                ).divide(3).divide(100).rename('bulk_density')).addBands(
                ee.Image(self.SOILGRIDS_CARBON).toFloat().select('soc_0-5cm_mean').add(
                    ee.Image(self.SOILGRIDS_CARBON).select('soc_5-15cm_mean')).add(
                    ee.Image(self.SOILGRIDS_CARBON).select('soc_15-30cm_mean')
                ).divide(3).divide(10).rename('soil_carbon')).toFloat().addBands(
                ee.Image(self.SOILGRIDS_STONECONTENT).select('cfvo_0-5cm_mean').add(
                    ee.Image(self.SOILGRIDS_STONECONTENT).select('cfvo_5-15cm_mean')).add(
                    ee.Image(self.SOILGRIDS_STONECONTENT).select('cfvo_15-30cm_mean')
                ).divide(3).divide(100).rename('stone_content')).addBands(
                ee.Image(self.SOILGRIDS_DEPTH).toFloat().select(0).toFloat().rename('depth')).unmask(0).clip(
                self.clip_aoi)

            return soil_properties

        except Exception as e:
            logger.error(
                f'Importing soil properties data from \
                {self.OPENLANDMAP_TEXTURE} and \
                {self.SOILGRIDS_SANDCONTENT} failed: {e}',
                exc_info=True)

    def get_soil_permeability(self):
        """
               Get the soil permeability (based on texture) Africa using OpenLandMap. This method
               should only be used when  AOI is NOT in Africa.
               :return: GEE imageCollection of soil permeability
        """

        logger.info(f'Importing and calculating soil permeability data from {self.OPENLANDMAP_TEXTURE}')
        try:

            texture = ee.Image(self.OPENLANDMAP_TEXTURE).select('b0').addBands(
                ee.Image(self.OPENLANDMAP_TEXTURE).select('b10')).addBands(
                ee.Image(self.OPENLANDMAP_TEXTURE).select('b30')).reduce(
                ee.Reducer.mode()).rename('texture')

            permeability_class = texture.eq(1).Or(texture.eq(2)).multiply(6).add(
                texture.eq(3).Or(texture.eq(5)).multiply(5)).add(
                texture.eq(4).Or(texture.eq(6)).multiply(4)).add(
                texture.eq(7).Or(texture.eq(8)).Or(texture.eq(10)).multiply(3)).add(
                texture.eq(9).Or(texture.eq(11)).multiply(2)).add(
                texture.eq(12).multiply(1)).clip(self.clip_aoi)

            return permeability_class

        except Exception as e:
            logger.warning(
                f'Importing and calculating soil permeability data from {self.OPENLANDMAP_TEXTURE} failed {e}',
                exc_info=True)

    def get_hydrosoil_group(self):
        """
               Get the hydrosoil group (based on texture) Africa using OpenLandMap. This method
               should only be used when  AOI is NOT in Africa.
               :return: GEE Image of shydrosoil group
        """
        logger.info(f'Importing and calculating hydrosoil group data from {self.OPENLANDMAP_TEXTURE}')
        try:

            texture = ee.Image(self.OPENLANDMAP_TEXTURE).select('b0').addBands(
                ee.Image(self.OPENLANDMAP_TEXTURE).select('b10')).addBands(
                ee.Image(self.OPENLANDMAP_TEXTURE).select('b30')).reduce(
                ee.Reducer.mode()).rename('texture')
            hydrosoil_group = ee.Image(1).where(
                texture.eq(12).Or(texture.eq(11)).Or(texture.eq(9)), 1).where(
                texture.eq(7).Or(texture.eq(8)).Or(texture.eq(10)), 2).where(
                texture.eq(6), 3).where(
                texture.eq(4).Or(texture.eq(5)).Or(texture.eq(3)).Or(texture.eq(1)).Or(texture.eq(2)),
                4).clip(self.clip_aoi)


            return hydrosoil_group

        except Exception as e:
            logger.error(f'Importing and calculating hydrosoil group data from {self.OPENLANDMAP_TEXTURE} failed: {e}',
                           exc_info=True)

    def get_bs_freq(self, smoothing_years=3):
        """
        Get Baresoil frequency based on cloud free Senitnel 2 data. Bare soil frequency is calculated as
        the number of bare soil observations divided by the number of cloud-free observations
        Note that some helper functions are stored in the BareSoil class at the bottom of this module. This method
        makes use of custom code written in JS by William OUelette and translated to Python by Koen Tieskens
        :param smoothing_years: number of years to use for calculating the percentage
        :return: GEE IMage with bare soil frequency as fraction (0-1).
        """
        logger.info(f'Importing and calculating bare soil frequency group data from {self.S2_COLLECTION}')
        try:

            COLLECTION = self.S2_COLLECTION
            start_year = self.end_year - smoothing_years + 1
            start_date = f'{start_year}-01-01'
            end_date = f'{self.end_year}-12-31'
            s2_date_range = ee.Dictionary({'start': start_date, 'end': end_date})
            s2_cl = s2_masks.loadImageCollection(COLLECTION, s2_date_range, self.clip_aoi)
            masked_collection = cloud_mask(s2_cl)

            def geos_3func(img):
                """
                little helper function
                :param img:
                :return:
                """
                geos3 = s2_masks.addGEOS3Mask(img)
                out = img.updateMask(geos3)
                return out
            bs_collection = masked_collection.map(geos_3func)
            bs_freq = bs_collection.select('B2').reduce(ee.Reducer.count(), 4).divide(
                masked_collection.select('B2').reduce(ee.Reducer.count(), 4)).unmask(0).rename(
                'bare_soil_frequency').toFloat().clip(self.clip_aoi)

            return bs_freq

        except Exception as e:
            logger.warning(
                f'Importing and calculating bare soil frequency group data from {self.S2_COLLECTION} failed: {e}',
                exc_info=True)


    def get_water_freq(self):
        """
        Calculates water frequency using historical monthly water observations and smooths the
        result over the past three years.

        The function imports historical monthly water observations from JRC water frequency data,
        performs calculations to determine water frequency over three years, and applies a
        smoothing function. The smoothed water frequency is clipped to the area of interest
        (clip_aoi) and returned as an Earth Engine image.

        Raises:
            Exception: If importing or calculating water frequency fails.

        Returns:
            ee.Image: Clipped smoothed water frequency image for the area of interest.
        """
        logger.info(
            f"Importing and calculating water frequency using monthly water observations from {self.JRC_HIST_WATER}")
        try:
            jrc_water_freq = ee.ImageCollection(self.JRC_HIST_WATER).filterDate(str(self.end_year-3),
                                                                                str(self.end_year))
            def smoother(img):
                yearCol = ee.ImageCollection.fromImages(img.get('year_matches'))
                out = yearCol.map(lambda x:  x.eq(2).rename('water').selfMask()).select('water').reduce(
                    ee.Reducer.count(),4).divide(
                        yearCol.map(lambda y: y.gt(0).selfMask()).select('water').reduce(
                            ee.Reducer.count(),4)).unmask(0).rename('water_frequency').set(
                    'system:time_start',
                    ee.Date(img.date().update(
                        year=ee.Number(img.date().get('year')),
                        month=7,
                        day=1,
                        hour=0,
                        minute=0,
                        second=0)).millis())
                return out

            water_ts_smooth = harmonize_ts(jrc_water_freq).map(smoother).map(
                    lambda z: z.addBands(z.metadata('system:time_start').divide(3.154e10).add(1970))
                )

            water_freq = water_ts_smooth.sort('system:time_start', False).first().select('water_frequency')
            clipped = water_freq.clip(self.clip_aoi)
            return clipped
        except Exception as e:
            logger.error(
                f"Importing and calculating water frequency using monthly water observations from \
                {self.JRC_HIST_WATER} failed: {e}", exc_info=True)

    def get_land_cover(self):
        """
        Retrieves and processes land cover data using ESA worldcover V1

        This method fetches land cover data from ESA Worldcover V1
        dataset, processes it, and  updates it with additional land cover
        change information from GLAD landcover change data before
        clipping it to the specified area of interest (AOI).

        Raises:
            Exception: If importing or processing land cover data fails.

        Returns:
            ee.Image: The processed and clipped land cover map.
        """
        logger.info(f"Importing land cover data from {self.ESA_WORLDCOVER}")
        try:
            land_cover = ee.ImageCollection(self.ESA_WORLDCOVER).filterBounds(
                self.clip_aoi).mosaic().remap(  # remap landcover
                [10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100],
                [3, 4, 5, 2, 1, 8, 9, 6, 7, 3, 5]).rename('constant')

            # update with glad landcover change
            landmask = ee.Image(self.GLAD_LAND).mask()
            m20 = ee.Image(self.GLAD_LULC).updateMask(landmask)
            land_cover = land_cover.where(m20.gt(100).And(m20.lt(150)), 7).toInt32()
            clipped = land_cover.clip(self.clip_aoi)
            return clipped

        except Exception as e:
            logger.warning(f"Importing land cover data from {self.ESA_WORLDCOVER} failed: {e}", exc_info=True)

    def get_lcc(self):
        """
        Land capability class, see https://www.jswconline.org/content/early/2020/06/09/jswc.2020.00023
        :return: land capability classification gee image
        """
        if self.africa:
            soil_properties = self.get_soil_properties_africa()
            permeability_class = self.get_soil_permeability_africa()
        else:
            soil_properties = self.get_soil_properties()
            permeability_class = self.get_soil_permeability()
        logger.info("Calculating land capability class")
        try:
            coarse_fragments = soil_properties.select('stone_content')
            slope_deg = self.get_slope()
            slope_perc = slope_deg.multiply(ee.Number(math.pi).divide(180)).tan().multiply(100)
            om = soil_properties.select('soil_carbon').divide(1000).multiply(1.72)
            K = RUSLE_factors.factorK(self.africa)

            wp1 = soil_properties.select('sand').divide(100).multiply(-0.024).add(
                soil_properties.select('clay').divide(100).multiply(0.487)).add(
                om.multiply(0.006)).add(
                om.multiply(soil_properties.select('sand').divide(100)).multiply(0.005)).subtract(
                soil_properties.select('clay').divide(100).multiply(om).multiply(0.013)).add(
                soil_properties.select('sand').divide(100).multiply(soil_properties.select('clay').divide(100)).multiply(
                0.068)).add(0.031)

            wp = wp1.add(wp1.multiply(0.14).subtract(0.02))

            fc1 = soil_properties.select('sand').divide(100).multiply(-0.251).add(
                soil_properties.select('clay').divide(100).multiply(0.195)).add(
                om.multiply(0.011)).add(
                om.multiply(soil_properties.select('sand').divide(100)).multiply(0.006)).subtract(
                soil_properties.select('clay').divide(100).multiply(om).multiply(0.027)).add(
                soil_properties.select('sand').divide(100).multiply(soil_properties.select('clay').divide(100)).multiply(
                0.452)).add(0.299)

            fc = fc1.add(fc1.pow(2).multiply(1.283).subtract(fc1.multiply(0.374)).subtract(0.015))

            awc = fc.subtract(wp).multiply(ee.Image(1).subtract(coarse_fragments.divide(100))).multiply(100)

            lcc = slope_perc.lte(2).And(K.multiply(100).divide(0.1317).gt(32)).Or(
                    slope_perc.lte(5).And(K.multiply(100).divide(0.1317).lte(32))).Or(
                    soil_properties.select('depth').gte(100)).Or(
                    soil_properties.select('texture').gte(4).And(soil_properties.select('texture').lte(10))).Or(
                    coarse_fragments.lt(0.1)).And(awc.gt(18)).And(permeability_class.gte(3)).rename(
                    'lcc1').addBands(
                slope_perc.gt(2).And(slope_perc.lte(5)).And(K.multiply(100).divide(0.1317).gt(32)).Or(
                    slope_perc.gt(5).And(slope_perc.lte(10)).And(K.multiply(100).divide(0.1317).lte(32))).Or(
                    soil_properties.select('texture').eq(2).Or(soil_properties.select('texture').eq(3)).Or(
                    soil_properties.select('texture').eq(11)).Or(soil_properties.select('texture').eq(12))).Or(
                    coarse_fragments.gte(0.1).And(coarse_fragments.lt(3))).Or(
                    awc.lte(18).And(awc.gt(12))).Or(
                    permeability_class.eq(2)).multiply(2).rename(
                    'lcc2')).addBands(
                slope_perc.gt(5).And(slope_perc.lte(10)).And(K.multiply(100).divide(0.1317).gt(32)).Or(
                    slope_perc.gt(10).And(slope_perc.lte(15)).And(K.multiply(100).divide(0.1317).lte(32))).Or(
                    soil_properties.select('depth').gte(70).And(soil_properties.select('depth').lt(100))).Or(
                    soil_properties.select('texture').eq(1)).Or(
                    coarse_fragments.gte(3).And(coarse_fragments.lt(15))).Or(
                    awc.lte(12).And(awc.gt(6))).Or(
                    permeability_class.eq(1)).multiply(3).rename(
                    'lcc3')).addBands(
                slope_perc.gt(10).And(slope_perc.lte(15)).And(K.multiply(100).divide(0.1317).gt(32)).Or(
                    slope_perc.gt(15).And(slope_perc.lte(30)).And(K.multiply(100).divide(0.1317).lte(32))).Or(
                    soil_properties.select('depth').gte(50).And(soil_properties.select('depth').lt(70))).Or(
                    awc.lt(6)).multiply(4).rename(
                    'lcc4')).addBands(
                coarse_fragments.gte(15).And(coarse_fragments.lt(50))
                      .multiply(5).rename('lcc5')).addBands(
                slope_perc.gt(15).And(slope_perc.lte(30)).And(K.multiply(100).divide(0.1317).gt(32)).Or(
                    K.multiply(100).divide(0.1317).lte(32)).Or(
                    soil_properties.select('depth').gt(20).And(soil_properties.select('depth').lt(50))).multiply(6).rename(
                    'lcc6')).addBands(
                slope_perc.gt(30).And(slope_perc.lte(60)).Or(
                    coarse_fragments.gte(50).And(coarse_fragments.lt(90))).multiply(7).rename(
                    'lcc7')).addBands(
                slope_perc.gt(60).Or(soil_properties.select('depth').lte(20)).Or(
                    coarse_fragments.gte(90)).multiply(8).rename('lcc8'))

            lcc = lcc.reduce(ee.Reducer.max()).rename('lcc')
            clipped = lcc.clip(self.clip_aoi)

            return clipped

        except Exception as e:
            logger.error(f'Calculating land capability class failed: {e}', exc_info=True)

    def get_google_buildings(self):
        """
        Fetches and processes Google Building dataset using Google Earth Engine (GEE).

        This function imports building data from the Google Open Building dataset,
        filters the data based on a confidence score of 0.6 or greater, applies a transformation
        to add properties to features, converts the features to an image, and then clips the
        result to the Area of Interest (AOI).

        Returns
        -------
        ee.Image
            An Earth Engine image object representing the processed and clipped building data.

        Raises
        ------
        Exception
            If the import or processing of the Google Building data fails.
        """

        def addNumberProperty(feature):
            return feature.set('n', ee.Number(1))
        logger.info(f'Importing Google Building data from {self.GOOGLE_BUILDING}')
        try:
            open_buildings = ee.FeatureCollection(self.GOOGLE_BUILDING).filterBounds(self.clip_aoi)
            open_buildings_filt = open_buildings.filter(ee.Filter.gte('confidence', 0.60))

            google_image = open_buildings_filt.map(addNumberProperty).reduceToImage(
                properties=['n'],
                reducer=ee.Reducer.first())

            mask_open_buildings = google_image.unmask(0)
            clipped = mask_open_buildings.clip(self.clip_aoi)

            return clipped
        except Exception as e:
            logger.error(f'Importing Google Building data from {self.GOOGLE_BUILDING} failed: {e}', exc_info=True)

    def get_wsf2019(self):
        logger.info(f'Importing World Settlement Footprint 2019 data from {self.WSF}')
        try:
            wsf2019 = ee.ImageCollection(self.WSF).mosaic()
            mask_wsf = wsf2019.unmask(0).eq(255).selfMask().unmask(0)
            clipped = mask_wsf.clip(self.clip_aoi)

            return clipped
        except Exception as e:
            logger.error(f'Importing World Settlement Footprint 2019 data from {self.WSF} failed: {e}', exc_info=True)

    def get_wsf_evo(self):
        """
        Imports and processes the World Settlement Footprint evolution data.

        This method attempts to import the World Settlement Footprint (WSF) evolution
        dataset using the Earth Engine API. The data is filtered to match the area of
        interest (AOI), then merged into a single image mosaic and clipped to the AOI
        geometry. If the operation fails, an error is logged.

        Returns:
            ee.Image: A clipped Earth Engine image representing the World Settlement
            Footprint evolution for the specified AOI.

        Raises:
            Exception: If there is any error during the import or processing of the
            dataset, it will log the issue and fail silently by re-raising the error.
        """
        logger.info(f'Importing World Settlement Footprint evolution data from {self.WSF_EVO}')

        try:
            wsf_evo = ee.ImageCollection(self.WSF_EVO).filterBounds(self.clip_aoi).mosaic()
            wsf_evo_clip = wsf_evo.clip(self.clip_aoi)
            return wsf_evo_clip
        except Exception as e:
            logger.error(f'Importing World Settlement Footprint evolution data from {self.WSF_EVO} failed: {e}', exc_info=True)

    def get_lppi(self, smoothing_years=3):
        """
            Calculate the Productivity performance indicator (LPPI)

            This function calculates the productivity performance indicator (LPPI) based on Sentinel-2
            imagery data and other geospatial data, such as land cover and soil properties.
            It derives a stratification unit from the combination of soil texture and land
            cover class. The function then uses cloud masking to process Sentinel-2 images
            and computes the mean Normalized Difference Vegetation Index (NDVI) over a
            specific time period. the productivity performance indicator aims to describe the vigor of existing
            vegetation in the city. The performance productivity indicator is calculated for each unique
            combination of landcover and soil texture within the AOI. For pixels within each unique
            landcover and soil texture combination, the mean Normalized Difference Vegetation Index (NDVI)
            is calculated for each pixel for the most recent year of Sentinel 2 observations and is divided by the
            95th percentile of mean NDVI across all pixels in the landcover/soil texture class such that higher
             values represent greener vegetation

            Args:
                smoothing_years (int): Number of years to smooth the NDVI data. Default is 3.

            Returns:
                ee.Image: An Earth Engine image representing the LPPI.

            Raises:
                Exception: If an error occurs during the calculation of the LPPI.
        """
        logger.info('Importing and calculating the LPPI index for clustering')
        try:
            land_cover = self.get_land_cover()
            if self.africa:
                soil_properties = self.get_soil_properties_africa()
            else:
                soil_properties = self.get_soil_properties()

            # Derive a stratification, which is a combination of soil texture and land cover class
            units = (soil_properties.select('texture').add(land_cover)).multiply(100).rename('units')

            # get sentinel imagecollection with cloud masking
            COLLECTION = self.S2_COLLECTION
            start_year = self.end_year - smoothing_years + 1
            start_date = f'{start_year}-01-01'
            end_date = f'{self.end_year}-12-31'
            s2_date_range = ee.Dictionary({'start': start_date, 'end': end_date})
            s2_cl = s2_masks.loadImageCollection(COLLECTION, s2_date_range, self.clip_aoi)
            masked_collection = cloud_mask(s2_cl)

            # get mean NDVI
            ndvi_mean = masked_collection.select('NDVI').reduce(ee.Reducer.mean(), 4).rename('NDVI')

            def stratified_mean(key, value):
                """
                Get stratified mean for key value pairs
                :param key:
                :param value:
                :return: the mean NDVI divided over the max for the key value stratus
                """
                ndvi_max = ee.Number(ee.Algorithms.If(
                    ee.Number(value).gte(50),
                    ndvi_mean.updateMask(units.eq(ee.Number.parse(key))).reduceRegion(
                        reducer=ee.Reducer.percentile([95]),
                        geometry=self.clip_aoi,
                        scale=100,
                        tileScale=16,
                        maxPixels=1e13).get('NDVI'),
                    0))
                return ndvi_mean.divide(ndvi_max).toFloat()

            ndvi_list = ee.Dictionary(units.reduceRegion(
                reducer=ee.Reducer.frequencyHistogram(),
                scale=100,
                maxPixels=1e13,
                tileScale=4,
                geometry=self.clip_aoi
            ).get('units')).map(lambda key, value: stratified_mean(key, value)).values()  # map over the gee dictionary key and value

            lppi = ee.ImageCollection.fromImages(ndvi_list).mosaic().rename('productivity_performance').clip(self.aoi.geometry())
            clipped = lppi.clip(self.clip_aoi)
            return clipped

        except Exception as e:
            logger.error(f'Importing and calculating the LPPI index for clustering failed: {e}', exc_info=True)

    def get_wri_trees(self):
        logger.info(f'Importing WRI tree percentage from {self.WRI_TREES}')
        try:
            trees = ee.ImageCollection(self.WRI_TREES).mosaic()
            clipped = trees.clip(self.clip_aoi)

            return clipped
        except Exception as e:
            logger.error(f'Importing WRI tree percentage from {self.WRI_TREES} failed: {e}', exc_info=True)

    def get_hrsl_pop(self):
        logger.info(f'Importing Population data from {self.HRSL_POP}')
        try:
            # population divided by 9 to account for the 30m pixel vs 10m pixel
            pop = ee.ImageCollection(self.HRSL_POP).median().divide(ee.Number(9))
            clipped = pop.clip(self.clip_aoi)

            return clipped
        except Exception as e:
            logger.error(f'Importing Population data from {self.HRSL_POP} failed: {e}', exc_info=True)

    def download_all(self, include_basic_indicators: bool = True, include_extended_indicators: bool = True, include_rusle: bool = True, include_LST: bool = True):
        images = []

        if include_basic_indicators:
            images.append({'image':self.get_dem(), 'name': 'dem', 'scale': 10, 'geometry': self.clip_aoi})
            images.append({'image': self.get_dem(), 'name': 'dem_30', 'scale': 30, 'geometry': self.clip_aoi})
            images.append({'image':self.get_slope(), 'name': 'slope', 'scale': 10, 'geometry': self.clip_aoi})
            images.append({'image': self.get_land_cover(), 'name': 'land_cover', 'scale': 10, 'geometry': self.clip_aoi})
            images.append({'image':self.get_bs_freq(), 'name': 'bs_freq', 'scale': 10, 'geometry': self.clip_aoi})
            images.append({'image': self.get_lppi(), 'name': 'lppi', 'scale': 10, 'geometry': self.clip_aoi})
            images.append({'image': self.get_hrsl_pop(), 'name': 'hrsl_pop', 'scale': 10, 'geometry': self.clip_aoi})
            images.append({'image': self.get_wsf2019(), 'name': 'wsf2019', 'scale': 10, 'geometry': self.clip_aoi})
            images.append({'image': self.get_wsf_evo(), 'name': 'wsf_evo', 'scale': 30, 'geometry': self.clip_aoi})
        
        if include_extended_indicators:
            images.append({'image':self.get_hand(), 'name': 'hand', 'scale': 10, 'geometry': self.clip_aoi})
            images.append({'image':self.get_annual_precipitation(), 'name': 'annual_prec', 'scale': 10, 'geometry': self.clip_aoi})
            images.append({'image':self.get_heavy_prec(), 'name': 'heavy_prec', 'scale': 10, 'geometry': self.clip_aoi})
            images.append({'image': self.get_lcc(), 'name': 'lcc', 'scale': 10, 'geometry': self.clip_aoi})
            images.append({'image': self.get_google_buildings(), 'name': 'google_buildings', 'scale': 10, 'geometry': self.clip_aoi})
            images.append({'image': self.get_water_freq(), 'name': 'water_freq', 'scale': 10, 'geometry': self.clip_aoi})
            images.append({'image': self.get_wri_trees(), 'name': 'wri_trees', 'scale': 10, 'geometry': self.clip_aoi})
            
            if self.africa:
                images.append({'image':self.get_soil_properties_africa(), 'name': 'soil_properties', 'scale': 10, 'geometry': self.clip_aoi})
                images.append({'image':self.get_soil_permeability_africa(), 'name': 'soil_perm', 'scale': 10, 'geometry': self.clip_aoi})
                images.append({'image':self.get_hydrosoil_group_africa(), 'name': 'hydrosoil', 'scale': 10, 'geometry': self.clip_aoi})
            else:
                images.append({'image': self.get_soil_properties(), 'name': 'soil_properties', 'scale': 10, 'geometry': self.clip_aoi})
                images.append({'image': self.get_hydrosoil_group(), 'name': 'hydrosoil', 'scale': 10, 'geometry': self.clip_aoi})
                images.append({'image': self.get_soil_permeability(), 'name': 'soil_perm', 'scale': 10, 'geometry': self.clip_aoi})
                
        # append rusle factors
        # request simplified K factor due to unresolved NaN issues
        if include_rusle:
            images.append({'image': RUSLE_factors.factorR(), 'name': 'rusle_r', 'scale': 10, 'geometry': self.clip_aoi})
            images.append({'image': self.get_factor_K_Benavidez(), 'name': 'rusle_k', 'scale': 10, 'geometry': self.clip_aoi})
            images.append({'image': RUSLE_factors.factorLS(self.get_slope(), self.get_slope_rad(), self.get_aspect()), 'name': 'rusle_ls', 'scale': 10, 'geometry': self.clip_aoi})
        
        if include_LST:
            lst = self.get_lst()
            images = images + lst

        return images
    
    
    def download_lppi(self):
        images = []        
        images.append({'image': self.get_lppi(), 'name': 'lppi', 'scale': 10, 'geometry': self.clip_aoi})       
        return images

    def check_exisiting(self, name: str = None) -> bool:
        """
        Check whether file already exists in the eo_indicators folder.
        Args:
            name: (str) Name of the file. should be storder in the image dictionary of each layer
        Returns:
            Boolean: True if the file  exists, False if it does not exist.
        """
        local_path = os.path.join(self.root, 'eo_indicators')
        current_files = os.listdir(local_path)
        if f'{name}.tif' in current_files:
            return True
        else:
            return False


    def export_ee(self, image, epsg, over_write=False, first_try=True):
        """
        Exports Earth Engine image data to a raster file in GeoTIFF format. This method downloads
        image data from Google Earth Engine, reprojects it to the specified EPSG, and saves the
        result into a file with LZW compression.

        Parameters:
        image: dict
            A dictionary containing the following keys:
            - name: str - Name of the output file (without extension).
            - image: str - An ID or reference to the Earth Engine Image to be processed.
            - geometry: dict - Geometry information specifying the area of interest.
            - scale: int - Desired resolution of the exported image.
        epsg: int
            EPSG code specifying the coordinate reference system to use for the exported
            image. This should be automatically created and refers to local UTM epsg
        over_write: bool, optional
            Whether to overwrite an existing file with the same name. Defaults to True.

        Raises:
        Exception
            An exception is raised if the image data could not be exported for any reason, such as
            invalid parameters or file IO issues.
        """
        try:

            if not over_write:
                # stop if file already exists
                cont = self.check_exisiting(image['name'])
                if cont:
                    tqauto.write(f' {image["name"]} already exists, skipping')
                    return
                
            ic = ee.ImageCollection(image['image']).map(
                # set the start time to avoid warning
                lambda img: img.set(
                    'system:time_start',
                    ee.Date.fromYMD(self.end_year, 1, 1).millis())
            )

            ds = xr.open_dataset(ic, engine='ee',
                                     crs=f'EPSG:{epsg}',
                                     geometry=image['geometry'],
                                     scale=image['scale'])
            imagr = ds.isel(time=0).rename({'Y': 'y', 'X': 'x'}).transpose('y', 'x').rio.write_crs(f"EPSG:{epsg}")
            output = os.path.join(self.root, 'eo_indicators', f'{image["name"]}.tif')
            imagr.rio.to_raster(output, compress='lzw')
            tqauto.write(f' Exported {image["name"]} to GeoTIFF format')

        #TODO paste the exception here
        except ValueError as e:
            if first_try:
                self.get_image_from_asset(image, epsg)
                tqauto.write(f' Exporting {image["name"]} failed: {e}, exporting from asset')
            else:
                tqauto.write(f'Could not export {image["name"]}: {e}')

        except Exception as e:
            tqauto.write(f'Could not export {image["name"]}: {e}')

    def export_to_asset(self, image, epsg,  max_pixels = 1e13):
        """
        Export an image to a Google Earth Engine (GEE) asset.

        This function facilitates exporting a specified image to a GEE asset using the Earth Engine API.
        It initializes and starts the export task. This method is used when direct export is failing.

        Args:
            image (dict): A dictionary containing details of the image to export. Expected keys include:
                - 'image': The Earth Engine image object
                - 'name': A string representing the name of the image/asset
                - 'scale': The resolution of the image in meters
                - 'geometry': The region to which the image will be exported
            epsg (int): EPSG code representing the coordinate reference system (CRS) of the export.
            max_pixels (int, optional): A threshold for the max number of pixels allowed in the export.
                Defaults to 1e13.

        Returns:
            ee.batch.Task.Task: The export task object initialized and started by the function.
        """
        asset_id = f'projects/{self.project}/assets/{image["name"]}'
        crs = f'EPSG:{epsg}'
        task = ee.batch.Export.image.toAsset(
            image=image['image'],
            description=image['name'],
            assetId=asset_id,
            scale=image['scale'],
            crs=crs,
            region=image['geometry'],
            maxPixels=max_pixels
        )

        task.start()

        logger.info(f"Started export task to asset: {asset_id}")
        return task


    def get_image_from_asset(self, image, epsg):
        """
        Get an image from an Earth Engine (EE) asset by exporting it and retrieving the processed
        result.

        This method repeatedly checks the status of the export task until it is either completed or
        encounters an error. If the task fails, it raises an exception with detailed information about
        the failure. Upon successful completion, the method retrieves the exported image exports it to
        a local geotif using the export.ee() method. This method is used when direct export is failing.

        Parameters:
        image (dict): A dictionary containing details of the image to export. Expected keys include:
                - 'image': The Earth Engine image object
                - 'name': A string representing the name of the image/asset
                - 'scale': The resolution of the image in meters
                - 'geometry': The region to which the image will be exported
        epsg: Coordinate reference system for the exported image.

        Raises:
        Exception: If the export task fails due to any error, with an error message describing
                   the failure in detail.
        Exception: If the export task encounters an unknown state, with information about the
                   unexpected task status.
        """
        task = self.export_to_asset(image, epsg)

        while True:
            task_status = task.status()['state']
            if task_status == 'READY':
                time.sleep(10)
            elif task_status == 'FAILED':
                raise Exception(f"Export task to asset {image['name']} failed: {task.status()['error_message']}")
            elif task_status == 'COMPLETED':
                break
            elif task_status == 'RUNNING':
                time.sleep(10)
            else:
                raise Exception(f"Export task to asset {image['name']} failed: {task.status()}")

        asset = task.config['assetExportOptions']['earthEngineDestination']['name']
        image2 = ee.Image(asset)
        image['image'] = image2

        self.export_ee(image, epsg, first_try=False)


def harmonize_ts(collection):
    """
    A function to group images per year for aggregation
    :param collection: ImageCollection
    :return: Grouped image collection
    """

    collection = collection.map(lambda img: img.set('year', ee.String(img.date().get('year'))))
    distinctYearCol = collection.distinct('year')
    filter = ee.Filter.equals(leftField='year', rightField='year')
    join = ee.Join.saveAll('year_matches')
    joinCol = ee.ImageCollection(join.apply(distinctYearCol, collection, filter))
    return joinCol

class ExportIndicators:

    def __init__(self, bucket, country, city, aoi):
        """

        :param bucket: (str)Name of GCS bucket
        :param country: (str)Name of Country
        :param city: (str)Name of City
        :param aoi: GEE Feature object of AOI
        """
        self.country = country
        self.city = city
        self.aoi = aoi
        self.bucket = bucket


    def export_image(self, image, description, region, folder='eo_indicators',
                         scale=10, crs='EPSG:4326'):
        """
        Create an Export GEE image to Google cloud batch task
        :param image: GEE image object
        :param description: Name for the file
        :param region: GEE geometry object to clip the image to
        :param folder: 'folder' in GCS city folder, or prefix
        :param scale: Reolution, default = 10m
        :param bucket: name of the bucket
        :param crs: crs, should be wgs84 (EPSG:4326)
        :return: Export task (run it by task.start())
        """
        gcs_folder = f"nbsos/{self.country}/{self.city}/{folder}"
        logger.info(f"Exporting {description} to GCS: {gcs_folder}/{description}")
        try:

            task = ee.batch.Export.image.toCloudStorage(image=image,
                                                        description=description,
                                                        fileNamePrefix=f'{gcs_folder}/{description}',
                                                        scale=scale,
                                                        bucket=self.bucket,
                                                        crs=crs,
                                                        region=region,
                                                        maxPixels=1e13)
            return task

        except Exception as e:
            logger.error(f'Exporting {description} to GCS: {gcs_folder}/{description} failed: {e}', exc_info=True)

    def export_image_to_drive(self, image, description, scale=10, crs='EPSG:4326', folder='eo_indicators'):
        """

        Exports an image to a Google Drive folder with specified parameters.

        This method facilitates exporting an Earth Engine image to a designated Google Drive
        folder. The export settings include image scale, projection information, and region
        of interest. The export parameters can be customized to suit specific requirements.

        Parameters:
            image: ee.Image
                The Earth Engine image to be exported.
            description: str
                A descriptive name for the export task.
            scale: int, optional
                The resolution of the exported image in meters per pixel. Defaults to 10.
            crs: str, optional
                The coordinate reference system (CRS) for the exported image. Defaults to 'EPSG:4326'.
            folder: str, optional
                The name of the Google Drive folder where the image will be saved. Defaults to 'eo_indicators'.

        Returns:
            ee.batch.Task
                The Earth Engine export task created for the export process.

        Raises:
            Exception
                If the export task fails due to issues during the process.
        """
        drive_folder = f"{self.city}_{folder}"
        logger.info(f"Exporting {description} to GCS: {drive_folder}/{description}")
        try:

            task = ee.batch.Export.image.toDrive(image=image,
                                                 folder=drive_folder,
                                                 description=description,
                                                 scale=scale,
                                                 crs=crs,
                                                 region=self.aoi.geometry(),
                                                 maxPixels=1e13)
            return task

        except Exception as e:
            logger.error(f'Exporting {description} to Google Drive: {drive_folder}/{description} failed: {e}', exc_info=True)


    def export_all(self, images, crs):
        tasks = []
        for image in images:
            try:
                task = self.export_image(image['image'], image['name'], region=image['geometry'], scale=image['scale'], crs=crs)
                task.start()
                tasks.append(task)
            except Exception as e:
                print(f'cannot export {image["name"]}')
        return tasks

def cloud_mask(col):
    """
    Filters and processes an Earth Engine ImageCollection to remove cloud and shadow artifacts, add vegetation
    indices, and select relevant spectral bands. The function applies a chain of necessary masks and
    transformations to ensure that the ImageCollection is ready for analysis.

    Parameters
    ----------
    col : ee.ImageCollection
        The input ImageCollection that includes Sentinel-2 imagery for further preprocessing.

    Returns
    -------
    ee.ImageCollection
        The processed and filtered ImageCollection with added vegetation indices, masks applied, and
        pre-selected spectral bands.
    """
    JRC_SURFACE_WATER = "JRC/GSW1_2/GlobalSurfaceWater"
    band_list = ['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B11', 'B12', 'fcover', 'NDVI']
    jrc_water = ee.Image(JRC_SURFACE_WATER).select('max_extent').eq(1)

    def mapper(img):
        out = img.addBands(
            img.normalizedDifference(
                ['B8', 'B4']).rename(
                'NDVI')).copyProperties(
            img, ['system:time_start'])
        return out

    out = col.filter(
        ee.Filter.notNull(['MEAN_INCIDENCE_AZIMUTH_ANGLE_B3',
                           'MEAN_INCIDENCE_AZIMUTH_ANGLE_B4',
                           'MEAN_INCIDENCE_AZIMUTH_ANGLE_B5',
                           'MEAN_INCIDENCE_AZIMUTH_ANGLE_B6',
                           'MEAN_INCIDENCE_AZIMUTH_ANGLE_B7',
                           'MEAN_INCIDENCE_AZIMUTH_ANGLE_B8A',
                           'MEAN_INCIDENCE_AZIMUTH_ANGLE_B11',
                           'MEAN_INCIDENCE_AZIMUTH_ANGLE_B12',
                           'MEAN_INCIDENCE_ZENITH_ANGLE_B3',
                           'MEAN_INCIDENCE_ZENITH_ANGLE_B4',
                           'MEAN_INCIDENCE_ZENITH_ANGLE_B5',
                           'MEAN_INCIDENCE_ZENITH_ANGLE_B6',
                           'MEAN_INCIDENCE_ZENITH_ANGLE_B7',
                           'MEAN_INCIDENCE_ZENITH_ANGLE_B8A',
                           'MEAN_INCIDENCE_ZENITH_ANGLE_B11',
                           'MEAN_INCIDENCE_ZENITH_ANGLE_B12',
                           'MEAN_SOLAR_AZIMUTH_ANGLE',
                           'MEAN_SOLAR_ZENITH_ANGLE'])
    ).map(
        s2_masks.addCloudShadowMask(jrc_water.Not(), 1e4)).map(
        s2_masks.applyCloudShadowMask).map(
        S2Fcover.fcover(1e4)).map(
        mapper).select(band_list)
    return out

