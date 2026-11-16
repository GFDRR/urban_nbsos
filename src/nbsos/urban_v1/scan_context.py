import geopandas as gpd
from typing import Tuple

class ScanContext:
    """This class acts as a container to set relevant variables for defining the spatial context of a scan.
    """
    _local_root: str = None
    _country_name: str = None
    _city_name: str = None
    _location_in_africa: bool = None
    _area_of_interest: str = None
    _aoi_centroid: any = None

    def __init__(self, local_root: str = None, country: str = None, city: str = None, location_in_africa: bool = False, area_of_interest: str = None):
        self.local_root = local_root
        self.country = country
        self.city = city
        self.location_in_africa = location_in_africa
        self.local_area_of_interest_file = area_of_interest

    @property
    def local_root(self) -> str:
        """The path to the local root folder, i.e., the data repository used by the scan.

        :getter: Returns the path to the local root folder.
        :setter: Sets the path to the local root folder.
        :type: str
        """
        return self._local_root
    
    @local_root.setter
    def local_root(self, value: str) -> None:
        self._local_root = value

    @property
    def country(self) -> str:
        """Name of the country of the NBSOS case to be assessed.

        :getter: Returns the name of the country of the NBSOS case to be assessed.
        :setter: Sets the name of the country of the NBSOS case to be assessed.
        :type: str
        """
        return self._country_name
    
    @country.setter
    def country(self, value: str) -> None:
        self._country_name = value

    @property
    def city(self) -> str:
        """Name of the city of the NBSOS case to be assessed.

        :getter: Returns the name of the city of the NBSOS case to be assessed.
        :setter: Sets the name of the city of the NBSOS case to be assessed.
        :type: str
        """
        return self._city_name
    
    @city.setter
    def city(self, value: str) -> None:
        self._city_name = value

    @property
    def location_in_africa(self) -> bool:
        """Defines whether the location of the NBSOS case to be assessed is within the African continent or not. It set to True, uses certain additional datasources in the preparation of indicators. By default, False. 

        :setter: Sets the indicator if the location is in Africa.
        :getter: Returns the indicator if the location is in Africa.
        :type: bool
        """
        return self._location_in_africa
    
    @location_in_africa.setter
    def location_in_africa(self, value: bool) -> None:
        self._location_in_africa = value

    @property
    def local_area_of_interest_file(self) -> str:
        """Path to the shapefile of the area of interest (AOI) of the NBSOS case to be assessed.

        :getter: Returns the path to the shapefile of the area of interest (AOI).
        :setter: Sets the path to the shapefile of the area of interest (AOI).
        :type: str
        """
        return self._area_of_interest
    
    @local_area_of_interest_file.setter
    def local_area_of_interest_file(self, value: str) -> None:
        # store the path to the shapefile of the AOI
        self._area_of_interest = value
    
    def get_centroid_wgs84(self) -> Tuple[float, float]:   
        """Return the x and y coordinates of the AOI centroid as tuple. The centroid is constructed for the union of all polygon geometries in the AOI shapefile.

        :return: WGS84 x and y coordinates of the AOI centroid.
        :rtype: Tuple[float, float]
        """

        # return the centroid of the area of interest in WGS84
        df_aoi = gpd.read_file(self.local_area_of_interest_file)
        df_aoi = df_aoi.to_crs(4326)
        centroid = df_aoi.geometry.union_all().centroid 
        return centroid.x, centroid.y

    @staticmethod
    def in_notebook() -> bool:
        """Identify if the environment the scan is executed in is a Jupyter notebook.

        :return: Indicator if the scan is executed in a Jupyter notebook.
        :rtype: bool
        """
        try:
            from IPython import get_ipython
            if 'IPKernelApp' not in get_ipython().config:  # pragma: no cover
                return False
        except ImportError:
            return False
        except AttributeError:
            return False
        return True