import logging
import rasterio as rio
import os
import numpy as np
from matplotlib import pyplot


class InvestFlood:
    """
    Calculates runoff, retention, runoff volume and retention volume. Input needed
    is storm depth in mm and a raster with CN-numbers. All equations are based on
    Urban Invest flood mitigation module, found at
    https://invest-userguide.readthedocs.io/en/latest/urban_flood_mitigation.html

    """

    def __init__(self, name, cn_path, storm_depth, output_folder):

        self.output_folder = output_folder
        self.name = name
        self.cn, self.meta = open_raster(cn_path)
        self.p = storm_depth
        self.smax = None
        self.runoff = None
        self.retention = None
        self.retention_volume = None
        self.runoff_volume = None
        lg = Logger()
        self.__logger = lg.get()


    def _calc_smax(self, p_1=25400, p_2=254):
        """
        Invest function to calculate Smax as a function of Curve Number,
        based on empirical equation found at
        https://invest-userguide.readthedocs.io/en/latest/urban_flood_mitigation.html
        :param cn: np array with Curve Number
        :param f: empirical parameter st at 25400 in Invest
        :return: np array with Smax
        """

        self.__logger.info(f"Calculating s(max)")

        try:
            c = zero_divide(p_1, self.cn)
            smax = c - p_2
            self.smax = RasterImage(smax, self.meta)

            return smax

        except Exception as err:
            self.__logger.error("could not calculate S(max)")
            raise err


    def _calc_runoff(self, lam=0.2):
        """
        Conditional function from InVest to calculate runoff as a function of
        Curve Number.
        https://invest-userguide.readthedocs.io/en/latest/urban_flood_mitigation.html
        If storm depth (p) > Smax * Lambda runoff is 0
        :param p: Storm depth in mm
        :param lam: Lambda parameter, set default to 0.2 by InVest
        :param smax: np array with Smax calculated by calc_smax
        :return: np array with runoff in mm
        """

        self.__logger.info(f"Calculating runoff, lambda = {lam}")

        try:
            smax = self._calc_smax()
            new_smax = smax * lam
            a = (self.p - new_smax) ** 2
            b = self.p + (1 - lam) * smax

            q = zero_divide(a,b)
            q2 = np.where((self.p > new_smax) & (q > 0), q, 0)
            self.runoff = RasterImage(q2, self.meta)

            return q2

        except Exception as err:
            self.__logger.error("Could not calculate runoff")
            raise err

    def _get_cell_size(self):
        """
        Get pixel size from meta data of the cn input raster. I'm not sure if this
        covers all or most, but it will work with data in meters and decimal degrees
        It does assume square raster cells
        :return: pixel size in square meters
        """

        self.__logger.info(f"Calculating pixel size")

        try:
            # use x-dimension
            pixel_size = abs(self.meta['transform'][0])
            

            if pixel_size < 0.5: # units are in degrees
                # get extent along y-dimension, but return absolute value 
                pixel_size = abs(self.meta['transform'][4])
                self.pixel_size = (pixel_size * 111.32 * 10**3) ** 2

            else: # units are in meters
                self.pixel_size = pixel_size

            self.__logger.info(f"Pixel size = {self.pixel_size} m^2")


            return pixel_size

        except Exception as err:
            self.__logger.error("Could not calculate pixel size")
            raise err



    def _calc_run_retention(self):
        """
        Calcuate the runoff retention per pixel
        :return: a raster image with values of runoff retention per pixel
        """

        self.__logger.info(f"Calculating runoff retention")

        try:
            temp = zero_divide(self.runoff, self.p)
            r = 1 - temp
            self.retention = RasterImage(r, self.meta)
            return r

        except Exception as err:
            self.__logger.error("Could not calculate runoff retention")
            raise err


    def _calc_run_retention_volume(self):
        """
        Calcuate the runoff retenion in volumne
        :return: raster image with runoff retention in cubic meters
        """

        self.__logger.info(f"Calculating runoff retention volume")

        try:
            r_m3 = self.retention * self.p * self.pixel_size * 10 ** -3
            self.retention_volume = RasterImage(r_m3, self.meta)
            return r_m3

        except Exception as err:
            self.__logger.error("Could not calculate runoff retention volume")
            raise err

    def _calc_runoff_volume(self):
        """
        Calculate the runoff volumne (flood volumne) in cubic meter
        :return:
        """

        self.__logger.info(f"Calculating runoff volume")

        try:
            q_m3 = self.runoff * self.pixel_size * 10 ** -3
            self.runoff_volume = RasterImage(q_m3, self.meta)

            return q_m3

        except Exception as err:
            self.__logger.error("Could not calculate runoff  volume")
            raise err


    def run(self):
        """
        Run the Urban Flood Invest model (without damage function)
        :return: -
        """

        self._calc_runoff()
        self._get_cell_size()
        self._calc_run_retention()
        self._calc_run_retention_volume()
        self._calc_runoff_volume()

    def export_layers(self, runoff=False, smax=False, retention=False,
                      runoff_volume=True, retention_volume=False):
        """
        Export tif layers of output files
        :param runoff: set to True to export runoff tif
        :param smax: set to True to export smax tif
        :param retention: set to True to export retention tif
        :param runoff_volume: set to True to export runoff_volume tif
        :param retention_volume: set to True to export retention_volume
        :return:
        """
        create_folder(self.output_folder)

        export_list = [key for key, value in locals().items() if value][1:]

        self.__logger.info(f"Exporting the following layers:"
                           f"{export_list}")

        for layer in export_list:
            filename = os.path.join(self.output_folder, f'{self.name}_{layer}.tif')

            try:
                self.__getattribute__(layer).export(filename)
                self.__logger.info(f"{filename} created")

            except Exception as err:
                self.__logger.warning(f"Could not create {filename}"
                                      f"{err}")
                continue


class RasterImage(np.ndarray):

    def __new__(cls, a, meta):
        """
        add plotting and exporting function to np.array
        :param a:
        :param meta: a copy of meta data from template raster tif
        :return: np.ndarray with additional functions
        """

        obj = np.asarray(a).view(cls)
        obj.meta = meta

        return obj

    def imshow(self, **kwargs):
        """
        plot raster layer
        :param kwargs: see pyplot.imshow for input parameters
        :return:
        """

        pyplot.imshow(self, **kwargs)

    def export(self, out_ras):
        """
        Use rasterio to export raster to tif
        :param out_ras: "path to output raster
        :return:
        """

        self.meta.update({'dtype': self.dtype}) # ascertain correct dtype
        with rio.open(out_ras, 'w', **self.meta) as dst:
            dst.write(self, 1)


def open_raster(in_ras, **kwargs):
    """
    Function to get np.ndarray from tif with meta data for template
    :param in_ras: input tif
    :param kwargs: any pararmters to update metadata
    :return:
    """

    with rio.open(in_ras, 'r+') as rst:
        data = rst.read(1)
        rst.nodata = 0
        meta = rst.meta.copy()
        meta.update(compress='lzw', **kwargs)

    return data, meta


def zero_divide(a,b):
    """
    Divide rasters and work araound zero divides. For now it gives 0 as result.
    :param a: numerator
    :param b:denominator
    :return: output raster with 0 as divide by 0 result
    """

    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(a, b)
        c[c == np.inf] = 0
        c = np.nan_to_num(c)

    return c


def create_folder(folder_name):
    try:
        os.makedirs(folder_name)
    except FileExistsError:
        pass


class Logger:

    def __init__(self, level=logging.INFO):

        logfmt = "[%(asctime)s - %(levelname)s] - %(message)s"
        dtfmt = "%Y-%m-%d %I:%M:%S"
        logging.basicConfig(level=logging.CRITICAL, format=logfmt, datefmt=dtfmt)
        self.logger = logging.getLogger()
        self.logger.setLevel(level=level)

    def get(self):

        return self.logger

