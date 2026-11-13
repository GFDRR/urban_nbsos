from . import gis_tools

from skimage.measure import label, regionprops_table
from skimage.util import map_array
from astropy.convolution import convolve
from skimage import draw
import rasterio
import numpy as np
import matplotlib.pyplot as plt
import os

import logging
logger = logging.getLogger('nbsos.scan')


def focal_stat(in_raster, radius, norm=False ,in_file = 'tif', nan_treatment='fill'):
    """
    Function to perform focal statistics over raster
    :param in_raster: path to input raster
    :param radius: radius (in cells)
    :param mode: choose sum or mean
    :return: numpy array with focal statistics output
    """
    if in_file == 'tif':
        # convert in_raster to np array
        with rasterio.open(in_raster) as in_ras:
            in_array = in_ras.read(1)
    else:
        in_array = in_raster

    kernel = create_kernel(radius)
    in_array = in_array.astype('float')
    in_array[in_array == 0] = np.nan
    # execute convolution with kernel
    out_array = convolve(in_array, kernel, normalize_kernel=norm, nan_treatment=nan_treatment)
    foc0 = np.nan_to_num(out_array)
    return foc0


def create_kernel(radius=4):
    """
    Function to be used by focal statistics to create the moving window kernel
    :param radius: radius in cells
    :return: kernel
    """
    size = (radius * 2) + 1
    arr = np.zeros((size, size))
    rr, cc = draw.disk((radius, radius), radius+1)
    arr[rr, cc] = 1
    return arr

def label_image(image, min_size):
    l = label(image)
    table = regionprops_table(
        l,
        properties=('label', 'area')
    )
    condition = (table['area'] >= min_size)

    # zero out labels not meeting condition
    input_labels = table['label']
    output_labels = input_labels * condition

    new_image = map_array(
        l, input_labels, output_labels
    )
    return new_image

def get_300_demand(local_root, country, city, nbs, weight, min_size=100):
    """
    Function to find the demand for greenspace based on the 300m rule from Konijnendijk
    :param status_mask: raster tif with value 1 for existing greenspace,
    :param pop: rater tif with population
    :param creation_mask: raster tif with 1 for creation opportunities of nbs
    :param weight: weight for recreation value of the nbs
    :param min_size: minimum size of existing greenspace patches to count in pixels. Konijnendijk has 1ha (100pixels)
    :return: demand for recreation and priority areas.
    """
    status_mask = os.path.join(local_root, country, city, 'suitability', 'nbs_protection.tif')
    pop = os.path.join(local_root, country, city, 'eo_indicators', 'hrsl_pop.tif')
    creation_mask = os.path.join(local_root, country, city, 'suitability', f'{nbs}_creation.tif')
    with rasterio.open(status_mask) as src:
        # protection raster, i.e. where there is current greenspace
        nbs_p = src.read(1)
    with rasterio.open(pop) as src:
        # population raster
        popul = src.read(1)
    with rasterio.open(creation_mask) as src:
        # where nbs is defined as creation opportunity
        nbs_c = src.read(1)

    protection = np.where(nbs_p == 1, 1, 0) # make sure that we have the correct value for protection
    new_protection = label_image(protection, min_size=min_size) # only keep those patches that are larger than the min_size in pixels
    prot_foc = focal_stat(new_protection, 29, in_file='np')
    no_access = np.where(prot_foc == 0, 1, 0) # where there is no access to 300m park

    """
    gis_tools.export_raster(no_access, pop, 'c:/tmp/temp_ras.tif')
    with rasterio.open('c:/tmp/temp_ras.tif') as src:
        no_access = src.read(1)
    """

    # population that needs park access
    pop_in_need = popul * no_access

    #How many people can be provided access to greensapce by creating nbs here
    # this will be the demand for recreation services in the nbsos ppt
    pop_300 = focal_stat(pop_in_need, 29, in_file='np')

    # the demand filtered for nbs creation opportunities multiplied byt the weight for the nbs
    nbs_potential = nbs_c * pop_300 * weight

    return nbs_potential, pop_300


def get_300_protection(local_root, country, city, out_path, protection_mask_type: str, min_size=100):
    """
    Analyze greenspace protection and population serviceability within a 300m radius.

    This function processes raster data to identify green spaces (protection)
    based on a given minimum size threshold and calculates the population
    serviced within a 300m radius by these green spaces. The output is exported
    as a new raster file depicting the population benefits of the identified
    protection areas.

    Parameters:
    ----------
    local_root : str
        The base directory of the input and output data structure.
    country : str
        The name of the country used to specify the folder structure.
    city : str
        The name of the city used to specify the folder structure.
    out_path : str
        The path of the output raster to be saved.
    protection_mask_type : str
        The type of protection mask to use, either "suitability" or "lppi", which determines the source raster for protection areas.
    min_size : int, optional
        The minimum patch size in pixels to consider as valid protection
        (default is 100).

    Raises:
    ------
    FileNotFoundError
        If the specified raster files do not exist in the provided directory
        structure.
    ValueError
        If the input raster data is malformed or incompatible with the process.
    """
    
    if protection_mask_type == "suitability":
        status_mask = os.path.join(local_root, country, city, 'suitability', 'nbs_protection.tif')
    elif protection_mask_type == "lppi":
        status_mask = os.path.join(local_root, country, city, 'suitability', 'nbs_protection_clustering.tif')
    
    
    pop = os.path.join(local_root, country, city, 'eo_indicators', 'hrsl_pop.tif')

    with rasterio.open(status_mask) as src:
        # protection raster, i.e. where there is current greenspace
        nbs_p = src.read(1)
    with rasterio.open(pop) as src:
        # population raster
        popul = src.read(1)

    if protection_mask_type == "suitability":
        protection = np.where(nbs_p == 1, 1, 0) # make sure that we have the correct value for protection
    elif protection_mask_type == "lppi":
        protection = np.where(nbs_p == 1, 1, 0)  # make sure that we have the correct value for protection
        
    
    new_protection = label_image(protection, min_size=min_size) # only keep those patches that are larger than the min_size in pixels
    new_protection = np.where(new_protection >0, 1, 0)
    #How many people are serviced within 300m
    pop_300 = focal_stat(popul, 29, in_file='np')
    pop_serviced = np.where(new_protection == 1, pop_300, 0)
    gis_tools.export_raster(pop_serviced, pop, out_path)


