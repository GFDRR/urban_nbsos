import rasterio as rio
from rasterio import features
import os
import geopandas as gpd
from scipy import sparse
from rasterio.warp import reproject, Resampling, calculate_default_transform
import numpy as np
from astropy.convolution import convolve
from skimage import draw
import math
import pandas as pd
from typing import Tuple, List, Dict
import exactextract as ee
from rasterio import Affine
import xarray as xr


def clip(rst_land_cover: rio.DatasetReader, clip_shape: gpd.GeoDataFrame, ) -> np.ndarray:
    """Clip an input raster to the extent of a polygon contained in clip_shape.

    :param rst_land_cover: The input raster to clip.
    :type rst_land_cover: rasterio.DatasetReader
    :param clip_shape: The shapefile containing the extent to clip as polygon feature.
    :type clip_shape: str
    :return: _description_
    :rtype: np.ndarray
    """
    clipped_land_use_matrix = rio.mask.mask(rst_land_cover, clip_shape['geometry'], indexes=1)
    return clipped_land_use_matrix

def make_exact_zonal_statistics(zone_shp: gpd.GeoDataFrame, ds_values: xr.DataArray, statistics: Dict[str,str]) -> gpd.GeoDataFrame:
    """Compute zonal statistics for array over features (that define zones) in a GeoDataFrame, and append the GeoDataFrame with computed statistics. 
       Computed indicators are weighted by area of intersection between zone and raster cell. 

    :param zone_shp: Features defining zones  
    :type zone_shp: gpd.GeoDataFrame
    :param value_mtx: Array of values for which statistics will be computed per zone
    :type value_mtx: np.ndarray
    :param transform: The Affine transform of the value array
    :type transform: Affine
    :param statistics: Dictionary of statistics to be computed (dictionary vale), and corresponding field names to be written to the GeoDataFrame (dictionary key)
    :type statistics: Dict[str,str]
    :return: GeoDataFrame of features with appended zonal statistics 
    :rtype: gpd.GeoDataFrame
    """

    stat_params = [f for f in statistics.values()]
    stats = ee.exact_extract(ds_values, zone_shp, stat_params, output='pandas')
    for field_name, param_name in statistics.items():
        zone_shp[field_name] = stats[param_name]

    return zone_shp


def get_area_crs(area_of_interest_filename: str) -> int:
        
    # Note: This method by Koen        
    # make sure the input is in wgs 84
    aoi = gpd.read_file(area_of_interest_filename).to_crs(4326)
    
    lon = aoi.geometry.centroid[0].x
    lat = aoi.geometry.centroid[0].y
    
    """Based on lat and lng, return best utm epsg-code"""
    utm_band = str((math.floor((lon + 180) / 6) % 60) + 1)
    
    if len(utm_band) == 1:
        utm_band = '0'+utm_band
    
    if lat >= 0:
        epsg_code = '326' + utm_band
        return int(epsg_code)
    
    epsg_code = '327' + utm_band        
    return int(epsg_code)   

def make_gaussian_kernel(d: float, sigma: float, rs: Tuple[int,int]) -> np.ndarray:
    """This function implements the R's terra::focalMat function of type Gauss to make a 
    Gaussian focal weights matrix. The size of the kernel that is returned is determined 
    by the specified kernel radius and resolution. The kernel is normalized so that the 
    sum of all elements is 1. 

    :param d: Radius of the kernel in meters
    :type d: float
    :param sigma: Sigma used in the Gaussian formula
    :type sigma: float
    :param rs: Raster resolution
    :type rs: Tuple[int,int]
    :return: Matrix with Gaussian kernel
    :rtype: np.ndarray
    """

    # determine kernel shape based on resolution
    # assert positive values
    nx = 1 + 2 * math.floor( d / abs(rs[0]))
    ny = 1 + 2 * math.floor( d / abs(rs[1]))

    # make result matrix: This is used to iterate over cells
    m = np.zeros((nx, ny), dtype=np.float32)

    # determine half resolution
    hrx = rs[0] / 2
    hry = rs[1] / 2

    # make coordinates of grid cell centroids with respect to total extent (d) of kernel. 
    # Start of iteration is at minus the half resolution, so we can add to the initial coordinates multiples of resolution  
    xr = (nx * rs[0]) / 2 * -1 - hrx
    yr = (ny * rs[1]) / 2 * -1 - hry

    # np has row,col indexing!
    # This works and appears to reproduce R nicely, but there might be better methods such as intializing an DataArray.
    # So there might be potential for improvement here.
    centroid_coordinates = []
    for index, x in np.ndenumerate(m):
    
        crow = index[0] + 1 # add one to offset zero-based indexing
        ccol = index[1] + 1 # add one to offset zero-based indexing^

        xcor = xr + (ccol * rs[0])
        ycor = yr + (crow * rs[1])        

        centroid_coordinates.append([xcor, ycor])

    # make a dataframe of centroid coordinates    
    df = pd.DataFrame(centroid_coordinates, columns=['x', 'y'])
    
    # square all columns
    dfsqr = df.apply(lambda x: x**2)
    
    # apply 2D Gaussian filter formula, that is also used in terra::focalMat.
    dfsqr['m'] = 1/(2 * math.pi * (sigma**2)) * np.exp(-(dfsqr['x'] + dfsqr['y'])/(2 * (sigma**2)))

    # reshape to matrix, and normalize so that sum of all elements is 1
    m = dfsqr['m'].values.reshape(nx, ny)
    m = m/np.sum(m)

    return m

def export_raster_dataset(ds: xr.Dataset, out_path: str):
    """Export xarray dataset to raster file.

    :param ds: Dataset to be exported
    :type ds: xr.Dataset
    :param out_path: Path to output raster
    :type out_path: str
    """
    
    # assert that target directory exists
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    # assert scale factor is 1
    attr = ds.attrs.copy()
    attr.update({'scale_factor' : 1})
    ds.rio.set_attrs(attr, inplace=True)
    ds.rio.to_raster(out_path)
    

def export_raster(in_band, mask, out_raster, meta=None):
    """
    Function to export numpy array to geotiff
    :param in_band: numpy array to be used as band 1
    :param mask: template raster mask
    :param out_raster: path to output raster
    :return: -
    """
    
    rst_meta = template(mask) if meta is None else meta
    rst_meta['compress'] = 'lzw' # asert also for custom meta
    with rio.open(out_raster,  'w', **rst_meta) as dst:
        dst.write(in_band, 1)

def template(template_raster, **kwargs):
    """
    return the meta data of a template raster to be used to create new rasters
    :param template_raster: path to template raster (geotif)
    :param kwargs: extra arguments to be passed to meta.update
    :return: meta data of template raster
    """
    with rio.open(template_raster, 'r+') as rst:
        #rst.nodata = 0
        meta = rst.meta.copy()
        meta.update(compress='lzw', **kwargs)
    return meta

def get_points(in_ras, key):

    # get geo information from flood pop raster
    with rio.open(in_ras) as src:
        gt = src.transform
        arr = src.read(1)

    # sparse matrix to exclude zeros (scipy)
    M = sparse.coo_matrix(arr.T)
    vals = M.data  # values
    ys = M.col  # x coordinates
    xs = M.row  # y_coorinates

    # get coordinates using gt values to transform to wgs_84
    coords = [
        {'x': x * gt[0] + gt[2],
         'y': y * gt[4] + gt[5],
         key: val}
        for x, y, val in zip(xs, ys, vals)
        if val > 0]

    return coords

def rasterize(shape_path, templ, out_ras, attr, fill_value: float = 0, nodata_value: float = np.nan):
    polys = gpd.read_file(shape_path)
    # assert projection 
    tmp = rio.open(templ)
    polys = polys.to_crs(tmp.crs)
    rasterize_df(polys, templ, out_ras, attr, fill_value, nodata_value)



def rasterize_geom(polys: gpd.GeoDataFrame, templ: str, out_ras: str, burn_value: float = 1, fill_value: float = 0, nodata_value: float = np.nan):
    meta = template(templ)
    
    # update meta with nodata value
    meta.update({
        'nodata' : nodata_value,
        'dtype' : 'float32'
    })

    with rio.open(out_ras, 'w+', **meta) as out:
        out_arr = out.read(1)

        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom, burn_value) for geom in polys.geometry)        
        burned = features.rasterize(shapes=shapes, fill=fill_value, out=out_arr, transform=out.transform)
        out.write_band(1, burned)

def rasterize_df(polys: gpd.GeoDataFrame, templ: str, out_ras: str, attr: str, fill_value: float = 0, nodata_value: float = np.nan):
    meta = template(templ)
    
    # update meta with nodata value
    meta.update({
        'nodata' : nodata_value,
        'dtype' : 'float32'
    })

    with rio.open(out_ras, 'w+', **meta) as out:
        out_arr = out.read(1)
        
        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom,value) for geom, value in zip(polys.geometry, polys[attr]))
        burned = features.rasterize(shapes=shapes, fill=fill_value, out=out_arr, transform=out.transform)
        out.write_band(1, burned)


def resample_raster(input, output, mask, temp_folder):
    temp_output = os.path.join(temp_folder, 'temp_ras.tif')

    # read in the mask variables
    with rio.open(mask) as src:
        xres = src.res[0]
        yres = src.res[1]
        left, bottom, right, top = src.bounds
        trans = src.transform
    window = rio.windows.from_bounds(left, bottom, right, top, trans)

    # read in the input raster
    with rio.open(input) as dataset:
        scale_factor_x = dataset.res[0]/xres
        scale_factor_y = dataset.res[1]/yres

        profile = dataset.profile.copy()
        # resample data to target shape
        data = dataset.read(
            out_shape=(dataset.count,
                       dataset.height * scale_factor_y,
                       dataset.width * scale_factor_x)
        )

        # scale image transform
        transform = dataset.transform * dataset.transform.scale(
            (1 / scale_factor_x),
            (1 / scale_factor_y)
        )
        profile.update({"height": data.shape[-2],
                        "width": data.shape[-1],
                        "transform": transform})

    with rio.open(temp_output, "w", **profile) as dataset:
        dataset.write(data)

    with rio.open(temp_output) as src:
        new_array = src.read(1, window=window)

    export_raster(new_array, mask, output)

"""
infile = os.path.join(local_root, country, 'fathom', 'P_1in5.tif')
outfile = os.path.join(local_root, country, 'fathom', 'P_1in_5_10m.tif')
match = os.path.join(local_root, country, city, 'eo_indicators', 'land_cover.tif')
output = os.path.join(local_root, country, 'fathom', 'P_1in_5_10m2.tif')

def resample_raster(input, mask, temp_folder):

    # read in the mask variables
    with rio.open(mask) as src:
        xres = src.res[0]
        yres = src.res[1]

    # read in the input raster
    with rio.open(input) as dataset:
        scale_factor_x = dataset.res[0]/xres
        scale_factor_y = dataset.res[1]/yres

        profile = dataset.profile.copy()
        # resample data to target shape
        data = dataset.read(
            out_shape=(dataset.count,
                       int(dataset.height * scale_factor_y),
                       int(dataset.width * scale_factor_x))
        )

        # scale image transform
        transform = dataset.transform * dataset.transform.scale(
            (1 / scale_factor_x),
            (1 / scale_factor_y)
        )
        profile.update({"height": data.shape[-2],
                        "width": data.shape[-1],
                        "transform": transform})

    with rio.open(out, "w", **profile) as dataset:
        dataset.write(data)

def mask_raster(input, mask, output):
    with rio.open(mask) as src:
        left, bottom, right, top = src.bounds
        trans = src.transform

    window = rio.windows.from_bounds(left, bottom, right, top, trans)
    with rio.open(input) as src:
        new_array = src.read(1, window=window)

    export_raster(new_array, mask, output)
"""
def reproj_match(infile, match, outfile, resampling_method: Resampling = Resampling.nearest):
    """Reproject a file to match the shape and projection of existing raster.
    SOURCE: https://pygis.io/docs/e_raster_resample.html
    Parameters
    ----------
    infile : (string) path to input file to reproject
    match : (string) path to raster with desired shape and projection
    outfile : (string) path to output file tif
    """
    # open input
    with rio.open(infile) as src:
        src_transform = src.transform

        # open input to match
        with rio.open(match) as match:
            dst_crs = match.crs

            # calculate the output transform matrix
            dst_transform, dst_width, dst_height = calculate_default_transform(
                src.crs,  # input CRS
                dst_crs,  # output CRS
                match.width,  # input width
                match.height,  # input height
                *match.bounds,  # unpacks input outer boundaries (left, bottom, right, top)
            )

        # set properties for output
        dst_kwargs = src.meta.copy()
        dst_kwargs.update({"crs": dst_crs,
                           "transform": dst_transform,
                           "width": dst_width,
                           "height": dst_height,
                           "nodata": 0})

        with rio.open(outfile, "w", **dst_kwargs) as dst:
            # iterate through bands and write using reproject function
            for i in range(1, src.count + 1):
                reproject(
                    source=rio.band(src, i),
                    destination=rio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=dst_transform,
                    dst_crs=dst_crs,
                    resampling=resampling_method)


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
        with rio.open(in_raster) as in_ras:
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
