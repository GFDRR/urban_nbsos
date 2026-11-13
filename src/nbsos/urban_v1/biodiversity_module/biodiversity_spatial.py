from .biodiversity_exports import HabitatThreat, ThreatProperties, ThreatDecayFunction, LandCoverClassVegetationProperties
from .biodiversity_logging import biodiversity_logger

import math
import os
from pathlib import Path
import xarray as xr 
from xrspatial import proximity  
import numpy as np
from typing import List, Tuple
import geopandas as gpd
import rasterio
import rasterio.features
from scipy import ndimage
import pandas as pd
from rich.progress import Progress, TaskProgressColumn, TimeElapsedColumn, MofNCompleteColumn, TextColumn, BarColumn
import concurrent.futures
from enum import Enum
import numpy.ma as ma
from skimage.measure import label, regionprops, regionprops_table
from rasterio.warp import calculate_default_transform, reproject, Resampling
from scipy import ndimage
import sys
from skimage.util import map_array
import rioxarray

class PatchContiguity(Enum):
        Queen = 0
        Rook = 1

class Spatial:
    """This class provides methods for the treatment of spatial raster data, including methods for the filtering of habitat,
       identification of connected patches, computation of proximity rasters, decay surfaces, etc. 
    """

    @staticmethod
    def get_distance_raster(mtx_input: np.ndarray, resolution: Tuple[int,int]) -> np.ndarray:        
        n, m = mtx_input.shape
        xr_rst = xr.DataArray(mtx_input, dims=['y', 'x'], name='raster')
        xr_rst['y'] = np.arange(n)[::-1]
        xr_rst['x'] = np.arange(m)
        mtx_res = proximity(xr_rst, target_values = [1]).to_numpy()
        
        # mtx_res is now in units of pixels. we scale by resolution to have a raster of distance
        mtx_res = mtx_res * resolution[0]
        
        return mtx_res

    @staticmethod
    def get_threat_decay_surface(threat: ThreatProperties, mtx_input: np.ndarray, resolution: Tuple[int,int]) -> np.ndarray:
        
        # handle empty matrix: if sum of input matrix is zero, the decay matrix is entirely zero
        if np.sum(mtx_input) == 0:
            return np.zeros(mtx_input.shape, np.float32) 

        mtx_proximity = Spatial.get_distance_raster(mtx_input, resolution)
        out_of_reach_mask = mtx_proximity > (threat.get_maximum_distance() * 1000)
        
        # now scale this matrix according to threat properties
        if threat.get_decay_function() == ThreatDecayFunction.Linear:
            # irxy = 1 - (dxy/dmax)
            np.divide(mtx_proximity, (threat.get_maximum_distance() * 1000), out=mtx_proximity)
            np.subtract(1, mtx_proximity, out=mtx_proximity)

        elif threat.get_decay_function() == ThreatDecayFunction.Exponential:
            # irxy = exp(-1*(2.99/dmax)*dxy)
            const_factor = -1 * (2.99/(threat.get_maximum_distance() * 1000))
            np.multiply(const_factor, mtx_proximity, out=mtx_proximity)
            np.exp(mtx_proximity, out=mtx_proximity)

        # mask negative values
        mtx_proximity[out_of_reach_mask] = 0
        return mtx_proximity

    @staticmethod
    def vector_to_raster(shape_geodf: gpd.GeoDataFrame, out_shape: Tuple[int,int], out_transform, fill_attribute: str = None, background_fill=0, default_value=1, all_touched: bool = False, dtype = np.float32) -> np.ndarray:
        
        # rasterize shapefiles passed as geodf
        # this will construct the required iterable for rasterization depending on whether a shape attrib is used for values; otherwise default value will be used
        if fill_attribute is not None:
            shapes = [(geom,value) for geom, value in zip(shape_geodf.geometry, shape_geodf[fill_attribute])]
        else: 
            shapes = [geom for geom in shape_geodf.geometry] 
        
        if len(shapes) > 0:
            mtx_rasterized_shapes = rasterio.features.rasterize(shapes=shapes, fill=background_fill, default_value=default_value, all_touched=all_touched, out_shape=out_shape, dtype=dtype, transform=out_transform)
        else:
            mtx_rasterized_shapes = np.full(shape=out_shape, fill_value=background_fill, dtype=dtype)
        
        return mtx_rasterized_shapes
    
    @staticmethod
    def make_road_layer(geodf_roads: gpd.GeoDataFrame, road_type: str, out_shape: Tuple[int,int], out_transform) -> np.ndarray:

        # make raster of roads to be re-used across scenarios
        roads_of_type = geodf_roads[(geodf_roads['highway'] == road_type)].copy()  
        mtx_roads = Spatial.vector_to_raster(roads_of_type, out_shape, out_transform, fill_attribute=None)
        return mtx_roads  
    
    @staticmethod
    def get_total_area(mtx_land_cover: np.ndarray, nodata_values: List[int], resolution: Tuple[float, float]) -> float:

        # to determine extent, we cannot simply use shape of matrix, as this may include nodata values.
        # therefore, determine the number of actual data pixels
        # then, multiply pixel count with specified resolution in x and y dim

        # determine nodata pixels by masking nodata values
        nodata_pixels = np.isin(mtx_land_cover, nodata_values, invert=False)
        # make masked array for further processing
        mtx_masked = ma.masked_array(mtx_land_cover, mask=nodata_pixels)
        return mtx_masked.count() * (resolution[0] * resolution[1])
    
    @staticmethod
    def mask_habitat_patches(mtx_land_cover: np.ndarray, habitat_values: List[LandCoverClassVegetationProperties]) -> np.ndarray:

        habitats = mtx_land_cover.copy()
        
        # get habitat class values        
        class_values = [habitat_type.get_land_cover_class() for habitat_type in habitat_values] 
            
        # mask pixels that are habitat, and return masked habitats as Bool layer
        habitat_mask = np.isin(habitats, class_values, invert=False)    
        habitats[~habitat_mask] = 0
        habitats[habitat_mask] = 1

        return habitats

    @staticmethod
    def _compute_patches(mtx_input: np.ndarray, contiguity: PatchContiguity) -> Tuple[int, np.ndarray]:
        
        # now determine patches. to do so:
        # define connectivity of patches
        # first, instantiate 3x3 mtx filled with 1.
        # this is equal to queens contiguity, and we derive rook if needed
        patch_connectivity = np.full((3,3), 1)
        if contiguity is PatchContiguity.Rook:
            # set corners to 0 for rook
            patch_connectivity[0][0] = 0
            patch_connectivity[0][2] = 0
            patch_connectivity[2][0] = 0
            patch_connectivity[2][2] = 0

        # then make patch layer from habitat layer
        mtx_patch_features = np.zeros(mtx_input.shape, dtype=np.int32)
        patch_count = ndimage.label(mtx_input, structure=patch_connectivity, output=mtx_patch_features)   

        return patch_count, mtx_patch_features

    @staticmethod
    def determine_patches(
        mtx_land_cover: np.ndarray, 
        habitat_values: List[LandCoverClassVegetationProperties], 
        resolution: Tuple[float,float], 
        contiguity = PatchContiguity.Rook,
        min_pixels: int = None) -> Tuple[int, float, np.ndarray]:
        
        # print(f"Determining patches with min_pixels={min_pixels}")

        def _get_labelled_habitat_mask(mtx_land_cover: np.ndarray, habitat_values: List[LandCoverClassVegetationProperties], resolution):

            # get mask of patches of current class
            mtx_masked_habitat = Spatial.mask_habitat_patches(mtx_land_cover, habitat_values)
        
            # determine scenario-specific habitat area atotal
            habitat_total_area = Spatial.get_total_area(mtx_masked_habitat, nodata_values=[0], resolution=resolution)
            patch_count, mtx_patch_features = Spatial._compute_patches(mtx_masked_habitat, contiguity=contiguity)

            return mtx_masked_habitat, habitat_total_area, patch_count, mtx_patch_features


        # get habitat mask in a first pass, unfiltered
        mtx_masked_habitat, habitat_total_area, patch_count, mtx_patch_features = _get_labelled_habitat_mask(mtx_land_cover=mtx_land_cover, habitat_values=habitat_values, resolution=resolution)
        # print(f"Identified {patch_count} in first pass with total area of {habitat_total_area} m²")

        # at this point, we have the following elements
        # mtx_masked_habitat: a binary mask of habitat pixels
        # habitat_total_area: the total area of habitat in m²
        # mtx_patch_features: a labelled matrix of habitat patches, where each patch has a unique integer label


        # for debug purposes
        # ds_test = rioxarray.open_rasterio('c:/users/sebsc/Documents/urban_nbsos_data/India/Chennai/eo_indicators/lppi.tif')
        # ds_out = ds_test.copy()
        # ds_out.values[0] = mtx_masked_habitat
        # ds_out.rio.to_raster('c:/users/sebsc/Documents/urban_nbsos_data/India/Chennai/test_masked_habitat.tif')
        # ds_out.values[0] = mtx_patch_features
        # ds_out.rio.to_raster('c:/users/sebsc/Documents/urban_nbsos_data/India/Chennai/test_patch_features.tif')


        # for each identified patch, determine properties including the number of pixels and area in m²
        patch_property_table = regionprops_table(label_image=mtx_patch_features, properties=('label', 'num_pixels', 'perimeter'))
        patch_property_table['area'] = patch_property_table['num_pixels'] * (resolution[0] * resolution[1])
        patch_property_table['perimeter'] = patch_property_table['perimeter'] * resolution[0]        

        # at this point, patch_property_table, respectively, patch_df contains the properties of all identified patches, including their area in m² and perimeter in m.
        # this is for all identified patches, irrespective of patch pixel count



        # if user requested filtering-down of habitat mask, in the following, we need to filter as follows:
        # mtx_patch_features: overwrite with 0 for all pixels that are part of patches that do not meet the minimum pixel count
        # mtx_masked_habitat: overwrite with 0 for all pixels that are part of patches that do not meet the minimum pixel count
        


        # filter region props table to only include patches with more than X pixels unless min_pixels is None
        if min_pixels is not None:

            # create a boolean mask to filter labelled habitat patches by pixel size 
            keep_mask = patch_property_table['num_pixels'] >= min_pixels

            # reconstruct dictionary
            patch_property_table = {
                key: values[keep_mask]
                for key, values in patch_property_table.items()
            }

  
        # map labelled features from input matrix to computed area 
        mtx_patch_area = np.zeros(mtx_patch_features.shape, dtype=np.float32)
        map_array(
            mtx_patch_features,
            patch_property_table['label'],
            patch_property_table['area'],
            out=mtx_patch_area
        )

        # overwrite the patch area matrix with 0 for all pixels that are not part of a patch
        mtx_patch_features[mtx_patch_area == 0] = 0
        mtx_masked_habitat[mtx_patch_area == 0] = 0

        # update patch count
        patch_count = len(patch_property_table['label'])
        habitat_total_area = Spatial.get_total_area(mtx_masked_habitat, nodata_values=[0], resolution=resolution)
        
        # print(f"Filtered {patch_count} patches with total area of {habitat_total_area} m²")
        

        # for debug purposes
        # ds_out.values[0] = mtx_masked_habitat
        # ds_out.rio.to_raster('c:/users/sebsc/Documents/urban_nbsos_data/India/Chennai/text_mtx_masked_habitat_repl.tif')
        # ds_out.values[0] = mtx_patch_area
        # ds_out.rio.to_raster('c:/users/sebsc/Documents/urban_nbsos_data/India/Chennai/test_patch_area_repl.tif')
        # ds_out.values[0] = mtx_patch_features
        # ds_out.rio.to_raster('c:/users/sebsc/Documents/urban_nbsos_data/India/Chennai/test_patch_features_repl.tif')

        patch_df = pd.DataFrame({
            'label' : patch_property_table['label'],
            'pixels' : patch_property_table['num_pixels'],
            'area' : patch_property_table['area'],
            'perimeter' : patch_property_table['perimeter']
        })

        # for debug purposes
        # patch_df.to_excel('c:/users/sebsc/Documents/urban_nbsos_data/India/Chennai/test_patch_df_excl.xlsx', index=False)
        
        return patch_count, habitat_total_area, mtx_patch_features, patch_df, mtx_patch_area

    @staticmethod
    def edge_dilation(mtx_patches: np.ndarray, mtx_roads: np.ndarray, iterations: int) -> Tuple[np.ndarray, int, np.ndarray]:
        
        # this should return a matrix with binary dilation applied to patch structures so that they are buffered by specified number of iterations 
        mtx_dilated_edges = ndimage.binary_dilation(mtx_patches, structure=None, iterations=iterations).astype(np.int32)

        # dilated edges result in a 0/1 matrix, where 1 indicates areas within connectivity threshold
        # burn in roads to model anthropogenic impact
        mtx_dilated_edges[((mtx_roads == 1) & (mtx_dilated_edges > 0))] = 0
        components_count, mtx_component_labels = Spatial._compute_patches(mtx_dilated_edges, contiguity=PatchContiguity.Rook)
         
        return mtx_dilated_edges, components_count, mtx_component_labels

    @staticmethod
    def clip_to_area_of_interest(rst_land_cover: rasterio.DatasetReader, shape_aoi, ) -> np.ndarray:
        """Clip an input raster to the extent of the area of interest polygon shapefile.

        :param rst_land_cover: The input raster to clip.
        :type rst_land_cover: rasterio.DatasetReader
        :param shape_aoi: The shapefile containing the area of interest as polygon feature.
        :type shape_aoi: str
        :return: _description_
        :rtype: np.ndarray
        """
        cropped_land_use_matrix = rasterio.mask.mask(rst_land_cover, shape_aoi, indexes=1)
        return cropped_land_use_matrix    

    
    @staticmethod
    def _get_patch_properties(patch_idx: int, patch_slices: List[any], mtx_patch_features: np.ndarray, resolution: Tuple[float,float]) -> Tuple[int, int, float, np.ndarray]:
        
        obj_slice = patch_slices[patch_idx]
        obj_label = patch_idx + 1

        # get the current habitat patch as slice 
        current_habitat_patch_slice = mtx_patch_features[obj_slice].copy()
        
        # cut out the current patch from the habitat matrix using the slice, then
        # properly mask out pixels of neighbouring patches that may be overlapping the slice
        current_patch_mask = np.isin(current_habitat_patch_slice, [obj_label], invert=False)
        current_habitat_patch_slice[current_patch_mask] = 1
        current_habitat_patch_slice[~current_patch_mask] = 0       

        # now that we have zeroed all non-clump pixels, determine patch area
        # the area in sqm is equal to the number of patch pixels * (resolution x by resolution y)
        current_pixel_count = np.sum(current_habitat_patch_slice)
        current_patch_size = (current_pixel_count * (resolution[0] * resolution[1])) 
        mtx_current_patch_size = np.multiply(current_habitat_patch_slice, current_patch_size)
        
        return patch_idx, current_pixel_count, current_patch_size, mtx_current_patch_size

    @staticmethod
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

    @staticmethod
    def match_reference_raster(src_filename: str, template_meta, out_path: str):

        # make proper outfile name
        infile_path = Path(src_filename)
        out_filename = infile_path.name
        outfile_path = os.path.join(out_path, out_filename)

        # source is the raster of origin
        ref_source = rasterio.open(src_filename)        
        
        # template is the destination raster properties
        # template meta is passed to the method and used as reference

        # out is the raster to be written
        out_ref = rasterio.open(outfile_path, 'w', **template_meta)
        rasterio.warp.reproject(
            source=ref_source.read(1),
            destination=rasterio.band(out_ref, 1),
            src_transform=ref_source.transform,
            src_crs=ref_source.crs,
            dst_transform=out_ref.transform,
            dst_crs=out_ref.crs,
            resampling=Resampling.nearest
        )
        out_ref.close()  
        return outfile_path      

    @staticmethod
    def project_raster_datasource(in_file: str, to_crs: int, output_folder: str, out_filename: str = None) -> str:
        """Reprojects a raster file to the recommended EPSG and returns the path of the reprojected file. 

        :param in_file: Raster to reproject.
        :type in_file: str
        :return: Path to the reprojected raster file.
        :rtype: str
        """

        dst_crs = f'EPSG:{to_crs}'

        infile_path = Path(in_file)       
        out_filename = infile_path.name if out_filename is None else out_filename
        outfile_path = os.path.join(output_folder, out_filename)

        # reproject src to dst
        # take care that it matches the reference layer, which should be land cover
        with rasterio.open(in_file) as src:
            
            biodiversity_logger.info(f"Reprojecting {in_file} from {src.crs} to {dst_crs}")
            transform, width, height = calculate_default_transform(
                src.crs, dst_crs, src.width, src.height, *src.bounds)
            kwargs = src.meta.copy()
            kwargs.update({
                'crs': dst_crs,
                'transform': transform,
                'width': width,
                'height': height
            })

            with rasterio.open(outfile_path, 'w', **kwargs) as dst:
                for i in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(src, i),
                        destination=rasterio.band(dst, i),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=transform,
                        dst_crs=dst_crs,
                        resampling=Resampling.nearest)

        return outfile_path