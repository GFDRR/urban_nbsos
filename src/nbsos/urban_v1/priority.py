from .nbsos_exports import ScanBenefit, Benefits, NbsTypes
from .scan_context import ScanContext
from .gis_tools import export_raster, rasterize, rasterize_df
from .biodiversity_module.biodiversity_exports import PhysiognomicLandscapeType
from .gis_tools import rasterize_geom       
from typing import List, Tuple, Dict
import rasterstats
import geopandas as gpd
import pandas as pd
from geojson import FeatureCollection
import os
import exactextract as ee
import numpy as np
import rasterio
import rioxarray
from pathlib import Path
import jenkspy
from tqdm.auto import tqdm as tqauto
import random
import shutil

import logging
logger = logging.getLogger('nbsos.scan')


def get_benefit_weight_combinations(focal_benefit: Benefits, list_of_benefits: List[Benefits], focal_weight: float = 0.8) -> Dict[Benefits, float]:
    
    # determine number of benefits other than focal benefit
    number_of_remainer_benefits = len(list_of_benefits) - 1 # this is the number of benefits to which we need to allocate 20% of weight equally; 1 benefit will be 80% 
    # determine weights for benefits other than focal benefit
    residual_total_weight = 1.0 - focal_weight
    if number_of_remainer_benefits > 0:
        residual_benefit_weight = residual_total_weight / number_of_remainer_benefits
    else:
        residual_benefit_weight = 0.0
    
    weights = []
    for b in list_of_benefits:
        if b is focal_benefit:
            weights.append(focal_weight)
        else:
            weights.append(residual_benefit_weight)
    return dict(zip(list_of_benefits, weights))


def get_weighted_raster(
        context: ScanContext,
        raster_path: str, 
        current_mtx, 
        weight: float,
        weighting_schema_name: str,
        normalize_values: bool = True,
        norm_min_value: float = None,
        norm_max_value: float = None,        
    ) -> np.ndarray:
    
    # assert min and max values for normalization are provided if normalization is requested
    if normalize_values:
        if norm_min_value is None or norm_max_value is None:
            raise ValueError('Normalization requested but min and max values for normalization not provided. Please provide norm_min_value and norm_max_value for normalization.')

    
    # if no normalization is requested, simply return the weighted raster combination.
    # if weight is one, this is equal to getting the original raster    
    if not normalize_values:
        return current_mtx * weight

    local_max_val = np.nanmax(current_mtx)
    local_min_val = np.nanmin(current_mtx)
    local_range = local_max_val - local_min_val

    # assert that the value range covers ALL raw benefit values, to normalize nbs performance relative to each other 
    abs_max_val = norm_max_value
    abs_min_val = norm_min_value
    abs_range =  abs_max_val - abs_min_val

    result = np.zeros(current_mtx.shape)  # Initialize result array with zeros
    if abs_range == 0:
        return result
    else:
        # apply the following function where current_mtx is not 0 or nan
        np.divide((current_mtx - abs_min_val), abs_range, where=(current_mtx != 0), out=result)

        # multiply by weight
        result = result * weight
        return result





def get_raster_paths_for_nbs_and_benefit(context: ScanContext, nbs: NbsTypes, current_benefit: ScanBenefit) -> Tuple[str, str, str]:
    """Returns nbs creation mask, path to raw-value benefit raster, and raw storage benefits for nbs with storage, if applicable.
    :param context: ScanContext object containing paths and parameters for the scan.
    :type context: ScanContext
    :param nbs: NbsType for which to get the raster paths.
    :type nbs: NbsTypes
    :param current_benefit: ScanBenefit object containing the benefit type and parameters for which to get the raster paths.
    :type current_benefit: ScanBenefit
    """
    raster_path = None
    storage_path = None
    nbs_creation_path = None
    
    nbs_creation_path = os.path.join(context.local_root, context.country, context.city, 'suitability', f'{nbs.get_name()}_creation.tif')
    current_benefit_type = current_benefit.benefit

    
    if current_benefit_type is Benefits.UrbanHeat:
        raster_path = os.path.join(context.local_root, context.country, context.city, 'Heat', 'outputs', 'heat_benefits', f'{nbs.get_name()}_heat_benefits.tif')

    elif current_benefit_type is Benefits.Recreation:
        raster_path = os.path.join(context.local_root, context.country, context.city, 'rec_benefit', f'recreation_benefit_{nbs.get_name()}.tif')            
    
    elif current_benefit_type is Benefits.AvoidedErosion:
        raster_path = os.path.join(context.local_root, context.country, context.city, 'Erosion', 'output', 'erosion_benefits', f'diff_{nbs.get_name()}_erosion.tif')            
    
    elif current_benefit_type is Benefits.SocialCohesion:
        raster_path = os.path.join(context.local_root, context.country, context.city, 'rec_benefit', f'social_cohesion_benefit_{nbs.get_name()}.tif')            
    
    elif current_benefit_type is Benefits.PluvialFlooding:
        raster_path = os.path.join(context.local_root, context.country, context.city, 'flood_benefit', f'flood_benefit_{nbs.get_name()}.tif')            
        
        # TODO : REVIEW NBS FAMILIES WITH STORAGE
        nbs_with_storage = [NbsTypes.BioretentionAreas, NbsTypes.InlandWetlands, NbsTypes.OpenGreenSpace]
        if nbs in nbs_with_storage:
            storage_path = os.path.join(context.local_root, context.country, context.city, 'flood_benefit', f'flood_benefit_{nbs.get_name()}_storage.tif')            
        
    elif current_benefit_type is Benefits.Biodiversity:        
        # get the raster, depending on landscape type to assess as defined in scan params, or Any by default if no        
        if current_benefit.params is not None:
            current_landscape_type = current_benefit.params.get('landscape_types_to_assess', PhysiognomicLandscapeType.Any)
        else:
            current_landscape_type = PhysiognomicLandscapeType.Any

        current_landscape_type_name = current_landscape_type.name
        raster_path = os.path.join(context.local_root, context.country, context.city, 'Biodiversity', 'outputs', 'connectivity', current_landscape_type_name, f'diffs_{nbs.get_name()}_creation_in_{current_landscape_type_name}.tif')
           
    elif current_benefit_type is Benefits.JobsCreation:
        raster_path = os.path.join(context.local_root, context.country, context.city, 'JobsCreation', f'{nbs.get_name()}_jobs_benefits.tif')

    return nbs_creation_path, raster_path, storage_path


def get_benefit_raster_for_nbs(
        context: ScanContext, 
        nbs: NbsTypes,
        nbs_list: List[NbsTypes], 
        benefit: ScanBenefit, 
        weight: float = 1.0,
        weighting_schema_name: str = None,
        normalize_values: bool = True,
        norm_min_value: float = None,
        norm_max_value: float = None,        
    ):

    if normalize_values:
        if norm_min_value is None or norm_max_value is None:
            raise ValueError('Normalization requested but min and max values for normalization not provided. Please provide norm_min_value and norm_max_value for normalization.')

    weighted_mtx = None
    nbs_creation_path, raster_path, storage_path = get_raster_paths_for_nbs_and_benefit(context, nbs, benefit)

    # intersect with creation opportunity mask
    nbs_creation_mask = rioxarray.open_rasterio(nbs_creation_path).values[0]    

    if raster_path is not None:
        # if valid path returned, open the creation benefits raster from raster_path
        benefit_value_mtx = rioxarray.open_rasterio(raster_path).values[0]
    else:
        # if None returned (this seems to potentially happen with carbon benefits), 
        # create a matrix of zeros, indicating no benefits for the nbs family
        # we use the creation mask for this purpose, which should always be valid
        benefit_value_mtx = np.zeros(nbs_creation_mask.shape, dtype=np.float32)


    # make intersects before the weighting
    benefit_value_mtx = nbs_creation_mask * benefit_value_mtx 

    weighted_mtx = get_weighted_raster(
        context, 
        raster_path,
        benefit_value_mtx, 
        weight, 
        weighting_schema_name, 
        normalize_values=normalize_values,
        norm_min_value=norm_min_value,
        norm_max_value=norm_max_value,
    )
    
    if storage_path is not None:

        storage_mtx = rioxarray.open_rasterio(storage_path).values[0]
        storage_mtx = nbs_creation_mask * storage_mtx

        weighted_mtx_storage = get_weighted_raster(
            context, 
            storage_path,
            storage_mtx, 
            weight, 
            weighting_schema_name, 
            normalize_values=normalize_values,
            norm_min_value=norm_min_value,
            norm_max_value=norm_max_value,
        )
        weighted_mtx = np.add(weighted_mtx, weighted_mtx_storage)


    return weighted_mtx

def get_weighted_benefits_for_nbs(
        context: ScanContext, 
        nbs: NbsTypes, 
        nbs_list: List[NbsTypes],
        benefits_list: List[ScanBenefit], 
        nbs_benefit_value_range: Dict[Benefits, Dict[str, float]],
        weights: Dict[Benefits, float], 
        weighting_schema_name: str,        
    ) -> np.ndarray:	
    
    benefits_raster = None
    for benefit in benefits_list:

        # obtain the absolute minimum and maximum values provided across all nbs, for the current benefit
        current_value_range = nbs_benefit_value_range[benefit.benefit]
        min_val = current_value_range['absolute_min']
        max_val = current_value_range['absolute_max']

        weighted_mtx = get_benefit_raster_for_nbs(
            context, 
            nbs, 
            nbs_list, 
            benefit, 
            weights[benefit.benefit], 
            weighting_schema_name, 
            normalize_values=True,
            norm_min_value=min_val,
            norm_max_value=max_val,
        )
        benefits_raster = np.add(benefits_raster, weighted_mtx) if benefits_raster is not None else weighted_mtx        
    
    return benefits_raster

def get_nbs_mca(
        context: ScanContext, 
        nbs_list: List[NbsTypes], 
        benefits_list: List[ScanBenefit],
        percentiles_for_binning: List[int] = [0, 25, 75, 100],
        heavy_focal_weight: float = 0.8,
        assess_equal_weights_scenario: bool = True,
        assess_benefit_heavy_scenario: bool = True,
        assess_benefit_only_scenario: bool = True,
        custom_scenarios: Dict[str, Dict[Benefits, float]] = None,
        highest_benefits_percentiles: List[int] = [80],
        focus_area_paths: List[str] = None,
        right: bool = True,
    ):
    
    logger.info('Performing multi-criteria analysis')

    # open hex raster to get metadata and mask AOI by zero values in the hex raster
    template_raster = os.path.join(context.local_root, context.country, context.city, 'ancillary', f'{context.city}_aoi_hex.tif')
    with rasterio.open(template_raster) as rst:
        out_meta = rst.meta.copy()
        out_shape = rst.shape
        
    out_meta.update({
        'nodata' : np.nan,
        'dtype' : np.float32
    })

    nbs_benefit_value_range = {}

    # determine percentiles based on maximum benefits across all nbs, per each benefit
    # bin each nbs-benefit raster using these bins, and write output rasters
    for benefit in tqauto(benefits_list, 'Binning nbs-benefit rasters based on overall percentiles'):

        out_path = os.path.join(context.local_root, context.country, context.city, 'output', 'nbs_compare', benefit.benefit.get_name())
        os.makedirs(out_path, exist_ok=True)

        # get the raster stack for the current benefit, that includes all nbs
        unweighted_raw_nbs_benefit_rasters = get_raw_benefit_raster_stack(context, nbs_list, benefit)
        
        # combine all unweighted creation benefit rasters for all nbs (of current benefit) into a single stack
        unweighted_stacked_rasters = np.stack(list(unweighted_raw_nbs_benefit_rasters.values()))
        
        # determine minimum and maximum benefit value, per pixel, for current benefit, across all nbs
        # determine min and max values
        mtx_max_value = np.full(out_shape, np.nan, dtype=np.float32)
        np.nanmax(unweighted_stacked_rasters, axis=0, out=mtx_max_value)
        mtx_min_value = np.full(out_shape, np.nan, dtype=np.float32)
        np.nanmin(unweighted_stacked_rasters, axis=0, out=mtx_min_value)

        # write out for possible later use
        out_raster = os.path.join(out_path, f'{benefit.benefit.get_name()}_max_vals.tif')
        export_raster(mtx_max_value, template_raster, out_raster, out_meta)
        out_raster = os.path.join(out_path, f'{benefit.benefit.get_name()}_min_vals.tif')
        export_raster(mtx_min_value, template_raster, out_raster, out_meta)
        
        # determine min and max absolute values across all nbs for current benefit
        max_val = np.nanmax(mtx_max_value)
        min_val = np.nanmin(mtx_min_value)

        # store these as reference for later
        nbs_benefit_value_range[benefit.benefit] = {
            'absolute_max' : max_val,
            'absolute_min' : min_val
        }

        # create bins based on max raster, but check if min value is smaller than first bin and if so, replace first bin with min value to avoid losing information in binning
        bins = get_bins(mtx_max_value, 'percentiles', percentiles_for_binning)       
        bins_replaced = bins.copy()         
        if min_val < bins[0]:
            bins_replaced[0] = min_val  

        # iterate over nbs and bin their benefit rasters according to overall bins
        for nbs in nbs_list:
            # nbs performance matrix
            raw_mtx = unweighted_raw_nbs_benefit_rasters[nbs.get_name()]
            # get mask of 0 values and nan for raw_mtx
            current_mask = np.isin(raw_mtx, [0], invert=False) | np.isnan(raw_mtx) 

            binned_mtx = apply_bins(raw_mtx, bins_replaced, right=True)
            if binned_mtx is not None:
                # write binned raster
                binned_mtx[current_mask] = np.nan
                
                out_raster = os.path.join(out_path, f'{nbs.get_name()}_{benefit.benefit.get_name()}_binned.tif')
                export_raster(binned_mtx, template_raster, out_raster, out_meta)



    # apply the following weighting schemes
    # equal weights
    list_of_benefits = [f.benefit for f in benefits_list]
    if assess_equal_weights_scenario:
        
        equal_weights = [1/len(list_of_benefits) for f in benefits_list]
        equal_weights = dict(zip(list_of_benefits, equal_weights))

        # conduct equal weighting
        out_path = os.path.join(context.local_root, context.country, context.city, 'output', 'mca', 'equal')
        write_weighted_combinations(
            context, 
            nbs_list, 
            benefits_list,
            nbs_benefit_value_range, 
            equal_weights, 
            'equal_weights', 
            out_path, 
            percentiles_for_binning=percentiles_for_binning,
            highest_benefits_percentiles_for_binning=highest_benefits_percentiles,
            focus_area_paths=focus_area_paths,
            right=right,
        )

    # assess additional weighting schemas
    # this only works for > 1 benefits    
    if assess_benefit_heavy_scenario:
        if len(list_of_benefits) > 1:
            
            # benefit-heavy (80% vs. other benefits)
            benefit_heavy_schemes = {}
            for b in list_of_benefits:
                # when we iterate over this list, each benefit b should be heavy once
                heavy_weight_set = get_benefit_weight_combinations(b, list_of_benefits, focal_weight=heavy_focal_weight)
                benefit_heavy_schemes[b] = heavy_weight_set

            # iterate over focal benefits in these schemas and conduct the weighting
            out_path = os.path.join(context.local_root, context.country, context.city, 'output', 'mca', 'benefit_heavy')
            for focal_benefit, weight_schema in tqauto(benefit_heavy_schemes.items(), 'Assessing benefit-heavy weighting schema'):
                # generate a schema name
                schema_name = f'{focal_benefit.value}_heavy'
                # get weighted combinations
                write_weighted_combinations(
                    context, 
                    nbs_list, 
                    benefits_list, 
                    nbs_benefit_value_range,
                    weight_schema, 
                    schema_name, 
                    out_path, 
                    percentiles_for_binning=percentiles_for_binning,
                    highest_benefits_percentiles_for_binning=highest_benefits_percentiles,
                    focus_area_paths=focus_area_paths,
                    right=right,
                )

        else:
            logger.warning('Skip benefit-heavy weighting scenarios: Only one benefit in model')

    if assess_benefit_only_scenario:
        # benefit-only (100%, neglecting other benefits) - this should also be done for single-benefit models
        benefit_only_schemas = {}
        for b in list_of_benefits:
            benefit_only_weight_set = get_benefit_weight_combinations(b, list_of_benefits, focal_weight=1.0)
            benefit_only_schemas[b] = benefit_only_weight_set

        # iterate over focal benefits in these schemas and conduct the weighting
        out_path = os.path.join(context.local_root, context.country, context.city, 'output', 'mca', 'benefit_only')
        for focal_benefit, weight_schema in tqauto(benefit_only_schemas.items(), 'Assessing benefit-only weighting schema'):
            # generate a schema name
            schema_name = f'{focal_benefit.value}_only'
            # get weighted combinations
            write_weighted_combinations(
                context, 
                nbs_list, 
                benefits_list,
                nbs_benefit_value_range, 
                weight_schema, 
                schema_name, 
                out_path, 
                percentiles_for_binning=percentiles_for_binning,
                highest_benefits_percentiles_for_binning=highest_benefits_percentiles,
                focus_area_paths=focus_area_paths,
                right=right,
            )

    if custom_scenarios is not None:
        for schema_name, weight_schema in tqauto(custom_scenarios.items(), 'Assessing custom weighting schemas'):

            # make weighted combination according to custom scenario
            out_path = os.path.join(context.local_root, context.country, context.city, 'output', 'mca', schema_name)

            required_paths = [
                out_path,
                os.path.join(out_path, 'nbs_weighted')
            ]
            # assert outpath exists
            for p in required_paths:
                os.makedirs(p, exist_ok=True)
            
            # get weighted combinations
            write_weighted_combinations(
                context, 
                nbs_list, 
                benefits_list,
                nbs_benefit_value_range, 
                weight_schema, 
                f'{schema_name}_weights', 
                out_path, 
                percentiles_for_binning=percentiles_for_binning,
                highest_benefits_percentiles_for_binning=highest_benefits_percentiles,
                focus_area_paths=focus_area_paths,
                right=right,
            )


def apply_bins(mtx: np.ndarray, bins: List[float], right: bool = True):
    
    mtx_binned = None
    try:
        # bin according to percentiles, right-inclusive                
        mtx_binned = np.digitize(mtx, bins, right=right)
        mtx_binned = mtx_binned.astype(np.float32)
        
        # conduct some masking on the percentiles raster.
        # Mask 0 values if needed. This is the case if the raster includes values lower than first bin's left side, but 0 is not included in the bins due to excluding the value in percentiles computation.
        mtx_binned[mtx_binned == 0] = np.nan
        # Mask values with len(bins). These are values outside the max value. That could be NaN values.
        mtx_binned[mtx_binned == len(bins)] = np.nan

    except:
        mtx_binned = None
    return mtx_binned

def get_bins(mtx: np.ndarray, binning_method: str, percentiles: List[int] = None) -> List[float]:

    # return bins, created from a raster
    if binning_method == 'jenks':
        mtx_all_vals = np.nan_to_num(mtx).flatten() 
        sampled_indices = random.sample(range(len(mtx_all_vals)), 10000) 
        sampled_values = mtx_all_vals[sampled_indices].tolist()

        # assert min and max inclusive sample
        sampled_values.append(np.min(mtx_all_vals))
        sampled_values.append(np.max(mtx_all_vals))

        bins = jenkspy.jenks_breaks(sampled_values, n_classes=5)

    elif binning_method == 'percentiles':
        if percentiles is None:
            raise ValueError('Percentiles for binning not provided. Please provide a list of percentiles.')        
        
        # this method specifically creates percentiles for values > 0,
        # as zero is considered background or not of interest
        tmp = mtx[mtx > 0].flatten()
        bins = np.nanpercentile(tmp, percentiles)
        
    return bins

def bin_raster(mtx: np.ndarray, binning_method: str, percentiles: List[int] = None, right: bool = True) -> np.ndarray:
   
    # result matrix
    mtx_binned = None

    # get bins
    bins = get_bins(mtx, binning_method, percentiles)

    # try to apply bins to raster 
    mtx_binned = apply_bins(mtx, bins, right=right)
    return mtx_binned

def get_raw_benefit_raster_stack(
        scan_context: ScanContext,
        nbs_list: List[NbsTypes],
        benefit: ScanBenefit,
    ) -> Dict[str, np.ndarray]:
    
    # store rasters for stacking
    raw_nbs_benefit_value_rasters = {}

    # iterate over nbs families
    for nbs in nbs_list:
        
        # simply return unweighted and not normalized raster for current nbs and benefit
        unweighted_nbs_benefit_raster = get_benefit_raster_for_nbs(
            scan_context, 
            nbs, 
            nbs_list, 
            benefit, 
            1.0, 
            None, 
            normalize_values=False, 
            norm_min_value=None, 
            norm_max_value=None,
        )
        
        # replace 0 values with nan as zero in the benefits is in all cases no impact == similar performance than baseline
        # replacing by NaN will also yield negative values as local maxima
        unweighted_nbs_benefit_raster[unweighted_nbs_benefit_raster == 0] = np.nan
        raw_nbs_benefit_value_rasters[nbs.get_name()] = unweighted_nbs_benefit_raster

    return raw_nbs_benefit_value_rasters


def get_weighted_raster_stack(
        context: ScanContext, 
        nbs_list: List[NbsTypes], 
        benefits_list: List[ScanBenefit], 
        nbs_benefit_value_range: Dict[Benefits, Dict[str, float]],
        weighting_schema: Dict[Benefits, float], 
        weighting_schema_name: str, 
        out_path: str,
        percentiles_for_binning: List[int] = [0, 25, 75, 100],
    ):

    
    # use hex raster as template raster
    template_raster = os.path.join(context.local_root, context.country, context.city, 'ancillary', f'{context.city}_aoi_hex.tif')
    with rasterio.open(template_raster) as rst:
        out_meta = rst.meta.copy()
        out_shape = rst.shape
        suitability_mask = rst.read(1)
        cluster_mask = rst.read(1)

    cluster_mask[cluster_mask > 0] = 1
    suitability_mask[suitability_mask > 0] = 1

    out_meta.update({
        'nodata' : np.nan,
        'dtype' : np.float32
    })
    
    # cluster_raster = os.path.join(context.local_root, context.country, context.city, 'cluster_maps', 'clusters.tif')
    # with rasterio.open(cluster_raster) as rst:
    #     cluster_mask = rst.read(1)


    weighted_rasters = [
        np.full(out_shape, fill_value=-9999, dtype=np.float32)
    ]

    for nbs in nbs_list:
        
        # get weighted benefits raster for current nbs
        mtx_weighted_benefits = get_weighted_benefits_for_nbs(
            context, 
            nbs, 
            nbs_list, 
            benefits_list, 
            nbs_benefit_value_range,
            weighting_schema, 
            weighting_schema_name,            
        )

        # write weighted benefit raster
        mtx_benefits_masked = mtx_weighted_benefits * suitability_mask        
        mtx_benefits_masked[np.isnan(cluster_mask)] = np.nan
        mtx_benefits_masked[mtx_benefits_masked == 0] = np.nan # mask out 0 values, as these are outside creation areas and not of interest        
        weighted_rasters.append(mtx_benefits_masked)

        out_raster = os.path.join(out_path, 'nbs_weighted', f'{nbs.get_name()}_{weighting_schema_name}.tif')
        export_raster(mtx_benefits_masked, template_raster, out_raster, out_meta)    
        
        write_binned_raster(
            out_raster, 
            template_raster, 
            os.path.join(out_path, 'nbs_weighted', f'{nbs.get_name()}_{weighting_schema_name}_binned.tif'), 
            'percentiles',
            percentiles_for_binning
        )
        # write highest quintile raster
        write_binned_raster(
            out_raster, 
            template_raster, 
            os.path.join(out_path, 'nbs_weighted', f'{nbs.get_name()}_{weighting_schema_name}_highest_quintile.tif'), 
            'percentiles',
            [0, 80, 100]
        )
        

    weighted_stack = np.stack(weighted_rasters)
    return weighted_stack

def write_weighted_combinations(
        scan_context: ScanContext, 
        nbs_list: List[NbsTypes], 
        benefits_list: List[ScanBenefit], 
        nbs_benefit_value_range: Dict[Benefits, Dict[str, float]],
        weighting_schema: Dict[Benefits, float], 
        weighting_schema_name: str, 
        out_path: str, 
        percentiles_for_binning: List[int] = [0, 25, 75, 100],
        highest_benefits_percentiles_for_binning: List[int] = [80],
        focus_area_paths: List[str] = None,
        right: bool = True):

    tqauto.write(f'Assessing {weighting_schema_name} weighting schema')

    # get hex raster
    # use hex path as template as well 
    hex_path = os.path.join(scan_context.local_root, scan_context.country, scan_context.city, 'ancillary', f'{scan_context.city}_aoi_hex.tif')
    with rasterio.open(hex_path) as rst:
        cluster_mask = rst.read(1)
        out_meta = rst.meta.copy()
        out_shape = rst.shape
        suitability_mask = rst.read(1)

    cluster_mask[cluster_mask > 0] = 1
    suitability_mask[suitability_mask > 0] = 1

    out_meta.update({
        'nodata' : np.nan,
        'dtype' : np.float32
    })

    # get weighted raster stack
    weighted_stack = get_weighted_raster_stack(
        scan_context, 
        nbs_list, 
        benefits_list, 
        nbs_benefit_value_range,
        weighting_schema, 
        weighting_schema_name, 
        out_path,
        percentiles_for_binning=percentiles_for_binning,
    )

    # determine max benefit per pixel across stack
    mtx_max_value = np.zeros(out_shape, dtype=np.float32)
    np.nanmax(weighted_stack, axis=0, out=mtx_max_value)
    
    # max value raster returns -9999 in nodata locations: mask these to nan
    mtx_max_value[mtx_max_value == -9999] = np.nan

    # multiply by suitability mask to automatically mask out non-creation areas
    mtx_max_value = mtx_max_value * suitability_mask
    mtx_max_value[np.isnan(cluster_mask)] = np.nan
    
    # export maximum benefit raster for manual inspection and use in binning
    max_benefits_raster_path = os.path.join(out_path, f'{weighting_schema_name}_max_benefit.tif')
    export_raster(mtx_max_value, hex_path, max_benefits_raster_path, out_meta)
    
    write_binned_raster(
        max_benefits_raster_path, 
        hex_path, 
        os.path.join(out_path, f'{weighting_schema_name}_max_benefit_binned.tif'), 
        'percentiles',
        percentiles_for_binning,
        right=right
    )

    # identify optimal nbs 
    mtx_optimal_nbs = np.zeros(out_shape, dtype=np.int64)
    # identify index of highest benefit value across stack to identify optimal nbs
    np.nanargmax(weighted_stack, axis=0, out=mtx_optimal_nbs)
    mtx_optimal_nbs = mtx_optimal_nbs.astype(np.float32)

    # mtx_optimal_nbs will return 0 where -9999 was the max value due to all other layers being nan; 
    # this value should be properly masked here
    mtx_optimal_nbs[mtx_optimal_nbs == 0] = np.nan

    # replace index with the land cover class value of the corresponding nbs
    for idx, nbs in enumerate(nbs_list):
        # get land cover class for nbs
        nbs_land_cover_class = nbs.get_land_cover_class_value()        
        mtx_optimal_nbs[mtx_optimal_nbs == idx+1] = nbs_land_cover_class # +1 due to zero-based index
    
    # mask by max benefit raster
    #mtx_optimal_nbs = mtx_optimal_nbs * suitability_mask
    mtx_optimal_nbs[np.isnan(mtx_max_value)] = np.nan    
    mtx_optimal_nbs[np.isnan(cluster_mask)] = np.nan

    optimum_nbs_raster_path = os.path.join(out_path, f'{weighting_schema_name}_optimum_nbs.tif')
    export_raster(mtx_optimal_nbs, hex_path, optimum_nbs_raster_path, out_meta)

    # identify nbs with highest benefits: 
    # write the nbs with highest benefits (by default, this is 20% or above 80th percentile)
    # iterate over list of highest benefits percentiles provided, which includes the default value if not specified otherwise by user
    for highest_benefit_percentile in highest_benefits_percentiles_for_binning:

        # create label for current highest benefits percentile
        highest_label = f"{100 - highest_benefit_percentile}"

        def get_masked_benefit_matrix(mtx_benefit_max_values: np.ndarray, focus_area_path: str) -> np.ndarray:                

            # read shapefile
            tqauto.write("Masking raster by focal area shapefile...")
            focus_area_gdf = gpd.read_file(focus_area_path)

            # get the target crs from out_meta
            target_crs = out_meta['crs']
            focus_area_gdf = focus_area_gdf.to_crs(target_crs)

            # intersect current raster with mask 
            # rasterize the focal_gdf to the same resolution and extent as the current raster
            # hex path is the template path                                         
            focal_mask = rasterize_geom(focus_area_gdf, hex_path, None)
            result_mtx = mtx_benefit_max_values.copy()
            result_mtx[focal_mask != 1] = np.nan
            return result_mtx, focal_mask

        def write_x_highest_benefits(highest_benefit_percentile: int, mtx_benefit_max_values: np.ndarray, out_path: str, mask: np.ndarray = None):

            #highest_benefits_raster_path = out_path + "highest_benefits.tif"
            mtx_max_benefit_quintiles = write_binned_raster(
                input_raster_path=None,
                template_raster_path=hex_path,
                out_path=None, # highest_benefits_raster_path,
                binning_method='percentiles',
                percentiles=[0, highest_benefit_percentile, 100],
                input_mtx=mtx_benefit_max_values,
                input_meta=out_meta
            )

            if mtx_max_benefit_quintiles is not None:                
                # make new matrix to overwrite values
                mtx_optimal_nbs_highest_benefits = mtx_optimal_nbs.copy()
                # mask out optimal nbs in pixels lower highest benefit percentile
                mtx_optimal_nbs_highest_benefits[mtx_max_benefit_quintiles < 2] = np.nan
                if mask is not None:
                    mtx_optimal_nbs_highest_benefits[mask != 1] = np.nan
                # export raster 
                export_raster(mtx_optimal_nbs_highest_benefits, hex_path, out_path, out_meta) 



        # write highest benefit optimal nbs for the complete aoi
        optimal_nbs_highest_benefits_raster_path = os.path.join(out_path, f'{weighting_schema_name}_optimum_nbs_highest_{highest_label}_benefits.tif')        
        write_x_highest_benefits(highest_benefit_percentile, mtx_max_value, optimal_nbs_highest_benefits_raster_path)



        # for this percentile, iterate over focus areas to determine highest-x-percent benefits within each area, and relative to each area, ignoring the aoi values outside of the focus area
        if focus_area_paths is not None:
            for focus_path in focus_area_paths: 

                # strip file ending and path from focus mask path
                focus_name = os.path.splitext(os.path.basename(focus_path))[0]

                masked_benefit_mtx, focal_mask = get_masked_benefit_matrix(mtx_max_value, focus_path)
                masked_raster_path = os.path.join(out_path, focus_name, f'{weighting_schema_name}_max_benefit_{focus_name}.tif')
                os.makedirs(os.path.dirname(masked_raster_path), exist_ok=True)
                export_raster(masked_benefit_mtx, hex_path, masked_raster_path, out_meta)
    
                optimal_nbs_highest_benefits_focus_area_raster_path = os.path.join(out_path, focus_name, f'{weighting_schema_name}_optimum_nbs_highest_{highest_label}_benefits_{focus_name}.tif')
                os.makedirs(os.path.dirname(optimal_nbs_highest_benefits_focus_area_raster_path), exist_ok=True)
                write_x_highest_benefits(highest_benefit_percentile, masked_benefit_mtx, optimal_nbs_highest_benefits_focus_area_raster_path, mask=focal_mask)
            

        
def write_binned_raster(
        input_raster_path: str = None, 
        template_raster_path: str = None, 
        out_path: str = None, 
        binning_method: str = 'percentiles', 
        percentiles: List[int] = None, 
        right: bool = True,
        input_mtx: np.ndarray = None,
        input_meta: dict = None

    ) -> np.ndarray:
    """Bin raster values based on Jenks or percentiles, and return the binned matrix. If no out_path is specified, the result of the operation is not written to disk.
    The input can either be a path to an input raster, or a np.ndarray and a corresponding metadata dictionary.

    :param input_raster_path: Path to the input raster to be binned
    :type input_raster_path: str
    :param template_raster_path: Path to the template raster for alignment
    :type template_raster_path: str
    :param out_path: Path to save the binned raster
    :type out_path: str
    :param binning_method: Method used for binning ('jenks' or 'percentiles')
    :type binning_method: str
    :param percentiles: Percentiles for binning, required if binning_method is 'percentiles'
    :type percentiles: List[int], optional
    :param right: Whether the bins include the right edge
    :type right: bool, optional
    :return: Binned raster as a numpy array
    :rtype: np.ndarray
    :param input_mtx: Input raster as a numpy array, optional if input_raster_path is provided
    :type input_mtx: np.ndarray, optional
    :param input_meta: Metadata dictionary corresponding to input_mtx, required if input_mtx is provided
    :type input_meta: dict, optional

    """

    if input_mtx is not None and input_meta is not None and input_raster_path is not None:
        raise ValueError("Provide either input_mtx and input_meta, or input_raster_path and template path, not both.")
    

    if input_mtx is not None and input_meta is not None:
        mtx = input_mtx
        meta = input_meta.copy()

    elif input_raster_path is not None:

        with rasterio.open(input_raster_path) as rst:
            mtx = rst.read(1)                
            meta = rst.meta.copy()
    
    # conduct actual binning
    binned_raster = bin_raster(mtx, binning_method, percentiles, right=right)
    
    if binned_raster is not None:

        if template_raster_path is not None and out_path is not None:

            meta.update({
                'nodata' : np.nan,
                'dtype' : np.float32
            })
            
            export_raster(
                binned_raster,
                template_raster_path,
                out_path,
                meta 
            )
    else:
        tqauto.write(f'Failed binning. Skip writing to file.')

    return binned_raster


def get_priority_scores(local_root, country, city, benefit_parameters: List[ScanBenefit], binning_method: str = 'jenks', percentiles: List[int] = [0,20,40,60,80,100]):
    """Determine priority areas for benefits of interest. The method will write priorities as shapefile, raster, and binned raster according to binning method used (jenks or percentiles). 

    :param local_root: Path to local root folder (data repository)
    :type local_root: str
    :param country: Country name of case
    :type country: str
    :param city: City name
    :type city: str
    :param benefit_parameters: List of benefits of interest 
    :type benefit_parameters: List[ScanBenefit]
    :param binning_method: Method used for binning (categorizing) priority scores, either 'jenks' or 'percentiles', defaults to 'jenks'.
    :type binning_method: str, optional
    :param percentiles: Percentiles for binning, defaults to [0,20,40,60,80,100]
    :type percentiles: List[int], optional
    """

    list_of_benefits = [f.benefit for f in benefit_parameters]

    shapes = os.path.join(local_root, country, city, 'ancillary', f'{city}_hex.shp')
    out_shape = os.path.join(local_root, country, city, 'output', 'priority', f'{city}_priorities.shp')

    template_raster = os.path.join(local_root, country, city, 'eo_indicators', 'lppi.tif')

    gdf = gpd.read_file(shapes)
    
    # columns to retain
    columns_to_retain = ['neigh_id', 'geometry']
    gdf = gdf.loc[:, columns_to_retain]
    gdf = gdf.to_crs(4326)

    for benefit in tqauto(list_of_benefits):

        if benefit is Benefits.PluvialFlooding:
            tqauto.write('Assessing flood priority')
            flood_demand = os.path.join(local_root, country, city, 'output', 'flood_risk.tif')  
            flood_raster = rioxarray.open_rasterio(flood_demand)
            gdf = gdf.to_crs(flood_raster.rio.crs)
            flood_stats = ee.exact_extract(flood_raster, gdf, 'mean', output='pandas')
            gdf['flood'] = flood_stats['mean']    
            
            out_path = os.path.join(local_root, country, city, 'output', 'priority', f'flood_priority.tif')
            rasterize_df(
                gdf,
                template_raster,
                out_path,
                'flood'
            )

            write_binned_raster(
                out_path,
                template_raster,
                os.path.join(local_root, country, city, 'output', 'priority', f'flood_priority_binned.tif'),
                binning_method,
                percentiles
            )

            tqauto.write('Assessing flood exposure')            
            in_path = os.path.join(local_root, country, city, 'temp', 'risk_neigh.tif')
            out_path = os.path.join(local_root, country, city, 'output', 'priority', f'flood_exposure.tif')
            shutil.copy(
                src=in_path,
                dst=out_path
            )

            write_binned_raster(
                out_path,
                template_raster,
                os.path.join(local_root, country, city, 'output', 'priority', f'flood_exposure_binned.tif'),
                binning_method,
                percentiles
            )

        if benefit is Benefits.Recreation:
            tqauto.write('Assessing recreation priority')
            recreation_demand = os.path.join(local_root, country, city, 'output', 'recreation_demand.tif')
            recreation_raster = rioxarray.open_rasterio(recreation_demand)
            gdf = gdf.to_crs(recreation_raster.rio.crs)
            recreation_stats = ee.exact_extract(recreation_raster, gdf, 'mean', output='pandas')
            gdf['rec'] = recreation_stats['mean']   

            out_path = os.path.join(local_root, country, city, 'output', 'priority', f'recreation_priority.tif')
            rasterize_df(
                gdf,
                template_raster,
                out_path,
                'rec'
            ) 

            write_binned_raster(
                out_path,
                template_raster,
                os.path.join(local_root, country, city, 'output', 'priority', f'recreation_priority_binned.tif'),
                binning_method,
                percentiles
            )

        if benefit is Benefits.SocialCohesion:
            tqauto.write('Assessing social cohesion priority')
            social_cohesion_demand = os.path.join(local_root, country, city, 'output', 'social_cohesion_demand.tif')
            social_cohesion_raster = rioxarray.open_rasterio(social_cohesion_demand)
            gdf = gdf.to_crs(social_cohesion_raster.rio.crs)
            social_cohesion_stats = ee.exact_extract(social_cohesion_raster, gdf, 'mean', output='pandas')
            gdf['cohes'] = social_cohesion_stats['mean']   

            out_path = os.path.join(local_root, country, city, 'output', 'priority', f'social_cohesion_priority.tif')
            rasterize_df(
                gdf,
                template_raster,
                out_path,
                'cohes'
            ) 

            write_binned_raster(
                out_path,
                template_raster,
                os.path.join(local_root, country, city, 'output', 'priority', f'social_cohesion_priority_binned.tif'),
                binning_method,                
                percentiles
            )


        if benefit is Benefits.Biodiversity:
            # for biodiversity, we have not yet a clear demand perspective
            pass

        if benefit is Benefits.AvoidedErosion:
            tqauto.write('Assessing erosion priority')
            erosion_demand = os.path.join(local_root, country, city, 'output', 'priority', 'erosion_demand.tif')
            erosion_raster = rioxarray.open_rasterio(erosion_demand)
            gdf = gdf.to_crs(erosion_raster.rio.crs)
            erosion_stats = ee.exact_extract(erosion_raster, gdf, 'mean', output='pandas')
            gdf['erosion'] = erosion_stats['mean']   

            out_path = os.path.join(local_root, country, city, 'output', 'priority', f'erosion_priority.tif')
            rasterize_df(
                gdf,
                template_raster,
                out_path,
                'erosion'
            ) 

            write_binned_raster(
                out_path,
                template_raster,
                os.path.join(local_root, country, city, 'output', 'priority', f'erosion_priority_binned.tif'),
                binning_method,                
                percentiles
            )

        if benefit is Benefits.UrbanHeat:
            tqauto.write('Assessing urban heat mitigation priority')
            heat_demand = os.path.join(local_root, country, city, 'output', 'heat_demand.tif')
            heat_raster = rioxarray.open_rasterio(heat_demand)
            gdf = gdf.to_crs(heat_raster.rio.crs)
            heat_stats = ee.exact_extract(heat_raster, gdf, 'mean', output='pandas')
            gdf['heat'] = heat_stats['mean']

            out_path = os.path.join(local_root, country, city, 'output', 'priority', f'heat_priority.tif')
            rasterize_df(
                gdf,
                template_raster,
                out_path,
                'heat'
            )

            write_binned_raster(
                out_path,
                template_raster,
                os.path.join(local_root, country, city, 'output', 'priority', f'heat_priority_binned.tif'),
                binning_method,
                percentiles
            )

            


    gdf.to_file(out_shape, driver='ESRI Shapefile')
    logger.info('Priority assessment completed')

