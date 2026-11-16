from .biodiversity_exports import PhysiognomicLandscapeType, VegetationLayer, LandCoverClassVegetationProperties, HabitatThreat, ThreatProperties, ThreatDecayFunction, ThreatSourceType, LandCover, LandCoverClassProperties, Metrics
from .biodiversity_spatial import Spatial, PatchContiguity
from .biodiversity_logging import biodiversity_logger

import numpy as np
from typing import Tuple, List
from skimage import data, util, measure
import concurrent.futures
import pandas as pd
import rioxarray as riox
from typing import Dict, List, Tuple
import sys
from scipy import ndimage

class BiodiversityIndicators:
    """This class encapsulates methods for the raster-based computation of indicators. 
    """
    @staticmethod
    def get_habitats_of_physiognomic_type(landscape_type: PhysiognomicLandscapeType, habitat_types: List[LandCoverClassVegetationProperties]) -> List[LandCoverClassVegetationProperties]:        

        if landscape_type is PhysiognomicLandscapeType.Any:
            # if we need any landscape type and do not differentiate between present layers, 
            # simply return the initial list.
            return habitat_types 

        habitats_of_interest = []
        for habitat_type in habitat_types:
            if PhysiognomicLandscapeType.get_landscape_type(habitat_type) is landscape_type:
                habitats_of_interest.append(habitat_type)
        return habitats_of_interest

    

    @staticmethod
    def mean_of_patch(patch_matrix: np.ndarray, in_matrix: np.ndarray) -> pd.DataFrame:
        region_stats = measure.regionprops_table(label_image=patch_matrix, intensity_image=in_matrix, properties=('label', 'area', 'intensity_mean'))
        stats_df = pd.DataFrame(region_stats)  
        stats_df['intensity_mean'] = stats_df['intensity_mean'] * stats_df['area']

        return stats_df
      

    @staticmethod
    def determine_metrics(df_patches: pd.DataFrame, reference_total_area: float) -> Tuple[pd.DataFrame, Dict[Metrics, float]]:
        
        dict_of_metrics = {}
        patch_count = df_patches.shape[0]
        
        df_patches['coh'] = ((df_patches['area'] / reference_total_area) ** 2)
        df_patches['as'] = (df_patches['area'] ** 2) 
        df_patches['reference_area'] = reference_total_area
        df_patches['m'] = ((df_patches['area'] ** 2) / df_patches['reference_area'])

        dict_of_metrics[Metrics.PatchCount] = float(patch_count)        
        dict_of_metrics[Metrics.Coherence] = float(df_patches["coh"].sum())                 
        dict_of_metrics[Metrics.LandscapeDivision] = float(1 - dict_of_metrics[Metrics.Coherence])
        dict_of_metrics[Metrics.SplittingIndex] = float(1 / dict_of_metrics[Metrics.Coherence])
        dict_of_metrics[Metrics.EffectiveMeshSize] = float(df_patches["m"].sum())

        return df_patches, dict_of_metrics

    @staticmethod
    def determine_metrics_for_connected_components(df_patches: pd.DataFrame, reference_total_area: float) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[Metrics, float]]:
        
        dict_of_metrics = {}

        # total number of patches in the landscape, irrespective of connected component membership
        patch_count = df_patches.shape[0]

        df_patches['reference_area'] = reference_total_area


        # to determine m for connected components, we have to group by component label and sum (component_size**2)/reference_area
        aggregated_df = df_patches.groupby('component_label').agg({'component_size': 'mean'}).reset_index()
        aggregated_df['as'] = (aggregated_df['component_size'] ** 2)
        aggregated_df['reference_area'] = reference_total_area
        aggregated_df['m'] = (aggregated_df['as'] / aggregated_df['reference_area'])

        # metrics such as coherence could theoretically be evaluated as before, although they wouldn't be well-defined anyway
        # we could determine C as m/At 
        # we could determine D as 1-C
        # we could determine S as 1/C

        # for comparison purposes, we could also quickly determine effective mesh size unconnected, to make differences more obvious
        df_patches['m'] = ((df_patches['area'] ** 2) / df_patches['reference_area'])

        # compute connected effective mesh size
        connected_effective_mesh_size = float(aggregated_df["m"].sum())

        dict_of_metrics[Metrics.PatchCount] = float(patch_count)
        dict_of_metrics[Metrics.Coherence] = connected_effective_mesh_size / reference_total_area
        dict_of_metrics[Metrics.LandscapeDivision] = float(1 - dict_of_metrics[Metrics.Coherence])
        dict_of_metrics[Metrics.SplittingIndex] = float(1 / dict_of_metrics[Metrics.Coherence])
        dict_of_metrics[Metrics.EffectiveMeshSize] = float(df_patches["m"].sum())
        dict_of_metrics[Metrics.ComponentCount] = float(aggregated_df.shape[0]) 
        dict_of_metrics[Metrics.ConnectedEffectiveMeshSize] = connected_effective_mesh_size

        return df_patches, aggregated_df, dict_of_metrics


    @staticmethod
    def coherence(mtx_patch_size: np.ndarray, total_area: float) -> np.ndarray:

        # determine actual COH at patch level
        # coherence COH for a given patch is defined as (Ai/At) ** 2
        mtx_coherence = np.zeros(mtx_patch_size.shape, dtype=np.float32)
        np.divide(mtx_patch_size, total_area, out=mtx_coherence)
        #np.pow(mtx_coherence, 2, out=mtx_coherence)        
        mtx_coherence = (mtx_coherence ** 2)
        return mtx_coherence
        
    @staticmethod
    def effective_mesh_size(mtx_patch_size: np.ndarray, total_area: float) -> np.ndarray:

        # determine patch contribution to effective mesh size m
        # patch contribution to m is 1/At * (patch size^^2 + ...)
        mtx_mesh_size = np.zeros(mtx_patch_size.shape, dtype=np.float32)
        #np.pow(mtx_patch_size, 2, out=mtx_mesh_size)        
        mtx_mesh_size = (mtx_patch_size ** 2)
        np.divide(mtx_mesh_size, total_area, out=mtx_mesh_size)
        return mtx_mesh_size
    
    @staticmethod
    def connected_effective_mesh_size(mtx_patch_labels: np.ndarray, df_patches: pd.DataFrame, mtx_patch_size: np.ndarray, mtx_component_labels: np.ndarray, total_area: float) -> Tuple[np.ndarray, np.ndarray]: 
        
        # result matrixes
        mtx_corrected_component_labels = mtx_component_labels.copy()
        mtx_mesh_size = np.zeros(mtx_patch_size.shape, dtype=np.float32)
        mtx_number_of_patches_in_component = np.zeros(mtx_patch_size.shape, dtype=np.float32)


        # add required columns in patch dataframe
        df_patches['component_label'] = 0
        df_patches['component_size'] = 0.0

        # Build fast label-indexed lookups to avoid repeated dataframe scans.
        patch_labels = df_patches['label'].to_numpy(dtype=np.int32)
        patch_areas = df_patches['area'].to_numpy(dtype=np.float32)
        max_patch_label = int(patch_labels.max()) if patch_labels.size > 0 else 0

        area_by_label = np.zeros(max_patch_label + 1, dtype=np.float32)
        area_by_label[patch_labels] = patch_areas

        component_by_label = np.zeros(max_patch_label + 1, dtype=np.int32)
        component_size_by_label = np.zeros(max_patch_label + 1, dtype=np.float32)

        # valid patch pixels are pixels that belong to a non-background patch
        valid_pixel_mask = mtx_patch_labels > 0
        if np.any(valid_pixel_mask):
            flat_component_labels = mtx_component_labels[valid_pixel_mask].astype(np.int32)
            flat_patch_labels = mtx_patch_labels[valid_pixel_mask].astype(np.int32)

            # unique (component, patch) pairs remove pixel-level duplication
            comp_patch_pairs = np.unique(
                np.stack((flat_component_labels, flat_patch_labels), axis=1),
                axis=0
            )

            # Component assignment for each patch label (max preserves prior loop semantics).
            np.maximum.at(component_by_label, comp_patch_pairs[:, 1], comp_patch_pairs[:, 0])

            # count distinct patches per component
            components_with_patch_refs = comp_patch_pairs[:, 0]
            max_component_label = int(mtx_component_labels.max()) if mtx_component_labels.size > 0 else 0
            patch_count_by_component = np.zeros(max_component_label + 1, dtype=np.float32)
            unique_components_with_refs, counts_per_component = np.unique(components_with_patch_refs, return_counts=True)
            patch_count_by_component[unique_components_with_refs] = counts_per_component.astype(np.float32)

            # sum patch areas by assigned component
            patch_component_refs = component_by_label[patch_labels]
            component_size_by_component = np.bincount(
                patch_component_refs,
                weights=patch_areas,
                minlength=max_component_label + 1
            ).astype(np.float32)
            component_size_by_label = component_size_by_component[component_by_label]

            # write rasters in one pass
            mtx_number_of_patches_in_component[valid_pixel_mask] = patch_count_by_component[flat_component_labels]
            mtx_mesh_size[valid_pixel_mask] = component_size_by_component[flat_component_labels]

            # components without associated patches are set to background (0)
            has_patches_by_component = np.zeros(max_component_label + 1, dtype=bool)
            has_patches_by_component[unique_components_with_refs] = True
            invalid_component_mask = (
                (mtx_corrected_component_labels > 0)
                & (~has_patches_by_component[mtx_corrected_component_labels])
            )
            mtx_corrected_component_labels[invalid_component_mask] = 0

        # write final component attributes back to patch dataframe in one vectorized pass
        df_patches['component_label'] = component_by_label[patch_labels]
        df_patches['component_size'] = component_size_by_label[patch_labels]
        
        # at this point, component sizes are written to the matrix
        # now, compute effective mesh size considering connected components
        mtx_mesh_size = (mtx_mesh_size ** 2)
        # now we need to divide by the total area of the component to get the patch-specific contribution to total effective mesh size
        np.divide(mtx_mesh_size, total_area, out=mtx_mesh_size)

        return mtx_corrected_component_labels, mtx_mesh_size, mtx_number_of_patches_in_component, df_patches
        

    @staticmethod
    def get_vegetation_layer_weight_mask(mtx_land_cover: np.ndarray, habitat_values: List[LandCoverClassVegetationProperties]) -> np.ndarray:
        
        scaling_mtx = np.zeros(mtx_land_cover.shape, np.float32)        
        for habitat_value in habitat_values:
            # result pixel values is input * (layer count/3)
            habitat_scale_factor = (float(habitat_value.get_vegetation_layer_count()) / 3)
            scaling_mtx[mtx_land_cover == habitat_value.get_land_cover_class()] = habitat_scale_factor 
        return scaling_mtx
    
    @staticmethod
    def vegetation_layer_scaling(mtx_parameter_to_scale: np.ndarray, mtx_scale: np.ndarray) -> np.ndarray:        
        return np.multiply(mtx_parameter_to_scale, mtx_scale)
    
    @staticmethod
    def habitat_quality(mtx_land_cover, geodf_roads, sensitivity: List[LandCover], threats: List[HabitatThreat], resolution:Tuple[int,int], transform, k: float = 0.05, nodata_values: List[int] = [0]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        
        # some predefined constants
        z = 2.5        
        k_scaled = k**z
        
        # sanity check
        if threats is None or sensitivity is None:
            return None
        
        else:
            
            # define sum of degradation result raster
            sum_of_degradation = np.zeros(mtx_land_cover.shape, np.float32)

            # make habitat quality value surface
            mtx_habitat_suitability = np.zeros(mtx_land_cover.shape, np.float32)
            for land_cover in sensitivity:
                current_mask = mtx_land_cover == land_cover.value.get_land_cover_class_value()
                mtx_habitat_suitability[current_mask] = land_cover.value.get_habitat_suitability()

            # sum of threat weights
            sum_of_threat_weights = sum([f.value.get_relative_weight() for f in threats])
           
            # iterate over all threats
            current_threat: HabitatThreat
            for current_threat in threats:
                biodiversity_logger.info(f"Assessing threat {current_threat.name}")

                # make threat source layer (presence/absence) derived from relevant input data sources
                mtx_threat = np.zeros(mtx_land_cover.shape, np.float32)
                if current_threat.value.get_source_type() == ThreatSourceType.Land_Cover:
                    threat_class_mask = np.isin(mtx_land_cover, [current_threat.value.get_value_for_selection()], invert=False)
                    mtx_threat[threat_class_mask] = 1

                elif current_threat.value.get_source_type() == ThreatSourceType.Road_Barriers:
                    mtx_road_type = Spatial.make_road_layer(geodf_roads, current_threat.value.get_value_for_selection(), mtx_land_cover.shape, transform)
                    mtx_threat += mtx_road_type

                # make decay surface (linear or exponential) depending on distance/proximity to current threat sources
                decay_surface = Spatial.get_threat_decay_surface(threat=current_threat.value, mtx_input=mtx_threat, resolution=resolution)                    

                # make sensitivity surface for the current threat
                sensitivity_surface = np.zeros(mtx_land_cover.shape, np.float32)                
                for land_cover in sensitivity:                
                    current_mask = mtx_land_cover == land_cover.value.get_land_cover_class_value()
                    sensitivity_surface[current_mask] = land_cover.value.get_sensitivity_to_threat(current_threat)

                # make normalized weight for current threat
                normalized_weight = current_threat.value.get_relative_weight() / sum_of_threat_weights

                sum_of_degradation += (decay_surface * sensitivity_surface * normalized_weight)
            
            # continue with the assessment: compute habitat quality
            habitat_quality = mtx_habitat_suitability * (1 - ((sum_of_degradation**z)/((sum_of_degradation**z) + k_scaled)))

        return mtx_habitat_suitability, sum_of_degradation, habitat_quality