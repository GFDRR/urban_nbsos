from .gis_tools import export_raster

from importlib import resources as impresources

import pandas as pd
import rasterio as rio
import numpy as np
import os

class CN:

    def __init__(self, local_root, country, city):
        
        # nbs_csv = nbs_cn.csv
        # scenario_csv = cn_sim.csv 
        
        my_resources = impresources.files("nbsos.urban_v1")
        
        nbs_csv = my_resources.joinpath("data", "nbs_cn.csv")
        scenario_csv = my_resources.joinpath("data", "cn_sim.csv")        

        self.nbs = pd.read_csv(nbs_csv)
        self.scenarios = pd.read_csv(scenario_csv)
        self.base = os.path.join(local_root, country, city)
        self.land_cover = self.open_raster(os.path.join(self.base, 'eo_indicators', 'land_cover.tif'))
        self.bs_freq = self.open_raster(os.path.join(self.base, 'eo_indicators', 'bs_freq.tif'))

        # get hydrosoil, with water as '4'
        hydrosoil = self.open_raster(os.path.join(self.base, 'eo_indicators', 'hydrosoil.tif'))
        try:
            water = self.open_raster(os.path.join(self.base, 'eo_indicators', 'water_mask.tif'))
        except Exception as e: # some areas have no water
            water = np.zeros(self.bs_freq.shape)
        self.hydrosoil = np.where((water == 1) | (hydrosoil == 0), 4, hydrosoil)

        # test nan-mask correction similar to jenks clustering based on hydrosoil layer
        self.nan_mask = np.isnan(self.hydrosoil)
        self.raster_shape = self.land_cover.shape

        self.land_cover = self.land_cover[~self.nan_mask]
        self.bs_freq = self.bs_freq[~self.nan_mask]
        self.hydrosoil = self.hydrosoil[~self.nan_mask]


    def open_raster(self, path):
        with rio.open(path) as src:
            ras = src.read(1)
        return ras

    def get_curve_number(self, nbs):
        creation_path = os.path.join(self.base, 'suitability', f'{nbs}_creation.tif')
        nbs_creation = self.open_raster(creation_path)

        nbs_creation = nbs_creation[~self.nan_mask]

        bs_freq = self.nbs.loc[self.nbs['nbs'] == nbs]['bs_freq'].item()
        land_cover = self.nbs.loc[self.nbs['nbs'] == nbs]['land_cover'].item()
        bs_freq_sim = np.where(nbs_creation > 0, bs_freq, self.bs_freq)
        land_cover_sim = np.where(nbs_creation > 0, land_cover, self.land_cover)
        bs_freq_sim_recl = self.reclass_bs_freq(bs_freq_sim, land_cover_sim)
        hydro_lc = (self.hydrosoil * 10) + land_cover_sim
        sim_class = (bs_freq_sim_recl * 100) + hydro_lc
        sim_class = sim_class.astype('int32')
        sim_class = np.where(sim_class < 0, 0, sim_class)  # fix for nan conversion
        lookup = [(row['sim_class'], row['cn']) for idx, row in self.scenarios.iterrows()]
        idx, val = np.asarray(lookup).T
        idx = idx.astype('int32')
        lookup_array = np.zeros(idx.max() + 1)
        lookup_array[idx] = val

        curve_number_sim = lookup_array[sim_class]

        return curve_number_sim

    def get_curve_number_baseline(self):

        bs_freq_sim = self.bs_freq
        land_cover_sim = self.land_cover
        bs_freq_sim_recl = self.reclass_bs_freq(bs_freq_sim, land_cover_sim)
        hydro_lc = (self.hydrosoil * 10) + land_cover_sim
        sim_class = (bs_freq_sim_recl * 100) + hydro_lc
        sim_class = sim_class.astype('int32')
        sim_class = np.where(sim_class < 0, 0, sim_class) # fix for nan conversion
        lookup = [(row['sim_class'], row['cn']) for idx, row in self.scenarios.iterrows()]
        idx, val = np.asarray(lookup).T
        idx = idx.astype('int32')
        lookup_array = np.zeros(idx.max() + 1)
        lookup_array[idx] = val

        curve_number_sim = lookup_array[sim_class]

        return curve_number_sim

    def get_cns(self, nbs_list):
        cn_baseline = self.get_curve_number_baseline()

        cn_baseline_full = np.zeros(self.raster_shape, dtype=np.float32)
        cn_baseline_full[~self.nan_mask] = cn_baseline

        mask = os.path.join(self.base, 'eo_indicators', 'land_cover.tif')
        output_folder = os.path.join(self.base, 'cn')
        export_raster(cn_baseline_full, mask, os.path.join(output_folder, 'cn_baseline.tif'))
        for nbs in nbs_list:
            cn = self.get_curve_number(nbs)
            cn_full = np.zeros(self.raster_shape, dtype=np.float32)
            cn_full[~self.nan_mask] = cn

            export_raster(cn_full, mask, os.path.join(output_folder, f'cn_{nbs}.tif'))

    def reclass_bs_freq(self, bs_freq, land_cover):
        scens = self.scenarios.loc[self.scenarios['Hydrosoil'] == 1]
        lc_classes = scens['land_cover'].unique()
        bs_freq_recl = bs_freq.copy()
        for lc_class in lc_classes:
            bs_mins = scens.loc[scens['land_cover'] == lc_class]['bs_freq_min'].unique()
            bs_maxs = scens.loc[scens['land_cover'] == lc_class]['bs_freq_max'].unique()
            bs_cats = scens.loc[scens['land_cover'] == lc_class]['bs_freq_cat'].unique()

            for bs_min, bs_max, bs_cat in zip(bs_mins, bs_maxs, bs_cats):
                bs_freq_recl = np.where((land_cover == lc_class) & (bs_freq >= bs_min) & (bs_freq <= bs_max),
                                   bs_cat, bs_freq_recl)
        bs_freq_recl = np.where(land_cover == 6, 0, bs_freq_recl)

        return bs_freq_recl