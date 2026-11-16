
from . import invest_flood, gis_tools

import rasterio as rio
import numpy as np
import os
from tqdm import tqdm
from tqdm.auto import tqdm as tqauto

class Runoff:

    def __init__(self, local_root, country, city):

        self.root = os.path.join(local_root, country, city)
        self.risk = os.path.join(self.root, 'output', 'flood_risk.tif')
        self.output_folder = os.path.join(self.root, 'flood_benefit')
        self.input_folder = os.path.join(self.root, 'cn')
        self.scratch_folder = os.path.join(self.root, 'temp')



    def nbs_invest(self, cn, storm_depth=150):
        cn_path = os.path.join(self.input_folder, cn)
        name = cn[3:-4] # remove the cn- and .tif from the filename for the NBS name
        invest = invest_flood.InvestFlood(name, cn_path, storm_depth, self.scratch_folder)
        invest.run()
        invest.export_layers(runoff=False,
                             retention=False,
                             runoff_volume=True,
                             retention_volume=False)

    def get_runoff(self):
        cns = [x for x in os.listdir(self.input_folder) if x.endswith(".tif")]
        baseline_cn = [cn for cn in cns if 'baseline' in cn][0]
        cns.remove(baseline_cn)

        self.nbs_invest(baseline_cn)
        baseline_runoff = os.path.join(self.scratch_folder, 'baseline_runoff_volume.tif')
        flooder = FloodPotential(baseline_runoff, self.scratch_folder, self.risk, self.risk)

        for cn in tqauto(cns):
            name = cn[3:-4] # remove the cn- and .tif from the filename for the NBS name
            print(name)
            self.nbs_invest(cn)
            flooder.get_potential(os.path.join(self.scratch_folder, f'{name}_runoff_volume.tif'),
                                  os.path.join(self.output_folder, f'flood_benefit_{name}.tif'))


class FloodPotential:

    def __init__(self, base_runoff, scratch, potential, mask):
        self.base_runoff = base_runoff
        self.potential = potential
        self.scratch = scratch
        self.mask = mask

    def get_dif(self, nbs_runoff):
        with rio.open(self.base_runoff) as src:
            base = src.read(1)
        with rio.open(nbs_runoff) as src:
            nbs = src.read(1)

        dif = base - nbs

        return dif

    def get_potential(self, nbs_runoff, out_ras):
        dif = self.get_dif(nbs_runoff)
        with rio.open(self.potential) as src:
            pot = src.read(1)
        pot = np.where(pot < 0, 0, pot)
        #pot = np.where(pot > 10000, 0, pot)

        effect = dif * pot

        effect1 = np.where(effect < 0, 0, effect)
        gis_tools.export_raster(effect1, self.mask, out_ras)
