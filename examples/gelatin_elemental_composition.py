# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 13:47:49 2026

@author: southan
"""


import numpy as np
import pandas as pd
from pyProtein import amino_acids as gelatin_amino_acids
from pyProtein import protein

# order of amino acid abundances:
# ['D', 'N', 'T', 'S', 'E',
# 'Q', 'G', 'A', 'C', 'V', 'M', 'I', 'L', 'Y', 'F', 'H', 'K',
# 'R', 'P', 'W', 'Hyp'm 'Hyl']

gel_b_handbook = protein(
    'res_per_1000',
    [46, 0, 18, 33, 72, 0, 335, 117, 0, 22, 3.9, 11, 24.3, 1.2, 14, 4.2, 28, 48, 124, 0, 93, 4.3])

gel_b_sewald = protein(
    'res_per_1000',
    [45, 0, 16, 32, 74, 0, 341, 116, 0, 21, 5, 12, 25, 1, 12, 4, 27, 48, 113, 0, 102, 6])

gel_b_sewald_mmol = protein(
    'mmol_g',
    np.array([0.4968341550984979, 0.0, 0.17665214403502147, 0.35330428807004294, 
     0.8170161661619744, 0.0, 3.7648988197463953, 1.2807280442539057, 0.0,
    0.23185593904596571, 0.05520379501094421, 0.13248910802626612,
    0.27601897505472106, 0.011040759002188842, 0.13248910802626612,
    0.04416303600875537, 0.29810049305909875, 0.5299564321050645,
    1.247605767247339, 0.0, 1.1261574182232619, 0.06624455401313306])*0.9796,
    mod_types=['methacryl'], mod_abundances=[0.3])

gel_b_sewald_mod = protein(
    'res_per_1000',
    [45, 0, 16, 32, 74, 0, 341, 116, 0, 21, 5, 12, 25, 1, 12, 4, 27, 48, 113, 0, 102, 6],
    mod_types=['methacryl'], mod_abundances=[33])

print('Gelatin type B (handbook)\n', gel_b_handbook.elemental_composition())
print('Gelatin type B (Sewald)\n', gel_b_sewald.elemental_composition())
A = gel_b_sewald_mod.elemental_composition()
print('Gelatin type B (Sewald) modified\n', A)