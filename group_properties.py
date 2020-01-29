# -*- coding: utf-8 -*-
"""
Collection of properties for different moieties present in polyampholytes.
"""

import numpy as np
import pandas as pd


amino_acids = pd.DataFrame(
        [], index=['D', 'N', 'T', 'S', 'E', 'Q', 'G', 'A', 'C', 'V', 'M', 'I',
        'L', 'Y', 'F', 'H', 'K', 'R', 'P', 'W', 'N_term', 'C_term'])

amino_acids['long_name'] = ['aspartic acid', 'asparagine', 'threonine',
           'serine', 'glutamic acid', 'glutamine', 'glycine', 'alanine',
           'cysteine', 'valine', 'methionine', 'isoleucine', 'leucine',
           'tyrosine', 'phenylalanine', 'histidine', 'lysine', 'arginine',
           'proline', 'tryptophan', 'N-terminus', 'C-terminus']

amino_acids['charge_indicator'] = [-1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, 0,
           -1, 0, 1, 1, 1, 0, 0, 1, -1]

amino_acids['molecular_weight'] = [133.1032, 132.1184, 119.1197, 105.0930,
           147.1299, 146.1451, 75.0669, 89.0935, 121.1590, 117.1469, 149.2124,
           131.1736, 131.1736, 181.1894, 165.1900, 155.1552, 146.1882,
           115.1310, 204.2262, 174.2017, np.nan, np.nan]

# according to Bjellqvist et al., Electrophoresis 1994, 15 (1), 529-539.
# DOI: 10.1002/elps.1150150171.
# pKa for N_term taken for glycine at N_term, other values are in the reference
# and could be added here.
amino_acids['pKa_bjellqvist'] = [4.05, np.nan, np.nan, np.nan, 4.45, np.nan,
           np.nan, np.nan, 9, np.nan, np.nan, np.nan, np.nan, 10, np.nan, 5.98,
           10, 12, np.nan, np.nan, 7.5, 3.55]

amino_acids['pKa_IPC_protein'] = [3.872, np.nan, np.nan, np.nan, 4.412, np.nan,
           np.nan, np.nan, 7.555, np.nan, np.nan, np.nan, np.nan, 10.85,
           np.nan, 5.637, 9.052, 11.84, np.nan, np.nan, 9.094, 2.869]
