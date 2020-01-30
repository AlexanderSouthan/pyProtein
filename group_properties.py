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

amino_acids['long_name'] = ['aspartic acid', 'asparagine', 'threonine',
           'serine', 'glutamic acid', 'glutamine', 'glycine', 'alanine',
           'cysteine', 'valine', 'methionine', 'isoleucine', 'leucine',
           'tyrosine', 'phenylalanine', 'histidine', 'lysine', 'arginine',
           'proline', 'tryptophan', 'N-terminus', 'C-terminus']

amino_acids['three_letter_code'] = ['Asp', 'Asn', 'Thr', 'Ser', 'Glu', 'Gln',
           'Gly', 'Ala', 'Cys', 'Val', 'Met', 'Ile', 'Leu', 'Tyr', 'Phe',
           'His', 'Lys', 'Arg', 'Pro', 'Trp', None, None]

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
amino_acids['pka_bjellqvist'] = np.nan
amino_acids.at['D', 'pka_bjellqvist'] = 4.05
amino_acids.at['E', 'pka_bjellqvist'] = 4.45
amino_acids.at['H', 'pka_bjellqvist'] = 5.98
amino_acids.at['C', 'pka_bjellqvist'] = 9
amino_acids.at['Y', 'pka_bjellqvist'] = 10
amino_acids.at['K', 'pka_bjellqvist'] = 10
amino_acids.at['R', 'pka_bjellqvist'] = 12
amino_acids.at['N_term', 'pka_bjellqvist'] = 7.5
amino_acids.at['C_term', 'pka_bjellqvist'] = 3.55


# according to Kozlowski, Biology Direct 2016, 11 (1), 55.
# DOI: 10.1186/s13062-016-0159-9.
amino_acids['pka_ipc_protein'] = np.nan
amino_acids.at['D', 'pka_ipc_protein'] = 3.872
amino_acids.at['E', 'pka_ipc_protein'] = 4.412
amino_acids.at['H', 'pka_ipc_protein'] = 5.637
amino_acids.at['C', 'pka_ipc_protein'] = 7.555
amino_acids.at['Y', 'pka_ipc_protein'] = 10.85
amino_acids.at['K', 'pka_ipc_protein'] = 9.052
amino_acids.at['R', 'pka_ipc_protein'] = 11.84
amino_acids.at['N_term', 'pka_ipc_protein'] = 9.094
amino_acids.at['C_term', 'pka_ipc_protein'] = 2.869


# according to http://emboss.sourceforge.net/
amino_acids['pka_emboss'] = np.nan
amino_acids.at['D', 'pka_emboss'] = 3.9
amino_acids.at['E', 'pka_emboss'] = 4.1
amino_acids.at['H', 'pka_emboss'] = 6.5
amino_acids.at['C', 'pka_emboss'] = 8.5
amino_acids.at['Y', 'pka_emboss'] = 10.1
amino_acids.at['K', 'pka_emboss'] = 10.8
amino_acids.at['R', 'pka_emboss'] = 12.5
amino_acids.at['N_term', 'pka_emboss'] = 8.6
amino_acids.at['C_term', 'pka_emboss'] = 3.6

[3.9, np.nan, np.nan, np.nan, 4.412, np.nan,
           np.nan, np.nan, 7.555, np.nan, np.nan, np.nan, np.nan, 10.85,
           np.nan, 5.637, 9.052, 11.84, np.nan, np.nan, 8.6, 3.6]