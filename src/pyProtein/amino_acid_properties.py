# -*- coding: utf-8 -*-
"""
Collection of properties for different moieties present in proteins.
"""

import numpy as np
import pandas as pd


# Start of amino acid data.
amino_acids = pd.DataFrame(
        [], index=['D', 'N', 'T', 'S', 'E', 'Q', 'G', 'A', 'C', 'V', 'M', 'I',
                   'L', 'Y', 'F', 'H', 'K', 'R', 'P', 'W', 'Hyp', 'Hyl'])

amino_acids['long_name'] = ['aspartic acid', 'asparagine', 'threonine',
                            'serine', 'glutamic acid', 'glutamine', 'glycine',
                            'alanine', 'cysteine', 'valine', 'methionine',
                            'isoleucine', 'leucine', 'tyrosine',
                            'phenylalanine', 'histidine', 'lysine', 'arginine',
                            'proline', 'tryptophan', 'hydroxyproline',
                            'hydroxylysine']

amino_acids['three_letter_code'] = ['Asp', 'Asn', 'Thr', 'Ser', 'Glu', 'Gln',
                                    'Gly', 'Ala', 'Cys', 'Val', 'Met', 'Ile',
                                    'Leu', 'Tyr', 'Phe', 'His', 'Lys', 'Arg',
                                    'Pro', 'Trp', 'Hyp', 'Hyl']

amino_acids['charge_indicator'] = [-1, 0, 0, 0, -1, 0, 0, 0, -1, 0, 0, 0, 0,
                                   -1, 0, 1, 1, 1, 0, 0, 0, 1]

amino_acids['C'] = [4, 4, 4, 3, 5, 5, 2, 3, 3, 5, 5, 6, 6, 9, 9, 6, 6, 6, 5,
                    11, 5, 6]
amino_acids['H'] = [7, 8, 9, 7, 9, 10, 5, 7, 7, 11, 11, 13, 13, 11, 11, 9, 14,
                    14, 9, 12, 9, 14]
amino_acids['N'] = [1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 2, 4, 1, 2,
                    1, 2]
amino_acids['O'] = [4, 3, 3, 3, 4, 3, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2,
                    3, 3]
amino_acids['S'] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0]

amino_acids['molar_mass'] = (
    12.0107 * amino_acids['C'] +
    1.0079 * amino_acids['H'] +
    14.0067 * amino_acids['N'] +
    15.9994 * amino_acids['O'] +
    32.065 * amino_acids['S'])

amino_acids['N_content'] = (
    14.0067 * amino_acids['N'] / amino_acids['molar_mass'])

amino_acids['molar_mass_residue'] = (amino_acids['molar_mass'] -
                                     2*1.0079 - 15.9994)

amino_acids['N_content_residue'] = (
    14.0067 * amino_acids['N'] / amino_acids['molar_mass_residue'])

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
amino_acids.at['Hyl', 'pka_bjellqvist'] = 10


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
amino_acids.at['Hyl', 'pka_ipc_protein'] = 9.052


# according to http://emboss.sourceforge.net/
amino_acids['pka_emboss'] = np.nan
amino_acids.at['D', 'pka_emboss'] = 3.9
amino_acids.at['E', 'pka_emboss'] = 4.1
amino_acids.at['H', 'pka_emboss'] = 6.5
amino_acids.at['C', 'pka_emboss'] = 8.5
amino_acids.at['Y', 'pka_emboss'] = 10.1
amino_acids.at['K', 'pka_emboss'] = 10.8
amino_acids.at['R', 'pka_emboss'] = 12.5
amino_acids.at['Hyl', 'pka_emboss'] = 10.8

# Start of data for chain modifications such as end groups (N terminus,
# C terminus) or side chain modifications
chain_modifications = pd.DataFrame([], index=['N_term', 'C_term'])  # ,
                                              # 'Methacryl', 'Acetyl'])
chain_modifications['long_name'] = ['N-terminus', 'C-terminus']  # ,
                                    # 'Methacryl group', 'Acetyl group']

chain_modifications['charge_indicator'] = [1, -1]  # , 0, 0]

chain_modifications['C'] = [0, 0]  # , 4, 2]
chain_modifications['H'] = [1, 1]  # , 4, 2]
chain_modifications['N'] = [0, 0]  # , 0, 0]
chain_modifications['O'] = [0, 1]  # , 1, 1]

chain_modifications['molar_mass_residue'] = (
    12.0107 * chain_modifications['C'] +
    1.0079 * chain_modifications['H'] +
    14.0067 * chain_modifications['N'] +
    15.9994 * chain_modifications['O'])

chain_modifications['N_content_residue'] = (
    14.0067 * chain_modifications['N'] / chain_modifications[
        'molar_mass_residue'])

chain_modifications['pka_bjellqvist'] = np.nan
chain_modifications.at['N_term', 'pka_bjellqvist'] = 7.5
chain_modifications.at['C_term', 'pka_bjellqvist'] = 3.55

chain_modifications['pka_ipc_protein'] = np.nan
chain_modifications.at['N_term', 'pka_ipc_protein'] = 9.094
chain_modifications.at['C_term', 'pka_ipc_protein'] = 2.869

chain_modifications['pka_emboss'] = np.nan
chain_modifications.at['N_term', 'pka_emboss'] = 8.6
chain_modifications.at['C_term', 'pka_emboss'] = 3.6

# This could hold the pKa values of acidic or basic modifications in the
# future. Currently it only makes sure that if no valid pKa scale is geven for
# C and N terminus, they are ignored.
chain_modifications['pka_other'] = np.nan
