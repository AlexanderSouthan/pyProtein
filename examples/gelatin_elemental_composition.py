# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 13:47:49 2026

@author: southan
"""


import pandas as pd
from pyProtein import amino_acids as gelatin_amino_acids

# write gelatin amino acid abundances into amino acid table
gelatin_amino_acids.at['A', 'per1000'] = 117 # alanine
gelatin_amino_acids.at['R', 'per1000'] = 48 # arginine
gelatin_amino_acids.at['N', 'per1000'] = 0 # asparagine
gelatin_amino_acids.at['D', 'per1000'] = 46 # aspartic acid
gelatin_amino_acids.at['C', 'per1000'] = 0 # cysteine
gelatin_amino_acids.at['E', 'per1000'] = 72 # glutamic acid
gelatin_amino_acids.at['Q', 'per1000'] = 0 # glutamine
gelatin_amino_acids.at['G', 'per1000'] = 335 # glycine
gelatin_amino_acids.at['H', 'per1000'] = 4.2 # histidine
gelatin_amino_acids.at['Hyp', 'per1000'] = 93 # hydroxyproline
gelatin_amino_acids.at['Hyl', 'per1000'] = 4.3 # hydroxylysine
gelatin_amino_acids.at['I', 'per1000'] = 11 # isoleucine
gelatin_amino_acids.at['L', 'per1000'] = 24.3 # leucine
gelatin_amino_acids.at['K', 'per1000'] = 28 # lysine
gelatin_amino_acids.at['M', 'per1000'] = 3.9 # methionine
gelatin_amino_acids.at['F', 'per1000'] = 14 # phenylalanine
gelatin_amino_acids.at['P', 'per1000'] = 124 # proline
gelatin_amino_acids.at['S', 'per1000'] = 33 # serine
gelatin_amino_acids.at['T', 'per1000'] = 18 # threonine
gelatin_amino_acids.at['W', 'per1000'] = 0 # tryptophan
gelatin_amino_acids.at['Y', 'per1000'] = 1.2 # tyrosine
gelatin_amino_acids.at['V', 'per1000'] = 22 # valine


# multiply amino acid abundances with number of the different atoms in amino acids
gelatin_amino_acids['total_c_per1000'] = gelatin_amino_acids['per1000'] * gelatin_amino_acids['C']
gelatin_amino_acids['total_h_per1000'] = gelatin_amino_acids['per1000'] * gelatin_amino_acids['H']
gelatin_amino_acids['total_n_per1000'] = gelatin_amino_acids['per1000'] * gelatin_amino_acids['N']
gelatin_amino_acids['total_o_per1000'] = gelatin_amino_acids['per1000'] * gelatin_amino_acids['O']
gelatin_amino_acids['total_s_per1000'] = gelatin_amino_acids['per1000'] * gelatin_amino_acids['S']


# sum up the number of different atoms for all amino acids
gelatin_elemental_composition = pd.DataFrame( # in total number per 1000 residues
    [gelatin_amino_acids['total_c_per1000'].sum(),
     gelatin_amino_acids['total_h_per1000'].sum(),
     gelatin_amino_acids['total_n_per1000'].sum(),
     gelatin_amino_acids['total_o_per1000'].sum(),
     gelatin_amino_acids['total_s_per1000'].sum()],
    index=['C', 'H', 'N', 'O', 'S'], columns=['atom count'])


gelatin_elemental_composition['molar mass [g/mol]'] = [12.011, 1.008, 14.007, 15.999, 32.07]
gelatin_elemental_composition['relative element masses'] = (
    gelatin_elemental_composition['atom count'] *
    gelatin_elemental_composition['molar mass [g/mol]'])
total_mass = gelatin_elemental_composition['relative element masses'].sum()