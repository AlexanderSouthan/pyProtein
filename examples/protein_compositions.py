# -*- coding: utf-8 -*-
"""
Amino acid analysis results for gelatin type A (pig skin) and type B (bovine
bone) taken from:
Sewald et al., Macromol. Biosci. 2018, 18 (12), 1800168.
DOI: 10.1002/mabi.201800168.

Correction factors for asparagine and glutamine abundances in gelatin type A
derived from:
R. Schrieber, H. Gareis, Gelatine Handbook: Theory and Industrial Practice,
John Wiley & Sons, Weinheim (Germany), 2007.
"""

import numpy as np
import pandas as pd
from pyPolyampholyte.group_properties import amino_acids

proteins = {'aa_data': amino_acids}
proteins['gelatin_type_a_porcine'] = pd.DataFrame([], index=proteins['aa_data'].index)
proteins['gelatin_type_b_bovine'] = pd.DataFrame([], index=proteins['aa_data'].index)
# These are real amino acid results, so values for glutamine and asparagine are
# incorrect. They were hydrolyzed during analysis.
proteins['gelatin_type_a_porcine']['wt%'] = [
    5.91, 0, 1.72, 3.45, 10.50, 0, 24.89, 9.24, 0, 2.42, 0.90, 1.23, 2.92,
    0.60, 2.08, 0.75, 3.79, 8.42, 13.58, 0, 11.82, 1.31]

proteins['gelatin_type_b_bovine']['wt%'] = [
    5.68, 0, 1.84, 3.22, 10.31, 0, 24.39, 9.84, 0, 2.35, 0.74, 1.46, 3.12,
    0.18, 1.92, 0.59, 3.70, 8.04, 12.37, 0, 12.73, 1.22]

# Conversion of mass-% to values in mmol/g.
proteins['gelatin_type_a_porcine']['mmol/g'] = (
    proteins['gelatin_type_a_porcine']['wt%']/proteins['aa_data']['molar_mass']*10)

proteins['gelatin_type_b_bovine']['mmol/g'] = (
    proteins['gelatin_type_b_bovine']['wt%']/proteins['aa_data']['molar_mass']*10)

# Correction of wrong asparagine and glutamine values with literature values
proteins['gelatin_type_a_porcine'].at['N', 'mmol/g'] = (
    proteins['gelatin_type_a_porcine'].at['D', 'mmol/g'] * 16/(16+29))
proteins['gelatin_type_a_porcine'].at['D', 'mmol/g'] -= (
    proteins['gelatin_type_a_porcine'].at['N', 'mmol/g'])
proteins['gelatin_type_a_porcine'].at['Q', 'mmol/g'] = (
    proteins['gelatin_type_a_porcine'].at['E', 'mmol/g'] * 25/(25+48))
proteins['gelatin_type_a_porcine'].at['E', 'mmol/g'] -= (
    proteins['gelatin_type_a_porcine'].at['Q', 'mmol/g'])

# Normalize gelatin amino acid abundances
proteins['gelatin_type_a_porcine']['mol% (norm)'] = (100 *
    proteins['gelatin_type_a_porcine']['mmol/g'] /
    np.sum(proteins['gelatin_type_a_porcine']['mmol/g']))
proteins['gelatin_type_b_bovine']['mol% (norm)'] = (100 *
    proteins['gelatin_type_b_bovine']['mmol/g']/
    np.sum(proteins['gelatin_type_b_bovine']['mmol/g']))

proteins['gelatin_type_a_porcine']['per 1000 residues'] = proteins[
    'gelatin_type_a_porcine']['mol% (norm)'] * 10
proteins['gelatin_type_b_bovine']['per 1000 residues'] = proteins[
    'gelatin_type_b_bovine']['mol% (norm)'] * 10

proteins['gelatin_type_a_porcine']['mmol/g (norm)'] = (1/np.sum(
    proteins['gelatin_type_a_porcine']['mol% (norm)'] *
    proteins['aa_data']['molar_mass']) *
    proteins['gelatin_type_a_porcine']['mol% (norm)'] * 1000)

proteins['gelatin_type_b_bovine']['mmol/g (norm)'] = (1/np.sum(
    proteins['gelatin_type_b_bovine']['mol% (norm)'] *
    proteins['aa_data']['molar_mass']) *
    proteins['gelatin_type_b_bovine']['mol% (norm)'] * 1000)

proteins['gelatin_type_a_porcine']['mmol/g (norm, per residue)'] = (1/np.sum(
    proteins['gelatin_type_a_porcine']['mol% (norm)'] *
    proteins['aa_data']['molar_mass_residue']) *
    proteins['gelatin_type_a_porcine']['mol% (norm)'] * 1000)

proteins['gelatin_type_b_bovine']['mmol/g (norm, per residue)'] = (1/np.sum(
    proteins['gelatin_type_b_bovine']['mol% (norm)'] *
    proteins['aa_data']['molar_mass_residue']) *
    proteins['gelatin_type_b_bovine']['mol% (norm)'] * 1000)

################################################################
