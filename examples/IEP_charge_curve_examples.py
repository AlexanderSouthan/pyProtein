# -*- coding: utf-8 -*-
"""
Demostrates the function by IEP calculation of gelatin type A, gelatin
type B, bovine serum albumin

Amino acid abundances for gelatin type A and type B as well as
experimentally determined results taken from
Sewald et al., Macromol. Biosci. 2018, 18 (12), 1800168.
DOI: 10.1002/mabi.201800168.

BSA sequence from https://www.uniprot.org/uniprot/P02769, experimentally
determined IEP from Salis et al., Langmuir 2011, 27 (18), 11597-11604.
DOI: 10.1021/la2024605. (Actually values were reported there between 4.7
and 5.6, so an intermediate value of 5.15 was used).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import pyProtein
from protein_compositions import proteins

# Define gelatin type A and gelatine type B via abundance of amino acids given
# in mmol/g. Could be any other unit, as long as the relative abundances are
# correct.
gelatin_type_a = pyProtein.protein(
    'mmol_g', proteins['gelatin_type_a_porcine']['mmol/g (norm, per residue)'].tolist())
gelatin_type_b = pyProtein.protein(
    'mmol_g', proteins['gelatin_type_b_bovine']['mmol/g'].tolist())

# GbM2 = pyProtein.protein(
#     'mmol_g', mmol_g=proteins['gelatin_type_b_bovine']['mmol/g'].tolist(),
#     mod_types=['Methacryl'], mod_abundances=[0.29],
#     mod_sites=['K'])

# GbM2A8 = pyProtein.protein(
#     'mmol_g', mmol_g=proteins['gelatin_type_b_bovine']['mmol/g'].tolist(),
#     mod_types=['Methacryl', 'Acetyl'], mod_abundances=[0.36, 0.418],
#     mod_sites=['K', 'K'])

# GbM10 = pyProtein.protein(
#     'mmol_g', mmol_g=proteins['gelatin_type_b_bovine']['mmol/g'].tolist(),
#     mod_types=['Methacryl'], mod_abundances=[0.95],
#     mod_sites=['K'])

# Define bovine serum albumin via its amino acid sequence.
# DTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA
bovine_serum_albumin_mature = (
        'DTHKSEIAHRFKDLGEEHFKGLVLI'
        'AFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDEL'
        'CKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLC'
        'DEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAED'
        'KGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQK'
        'FPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTI'
        'SSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQ'
        'EAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHAC'
        'YSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQV'
        'STPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKT'
        'PVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTL'
        'PDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKE'
        'ACFAVEGPKLVVSTQTALA')
bovine_serum_albumin = pyProtein.protein(
        'sequence', bovine_serum_albumin_mature,
        pka_data='pka_bjellqvist', mod_types=['N_term', 'C_term'],
        mod_abundances=[1, 1], pka_scales=['pka_bjellqvist', 'pka_bjellqvist'])
# bovine_serum_albumin_mod = pyProtein.protein(
#         'sequence', sequence=bovine_serum_albumin_mature,
#         pka_data='pka_bjellqvist', mod_types=['N_term', 'C_term', 'Methacryl'],
#         mod_abundances=[1, 1, 10], mod_sites=['any', 'any', 'K'])

# Define the pH range used for calculations. First number is the lower limit,
# second number is the upper limit.
pH_range = [0, 14]

# Prepare a DataFrame that will hold the results of calculations.
results = pd.DataFrame([], index=['IEP_pka_bjellqvist', 'IEP_pka_ipc_protein',
                    'IEP_pka_emboss', 'IEP_experimental',
                    'avg_residue_molar_mass', 'total_N'], columns=[
                        'Gelatin type A', 'Gelatin type B',
                        'Bovine serum albumin'])

# Calculate average amino acid residue molar masses and store them in results
# DataFrame
results.at['avg_residue_molar_mass', 'Gelatin type A'] = (
    gelatin_type_a.mean_residue_molar_mass())
results.at['avg_residue_molar_mass', 'Gelatin type B'] = (
    gelatin_type_b.mean_residue_molar_mass())
results.at['avg_residue_molar_mass', 'Bovine serum albumin'] = (
    bovine_serum_albumin.mean_residue_molar_mass())

# Calculate total N content and store them in results DataFrame
results.at['total_N', 'Gelatin type A'] = (
    gelatin_type_a.n_content())
results.at['total_N', 'Gelatin type B'] = (
    gelatin_type_b.n_content())
results.at['total_N', 'Bovine serum albumin'] = (
    bovine_serum_albumin.n_content())

# Put experimental IEP values into results DataFrame.
results.at['IEP_experimental', 'Gelatin type A'] = 8.8
results.at['IEP_experimental', 'Gelatin type B'] = 4.9
results.at['IEP_experimental', 'Bovine serum albumin'] = 5.15

# Calculate the results of the three proteins and write into results DataFrame.
# Note that in between the pka_data property is changed to do the calculations
# based on the different pKa value tables.
results.at['IEP_pka_bjellqvist', 'Gelatin type A'] = round(
        gelatin_type_a.IEP(ph_range=pH_range), 2)
results.at['IEP_pka_bjellqvist', 'Gelatin type B'] = round(
        gelatin_type_b.IEP(ph_range=pH_range), 2)
results.at['IEP_pka_bjellqvist', 'Bovine serum albumin'] = round(
        bovine_serum_albumin.IEP(ph_range=pH_range), 2)
gelatin_type_a.pka_data = 'pka_emboss'
gelatin_type_a.initialize_pka_dataset()
gelatin_type_b.pka_data = 'pka_emboss'
gelatin_type_b.initialize_pka_dataset()
bovine_serum_albumin.pka_data = 'pka_emboss'
bovine_serum_albumin.initialize_pka_dataset()
results.at['IEP_pka_emboss', 'Gelatin type A'] = round(
        gelatin_type_a.IEP(ph_range=pH_range), 2)
results.at['IEP_pka_emboss', 'Gelatin type B'] = round(
        gelatin_type_b.IEP(ph_range=pH_range), 2)
results.at['IEP_pka_emboss', 'Bovine serum albumin'] = round(
        bovine_serum_albumin.IEP(ph_range=pH_range), 2)
gelatin_type_a.pka_data = 'pka_ipc_protein'
gelatin_type_a.initialize_pka_dataset()
gelatin_type_b.pka_data = 'pka_ipc_protein'
gelatin_type_b.initialize_pka_dataset()
bovine_serum_albumin.pka_data = 'pka_ipc_protein'
bovine_serum_albumin.initialize_pka_dataset()
results.at['IEP_pka_ipc_protein', 'Gelatin type A'] = round(
        gelatin_type_a.IEP(ph_range=pH_range), 2)
results.at['IEP_pka_ipc_protein', 'Gelatin type B'] = round(
        gelatin_type_b.IEP(ph_range=pH_range), 2)
results.at['IEP_pka_ipc_protein', 'Bovine serum albumin'] = round(
        bovine_serum_albumin.IEP(ph_range=pH_range), 2)

# The rest is just for plotting the results

# Calculate positions of bars grouped by columns of results
positions = []
pos_counter = 1
for ii in range(len(results.columns)):
    for jj in range(len(results.index)-2):
        positions.append(pos_counter)
        pos_counter += 1
    pos_counter += 1

plt.figure(0)
plt.bar(positions, results.values[:-2].T.flatten(), ec='k', fc='skyblue')
plt.ylabel('IEP')

tick_positions = np.mean(np.array(positions).reshape(
        len(results.columns), len(results.index)-2), axis=1)
plt.xticks(tick_positions, results.columns)

# Write used pKa values as text into bars
for kk, curr_pos in enumerate(positions):
    plt.text(curr_pos, 0.5, results.index[
            kk % (len(results.index)-2)], {'ha': 'center', 'va': 'bottom'},
        rotation=90)

vline_positions = tick_positions[0:-1] + np.diff(tick_positions)/2
plt.vlines(vline_positions, plt.ylim()[0], plt.ylim()[1], linestyles='--',
            linewidths=1)

charge_curve_typeA = gelatin_type_a.charge_curve(ph_range=pH_range)
charge_curve_typeB = gelatin_type_b.charge_curve(ph_range=pH_range)
charge_curve_BSA = bovine_serum_albumin.charge_curve(
        ph_range=pH_range)

plt.figure(1)
plt.plot(charge_curve_typeA[0], charge_curve_typeA[1])
plt.plot(charge_curve_typeB[0], charge_curve_typeB[1])
plt.plot(charge_curve_BSA[0], charge_curve_BSA[1])
plt.hlines(0, pH_range[0], pH_range[1], linestyles='--', linewidths=1)
plt.xlabel('pH')
plt.ylabel('relative charge')
plt.xlim((pH_range[0], pH_range[1]))
plt.legend(['Gelatin type A', 'Gelatin type B', 'Bovine serum albumin'])
