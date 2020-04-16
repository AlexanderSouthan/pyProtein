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

import pyPolyampholyte

# Define gelatin type A and gelatine type B via abundance of amino acids given
# in mmol/g. Could be any other unit, as long as the relative abundances are
# correct.
gelatin_type_a = pyPolyampholyte.polyampholyte(
    'protein', abundance=[0.286, 0.158, 0.144, 0.328, 0.469, 0.245, 3.314,
                          1.037, 0, 0.206, 0.06, 0.094, 0.223, 0.033,
                          0.126, 0.048, 0.259, 0.483, 1.180, 0, 0.902, 0.066])
gelatin_type_b = pyPolyampholyte.polyampholyte(
    'protein', abundance=[0.427, 0, 0.154, 0.306, 0.701, 0, 3.248, 1.104,
                          0, 0.201, 0.05, 0.111, 0.238, 0.01, 0.116, 0.038,
                          0.253, 0.462, 1.075, 0, 0.971, 0.061])

# Define bovine serum albumin via its amino acid sequence.
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
bovine_serum_albumin = pyPolyampholyte.polyampholyte(
        'protein', sequence=bovine_serum_albumin_mature,
        pka_data='pka_bjellqvist')

# Define the pH range used for calculations. First number is the lower limit,
# second number is the upper limit.
pH_range = [0, 14]

# Prepare a DataFrame that will hold the results of calculations.
results = pd.DataFrame([], index=['IEP_pka_bjellqvist', 'IEP_pka_ipc_protein',
                    'IEP_pka_emboss', 'IEP_experimental',
                    'avg_residue_molar_mass'], columns=[
                        'Gelatin type A', 'Gelatin type B',
                        'Bovine serum albumin'])

# Calculate average amino acid residue molar masses and store them in results
# DataFrame
results.at['avg_residue_molar_mass', 'Gelatin type A'] = (
    gelatin_type_a.calc_mean_residue_mw())
results.at['avg_residue_molar_mass', 'Gelatin type B'] = (
    gelatin_type_b.calc_mean_residue_mw())
results.at['avg_residue_molar_mass', 'Bovine serum albumin'] = (
    bovine_serum_albumin.calc_mean_residue_mw())
                        
# Put experimental IEP values into results DataFrame.
results.at['IEP_experimental', 'Gelatin type A'] = 8.8
results.at['IEP_experimental', 'Gelatin type B'] = 4.9
results.at['IEP_experimental', 'Bovine serum albumin'] = 5.15

# Calculate the results of the three proteins and write into results DataFrame.
# Note that in between the pka_data property is changed to do the calculations
# based on the different pKa value tables.
results.at['IEP_pka_bjellqvist', 'Gelatin type A'] = round(
        gelatin_type_a.calc_IEP(ph_range=pH_range), 2)
results.at['IEP_pka_bjellqvist', 'Gelatin type B'] = round(
        gelatin_type_b.calc_IEP(ph_range=pH_range), 2)
results.at['IEP_pka_bjellqvist', 'Bovine serum albumin'] = round(
        bovine_serum_albumin.calc_IEP(ph_range=pH_range), 2)
gelatin_type_a.pka_data = 'pka_emboss'
gelatin_type_b.pka_data = 'pka_emboss'
bovine_serum_albumin.pka_data = 'pka_emboss'
results.at['IEP_pka_emboss', 'Gelatin type A'] = round(
        gelatin_type_a.calc_IEP(ph_range=pH_range), 2)
results.at['IEP_pka_emboss', 'Gelatin type B'] = round(
        gelatin_type_b.calc_IEP(ph_range=pH_range), 2)
results.at['IEP_pka_emboss', 'Bovine serum albumin'] = round(
        bovine_serum_albumin.calc_IEP(ph_range=pH_range), 2)
gelatin_type_a.pka_data = 'pka_ipc_protein'
gelatin_type_b.pka_data = 'pka_ipc_protein'
bovine_serum_albumin.pka_data = 'pka_ipc_protein'
results.at['IEP_pka_ipc_protein', 'Gelatin type A'] = round(
        gelatin_type_a.calc_IEP(ph_range=pH_range), 2)
results.at['IEP_pka_ipc_protein', 'Gelatin type B'] = round(
        gelatin_type_b.calc_IEP(ph_range=pH_range), 2)
results.at['IEP_pka_ipc_protein', 'Bovine serum albumin'] = round(
        bovine_serum_albumin.calc_IEP(ph_range=pH_range), 2)

# The rest is just for plotting the results

# Calculate positions of bars grouped by columns of results
positions = []
pos_counter = 1
for ii in range(len(results.columns)):
    for jj in range(len(results.index)):
        positions.append(pos_counter)
        pos_counter += 1
    pos_counter += 1

plt.figure(0)
plt.bar(positions, results.values.T.flatten(), ec='k', fc='skyblue')
plt.ylabel('IEP')

tick_positions = np.mean(np.array(positions).reshape(
        len(results.columns), len(results.index)), axis=1)
plt.xticks(tick_positions, results.columns)

# Write used pKa values as text into bars
for kk, curr_pos in enumerate(positions):
    plt.text(curr_pos, 0.5, results.index[
            kk % len(results.index)], {'ha': 'center', 'va': 'bottom'},
        rotation=90)

vline_positions = tick_positions[0:-1] + np.diff(tick_positions)/2
plt.vlines(vline_positions, plt.ylim()[0], plt.ylim()[1], linestyles='--',
           linewidths=1)

charge_curve_typeA = gelatin_type_a.calc_charge_curve(ph_range=pH_range)
charge_curve_typeB = gelatin_type_b.calc_charge_curve(ph_range=pH_range)
charge_curve_BSA = bovine_serum_albumin.calc_charge_curve(
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
