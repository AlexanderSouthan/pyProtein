# -*- coding: utf-8 -*-
"""
Does basic calculations on polyampholytes such as polypeptides, e.g. the
isoelectric point similar to the "ExPASy Compute pI/Mw tool".
"""

import numpy as np
import pandas as pd
from scipy.optimize import brentq
import matplotlib.pyplot as plt

# import from own files
import group_properties
#############

class polyampholyte:
    def __init__(self, mode, **kwargs):
        """
        Parameters
        ----------
        mode : str
            Mode defining which input parameters are required for
            initialize_dataset. Currently only 'protein'.
        **kwargs : datatypes depend on mode
            abundance : list of float
                For mode 'protein', contains abundance of amino acids
                and C- and N-terminus in the order ['D', 'N', 'T', 'S', 'E',
                'Q', 'G', 'A', 'C', 'V', 'M', 'I', 'L', 'Y', 'F', 'H', 'K',
                'R', 'P', 'W', 'N_term', 'C_term']

        Returns
        -------
        None.

        """
        self.mode = mode
        self.kwargs = kwargs
        self.initialize_dataset()

    def initialize_dataset(self):
        """
        Should transform input data into DataFrame containing amino acid
        abundances, pKa values and charge_indicator of relevant groups.
        Other modes for example based on the amino acid sequence need to
        be added.

        Currently the only mode is 'protein' with options 'abundance' and
        'sequence'.
        """

        if self.mode == 'protein':
            self.dataset = group_properties.amino_acids.copy()
            
            if 'abundance' in self.kwargs:
                abundance = self.kwargs.get('abundance')
                # if less abundance values than entries in amino acid table are
                # given, remaining abundances are set to zero
                abundance.extend(
                        (len(self.dataset.index) - len(abundance)) * [0])
                self.dataset['abundance_input'] = abundance
            elif 'sequence' in self.kwargs:
                sequence = self.kwargs.get('sequence')
                abundance = []
                # occurences of the first 20 amino acids in the sequence are
                # counted, needs to be adapted if more than 20 amino acids such
                # as hydroxylysine are added
                for key in self.dataset.index[:20]:
                    abundance.append(sequence.count(key))
                # N_term and C_term abundances
                abundance.extend([1, 1])
                self.dataset['abundance_input'] = abundance
            else:
                raise ValueError(
                        'Unknown or missing input for protein composition.')

            # normalize abundances to sum to make different inputs comparable
            self.dataset['abundance_norm'] = (
                    self.dataset['abundance_input']/
                    self.dataset['abundance_input'].sum())

            # look for user input for pKa data
            # else default to 'pka_bjellqvist'
            if 'pka_data' in self.kwargs:
                self.pka_data = self.kwargs.get('pka_data')
            else:
                self.pka_data = 'pka_bjellqvist'

            # identify IEP relevant rows in dataset
            # and select for IEP calculations
            data_mask = ~self.dataset[self.pka_data].isna().values
            self.IEP_dataset = self.dataset.iloc[data_mask]
        else:
            raise ValueError('Unknown mode for IEP calculation.')

    def calc_charge(self, pH):
        """
        Calculates the net charge of the polyampholyte at a given pH.

        Parameters
        ----------
        pH : float
            The pH used for net charge calculation.

        Returns
        -------
        charge : float
            The net charge of the polyampholyte at the given pH.

        """
        charge = np.sum(self.IEP_dataset['charge_indicator'].values *
                        self.IEP_dataset['abundance_norm'].values /
                        (1+10**(self.IEP_dataset['charge_indicator'].values *
                                (pH-self.IEP_dataset[
                                        self.pka_data].values))))
        return charge

    def calc_charge_curve(self, ph_range=[0, 14], data_points=100):
        """
        Calculates the charge curve of the polyampholyte in a given
        pH range.

        Parameters
        ----------
        ph_range : list of floats, optional
            First value is the lower pH limit and second value the upper pH
            limit used for calculations. The default is [0, 14].
        data_points : int, optional
            Number of data points calculated. The default is 100.

        Returns
        -------
        curve : ndarray
            2D array with shape (2,data_points) containing the pH values
            in the first row and the net charges in the second row.

        """
        pH = np.linspace(ph_range[0], ph_range[1], data_points)

        curve = np.sum(self.IEP_dataset['charge_indicator'].values *
                       self.IEP_dataset['abundance_norm'].values /
                       (1+10**(self.IEP_dataset['charge_indicator'].values *
                               (pH[:, np.newaxis] -
                                self.IEP_dataset[self.pka_data].values
                                ))), axis=1)
        return np.array([pH, curve])

    def calc_IEP(self, ph_range=[0, 14]):
        """
        Calculates the isoelectric point (IEP) of the polyampholyte, i.e. the
        pH value with a net charge of zero.

        Parameters
        ----------
        ph_range : list of floats, optional
            First value is the lower pH limit and second value the upper pH
            limit used for calculations. The default is [0, 14].

        Returns
        -------
        IEP : float or np.nan
            The calculated IEP is returned as float. If no IEP was found np.nan
            is returned.

        """
        try:
            IEP = brentq(self.calc_charge, ph_range[0], ph_range[1])
            return IEP
        except:
            return np.nan


if __name__ == "__main__":
    # demostrates the function by IEP calculation of gelatin type A, gelatin
    # type B, bovine serum albumin
    
    # amino acid abundances for gelatin type A and type B as well as
    # experimentally determined IEPs taken from
    # Sewald et al., Macromol. Biosci. 2018, 18 (12), 1800168.
    # DOI: 10.1002/mabi.201800168.
    gelatin_type_a = polyampholyte(
        'protein', abundance=[0.286, 0.158, 0.144, 0.328, 0.469, 0.245, 3.314,
                              1.037, 0, 0.206, 0.06, 0.094, 0.223, 0.033,
                              0.126, 0.048, 0.325, 0.483, 2.082, 0])

    gelatin_type_b = polyampholyte(
        'protein', abundance=[0.427, 0, 0.154, 0.306, 0.701, 0, 3.248, 1.104,
                              0, 0.201, 0.05, 0.111, 0.238, 0.01, 0.116, 0.038,
                              0.314, 0.462, 2.046, 0, 0, 0])
    
    # BSA sequence from https://www.uniprot.org/uniprot/P02769, experimentally
    # determined IEP from Salis et al., Langmuir 2011, 27 (18), 11597-11604.
    # DOI: 10.1021/la2024605. (Actually values were reported there between 4.7
    # and 5.6, so an intermediate value of 5.15 was used). 
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
    
    bovine_serum_albumin = polyampholyte(
            'protein',sequence=bovine_serum_albumin_mature,
            pka_data='pka_bjellqvist')

    pH_range = [0, 14]
    
    IEPs = pd.DataFrame([],index=['pka_bjellqvist', 'pka_ipc_protein',
                        'pka_emboss','experimental'],columns=['Gelatin type A',
                                    'Gelatin type B', 'Bovine serum albumin'])
    
    IEPs.at['experimental','Gelatin type A'] = 8.8
    IEPs.at['experimental','Gelatin type B'] = 4.9
    IEPs.at['experimental','Bovine serum albumin'] = 5.15
    
    IEPs.at['pka_bjellqvist', 'Gelatin type A'] = round(
            gelatin_type_a.calc_IEP(ph_range=pH_range), 2)
    IEPs.at['pka_bjellqvist', 'Gelatin type B'] = round(
            gelatin_type_b.calc_IEP(ph_range=pH_range), 2)
    IEPs.at['pka_bjellqvist', 'Bovine serum albumin'] = round(
            bovine_serum_albumin.calc_IEP(ph_range=pH_range), 2)
    gelatin_type_a.pka_data = 'pka_emboss'
    gelatin_type_b.pka_data = 'pka_emboss'
    bovine_serum_albumin.pka_data = 'pka_emboss'
    IEPs.at['pka_emboss', 'Gelatin type A'] = round(
            gelatin_type_a.calc_IEP(ph_range=pH_range), 2)
    IEPs.at['pka_emboss', 'Gelatin type B'] = round(
            gelatin_type_b.calc_IEP(ph_range=pH_range), 2)
    IEPs.at['pka_emboss', 'Bovine serum albumin'] = round(
            bovine_serum_albumin.calc_IEP(ph_range=pH_range), 2)
    gelatin_type_a.pka_data = 'pka_ipc_protein'
    gelatin_type_b.pka_data = 'pka_ipc_protein'
    bovine_serum_albumin.pka_data = 'pka_ipc_protein'
    IEPs.at['pka_ipc_protein', 'Gelatin type A'] = round(
            gelatin_type_a.calc_IEP(ph_range=pH_range), 2)
    IEPs.at['pka_ipc_protein', 'Gelatin type B'] = round(
            gelatin_type_b.calc_IEP(ph_range=pH_range), 2)
    IEPs.at['pka_ipc_protein', 'Bovine serum albumin'] = round(
            bovine_serum_albumin.calc_IEP(ph_range=pH_range), 2)

    # calculate positions of bars grouped by columns of IEPs
    positions = []
    pos_counter = 1
    for ii in range(len(IEPs.columns)):
        for jj in range(len(IEPs.index)):
            positions.append(pos_counter)
            pos_counter += 1
        pos_counter += 1

    plt.figure(0)
    plt.bar(positions, IEPs.values.T.flatten(), ec='k', fc='skyblue')
    plt.ylabel('IEP')

    tick_positions = np.mean(np.array(positions).reshape(
            len(IEPs.columns), len(IEPs.index)), axis=1)
    plt.xticks(tick_positions, IEPs.columns)

    # write used pKa values as text into bars
    for kk, curr_pos in enumerate(positions):
        plt.text(curr_pos, 0.5, IEPs.index[
                kk % len(IEPs.index)],{'ha': 'center', 'va': 'bottom'},
        rotation=90)
        
    vline_positions = tick_positions[0:-1] + np.diff(tick_positions)/2
    plt.vlines(vline_positions,plt.ylim()[0],plt.ylim()[1], linestyles='--',
               linewidths=1)

    charge_curve_typeA = gelatin_type_a.calc_charge_curve(ph_range=pH_range)
    charge_curve_typeB = gelatin_type_b.calc_charge_curve(ph_range=pH_range)
    charge_curve_BSA = bovine_serum_albumin.calc_charge_curve(
            ph_range=pH_range)

    plt.figure(1)
    plt.plot(charge_curve_typeA[0], charge_curve_typeA[1])
    plt.plot(charge_curve_typeB[0], charge_curve_typeB[1])
    plt.plot(charge_curve_BSA[0], charge_curve_BSA[1])
    plt.hlines(0, pH_range[0], pH_range[1], linestyles='--',
                            linewidths=1)
    plt.xlabel('pH')
    plt.ylabel('relative charge')
    plt.xlim((pH_range[0], pH_range[1]))
    plt.legend(['Gelatin type A', 'Gelatin type B', 'Bovine serum albumin'])
