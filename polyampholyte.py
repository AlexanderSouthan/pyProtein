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
            abundance : ndarray
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
            self.dataset = group_properties.amino_acids
            
            if 'abundance' in self.kwargs:
                self.dataset['abundance_input'] = self.kwargs.get('abundance')
            elif 'sequence' in self.kwargs:
                sequence = self.kwargs.get('sequence')
                abundance = []
                for key in self.dataset.index:
                    abundance.append(sequence.count(key))
                self.dataset['abundance_input'] = abundance
            else:
                raise ValueError(
                        'Unknown or missing input for protein composition.')
            
            self.dataset['abundance_norm'] = ( #  normalize to sum to make different inputs comparable
                    self.dataset['abundance_input']/self.dataset['abundance_input'].sum())
            if 'pKa_data' in self.kwargs:
                self.pKa_data = self.kwargs.get('pKa_data')
            else: #  defaults to the following value
                self.pKa_data = 'pKa_bjellqvist'

            data_mask = ~self.dataset[self.pKa_data].isna().values
            self.IEP_dataset = self.dataset.iloc[data_mask]
        else:
            raise ValueError('Unknown mode for IEP calculation')

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
                                        self.pKa_data].values))))
        return charge

    def calc_charge_curve(self, range=[0, 14], data_points=100):
        """
        Calculates the charge curve of the polyampholyte in a given
        pH range.

        Parameters
        ----------
        range : list of floats, optional
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
        pH = np.linspace(range[0], range[1], data_points)

        curve = np.sum(self.IEP_dataset['charge_indicator'].values *
                       self.IEP_dataset['abundance_norm'].values /
                       (1+10**(self.IEP_dataset['charge_indicator'].values *
                               (pH[:, np.newaxis] -
                                self.IEP_dataset[self.pKa_data].values
                                ))), axis=1)
        return np.array([pH, curve])

    def calc_IEP(self, range=[0, 14]):
        """
        Calculates the isoelectric point (IEP) of the polyampholyte, i.e. the
        pH value with a net charge of zero.

        Parameters
        ----------
        range : list of floats, optional
            First value is the lower pH limit and second value the upper pH
            limit used for calculations. The default is [0, 14].

        Returns
        -------
        IEP : float or np.nan
            The calculated IEP is returned as float. If no IEP was found np.nan
            is returned.

        """
        try:
            IEP = brentq(self.calc_charge, range[0], range[1])
            return IEP
        except:
            return np.nan


if __name__ == "__main__":
    # demostrates the function by IEP calculation of gelatin type A, gelatin
    # type B, bovine serum albumin
    gelatin_type_a = polyampholyte(
        'protein', abundance=np.array(
            [0.286, 0.158, 0.144, 0.328, 0.469, 0.245, 3.314, 1.037, 0, 0.206,
             0.06, 0.094, 0.223, 0.033, 0.126, 0.048, 0.325, 0.483, 2.082, 0,
             0, 0]))

    gelatin_type_b = polyampholyte(
        'protein', abundance=np.array(
            [0.427, 0, 0.154, 0.306, 0.701, 0, 3.248, 1.104, 0, 0.201, 0.05,
             0.111, 0.238, 0.01, 0.116, 0.038, 0.314, 0.462, 2.046, 0, 0, 0]))
    
    # BSA sequence from https://www.uniprot.org/uniprot/P02769
    bovine_serum_albumin_sequence = (
            'MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLI'
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
            'protein',sequence=bovine_serum_albumin_sequence,
            pKa_data='pKa_IPC_protein')

    pH_range = [0, 14]
    
    print('pH range used for calculations:', pH_range)
    print('IEP gelatin type A:',
          gelatin_type_a.calc_IEP(range=pH_range))
    print('IEP gelatin type B:',
          gelatin_type_b.calc_IEP(range=pH_range))
    print('IEP bovine serum albumin:',
          bovine_serum_albumin.calc_IEP(range=pH_range))

    charge_curve_typeA = gelatin_type_a.calc_charge_curve(range=pH_range)
    charge_curve_typeB = gelatin_type_b.calc_charge_curve(range=pH_range)
    charge_curve_BSA = bovine_serum_albumin.calc_charge_curve(range=pH_range)

    plt.plot(charge_curve_typeA[0], charge_curve_typeA[1])
    plt.plot(charge_curve_typeB[0], charge_curve_typeB[1])
    plt.plot(charge_curve_BSA[0], charge_curve_BSA[1])
    plt.hlines(0, pH_range[0], pH_range[1], linestyles='--',
                            linewidths=1)
    plt.xlabel('pH')
    plt.ylabel('relative charge')
    plt.xlim((pH_range[0], pH_range[1]))
    plt.legend(['Gelatin type A', 'Gelatin type B'])
