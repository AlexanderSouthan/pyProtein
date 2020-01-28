# -*- coding: utf-8 -*-
"""
Does basic calculations on polyampholytes such as polypeptides, e.g. the
isoelectric point similar to the "ExPASy Compute pI/Mw tool".
"""

import numpy as np
import pandas as pd
from scipy.optimize import brentq
import matplotlib.pyplot as plt


class polyampholyte:
    def __init__(self, mode, **kwargs):
        """
        Parameters
        ----------
        mode : str
            Mode defining which input parameters are required for
            initialize_dataset. Currently only 'protein_residues'.
        **kwargs : datatypes depend on mode
            abundance : ndarray
                For mode 'protein_residues', contains abundance of amino acids
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

        Currently the only mode is 'protein_residues'.
        """

        if self.mode == 'protein_residues':
            # abundance = kwargs.get(abundance)

            self.dataset = pd.DataFrame(
                [], index=['D', 'N', 'T', 'S', 'E', 'Q', 'G', 'A', 'C', 'V',
                           'M', 'I', 'L', 'Y', 'F', 'H', 'K', 'R', 'P', 'W',
                           'N_term', 'C_term'])
            self.dataset['long_name'] = ['aspartic acid', 'asparagine',
                                         'threonine', 'serine',
                                         'glutamic acid', 'glutamine',
                                         'glycine', 'alanine', 'cysteine',
                                         'valine', 'methionine', 'isoleucine',
                                         'leucine', 'tyrosine',
                                         'phenylalanine', 'histidine',
                                         'lysine', 'arginine', 'proline',
                                         'tryptophan', 'N-terminus',
                                         'C-terminus']
            self.dataset['pKa'] = [4.05, np.nan, np.nan, np.nan, 4.45, np.nan,
                                   np.nan, np.nan, 9, np.nan, np.nan, np.nan,
                                   np.nan, 10, np.nan, 5.98, 10, 12, np.nan,
                                   np.nan, 7.5, 3.55]
            # pKa for glycine at N_term according to
            # doi 10.1002/elps.1150150171
            self.dataset['charge_indicator'] = [-1, 0, 0, 0, -1, 0, 0, 0, -1,
                                                0, 0, 0, 0, -1, 0, 1, 1, 1, 0,
                                                0, 1, -1]
            self.dataset['abundance'] = self.kwargs.get('abundance')
            data_mask = ~self.dataset['pKa'].isna().values
            self.active_data = self.dataset.iloc[data_mask]
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
        charge = np.sum(self.active_data['charge_indicator'].values *
                        self.active_data['abundance'].values /
                        (1+10**(self.active_data['charge_indicator'].values *
                                (pH-self.active_data['pKa'].values))))
        return charge

    def calc_titration_curve(self, range=[0, 14], data_points=100):
        """
        Calculates the titration curve of the polyampholyte in a given
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

        curve = np.sum(self.active_data['charge_indicator'].values *
                       self.active_data['abundance'].values /
                       (1+10**(self.active_data['charge_indicator'].values *
                               (pH[:, np.newaxis] -
                                self.active_data['pKa'].values))), axis=1)
        return np.array([pH, curve])

    def calc_IEP(self):
        """
        Calculates the isoelectric point (IEP) of the polyampholyte, i.e. the
        pH value with a net charge of zero.

        Returns
        -------
        IEP : float
            The calculated IEP.

        """
        IEP = brentq(self.calc_charge, 0, 14)

        return IEP


if __name__ == "__main__":
    # demostrates the function by IEP calculation of gelatin type A and
    # gelatin type B
    gelatin_type_a = polyampholyte(
        'protein_residues', abundance=np.array(
            [0.286, 0.158, 0.144, 0.328, 0.469, 0.245, 3.314, 1.037, 0, 0.206,
             0.06, 0.094, 0.223, 0.033, 0.126, 0.048, 0.325, 0.483, 2.082, 0,
             0, 0]))

    gelatin_type_b = polyampholyte(
        'protein_residues', abundance=np.array(
            [0.427, 0, 0.154, 0.306, 0.701, 0, 3.248, 1.104, 0, 0.201, 0.05,
             0.111, 0.238, 0.01, 0.116, 0.038, 0.314, 0.462, 2.046, 0, 0, 0]))

    print('IEP gelatin type A: ',
          gelatin_type_a.calc_IEP())
    print('IEP gelatin type B: ',
          gelatin_type_b.calc_IEP())

    titration_curve_typeA = gelatin_type_a.calc_titration_curve()
    titration_curve_typeB = gelatin_type_b.calc_titration_curve()

    plt.plot(titration_curve_typeA[0], titration_curve_typeA[1],
             titration_curve_typeB[0], titration_curve_typeB[1],
             titration_curve_typeA[0], np.zeros_like(titration_curve_typeA[0]))
