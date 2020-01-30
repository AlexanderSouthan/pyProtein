# -*- coding: utf-8 -*-
"""
Does basic calculations on polyampholytes such as polypeptides, e.g. the
isoelectric point similar to the "ExPASy Compute pI/Mw tool".
"""

import numpy as np
from scipy.optimize import brentq

from . import group_properties


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
                    self.dataset['abundance_input'] /
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
        except ValueError:
            IEP = np.nan
        return IEP
