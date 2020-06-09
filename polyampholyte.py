# -*- coding: utf-8 -*-
"""
Does basic calculations on polyampholytes such as polypeptides, e.g. the
isoelectric point similar to the "ExPASy Compute pI/Mw tool".
"""

import numpy as np
import pandas as pd
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
            mmol_g : list of float
                For mode 'protein', contains abundance of amino acids
                and C- and N-terminus in the order ['D', 'N', 'T', 'S', 'E',
                'Q', 'G', 'A', 'C', 'V', 'M', 'I', 'L', 'Y', 'F', 'H', 'K',
                'R', 'P', 'W', 'Hyp'm 'Hyl']
            sequence : str
                For mode 'protein', alternative to mmol_g. Has to be a
                string consisting of one letter codes for amino acids.
            absolute : list of float
                For mode 'protein', alternative to mmol_g and sequence. Is a
                list in the same order as for mmol_g, but with absolute
                abundances of amino acids, e.g. deducted from the sequence.
            res_per_1000 : list of float
                For mode 'protein', list in the same order as mmol_g, but
                abundances normed to 1000 amino acid residues.

        Returns
        -------
        None.

        """
        self.mode = mode
        self.kwargs = kwargs
        self.initialize_main_chain()
        self.initialize_modifications()
        self.initialize_pka_dataset()

        self.mass_calc_dataset = pd.concat(
            [self.main_chain, self.modifications])[
                ['molar_mass_residue', 'N_content_residue', 'abundance_input',
                 'abundance_norm']]

        # the following calculation gives mass percents on residue basis
        self.mass_calc_dataset['mass_percent_residue'] = (
            self.mass_calc_dataset['abundance_input'] *
            self.mass_calc_dataset['molar_mass_residue'] /
            np.sum(
                self.mass_calc_dataset['abundance_input'] *
                self.mass_calc_dataset['molar_mass_residue'])
            )

    def initialize_main_chain(self):
        """
        Transforms input data into DataFrame containing amino acid
        abundances, pKa values and charge_indicator of relevant groups.

        Currently the only mode is 'protein' with options 'abundance',
        'absolute', and 'sequence'.
        """
        if self.mode == 'protein':
            self.main_chain = group_properties.amino_acids.copy()
            self.amino_acid_sequence = None
            self.number_of_residues = None

            if 'mmol_g' in self.kwargs:
                self.abundance_unit = 'mmol_g'
                abundance = self.kwargs.get('mmol_g')
                # if less abundance values than entries in amino acid table are
                # given, remaining abundances are set to zero
                abundance.extend(
                        (len(self.main_chain.index) - len(abundance)) * [0])
                self.main_chain['abundance_input'] = abundance
            elif 'absolute' in self.kwargs:
                self.abundance_unit = 'abs'
                abundance = self.kwargs.get('absolute')
                # if less abundance values than entries in amino acid table are
                # given, remaining abundances are set to zero
                abundance.extend(
                        (len(self.main_chain.index) - len(abundance)) * [0])
                self.main_chain['abundance_input'] = abundance
            elif 'res_per_1000' in self.kwargs:
                self.abundance_unit = 'res_per_1000'
                abundance = self.kwargs.get('res_per_1000')
                # if less abundance values than entries in amino acid table are
                # given, remaining abundances are set to zero
                abundance.extend(
                        (len(self.main_chain.index) - len(abundance)) * [0])
                self.main_chain['abundance_input'] = abundance
            elif 'sequence' in self.kwargs:
                self.abundance_unit = 'abs'
                self.amino_acid_sequence = self.kwargs.get('sequence')
                abundance = []
                # occurences of the first 20 amino acids in the sequence are
                # counted, needs to be adapted if more than 20 amino acids such
                # as hydroxylysine are added
                self.number_of_residues = len(self.amino_acid_sequence)
                for key in self.main_chain.index[:20]:
                    abundance.append(self.amino_acid_sequence.count(key))
                # Hyp, Hyl abundances
                abundance.extend(
                    (len(self.main_chain.index) - len(abundance)) * [0])
                self.main_chain['abundance_input'] = abundance
            else:
                raise ValueError(
                        'Unknown or missing input for protein composition.')

            # Normalize abundances to sum to make different inputs comparable.
            self.main_chain['abundance_norm'] = (
                self.main_chain['abundance_input'] /
                np.sum(self.main_chain['abundance_input'].values))

        else:
            raise ValueError('Unknown mode for IEP calculation.')

    def initialize_modifications(self):
        types = self.kwargs.get('mod_types', [])
        abundances = self.kwargs.get('mod_abundances', [])
        sites = self.kwargs.get('mod_sites', ['any']*len(types))

        self.modifications = group_properties.chain_modifications.copy()

        self.modifications['abundance_input'] = np.nan
        for mod_type, abundance, site in zip(types, abundances, sites):
            self.modifications.at[mod_type, 'abundance_input'] = abundance
            self.modifications.at[mod_type, 'modified_residues'] = site
            
        self.modifications.dropna(subset=['abundance_input'], inplace=True)

        # Absolute modification abundances are transformed to abundance
        # per amino acid residue
        self.modifications['abundance_norm'] = (
            self.modifications['abundance_input'].values /
            np.sum(self.main_chain['abundance_input'].values))

        if self.abundance_unit == 'mmol_g':
            main_chain_fraction = 1 - np.sum(
                self.modifications['abundance_input']/1000 *
                self.modifications['molar_mass_residue'])

            self.main_chain['abundance_input'] *= main_chain_fraction

    def initialize_pka_dataset(self):
        # This method needs heavy updating. It must collect relevant entries from
        # self.modifications and self.main_chain and recalculate the amino
        # acid abundances of relevant amino acids. Then, alle other functions
        # dealing with charges need updating. All other functions dealing with
        # masses are already up to date and work also with modifications.

        # look for user input for pKa data
        # else default to 'pka_bjellqvist'
        if 'pka_data' in self.kwargs:
            self.pka_data = self.kwargs.get('pka_data')
        else:
            self.pka_data = 'pka_bjellqvist'

        # identify IEP relevant rows in dataset
        # and select for IEP calculations
        data_mask = ~self.main_chain[self.pka_data].isna().values
        self.IEP_dataset = self.main_chain.iloc[data_mask]

    def charge(self, pH):
        """
        Calculate the net charge of the polyampholyte at a given pH.

        Parameters
        ----------
        pH : float
            The pH used for net charge calculation.

        Returns
        -------
        charge : float
            The net charge of the polyampholyte at the given pH.

        """
        
        # There might be a problem with the normalized abundances, because
        # C-term and N_term abundances are not included in the denominator
        # of normalization. Therefore, the sum of normalized abundances of all
        # amino acids and the termini is a little larger than one. Must be
        # checked in the future, also relevant for charge_curve. Probably
        # especially important for small peptides. Is irrelevant if no amino
        # acid sequence is given because then the termini are and must be
        # ignored anyway in the calculations.
        charge = np.sum(self.IEP_dataset['charge_indicator'].values *
                        self.IEP_dataset['abundance_norm'].values /
                        (1+10**(self.IEP_dataset['charge_indicator'].values *
                                (pH-self.IEP_dataset[
                                        self.pka_data].values))))
        return charge

    def charge_curve(self, ph_range=[0, 14], data_points=100):
        """
        Calculate the charge curve of the polyampholyte in a given
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

    def IEP(self, ph_range=[0, 14]):
        """
        Calculate the isoelectric point (IEP) of the polyampholyte, i.e. the
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
            IEP = brentq(self.charge, ph_range[0], ph_range[1])
        except ValueError:
            IEP = np.nan
        return IEP

    def molar_mass(self):
        # assert self.abundance_unit == 'abs', 'Amino acid abundances are not absolute values.'

        if self.abundance_unit == 'abs':
            molar_mass = np.sum(
                self.mass_calc_dataset['abundance_input'].values *
                self.mass_calc_dataset['molar_mass_residue'].values)
            return molar_mass
        else:
            return None

    def mean_residue_molar_mass(self):
        """
        Calculate the mean residue molecular weight of proteins.

        Calculation is done without taking C and N terminus into account if
        only relative abundances are known and not the entire sequence.

        Returns
        -------
        mean_residue_mw : float
            The mean residue molecular weight based on the abundances of amino
            acids.

        """
        # Calculate average molar masses without formally condensed water
        # during protein formation for each amino acid.
        mean_residue_molar_mass = np.sum(
            self.mass_calc_dataset['molar_mass_residue'].values *
            self.mass_calc_dataset['abundance_norm'].values)

        return mean_residue_molar_mass

    def n_content(self):
        n_content = (
            np.sum(
                self.mass_calc_dataset['N_content_residue'] *
                self.mass_calc_dataset['abundance_norm'] *
                self.mass_calc_dataset['molar_mass_residue']) /
            np.sum(
                self.mass_calc_dataset['abundance_norm'] *
                self.mass_calc_dataset['molar_mass_residue']
                )
            )
        return n_content
