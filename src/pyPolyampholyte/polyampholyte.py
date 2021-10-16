# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from . import group_properties


class polyampholyte:
    """
    Does basic calculations on polyampholytes such as polypeptides, e.g.
    the isoelectric point similar to the "ExPASy Compute pI/Mw tool".
    """
    def __init__(self, mode, **kwargs):
        """
        Parameters
        ----------
        mode : str
            Mode defining which input parameters are required for
            initialize_dataset. Currently only 'protein'.
        **kwargs
            pka_data : string
                Gives the pKa dataset used for charge and IEP calculations.
                Allowed values are 'pka_ipc_protein', 'pka_emboss' or
                'pka_bjellqvist'. Default is 'pka_bjellqvist'.
        **kwargs for mode=='protein'
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
                abundances normalized to 1000 amino acid residues.
        **kwargs in case of modifications are present
            mod_types : list of strings
                Allowed values are the index values of
                group_properties.chain_modifications, defaults to empty list.
            mod_abundances : list of float
                The abundances of the modifications given by mod_types. Must be
                given in the same units as the input data for the main chain
                andcontain the same amount of elements like mod_types. Default
                is an empty list.
            mod_sites : list of strings
                The sites on the main chain which are modified, so allowed
                values are the index values in self.main_chain or 'any'. 'any'
                means that the modification site is undefined. In this case,
                the modification is only important e.g. for molar mass
                calculations, but irrelevant e.g. for IEP calculations. Must
                contain the same amount of elements like mod_types. Default is
                a list the same size as mod_types with only 'any' entries.
            pka_scales : list of strings
                The pKa sclaes used for the different modifications. Must
                contains the same amount of entires like mod_types. Default is
                a list with the same size as mod_types with only
                'pka_other' entries.

        Returns
        -------
        None.

        """
        self.mode = mode

        if self.mode == 'protein':
            self.amino_acid_sequence = None
            self.number_of_residues = None
            if 'mmol_g' in kwargs:
                self.abundance_unit = 'mmol_g'
                self.abundance = kwargs.get('mmol_g')
            elif 'absolute' in kwargs:
                self.abundance_unit = 'absolute'
                self.abundance = kwargs.get('absolute')
            elif 'res_per_1000' in kwargs:
                self.abundance_unit = 'res_per_1000'
                self.abundance = kwargs.get('res_per_1000')
            elif 'sequence' in kwargs:
                self.abundance_unit = 'sequence'
                self.amino_acid_sequence = kwargs.get('sequence')
                self.abundance = []
            else:
                raise ValueError('Mode is protein, but no valid option for'
                                 ' amino acid composition given.')

            self.mod_types = kwargs.get('mod_types', [])
            self.mod_abundances = kwargs.get('mod_abundances', [])
            self.mod_sites = kwargs.get('mod_sites',
                                        ['any']*len(self.mod_types))
            self.pka_scales = kwargs.get(
                'pka_scales', ['pka_other']*len(self.mod_types))

            # look for user input for pKa data else default to 'pka_bjellqvist'
            self.pka_data = kwargs.get('pka_data', 'pka_bjellqvist')

        else:
            raise ValueError('No valid value for mode given.')

        self.initialize_main_chain()
        self.initialize_modifications()
        # self.combine_main_chain_modifications()
        self.initialize_pka_dataset()

    def initialize_main_chain(self):
        """
        Transforms input data into DataFrame containing amino acid
        abundances, pKa values and charge_indicator of relevant groups.
        """
        if self.mode == 'protein':
            self.main_chain = pd.DataFrame(
                [], index=group_properties.amino_acids.index)  # group_properties.amino_acids.copy()
            self.main_chain['molar_mass_residue'] = (
                group_properties.amino_acids['molar_mass_residue'])
            self.main_chain['N_content_residue'] = (
                group_properties.amino_acids['N_content_residue'])

            if self.abundance_unit in ['mmol_g', 'absolute', 'res_per_1000']:
                # if less abundance values than entries in amino acid table are
                # given, remaining abundances are set to zero
                self.abundance.extend(
                        (len(self.main_chain.index) -
                         len(self.abundance)) * [0])
                self.main_chain[
                    'abundance_' + self.abundance_unit] = self.abundance
            elif self.abundance_unit == 'sequence':
                # occurences of the first 20 amino acids in the sequence are
                # counted, needs to be adapted if more than 20 amino acids such
                # as hydroxylysine are added
                self.number_of_residues = len(self.amino_acid_sequence)
                for key in self.main_chain.index[:20]:
                    self.abundance.append(self.amino_acid_sequence.count(key))
                # Hyp, Hyl abundances
                self.abundance.extend(
                    (len(self.main_chain.index) - len(self.abundance)) * [0])
                self.main_chain[
                    'abundance_' + self.abundance_unit] = self.abundance

            # Normalize abundances to sum to make different inputs comparable.
            self.main_chain['abundance_norm'] = (
                self.main_chain['abundance_' + self.abundance_unit].values /
                np.sum(self.main_chain['abundance_' +
                                        self.abundance_unit].values))

    def initialize_modifications(self):
        self.modifications = pd.DataFrame(
            [], index=group_properties.chain_modifications.index)  # group_properties.chain_modifications.copy()
        self.modifications['molar_mass_residue'] = (
            group_properties.chain_modifications['molar_mass_residue'])
        self.modifications['N_content_residue'] = (
            group_properties.chain_modifications['N_content_residue'])

        self.modifications['abundance_' + self.abundance_unit] = np.nan
        for mod_type, abundance, site, pka_scale in zip(self.mod_types,
                                                        self.mod_abundances,
                                                        self.mod_sites,
                                                        self.pka_scales):
            self.modifications.at[
                mod_type, 'abundance_' + self.abundance_unit] = abundance
            self.modifications.at[mod_type, 'modified_residues'] = site
            self.modifications.at[mod_type, 'pka_scale'] = pka_scale

        self.modifications.dropna(subset=['abundance_' + self.abundance_unit],
                                  inplace=True)

        # Modification abundances are normalized on the basis of the input data
        # of the main chain because the modifications do not introduce any
        # additional residues/repeating units
        self.modifications['abundance_norm'] = (
            self.modifications['abundance_' + self.abundance_unit].values /
            np.sum(self.main_chain['abundance_' + self.abundance_unit].values))

    def initialize_pka_dataset(self):
        self.IEP_dataset = pd.DataFrame([], index=self.main_chain.index)
        self.IEP_dataset['abundance_norm'] = self.main_chain['abundance_norm']
        self.IEP_dataset['charge_indicator'] = group_properties.amino_acids[
            'charge_indicator']
        self.IEP_dataset['pka_data'] = group_properties.amino_acids[
            self.pka_data]
        
        for idx in self.modifications.index:
            mod_idx = self.modifications.at[idx, 'modified_residues']
            mod_pka_scale = self.modifications.at[idx, 'pka_scale']
            if mod_idx != 'any':
                self.IEP_dataset.at[mod_idx, 'abundance_norm'] -= (
                    self.modifications.at[idx, 'abundance_norm'])
            self.IEP_dataset.at[idx, 'abundance_norm'] = self.modifications.at[
                idx, 'abundance_norm']
            self.IEP_dataset.at[idx, 'charge_indicator'] = (
                group_properties.chain_modifications.at[idx,
                                                        'charge_indicator'])
            self.IEP_dataset.at[idx, 'pka_data'] = (
                group_properties.chain_modifications.at[idx, mod_pka_scale])
            
        self.IEP_dataset.dropna(subset=['pka_data'], inplace=True)

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
                                (pH-self.IEP_dataset['pka_data'].values))))
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
                                self.IEP_dataset['pka_data'].values
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

    def molar_mass(self, molecule_part='all'):
        """
        Calculate the molar mass of the polyelectrolyte.
        
        Works only in case absolute abundances of the residues/repeating units
        are known, so when the sequence or absolute abundances are given.

        Parameters
        ----------
        molecule_part : str, optional
            Defines if the molar mass is calculated for the full modified
            polymer ('all'), only the main chain ('main_chain'), or only the
            modifications ('mods'). Default is 'all'.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        if self.abundance_unit in ['absolute', 'sequence']:
            main_chain_contribution = np.sum(
                self.main_chain['abundance_' + self.abundance_unit].values *
                self.main_chain['molar_mass_residue'].values)
            modification_contribution = np.sum(
                self.modifications['abundance_' + self.abundance_unit].values *
                self.modifications['molar_mass_residue'].values)
            molar_mass = main_chain_contribution + modification_contribution

            if molecule_part == 'all':
                return molar_mass
            elif molecule_part == 'main_chain':
                return main_chain_contribution
            elif molecule_part == 'mods':
                return modification_contribution
        else:
            return None

    def mean_residue_molar_mass(self, molecule_part='all'):
        """
        Calculate the mean residue molecular weight of proteins.

        Calculation is done without taking C and N terminus into account if
        only relative abundances are known and not the entire sequence.

        Parameters
        ----------
        molecule_part : str, optional
            Defines if the mean residue molecular mass is calculated for the
            full modified polymer ('all'), only the main chain ('main_chain'),
            or only the modifications ('mods'). Default is 'all'.

        Returns
        -------
        mean_residue_mass : float
            The mean residue molecular weight based on the abundances of amino
            acids of the molecule part given with molecule_part.

        """
        main_chain_contribution = np.sum(
            self.main_chain['molar_mass_residue'].values *
            self.main_chain['abundance_norm'].values)
        modification_contribution = np.sum(
            self.modifications['molar_mass_residue'].values *
            self.modifications['abundance_norm'].values)
        mean_residue_molar_mass = (
            main_chain_contribution + modification_contribution)

        if molecule_part == 'all':
            return mean_residue_molar_mass
        elif molecule_part == 'main_chain':
            return main_chain_contribution
        elif molecule_part == 'mods':
            return modification_contribution

    def n_content(self, molecule_part='all'):
        """
        Calculate the nitrogen mass fraction of polymers.

        Parameters
        ----------
        molecule_part : str, optional
            Defines if the n content is calculated for the full modified
            polymer ('all'), only the main chain ('main_chain'), or only the
            modifications ('mods'). Default is 'all'.

        Returns
        -------
        n_content : float
            The n_content based on the abundances of amino acids of the
            molecule part given with molecule_part.

        """
        mean_chain = self.mean_residue_molar_mass(molecule_part='main_chain')
        mean_mods = self.mean_residue_molar_mass(molecule_part='mods')

        main_chain_contribution = (
            np.sum(
                self.main_chain['N_content_residue'] *
                self.main_chain['abundance_norm'] *
                self.main_chain['molar_mass_residue']) /
            mean_chain
            )
        modification_contribution = (
            np.sum(
                self.modifications['N_content_residue'] *
                self.modifications['abundance_norm'] *
                self.modifications['molar_mass_residue']) /
            mean_mods
            ) if mean_mods != 0 else 0
        n_content = main_chain_contribution + modification_contribution

        if molecule_part == 'all':
            return n_content
        elif molecule_part == 'main_chain':
            return main_chain_contribution
        elif molecule_part == 'mods':
            return modification_contribution
