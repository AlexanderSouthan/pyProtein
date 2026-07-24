# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from . import amino_acid_properties
from little_helpers import table_of_elements


class protein:
    """
    Does basic calculations on protein, enzymes and polypeptides.

    Calculations include e.g. the isoelectric point similar to the "ExPASy
    Compute pI/Mw tool".
    """

    def __init__(self, mode, aa_abundances, pka_data='pka_bjellqvist',
                 **kwargs):
        """
        Initialize a protein instance.

        Parameters
        ----------
        mode : str
            Mode defining how amino acid abundances are passed to the class.
            Currently can be 'mmol_g', 'sequence', 'absolute' or
            'res_per_1000'.
        aa_abundances : list or string
            Gives the amino acid abundances in the unit defined by mode.
            Format is as follows:
            mode == 'mmol_g' : list of float or ndarray
                For mode 'protein', contains abundance of amino acids
                and C- and N-terminus in the order ['D', 'N', 'T', 'S', 'E',
                'Q', 'G', 'A', 'C', 'V', 'M', 'I', 'L', 'Y', 'F', 'H', 'K',
                'R', 'P', 'W', 'Hyp'm 'Hyl']. Because this is normalized to the
                mass, pay attention to give the data for the modified protein
                in case you want to specify modifications with the mod kwargs
                described below. In the future, the calculations should take
                the masses of the modifications into account, and then it
                should be possible to also give the value normalized to the
                mass of the unmodified protein, but this is not implemented
                yet.
            mode == 'sequence' : str
                Alternative to mmol_g. Has to be a string consisting of one
                letter codes for amino acids.
            mode == 'absolute' : list of float or ndarray
                Alternative to mmol_g and sequence. Is a list in the same order
                as for mmol_g, but with absolute abundances of amino acids,
                e.g. deducted from the sequence.
            mode == 'res_per_1000' : list of float or ndarray
                List in the same order as mmol_g, but abundances normalized to
                1000 amino acid residues.
        pka_data : string, optional
            Gives the pKa dataset used for charge and IEP calculations. Allowed
            values are 'pka_ipc_protein', 'pka_ipc2_protein', 'pka_emboss' or
            'pka_bjellqvist'. Default is 'pka_bjellqvist'.
        **kwargs in case of modifications are present
            mod_types : list of strings
                Allowed values are the index values of
                amino_acid_properties.chain_modifications, defaults to empty
                list. Currently only 'N_term' and 'C_term' are allowed in order
                to account for the N and C terminus. Will be expanded in the
                future.
            mod_abundances : list of float
                The abundances of the modifications given by mod_types. Must be
                given in the same units as the input data for the main chain
                and contain the same amount of elements like mod_types. Default
                is an empty list. Usually, is [1, 1] to account for one C
                terminus and one N terminus. If it is given in mmol/g, then pay
                attention to give the number normalized to the mass of the 
                modified version of the protein. In the future, the
                calculations should take the masses of the modifications into
                account, and then it should be possible to also give the value
                normalized to the mass of the unmodified protein, but this is
                not implemented yet. 
            mod_sites : list of strings
                The sites on the main chain which are modified, so allowed
                values are the index values in self.main_chain or 'any'. 'any'
                means that the modification site is undefined. In this case,
                the modification is only important e.g. for molar mass
                calculations, but irrelevant e.g. for IEP calculations. Must
                contain the same amount of elements like mod_types. Default is
                a list the same size as mod_types with only 'any' entries.
            mod_distrib : list of strings
                Defines how the modifications are distributed on the mod_sites.
                Allowed values are 'equal' (all mod_sites are modified equally,
                if abundances allow) or 'sequential' (first the first entry in
                mod_sites is modifies completely, then the second, and so on).
                Default is 'equal'.
            pka_mods : list of strings
                The pKa sclaes used for the different modifications. Must
                contains the same amount of entires like mod_types. Default is
                a list with the same size as mod_types with only
                'pka_other' entries.

        Returns
        -------
        None.

        """
        self.mode = mode
        self.amino_acid_sequence = None
        self.number_of_residues = None

        allowed_modes = ['mmol_g', 'sequence', 'absolute', 'res_per_1000']

        if self.mode == allowed_modes[0]:  # 'mmol_g'
            self.abundance_unit = 'mmol_g'
            abundance = aa_abundances
        elif self.mode == allowed_modes[1]:  # 'sequence'
            self.abundance_unit = 'sequence'
            self.amino_acid_sequence = aa_abundances
            abundance = []
        elif self.mode == allowed_modes[2]:  # 'absolute'
            self.abundance_unit = 'absolute'
            abundance = aa_abundances
        elif self.mode == allowed_modes[3]:  # 'res_per_1000'
            self.abundance_unit = 'res_per_1000'
            abundance = aa_abundances
        else:
            raise ValueError('No valid option for mode given. Allowed values '
                             'are {}'.format(allowed_modes))

        if isinstance(abundance, np.ndarray):
            abundance = abundance.tolist()
        abundance = [float(i) for i in abundance]

        mod_types = kwargs.get('mod_types', [])
        mod_abundances = kwargs.get('mod_abundances', [])
        mod_sites = kwargs.get('mod_sites',
                                    ['any']*len(mod_types))
        pka_mods = kwargs.get(
            'pka_mods', ['pka_other']*len(mod_types))
        mod_distrib = kwargs.get(
            'mod_distrib', ['equal']*len(mod_types))

        # Make sure that kwargs contain equal numbers of elements.
        assert (
            len(mod_types) == len(mod_abundances) ==
            len(pka_mods) == len(mod_sites) == 
            len(mod_distrib)), (
            'Length mismatch, mod_types, mod_abundances, pka_mods, '
            'mod_sites, and mod_distrib must contain an equal number of '
            'elements, but contain '
            '{}, {}, {}, {}, and {} elements.'.format(
                len(mod_types), len(mod_abundances),
                len(pka_mods), len(mod_sites),
                len(mod_distrib)))

        self.pka_data = pka_data

        self.initialize_main_chain(abundance)
        self.initialize_modifications(
            mod_types, mod_abundances, mod_sites, pka_mods, mod_distrib)
        self.normalize_abundances()        
        self.initialize_pka_dataset()

    def initialize_main_chain(self, abundance):
        """
        Transform input data into DataFrame.

        Resulting DataFrame contains amino acid abundances (both as given and
        normalized), N contents and molar masses of resisues.
        """
        self.main_chain = pd.DataFrame(
            [], index=amino_acid_properties.amino_acids.index)

        if self.abundance_unit in ['mmol_g', 'absolute', 'res_per_1000']:
            # if less abundance values than entries in amino acid table are
            # given, remaining abundances are set to zero
            abundance.extend(
                    (len(self.main_chain.index) -
                     len(abundance)) * [0])
            self.main_chain[
                'abundance_' + self.abundance_unit] = abundance

        elif self.abundance_unit == 'sequence':
            # occurences of the first 20 amino acids in the sequence are
            # counted, needs to be adapted if more than 20 amino acids such
            # as hydroxylysine are added
            self.number_of_residues = len(self.amino_acid_sequence)
            for key in self.main_chain.index[:20]:
                abundance.append(self.amino_acid_sequence.count(key))
            # Hyp, Hyl abundances
            abundance.extend(
                (len(self.main_chain.index) - len(abundance)) * [0])
            self.main_chain[
                'abundance_' + self.abundance_unit] = abundance

        self.norm_factor = np.sum(abundance)

    def initialize_modifications(self, mod_types, mod_abundances, mod_sites,
                                 pka_mods, mod_distrib):
        self.modifications = pd.DataFrame(
            [], index=amino_acid_properties.chain_modifications.index,
            columns=['abundance_' + self.abundance_unit, 'modified_residues',
                     'pka_scale', 'theo_max'])

        self.modifications['abundance_' + self.abundance_unit] = np.nan
        self.modifications['modified_residues'] = self.modifications['modified_residues'].astype('object')
        for mod_type, abundance, site, pka_scale, distrib in zip(
                mod_types, mod_abundances, mod_sites,
                pka_mods, mod_distrib):
            self.modifications.at[
                mod_type, 'abundance_' + self.abundance_unit] = abundance
            self.modifications.at[mod_type, 'pka_scale'] = pka_scale
            self.modifications.at[mod_type, 'mod_distrib'] = distrib
            if site == 'any':
                self.modifications.at[mod_type, 'modified_residues'] = []
                self.modifications.at[mod_type, 'theo_max'] = 'undefined'
                mod_keys = self.main_chain.index.values
            else:
                self.modifications.at[mod_type, 'modified_residues'] = site
                mod_keys = self.modifications.at[mod_type, 'modified_residues'] 

            theo_max = self.main_chain.loc[
                mod_keys, 'abundance_' + self.abundance_unit].sum()
            self.modifications.at[mod_type, 'theo_max'] = theo_max
            if theo_max < abundance:
                raise ValueError(
                    'Theoretical upper limit for modification is {}, '
                    'but {} was given'.format(theo_max, abundance))
        self.modifications.dropna(
            subset=['abundance_' + self.abundance_unit], inplace=True)

        # Subtract the modfication abundances from main chain abundances. This
        # is important for the IEP calculations if IEP relevant residues are
        # modified.
        self.main_chain['abundance_iep_' + self.abundance_unit] = (
            self.main_chain['abundance_' + self.abundance_unit])
        for curr_mod, distrib in zip(mod_types, mod_distrib):
            remaining = self.modifications.at[curr_mod, 'abundance_' + self.abundance_unit]
            mod_subset = self.main_chain.loc[
                self.modifications.at[curr_mod, 'modified_residues'],
                'abundance_iep_' + self.abundance_unit]
            # the modification abundance is subtracted evenly for all
            # modification sites, if abundances allow (negative abundances are
            # not allowed).
            if distrib == 'equal':
                while remaining > 1E-12:
                    active = mod_subset[mod_subset > 0]
                    if active.empty:
                        break
                    share = remaining/len(active)
                    deduction = active.clip(upper=share)
                    mod_subset.loc[active.index] -= deduction
                    remaining -= deduction.sum()
            # the modification abundance is subtrated from one after the other
            # modification sites
            elif distrib == 'sequential':
                for idx in self.modifications.at[curr_mod, 'modified_residues']:
                    if remaining <= 0:
                        break
                    deduction = min(mod_subset.loc[idx], remaining)
                    mod_subset.loc[idx] -= deduction
                    remaining -= deduction
            else:
                raise ValueError(
                    'Allowed values for mod_distrib are \'equal\' and '
                    '\'sequential\', but mod_distrib is {}.'.format(
                        mod_distrib))
            # write the unmodified amino acid abundances into main chain table
            self.main_chain.loc[mod_subset.index, 'abundance_iep_' + self.abundance_unit] = mod_subset

    def normalize_abundances(self):
        # Normalize abundances to sum to make different inputs comparable.
        self.main_chain['abundance_norm'] = (
            self.main_chain['abundance_' + self.abundance_unit].values /
            self.norm_factor)
        self.main_chain['abundance_iep_norm'] = (
            self.main_chain['abundance_iep_' + self.abundance_unit].values /
            self.norm_factor)

        # Modification abundances are normalized on the basis of the input data
        # of the main chain because the modifications do not introduce any
        # additional residues/repeating units
        self.modifications['abundance_norm'] = (
            self.modifications['abundance_' + self.abundance_unit].values /
            self.norm_factor)

    def initialize_pka_dataset(self):
        # Write the amino acid characteristics important for IEP calculation
        # into the IEP dataset.
        self.IEP_dataset = pd.DataFrame([], index=self.main_chain.index)
        self.IEP_dataset['abundance_norm'] = self.main_chain['abundance_iep_norm']
        self.IEP_dataset['charge_indicator'] = (
            amino_acid_properties.amino_acids['charge_indicator'])
        self.IEP_dataset['pka_data'] = amino_acid_properties.amino_acids[
            self.pka_data]

        # Write the modification characteristics important for IEP calculation
        # into the IEP dataset
        for idx in self.modifications.index:
            mod_pka_scale = self.modifications.at[idx, 'pka_scale']
            self.IEP_dataset.at[idx, 'abundance_norm'] = self.modifications.at[
                idx, 'abundance_norm']
            self.IEP_dataset.at[idx, 'charge_indicator'] = (
                amino_acid_properties.chain_modifications.at[
                    idx, 'charge_indicator'])
            self.IEP_dataset.at[idx, 'pka_data'] = (
                amino_acid_properties.chain_modifications.at[
                    idx, mod_pka_scale])

        # Drop all entries that are not relevant for IEP because they have no
        # pKa values.
        self.IEP_dataset.dropna(subset=['pka_data'], inplace=True)

    def charge(self, pH, mode='overall', output='norm', values_only=False):
        """
        Calculate the net charge of the protein at a given pH.

        Parameters
        ----------
        pH : float
            The pH used for net charge calculation.
        mode : str, optional
            What kind of charge is calculated. Allowed values are 'overall', so
            an overall charge of the protein is calculated, or 'individual', so
            that the charges of the individual groups are calculated, or
            'aggregate', so that the positive and negative charges are summed
            up and are both returned. Default
            is 'overall'.
        output : str, optional
            Define if the output is normalized ('norm') or gives the charge
            value in the unit of the initially given abundances ('absolute').
            Default is 'norm'.
        values_only : boolean, optional
            If True, only the numbers are returned as an array without the
            context of a pandas DataFrame or Series. Default is False.

        Returns
        -------
        charge : pandas Series or DataFrame
            A pandas object containing the pH values as index and the charges
            in the columns, depending on the mode. With mode = 'overall', a
            Series is returned with the overall charge. With
            mode = 'individual', a DataFrame is returned with the charges
            caused by the individual amino acids in the columns. With
            mode = 'aggregate', a DataFrame is returned with the aggregated
            positive and negative charges in the columns.

        """
        charge = self.charge_curve(
            ph_range=[pH, pH], data_points=1, mode=mode, output=output)
        if mode == 'overall' and values_only:
            charge = charge.item()
        elif (mode in ['individual', 'aggegate']) and values_only:
            charge = charge.to_numpy()

        return charge

    def charge_curve(
            self, ph_range=[0, 14], data_points=100, mode='overall', 
            output='norm'):
        """
        Calculate the charge curve of the protein in a given pH range.

        Parameters
        ----------
        ph_range : list of floats, optional
            First value is the lower pH limit and second value the upper pH
            limit used for calculations. The default is [0, 14].
        data_points : int, optional
            Number of data points calculated. The default is 100.
        mode : str, optional
            What kind of charge is calculated. Allowed values are 'overall', so
            an overall charge of the protein is calculated, or 'individual', so
            that the charges of the individual groups are calculated, or
            'aggregate', so that the positive and negative charges are summed
            up and are both returned. Default
            is 'overall'.
        output : str, optional
            Define if the output is normalized ('norm') or gives the charge
            value in the unit of the initially given abundances ('absolute').
            Default is 'norm'.

        Returns
        -------
        curve : pandas Series or DataFrame
            A pandas object containing the pH values as index and the charges
            in the columns, depending on the mode. With mode = 'overall', a
            Series is returned with the overall charge. With
            mode = 'individual', a DataFrame is returned with the charges
            caused by the individual amino acids in the columns. With
            mode = 'aggregate', a DataFrame is returned with the aggregated
            positive and negative charges in the columns.

        """
        pH = np.linspace(ph_range[0], ph_range[1], data_points)
        charges = (
            self.IEP_dataset['charge_indicator'].values *
            self.IEP_dataset['abundance_norm'].values /
            (1+10**(self.IEP_dataset['charge_indicator'].values *
                    (pH[:, np.newaxis] -
                     self.IEP_dataset['pka_data'].values))))

        if mode == 'overall':
            curve = pd.Series(
                np.sum(charges, axis=1), index=pd.Index(pH, name='pH value'),
                name='overall charge')
        elif mode == 'individual':
            curve = pd.DataFrame(
                charges, columns = self.IEP_dataset.index,
                index=pd.Index(pH, name='pH value'))
        elif mode == 'aggregate':
            curve = pd.DataFrame({
                'negative charges': np.clip(charges, None, 0).sum(axis=1),
                'positive charges': np.clip(charges, 0, None).sum(axis=1)
                }, index=pd.Index(pH, name='pH value'))
        else:
            raise ValueError(
                'No valid mode given. Allowed values are \'overall\', '
                '\'aggregate\' or \'individual\', but {} was given.'.format(
                    mode))

        if output == 'absolute':
            corr_factor = self.norm_factor
        elif output == 'norm':
            corr_factor = 1
        else:
            raise ValueError(
                'No valid output given. Allowed values are \'norm\' or '
                '\'absolute\', but {} was given.'.format(output))
        
        curve = curve * corr_factor
        return curve

    def IEP(self, ph_range=[0, 14]):
        """
        Calculate the isoelectric point (IEP) of the protein.

        The IEP is the pH value at which the protein has net charge of zero.

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
            IEP = brentq(
                self.charge, ph_range[0], ph_range[1],
                args=('overall', 'norm', True))
        except ValueError:
            IEP = np.nan
        return IEP

    def molar_mass(self, molecule_part='all'):
        """
        Calculate the molar mass of the protein.

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
        float
            The molar mass of the protein.

        """
        if self.abundance_unit in ['absolute', 'sequence']:
            main_chain_contribution = np.sum(
                self.main_chain['abundance_' + self.abundance_unit].values *
                amino_acid_properties.amino_acids['molar_mass_residue'].values)
            modification_contribution = np.sum(
                self.modifications['abundance_' + self.abundance_unit].values *
                amino_acid_properties.chain_modifications.loc[
                    self.modifications.index, 'molar_mass_residue'].values)
            molar_mass = main_chain_contribution + modification_contribution

            allowed_parts = ['all', 'main_chain', 'mods']
            if molecule_part == allowed_parts[0]:  # 'all'
                return molar_mass
            elif molecule_part == allowed_parts[1]:  # 'main_chain'
                return main_chain_contribution
            elif molecule_part == allowed_parts[2]:  # 'mods'
                return modification_contribution
            else:
                raise ValueError('No valid value given for molecule_part. '
                                 'Allowed values are {}, but \'{}\' was given'
                                 '.'.format(allowed_parts, molecule_part))
        else:
            raise Exception('Molar mass can only be calculated if mode is '
                            '\'absolute\' or \'sequence\', however mode is '
                            '\'{}\'.'.format(self.abundance_unit))

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
            amino_acid_properties.amino_acids['molar_mass_residue'].values *
            self.main_chain['abundance_norm'].values)
        modification_contribution = np.sum(
            amino_acid_properties.chain_modifications.loc[
                self.modifications.index, 'molar_mass_residue'].values *
            self.modifications['abundance_norm'].values)
        mean_residue_molar_mass = (
            main_chain_contribution + modification_contribution)

        allowed_parts = ['all', 'main_chain', 'mods']
        if molecule_part == allowed_parts[0]:  # 'all'
            return mean_residue_molar_mass
        elif molecule_part == allowed_parts[1]:  # 'main_chain'
            return main_chain_contribution
        elif molecule_part == allowed_parts[2]:  # 'mods'
            return modification_contribution
        else:
            raise ValueError('No valid value given for molecule_part. '
                             'Allowed values are {}, but \'{}\' was given'
                             '.'.format(allowed_parts, molecule_part))

    def n_content(self):
        """
        Calculate the nitrogen mass fraction of polymers. This is a legacy
        function which calls the more versatile function
        self.elemental_composition.

        Returns
        -------
        n_content : pandas Series
            The n_content of the main chain, the modifications and the
            complete modified protein based on the abundances of amino acids
            and modifications.

        """
        n_content = self.elemental_composition(atom_labels=['N'])
        
        return n_content/100

    def elemental_composition(self, atom_labels=['C', 'H', 'N', 'O', 'S']):
        """
        Calculates the elemental composition of the protein.

        Returns
        -------
        element_comp : pandas DataFrame
            Contains information about the normalized atom counts of the main
            chain, the modifications, the complete modified protein, as well as
            the mass percentage of the different atoms for the three parts.
            The index contains all elements specified in atom_labels.

        """
        element_count = {}
        element_count['main_chain'] = pd.DataFrame(
            [], index=amino_acid_properties.amino_acids.index)
        element_count['mods'] = pd.DataFrame(
            [], index=amino_acid_properties.chain_modifications.index)
        
        all_atoms = ['C', 'H', 'N', 'O', 'S']
        for curr_atom in all_atoms:
            # calculate the normalized atom numbers for the main chain residues
            element_count['main_chain'][curr_atom] = (
                self.main_chain['abundance_norm'] *
                amino_acid_properties.amino_acids[curr_atom + '_residue'])
            # calculate the normalized atom numbers for the modification types
            element_count['mods'][curr_atom] = (
                self.modifications['abundance_norm'] *
                amino_acid_properties.chain_modifications[curr_atom])
        element_count['mods'].dropna(inplace=True)

        # calculate the aggregates for main chain and modification atom numbers
        element_comp = element_count['main_chain'].sum().to_frame(
            name='atom count main chain (norm)')
        element_comp['atom count mods (norm)'] = element_count['mods'].sum()
        element_comp['atom count all (norm)'] = (
            element_comp['atom count mods (norm)'].add(
            element_comp['atom count main chain (norm)'], fill_value=0))

        # calculate the mass percentages of the different atoms
        mass_norm = {}
        atomic_masses = table_of_elements.loc[all_atoms, 'atomic weight']
        mass_norm['main_chain'] = (
            element_comp['atom count main chain (norm)'].multiply(
                atomic_masses))
        mass_norm['mods'] = (
            element_comp['atom count mods (norm)'].multiply(
                atomic_masses))
        mass_norm['all'] = (
            element_comp['atom count all (norm)'].multiply(
                atomic_masses))
        
        # write mass percentages into output table
        element_comp['mass % main chain'] = (
            mass_norm['main_chain']/mass_norm['main_chain'].sum()*100)
        element_comp['mass % mods'] = (
            mass_norm['mods']/mass_norm['mods'].sum()*100)
        element_comp['mass % all'] = (
            mass_norm['all']/mass_norm['all'].sum()*100)
        
        return element_comp.loc[atom_labels]

    def check_abundances(self):
        """
        This is just a diagnosis function to check if the input data makes sense.
        
        For input data in mmol/g, the function calculates the mass of all main
        chain elements and modifications and outputs the total mass which
        should be equal to 1 g. 
        For input data in redidues per 1000, the function calculates the total
        number of residues and outputs that number, which should be equal to
        1000.

        Returns
        -------
        None.

        """
        if self.abundance_unit == 'mmol_g':
            main_chain_mass = self.main_chain['abundance_mmol_g'].mul(
                amino_acid_properties.amino_acids['molar_mass_residue']).sum()/1000
            mod_mass = self.modifications['abundance_mmol_g'].mul(
                amino_acid_properties.chain_modifications['molar_mass_residue']).sum()/1000
            print('The data in mmol/g result in a main chain mass of {} g '
                  'and a mass of the modifications of {} g. The total mass '
                  'thus is {}, while a mass of 1.000 g is expected.'.format(
                      round(main_chain_mass, 4), round(mod_mass, 4),
                      round(main_chain_mass+mod_mass, 4)))
        elif self.abundance_unit == 'res_per_1000':
            number_of_residues = self.main_chain['abundance_res_per_1000'].sum()
            number_of_mods = self.modifications['abundance_res_per_1000'].sum()
            print(
                'Total number of residues is {}, while 1000 is '
                'expected.'.format(number_of_residues))
            print(
                'Total number of modifications is {}.'.format(number_of_mods))
        else:
            
            print('No check possible for input of type {}.'.format(
                self.abundance_unit))
