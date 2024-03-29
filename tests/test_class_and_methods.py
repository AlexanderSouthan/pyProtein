#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 16 19:59:07 2021

@author: Alexander Southan
"""

import numpy as np
import pandas as pd
import unittest

from src.pyProtein import protein
from src.pyProtein import amino_acids


class TestInitialize(unittest.TestCase):

    def test_mode_sequence(self):

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
        bsa = protein(
                'sequence', bovine_serum_albumin_mature,
                pka_data='pka_bjellqvist', mod_types=['N_term', 'C_term'],
                mod_abundances=[1, 1],
                pka_scales=['pka_bjellqvist', 'pka_bjellqvist'])

        # IEP tests in the following. This also tests the charge method as this
        # is called by the IEP method.
        bsa_iep_bjellqvist = bsa.IEP()

        # Change pKa scale to 'pka_emboss'
        bsa.pka_data = 'pka_emboss'
        bsa.initialize_pka_dataset()
        bsa_iep_emboss = bsa.IEP()

        # Change pKa scale to 'pka_ipc_protein'
        bsa.pka_data = 'pka_ipc_protein'
        bsa.initialize_pka_dataset()
        bsa_iep_ipc_protein = bsa.IEP()

        self.assertAlmostEqual(bsa_iep_bjellqvist, 5.59867756, 4)
        self.assertAlmostEqual(bsa_iep_emboss, 5.6662233, 4)
        self.assertAlmostEqual(bsa_iep_ipc_protein, 5.40694206, 4)

        IEP_nan = bsa.IEP(ph_range=[0, 1])
        self.assertTrue(np.isnan(IEP_nan))

        # Charge curve calculation test in the following.
        ph_range = [3, 12]
        data_points = 121
        bsa_charge_curve = bsa.charge_curve(ph_range=ph_range,
                                            data_points=data_points)
        self.assertEqual(bsa_charge_curve[0, 0], ph_range[0])
        self.assertEqual(bsa_charge_curve[0, -1], ph_range[1])
        self.assertAlmostEqual(bsa_charge_curve[1, 0], 0.15918538, 4)
        self.assertAlmostEqual(bsa_charge_curve[1, -1], -0.24734945, 4)
        self.assertEqual(bsa_charge_curve.shape[0], 2)
        self.assertEqual(bsa_charge_curve.shape[1], data_points)

        # Molar mass calculation test in the following
        bsa_m_full = bsa.molar_mass()
        bsa_m_main_chain = bsa.molar_mass(molecule_part='main_chain')
        bsa_m_end_groups = bsa.molar_mass(molecule_part='mods')

        self.assertAlmostEqual(bsa_m_full, 66432.081799, 4)
        self.assertAlmostEqual(bsa_m_main_chain, 66414.066599, 4)
        self.assertAlmostEqual(bsa_m_end_groups, 18.0152, 4)

        # Mean residue molar mass calculation test in the following
        bsa_mrm_full = bsa.mean_residue_molar_mass()
        bsa_mrm_main_chain = bsa.mean_residue_molar_mass(
            molecule_part='main_chain')
        bsa_mrm_end_groups = bsa.mean_residue_molar_mass(molecule_part='mods')

        self.assertAlmostEqual(bsa_mrm_full, 113.9486823, 4)
        self.assertAlmostEqual(bsa_mrm_main_chain, 113.9177814, 4)
        self.assertAlmostEqual(bsa_mrm_end_groups, 0.03090085, 4)

        # Nitrogen content test in the following
        bsa_n_full = bsa.n_content()
        bsa_n_main_chain = bsa.n_content(molecule_part='main_chain')
        bsa_n_end_groups = bsa.n_content(molecule_part='mods')

        self.assertAlmostEqual(bsa_n_full, 0.1647125, 4)
        self.assertAlmostEqual(bsa_n_main_chain, 0.16471258, 4)
        self.assertAlmostEqual(bsa_n_end_groups, 0, 4)

        # test errors
        self.assertRaises(ValueError, protein, 'abs', [])
        self.assertRaises(ValueError, bsa.molar_mass, 'al')
        self.assertRaises(ValueError, bsa.mean_residue_molar_mass, 'al')
        self.assertRaises(ValueError, bsa.n_content, 'al')

    def test_mode_mmol_g(self):

        # test with a list as input
        gelatin_type_a_porcine = [0.42673911, 0, 0.15446764, 0.30639733,
                                  0.70074608, 0, 3.24912344, 1.10446511, 0,
                                  0.20060455, 0.04959423, 0.1113039,
                                  0.23785491, 0.00993443, 0.11623071,
                                  0.03802669, 0.25310048, 0.4615374,
                                  1.07443666, 0, 0.97079605, 0.07522209]
        gel_a = protein(
            'mmol_g', gelatin_type_a_porcine)

        # test with a numpy array as input
        gelatin_type_a_porcine_np = np.array(gelatin_type_a_porcine)
        gel_a_np = protein(
            'mmol_g', gelatin_type_a_porcine_np)

        # test errors
        self.assertRaises(Exception, gel_a.molar_mass)

    def test_absolute(self):
        col_1_a1_sequence = (
            'QLSYGYDEKSTGGISVPGPMGPSGPRGLPGPPGAPGPQGFQGPPGEPGEPGASGPMGPRGPPGP'
            'PGKNGDDGEAGKPGRPGERGPPGPQGARGLPGTAGLPGMKGHRGFSGLDGAKGDAGPAGPKGEP'
            'GSPGENGAPGQMGPRGLPGERGRPGAPGPAGARGNDGATGAAGPPGPTGPAGPPGFPGAVGAKG'
            'EAGPQGPRGSEGPQGVRGEPGPPGPAGAAGPAGNPGADGQPGAKGANGAPGIAGAPGFPGARGP'
            'SGPQGPGGPPGPKGNSGEPGAPGSKGDTGAKGEPGPVGVQGPPGPAGEEGKRGARGEPGPTGLP'
            'GPPGERGGPGSRGFPGADGVAGPKGPAGERGSPGPAGPKGSPGEAGRPGEAGLPGAKGLTGSPG'
            'SPGPDGKTGPPGPAGQDGRPGPPGPPGARGQAGVMGFPGPKGAAGEPGKAGERGVPGPPGAVGP'
            'AGKDGEAGAQGPPGPAGPAGERGEQGPAGSPGFQGLPGPAGPPGEAGKPGEQGVPGDLGAPGPS'
            'GARGERGFPGERGVQGPPGPAGPRGANGAPGNDGAKGDAGAPGAPGSQGAPGLQGMPGERGAAG'
            'LPGPKGDRGDAGPKGADGSPGKDGVRGLTGPIGPPGPAGAPGDKGESGPSGPAGPTGARGAPGD'
            'RGEPGPPGPAGFAGPPGADGQPGAKGEPGDAGAKGDAGPPGPAGPAGPPGPIGNVGAPGAKGAR'
            'GSAGPPGATGFPGAAGRVGPPGPSGNAGPPGPPGPAGKEGGKGPRGETGPAGRPGEVGPPGPPG'
            'PAGEKGSPGADGPAGAPGTPGPQGIAGQRGVVGLPGQRGERGFPGLPGPSGEPGKQGPSGASGE'
            'RGPPGPMGPPGLAGPPGESGREGAPGAEGSPGRDGSPGAKGDRGETGPAGPPGAPGAPGAPGPV'
            'GPAGKSGDRGETGPAGPTGPVGPVGARGPAGPQGPRGDKGETGEQGDRGIKGHRGFSGLQGPPG'
            'PPGSPGEQGPSGASGPAGPRGPPGSAGAPGKDGLNGLPGPIGPPGPRGRTGDAGPVGPPGPPGP'
            'PGPPGPPSAGFDFSFLPQPPQEKAHDGGRYYRA')
        col_1_a1_abundances = analyze_sequence(col_1_a1_sequence)
        col_1_a1_abundances['K'] -= 10
        col_1_a1_abundances['P'] -= 111
        col_1_a1_abundances['Hyl'] += 10
        col_1_a1_abundances['Hyp'] += 111

        # Collagen type I human uniprot.org
        col_1_a1 = protein(
            'absolute', col_1_a1_abundances.to_list(),
            pka_data='pka_bjellqvist', mod_types=['N_term', 'C_term'],
            mod_abundances=[1, 1])

        col_mm = col_1_a1.molar_mass()

    def test_res_per_1000(self):
        gel_a_wt_perc = np.array([
            5.91, 0, 1.72, 3.45, 10.50, 0, 24.89, 9.24, 0, 2.42, 0.90, 1.23,
            2.92, 0.60, 2.08, 0.75, 3.79, 8.42, 13.58, 0, 11.82, 1.31])
        gel_a_mmol_g = gel_a_wt_perc/amino_acids['molar_mass']*10
        gel_a_mmol_g_norm = (100 * gel_a_mmol_g / np.sum(gel_a_mmol_g))
        gel_a_per_1000 = gel_a_mmol_g_norm * 10

        gel_a = protein('res_per_1000', gel_a_per_1000.values)
        print(gel_a.mean_residue_molar_mass())

        # test errors
        self.assertRaises(Exception, gel_a.molar_mass)

def analyze_sequence(sequence):
    abundance = []
    # occurences of the first 20 amino acids in the sequence are
    # counted, needs to be adapted if more than 20 amino acids such
    # as hydroxylysine are added
    for key in amino_acids.index[:20]:
        abundance.append(sequence.count(key))
    # Hyp, Hyl abundances
    abundance.extend(
        (len(amino_acids.index) - len(abundance)) * [0])
    return pd.Series(abundance, index=amino_acids.index)
