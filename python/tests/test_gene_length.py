"""Validate spatial factorization against the independently assembled full operator."""
from dataclasses import replace
import unittest

import numpy as np
from scipy.sparse.linalg import spsolve

from cpadynamics.full_model import FullModel
from cpadynamics.gene_length import GeneLengthModel, length_population
from cpadynamics.observables import cleavage_profile
from cpadynamics.parameters import working_parameters


class GeneLengthTests(unittest.TestCase):
    def test_matches_full_sparse_model_including_pas_at_first_node(self):
        p = replace(working_parameters(), k_in=0.0002)
        lengths = [100, 300, 1000]
        for m in (1, 5):
            for ef in (0., 42000.):
                with self.subTest(m=m, ef=ef):
                    profiles = GeneLengthModel(p, m, lengths, 400).evaluate(ef)
                    for i, length in enumerate(lengths):
                        local_p = replace(p, PASposition=length, geneLength_bp=length+400)
                        model = FullModel(local_p, m)
                        occupation = spsolve(-(model.a0+ef*model.ae), model.unit_source)
                        solution = model.summarize(np.r_[occupation, ef, 1/p.k_in])
                        np.testing.assert_allclose(profiles.pol_per_flux[i], occupation.sum(), rtol=1e-11)
                        np.testing.assert_allclose(profiles.e_per_flux[i], model.e_counts @ occupation, atol=1e-11)
                        reference = cleavage_profile(solution.R, solution.RHE, local_p)
                        np.testing.assert_allclose(profiles.cleavage_cdf[i], reference.cdf, atol=1e-12)
                        np.testing.assert_allclose(profiles.cad50_bp[i], reference.cad50_bp, rtol=1e-11)

    def test_zero_e_and_no_cleavage_do_not_fabricate_cad(self):
        p = working_parameters()
        for parameters, ef in ((p, 0.), (replace(p, kc=0), 40000.)):
            profile = GeneLengthModel(parameters, 5, [2500, 5000], 500).evaluate(ef)
            self.assertTrue(np.isnan(profile.cad50_bp).all())
            np.testing.assert_array_equal(profile.cleavage_cdf, 0.)

    def test_lattice_bin_probabilities_and_shared_conservation(self):
        lengths, weights, mass = length_population(minimum=100, maximum=1000)
        self.assertAlmostEqual(weights.sum(), 1.)
        self.assertEqual(weights[-1], 0.)
        self.assertTrue(0 < mass < 1)
        p = replace(working_parameters(), k_in=0.0002)
        pools, profile = GeneLengthModel(p, 1, lengths, 400).solve_shared_pools(weights)
        self.assertTrue(0 < pools['Pol_free'] < p.Pol_total)
        self.assertTrue(0 < pools['E_free'] < p.E_total)
        flux = p.k_in*pools['Pol_free']
        self.assertAlmostEqual(pools['Pol_free']+10000*flux*(weights @ profile.pol_per_flux), p.Pol_total, places=6)
        self.assertAlmostEqual(pools['E_free']+10000*flux*(weights @ profile.e_per_flux), p.E_total, places=6)


if __name__ == '__main__':
    unittest.main()
