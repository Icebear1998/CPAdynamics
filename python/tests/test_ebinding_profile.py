"""Check the profile averaging convention independently of a steady solve."""
from dataclasses import replace
import unittest

import numpy as np

from cpadynamics.full_model import FullModel
from plot_ebinding_profile import binding_profiles, profile_parameters


class BindingProfileTests(unittest.TestCase):
    def test_recognized_and_unrecognized_populations_are_weighted_together(self):
        model = FullModel(replace(profile_parameters(), geneLength_bp=300, PASposition=200), 2)
        state = model.initial_state()
        # First node: two R(1,0) and one R(2,2), giving E=2/3, p=4/3.
        state[model.r_states.index((1, 0))] = 2
        state[model.r_states.index((2, 2))] = 1
        # PAS: three R(0,0) and one RHE(2,1), giving E=1/4, p=1/2.
        state[model.sr + model.r_states.index((0, 0))] = 3
        state[model.r_size + model.rhe_states.index((2, 1))] = 1
        avg_e, avg_p = binding_profiles(model, model.summarize(state))
        np.testing.assert_allclose(avg_e[:2], [2/3, 1/4])
        np.testing.assert_allclose(avg_p[:2], [4/3, 1/2])
        self.assertTrue(np.isnan(avg_e[2]))
        self.assertTrue(np.isnan(avg_p[2]))

    def test_weighted_profile_recovers_total_bound_e(self):
        parameters = replace(profile_parameters(), geneLength_bp=1000, PASposition=600)
        model = FullModel(parameters, 5)
        solution = model.solve_steady_state()
        avg_e, avg_p = binding_profiles(model, solution)
        polymerases = solution.R.copy()
        polymerases[model.pas:] += solution.RHE
        self.assertAlmostEqual(avg_e @ polymerases + solution.E_free, parameters.E_total, places=6)
        self.assertTrue(np.all((0 <= avg_e) & (avg_e <= avg_p) & (avg_p <= 5)))


if __name__ == '__main__':
    unittest.main()
