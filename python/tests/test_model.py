"""Scientific regression tests: run with PYTHONPATH=python python3 -m unittest discover -s python/tests -v."""
import unittest
from dataclasses import replace

import numpy as np

from cpadynamics import Parameters, FullModel, cleavage_profile, solve_equilibrium
from cpadynamics.states import internal_generator, state_map


class MatlabParityTests(unittest.TestCase):
    def test_supplied_matlab_cad_values(self):
        for m, expected in enumerate([1191, 728, 571, 490, 442, 407], 1):
            with self.subTest(M=m):
                p = replace(Parameters(), kHon=4, kHoff=2, kc=0.13, kEoff_engaged=0)
                solution = solve_equilibrium(p, m)
                profile = cleavage_profile(solution.R, solution.RHE, p)
                self.assertLess(abs(profile.cad50_bp - expected), 0.5)
                self.assertLess(abs(solution.e_residual), 1e-6)

    def test_first_bin_includes_cleavage_and_uses_right_edge(self):
        p = replace(Parameters(), k_e=0, k_e2=0, kc=2)
        profile = cleavage_profile(np.zeros(2), np.array([3., 0.]), p)
        np.testing.assert_equal(profile.distances_bp, [0, 100, 200])
        np.testing.assert_allclose(profile.cdf, [0, 1, 1])
        self.assertEqual(profile.cad50_bp, 50)

    def test_unreached_threshold_is_censored(self):
        p = replace(Parameters(), kc=0)
        profile = cleavage_profile(np.ones(2), np.zeros(2), p)
        self.assertTrue(np.isnan(profile.cad50_bp))
        self.assertEqual(profile.max_exit_cdf, 0)


class FullKineticsTests(unittest.TestCase):
    def test_m1_generator_matches_three_state_equations(self):
        p = Parameters()
        a = internal_generator(1, 0.7, 40000, p)
        expected = [[-.7, 1, 0], [.7, -1.1, .5], [0, .1, -.5]]
        np.testing.assert_allclose(a, expected, atol=1e-14)

    def test_engaged_factor_excluded_from_ordinary_unbinding(self):
        p = replace(Parameters(), kPon_min=0, kPon_slope=0, kPoff=0, kEon=0)
        states = state_map(5, recognized=True)
        a = internal_generator(5, 0, 0, p, recognized=True)
        index = {state: i for i, state in enumerate(states)}
        np.testing.assert_equal(a[:, index[(1, 1)]], 0)
        self.assertEqual(a[index[(3, 2)], index[(3, 3)]], 2*p.kEoff)
        np.testing.assert_allclose(a.sum(axis=0), 0)

    def test_h_dissociation_preserves_bound_e(self):
        p = replace(Parameters(), geneLength_bp=200, PASposition=100,
                    k_e=0, k_e2=0, kc=0, kEon=0, kEoff=0,
                    kPon_min=0, kPon_slope=0, kPoff=0, kHon=0, k_in=0, kEoff_engaged=0)
        model = FullModel(p, 1)
        y = np.zeros(model.size + 2)
        y[model.r_size] = 1
        derivative = model.rhs(0, y)
        self.assertEqual(derivative[2], p.kHoff)
        self.assertEqual(derivative[model.r_size], -p.kHoff)
        self.assertEqual(derivative[-2], 0)  # no E is released by H dissociation
        self.assertEqual(derivative[-1], 0)

    def test_recognition_selects_one_of_e_factors_only_after_pas(self):
        p = replace(Parameters(), geneLength_bp=300, PASposition=200,
                    k_e=0, k_e2=0, kc=0, kEon=0, kEoff=0, kHoff=0,
                    kPon_min=0, kPon_slope=0, kPoff=0, k_in=0, kEoff_engaged=0)
        model = FullModel(p, 5)
        local_r = model.r_states.index((3, 2))
        local_rhe = model.rhe_states.index((3, 2))
        y = np.zeros(model.size+2)
        y[local_r] = 1  # upstream of PAS: no recognition
        np.testing.assert_equal(model.rhs(0, y), 0)
        y[local_r] = 0
        y[model.sr+local_r] = 1  # PAS: counts are preserved
        expected = np.zeros_like(y)
        expected[model.sr+local_r] = -2*p.kHon
        expected[model.r_size+local_rhe] = 2*p.kHon
        np.testing.assert_allclose(model.rhs(0, y), expected, atol=1e-14)

    def test_elongation_preserves_internal_state_and_runoff_releases_e(self):
        p = replace(Parameters(), geneLength_bp=200, PASposition=100,
                    kc=0, kEon=0, kEoff=0, kHoff=0, kHon=0,
                    kPon_min=0, kPon_slope=0, kPoff=0, k_in=0, kEoff_engaged=0)
        model = FullModel(p, 5)
        local = model.rhe_states.index((4, 3))
        y = np.zeros(model.size+2)
        y[model.r_size+local] = 1
        expected = np.zeros_like(y)
        expected[model.r_size+local] = -p.k_e2
        expected[model.r_size+model.sh+local] = p.k_e2
        np.testing.assert_allclose(model.rhs(0, y), expected, atol=1e-14)
        y[model.r_size+local] = 0
        y[model.r_size+model.sh+local] = 1
        derivative = model.rhs(0, y)
        self.assertAlmostEqual(derivative[-2], 3*p.k_e2)
        self.assertAlmostEqual(derivative[-1], p.k_e2)

    def test_cleavage_at_pas_returns_all_resources(self):
        p = replace(Parameters(), geneLength_bp=200, PASposition=100,
                    k_e=0, k_e2=0, kEon=0, kEoff=0, kHoff=0,
                    kPon_min=0, kPon_slope=0, kPoff=0, kHon=0, k_in=0, kEoff_engaged=0)
        model = FullModel(p, 5)
        y = np.zeros(model.size + 2)
        source = model.r_size + model.rhe_states.index((3, 3))
        y[source] = 2
        derivative = model.rhs(0, y)
        self.assertAlmostEqual(derivative[source], -2*p.kc)
        self.assertAlmostEqual(derivative[-2], 6*p.kc)
        self.assertAlmostEqual(derivative[-1], 2*p.kc)

    def test_full_steady_state_is_physical_and_satisfies_dynamic_equations(self):
        for m in range(1, 8):
            with self.subTest(M=m):
                model = FullModel(Parameters(), m)
                s = model.solve_steady_state()
                self.assertGreaterEqual(s.state.min(), -1e-10)
                self.assertGreaterEqual(s.E_free, 0)
                self.assertGreaterEqual(s.Pol_free, 0)
                self.assertLess(abs(s.e_residual), 1e-6)
                self.assertLess(abs(s.pol_residual), 1e-6)
                self.assertLess(s.rhs_max_abs, 1e-7)
                np.testing.assert_allclose(model.rhs(0, s.state), 0, atol=1e-7)

    def test_zero_e_yields_no_recognition(self):
        p = replace(Parameters(), E_total=0)
        s = FullModel(p, 5).solve_steady_state()
        np.testing.assert_allclose(s.RHE, 0, atol=1e-12)
        self.assertEqual(s.E_free, 0)

    def test_analytic_jacobian_and_pool_conservation_away_from_steady_state(self):
        p = replace(Parameters(), geneLength_bp=500, PASposition=300, kEoff_engaged=.05)
        model = FullModel(p, 3)
        rng = np.random.default_rng(918)
        x = rng.uniform(.1, 1, model.size)
        y = np.r_[x, p.E_total-model.e_counts@x, p.Pol_total-x.sum()]
        derivative = model.rhs(0, y)
        self.assertAlmostEqual(model.e_counts@derivative[:-2]+derivative[-2], 0, places=8)
        self.assertAlmostEqual(derivative[:-2].sum()+derivative[-1], 0, places=8)
        direction = rng.normal(size=y.size)
        step = .001
        numerical = (model.rhs(0, y+step*direction)-model.rhs(0, y-step*direction))/(2*step)
        np.testing.assert_allclose(model.jacobian(0, y)@direction, numerical, rtol=1e-7, atol=1e-7)

    def test_time_integration_from_empty_gene_agrees_with_steady_solver(self):
        p = replace(Parameters(), geneLength_bp=1000, PASposition=600,
                    Pol_total=100, E_total=200, kEon=.001)
        for m in (1, 5):
            with self.subTest(M=m):
                model = FullModel(p, m)
                steady = model.solve_steady_state()
                result = model.integrate(1000, t_eval=[0, 1000], rtol=1e-8, atol=1e-10)
                final = model.summarize(result.y[:, -1])
                np.testing.assert_allclose(final.state, steady.state, atol=1e-6, rtol=1e-6)
                self.assertLess(final.rhs_max_abs, 1e-7)
                self.assertLess(abs(final.e_residual), 1e-6)
                self.assertLess(abs(final.pol_residual), 1e-6)


if __name__ == '__main__':
    unittest.main()
