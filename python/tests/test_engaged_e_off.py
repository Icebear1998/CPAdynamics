"""The H-engaged E can detach, followed by rapid EH disassembly."""
from dataclasses import replace
import unittest

import numpy as np

from cpadynamics import Parameters, FullModel, cleavage_profile


class EngagedEOffTests(unittest.TestCase):
    def test_each_anchor_loss_releases_one_e_and_preserves_polymerase_and_p(self):
        base = replace(Parameters(), geneLength_bp=300, PASposition=200, kEoff_engaged=0)
        rate = .07
        for m in (1, 5):
            old = FullModel(base, m)
            new = FullModel(replace(base, kEoff_engaged=rate), m)
            for j in range(new.n_pas):
                for h_local, (p, e) in enumerate(new.rhe_states):
                    with self.subTest(M=m, node=j, p=p, e=e):
                        source = new.r_size+j*new.sh+h_local
                        target = (new.pas+j)*new.sr+new.r_states.index((p, e-1))
                        y = np.zeros(new.size+2)
                        y[source] = 3
                        difference = new.rhs(0, y)-old.rhs(0, y)
                        expected = np.zeros_like(y)
                        expected[source] = -3*rate
                        expected[target] = 3*rate
                        expected[-2] = 3*rate  # exactly one E per event
                        np.testing.assert_allclose(difference, expected, atol=1e-13)

    def test_h_off_and_anchor_off_are_separate_competing_routes(self):
        p = replace(Parameters(), geneLength_bp=100, PASposition=100,
                    k_e=0, k_e2=0, kc=0, kEon=0, kEoff=0,
                    kPon_min=0, kPon_slope=0, kPoff=0, kHon=0, k_in=0,
                    kHoff=.2, kEoff_engaged=.05)
        model = FullModel(p, 1)
        y = np.zeros(model.size+2)
        y[model.r_size] = 1  # RHE(1,1)
        derivative = model.rhs(0, y)
        self.assertAlmostEqual(derivative[model.r_states.index((1, 1))], .2)
        self.assertAlmostEqual(derivative[model.r_states.index((1, 0))], .05)
        self.assertAlmostEqual(derivative[model.r_size], -.25)
        self.assertAlmostEqual(derivative[-2], .05)
        self.assertAlmostEqual(derivative[-1], 0)

    def test_zero_rate_preserves_previous_full_model_results(self):
        for m, expected in [(1, 1535.6535807633459), (5, 551.4944450683093)]:
            p = replace(Parameters(), kHon=4, kHoff=2, kc=0.13, kEoff_engaged=0)
            s = FullModel(p, m).solve_steady_state()
            self.assertAlmostEqual(cleavage_profile(s.R, s.RHE, p).cad50_bp, expected, places=6)

    def test_nonzero_rate_conserves_resources_at_steady_state(self):
        for m in range(1, 7):
            for rate in (.05, .5):
                with self.subTest(M=m, rate=rate):
                    s = FullModel(replace(Parameters(), kEoff_engaged=rate), m).solve_steady_state()
                    self.assertGreaterEqual(s.state.min(), -1e-10)
                    self.assertLess(abs(s.e_residual), 1e-6)
                    self.assertLess(abs(s.pol_residual), 1e-6)
                    self.assertLess(s.rhs_max_abs, 1e-7)

    def test_nonzero_rate_time_integration_agrees_with_steady_solution(self):
        p = replace(Parameters(), geneLength_bp=1000, PASposition=600,
                    Pol_total=100, E_total=200, kEon=.001, kEoff_engaged=.05)
        for m in (1, 5):
            model = FullModel(p, m)
            steady = model.solve_steady_state()
            result = model.integrate(1000, t_eval=[0, 1000], rtol=1e-8, atol=1e-10)
            final = model.summarize(result.y[:, -1])
            np.testing.assert_allclose(final.state, steady.state, rtol=1e-6, atol=1e-6)
            self.assertLess(final.rhs_max_abs, 1e-7)
            self.assertLess(abs(final.e_residual), 1e-6)
            self.assertLess(abs(final.pol_residual), 1e-6)


if __name__ == '__main__':
    unittest.main()
