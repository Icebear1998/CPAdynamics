"""Reference equilibrium closure for regression against CURRENT MATLAB.

Uses the same 100-point free-E grid, linear interpolation, PAS-frozen recognition
multiplier and flux equations. At fixed free E, the aggregate steady equations
are linear; a forward solve replaces fsolve without changing those equations.
The outer root is solved more tightly than the historical Python port.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import solve
from scipy.optimize import brentq

from .states import internal_components, state_map


@dataclass
class EquilibriumSolution:
    R: np.ndarray
    RHE: np.ndarray
    E_free: float
    Pol_free: float
    avg_E: np.ndarray
    avg_P: np.ndarray
    kHon_effective: float
    e_residual: float


def solve_equilibrium(parameters, m, num_grid_points=100):
    p = parameters
    if p.k_e <= 0 or p.k_e2 <= 0 or p.k_in <= 0:
        raise ValueError('Steady solver requires positive initiation and elongation')
    if num_grid_points < 2:
        raise ValueError('At least two interpolation points are required')
    states = np.asarray(state_map(m))
    n, pas, n_pas = p.geometry()
    grid = np.linspace(0, p.E_total, num_grid_points)
    means = np.zeros((num_grid_points, n, 2))
    rhs = np.zeros(len(states)); rhs[-1] = 1
    for i in range(n):
        a0, ae = internal_components(m, p.kPon_min+p.kPon_slope*i, p)
        for j, ef in enumerate(grid):
            a = a0 + ef*ae
            a[-1, :] = 1
            distribution = solve(a, rhs, check_finite=False)
            means[j, i] = distribution @ states  # P then E

    def profile(ef):
        if p.E_total == 0:
            return means[0]
        j = min(max(np.searchsorted(grid, ef, side='right') - 1, 0), len(grid)-2)
        fraction = (ef-grid[j])/(grid[j+1]-grid[j])
        return (1-fraction)*means[j] + fraction*means[j+1]

    def populations(ef):
        averages = profile(ef)
        h = p.kHon*averages[pas, 1]
        r, rhe = np.zeros(n), np.zeros(n_pas)
        r[:pas] = 1/p.k_e  # unit initiation flux
        block = np.array([[p.k_e+h, -p.kHoff], [-h, p.k_e2+p.kHoff+p.kc]])
        for j, i in enumerate(range(pas, n)):
            incoming = [p.k_e*r[i-1] if i else 1,
                        p.k_e2*rhe[j-1] if j else 0]
            r[i], rhe[j] = solve(block, incoming, check_finite=False)
        flux = p.Pol_total/(1/p.k_in+r.sum()+rhe.sum())
        return r*flux, rhe*flux, averages, h, flux/p.k_in

    def residual(ef):
        r, rhe, averages, _, _ = populations(ef)
        return ef + r@averages[:, 1] + rhe@averages[pas:, 1] - p.E_total

    ef = brentq(residual, 0, p.E_total, xtol=1e-9) if p.E_total else 0.
    r, rhe, averages, h, pf = populations(ef)
    return EquilibriumSolution(r, rhe, ef, pf, averages[:, 1], averages[:, 0], h, residual(ef))
