"""Full finite-rate gene-length analysis with shared genome-wide resources.

Forward transport makes the steady operator block triangular in position.
Solving each local R/RHE block in order is exactly the full sparse solve;
transport remains in every block. Upstream profiles can be shared between
genes with different PAS positions. No binding equilibrium or occupancy
interpolation is used.
"""
from dataclasses import dataclass

import numpy as np
from scipy.linalg import block_diag, lu_factor, lu_solve
from scipy.optimize import brentq
from scipy.special import ndtr

from .states import state_map, internal_components


@dataclass
class LengthProfiles:
    pol_per_flux: np.ndarray
    e_per_flux: np.ndarray
    cleavage_cdf: np.ndarray
    distances_bp: np.ndarray
    cad50_bp: np.ndarray
    max_flux_residual: float


def length_population(spacing=100, minimum=2500, maximum=200000,
                      median=22000, log10_sigma=0.68):
    """Integrate the truncated lognormal exactly over each lattice-length bin.

    With a downstream extension that is a multiple of spacing, all lengths
    in [j*spacing,(j+1)*spacing) have identical discretized geometry. Include
    the upper endpoint for plotting, with zero probability mass.
    """
    if (spacing <= 0 or minimum < spacing or maximum <= minimum or median <= 0
            or log10_sigma <= 0 or minimum % spacing or maximum % spacing):
        raise ValueError('Require positive, lattice-aligned length limits and distribution')
    lengths = np.arange(minimum, maximum + spacing, spacing, dtype=float)
    cdf = ndtr(np.log(lengths / median) / (log10_sigma * np.log(10)))
    mass = cdf[-1] - cdf[0]
    weights = np.r_[np.diff(cdf) / mass, 0.]
    return lengths, weights, float(mass)


class GeneLengthModel:
    """Conditional resource demand per unit initiation flux for many lengths."""

    def __init__(self, parameters, m, lengths, downstream_bp=5000):
        self.p = p = parameters
        self.lengths = np.asarray(lengths, dtype=float)
        if (self.lengths.ndim != 1 or not self.lengths.size
                or not np.isfinite(self.lengths).all()
                or np.any(self.lengths < p.L_a)
                or np.any(self.lengths % p.L_a)
                or downstream_bp < 0 or downstream_bp % p.L_a
                or min(p.k_in, p.k_e, p.k_e2) <= 0):
            raise ValueError('Require lattice-aligned lengths and positive transport/initiation')
        self.pas = (self.lengths / p.L_a).astype(int) - 1
        self.post_nodes = int(downstream_bp / p.L_a) + 1
        nodes = int(self.pas.max()) + self.post_nodes
        r_states, h_states = state_map(m), state_map(m, True)
        self.sr, self.sh = len(r_states), len(h_states)
        self.er = np.array(r_states)[:, 1]
        self.eh = np.array(h_states)[:, 1]
        self.transport = np.r_[np.full(self.sr, p.k_e), np.full(self.sh, p.k_e2)]
        cross = np.zeros((self.sr + self.sh,) * 2)
        lookup = {state: i for i, state in enumerate(r_states)}
        for j, (phos, e) in enumerate(h_states):
            r, h = lookup[(phos, e)], self.sr + j
            for source, target, rate in (
                    (r, h, e * p.kHon), (h, r, p.kHoff),
                    (h, lookup[(phos, e-1)], p.kEoff_engaged)):
                cross[target, source] += rate
                cross[source, source] -= rate
        self.r0, self.h0 = [], []
        for node in range(nodes):
            kpon = p.kPon_min + p.kPon_slope * node
            r0, re = internal_components(m, kpon, p)
            h0, he = internal_components(m, kpon, p, True)
            self.r0.append(p.k_e * np.eye(self.sr) - r0)
            self.h0.append(np.diag(self.transport) - block_diag(r0, h0) - cross
                           + np.diag(np.r_[np.zeros(self.sr), np.full(self.sh, p.kc)]))
        self.re, self.he = re, block_diag(re, he)

    def evaluate(self, e_free):
        """Return finite-rate occupancies and CDFs at prescribed free E."""
        if not np.isfinite(e_free) or e_free < 0:
            raise ValueError('Free E must be finite and nonnegative')
        p = self.p
        last_pas = int(self.pas.max())
        upstream = np.zeros((last_pas + 1, self.sr))
        incoming = np.zeros(self.sr)
        incoming[0] = 1.
        for node in range(last_pas):
            upstream[node + 1] = lu_solve(lu_factor(self.r0[node] - e_free*self.re), incoming)
            incoming = p.k_e * upstream[node + 1]
        # Row j holds the occupancy strictly before PAS index j.
        prefix_pol = np.cumsum(upstream.sum(axis=1))
        prefix_e = np.cumsum(upstream @ self.er)
        pol, bound_e = prefix_pol[self.pas].copy(), prefix_e[self.pas].copy()
        cleavage = np.zeros((len(self.lengths), self.post_nodes))
        terminal = np.zeros(len(self.lengths))
        factors = {}
        for i, pas in enumerate(self.pas):
            incoming = np.r_[p.k_e * upstream[pas], np.zeros(self.sh)]
            if pas == 0:
                incoming[0] = 1.
            for j in range(self.post_nodes):
                node = int(pas + j)
                if node not in factors:
                    factors[node] = lu_factor(self.h0[node] - e_free*self.he)
                local = lu_solve(factors[node], incoming)
                if not np.isfinite(local).all() or local.min() < -1e-9:
                    raise RuntimeError('Invalid conditional microstate population')
                pol[i] += local.sum()
                bound_e[i] += self.er @ local[:self.sr] + self.eh @ local[self.sr:]
                cleavage[i, j] = p.kc * local[self.sr:].sum()
                incoming = self.transport * local
            terminal[i] = incoming.sum()
        total_out = cleavage.sum(axis=1) + terminal
        residual = float(np.max(np.abs(total_out - 1)))
        if residual > 1e-8:
            raise RuntimeError(f'Conditional flux balance failed: {residual}')
        cdf = np.column_stack((np.zeros(len(pol)), np.cumsum(cleavage, axis=1) / total_out[:, None]))
        distances = np.arange(self.post_nodes + 1) * p.L_a
        cad = np.full(len(pol), np.nan)
        for i in np.flatnonzero(cdf[:, -1] >= 0.5):
            j = np.searchsorted(cdf[i], 0.5, side='left')
            cad[i] = distances[j-1] + p.L_a * (0.5-cdf[i, j-1]) / (cdf[i, j]-cdf[i, j-1])
        return LengthProfiles(pol, bound_e, cdf, distances, cad, residual)

    def solve_shared_pools(self, weights, num_genes=10000):
        """Conserve global pools; local initiation is k_in*R_free per gene."""
        weights = np.asarray(weights, dtype=float)
        if (weights.shape != self.lengths.shape or not np.isfinite(weights).all()
                or np.any(weights < 0) or not np.isclose(weights.sum(), 1, atol=1e-12)
                or num_genes < 1 or int(num_genes) != num_genes):
            raise ValueError('Require normalized nonnegative weights and positive gene count')
        p = self.p
        calls = 0

        def evaluate(ef):
            nonlocal calls
            calls += 1
            profile = self.evaluate(ef)
            rf = p.Pol_total / (1 + num_genes*p.k_in*(weights @ profile.pol_per_flux))
            e_bound = num_genes*p.k_in*rf*(weights @ profile.e_per_flux)
            return ef + e_bound - p.E_total, rf, profile

        ef = brentq(lambda e: evaluate(e)[0], 0, p.E_total, xtol=1e-8) if p.E_total else 0.
        e_residual, rf, profile = evaluate(ef)
        pol_bound = num_genes*p.k_in*rf*(weights @ profile.pol_per_flux)
        pol_residual = rf + pol_bound - p.Pol_total
        if max(abs(e_residual), abs(pol_residual)) > 1e-6:
            raise RuntimeError('Shared genome pools failed conservation')
        return dict(E_free=float(ef), Pol_free=float(rf),
                    E_conservation_residual=float(e_residual),
                    Pol_conservation_residual=float(pol_residual),
                    conditional_evaluations=calls), profile
