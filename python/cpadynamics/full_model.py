"""Full finite-rate R/RHE kinetics, with one H-engaged E per RHE complex.

The engaged E is excluded from ordinary E unbinding. Its optional, separate
CTD-detachment route (kEoff_engaged) immediately disassembles EH, returns one E
to the free pool, and leaves an unrecognized R with one fewer bound E.

State ordering: all R nodes (node-major), all RHE nodes (node-major), free E,
free Pol II. Internal state maps are ordered by p, then e. No local equilibrium
distribution or occupancy interpolant enters these equations.
"""
from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.sparse.linalg import spsolve

from .states import state_map, internal_components


@dataclass
class FullSolution:
    state: np.ndarray
    R_micro: np.ndarray
    RHE_micro: np.ndarray
    R: np.ndarray
    RHE: np.ndarray
    E_free: float
    Pol_free: float
    e_residual: float
    pol_residual: float
    rhs_max_abs: float


class FullModel:
    def __init__(self, parameters, m):
        self.parameters = p = parameters
        self.m = m
        self.r_states = state_map(m)
        self.rhe_states = state_map(m, recognized=True)
        self.n, self.pas, self.n_pas = p.geometry()
        self.sr, self.sh = len(self.r_states), len(self.rhe_states)
        self.r_size = self.n*self.sr
        self.size = self.r_size + self.n_pas*self.sh
        self.e_counts = np.r_[np.tile(np.array(self.r_states)[:, 1], self.n),
                              np.tile(np.array(self.rhe_states)[:, 1], self.n_pas)]
        self.unit_source = np.zeros(self.size)
        self.unit_source[0] = 1  # initiation is exclusively into R(0,0)

        components_r = [internal_components(m, p.kPon_min+p.kPon_slope*i, p)
                        for i in range(self.n)]
        components_rhe = [internal_components(m, p.kPon_min+p.kPon_slope*i, p, True)
                          for i in range(self.pas, self.n)]
        all_components = components_r + components_rhe
        a0 = sparse.block_diag([sparse.csr_matrix(c[0]) for c in all_components], format='csc')
        self.ae = sparse.block_diag([sparse.csr_matrix(c[1]) for c in all_components], format='csc')

        # Transport preserves (p,e). Final-node outflow leaves the lattice.
        r_transport = self._transport(self.n, self.sr, p.k_e)
        rhe_transport = self._transport(self.n_pas, self.sh, p.k_e2)
        rhe_transport -= p.kc*sparse.eye(self.n_pas*self.sh, format='csc')
        a0 += sparse.block_diag([r_transport, rhe_transport], format='csc')

        # R(p,e) -> RHE(p,e) selects one of e bound E; H off preserves counts.
        r_lookup = {state: i for i, state in enumerate(self.r_states)}
        rows, cols, rates = [], [], []
        for j in range(self.n_pas):
            for h_index, state in enumerate(self.rhe_states):
                r = (self.pas+j)*self.sr+r_lookup[state]
                h = self.r_size+j*self.sh+h_index
                recognition = state[1]*p.kHon
                rows.extend([h, r, r, h])
                cols.extend([r, r, h, h])
                rates.extend([recognition, -recognition, p.kHoff, -p.kHoff])
                # One engaged E, hence NO e multiplicity on this route.
                # Phosphorylation p is preserved; recognition is lost even if
                # other E remain. Those factors can subsequently recognize H.
                r_after_loss = (self.pas+j)*self.sr+r_lookup[(state[0], state[1]-1)]
                rows.extend([r_after_loss, h])
                cols.extend([h, h])
                rates.extend([p.kEoff_engaged, -p.kEoff_engaged])
        a0 += sparse.coo_matrix((rates, (rows, cols)), shape=(self.size, self.size)).tocsc()
        self.a0 = a0.tocsc()
        self.a0.eliminate_zeros()
        self.ae.eliminate_zeros()

        # Pool rows account for association, dissociation, cleavage and run-off.
        weights = sparse.csr_matrix(self.e_counts[None, :])
        self.e_release0 = -np.asarray((weights@self.a0).toarray()).ravel()
        self.e_release1 = -np.asarray((weights@self.ae).toarray()).ravel()
        self.pol_release0 = -np.asarray(self.a0.sum(axis=0)).ravel()
        self.pol_release1 = -np.asarray(self.ae.sum(axis=0)).ravel()
        self.initiation = p.k_in*self.unit_source

    @staticmethod
    def _transport(nodes, states, rate):
        size = nodes*states
        out = -rate*sparse.eye(size, format='csc')
        if nodes > 1:
            out += sparse.diags(rate*np.ones(size-states), -states, shape=(size, size), format='csc')
        return out

    def initial_state(self):
        """Empty gene; all E and polymerase initially free."""
        return np.r_[np.zeros(self.size), self.parameters.E_total, self.parameters.Pol_total]

    def rhs(self, time, state):
        """Simultaneous mass-action kinetics including explicit free pools."""
        x, ef, pf = state[:self.size], state[-2], state[-1]
        dx = self.a0@x + ef*(self.ae@x) + pf*self.initiation
        de = self.e_release0@x + ef*(self.e_release1@x)
        dp = self.pol_release0@x + ef*(self.pol_release1@x) - self.parameters.k_in*pf
        return np.r_[dx, de, dp]

    def jacobian(self, time, state):
        """Analytic sparse Jacobian; explicit pools avoid a dense rank-one update."""
        x, ef = state[:self.size], state[-2]
        return sparse.bmat([
            [self.a0+ef*self.ae, sparse.csc_matrix((self.ae@x)[:, None]),
             sparse.csc_matrix(self.initiation[:, None])],
            [sparse.csc_matrix((self.e_release0+ef*self.e_release1)[None, :]),
             sparse.csc_matrix([[self.e_release1@x]]), sparse.csc_matrix((1, 1))],
            [sparse.csc_matrix((self.pol_release0+ef*self.pol_release1)[None, :]),
             sparse.csc_matrix([[self.pol_release1@x]]),
             sparse.csc_matrix([[-self.parameters.k_in]])],
        ], format='csc')

    def summarize(self, state):
        y = np.asarray(state, dtype=float)
        if y.shape != (self.size+2,) or not np.isfinite(y).all():
            raise ValueError('State has incorrect shape or nonfinite entries')
        r = y[:self.r_size].reshape(self.n, self.sr)
        rhe = y[self.r_size:self.size].reshape(self.n_pas, self.sh)
        return FullSolution(
            y.copy(), r.copy(), rhe.copy(), r.sum(axis=1), rhe.sum(axis=1),
            float(y[-2]), float(y[-1]),
            float(y[-2]+self.e_counts@y[:self.size]-self.parameters.E_total),
            float(y[-1]+y[:self.size].sum()-self.parameters.Pol_total),
            float(np.max(np.abs(self.rhs(0, y)))))

    def solve_steady_state(self):
        """Solve ALL finite-rate equations at steady state, without time separation.

        At fixed free E, bound-state equations are linear. Solve them per unit
        initiation flux, enforce Pol conservation by scaling that flux, then
        find free E from E conservation. This is an exact algebraic reduction
        of the full steady equations, not a binding-equilibrium approximation.
        """
        p = self.parameters
        if min(p.k_e, p.k_e2, p.k_in) <= 0:
            raise ValueError('Steady solver requires positive initiation and elongation')

        def conditional(ef):
            occupation = spsolve(-(self.a0+ef*self.ae), self.unit_source)
            if not np.isfinite(occupation).all() or occupation.min() < -1e-9:
                raise RuntimeError('Conditional steady solve produced an invalid population')
            flux = p.Pol_total/(1/p.k_in+occupation.sum())
            return occupation*flux, flux/p.k_in

        def residual(ef):
            x, _ = conditional(ef)
            return ef+self.e_counts@x-p.E_total

        ef = brentq(residual, 0, p.E_total, xtol=1e-9) if p.E_total else 0.
        x, pf = conditional(ef)
        solution = self.summarize(np.r_[x, ef, pf])
        if max(abs(solution.e_residual), abs(solution.pol_residual), solution.rhs_max_abs) > 1e-6:
            raise RuntimeError('Full steady solution failed conservation or kinetic residual checks')
        return solution

    def integrate(self, t_end, initial=None, rtol=1e-7, atol=1e-9, t_eval=None):
        """Integrate with SciPy BDF; returns its OdeResult (time is in seconds).

        A successful integration need not have reached steady state; inspect
        summarize(result.y[:,-1]).rhs_max_abs before calling it converged.
        Provide t_eval to retain only requested output times for large systems.
        """
        if not np.isfinite(t_end) or t_end <= 0:
            raise ValueError('t_end must be positive and finite')
        y0 = self.initial_state() if initial is None else np.asarray(initial, dtype=float)
        summary = self.summarize(y0)
        if y0.min() < 0:
            raise ValueError('Initial populations must be nonnegative')
        if max(abs(summary.e_residual), abs(summary.pol_residual)) > 1e-6:
            raise ValueError('Initial state does not conserve the configured pools')
        result = solve_ivp(self.rhs, (0, t_end), y0, method='BDF', jac=self.jacobian,
                           rtol=rtol, atol=atol, t_eval=t_eval)
        if not result.success:
            raise RuntimeError(result.message)
        if not np.isfinite(result.y).all() or result.y.min() < -max(100*atol, 1e-7):
            raise RuntimeError('Integration produced materially negative or nonfinite populations')
        return result
