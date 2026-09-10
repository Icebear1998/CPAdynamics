"""Generate fixed (p,e) state maps and column-oriented kinetic matrices.

The R transitions match MATLAB at positive rates. Unlike the old scratch-matrix
pruning, state dimensions remain fixed when rates or free E vanish.
"""
from numbers import Integral
import numpy as np


def state_map(m, recognized=False):
    if isinstance(m, bool) or not isinstance(m, Integral) or m < 1:
        raise ValueError('M must be a positive integer')
    return tuple((p, e) for p in range(m + 1)
                 for e in range(1 if recognized else 0, p + 1))


def internal_components(m, kpon, parameters, recognized=False):
    """Return A0, AE such that A(E_free) = A0 + E_free * AE.

    Phosphorylation rates have no site multiplicity, matching current MATLAB.
    In RHE only e-1 bound E factors use ordinary kEoff. The engaged E's separate
    detachment/disassembly route changes RHE to R and is built in full_model.py.
    """
    states = state_map(m, recognized)
    index = {state: i for i, state in enumerate(states)}
    a0 = np.zeros((len(states), len(states)))
    ae = np.zeros_like(a0)

    def add(matrix, source, target, rate):
        j, i = index[source], index[target]
        matrix[i, j] += rate
        matrix[j, j] -= rate

    for p, e in states:
        source = (p, e)
        if p < m:
            add(a0, source, (p + 1, e), kpon)
        if p > e:
            add(a0, source, (p - 1, e), parameters.kPoff)
            add(ae, source, (p, e + 1), (p - e) * parameters.kEon)
        exchangeable = e - int(recognized)
        if exchangeable > 0:
            add(a0, source, (p, e - 1), exchangeable * parameters.kEoff)
    return a0, ae


def internal_generator(m, kpon, e_free, parameters, recognized=False):
    a0, ae = internal_components(m, kpon, parameters, recognized)
    return a0 + e_free * ae
