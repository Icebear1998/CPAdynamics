"""Flux-based cleavage CDF, including PAS cleavage and terminal run-off."""
from dataclasses import dataclass
import numpy as np


@dataclass
class CleavageProfile:
    distances_bp: np.ndarray  # includes origin
    cdf: np.ndarray           # includes zero at origin
    cad50_bp: float
    max_exit_cdf: float
    total_outflux: float


def cleavage_profile(R, RHE, parameters):
    """First RHE bin is (0,L_a]; never extrapolate an unreached median.

    This matches the corrected MATLAB coordinates and fluxes. For a plateau,
    use the earliest crossing with interpolation from the preceding point;
    MATLAB's legacy plateau branch instead snaps to a bin edge.
    """
    r, rhe = np.asarray(R, dtype=float), np.asarray(RHE, dtype=float)
    if r.ndim != 1 or rhe.ndim != 1 or min(r.size, rhe.size) < 1:
        raise ValueError('R and RHE must be nonempty one-dimensional profiles')
    if not (np.isfinite(r).all() and np.isfinite(rhe).all()):
        raise ValueError('Profiles must be finite')
    if min(r.min(), rhe.min()) < -1e-9:
        raise ValueError('Profiles must be nonnegative')
    flux = parameters.kc * rhe
    total = float(flux.sum() + parameters.k_e*r[-1] + parameters.k_e2*rhe[-1])
    cdf = np.r_[0., np.cumsum(flux)/total if total > 1e-9 else np.zeros_like(rhe)]
    distances = np.arange(rhe.size + 1, dtype=float) * parameters.L_a
    cad = float('nan')
    crossing = np.flatnonzero(cdf >= .5)
    if crossing.size:
        j = crossing[0]
        cad = float(distances[j-1] + (distances[j]-distances[j-1])
                    * (.5-cdf[j-1])/(cdf[j]-cdf[j-1]))
    return CleavageProfile(distances, cdf, cad, float(cdf[-1]), total)
