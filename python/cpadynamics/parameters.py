"""Parameters matching the MATLAB baseline supplied on 2026-09-15."""
from dataclasses import dataclass, fields
import math


@dataclass(frozen=True)
class Parameters:
    L_a: float = 100.0
    geneLength_bp: float = 25000.0
    PASposition: float = 20000.0
    k_in: float = 2.0
    k_e: float = 0.65
    k_e2: float = 0.30
    E_total: float = 100000.0
    Pol_total: float = 70000.0
    kEon: float = 2.5e-6
    kEoff: float = 0.5
    kHon: float = 7.0  # base rate PER E; never overwritten by an effective rate
    kHoff: float = 1.0
    kc: float = 0.15
    kPon_min: float = 0.01
    kPon_slope: float = 0.005
    kPoff: float = 1.0
    # Full model only: CTD detachment of the H-engaged E, followed by rapid EH
    # disassembly.
    kEoff_engaged: float = 0.05

    def __post_init__(self):
        for field in fields(self):
            value = getattr(self, field.name)
            if not math.isfinite(value) or value < 0:
                raise ValueError(f'{field.name} must be finite and nonnegative')
        if self.L_a <= 0:
            raise ValueError('L_a must be positive')
        n, pas, _ = self.geometry()
        if n < 1 or not 0 <= pas < n:
            raise ValueError('Geometry must contain at least one node and a PAS inside the gene')

    def geometry(self):
        """MATLAB geometry, with PAS converted to a ZERO-based index."""
        n = math.floor(self.geneLength_bp / self.L_a)
        pas = math.floor(self.PASposition / self.L_a) - 1
        return n, pas, n - pas


def working_parameters():
    """Return the shared MATLAB/Python default parameter set."""
    return Parameters()
