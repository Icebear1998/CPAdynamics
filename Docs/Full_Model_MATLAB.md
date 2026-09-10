# Full finite-rate MATLAB model

`EBindingNumberVsCad.m` now runs the full R/RHE model with
`P.kEoff_engaged = 0.05` s⁻¹. This is a trial parameter, not a fitted estimate.
The sweep uses M = 1–6; change `EBindingNumber_values` to `1:7` or `[1,5]`
to change the cases. The state generator changes both maximum phosphorylation
and maximum E occupancy with M, without requiring new equations by hand.

## Reactions and states

At each node, R has states `(p,e)` with `0 <= e <= p <= M`. From the PAS onward,
RHE has states with `1 <= e <= p <= M`; exactly one of its E factors is engaged
with H. Here e counts all attached E factors, including the engaged one.

| Reaction | Rate per polymerase |
| --- | --- |
| Phosphorylation, p → p+1, if p < M | kPon at that node |
| Dephosphorylation, p → p−1, if p > e | kPoff |
| E binding, e → e+1 | (p−e) kEon E_free |
| Ordinary E dissociation from R | e kEoff |
| Ordinary E dissociation from RHE | (e−1) kEoff |
| R(p,e) → RHE(p,e), from PAS onward | e kHon |
| RHE(p,e) → R(p,e), H dissociates first | kHoff |
| RHE(p,e) → R(p,e−1) + free E, engaged E detaches first | kEoff_engaged |
| Cleavage from every RHE node, including the first | kc |

The last dissociation route assumes that detached EH subsequently disassembles
rapidly. We collapse that short intermediate into immediate E return to the
free pool; we do not track detached EH or a free H pool. This is a specific
approximation for disassembly. Phosphorylation, E binding, recognition and
elongation all retain their finite rates.

Phosphorylation rates have no site multiplicity, matching the existing model.
R and RHE elongate at k_e and k_e2, preserving p and e. Initiation enters
R(p=0,e=0). Cleavage and terminal run-off return polymerase and all its E.
After H dissociates first, all remaining E use ordinary kEoff and can recognize
H again. The engaged-E detachment rate is independent of kHoff and has no e
multiplicity.

## Files and solver

| File | Purpose |
| --- | --- |
| `build_full_state_map.m` | Generate R/RHE microstates for any positive integer M |
| `build_full_internal_rate_matrices.m` | Local phosphorylation and ordinary E exchange |
| `build_full_rate_matrices.m` | Assemble sparse transport, recognition, dissociation and cleavage operators |
| `ode_dynamics_full_multipleE.m` | Full ODE, including free E and free polymerase |
| `full_model_jacobian.m` | Sparse analytic Jacobian for optional ode15s integration |
| `run_full_termination_simulation.m` | Solve the full steady state with resource conservation |
| `calculate_full_pas_cleavage_profile.m` | Flux CDF and CAD, including cleavage in the first bin |

The steady-state solver fixes E_free temporarily, solves all microstate and
transport equations as one sparse linear system per unit initiation flux,
scales the solution to conserve polymerase, and uses fzero to conserve E.
It checks nonnegative populations, both resource balances, kinetic residuals
and total polymerase outflux. This is a steady-state solution of the full
ODE, not a local equilibrium closure. kHon stays the base per-E rate.

If `kEoff_engaged` is omitted, the helper defaults to zero to reproduce the
earlier full model with protected engaged E. The sweep explicitly sets 0.05.
Zero does **not** restore the original time-separation model.

The first RHE bin represents `(0,100]` bp with the default spacing. Its
cleavage contributes to the CDF at 100 bp, with an explicit `(0,0)` origin
used for CAD interpolation. An unreached threshold gives NaN and a warning;
`diagnostics.max_exit_cdf` reports the fraction reached within the window.

## Run from the repository root

```matlab
EBindingNumberVsCad

% Or run one case:
P = default_parameters();
P.kEoff_engaged = 0.05;
[R, RHE, P_out, D] = run_full_termination_simulation(P, 5);
[cdf, distance_bp, CAD50, diagnostics] = ...
    calculate_full_pas_cleavage_profile(R, RHE, P_out);
```

Optional transient integration, starting with free resources and no bound
polymerases (choose duration and tolerances for the intended experiment):

```matlab
model = D.model;
X0 = [zeros(model.size,1); P.E_total; P.Pol_total];
options = odeset('RelTol',1e-8,'AbsTol',1e-10, ...
    'Jacobian',@(t,X) full_model_jacobian(t,X,model));
[t,X] = ode15s(@(t,X) ode_dynamics_full_multipleE(t,X,model), ...
    [0,50000], X0, options);
```

`D.R_micro` and `D.RHE_micro` store node-by-state populations, and
`D.R_states` and `D.RHE_states` identify their columns. The ODE vector instead
stores each node's microstates contiguously, followed by E_free and Pol_free.

Set `CPAD_FORCE_SAVE=true` to save the sweep's MAT, CSV and PNG files under
`SecondVersionResults/Full_EBindingNumber_vs_CAD/`. The existing
`CPAD_RESULTS_ROOT` environment override is respected. Other analysis scripts
continue to use their existing helpers.

## Verification

Python reference CAD50 targets with the current defaults and engaged-E
off-rate 0.05 s⁻¹:

| M | CAD50 (bp) |
| --- | ---: |
| 1 | 1698.670 |
| 2 | 996.955 |
| 3 | 765.393 |
| 4 | 650.695 |
| 5 | 581.572 |
| 6 | 536.227 |

Run the MATLAB numerical checks with:

```matlab
results = runtests('tests/test_full_termination_model.m');
assertSuccess(results);
```

Tests cover independent dissociation routes, E accounting, first-bin cleavage,
unreached thresholds, zero E, state sizes, analytic Jacobians, conservation,
Python CAD targets, and agreement between ode15s and the direct steady solve.
For M=1 and M=5, sparse operators, derivatives and Jacobians are also checked
against `tests/fixtures/full_model_python_reference.mat`, exported independently
from the Python implementation by `python/export_matlab_full_reference.py`.

MATLAB and Octave were unavailable during implementation. MATLAB syntax checks
passed with MISS_HIT, but the MATLAB numerical tests have not been executed.
The table above contains Python reference results, not measured MATLAB output.
