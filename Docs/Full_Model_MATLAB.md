# Full finite-rate MATLAB model

`EBindingNumberVsCad.m` compares the full finite-rate R/RHE model at
`engaged_E_off_rates = [0, 0.05, 0.5]` s⁻¹. These are trial
parameters, not fitted estimates. All three curves are recomputed from the
current defaults, with a common flux-based cleavage observable.
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
| `build_full_rate_matrices.m` | Assemble sparse transport, recognition, dissociation and cleavage operators |
| `ode_dynamics_full_multipleE.m` | Full ODE, including free E and free polymerase |
| `run_full_termination_simulation.m` | Solve the full steady state with resource conservation |
| `calculate_full_pas_cleavage_profile.m` | Flux CDF and CAD, including cleavage in the first bin |

The steady-state solver fixes E_free temporarily, solves all microstate and
transport equations as one sparse linear system per unit initiation flux,
scales the solution to conserve polymerase, and uses fzero to conserve E.
It checks nonnegative populations, both resource balances, kinetic residuals
and total polymerase outflux. This is a steady-state solution of the full
ODE, not a local equilibrium closure. kHon stays the base per-E rate.

The free-E root solve uses `TolX = eps`. A looser `TolX = 1e-9` can stop
before the absolute E-conservation tolerance of `1e-6` is reached at pool
sizes near 100,000. Validation errors report each residual and the minimum
population so numerical stopping errors can be distinguished from other failures.

If `kEoff_engaged` is omitted, the helper defaults to zero to reproduce the
earlier full model with protected engaged E. The sweep explicitly sets each rate.
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

Optional transient validation using the test Jacobian helper, starting with
free resources and no bound polymerases (run from the repository root):

```matlab
addpath('tests'); % Test-only analytic Jacobian
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
`SecondVersionResults/Full_EBindingNumber_vs_CAD/`. Saved matrices have one row
per M and columns ordered as equilibrium, then each engaged-E off-rate. The
CSV has one row per (model, M, rate), including diagnostics and error messages.
The equilibrium reference uses the existing fsolve solver and therefore
requires Optimization Toolbox. The existing
`CPAD_RESULTS_ROOT` environment override is respected. Other analysis scripts
continue to use their existing helpers.

## Verification

Historical Python reference CAD50 targets with engaged-E off-rate 0.05 s⁻¹
and the pre-appendix parameters (`k_e = 0.65`, `kEon = 2.5e-6`, `kHon = 4`,
`kHoff = 2`, `kc = 0.13`; other inputs unchanged):

| M | CAD50 (bp) |
| --- | ---: |
| 1 | 1698.670 |
| 2 | 996.955 |
| 3 | 765.393 |
| 4 | 650.695 |
| 5 | 581.572 |
| 6 | 536.227 |

These values are retained as Python regression fixtures with explicit parameter inputs.
They are not predictions for the current defaults returned by
`default_parameters.m`; rerun the sweep to obtain those predictions.

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

## Migrated MATLAB analyses

All active analysis scripts now use the full finite-rate model. The common
working default is `kEoff_engaged = 0.05` s⁻¹; the binding-capacity comparison
retains its explicit three-rate sweep. Full-model details provide
`avg_E_bound`, `avg_Ser2P`, `E_bound`, and `Pol_bound` from the microstates.

For shared-resource gene-length calculations, call
`run_full_termination_simulation(P,M,'FreePools',[R_free E_free])`.
This solves a gene at externally fixed free pools without imposing a separate
copy of genome-wide resource conservation on that gene. The local steady RHS
and flux are validated; residuals against unused local totals are NaN.

Rebuild the lookup pipeline in order: `GeneLengthGenerateGrid`,
`GeneLengthBuildInterpolation`, `GeneLengthAnalyze`. New artifacts carry
`full_gene_length_` prefixes and model/pool metadata. Old equilibrium grids
are rejected. the local `solve_full_genome_pools` function in `GeneLengthAnalyze.m` exploits linear occupancy scaling in
free Pol II and solves free E on its physical interval. Interpolation is
restricted to the grid and the length population is explicitly truncated to
that grid's interval. CAD thresholds not reached within the downstream window
remain NaN and retain their within-window cleavage diagnostics.

State indexing and internal reaction generators are local functions in
`build_full_rate_matrices.m`. Gene-grid validation is local to
`GeneLengthBuildInterpolation.m`; the shared-pool solver is local to
`GeneLengthAnalyze.m`. Run the three gene-length scripts in their usual order.
Rebuild the grid and interpolation after changing defaults; saved grids retain
the parameter set used to generate them.

The analytic Jacobian lives in `tests/full_model_jacobian.m` for reference,
finite-difference and transient-integration tests. Production model builders
and steady-state analyses do not include or depend on it.

`cpad_analysis_output_dir.m` remains shared by saving and analysis entry points
so `CPAD_RESULTS_ROOT` and default output locations are handled consistently.
