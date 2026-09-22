# Python full R/RHE model

This directory implements the full finite-rate model discussed in **Assess time
separation assumption**. It is independent of the older, git-ignored `pyCPA/`
package and does not modify MATLAB or the legacy analyses.

## Run

To reproduce the full-model figures and diagnostics from `CpaMultipleEMain.m`:

```bash
python3 python/run_cpa_multiple_e_main.py
```

This uses `EBindingNumber=5` and writes E/Ser2P and R/RHE plots, CAD50,
cleavage CDF, raw profiles, microstates and a manifest under
`SecondVersionResults/python_CPA_multipleE_main/<timestamp>/`. The profiles
span the full simulated gene and use position relative to PAS, matching the
MATLAB script. Pass `--m N` or `--output-dir PATH` to change those settings.

To run gene length versus CAD with the full finite-rate model:

```bash
python3 python/plot_gene_length_vs_cad.py
```

This implements the shared-pool analysis in `GeneLengthGenerateGrid.m`,
`GeneLengthBuildInterpolation.m`, and `GeneLengthAnalyze.m`. It uses the current
working parameters, capacity 5, 10,000 active genes, per-gene `k_in=2/10000`,
TSS-to-PAS lengths 2.5–200 kb, and a fixed 5 kb downstream extension. The
lognormal length distribution has median 22 kb and log10 standard deviation
0.68, normalized over the simulated range.

The Python implementation directly solves the full kinetic equations at every
100 bp lattice length. Forward block substitution reuses upstream states
without changing any reactions or transport rates. Exact lognormal probability
differences integrate each lattice-length bin, avoiding the MATLAB pipeline's
occupancy interpolation and coarse quadrature. It solves one shared pair of
free pools and checks representative genes against the complete sparse model.
This is equation parity with MATLAB, not an expectation of identical output
from a coarse MATLAB interpolation grid. MATLAB execution is not required.

The command saves PNG/SVG figures, all 1,976 CAD/occupancy rows, a selected-length
CSV, cleavage CDFs, a report, and a manifest with parameters, source hashes and
validation under `SecondVersionResults/python_GeneLengthAnalysis/<timestamp>/`.
Use `--m 1` to change capacity, `--output-dir PATH` for a chosen directory, or
`CPAD_RESULTS_ROOT` to relocate results. CAD is never extrapolated. Matching
MATLAB's PAS-bin convention, the 5 kb extension contains 51 cleavage bins with
CDF coordinates through 5.1 kb. The upper endpoint at 200 kb is plotted but has
zero mass under the continuous length distribution.

To reproduce the current `PlotEBindingProfile.m` with the full finite-rate model:

```bash
python3 python/plot_ebinding_profile.py
```

This runs capacities 1, 5 and 10 and plots average E (solid) and Ser2P (dashed)
from TSS to PAS (inclusive), with distance from TSS in kb on the x-axis.
Both averages include R and RHE populations. The full gene is still simulated;
raw CSVs retain the full profiles and their coordinates relative to PAS.
All analysis entry points use the same defaults as `default_parameters.m`:
`kEon=2.5e-6`, `kEoff=0.5`, `kEoff_engaged=0.05`, `kHon=7`,
`kHoff=1`, `kc=0.15`, `k_e=0.65`, and `k_e2=0.3`. The profile command saves
PNG/SVG figures, profile CSVs, microstates and a parameter/diagnostics manifest
under `SecondVersionResults/python_EBindingProfile/<timestamp>/` and prints the
absolute path. Use `--output-dir PATH` or `CPAD_RESULTS_ROOT` to change its location.

From the repository root, using Python 3.10 or later:

```bash
python3 -m pip install -r python/requirements.txt
python3 python/run_comparison.py --m 1 5
PYTHONPATH=python python3 -m unittest discover -s python/tests -v
```

Using explicit historical inputs (`kHon=4`, `kHoff=2`, `kc=0.13`,
`kEoff_engaged=0`), the comparison first verifies the equilibrium reference against all six supplied
MATLAB CAD50 values: **1191, 728, 571, 490, 442, 407 bp** for M=1..6. Each must
agree within 0.5 bp, the precision of the supplied numbers, before full-model
calculations run. MATLAB itself was unavailable during implementation; this is
regression against supplied MATLAB outputs, not a new MATLAB execution.

Results go to a new timestamped directory under
`SecondVersionResults/python_full_model/`. That generated-results directory is
git-ignored by the existing repository configuration; the code here is not.
The run prints the absolute output directory. To plot its saved CSVs:

```bash
python3 python/plot_comparison.py /absolute/path/to/the/results-directory
```

Generate all seven cases with `--m 1 2 3 4 5 6 7`. The same matrix generator
constructs each system. The default comparison focuses on M=1 and M=5.

To reproduce the six-point `EBindingNumberVsCad.m` sweep and compare it with
the full model:

```bash
python3 python/run_comparison.py --m 1 2 3 4 5 6
python3 python/plot_ebinding_sweep.py /absolute/path/to/the/results-directory
```

This creates `ebinding_number_vs_cad.png` and `.svg` from `comparison.csv`,
with maximum capacity on the x axis and CAD50 on the y axis. It uses the current
MATLAB script's six capacities, default parameters, and 50% threshold.

For an additional time-integration comparison starting with an empty gene:

```bash
python3 python/run_comparison.py --m 1 5 --transient-seconds 30000
```

This is slower than the direct steady solve. A successful integration is not
automatically a steady state: `dynamic_validation.csv` reports both the final
kinetic residual and agreement with the direct solve. The closed polymerase pool
can produce long transients after synchronized loading of an initially empty
gene. A 5,000-second run did not converge for the default gene.

During implementation, both default cases were independently verified from an
empty gene: M=1 at 30,000 s, and M=5 at 50,000 s (continuing its saved 30,000 s
state for another 20,000 s). Final maximum species-scaled deviations from the
direct solution were below 1e-6, and kinetic residuals were below 1e-4.
The M=5 dynamic CAD50 was 551.4948 bp versus 551.4944 bp from the direct solve.
Use `--transient-seconds 50000` to run a fresh trajectory over that longer interval.

## API

Set `PYTHONPATH=python`, or run from the `python` directory:

```python
from dataclasses import replace
from cpadynamics import Parameters, FullModel, cleavage_profile

parameters = Parameters()
# For a different H off-rate, for example:
# parameters = replace(parameters, kHoff=0.05)
# Independently allow engaged-E CTD detachment followed by rapid EH disassembly:
# parameters = replace(parameters, kEoff_engaged=0.05)
model = FullModel(parameters, m=5)
solution = model.solve_steady_state()
profile = cleavage_profile(solution.R, solution.RHE, parameters)
print(profile.cad50_bp, solution.E_free, solution.rhs_max_abs)

# Optional transient simulation from an empty gene:
trajectory = model.integrate(30000, t_eval=[0, 100, 1000, 10000, 30000])
final = model.summarize(trajectory.y[:, -1])
```

`Parameters.kHon` is always the base per-E recognition rate. Unlike the MATLAB
return struct it is never overwritten with an occupancy-dependent effective rate.
If changing `L_a`, also set `k_e=65/L_a` and `k_e2=30/L_a` to retain the default
physical elongation velocities. Parameters are immutable; use `replace` for
changes. The comparison CLI intentionally uses the exact default parameters.

## State definitions and reactions

- `R[i,p,e]`: unrecognized polymerase, `0 <= e <= p <= M`.
- `RHE[i,p,e]`: recognized polymerase, `1 <= e <= p <= M`.
- One of the `e` E factors in each RHE is H-engaged and excluded from ordinary
  E unbinding. A separate `kEoff_engaged` can allow its CTD detachment.
- M changes both the maximum phosphorylation level and maximum E occupancy.
- Initiation enters `R[0,0,0]` at `k_in * Pol_free`.

Internal transitions follow these rates per polymerase:

| Transition | R | RHE |
|---|---:|---:|
| Phosphorylation, p -> p+1 if p<M | kPon(i) | kPon(i) |
| Dephosphorylation, p -> p-1 if p>e | kPoff | kPoff |
| E binding, e -> e+1 if e<p | (p-e) kEon E_free | (p-e) kEon E_free |
| Ordinary E unbinding | e kEoff | (e-1) kEoff |

The phosphorylation transitions have **no site multiplicity**, matching the
current MATLAB model. Occupied phosphorylation sites are protected from
dephosphorylation; the H-engaged site is therefore protected too.

At every node from the PAS onward:

```text
R(p,e)   -> RHE(p,e)    at e*kHon
RHE(p,e) -> R(p,e)      at kHoff
RHE(p,e) -> R(p,e-1) + free E at kEoff_engaged
RHE(p,e) -> free pools at kc
```

Recognition selects one of the `e` factors to engage H and preserves `(p,e)`.
H dissociation preserves `(p,e)` too: it removes protection but does not itself
release E. Independently, the H-engaged E may detach from the CTD at
`kEoff_engaged`. The detached EH is assumed to disassemble rapidly, so the net
event returns exactly one E to the free pool, preserves p, and returns the
polymerase to unrecognized R(p,e-1). This applies even when other E remain:
those E can subsequently bind H again. There is only one engaged E, so this
new rate has **no e multiplier**. H is the PAS and does not become a separate
free-protein pool; the event removes recognition, not the transcript.

There is no RHE state with e=0. When the last E detaches, the destination is
R(p,0), not RHE(p,0). Additional E can bind and unbind normally while the complex
is recognized. The default `kEoff_engaged=0.05` enables this route in both
the steady and time solvers; setting it to zero disables engaged-E loss.

R elongates at `k_e`, RHE at `k_e2`, preserving `(p,e)`. Cleavage and final-node
run-off return one polymerase and all its bound E to the free pools. H denotes
PAS engagement; this implementation does not introduce a separate finite H pool.

At position i (zero-based), `kPon(i) = kPon_min + kPon_slope*i`.
Geometry matches MATLAB exactly: `N=floor(geneLength_bp/L_a)`,
`PAS=floor(PASposition/L_a)-1`, and `N_PAS=N-PAS`.

## Matrix formulation and solvers

`states.py` generates fixed state maps and local matrices. It keeps all allowed
states even when free E or a rate becomes zero, avoiding the old builder's
rate-dependent row/column pruning. For positive rates, its R matrix matches the
literal MATLAB builder (checked for M=1..7 at four phosphorylation rates).

`full_model.py` assembles two sparse matrices:

```text
A(E_free) = A0 + E_free * AE
dx/dt = A(E_free) x + k_in * Pol_free * b
dE_free/dt = -e_counts @ dx/dt
dPol_free/dt = -sum(dx/dt)
```

Here `x` contains all bound polymerase microstates and `b` injects into R(0,0)
at the first node. The actual pool rows are precomputed from these conservation
identities. **No local equilibrium distribution appears in this system.**

The state vector is all R nodes (node-major), then all RHE nodes (node-major),
then E_free and Pol_free. Within each node, `(p,e)` is ordered by p, then e.
`R_micro` and `RHE_micro` return these arrays as node-by-state matrices.

The number of bound-state variables is

```text
N*(M+1)*(M+2)/2 + N_PAS*M*(M+1)/2.
```

For the default geometry, M=1 has 801 variables, M=5 has 6,015, and M=7 has
10,428. Time integration adds two explicit free-pool variables.

Two solvers evaluate the **same full equations**:

1. **Direct steady state:** at a candidate E_free, solve `A u = -b` with sparse
   linear algebra. This is the occupation per unit initiation flux. Set
   `J=Pol_total/(1/k_in + sum(u))`, `x=J*u`, `Pol_free=J/k_in`; solve the scalar
   equation `E_free + e_counts @ x = E_total`. This is exact steady-state
   elimination, not a rapid-binding approximation. Every finite-rate transport
   and microstate equation remains in A. Residual checks reject invalid results.
2. **Time integration:** SciPy BDF with an analytic sparse Jacobian evolves
   all microstates and both pools together. By default, the gene starts empty,
   with both pools entirely free. Supplying `t_eval` limits stored output.

The equations are deterministic mass-action population equations. They remove
the fast-binding closure; they do not model stochastic fluctuations of finite
global pools, polymerase exclusion, or additional biochemical intermediates.

Solver references: [SciPy BDF](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.BDF.html)
and [SciPy sparse linear solve](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.spsolve.html).

## Cleavage convention and baseline parity

Cleavage is active in the first PAS-associated state and counted in its loss
term and its flux. That state represents the interval `(0,100] bp` at default
spacing; its cumulative flux is reported at 100 bp. CDF arrays explicitly
include `(0 bp, 0%)`. Normalization includes cleavage and both terminal run-off
fluxes, so cleavage fractions need not reach 100% inside the window.

An unreached 50% threshold returns NaN, never an extrapolated CAD. At a plateau,
the Python helper consistently interpolates the first threshold crossing.
The current MATLAB helper has a legacy branch that snaps to a bin edge when
any CDF increment is exactly zero. This edge-case difference does not affect
the six default reference results.

`equilibrium.py` retains the MATLAB model's 100-point E_free grid, linear
interpolation, PAS-frozen `kHon*mean_E(PAS)`, aggregate transport, and E pool
conservation. It solves stationary distributions by a normalized linear solve
instead of SVD, and the linear transport equations directly instead of fsolve.
These are equivalent equations at the default positive rates. The outer root
is solved to tighter tolerance than the older Python port.

With the historical regression inputs above, the full model predicts approximately **1535.654 bp (M=1)** and
**551.494 bp (M=5)**, versus **1191.205** and **441.784 bp** from the equilibrium
reference. This comparison changes finite-rate loading, protection in RHE, and
recognition based on actual occupancy. It does **not** isolate the error caused
solely by the timescale approximation.

## Validation and artifacts

Tests cover all six supplied MATLAB values, M=1 analytic internal equations,
protection of the last E, H dissociation without E release, cleavage/resource
release in the PAS bin, zero E, conservation away from steady state, the analytic
Jacobian against finite differences, physical steady states for M=1..7, and
agreement between BDF integration and the direct solver on a smaller gene for
M=1 and M=5. No MATLAB, FirstVersion, SparedCodes, or pyCPA source is modified.

The engaged-E extension adds tests of every RHE-to-R detachment mapping for
M=1 and M=5, the independent H-first and E-first pathways, zero-rate backward
compatibility, conservation at nonzero rates for M=1..6, and time integration
at a nonzero rate. The complete suite contains 18 tests.

Each comparison saves:

- `matlab_parity.csv`, `comparison.csv`, and `REPORT.md`.
- `cleavage_M*.csv` with origin-inclusive CDFs.
- `profiles_M*.csv` with aggregate populations and conditional occupancies.
- `microstates_M*.npz` with the complete state populations and maps.
- `manifest.json` with exact parameters, versions, and source SHA-256 hashes.
- Optional `dynamic_validation.csv`, `transient_M*.csv`, and final transient states.

Plotting is separate: `plot_comparison.py` reads the saved CSVs and writes PNG
and SVG figures without rerunning the model.

## Engaged-E detachment sensitivity

```bash
python3 python/run_engaged_e_sweep.py
python3 python/plot_engaged_e_sweep.py /absolute/path/to/the/results-directory
```

By default this runs M=1..6 at `kEoff_engaged/kEoff = 0, 0.1, 1`, corresponding
to 0, 0.05 and 0.5 s^-1 with the current defaults. These are illustrative
sensitivity values, not experimentally estimated rates. Adjust the sweep with
`--m 1 5` or `--relative-off-rates 0 0.01 0.1 1`.

Results are saved under `SecondVersionResults/python_engaged_E_dissociation/`
with the parameters, source hashes, CSVs, microstate NPZs, and a report. The
original MATLAB-equilibrium reference does not include the new reaction and is
calculated with the original parameters. All six supplied MATLAB values are
checked before running the new sensitivity calculation.

At 0.05 s^-1, CAD50 is approximately 1698.670 bp for M=1 and 581.572 bp for M=5,
versus 1535.654 and 551.494 bp when the new route is disabled. In this tested
regime the extra loss pathway lengthens CAD, moving it farther from the
original equilibrium reference. This reflects the chosen rapid-disassembly
mechanism and does not establish a measured biological off-rate.

The current default full model gives **CAD50 = 355.4668 bp at M=5**,
matching the MATLAB output supplied on 2026-09-15. Historical numerical results
above retain their original parameter inputs and are not current-default predictions.
