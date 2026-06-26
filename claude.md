# CPAdynamics — Project Context for AI Assistants

## What This Project Is

A MATLAB research codebase modeling **RNA Polymerase II (Pol II) transcription termination** in mammalian cells. The model simulates the kinetics of **Cleavage and Polyadenylation (CPA)** to predict:

- **Alternative Polyadenylation (APA)** site choice
- **Termination Commitment Distance (TCD)**: how far downstream of a poly(A) signal (PAS) polymerases travel before terminating
- Genome-wide **resource competition** effects (how total CPA factor concentrations affect individual gene behavior)

This is Version 2.0 of the model. The key advance over Version 1.0 is that each polymerase can now bind **multiple E factors simultaneously** (controlled by `EBindingNumber`).

## Key Biology

- **Pol II** transcribes along a gene, modeled as a 1D lattice with node spacing `L_a = 100 bp`
- **R**: elongating Pol II (can bind E factors along the gene)
- **REH**: terminating Pol II that has recognized the PAS and is committed to cleavage
- **E factor**: a CPA assembly factor that binds Pol II CTD (Ser2P-phosphorylated); binding is modeled as rapid equilibrium
- **PAS (poly(A) signal)**: the point on the gene where REH complexes begin forming; located at `PASposition` bp from TSS
- **Ser2P**: CTD phosphorylation state, modeled as increasing linearly from TSS with slope `kPon_slope`
- Termination occurs when REH complexes cleave RNA at rate `kc`

## Model Architecture

### Two-timescale hybrid:

1. **Slow (ODE system)**: Pol II elongation and termination. Solved numerically with `fsolve` at steady state.
2. **Fast (symbolic/numerical equilibrium)**: E factor binding to Pol II CTD. Pre-computed and stored as a function handle `P.RE_val_bind_E(Ef)`.

### Gene discretization:

- Gene split into `N = geneLength_bp / L_a` nodes
- PAS at node `PAS = PASposition / L_a`
- `N_PAS = N - PAS + 1` nodes after PAS (where REH exists)
- Geometry and steady-state free E are stored in the parameter struct: `P.N`, `P.PAS`, `P.N_PAS`, `P.Ef_ss` — no global variables are used

### Rate matrix for E binding:

States are indexed in a 2D (P-level, E-level) block structure. The rate matrix is built **numerically** via `build_rate_matrix_numerical.m`, which is shared by both `compute_steady_states_numerical.m` and `compute_avg_E_bound_numerical.m`. The old symbolic path (`construct_rate_matrix.m` + `compute_steady_states.m`) is no longer called by any active script.

## File Map

### Core simulation

| File                                | Role                                                                                    |
| ----------------------------------- | --------------------------------------------------------------------------------------- |
| `run_termination_simulation.m`      | Top-level function: sets up geometry (`P.N/PAS/N_PAS`), pre-computes E-binding grid, solves ODE in two steps, stores `P.Ef_ss` |
| `ode_dynamics_multipleE.m`          | ODE RHS; reads geometry and rates from `P`; no globals, no nested solver               |
| `build_rate_matrix_numerical.m`     | **Shared** numerical rate matrix builder used by both steady-state functions            |
| `compute_steady_states_numerical.m` | Numerical null-space (SVD) steady-state distributions for E-binding and Ser2P           |
| `compute_avg_E_bound_numerical.m`   | Computes average E bound at each gene position numerically                              |
| `calculate_pas_cleavage_profile.m`  | Flux-based CDF of termination events downstream of PAS                                  |

### Analysis scripts

| File                        | Role                                                                |
| --------------------------- | ------------------------------------------------------------------- |
| `CpaMultipleEMain.m`        | Single comprehensive run; plots R/REH profiles and Ser2P/AvgE       |
| `ParameterSweep1D.m`        | Sweep one parameter, plot TCD (50% cutoff position)                 |
| `ParameterSweep2D.m`        | Sweep two parameters simultaneously                                 |
| `SweepParameterPasUsage.m`  | Proximal PAS usage vs inter-PAS distance (parallel, `parfor`)       |
| `EBindingNumberVsCad.m`     | TCD vs EBindingNumber sweep                                         |
| `Sweep2DkHdkEdCad.m`        | 2D sweep of kHoff/kHon × kEoff/kEon ratios, plots CAD              |
| `PlotEBindingProfile.m`     | Avg E binding + Ser2P profiles for multiple EBindingNumbers         |
| `SimulateCpaAssembly.m`     | CPA assembly CDF vs distance (analogous to Chao et al. 1999 Fig. 8) |
| `SanityCheckMultipleE.m`    | Physical edge-case assertions; E and Pol II conservation checks     |

### Gene length analysis pipeline (run in order)

| File                            | Role                                                                                 |
| ------------------------------- | ------------------------------------------------------------------------------------ |
| `GeneLengthGenerateGrid.m`      | Step 1: Create 3D grid of (R_free, E_free, L) → R_occupied, E_occupied               |
| `GeneLengthBuildInterpolation.m`| Step 2: Fit interpolation functions; define log-normal gene length PDF               |
| `GeneLengthAnalyze.m`           | Step 3–5: Solve conservation equations self-consistently; compute TCD vs gene length |

### Utilities

| File                      | Role                                                       |
| ------------------------- | ---------------------------------------------------------- |
| `default_parameters.m`    | Standard parameter set; called by all analysis scripts     |
| `save_analysis_results.m` | Standardized result/plot saving to `SecondVersionResults/` |

## Standard Parameter Set

These are defined in `default_parameters.m`. All scripts call `default_parameters()` and override only what they need.

```matlab
P.L_a        = 100;        % bp per node
P.k_in       = 2;          % Pol II initiation rate
P.k_e        = 65/100;     % Elongation rate (before PAS)
P.k_e2       = 30/100;     % Elongation rate (after PAS, in REH)
P.E_total    = 100000;     % Total E factor pool
P.Pol_total  = 70000;      % Total Pol II pool
P.kEon       = 0.000001;  % E factor on-rate
P.kEoff      = 0.2;        % E factor off-rate
P.kHon       = 1;          % PAS recognition (hexamer) on-rate
P.kHoff      = 0.5;          % Hexamer off-rate
P.kc         = 0.13;        % Cleavage rate
P.kPon_min   = 0.01;       % Min Ser2P phosphorylation rate (at TSS)
P.kPon_slope = 0.005;      % Linear slope of kPon along gene
P.kPoff      = 1;          % Ser2P dephosphorylation rate
P.geneLength_bp = 25000;   % Total gene length
P.PASposition   = 20000;   % PAS position from TSS
EBindingNumber  = 5;       % Max E factors per polymerase
```

## Two-Step Solution Strategy in `run_termination_simulation.m`

0. **Pre-computation**: `P.RE_val_bind_E` is built as a fast interpolant over a pre-computed grid of 100 `Ef_val` points (avoids repeated SVD calls inside the solver):
   ```matlab
   P.RE_val_bind_E = @(Ef_val) interpolate_E_bound(Ef_val, Ef_grid, avg_E_bound_grid, avg_Ser2P_grid);
   ```
1. **Step 1**: Solve the ODE with initial `kHon`. After convergence, compute `P.Ef_ss` once using the internal helper `solve_Efree_steady_state(R, REH, P)` which calls `fzero` on the E-conservation constraint.
2. **Step 2**: Update `P.kHon = kHon_base * avg_E_bound(P.PAS)` (scale by average E bound at PAS), then re-solve ODE from the Step 1 solution. Update `P.Ef_ss` again from the final solution.

Geometry (`P.N`, `P.PAS`, `P.N_PAS`) and `P.Ef_ss` are stored in the returned struct `P`. The returned `r_E_BeforePas` and `r_P` outputs are always `[]` — they exist only for call-site compatibility.

**There is no symbolic cache.** The `SymbolicCache/` folder and `compute_steady_states.m` are no longer used.

## Output Organization

All results go to `SecondVersionResults/<analysis_type>/` with timestamped filenames. Subdirectories include:

- `CPA_multipleE_main/` (`CpaMultipleEMain`)
- `parameter_sweep_1D/` (`ParameterSweep1D`)
- `ProximalPASUsage_ParameterSweep/` (`SweepParameterPasUsage`)
- `Sanity_check_multipleE/` (`SanityCheckMultipleE`)
- `Sweep1DEbindingnumber/` (`EBindingNumberVsCad`)
- `Ser2P_Eaverage_Profile/` (`PlotEBindingProfile`)
- `GeneLengthAnalysis/` (`GeneLengthGenerateGrid`, `GeneLengthBuildInterpolation`, `GeneLengthAnalyze` — grid data and interpolation .mat files)

## Conservation Equations (Gene Length Analysis)

The genome-wide self-consistency condition:

$$R_{\text{total}} = R_{\text{free}} + N_{\text{genes}} \int R_{\text{occupied}}(R_f, E_f, L)\, f(L)\, dL$$
$$E_{\text{total}} = E_{\text{free}} + N_{\text{genes}} \int E_{\text{occupied}}(R_f, E_f, L)\, f(L)\, dL$$

where $f(L)$ is a log-normal gene length distribution fit to human genomic data. Solved with `fsolve` in `GeneLengthAnalyze.m`.

## Common Gotchas

- **No global variables**: geometry (`P.N`, `P.PAS`, `P.N_PAS`) and free E concentration (`P.Ef_ss`) are stored in the parameter struct `P` and passed explicitly. Do not reintroduce `global` declarations.
- **`P.Ef_ss` is set by `run_termination_simulation`** after each ODE solve via `fzero`; it is available on the returned `P` struct. Scripts that need `Ef_ss` should read `P_out.Ef_ss`.
- **No symbolic cache**: the old `SymbolicCache/` folder and `compute_steady_states.m` are no longer part of the active pipeline. Do not reintroduce them.
- **Parallel workers** (`parfor`): `SweepParameterPasUsage.m` and `GeneLengthAnalyze.m` use `parfor`; each worker calls `run_termination_simulation` independently and receives its own `P` struct with `P.Ef_ss` — no shared state needed.
- **`P.RE_val_bind_E`** is a function handle `@(Ef_val) ...` returning a `1×N` vector of average E bound at each node. It is set inside `run_termination_simulation.m` and must be present in `P` before calling `ode_dynamics_multipleE`.
- **Rate matrix**: all active code builds the rate matrix numerically via `build_rate_matrix_numerical.m`. `construct_rate_matrix.m` is kept in `SparedCodes/` for reference only.
- **File naming**: main analysis scripts use PascalCase (e.g., `CpaMultipleEMain.m`); helper functions use snake_case (e.g., `run_termination_simulation.m`, `build_rate_matrix_numerical.m`).
- The `FirstVersion/` folder contains legacy code (single E binding). Do not modify; use for reference only.
- `SparedCodes/` contains archived/experimental variants. Do not call any file there from active scripts.
