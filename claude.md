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
- Global variables: `N`, `PAS`, `N_PAS`, `Ef_ss`

### Rate matrix for E binding:

States are indexed in a 2D (P-level, E-level) block structure. The rate matrix is built **numerically** via `build_rate_matrix_numerical.m`, which is shared by both `compute_steady_states_numerical.m` and `compute_avg_E_bound_numerical.m`. The old symbolic path (`construct_rate_matrix.m` + `compute_steady_states.m`) is no longer called by any active script.

## File Map

### Core simulation

| File                                | Role                                                                                    |
| ----------------------------------- | --------------------------------------------------------------------------------------- |
| `run_termination_simulation.m`      | Top-level function: sets up geometry, builds `P.RE_val_bind_E`, solves ODE in two steps |
| `ode_dynamics_multipleE.m`          | ODE RHS; computes self-consistent `Ef_ss` via `fsolve` on first call                    |
| `construct_rate_matrix.m`           | Symbolic rate matrix — kept for reference; not called by active scripts                 |
| `build_rate_matrix_numerical.m`     | **Shared** numerical rate matrix builder used by both steady-state functions            |
| `compute_steady_states_numerical.m` | Numerical null-space (SVD) steady-state distributions for E-binding and Ser2P           |
| `compute_avg_E_bound_numerical.m`   | Computes average E bound at each gene position numerically                              |
| `calculate_pas_usage_profile.m`     | Flux-based CDF of termination events downstream of PAS                                  |

### Analysis scripts

| File                       | Role                                                                |
| -------------------------- | ------------------------------------------------------------------- |
| `CPA_multipleE_main.m`     | Single comprehensive run; plots R/REH profiles and Ser2P/AvgE       |
| `parameter_sweep_1D.m`     | Sweep one parameter, plot TCD (50% cutoff position)                 |
| `parameter_sweep_2D.m`     | Sweep two parameters simultaneously                                 |
| `APA_sweep_2PAS.m`         | Two-PAS model: proximal vs distal site usage vs distance/strength   |
| `SweepParameterPASusage.m` | Proximal PAS usage vs inter-PAS distance (parallel, `parfor`)       |
| `EBindingNumberVsCAD.m`    | TCD vs EBindingNumber sweep                                         |
| `sweep_2D_kEoff_kPon.m`    | 2D sweep of kEoff × kPon_slope, plots avg E at PAS                  |
| `PlotEbindingProfile.m`    | Avg E binding + Ser2P profiles for multiple EBindingNumbers         |
| `Simulate_CPAAssembly.m`   | CPA assembly CDF vs distance (analogous to Chao et al. 1999 Fig. 8) |

### Gene length analysis pipeline (run in order)

| File                                | Role                                                                                 |
| ----------------------------------- | ------------------------------------------------------------------------------------ |
| `generate_gene_length_grid.m`       | Step 1: Create 3D grid of (R_free, E_free, L) → R_occupied, E_occupied               |
| `build_gene_length_interpolation.m` | Step 2: Fit interpolation functions; define log-normal gene length PDF               |
| `analyze_gene_length_TCD.m`         | Step 3–5: Solve conservation equations self-consistently; compute TCD vs gene length |

### Utilities

| File                      | Role                                                       |
| ------------------------- | ---------------------------------------------------------- |
| `save_analysis_results.m` | Standardized result/plot saving to `SecondVersionResults/` |

## Standard Parameter Set

```matlab
P.L_a        = 100;        % bp per node
P.k_in       = 2;          % Pol II initiation rate
P.k_e        = 65/100;     % Elongation rate (before PAS)
P.k_e2       = 30/100;     % Elongation rate (after PAS, in REH)
P.E_total    = 100000;     % Total E factor pool
P.Pol_total  = 70000;      % Total Pol II pool
P.kEon       = 0.0000025;  % E factor on-rate
P.kEoff      = 0.1;        % E factor off-rate
P.kHon       = 2;          % PAS recognition (hexamer) on-rate
P.kHoff      = 1;          % Hexamer off-rate
P.kc         = 0.1;        % Cleavage rate
P.kPon_min   = 0.01;       % Min Ser2P phosphorylation rate (at TSS)
P.kPon_slope = 0.005;      % Linear slope of kPon along gene
P.kPoff      = 1;          % Ser2P dephosphorylation rate
P.geneLength_bp = 25000;   % Total gene length
P.PASposition   = 20000;   % PAS position from TSS
EBindingNumber  = 5;       % Max E factors per polymerase
```

## Two-Step Solution Strategy in `run_termination_simulation.m`

1. **Step 1** (`P.FirstRun = true`): Solve ODE with initial `kHon`, simultaneously solving for self-consistent `Ef_ss` (free E concentration) using a nested `fsolve`.
2. **Step 2** (`P.FirstRun = false`): Update `P.kHon = kHon_base * avg_E_bound(PAS)` (scale by average E bound at PAS), then re-solve ODE from Step 1 solution.

`P.RE_val_bind_E` is set inside `run_termination_simulation.m` as:

```matlab
P.RE_val_bind_E = @(Ef_val) compute_avg_E_bound_numerical(Ef_val, kPon_vals, P.kPoff, P.kEon, P.kEoff, n_states);
```

The returned `r_E_BeforePas` and `r_P` outputs are always `[]` — they exist only for call-site compatibility.

**There is no symbolic cache.** The `SymbolicCache/` folder and `compute_steady_states.m` are no longer used.

## Output Organization

All results go to `SecondVersionResults/<analysis_type>/` with timestamped filenames. Subdirectories include:

- `APA_sweep_2PAS/`
- `CPA_multipleE_main/`
- `parameter_sweep_1D/`
- `PASUsageAnalysis/`
- `PASUsagevsInterPASDistance/`
- `Sanity_check_multipleE/`
- `Sweep1DEbindingnumber/`
- `GeneLengthAnalysis/` (grid data and interpolation .mat files)

## Conservation Equations (Gene Length Analysis)

The genome-wide self-consistency condition:

$$R_{\text{total}} = R_{\text{free}} + N_{\text{genes}} \int R_{\text{occupied}}(R_f, E_f, L)\, f(L)\, dL$$
$$E_{\text{total}} = E_{\text{free}} + N_{\text{genes}} \int E_{\text{occupied}}(R_f, E_f, L)\, f(L)\, dL$$

where $f(L)$ is a log-normal gene length distribution fit to human genomic data. Solved with `fsolve` in `analyze_gene_length_TCD.m`.

## Common Gotchas

- **Global variables** (`N`, `PAS`, `N_PAS`, `Ef_ss`) must be declared `global` in any function that uses them and in the calling script.
- **No symbolic cache**: the old `SymbolicCache/` folder and `compute_steady_states.m` are no longer part of the active pipeline. Do not reintroduce them.
- **Parallel workers** (`parfor`): `SweepParameterPASusage.m` and `analyze_gene_length_TCD.m` use `parfor`; global variables are not shared across workers, so each worker must recompute `Ef_ss` locally.
- **`P.RE_val_bind_E`** is a function handle `@(Ef_val) ...` returning a `1×N` vector of average E bound at each node. It is set inside `run_termination_simulation.m` and must be present in `P` before calling `ode_dynamics_multipleE`.
- **Rate matrix**: all active code builds the rate matrix numerically via `build_rate_matrix_numerical.m`. `construct_rate_matrix.m` is kept for reference only.
- The `FirstVersion/` folder contains legacy code (single E binding). Do not modify; use for reference only.
- `SparedCodes/` contains archived/experimental variants. `compute_steady_states.m` lives there but is no longer on the active call path.
