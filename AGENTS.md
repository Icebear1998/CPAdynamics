# CPAdynamics

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
- **E factor**: a CPA assembly factor that binds Pol II CTD (Ser2P-phosphorylated); binding is modeled with explicit finite rates
- **PAS (poly(A) signal)**: the point on the gene where REH complexes begin forming; located at `PASposition` bp from TSS
- **Ser2P**: CTD phosphorylation state, modeled as increasing linearly from TSS with slope `kPon_slope`
- Termination occurs when REH complexes cleave RNA at rate `kc`

## Model Architecture

### Full finite-rate model

All active MATLAB analyses use `run_full_termination_simulation`. Phosphorylation,
E binding, PAS recognition, elongation and cleavage are explicit finite-rate
reactions. There is no rapid-equilibrium closure and `P.kHon` is never rescaled.
R states have `0 <= e <= p <= M`; recognized RHE states have `1 <= e <= p <= M`.
Exactly one E in each RHE state is engaged with H. `kEoff_engaged` controls its
independent detachment and rapid EH disassembly; the working default is 0.05 s^-1.

Geometry is `P.N = floor(geneLength_bp/L_a)`, `P.PAS = floor(PASposition/L_a)`
and `P.N_PAS = P.N-P.PAS+1`. All geometry and free pools are stored in `P`.

All active analyses and tests use the full finite-rate helpers in the project root.

## File Map

### Core simulation

| File | Role |
| --- | --- |
| `run_full_termination_simulation.m` | Full finite-rate steady state; closed or fixed free pools; microstate-derived averages and resource demand |
| `build_full_rate_matrices.m` | Sparse transport/reaction operator with local state-map and internal-generator helpers |
| `ode_dynamics_full_multipleE.m` | Full ODE RHS with explicit E and Pol II pools |
| `full_model_jacobian.m` | Sparse full ODE Jacobian |
| `calculate_full_pas_cleavage_profile.m` | Flux-based cleavage CDF and within-window CAD |

### Analysis scripts

| File                       | Role                                                                |
| -------------------------- | ------------------------------------------------------------------- |
| `CpaMultipleEMain.m`       | Single comprehensive run; plots R/REH profiles and Ser2P/AvgE       |
| `ParameterSweep1D.m`       | Sweep one parameter, plot TCD (50% cutoff position)                 |
| `ParameterSweep2D.m`       | Sweep two parameters simultaneously                                 |
| `SweepParameterPasUsage.m` | Proximal PAS usage vs inter-PAS distance (parallel, `parfor`)       |
| `EBindingNumberVsCad.m`    | TCD vs EBindingNumber sweep                                         |
| `Sweep2DCad.m`       | Pairwise CAD maps: kHd-kEd, kHd-kc and kEd-kc               |
| `PlotEBindingProfile.m`    | Avg E binding + Ser2P profiles for multiple EBindingNumbers         |
| `SimulateCpaAssembly.m`    | CPA assembly CDF vs distance (analogous to Chao et al. 1999 Fig. 8) |
| `SanityCheckMultipleE.m`   | Physical edge-case assertions; E and Pol II conservation checks     |

### Gene length analysis pipeline (run in order)

| File                             | Role                                                                                 |
| -------------------------------- | ------------------------------------------------------------------------------------ |
| `GeneLengthGenerateGrid.m`       | Step 1: Create 3D grid of (R_free, E_free, L) → R_occupied, E_occupied               |
| `GeneLengthBuildInterpolation.m` | Step 2: Fit interpolation functions; define log-normal gene length PDF               |
| `GeneLengthAnalyze.m`            | Step 3–5: Solve conservation equations self-consistently; compute TCD vs gene length |

### Utilities

| File                      | Role                                                       |
| ------------------------- | ---------------------------------------------------------- |
| `default_parameters.m`    | Standard parameter set; called by all analysis scripts     |
| `save_analysis_results.m` | Standardized result/plot saving to `Results/` |

## Standard Parameter Set

These are defined in `default_parameters.m`. All scripts call `default_parameters()` and override only what they need.

```matlab
P.L_a        = 100;        % bp per node
P.k_in       = 2;          % Pol II initiation rate
P.k_e        = 65/100;     % Elongation rate (before PAS)
P.k_e2       = 30/100;     % Elongation rate (after PAS, in REH)
P.E_total    = 100000;     % Total E factor pool
P.Pol_total  = 70000;      % Total Pol II pool
P.kEon       = 2.5e-6;  % E factor on-rate
P.kEoff      = 0.5;        % E factor off-rate
P.kEoff_engaged = 0.05;    % Independent engaged-E off-rate
P.kHon       = 7;          % PAS recognition (hexamer) on-rate
P.kHoff      = 1;          % Hexamer off-rate
P.kc         = 0.15;        % Cleavage rate
P.kPon_min   = 0.01;       % Min Ser2P phosphorylation rate (at TSS)
P.kPon_slope = 0.005;      % Linear slope of kPon along gene
P.kPoff      = 1;          % Ser2P dephosphorylation rate
P.geneLength_bp = 25000;   % Total gene length
P.PASposition   = 20000;   % PAS position from TSS
EBindingNumber  = 5;       % Max E factors per polymerase
```

## Full-model solution strategy

`run_full_termination_simulation(P, M)` builds one sparse operator, solves its
linear microstate equations at each trial free E, scales occupancy to conserve
Pol II, and uses `fzero` to conserve E. Returned `P.Ef_ss` and `P.Pol_free_ss`
are consistent with the returned R and RHE profiles.

`run_full_termination_simulation(P, M, 'FreePools', [R_free E_free])` holds both
free pools fixed. It is used for the gene-length lookup grid and profiles after
the shared genome-wide pools have been solved. Local totals must not be enforced
again in this mode. Read `details.Pol_bound` and `details.E_bound` for resource
demand. The full kinetic RHS and flux balance are still validated.

Read `details.avg_E_bound` and `details.avg_Ser2P` for averages over all polymerases
at each node, including RHE after PAS. Empty nodes have NaN averages and zero
resource demand. No binding equilibrium interpolant or symbolic cache is used.

The gene-length pipeline writes `full_gene_length_grid_data_*.mat`, then
`full_gene_length_interpolation_*.mat`, then `full_gene_length_TCD_analysis_*.mat`.
Rebuild both prerequisites after a model or parameter change. Old equilibrium
artifacts are rejected. The grid includes zero through the global pool totals;
interpolation never extrapolates and the length PDF is normalized over the grid's
finite length interval. The saved metadata records that interval and PDF mass.

## Output Organization

The standard saver writes to `Results/<analysis_type>/`; full binding-capacity and gene-length analyses use `SecondVersionResults/`. `CPAD_RESULTS_ROOT` overrides either root. The latter analyses use timestamped outputs. Subdirectories include:

- `CPA_multipleE_main/` (`CpaMultipleEMain`)
- `parameter_sweep_1D/` (`ParameterSweep1D`)
- `ProximalPASUsage_ParameterSweep/` (`SweepParameterPasUsage`)
- `Full_EBindingNumber_vs_CAD/` (`EBindingNumberVsCad`)
- `Ser2P_Eaverage_Profile/` (`PlotEBindingProfile`)
- `GeneLengthAnalysis/` (`GeneLengthGenerateGrid`, `GeneLengthBuildInterpolation`, `GeneLengthAnalyze` — grid data and interpolation .mat files)

## Conservation Equations (Gene Length Analysis)

The genome-wide self-consistency condition:

$$R_{\text{total}} = R_{\text{free}} + N_{\text{genes}} \int R_{\text{occupied}}(R_f, E_f, L)\, f(L)\, dL$$
$$E_{\text{total}} = E_{\text{free}} + N_{\text{genes}} \int E_{\text{occupied}}(R_f, E_f, L)\, f(L)\, dL$$

where $f(L)$ is a log-normal gene length distribution fit to human genomic data. Solved by the local `solve_full_genome_pools` function in `GeneLengthAnalyze.m`: analytically eliminate R_free using linear occupancy scaling and use bracketed `fzero` for E_free.

## Common Gotchas

- Never introduce globals, symbolic caches or equilibrium binding handles into active code.
- `P.kHon` remains the base per-E recognition rate in every full-model simulation.
- Use `calculate_full_pas_cleavage_profile`: unreached CAD thresholds return NaN. Report `diagnostics.max_exit_cdf` alongside them.
- `run_full_termination_simulation` requires positive initiation and elongation rates; do not turn solver errors into successful zero-population cases.
- Shared-pool gene-length runs must pass both pools explicitly through `'FreePools'` and reuse the grid's parameters and downstream length.
- `parfor` workers each receive their own parameter struct; no shared mutable solver state is needed.
- Average E and Ser2P profiles come from full microstates, not `P.RE_val_bind_E`.
- Main analysis scripts use PascalCase; helpers use snake_case.
- `FirstVersion/` and `SparedCodes/` contain legacy/reference code. Do not modify or call them from active scripts.
