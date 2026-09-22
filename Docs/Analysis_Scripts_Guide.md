# CPA Dynamics Analysis Scripts - User Guide

**Version 2.0**  
**Last Updated: September 2026**

## Table of Contents

1. [Overview](#overview)
2. [Core Analysis Scripts](#core-analysis-scripts)
3. [Parameter Sweep Scripts](#parameter-sweep-scripts)
4. [Specialized Analysis Scripts](#specialized-analysis-scripts)
5. [Output and Results](#output-and-results)
6. [Usage Guidelines](#usage-guidelines)
7. [Troubleshooting](#troubleshooting)

---

## Overview

This guide describes the complete suite of analysis scripts for the CPA (Cleavage and Polyadenylation) Dynamics model. The scripts investigate transcription termination and alternative polyadenylation (APA) through various computational approaches.

All scripts support an optional `saveData` (or `save_result`) flag at the top of the file. Set it to `true` to save results via the unified `save_analysis_results.m` utility. Results are written to `Results/<analysis_type>/`.

All active MATLAB analyses use the **full finite-rate R/RHE model** via `run_full_termination_simulation`. Binding and phosphorylation profiles are computed from microstates, and `kHon` stays at its base per-E value. The working engaged-E off-rate is `P.kEoff_engaged = 0.05` s⁻¹ in `default_parameters.m`.

No Symbolic Math or Optimization Toolbox is needed for the active steady solvers. Scripts using `parfor`/`parpool` require Parallel Computing Toolbox.

---

## Core Analysis Scripts

### 1. CpaMultipleEMain.m

**Purpose**: Single comprehensive CPA simulation run

**What it does**:

- Solves the full ODE system for polymerase dynamics
- Computes steady-state concentrations of R (elongating) and REH (terminating) polymerases
- Calculates Ser2P and average E binding from full microstate populations, including RHE downstream of PAS
- Generates two plots: Ser2P/AverageE profiles, and R/REH concentration profiles
- Reports polymerase distribution at steady state

**Key Parameters**:

- `EBindingNumber`: Max E factors per polymerase (default: 5)
- `geneLength_bp`: Total gene length (default: 25,000 bp)
- `PASposition`: Position of poly(A) signal (default: 20,000 bp)
- `saveData`: Set `true` to save results (default: `false`)

**Expected Results**:

- Plot 1: Ser2P and AverageE vs position relative to PAS
- Plot 2: R(l) and REH(l) concentrations with 50% TCD marker
- Console output: free/bound Pol II distribution and free E concentration

**When to use**:

- Initial model validation and baseline results
- Visualizing polymerase and E factor distributions

---

## Parameter Sweep Scripts

### 2. ParameterSweep1D.m

**Purpose**: Systematic exploration of how individual parameters affect the Cleavage Arrest Distance (CAD)

**What it does**:

- Iterates over a list of parameters, sweeping each across a biologically relevant range
- Calculates the position where 50% of polymerases have terminated (CAD₅₀)
- Generates a line plot per parameter showing CAD vs parameter value
- Marks the default parameter value on each plot

**Key Features**:

- Covers 11 parameters: `k_e`, `k_e2`, `E_total`, `Pol_total`, `kc`, `kEon`, `kEoff`, `k_in`, `kHoff`, `kHon`, `kPon_slope`
- Automatic linear/log range selection per parameter type
- `save_result` flag (default: `false`) to enable saving via `save_analysis_results`

**Expected Results**:

- One plot per parameter showing CAD₅₀ (bp) vs parameter value
- Default parameter value highlighted

**When to use**:

- Parameter sensitivity analysis
- Identifying which parameters most strongly control termination distance

### 3. ParameterSweep2D.m

**Purpose**: Explores interactions between pairs of parameters

**What it does**:

- Performs a 2D grid sweep over a configurable pair of parameters
- Calculates CAD₅₀ at each (param1, param2) combination
- Generates multi-line plots showing how one parameter modulates the other's effect

**Key Features**:

- Configurable parameter pairs via `param_pairs`
- `save_result` flag (default: `false`)
- Robust NaN handling for failed simulations

**Expected Results**:

- Multi-line plot: CAD₅₀ vs param1, with separate lines for each param2 value

**When to use**:

- Understanding parameter interactions
- Identifying combinations that produce specific termination behaviors

### 4. Sweep2DCad.m

**Purpose**: Full finite-rate CAD₅₀ contour maps for all three pairs of PAS recognition affinity, E factor affinity and cleavage rate.

**What it does**:

- Sweeps (x, y) = (kHd, kEd), (kc, kHd) and (kc, kEd), with `kHd = kHoff/kHon` and `kEd = kEoff/kEon`.
- Uses log-spaced ranges from `P.ranges` in `default_parameters.m`, including `P.ranges.kc`.
- Fixes on-rates at reference values and reconstructs off-rates from dissociation constants. The unswept parameter and `kEoff_engaged` stay at baseline.
- Uses `run_full_termination_simulation`; CAD thresholds outside the simulated window remain `NaN`.
- Computes the baseline once per binding capacity and marks it on each map.

**Key Features**:

- Set `nPoints`, `percent_cleavage` and `EBindingNumbers` in the script. Default binding capacity is 5; use `[1 5]` for comparison.
- `save_result` flag (default: `true`)
- Version-safe colormap (falls back from `turbo` to `jet`)
- All sweep and plotting helpers are local functions in `Sweep2DCad.m`; run that script to compute, plot and save all three pairs.
- Raw CAD, within-window cleavage fractions, ODE residuals and per-cell errors are saved with explicit matrix orientation (rows = y, columns = x).

**Expected Results**:

- Three figures per binding capacity, with logarithmic axes and a shared 0–1600 bp display scale. Saved CAD values are uncapped.
- Red dashed contours at 400 bp and 800 bp when crossed by the data
- Star marker at the base-parameter point, with a text-only `Base parameters CAD_50 = ... bp` legend in the upper-left corner
- Open circles for unreached thresholds and crosses for failed simulations; missing CAD regions stay gray.
- Separate output folders: `sweep_2D_kHd_kEd_CAD`, `sweep_2D_kc_kHd_CAD` and `sweep_2D_kc_kEd_CAD` under `Results/` (or `CPAD_RESULTS_ROOT` when configured).

**When to use**:

- Parameter robustness analysis for kHd, kEd and kc
- Demonstrating that the model is consistent with experimental CAD ranges

---

## Specialized Analysis Scripts

### 5. SweepParameterPasUsage.m

**Purpose**: Predicts proximal PAS usage as a function of inter-PAS distance, swept over a chosen global parameter

**What it does**:

- Runs parallel simulations (`parfor`) over a grid of inter-PAS distances
- For each simulation, computes the termination CDF; the CDF value at a given distance equals proximal site usage
- Sweeps a configurable parameter (`kHoff`, `E_total`, `k_e`, or `kc`) with preset value sets
- Plots proximal usage (%) vs inter-PAS distance on a log x-axis

**Key Features**:

- Parallel processing via `parfor`
- Switch-case parameter selection at the top of the script
- `save_results` flag (default: `true`)

**Expected Results**:

- Multi-line semilog plot: proximal usage (%) vs inter-PAS distance
- One line per sweep parameter value

**When to use**:

- Predicting APA site choice under different conditions
- Understanding how global factors modulate distance-dependent competition

### 6. EBindingNumberVsCad.m

**Purpose**: Investigates how the maximum number of E factors per polymerase affects the CAD

**What it does**:

- Sweeps `EBindingNumber` from 1 to 6
- Computes CAD₅₀ using the full finite-rate model for engaged-E off-rates of 0, 0.05 and 0.5 s⁻¹
- Plots one CAD (bp) vs EBindingNumber curve per engaged-E off-rate; zero is the protected-E limit

**Key Features**:

- `saveData` flag (default: `false`)
- Robust error handling; failed runs recorded as NaN

**Expected Results**:

- Scatter/line plot showing how CAD₅₀ changes with maximum E factor occupancy

**When to use**:

- Calibrating `EBindingNumber`
- Assessing sensitivity of CAD to E factor cooperativity

### 7. PlotEBindingProfile.m

**Purpose**: Generates combined average E binding and Ser2P spatial profiles for multiple EBindingNumbers

**What it does**:

- Runs simulations for `EBindingNumber` = 1, 3, 5
- Plots average E bound (solid) and Ser2P (dashed) profiles on a single figure, color-coded by EBindingNumber
- Optionally saves per-profile data as tab-delimited `.txt` files and the figure as `.png`

**Key Features**:

- `saveData` flag (default: `false`); when `true`, saves to `Results/Ser2P_Eaverage_Profile/` and `Results/SupportFigures/`
- x-axis in kb relative to PAS

**Expected Results**:

- Combined profile plot comparing E binding and phosphorylation across binding numbers

**When to use**:

- Visualizing how E factor occupancy and CTD phosphorylation vary along the gene
- Generating support figures

### 8. SimulateCpaAssembly.m

**Purpose**: Simulates CPA complex assembly kinetics as a function of distance downstream of the PAS

**What it does**:

- Runs the full finite-rate simulation and computes the cleavage-flux CDF
- Interpolates the CDF from the physical origin to a fine distance grid (0–1000 bp, 10 bp steps); queries beyond the simulated window return NaN
- Plots "% CPA Assembly Completion" vs poly(A)–anti-poly(A) separation distance
- Analogous to Figure 8 from Chao et al. (1999)

**Key Features**:

- `saveData` flag (default: `false`)
- Prints key values at 100, 200, …, 1000 bp

**Expected Results**:

- Sigmoidal CDF curve showing assembly completion vs distance
- Quantitative rescue percentages at biologically relevant separations

**When to use**:

- Comparing model predictions to cis-antisense rescue assay data

### 9. SanityCheckMultipleE.m

**Purpose**: Verifies model correctness through conservation checks and edge case tests

**What it does**:

- Checks Pol II and E factor conservation at steady state
- Verifies non-negativity of all state variables
- Tests no E supply, no E binding, no recognition and no cleavage; verifies that zero initiation/elongation is rejected by the algebraic steady solver

**Expected Results**:

- Assertions stop on failures; successful checks print PASS
- Baseline R/RHE spatial plot

**When to use**:

- After any model change, to confirm physical consistency
- Debugging unexpected simulation behavior

---

## Gene-length analysis pipeline

Run these in order after migrating or changing the model parameters:

```matlab
GeneLengthGenerateGrid
GeneLengthBuildInterpolation
GeneLengthAnalyze
```

The grid solves full finite-rate kinetics at prescribed `[R_free, E_free]` using `'FreePools'`. It records Pol II and E demand directly from microstates and covers zero through each global pool total. The interpolator requires a complete Cartesian grid, preserves linear scaling in free Pol II, and disables extrapolation.

`GeneLengthAnalyze` uses the grid's parameters and downstream window. It normalizes the length distribution over the simulated interval, eliminates free Pol II analytically, and solves E conservation with bounded `fzero` through its local `solve_full_genome_pools` function. It then computes full-model CAD profiles at the shared pools. Unreached thresholds remain NaN, with cleavage fractions, solver errors and residuals saved alongside them.

Artifacts use `full_gene_length_*` filenames in `SecondVersionResults/GeneLengthAnalysis/` (or `CPAD_RESULTS_ROOT`). Old equilibrium grid/interpolation files are not reused. Grid-node reconstruction errors do not establish interpolation accuracy between nodes; increase grid resolution to assess convergence.

## Verification

Run `runtests('tests/test_full_termination_model.m')`, `runtests('tests/test_full_model_analysis.m')` and `SanityCheckMultipleE` in MATLAB. The migration tests cover fixed-pool equivalence, resource accounting, empty states, analytical genome conservation and recovery of a single-gene closed-pool solution.

## Output and Results

### Saving Results

All scripts have a `saveData` (or `save_result` / `save_results`) flag at the top. Set it to `true` to save output via `save_analysis_results.m`.

**File Organization**:

- Results saved to `Results/<analysis_type>/`
- Subfolders created automatically

**File Types Generated**:

1. **Plot Files** (`.png`): High-quality visualizations
2. **Data Files** (`_data.txt`): Numerical results with parameter metadata

**File Naming Convention** (no timestamp):

```
CPA_main_EBinding5.png
CPA_main_EBinding5_data.txt
Sweep1D_kHoff_EBinding5.png
sweep2D_kHd_kEd_CAD_EBinding1_data.txt
ProxPASUsage_kHoff_data.txt
EBinding_vs_CAD_N1-6_data.txt
CPA_assembly_EBinding5_data.txt
```

---

## Usage Guidelines

### Prerequisites

- MATLAB (no Symbolic Math Toolbox required)
- Parallel Computing Toolbox (optional, for `SweepParameterPasUsage.m` and `GeneLengthAnalyze.m`)

### Recommended Workflow

1. **Start with `CpaMultipleEMain.m`** to understand basic model behavior
2. **Use `ParameterSweep1D.m`** to identify sensitive parameters
3. **Use `Sweep2DCad.m`** for pairwise robustness analysis over kHd, kEd and kc
4. **Use `SweepParameterPasUsage.m`** for APA site choice predictions
5. **Run `SanityCheckMultipleE.m`** after any model modifications

### Saving Output

Set the flag at the top of each script:

```matlab
saveData = true;   % or save_result / save_results depending on the script
```

### Computational Considerations

- **Parallel scripts** (`SweepParameterPasUsage.m`, `GeneLengthAnalyze.m`): Use `parfor`; a parallel pool starts automatically
- **Runtime**: Single runs ~seconds; large 2D sweeps ~minutes to hours
- **Memory**: 2D sweeps can require several GB of RAM

---

## Troubleshooting

### Common Issues

**1. "Solver failed" errors**:

- Check for unphysical parameter combinations (e.g., rates near zero)
- Verify the active full-model helpers (`run_full_termination_simulation`, `calculate_full_pas_cleavage_profile`, etc.) are on the MATLAB path

**2. Negative or very large concentrations**:

- Indicates parameter combinations outside physically reasonable range
- Start from default parameters and modify gradually

**3. Slow performance**:

- Reduce the number of sweep points
- Close unused parallel pools: `delete(gcp('nocreate'))`

**4. Missing plots or save errors**:

- Ensure write permissions for the `Results/` directory
- Check that `save_analysis_results.m` is on the MATLAB path

### Parameter Guidelines

**Biologically Reasonable Ranges**:

- `kHoff`: 0.1 – 10 (PAS recognition off-rate)
- `E_total`: 10,000 – 200,000 (CPA factor pool size)
- `kc`: 0.01 – 1.0 (cleavage rate)
- `EBindingNumber`: 1 – 8 (max E factors per Pol II)

**Numerical Stability**:

- Avoid extreme parameter ratios (>1000× differences)
- Use logarithmic spacing for rate constants in sweeps
- Start with default parameters (`default_parameters.m`) and modify one at a time

---

_This guide covers the core functionality of the CPA Dynamics analysis suite. For specific biological interpretations or advanced customizations, consult the accompanying research documentation._
