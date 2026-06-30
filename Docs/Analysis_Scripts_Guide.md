# CPA Dynamics Analysis Scripts - User Guide

**Version 2.0**  
**Last Updated: June 2026**

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

The model uses a fully **numerical** approach — no Symbolic Math Toolbox is required.

---

## Core Analysis Scripts

### 1. CpaMultipleEMain.m

**Purpose**: Single comprehensive CPA simulation run

**What it does**:

- Solves the full ODE system for polymerase dynamics
- Computes steady-state concentrations of R (elongating) and REH (terminating) polymerases
- Calculates Ser2P phosphorylation levels and average E factor binding profiles
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

### 4. Sweep2DkHdkEdCad.m

**Purpose**: 2D contour map of CAD₅₀ as a function of PAS recognition and E factor binding dissociation constants

**What it does**:

- Sweeps `kHd = kHoff/kHon` and `kEd = kEoff/kEon` over literature-plausible ranges
- Fixes on-rates at reference values; varies off-rates to achieve target dissociation constants
- Generates a filled contour plot with experimental target band (400–800 bp) overlaid
- Computes and marks the CAD at base parameters

**Key Features**:

- Configurable literature ranges for individual rates
- `save_result` flag (default: `true`)
- Version-safe colormap (falls back from `turbo` to `jet`)

**Expected Results**:

- Contour map of CAD₅₀ across (kHd, kEd) space
- Red dashed contours at 400 bp and 800 bp (experimental window)
- Star marker at base-parameter point

**When to use**:

- Parameter robustness analysis for kHon/kHoff and kEon/kEoff
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
- Computes CAD₅₀ for each value
- Plots CAD (bp) vs EBindingNumber with data labels

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

- Runs a single simulation and computes the termination CDF
- Interpolates the CDF to a fine grid of distances (0–1000 bp, 10 bp steps)
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
- Tests 6 edge cases (k_in = 0, k_e = 0, k_e2 = 0, kc = 0, E_total = 0, large kHon/kc)

**Expected Results**:

- Console output with PASS/FAIL for each check
- Diagnostic plots generated only for FAIL cases

**When to use**:

- After any model change, to confirm physical consistency
- Debugging unexpected simulation behavior

---

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
3. **Use `Sweep2DkHdkEdCad.m`** for robustness analysis over kHon/kEon space
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
- Verify all helper functions (`run_termination_simulation`, `calculate_pas_cleavage_profile`, etc.) are on the MATLAB path

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
