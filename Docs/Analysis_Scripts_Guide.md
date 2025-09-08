# CPA Dynamics Analysis Scripts - User Guide

**Version 2.0**  
**Last Updated: December 2024**

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

This guide describes the complete suite of analysis scripts for the CPA (Cleavage and Polyadenylation) Dynamics model. The scripts are designed to investigate different aspects of transcription termination and alternative polyadenylation (APA) through various computational approaches.

All scripts automatically save their results to organized folders within `SecondVersionResults/` using the unified `save_analysis_results.m` utility function.

---

## Core Analysis Scripts

### 1. CPA_multipleE_main.m

**Purpose**: Main simulation script that runs a single comprehensive CPA analysis

**What it does**:

- Solves the full ODE system for polymerase dynamics
- Computes steady-state concentrations of R (elongating) and REH (paused) polymerases
- Calculates Ser2P phosphorylation levels and average E factor binding
- Generates visualization plots of polymerase distributions

**Key Parameters**:

- `EBindingNumber`: Number of E factors that can bind (default: 3)
- `L_a`: Node size in base pairs (default: 100)
- `geneLength_bp`: Total gene length (default: 25,000 bp)
- `PASposition`: Position of polyadenylation site (default: 20,000 bp)

**Expected Results**:

- Two plots: Ser2P and AverageE vs position, and R/REH concentrations vs position
- Console output showing total polymerase distribution
- Saved data includes all concentration profiles and key metrics

**When to use**:

- Initial model validation
- Understanding basic polymerase dynamics
- Generating baseline results for comparison

---

## Parameter Sweep Scripts

### 2. parameter_sweep_1D.m

**Purpose**: Systematic exploration of how individual parameters affect termination cutoff positions

**What it does**:

- Sweeps through ranges of single parameters (kHoff, E_total, kc, etc.)
- Calculates the position where 25% of polymerases have terminated
- Generates plots showing parameter sensitivity
- Supports both linear and logarithmic parameter ranges

**Key Features**:

- Configurable parameter selection via `param_list`
- Automatic range determination based on parameter type
- Robust error handling for failed simulations
- High-resolution plot generation with timestamps

**Parameters Available**:

- `k_e`, `k_e2`: Elongation rates
- `E_total`: Total E factor concentration
- `kc`: Cleavage rate
- `kHon`, `kHoff`: Hexamer binding/unbinding rates
- `kEon`, `kEoff`: E factor binding rates

**Expected Results**:

- Line plots showing cutoff position vs parameter value
- Clear identification of parameter sensitivity regions
- Default parameter values highlighted on plots

**When to use**:

- Parameter sensitivity analysis
- Identifying critical parameter ranges
- Model calibration and validation

### 3. parameter_sweep_2D.m

**Purpose**: Explores interactions between pairs of parameters

**What it does**:

- Performs 2D parameter sweeps over parameter pairs
- Uses intelligent initial guess extrapolation for numerical stability
- Generates multi-line plots showing parameter interactions
- Optimized for computational efficiency

**Key Features**:

- Extrapolation-based initial guessing for faster convergence
- Configurable parameter pairs via `param_pairs`
- Automatic history tracking for solution continuation
- Robust handling of numerical failures

**Expected Results**:

- Multi-line plots showing cutoff positions across 2D parameter space
- Clear visualization of parameter interaction effects
- Legend showing different parameter values

**When to use**:

- Understanding parameter interactions
- Identifying parameter combinations that produce specific behaviors
- Advanced model exploration

---

## Specialized Analysis Scripts

### 4. PASUsageAnalysis.m

**Purpose**: Comprehensive 2D analysis of PAS usage as a function of global and local factors

**What it does**:

- Performs 2D sweep over kHoff (PAS strength) and E_total (global factor)
- Calculates cumulative polymerase exit at a fixed distance (300 bp)
- Uses parallel processing for computational efficiency
- Generates contour plots showing phase diagrams

**Key Features**:

- Parallel processing with `parfor` loops
- Flux-based calculation of termination probabilities
- Interpolation to specific distances for biological relevance
- High-resolution contour visualization

**Expected Results**:

- Contour plot showing cumulative exit percentage
- Clear phase boundaries between different termination regimes
- Quantitative metrics at biologically relevant distances

**When to use**:

- Understanding global vs local factor balance
- Predicting APA site usage under different conditions
- Generating phase diagrams for publication

### 5. PASUsagevsInterPASDistance.m

**Purpose**: Predicts how proximal PAS usage depends on the distance to downstream PAS sites

**What it does**:

- Simulates termination profile for a single parameter set
- Calculates proximal site usage probability vs inter-PAS distance
- Includes biological context (300 bp median distance marker)
- Provides direct predictions for experimental comparison

**Key Features**:

- Biologically motivated distance ranges (0-2500 bp)
- Direct calculation of site choice probabilities
- Visualization with experimental context markers
- Robust solver with error handling

**Expected Results**:

- Plot showing proximal usage % vs inter-PAS distance
- Exponential-like decay typical of competition models
- Reference line at median biological distance (300 bp)

**When to use**:

- Predicting experimental outcomes
- Understanding distance-dependent competition
- Validating model against APA-seq data

### 6. Sweep1DEbindingnumber.m

**Purpose**: Investigates how the maximum number of E factors affects average binding

**What it does**:

- Sweeps through different EBindingNumber values (1-3)
- Calculates average E factor binding at the PAS
- Generates plots showing saturation behavior
- Focuses on E factor binding dynamics

**Key Features**:

- Fixed parameter set for controlled comparison
- Focus on E factor binding saturation
- Clear visualization of binding vs capacity
- Robust error handling for edge cases

**Expected Results**:

- Plot showing average E binding vs EBindingNumber
- Saturation behavior at high EBindingNumber values
- Quantitative binding coefficients

**When to use**:

- Understanding E factor binding saturation
- Calibrating EBindingNumber parameter
- Investigating cooperative binding effects

### 7. SweepParameterPASusage.m

**Purpose**: Flexible parameter sweep tool for analyzing PAS usage across distances

**What it does**:

- Configurable parameter sweeps (E_total, k_e, kc, kHoff)
- Calculates cumulative exit CDF across distance ranges
- Uses parallel processing for efficiency
- Generates multi-line plots showing parameter effects

**Key Features**:

- Switch-case parameter selection system
- Parallel processing with `parfor`
- Flux-based CDF calculations
- Comprehensive distance range analysis (0-2500 bp)

**Expected Results**:

- Multi-line plots showing CDF vs distance for different parameter values
- Clear parameter-dependent differences in termination profiles
- Quantitative comparison across parameter ranges

**When to use**:

- Systematic parameter exploration
- Understanding parameter effects on termination profiles
- Generating data for meta-analysis

---

## Output and Results

### Automatic Data Saving

All scripts use the unified `save_analysis_results.m` function that provides:

**File Organization**:

- Results saved to `SecondVersionResults/[ScriptName]/`
- Automatic folder creation
- Timestamped filenames with parameter information

**File Types Generated**:

1. **Plot Files** (`.png`): High-quality visualizations
2. **Data Files** (`_data.txt`): Complete numerical results with metadata

**Data File Contents**:

- Complete parameter sets used
- Timestamp and analysis type
- Numerical results in tab-delimited format
- Detailed comments explaining data structure

### Typical File Naming Convention

```
CPA_main_EBinding3_20241220_143022.png
CPA_main_EBinding3_20241220_143022_data.txt
Sweep1D_kHoff_EBinding5_20241220_143022.png
PASUsage_kHoff1e-03-1e-01_Etotal60000-140000_20241220_143022.png
```

---

## Usage Guidelines

### Getting Started

1. **Prerequisites**: MATLAB with Symbolic Math Toolbox
2. **Dependencies**: Ensure `compute_steady_states.m` and `ode_dynamics_multipleE.m` are in path
3. **Memory**: Some scripts use parallel processing - ensure adequate RAM

### Recommended Workflow

1. **Start with `CPA_multipleE_main.m`** to understand basic model behavior
2. **Use `parameter_sweep_1D.m`** to identify sensitive parameters
3. **Apply `PASUsageAnalysis.m`** for comprehensive 2D analysis
4. **Use specialized scripts** for specific research questions

### Parameter Modification

- Modify parameter values in the script headers
- Key parameters: `EBindingNumber`, `kHoff`, `E_total`, `kc`
- Biological ranges provided in comments

### Computational Considerations

- **Parallel scripts**: Will automatically use available CPU cores
- **Runtime**: Ranges from minutes (single runs) to hours (large sweeps)
- **Memory**: 2D sweeps can require several GB of RAM

---

## Troubleshooting

### Common Issues

**1. "Solver failed" errors**:

- Reduce parameter ranges or increase numerical tolerances
- Check for unphysical parameter combinations
- Verify all dependencies are in MATLAB path

**2. "Negative E_f" warnings**:

- Indicates parameter combinations that deplete E factors
- Reduce E_total or increase other parameters
- These combinations are automatically flagged as unphysical

**3. Slow performance**:

- Reduce parameter sweep ranges
- Close parallel pools between runs: `delete(gcp('nocreate'))`
- Increase MATLAB memory allocation

**4. Missing plots**:

- Check that figures aren't being closed prematurely
- Verify save directories have write permissions
- Ensure sufficient disk space

### Parameter Guidelines

**Biologically Reasonable Ranges**:

- `kHoff`: 0.001 - 0.1 (PAS strength)
- `E_total`: 10,000 - 200,000 (protein concentration)
- `kc`: 0.01 - 1.0 (cleavage rate)
- `EBindingNumber`: 1 - 8 (binding sites)

**Numerical Stability**:

- Avoid extreme parameter ratios (>1000x differences)
- Use logarithmic spacing for rate constants
- Start with default parameters and modify gradually

### Getting Help

- Check MATLAB console for detailed error messages
- Review saved data files for diagnostic information
- Verify model assumptions match your biological system

---

## Advanced Usage

### Customizing Analysis

- Modify distance ranges in PAS usage scripts
- Add new parameters to sweep lists
- Customize plot aesthetics and output formats

### Integration with Experiments

- Use `PASUsagevsInterPASDistance.m` results to predict APA-seq outcomes
- Compare `PASUsageAnalysis.m` phase diagrams with condition-dependent data
- Validate parameter ranges using experimental measurements

### Extending the Model

- Add new analysis types to `save_analysis_results.m`
- Create custom parameter sweep combinations
- Implement new metrics for termination analysis

---

_This guide covers the core functionality of the CPA Dynamics analysis suite. For specific biological interpretations or advanced customizations, consult the accompanying research papers and model documentation._
