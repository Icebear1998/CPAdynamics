# CPA Dynamics: A Kinetic Model of Transcription Termination and Alternative Polyadenylation

**Version 2.0**  
**Last Updated: September 17, 2025**

## 1. Project Overview

This project provides a quantitative, multi-scale model of RNA Polymerase II (Pol II) transcription termination in mammalian cells. The model is implemented in MATLAB and simulates the kinetic competition that governs cleavage and polyadenylation (CPA), with a focus on predicting Alternative Polyadenylation (APA) site choice and termination commitment distance (TCD).

### Primary Goals

- Understand how local features (poly(A) signal strength) and global factors (CPA protein concentrations) regulate gene expression
- Predict mRNA 3'-UTR length determination through APA site choice
- Analyze termination commitment distance as a function of gene length and resource availability
- Model genome-wide resource competition effects on individual gene termination

## 2. Model Architecture

The simulation uses a hybrid, multi-scale framework separating processes on different timescales:

### a. Slow Process: Pol II Elongation (Main ODE System)

Core simulation modeled as Ordinary Differential Equations describing Pol II movement along discretized gene templates.

**Key State Variables:**

- **R**: Active, elongating Pol II complexes (can bind multiple E factors)
- **REH**: Terminating Pol II complexes that have recognized poly(A) signals and are committed to cleavage

**Critical Innovation: Multiple E Factor Binding**  
Unlike Version 1.0's single E binding limitation, **Version 2.0 allows each polymerase to bind multiple E factors simultaneously**, enabling cooperative binding effects and more realistic termination kinetics.

### b. Fast Process: CPA Factor Binding (Symbolic Equilibrium Model)

Early CPA factor ('E') binding to Pol II CTD assumed to be rapid equilibrium, pre-calculated using MATLAB's Symbolic Math Toolbox.

**Output**: `P.RE_val_bind_E` - analytical function for average 'E' factors bound to Pol II at any gene position, dependent on local CTD phosphorylation state.

## 3. Key Features & Capabilities

### Core Analysis Types

- **Parameter Sweeps**: Large-scale sensitivity analysis across biochemical rates and protein concentrations
- **APA Prediction**: Probability calculations for proximal vs. distal poly(A) site usage
- **Termination Commitment Distance**: Novel metric quantifying termination efficiency as function of gene length
- **Resource Competition Modeling**: Genome-wide effects without explicit multi-gene simulation

### Mathematical Framework

- **Flux-Based Termination Analysis**: CDF calculation based on actual termination events
- **Conservation Equations**: Self-consistent solutions for free vs. bound resource pools
- **Interpolation-Based Scaling**: 3D interpolation functions for efficient parameter space exploration
- **Biophysically Grounded Parameters**: First-principles estimation of key kinetic rates

## 4. Analysis Scripts

### Primary Analysis Scripts

- **`CPA_multipleE_main.m`**: Core multi-E binding analysis with termination profiles
- **`parameter_sweep_1D.m`**: Single-parameter sensitivity analysis with flux-based cutoff calculation
- **`parameter_sweep_2D.m`**: Two-parameter interaction analysis
- **`PASUsageAnalysis.m`**: Detailed PAS usage analysis across parameter ranges
- **`SweepParameterPASusage.m`**: Parameter sweeps focused on PAS usage patterns
- **`PASUsagevsInterPASDistance.m`**: Analysis of inter-PAS distance effects
- **`Sweep1DEbindingnumber.m`**: E-binding number sensitivity analysis

### Gene Length Analysis Pipeline

- **`generate_gene_length_grid.m`**: Creates 3D grid data for R_occupied and E_occupied as functions of (R_free, E_free, L)
- **`build_gene_length_interpolation.m`**: Builds interpolation functions and visualizes resource occupation patterns
- **`analyze_gene_length_TCD.m`**: Solves conservation equations and calculates termination commitment distance

### Utility Functions

- **`calculate_pas_usage_profile.m`**: Modular function for flux-based termination profile calculation
- **`save_analysis_results.m`**: Standardized data and plot saving across all analysis scripts
- **`ode_dynamics_multipleE.m`**: Core ODE system definition
- **`compute_steady_states.m`**: Symbolic calculation of steady-state solutions

## 5. Key Innovations

### Version 2.0 Updates

- **Multiple E Factor Binding**: **Major advancement from Version 1.0** - each polymerase can now bind multiple E factors simultaneously, enabling cooperative binding effects
- **Modular Architecture**: Extracted common calculations into reusable utility functions
- **Flux-Based Analysis**: Replaced concentration ratios with biologically accurate flux calculations
- **Gene Length Modeling**: Novel approach to model genome-wide competition without multi-gene simulation
- **Standardized Output**: Consistent data saving and organization across all analysis types
- **Resource Conservation**: Mathematical framework for self-consistent free/bound resource determination

### Key Difference from Version 1.0: Multiple E Binding

**Version 1.0 Limitation**: R → RE (single E) → REL pathway  
**Version 2.0 Innovation**: R can bind multiple E factors (RE₁, RE₂, RE₃, ..., REₙ), allowing:

- **Cooperative Effects**: Enhanced binding affinity with multiple factors
- **Realistic Kinetics**: Better match to experimental CPA complex formation
- **Tunable Termination**: Variable E binding numbers (EBindingNumber parameter)

### Termination Commitment Distance (TCD)

New metric quantifying the distance downstream of PAS where a specified fraction of polymerases have terminated. Enables:

- Comparison across different gene lengths
- Resource competition effects quantification
- Genome-wide termination efficiency prediction

## 6. File Organization

```
CPAdynamics/
├── Analysis Scripts/          # Primary analysis workflows
├── FirstVersion/             # Legacy implementations
├── Results/                  # Generated plots and data
├── SecondVersionResults/     # Current version outputs
├── Docs/                    # Documentation
└── Utility Functions/        # Modular helper functions
```

### Data Management

- **Automated Saving**: All scripts use `save_analysis_results.m` for consistent output
- **Organized Results**: Analysis-specific folders with timestamped files
- **Raw Data Export**: Text files with parameter values and numerical results
- **Plot Generation**: High-quality figures with metadata and timestamps
