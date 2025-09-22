# CPA Dynamics: A Kinetic Model of Transcription Termination and Alternative Polyadenylation

**Version 1.0**  
**Original Implementation**  
**Date: September 1, 2025**

## 1. Project Overview

This project provides a quantitative, kinetic model of RNA Polymerase II (Pol II) transcription termination in mammalian cells. The model is implemented in MATLAB and simulates the molecular competition that governs cleavage and polyadenylation (CPA), with a focus on predicting Alternative Polyadenylation (APA) site choice.

### Primary Goals

- Model the kinetic competition between elongation and termination at poly(A) sites
- Predict probability of proximal vs. distal poly(A) site usage
- Understand how CPA factor concentrations affect termination efficiency
- Provide quantitative framework for APA regulation

## 2. Model Architecture

The simulation is built on a detailed kinetic framework tracking multiple molecular species and their interactions along a discretized gene template.

### Core Molecular Species

- **R**: Active, elongating Pol II complexes (unbound)
- **RE**: Pol II complexes bound to **single** early CPA factor (E)
- **REL**: Pol II complexes bound to both early (E) and late (L) CPA factors, committed to termination

### Key Processes Modeled

1. **Pol II Initiation**: Constant rate initiation from promoter
2. **Elongation**: Position-dependent elongation rates (faster before PAS, slower after)
3. **Single E Factor Binding**: Position-dependent recruitment of **one E factor per polymerase**
4. **L Factor Binding**: Late CPA factor binding to RE complexes to form termination-competent REL
5. **Cleavage**: Final termination step from REL complexes

### **Critical Limitation: Single E Binding**

**Version 1.0 assumes each polymerase can bind at most one E factor**, creating a simple R → RE → REL pathway. This represents the minimal kinetic model for CPA factor recruitment but does not capture the cooperative binding effects observed experimentally.

## 3. Mathematical Framework

### ODE System Structure

The model uses a system of ordinary differential equations with three coupled species at each genomic position:

```matlab
% For each position n:
dR/dt = (flux_in) - (flux_out) - (E_binding) + (E_unbinding)
dRE/dt = (E_binding) - (E_unbinding) - (L_binding) + (L_unbinding) - (flux_out)
dREL/dt = (L_binding) - (L_unbinding) - (cleavage) - (flux_out)
```

### Conservation Laws

- **Total Pol II**: `Pol_total = sum(R) + sum(RE) + sum(REL) + Pol_free`
- **Total E Factors**: `E_total = sum(RE) + sum(REL) + E_free`
- **Total L Factors**: `L_total = sum(REL) + L_free`

### Position-Dependent Parameters

- **E Binding Rate**: `kE_on(n)` - increases linearly toward PAS, constant after PAS
- **Elongation Rate**: `k_e` before PAS, `k_e2` (slower) after PAS

## 4. Core MATLAB Files

### Primary Simulation Scripts

- **`CPA_odeSolver_main.m`**: Main simulation script for single-gene analysis
- **`CPA_odeSolver_multipleE.m`**: Extended version with multiple E binding sites
- **`ode_system.m`**: Core ODE system definition (R, RE, REL species)
- **`ode_system_multipleE.m`**: Extended ODE system with multiple E binding

### Analysis and Visualization

- **`CPA_backbone.m`**: Analysis of backbone termination profiles
- **`CPA_matrix.m`**: Matrix-based analysis methods
- **`MultipleEplot.m`**: Visualization of multiple E binding results
- **`projectionMiltipleE.m`**: Projection analysis for multiple E scenarios

### Testing and Validation

- **`sanity_check_model.m`**: Model validation and consistency checks
- **`Test2E2.m`**: Two-E binding scenario testing
- **`Test4E.m`**: Four-E binding scenario testing
- **`TestEf.m`**: E factor concentration testing

### Utility Functions

- **`ode_backbone.m`**: Backbone calculation utilities
- **`compute_normalized_ratios.m`**: Ratio analysis functions
- **`OneEFsolve.m`**: Single E factor solving routines

## 5. Key Parameters

### Kinetic Rate Constants

- **`k_in`**: Pol II initiation rate (molecules/time)
- **`k_c`**: Cleavage rate from REL complexes (1/time)
- **`k_e`**: Elongation rate before PAS (1/time)
- **`k_e2`**: Elongation rate after PAS (1/time, typically < k_e)
- **`kE_on`**: E factor binding rate (position-dependent)
- **`kE_off`**: E factor unbinding rate (1/time)
- **`kL_on`**: L factor binding rate (1/time)
- **`kL_off`**: L factor unbinding rate (1/time)

### Molecular Concentrations

- **`Pol_total`**: Total Pol II concentration (molecules)
- **`E_total`**: Total early CPA factor concentration (molecules)
- **`L_total`**: Total late CPA factor concentration (molecules)

### Geometric Parameters

- **`N`**: Total number of genomic positions (gene length / L_a)
- **`PAS`**: Position of poly(A) site (bp / L_a)
- **`L_a`**: Length per discretization unit (typically 100 bp)

## 6. Simulation Workflow

### 1. Parameter Setup

```matlab
% Define rate constants and concentrations
P.k_in = 2;
P.k_c = 0.2;
P.E_total = 70000;
% ... other parameters
```

### 2. Initial Conditions

- Start with zero or small positive values for all species
- Use ODE45 for transient simulation to approach steady state
- Refine with fsolve for exact steady-state solution

### 3. Steady-State Solution

```matlab
% Solve for steady state
options = optimoptions('fsolve', 'Display', 'iter');
[X, fval, exitflag] = fsolve(@(x) ode_system(0, x, P), X_init, options);
```

### 4. Analysis

- Extract R, RE, REL profiles along gene
- Calculate termination probability and read-through
- Analyze sensitivity to parameter changes

## 7. Key Features

### Biological Realism

- **Position-Dependent Rates**: E binding increases toward PAS
- **Multiple Factor Types**: Early (E) and late (L) CPA factors
- **Conservation Laws**: Explicit tracking of molecular conservation
- **Elongation Rate Changes**: Slower elongation after PAS recognition

### Mathematical Rigor

- **Steady-State Solutions**: Exact equilibrium calculations using fsolve
- **Transient Dynamics**: ODE45 integration for time-dependent behavior
- **Parameter Sensitivity**: Systematic exploration of parameter space
- **Convergence Checking**: Validation of numerical solutions

### Modular Design

- **Flexible Parameter Structure**: Easy modification of rate constants
- **Multiple E Binding**: Extension to complex binding scenarios
- **Testing Framework**: Comprehensive validation scripts
- **Visualization Tools**: Built-in plotting and analysis functions

## 8. Typical Results

### Termination Profiles

- **Sharp Termination**: High L factor concentrations lead to efficient termination shortly after PAS
- **Gradual Termination**: Low L factors result in extended termination regions
- **Read-Through**: Balance between termination and continued elongation

### Parameter Sensitivity

- **E Factor Concentration**: Higher E increases PAS recognition and termination preparation
- **L Factor Concentration**: Higher L increases termination efficiency
- **Cleavage Rate**: Higher k_c increases termination sharpness
- **Elongation Rates**: Ratio k_e/k_e2 affects termination window size

## 9. Limitations and Assumptions

### Model Simplifications

- **Single Gene**: No competition between multiple genes for factors
- **Single E Binding**: Each polymerase can bind at most one E factor (major limitation)
- **Homogeneous Factors**: All E factors assumed identical
- **Linear Geometry**: One-dimensional gene representation
- **Constant Initiation**: Steady promoter activity

### Key Limitation: Single E Factor Binding

The most significant limitation of Version 1.0 is the **single E binding assumption**. Experimental evidence suggests:

- **Cooperative Binding**: Multiple CPA factors can bind to single polymerases
- **Enhanced Termination**: Multiple E factors increase termination efficiency
- **Biological Realism**: Real CPA complexes involve multiple protein factors

This limitation was addressed in Version 2.0 with the **multiple E binding framework**.

### Numerical Considerations

- **Discretization**: Gene divided into finite segments (L_a = 100 bp typical)
- **Steady State**: Focus on equilibrium behavior
- **Convergence**: Requires careful initial guesses for complex parameter regimes

## 10. Future Extensions

### Potential Enhancements

- **Multiple Genes**: Competition for shared CPA factors
- **Stochastic Effects**: Noise and fluctuations in molecular counts
- **Chromatin Context**: Position-dependent accessibility and modifications
- **Alternative PAS**: Multiple poly(A) sites within single genes

### Technical Improvements

- **Parallel Computing**: Parameter sweeps with parfor
- **Advanced Solvers**: Specialized ODE solvers for stiff systems
- **Optimization**: Faster steady-state finding algorithms
- **Visualization**: Interactive parameter exploration tools

---

## File Organization (Version 1.0)

```
CPAdynamics/FirstVersion/
├── CPA_odeSolver_main.m          # Main single-gene simulation
├── CPA_odeSolver_multipleE.m     # Multiple E binding extension
├── ode_system.m                  # Core ODE system (R, RE, REL)
├── ode_system_multipleE.m        # Extended ODE system
├── CPA_backbone.m                # Backbone analysis
├── CPA_matrix.m                  # Matrix-based methods
├── MultipleEplot.m               # Multiple E visualization
├── projectionMiltipleE.m         # Projection analysis
├── sanity_check_model.m          # Model validation
├── Test2E2.m                     # Two E binding tests
├── Test4E.m                      # Four E binding tests
├── TestEf.m                      # E factor tests
├── ode_backbone.m                # Backbone utilities
├── compute_normalized_ratios.m   # Ratio calculations
└── OneEFsolve.m                  # Single E solver
```

This version represents the foundational implementation of the CPA dynamics model, providing a solid mathematical framework for understanding transcription termination kinetics and serving as the basis for subsequent enhancements and extensions.
