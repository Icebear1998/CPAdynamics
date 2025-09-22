# CPA Dynamics: A Kinetic Model of Transcription Termination and Alternative Polyadenylation

## Introduction

### Background

Transcription termination in mammalian cells is a highly regulated process that determines the length of mRNA 3'-untranslated regions (UTRs) and consequently affects gene expression, mRNA stability, and protein function. The cleavage and polyadenylation (CPA) machinery recognizes specific sequence elements downstream of protein-coding regions to terminate transcription and process pre-mRNA into mature mRNA molecules. Alternative polyadenylation (APA), the use of multiple polyadenylation sites within a single gene, occurs in over 70% of human genes and serves as a critical mechanism for post-transcriptional gene regulation.

The CPA process involves the coordinated action of multiple protein complexes that recognize the polyadenylation signal (PAS) and associated regulatory elements. Early CPA factors, including the cleavage and polyadenylation specificity factor (CPSF) and cleavage stimulation factor (CstF), bind to the elongating RNA polymerase II (Pol II) complex and recognize the PAS hexamer sequence (typically AAUAAA). Subsequently, late CPA factors, including the cleavage factors I and II (CFI and CFII), are recruited to form a complete termination-competent complex that cleaves the nascent RNA and releases the polymerase from the DNA template.

### The Challenge of Quantitative Modeling

Despite extensive biochemical characterization of individual CPA components, the kinetic competition that determines APA site choice remains poorly understood. Traditional approaches have focused on sequence-based prediction models or steady-state measurements, but these fail to capture the dynamic nature of the termination process. The challenge lies in modeling the multi-step kinetic pathway from PAS recognition to polymerase release, while accounting for the competition between continued elongation and termination commitment.

Several key factors complicate quantitative modeling of CPA dynamics:

1. **Multi-step kinetics**: The termination process involves sequential binding of multiple protein factors, each with distinct kinetic parameters.

2. **Spatial organization**: The efficiency of factor recruitment depends on the position along the gene and the local concentration of CPA proteins.

3. **Resource competition**: CPA factors are shared among all actively transcribing genes, creating genome-wide competition effects.

4. **Cooperative binding**: Multiple CPA factors can bind simultaneously to single polymerase complexes, leading to non-linear kinetic behavior.

### Termination Commitment Distance

A critical but poorly characterized parameter in transcription termination is the "termination commitment distance" (TCD) - the genomic distance downstream of a PAS at which a specified fraction of polymerases have terminated. TCD varies significantly among genes and cellular conditions, but its determinants and biological significance remain unclear. Understanding TCD is essential for predicting APA outcomes and designing therapeutic strategies that target 3'-UTR length.

### Objectives

This work presents a comprehensive kinetic model of CPA dynamics that addresses these challenges through several key innovations:

1. **Multi-scale modeling**: Separation of fast (factor binding) and slow (polymerase elongation) processes enables computational tractability while maintaining biological realism.

2. **Resource conservation**: Explicit tracking of free and bound CPA factor pools enables modeling of genome-wide competition effects.

3. **Multiple factor binding**: Unlike previous models limited to single factor interactions, our framework accommodates cooperative binding of multiple CPA factors to individual polymerases.

4. **Flux-based analysis**: Termination profiles are calculated based on actual molecular flux rather than steady-state concentrations, providing more accurate biological predictions.

5. **Gene length effects**: Novel mathematical framework for analyzing how gene length affects termination efficiency without explicit multi-gene simulation.

---

## Methods

### Model Architecture

#### Overview

The CPA dynamics model employs a hybrid, multi-scale approach that separates molecular processes occurring on different timescales. The core simulation tracks the movement and state transitions of RNA polymerase II complexes along discretized gene templates, while factor binding equilibria are pre-calculated using symbolic mathematics.

#### Spatial Discretization

Genes are discretized into nodes of length L_a (typically 100 bp), with node indices running from the transcription start site (TSS) to a specified distance beyond the polyadenylation site (PAS). The total number of nodes N is determined by the gene length, and the PAS position is specified as node index PAS.

#### Molecular Species

**Version 1.0 (Single E Binding Model)**:

- **R(n)**: Concentration of elongating polymerases at node n (unbound)
- **RE(n)**: Concentration of polymerases bound to single early CPA factor (E) at node n
- **REL(n)**: Concentration of polymerases bound to both early (E) and late (L) CPA factors at node n

**Version 2.0 (Multiple E Binding Model)**:

- **R(n)**: Concentration of elongating polymerases at node n (can bind multiple E factors)
- **REH(n)**: Concentration of polymerases committed to termination at node n (downstream of PAS only)

### Mathematical Framework

#### Version 1.0: Explicit Multi-Species Model

The single E binding model tracks three coupled species at each genomic position through the following system of ordinary differential equations:

**For positions before PAS (n ≤ PAS)**:

```
dR(n)/dt = k_e × R(n-1) - k_e × R(n) - k_E_on(n) × R(n) × E_free + k_E_off × RE(n)

dRE(n)/dt = k_E_on(n) × R(n) × E_free - k_E_off × RE(n) - k_e × RE(n) + k_e × RE(n-1)
           - k_L_on × RE(n) × L_free + k_L_off × REL(n)

dREL(n)/dt = k_L_on × RE(n) × L_free - k_L_off × REL(n) - k_e × REL(n) + k_e × REL(n-1)
```

**For positions after PAS (n > PAS)**:

```
dREL(n)/dt = ... - k_c × REL(n)  [includes cleavage term]
```

**Conservation Equations**:

```
E_free = E_total - Σ RE(n) - Σ REL(n)
L_free = L_total - Σ REL(n)
Pol_free = Pol_total - Σ R(n) - Σ RE(n) - Σ REL(n)
```

#### Version 2.0: Multiple E Binding with Symbolic Equilibrium

The multiple E binding model simplifies the molecular species while enabling cooperative factor binding through pre-calculated symbolic functions.

**Core ODE System**:

```
dR(n)/dt = k_in × δ(n,1) + k_e × R(n-1) - k_e × R(n) - k_H_on(n) × R(n) + k_H_off × REH(n)

dREH(n)/dt = k_H_on(n) × R(n) - k_H_off × REH(n) - k_c × REH(n) - k_e2 × REH(n) + k_e2 × REH(n-1)
```

**Multiple E Binding Framework**:
The average number of E factors bound to polymerases at position n is calculated using symbolic mathematics:

```matlab
P.RE_val_bind_E = matlabFunction(symbolic_expression_for_E_binding)
```

This function depends on:

- Local E factor concentration (E_free)
- Position-dependent binding affinity
- Cooperative binding parameters
- Maximum binding number (EBindingNumber)

#### Boundary Conditions and Initiation

**Initiation**: Polymerases enter at the TSS (n=1) with rate k_in proportional to free polymerase concentration.

**Elongation**: Polymerases move between adjacent nodes with rates k_e (before PAS) and k_e2 (after PAS, typically slower).

**Termination**: REH complexes undergo cleavage with rate k_c, removing them from the system.

**Gene End**: Polymerases reaching the gene end are removed with rates k_e (R species) and k_e2 (REH species).

### Parameter Estimation

#### Kinetic Rate Constants

**Elongation Rates**:

- k_e = 65/L_a s⁻¹ (before PAS)
- k_e2 = 30/L_a s⁻¹ (after PAS)

**CPA Factor Binding**:

- k_E_on: Position-dependent, increases linearly toward PAS
- k_E_off = 10 s⁻¹ (unbinding rate)
- k_H_on: Effective binding rate incorporating multiple E factors
- k_H_off = variable (depends on binding strength)

**Cleavage and Termination**:

- k_c = 0.05 s⁻¹ (cleavage rate)
- k_in = 2 s⁻¹ (initiation rate)

#### Molecular Concentrations

**Standard Conditions**:

- Pol_total = 70,000 molecules
- E_total = 70,000 molecules
- L_total = 100,000 molecules (Version 1.0 only)

**Parameter Ranges for Sensitivity Analysis**:

- E_total: 20,000 - 200,000 molecules
- k_c: 0.02 - 1.0 s⁻¹
- k_H_on: 0.01 - 10 s⁻¹

### Numerical Methods

#### Steady-State Solution

The model focuses on steady-state behavior, solved using MATLAB's `fsolve` function:

```matlab
options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
X = fsolve(@(x) ode_dynamics_multipleE(x, P), X0, options);
```

**Initial Conditions**: Small positive values (1e-6) for all species to avoid numerical singularities.

**Convergence Criteria**: Function tolerance of 1e-8 with maximum 1000 iterations.

#### Symbolic Mathematics

Version 2.0 employs MATLAB's Symbolic Math Toolbox to pre-calculate equilibrium binding expressions:

```matlab
syms E_free R_local real positive;
% Define binding equilibrium equations
% Solve for average E binding as function of E_free and R_local
P.RE_val_bind_E = matlabFunction(simplified_expression);
```

### Analysis Methods

#### Flux-Based Termination Analysis

Termination profiles are calculated based on molecular flux rather than concentrations:

**Total Outflux Calculation**:

```
flux_cleavage(n) = k_c × REH(n)
flux_R_exit = k_e × R(end)
flux_REH_exit = k_e2 × REH(end)
total_outflux = Σ flux_cleavage + flux_R_exit + flux_REH_exit
```

**Cumulative Distribution Function**:

```
exit_cdf(n) = Σ(i=1 to n) flux_cleavage(i) / total_outflux
```

This CDF represents the probability that a polymerase passing the PAS will terminate by position n.

#### Termination Commitment Distance (TCD)

TCD is defined as the distance downstream of PAS where a specified fraction (typically 75%) of polymerases have terminated:

```
TCD = interp1(exit_cdf, distances_bp, 0.75)
```

#### Parameter Sensitivity Analysis

**1D Parameter Sweeps**: Systematic variation of single parameters while holding others constant.

**2D Parameter Sweeps**: Exploration of parameter interactions through grid-based analysis.

**Parallel Computing**: Use of MATLAB's `parfor` for computational efficiency in parameter sweeps.

### Gene Length Analysis Framework

#### Motivation

To understand how gene length affects termination without explicit multi-gene simulation, we developed a novel interpolation-based approach that decouples local gene behavior from global resource pools.

#### Grid Generation

For various combinations of (R_free, E_free, L), where L represents TSS-to-PAS distance:

1. **Parameter Ranges**:

   - R_free: 100 - 10,000 molecules
   - E_free: 100 - 10,000 molecules
   - L: 2,500 - 200,000 bp

2. **Simulation**: Run single-gene model for each (R_free, E_free, L) combination

3. **Output**: Calculate R_occupied and E_occupied for each condition

#### Interpolation Functions

Build 3D interpolation functions using MATLAB's `scatteredInterpolant`:

```matlab
F_R_occupied = scatteredInterpolant(R_free_grid, E_free_grid, L_grid, R_occupied_grid);
F_E_occupied = scatteredInterpolant(R_free_grid, E_free_grid, L_grid, E_occupied_grid);
```

#### Conservation Equations

For a genome with gene length distribution f(L), the conservation equations become:

```
R_total = R_free + ∫ R_occupied(R_free, E_free, L) × f(L) dL
E_total = E_free + ∫ E_occupied(R_free, E_free, L) × f(L) dL
```

#### Self-Consistent Solution

Solve for (R_free, E_free) using `fsolve`:

```matlab
conservation_residual = @(x) [
    R_total - (x(1) + integral_R_occupied(x(1), x(2)));
    E_total - (x(2) + integral_E_occupied(x(1), x(2)))
];
[R_free_solution, E_free_solution] = fsolve(conservation_residual, initial_guess);
```

### Computational Implementation

#### Software Environment

- **Platform**: MATLAB R2021a or later
- **Required Toolboxes**: Symbolic Math Toolbox, Optimization Toolbox
- **Parallel Computing**: MATLAB Parallel Computing Toolbox (optional)

#### Code Organization

**Core Simulation**:

- `ode_dynamics_multipleE.m`: Main ODE system definition
- `compute_steady_states.m`: Symbolic equilibrium calculations

**Analysis Scripts**:

- `CPA_multipleE_main.m`: Single-gene analysis
- `parameter_sweep_1D.m`: Single-parameter sensitivity
- `parameter_sweep_2D.m`: Two-parameter interactions
- `PASUsageAnalysis.m`: Detailed termination analysis

**Gene Length Pipeline**:

- `generate_gene_length_grid.m`: Grid data generation
- `build_gene_length_interpolation.m`: Interpolation function construction
- `analyze_gene_length_TCD.m`: TCD analysis across gene lengths

**Utility Functions**:

- `calculate_pas_usage_profile.m`: Modular flux-based analysis
- `save_analysis_results.m`: Standardized data export

#### Performance Optimization

**Vectorization**: All calculations vectorized for MATLAB efficiency.

**Parallel Processing**: Parameter sweeps use `parfor` loops for multi-core execution.

**Memory Management**: Large datasets stored efficiently using sparse matrices where appropriate.

**Numerical Stability**: Careful handling of small concentrations and flux calculations to avoid numerical artifacts.

### Model Validation

#### Conservation Law Verification

All simulations verify mass conservation:

```matlab
total_check = sum(R_sol) + sum(REH_sol) + calculated_free_pools;
assert(abs(total_check - expected_total) < tolerance);
```

#### Parameter Sensitivity Bounds

Model behavior tested across physiologically relevant parameter ranges to ensure:

- Numerical stability
- Biologically reasonable predictions
- Convergence across parameter space

#### Comparison with Experimental Data

Model predictions compared against:

- Published APA site usage data
- Measured termination distances from RNA-seq
- CPA factor knockout/overexpression effects

This mathematical framework provides a comprehensive, computationally tractable approach to modeling CPA dynamics while maintaining the biological realism necessary for quantitative predictions of alternative polyadenylation outcomes.
