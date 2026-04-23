# CPA Dynamics Project Summary

**CPA Dynamics** is a kinetic, multi-scale modeling framework implemented in MATLAB designed to simulate RNA Polymerase II (Pol II) transcription termination and Alternative Polyadenylation (APA) site choice in mammalian cells.

## 🎯 Core Objectives
The model fundamentally aims to describe how global resource availability (total Pol II, CPA complex proteins) and local binding dynamics (Ser2P phosphorylation, sequential polyadenylation sites) alter gene expression boundary events. Specifically, it predicts:
1. **Alternative Polyadenylation (APA) Choice:** Competition between Proximal and Distal Poly(A) sites for factor binding and cleavage.
2. **Termination Commitment Distance / CAD:** The distance downstream of a PAS required to drop active elongation flux by a measurable threshold (e.g., 50%).

## 🔬 Mathematical Approach
The framework uses a hybrid equilibrium-kinetic solver to process two vastly different timescales simultaneously:
- **Rapid Equilibrium:** Solves for pre-binding of global Cleavage and Polyadenylation (CPA) factors (dubbed the **`E`** factor) onto the Pol II CTD using pre-computed Symbolic math (`compute_steady_states.m`).
- **Kinetic Flux Evaluation:** Numerically solves the spatial, one-dimensional flux of Active Polymerase (`R`) converting into Cleaving/Terminating Polymerase (`REH`) moving linearly across discretely binned base pairs of the gene using `fsolve` nonlinear convergence constraints (`ode_dynamics_multipleE.m`).

### Version 2.0 Key Innovation: Multiple `E` Binding
While earlier models restricted Pol II to binding just one single `E` factor, the V2.0 architecture supports modeling multiple, cooperative binding states. Polymerases can now simultaneously bind `1` to `N` factors, allowing for more realistic response curves natively tied to substrate concentrations.

## 📁 Key Architecture & Files

### Core Engine
* **`run_termination_simulation.m`:** The central nervous system of the model. Orchestrates the 2-step solution: solving initial flux profiles given a free `E` concentration, and then forcing a subsequent re-evaluation balancing global protein conservation.
* **`ode_dynamics_multipleE.m`:** Stores the differential equations controlling elongation and transition rates across each spatial gene "node" (discretized lengths defined by `L_a`).
* **`compute_steady_states.m`:** Performs initial symbolic algebra to map rate dependencies into usable MATLAB function handles to speed up simulation times significantly.

### Analysis Workflows
* **`CPA_multipleE_main.m`:** Reconstructs the 1D gradient view of the model predicting Ser2P loading over transcription distances against bound E factors.
* **`PASUsagevsInterPASDistance.m`:** Maps out APA site choice (Proximal vs Distal usage vectors) as a function of the genomic base-pair distance spacing the poly(A) sites apart.
* **`parameter_sweep_1D.m` / `PASUsageAnalysis.m`:** Robust scripts for plotting sensitivity boundaries. Often configured to map shifts in termination distances against modifications to dissociation constants (`kEoff`, `kHoff`) or base pool concentrations (`E_total`, `Pol_total`).

## ⚙️ Standard Kinetic Parameters
Simulations within this dataset typically define biological interactions through standard rate profiles initialized via a global `P` struct:
* **Global Pooling:** e.g., `P.E_total = 100000;`, `P.Pol_total = 70000;`
* **Translocation & Processivity:** `P.k_e` (elongation rate), `k_in` (initiation rate).
* **Phosphorylation Scaling:** Modeled linearly as increasing downstream binding capacity defined by `P.kPon_min = 0.01;` extending via a slope `P.kPon_slope = 0.005;`.
* **Binding/Unbinding constants:** `kEon`, `kEoff` for CPA factor scaffolding, and `kHon`/`kHoff` for complex recognition at sequence signals.
