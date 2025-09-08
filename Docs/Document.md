A Kinetic Model of Transcription Termination and Alternative Polyadenylation (APA)
Version 1.0
Last Updated: September 1, 2025

1. Project Overview
   This project provides a quantitative, multi-scale model of RNA Polymerase II (Pol II) transcription termination in mammalian cells. The model is written in MATLAB and is designed to simulate the kinetic competition that governs cleavage and polyadenylation (CPA), with a specific focus on predicting Alternative Polyadenylation (APA) site choice.

The primary goal is to understand how both local features (e.g., the strength of a poly(A) signal hexamer) and global cellular factors (e.g., the concentration of CPA proteins) interact to regulate gene expression by determining the length of mRNA 3'-UTRs.

2. Model Architecture
   The simulation is built on a hybrid, multi-scale framework that separates processes occurring on different timescales.

a. Slow Process: Pol II Elongation (Main ODE)
This is the core of the simulation, modeled as a system of Ordinary Differential Equations (ODEs). It describes the movement of Pol II along a discretized gene template. The key states tracked are:

R: The population of active, elongating Pol II complexes.

REH: The population of terminating Pol II complexes that have recognized the poly(A) signal (Hexamer) and are committed to cleavage.

b. Fast Process: CPA Factor Binding (Embedded Symbolic Model)
The binding of early CPA factors (denoted 'E') to the Pol II CTD is assumed to be a rapid equilibrium. Instead of simulating this in the main ODE, it is pre-calculated.

The model uses MATLAB's Symbolic Math Toolbox to derive an analytical function, P.RE_val_bind_E, that gives the average number of 'E' factors bound to Pol II at any position on the gene.

This binding is dependent on the local phosphorylation state of the Pol II CTD, which changes as the polymerase moves along the gene.

3. Key Features & Capabilities
   Parameter Sweeps: The code is structured to perform large-scale parameter sweeps to test the sensitivity of the system to different biochemical rates and protein concentrations.

APA Prediction: The model can predict the probability of using a proximal vs. a distal poly(A) site by analyzing the flux of polymerases that terminate versus those that read through the first site.

Biophysically Grounded: Key parameters like k
Hon
​
and k
Hoff
​
are estimated using first-principles biophysical calculations, grounding the model in physical reality.

4. Core MATLAB Files
   SweepParameterPASusage_parallel.m: The main script for running parameter sweeps. This is the primary entry point for simulations.

run_termination_simulation_parallel.m: A helper function that runs a single simulation for a given set of parameters.

ode_dynamics_multipleE_parallel.m: Defines the system of Ordinary Differential Equations for Pol II elongation and termination.

compute_steady_states.m: The script that performs the symbolic
