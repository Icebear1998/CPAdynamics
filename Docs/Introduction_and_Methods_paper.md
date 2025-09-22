## Introduction

### Background

Cleavage and polyadenylation (CPA) is the terminal step of eukaryotic mRNA biogenesis and a central regulator of gene expression. By recognizing sequence elements downstream of protein-coding regions, the CPA machinery cleaves nascent transcripts and adds poly(A) tails, enabling release of RNA Polymerase II (Pol II) and completion of transcription. A majority of mammalian genes (>70%) harbor multiple polyadenylation sites (PASs), leading to alternative polyadenylation (APA) and variable 3′-untranslated region (3′-UTR) lengths that influence mRNA stability, localization, translation, and protein function.

Despite detailed biochemical characterization of CPA components, a quantitative understanding of how kinetic competition between elongation and termination shapes APA site choice remains incomplete. In particular, the genome-wide sharing of CPA factors, the spatial organization of recruitment along genes, and the potential for cooperative binding create complex, non-linear behaviors that are not captured by sequence-based predictors alone.

### Problem statement

Existing models often treat PAS usage as a static, sequence-determined choice or rely on steady-state concentration ratios without explicitly modeling molecular flux. These approaches miss the inherently dynamic nature of termination: Pol II must both recognize the PAS and commit to cleavage while competing with continued elongation. Furthermore, the same pool of CPA factors services all actively transcribing genes, creating resource competition that couples termination efficiency across the genome.

### Key concept: termination commitment distance

We adopt the termination commitment distance (TCD) as a central quantitative descriptor. TCD is defined as the genomic distance downstream of a PAS at which a specified fraction (e.g., 75%) of polymerases have terminated. TCD captures both molecular kinetics and spatial organization, enabling comparisons across genes, conditions, and parameter regimes.

### Objective and contributions

This work presents a kinetic, multi-scale model of CPA dynamics capable of predicting APA outcomes and TCD. The main contributions are:

- **Multi-scale architecture**: Separation of slow Pol II motion (ODEs) from fast CPA factor binding (symbolic equilibrium) for tractable yet realistic simulations.
- **Cooperative binding**: A framework that allows multiple early CPA factors to bind a single polymerase, capturing cooperative effects absent from single-binding models.
- **Resource conservation**: Explicit accounting of free and bound molecular pools to model genome-wide resource competition.
- **Flux-based analysis**: Termination is quantified by actual molecular fluxes rather than steady-state concentrations.
- **Gene-length scaling**: A novel interpolation approach predicts gene-length effects on termination without simulating all genes simultaneously.

---

## Methods

### Model overview

We model transcription along a discretized gene with spatial step size `L_a` (typically 100 bp). The poly(A) site (PAS) resides at a specified index, and the domain extends downstream to capture post-PAS termination. The model comprises:

- **State variables**
  - `R(n)`: Concentration of elongating Pol II at position `n` (capable of binding multiple early CPA factors, E).
  - `REH(n)`: Concentration of Pol II complexes that have recognized PAS features and are committed to cleavage/termination (post-PAS dominant state).
- **Processes**
  - Initiation at the transcription start site (rate `k_in`).
  - Elongation before and after the PAS (rates `k_e` and `k_e2`, respectively, with `k_e2 ≤ k_e`).
  - Effective recognition/commitment to termination (conversion `R → REH` with position-dependent rate incorporating E binding).
  - Cleavage and release of `REH` (rate `k_c`).

This “Version 2.0” framework generalizes an earlier single-E binding model by allowing multiple E factors to bind Pol II, enabling cooperative recruitment and more realistic kinetics.

### Spatial discretization and boundary conditions

- Genes are divided into `N` nodes: `n = 1..N` with node length `L_a`.
- The PAS index `n = PAS` splits pre- and post-PAS regions.
- Boundary conditions:
  - Initiation: influx to `R(1)` at rate `k_in`.
  - Upstream boundary: no influx to `R(0)`.
  - Downstream removal: `R(N)` exits with rate `k_e`; `REH(N)` exits with rate `k_e2`.

### Multiple E-factor binding (fast equilibrium)

Early CPA factor binding is assumed to equilibrate rapidly relative to elongation. We pre-compute the average number of bound E factors at each position via symbolic analysis, producing an analytical function `P.RE_val_bind_E` that depends on:

- Free E concentration (`E_free`).
- Position-dependent affinity (reflecting CTD phosphorylation, local sequence, and chromatin context as modeled).
- Cooperative parameters and the maximum binding capacity (`EBindingNumber`).

This equilibrium contributes to the effective conversion rate from `R` to `REH` (PAS recognition/commitment), parameterized in the ODE system.

### ODE system (slow dynamics)

At each position `n`:

- Elongation transports `R` and `REH` between neighboring nodes with rates `k_e` (pre-PAS) and `k_e2` (post-PAS).
- Position-dependent recognition converts `R` to `REH` with an effective rate incorporating `P.RE_val_bind_E`.
- `REH` is cleaved at rate `k_c` and removed from the system.

Conservation equations track free and bound molecular pools, enforcing mass balance for Pol II and E factors:

- `Pol_total = Pol_free + Σ R(n) + Σ REH(n)`
- `E_total = E_free + E_bound(R, position)`

### Parameterization

Default parameters were chosen for biological plausibility and numerical stability and are tunable in analysis scripts:

- **Elongation**: `k_e = 65/L_a s⁻¹` pre-PAS; `k_e2 = 30/L_a s⁻¹` post-PAS.
- **Initiation and cleavage**: `k_in` (context dependent); `k_c ≈ 0.05 s⁻¹` baseline.
- **E binding**: Position-dependent on-rate embedded in the symbolic equilibrium; off-rates set by cooperative model parameters.
- **Pools**: Representative totals `Pol_total ≈ 7×10⁴`, `E_total ≈ 7×10⁴` (adjustable for sweeps).

Parameter ranges for sensitivity analyses typically span orders of magnitude (e.g., `E_total` 2×10⁴–2×10⁵; `k_c` 0.02–1.0 s⁻¹; effective binding rates 0.01–10 s⁻¹).

### Numerical solution strategy

We target steady-state solutions of the discretized ODEs:

- Initialize all state variables to small positive values to avoid singularities.
- Use MATLAB’s `fsolve` to solve for stationarity of the ODE right-hand sides defined in `ode_dynamics_multipleE.m`.
- Tight function tolerances (e.g., 1e−8) and bounded iteration counts ensure convergence.
- For difficult parameter regimes, transient integration can provide improved initial guesses before `fsolve` refinement.

All computations are vectorized where possible; parallelization (MATLAB Parallel Computing Toolbox) accelerates parameter sweeps.

### Flux-based termination and APA usage metrics

We quantify termination using fluxes rather than concentrations:

- Per-position cleavage flux: `flux_cleave(n) = k_c × REH(n)`.
- Exit fluxes at the gene end: `flux_R_exit = k_e × R(N)`, `flux_REH_exit = k_e2 × REH(N)`.
- Total outflux: `F_total = Σ flux_cleave + flux_R_exit + flux_REH_exit`.

The cumulative distribution of termination events downstream of PAS is the CDF:

- `CDF_terminate(n) = (Σ_{i=PAS..n} flux_cleave(i)) / F_total`.

- **Termination Commitment Distance (TCD)**: distance at which the CDF first reaches a preset threshold (e.g., 0.75). We compute TCD by interpolation over genomic distance.

For APA predictions with two PASs, proximal site usage at a given inter-PAS distance is the CDF evaluated at that distance downstream of the proximal PAS.

### Gene-length effects without multi-gene simulation

To capture genome-wide resource competition and gene-length dependence without explicitly simulating all genes simultaneously, we build interpolation functions from a grid of single-gene simulations across (`R_free`, `E_free`, `L`):

- Generate grids of occupied resources `R_occupied` and `E_occupied` with `generate_gene_length_grid.m`.
- Fit `scatteredInterpolant` functions mapping (`R_free`, `E_free`, `L`) to occupied pools with `build_gene_length_interpolation.m`.
- Enforce global conservation by solving for (`R_free`, `E_free`) that satisfy:
  - `R_total = R_free + ∫ R_occupied(R_free, E_free, L) f(L) dL`
  - `E_total = E_free + ∫ E_occupied(R_free, E_free, L) f(L) dL`

This self-consistent solution yields genome-wide free pools and enables prediction of gene-length–dependent TCD and APA profiles.

### Software and implementation

- **Language**: MATLAB (R2021a or later).
- **Core**: `ode_dynamics_multipleE.m`, `compute_steady_states.m`.
- **Analyses**: `CPA_multipleE_main.m`, `parameter_sweep_1D.m`, `parameter_sweep_2D.m`, `PASUsageAnalysis.m`, `PASUsagevsInterPASDistance.m`, `Sweep1DEbindingnumber.m`.
- **Utilities**: `calculate_pas_usage_profile.m`, `save_analysis_results.m`, grid/interpolation scripts for gene-length analysis.
- **Performance**: Vectorization and optional parallelization (`parfor`) are used to accelerate sweeps; numerical safeguards prevent negative pools and enforce conservation.

### Validation and robustness checks

- **Conservation**: Verify that total Pol II and E are conserved within numerical tolerance at steady state.
- **Sensitivity**: 1D and 2D parameter sweeps identify regimes of stable behavior and highlight nonlinearities due to cooperative binding and resource competition.
- **Comparisons to data**: Where available, compare TCD scales and APA usage trends to published datasets (e.g., termination distances inferred from RNA-seq; APA-seq proximal usage vs inter-PAS distance).
- **Numerical stability**: Diagnostics flag unphysical parameter combinations (e.g., negative `E_free`) and guide parameter range selection.

---

## Relevance and expected insights

By integrating cooperative E-factor binding, explicit resource conservation, and flux-based termination metrics, this framework connects molecular biophysics to gene-level APA outcomes and genome-wide competition. The TCD provides a compact, comparable measure of termination efficiency, while the interpolation-based scaling links gene length distributions to global CPA factor availability. Together, these components enable quantitative, testable predictions about how CPA regulation reshapes mRNA 3′-UTR landscapes across conditions.
