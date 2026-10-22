# **Methods**

## **A finite rate model of cotranscriptional CPA assembly**

We modeled how RNA polymerase II (Pol II) elongation, CTD phosphorylation, CPA-factor loading, PAS recognition, and RNA cleavage jointly determine the distance from a poly(A) signal (PAS) to completed assembly. We represented these processes as a coupled system of finite-rate reactions and transport along a discretized gene. 

Recognition was reversible … . Cleavage removed the polymerase from the modeled template and returned its bound resources to the free pools. We defined the primary readout, CPA assembly distance at 50% cleavage (${CAD}_{50}$),....

## **Gene discretization and polymerase microstates**

We represented each gene as a one-dimensional lattice with node spacing ${L}_{a}=100 bp$. For template length $G$ and the PAS coordinate is at ${x}_{PAS}$, we used $N=⌊G/{L}_{a}⌋$ nodes and PAS node $q=⌊{x}_{PAS}/{L}_{a}⌋$. Recognition and cleavage were permitted from node $q$ through $N$. The default single-gene template was 25 kb long, with its PAS 20 kb downstream of the transcription start site (TSS).

At node $i$, ${R}_{i}(p,e)$ denotes the population of unrecognized polymerases carrying $p$ phosphorylated Ser2P  sites and $e$ bound E factors, with $0\leq e\leq p\leq M$. Recognized polymerases were represented by ${RHE}_{i}(p,e)$, with $1\leq e\leq p\leq M$. Exactly one of the $e$ factors in an RHE state was engaged with H. 

## **CTD phosphorylation and E factor exchange**

Phosphorylation rate ${k}_{Pon}$ increased linearly along the gene with rate ${k}_{Pon,slope}$,  ${k}_{Pon}(i)={k}_{Pon,min}+{k}_{Pon,slope}(i-1)$. 

Dephosphorylation decreased $p$ by one at ${k}_{Poff}$ only when $p>e$, protecting E-occupied sites.

E associated with either R or RHE at $(p-e){k}_{Eon}{E}_{f}$, reflecting the number of unoccupied phosphorylated sites. Ordinary E \- E factor that does not in the recognition complex REH; dissociation occurred at $e{k}_{Eoff}$ in R and $(e-1){k}_{Eoff}$ in RHE. The recognition complex RHE dissociates in two ways described below.

## **PAS recognition and alternative dissociation routes**

PAS recognition was modeled as the transition from an elongating Pol II state (RE) to a PAS-bound REH state at the PAS. The effective PAS binding rate was obtained by scaling the basal rate kHon by the mean number of E factors loaded at the PAS. This coupling represents the assumption that Pol II molecules arriving with higher early-factor occupancy are more likely to enter a cleavage-competent, post-PAS complex. At and downstream of the PAS node, each bound E could initiate recognition, giving a transition $R(p,e)\rightarrow RHE(p,e)$ at $e{k}_{Hon}$. 

Recognition complex REH could dissociate via two distinct routes. H dissociation at ${k}_{Hoff}$ produced $RHE(p,e)\rightarrow R(p,e)$, retaining every CTD-bound E factor. Alternatively, the engaged E detached at ${k}_{Eoff,engaged}$, returning the polymerase to $R(p,e-1)$ and releasing one E. This second route assumed rapid disassembly of the detached EH complex. Because exactly one E was engaged, ${k}_{Eoff,engaged}$ had no $e$ multiplier. Remaining bound E factors could subsequently support recognition again.

| Event | Transition and rate per polymerase |
| :---- | :---- |
| Phosphorylation | $(p,e)\rightarrow (p+1,e)$ at ${k}_{Pon}(i)$, if $p<M$; in R and RHE |
| Dephosphorylation | $(p,e)\rightarrow (p-1,e)$ at ${k}_{Poff}$, if $p>e$; in R and RHE |
| E association | $(p,e)\rightarrow (p,e+1)$ at $(p-e){k}_{Eon}{E}_{f}$; in R and RHE |
| Ordinary E loss from R | $R(p,e)\rightarrow R(p,e-1)+{E}_{f}$ at $e{k}_{Eoff}$ |
| Ordinary E loss from RHE | $RHE(p,e)\rightarrow RHE(p,e-1)+{E}_{f}$ at $(e-1){k}_{Eoff}$ |
| PAS recognition | $R(p,e)\rightarrow RHE(p,e)$ at $e{k}_{Hon}$, for $i\geq q$ |
| H dissociation | $RHE(p,e)\rightarrow R(p,e)$ at ${k}_{Hoff}$ |
| Engaged E detachment | $RHE(p,e)\rightarrow R(p,e-1)+{E}_{f}$ at ${k}_{Eoff,engaged}$ |
| Cleavage | $RHE(p,e)\rightarrow {P}_{f}+e{E}_{f}$ at ${k}_{c}$ |

## **Initiation elongation and removal**

New polymerases entered node 1 in state $R(0,0)$ at flux ${k}_{in}{P}_{f}$. Elongation advanced a polymerase by one node without changing $p$ or $e$. R states moved at ${k}_{e}=(65 {bp s}^{-1})/{L}_{a}$, corresponding to 65 bp/s, and RHE states at ${k}_{e2}=(30 {bp s}^{-1})/{L}_{a}$, corresponding to 30 bp/s. Reduced velocity was associated with recognition, and unrecognized polymerases retained ${k}_{e}$ downstream of the PAS.

Cleavage removed RHE polymerases at ${k}_{c}$. At the final node, R and RHE polymerases could also leave by terminal run-off at their respective transport rates. Both cleavage and run-off returned one Pol II and all $e$ attached E factors to the free pools. Run-off represented escape from the finite simulated window and was not counted as cleavage.

## **Coupled kinetic equations and resource conservation**

Let $X$ collect all R and RHE microstate populations, ${E}_{f}$ denote free E, and ${P}_{f}$ denote free Pol II. The bound-state dynamics were written as

$$\frac{dX}{dt}=({A}_{0}+{E}_{f}{A}_{E})X+{P}_{f}b.$$

${A}_{0}$ contains phosphorylation, dephosphorylation, dissociation, recognition, cleavage, and transport. ${A}_{E}$ contains E-association transitions per unit free E. The source $b={k}_{in}u$ inflows polymerases into ${R}_{1}(0,0)$, where $u$ is a unit vector.

Free E was determined by enforcing conservation of total E. At steady state, the total E pool equals free E plus E bound across all R, RE, and REH states. Thus, the free E concentration used by the fast module was updated until local binding equilibria and global factor conservation were consistent. If $w$ contains the bound-E count $e$ of each microstate and $1$ denotes a vector of ones, conservation requires

$${E}_{total}={E}_{f}+{w}^{T}X,  {P}_{total}={P}_{f}+{1}^{T}X.$$

Free-pool derivatives were derived from the same reaction operator: $d{E}_{f}/dt=-{w}^{T}A({E}_{f})X$ and $d{P}_{f}/dt=-{1}^{T}A({E}_{f})X-{k}_{in}{P}_{f}$. Consequently, binding consumes free E, E detachment restores it, initiation consumes free Pol II, and cleavage and run-off recycle both resources.

## **CPA Assembly Distance**

We quantified CPA assembly efficiency using the CPA Assembly Distance (CAD), defined as the distance downstream of the PAS at which the cumulative commitment probability reaches 50%. Conceptually, CAD reports the downstream distance required for half of polymerases that encounter the PAS to commit to CPA. This definition was chosen because it maps naturally onto cis-antisense rescue assays, in which protection increases when CPA commitment occurs before downstream interference.

![][image1]

*Figure M3. CAD calculation from the commitment flux. The cleavage/commitment flux at downstream node l is k\_c · REH(l). CAD is obtained by normalizing the cumulative exit-flux distribution and interpolating the distance at which the cumulative probability reaches 50%.*

For each steady-state solution, the local commitment flux was computed as k\_c · REH(l). Run-off from the terminal node was treated as a terminal flux to ensure that polymerases reaching the end of the simulated template without commitment were included in the output distribution. The normalized cumulative sum of downstream commitment and terminal fluxes defined the termination CDF F(l), and CAD was obtained by linear interpolation at F(l) \= 0.5.

## **Gene-length and alternative polyadenylation analyses**

To test whether gene architecture affects CPA efficiency, we varied the TSS-to-PAS distance while holding the downstream model structure fixed and computed CAD for each simulated gene. To evaluate genome-wide resource competition, we precomputed a three-dimensional lookup grid for R\_occupied and E\_occupied as functions of free Pol II, free E, and gene length. These functions were interpolated and combined with a log-normal gene-length distribution to solve conservation equations for total Pol II and total E. This procedure links the single-gene assembly model to genome-wide competition for limiting transcription and CPA factors.

To model alternative polyadenylation (APA), we treated the termination CDF downstream of a proximal PAS as the probability that a polymerase commits to that PAS before reaching a distal PAS. Proximal PAS usage was estimated by evaluating the CDF at the inter-PAS distance. PAS strength was varied by sweeping the PAS dissociation parameter, allowing us to ask how intrinsic PAS strength and spacing jointly control proximal-site usage.

## **Implementation and parameters**

Simulations were implemented in MATLAB. The pipeline constructs the E-binding equilibrium rate matrix, solves the steady-state gene-scale ODE system with fsolve, and uses fzero to satisfy E-conservation after each ODE solve. Unless otherwise stated, figures were generated using the default parameter set described in the supplemental methods, which provide full parameter values and script-level reproducibility details.

[image1]: ![alt text](image.png)