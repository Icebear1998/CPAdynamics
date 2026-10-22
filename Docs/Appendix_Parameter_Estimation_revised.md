# **Appendix: Parameter Estimation**

## **A. Molecular Association Rate ${K}_{on}$**

Throughout this estimation, we use a basal molecular association rate of ${K}_{on}$ \= $5 \times 10^{5}$ M⁻¹ s⁻¹ for bimolecular interactions. This is a basal, orientation-limited association rate: protein–protein association without electrostatic steering typically proceeds at \~$5 \times 10^{5}$ M⁻¹s⁻¹, orders of magnitude below the Smoluchowski diffusion limit, because the partners must meet in a restricted set of orientations (Dhusia & Wu, 2021; Northrup & Erickson, 1992; Schreiber et al., 2009; Zhou, 1993).

The exact value of ${K}_{on}$ cancels from the dimensionless equilibrium ratios. For the hexamer-binding module, ${k}_{Hd}$ \= ${k}_{Hoff}$ / ${k}_{Hon}$, and for the E-factor CTD-binding module, ${k}_{Ed}$ \= ${k}_{Eoff}$ / ${k}_{Eon}$. Because the corresponding forward and reverse microscopic rates are scaled consistently from the same ${K}_{on}$ convention, ${K}_{on}$ cancels from these ratios, so equilibrium occupancies do not depend on it. It does not cancel from the absolute rates. ${k}_{Hon}$, ${k}_{Hoff}$, ${k}_{Eon}$, ${k}_{Eoff}$ and ${k}_{c}$ all scale with ${K}_{on}$ and compete with the elongation rates ${k}_{e}$ and ${k}_{e2}$, so ${K}_{on}$ sets how far Pol II travels before recognition and commitment complete. Simulation behavior is therefore controlled by ${k}_{Hd}$ and ${k}_{Ed}$ together with the absolute rates relative to elongation, including the commitment rate ${k}_{c}$, all of which are tied to this ${K}_{on}$ convention.

## **B. First-Principles Estimation of the Hexamer Binding Rate**

### **E-factor PAS Recognition Complex — Hexamer Binding (${k}_{Hon}$)**

The PAS recognition rate (${k}_{Hon}$) represents the tethered encounter between early CPA factors (E factors) pre-loaded on the Pol II CTD and the newly transcribed poly(A) signal hexamer. The relevant geometry includes two flexible tethers whose proximal ends attach to the Pol II body: the CTD arm and the nascent-RNA arm.

Because the CTD can be represented as an intrinsically disordered polymer, the worm-like chain (WLC) model approximates this problem well (Caliskan et al., 2005; Hyeon et al., 2006; Phillips et al., 2013).

**Physical parameterization.** We parameterize the model with a representative baseline and a physically motivated range:

* **Effective Contour Length (${L}_{eff}$):** Baseline 107 nm. The human Pol II CTD repeat array contains 52 heptad repeats, corresponding nominally to 364 residues and a fully extended contour length of approximately 131–138 nm (Corden et al., 1985; Rosonina & Blencowe, 2004; Hsin & Manley, 2012). The CTD is joined to the Pol II core by a flexible, unstructured linker of \~86 residues in yeast (Sdano et al., 2017). The human counterpart is the RPB1 stretch between the core and the repeats, which Vos et al. (2018) expressed as residues 1488–1592 (105 residues). That paper did not measure its length or flexibility, and part of it may be structured: UniProt (P24928) annotates predicted disorder only from residue 1546. We therefore treat \~105 residues as an upper bound on the flexible stretch. Treating it as a worm-like chain with \~0.36–0.38 nm per residue (our assumption, not a value from Vos et al.) gives a contour length of \~38–40 nm. If only residues 1546–1592 were flexible, the contour length would be \~17–18 nm, ${L}_{eff}$ would be \~86 nm, and ${k}_{Hon}$ would rise from 7.0 to 7.9 s⁻¹, within the range below. Assuming that E factors bind uniformly along the repeat array, their mean contour coordinate is therefore approximately

$${L}_{eff} \approx  38 + \frac{138}{2} \approx  107 nm$$

* We consequently use ${L}_{eff} = 107$ nm as the representative effective contour length

* **Adjusted Persistence Length (${l}_{p}$):** Baseline 2.0 nm. Range: \[1.0, 3.0\] nm. A naked, highly flexible IDP may have a short persistence length (\< 1.0 nm), whereas the phosphorylated CTD is expected to be more extended and stiffer due to its high negative charge (Portz et al., 2017).

* **Anchor Separation (d):** The distance between two anchors is estimated at d \= 10 nm, the distance from the CTD linker to the exit channel of the RNA.

**Effective concentration.** For independent, long, flexible chains, we approximate the two endpoint distributions as Gaussians. Their convolution gives a relative-coordinate Gaussian with mean-square span ${⟨{R}^{2}⟩}_{tot}$ \= ${⟨{R}^{2}⟩}_{CTD}$ \+ ${⟨{R}^{2}⟩}_{RNA}$. Using the 3D Gaussian probability density for a WLC end-to-end vector, the effective local concentration, J-density, at distance d from the anchor (Phillips et al., 2013\) is:

$${J}_{density}={\left({\frac{3}{2\pi ⟨{R}^{2}{⟩}_{tot}}}\right)}^{3/2}exp\left({\frac{-3{d}^{2}}{2⟨{R}^{2}{⟩}_{tot}}}\right)$$

**CTD arm:** For the ideal WLC, ${⟨{R}^{2}⟩}_{CTD}$ ≈ 2 × ${L}_{eff}$ × ${l}_{p}$ \= 428 nm², corresponding to R ≈ 20.69 nm. This is of the same order as the SAXS measurements of the phosphorylated Drosophila CTD (Portz et al., 2017). In this report, the radius of gyration, ${R}_{g}$ \= 6.06 nm, is for the 42-repeat fly CTD. ${R}_{g}$ is the root-mean-square distance of all chain segments from the chain’s center of mass, and R is the end-to-end distance. Using the relationship R ≈ 2.5 ${R}_{g}$, which holds approximately for unfolded or disordered protein ensembles (Riback et al., 2019), the phosphorylated fly CTD gives R ≈ 15.2 nm. Scaling this value from the 42-repeat fly CTD to the 52-repeat human CTD as √(52/42) yields:

$$R \approx  16.9 nm (SAXS-scaled); baseline WLC estimate: 20.7 nm$$

The 16.9 nm SAXS value is for the whole CTD without the linker, whereas the baseline WLC estimate uses the mean contour coordinate of an E factor, which includes the linker. Without the linker, the mean contour coordinate would be 138/2 = 69 nm and the WLC estimate would be R ≈ √(2 × 69 nm × 2 nm) ≈ 16.6 nm, similar to the SAXS-scaled value.

**RNA arm:** For 200–600 nt, the contour-length range is ${L}_{RNA}$ \= 200 × 0.61 to 600 × 0.61 \= 122–366 nm. The baseline 400-nt segment has ${L}_{RNA}$ \= 244 nm. With ${l}_{p,RNA}$ \= 1.26 nm, ${⟨{R}^{2}⟩}_{RNA}$ ≈ 2 ${l}_{p}$ L \= 614.88 nm² at baseline.

For our baseline ${{R}^{2}}_{tot}$ \= 428 \+ 614.88 \= 1042.88 nm² and d \= 10 nm, this yields ${J}_{density}$ ≈ 8.484 × 10⁻⁶ nm⁻³. Using the unit conversion 1 nm⁻³ = 1/(N_A × 10⁻²⁴ L) ≈ 1.66 M:

$${J}_{density} \approx  1.409 \times  10⁻⁵ M = 14.09 \mu M$$

Applying our baseline association rate (${K}_{on}$ \= 5 × 10⁵ M⁻¹s⁻¹):

$${k}_{Hon} = {K}_{on} · {J}_{eff}$$

$${k}_{Hon} \approx  5 \times  {10}^{5} M⁻¹ s⁻¹ \times  1.409 \times  {10}^{-5} M = 7 s⁻¹$$

**Range evaluation:**

$${k}_{Hon} \in  [4, 17] s⁻¹$$

### **E-factor PAS Recognition Complex — Hexamer Dissociation (${k}_{Hoff}$)**

The equilibrium dissociation constant ${K}_{d}$ for polyadenylation signal (PAS) hexamer recognition by mPSF (comprising CPSF-160, WDR33, CPSF-30, and Fip1) is approximately 3–40 nM for canonical AAUAAA with flanking sequences (Hamilton et al., 2019; Huang et al., 2025\) and roughly 2,000 nM for the isolated hexamer (Hamilton et al., 2019). Variant hexamers exhibit weaker binding, ranging from 79 nM to 13 μM (Huang et al., 2025). We used two distinct values for ${k}_{Hoff}$: a strong PAS in the 3–40 nM range and a weak PAS near 2,000 nM, consistent with Chao et al. (1999), who compared strong (SV40 late) and weak (SV40 early) PASs. The microscopic dissociation rates derived using our nuclear association rate are as follows:

${k}_{Hoff} = {K}_{on} · {K}_{d} = 5 \times  {10}^{5} M⁻¹ s⁻¹ \times  20 \times  {10}^{-9} M = 0.01 s⁻¹$ (strong PAS, middle of 3–40 nM)

${k}_{Hoff} = {K}_{on} · {K}_{d} = 5 \times  {10}^{5} M⁻¹ s⁻¹ \times  2000 \times  {10}^{-9} M = 1 s⁻¹$ (weak PAS)

**Range evaluation:**

$${k}_{Hoff} \in  [0.01, 1] s⁻¹$$

## **C. E-factor Binding Parameters**

### **Early Factor Association Rate (${k}_{Eon}$)**

E-factor binding to the Pol II CTD is assumed to proceed at the basal association rate ${K}_{on}$. We convert the molar rate constant into a per-molecule rate suitable for physiologically scaled simulations using:

$${k}_{Eon}=\frac{{K}_{on}}{{N}_{A}\cdot {V}_{nucleus}}$$

$${k}_{Eon} = \frac{5 \times  {10}^{5} M⁻¹ s⁻¹}{6.022 \times  {10}^{23} \times  3.74 \times  {10}^{-13} L} \approx  2.2 \times  {10}^{-6} molecule⁻¹ s⁻¹$$

where ${N}_{A}$ \= 6.022 × 10²³ molecules/mol and ${V}_{nucleus}$ \= 3.74 × 10⁻¹³ L corresponds to a representative 374 μm³ mammalian nucleus (BioNumbers). The formula gives 2.2 × 10⁻⁶; the simulations use the rounded working value ${k}_{Eon}$ = 2.5 × 10⁻⁶ molecule⁻¹ s⁻¹ (equivalent to a 332 μm³ nucleus), which lies inside the range below.
**Range evaluation:** Over the biological range of nuclear volumes \[300, 500\] μm³, ${k}_{Eon}$ ∈ \[1.66 × 10⁻⁶, 2.77 × 10⁻⁶\] molecule⁻¹s⁻¹.

### **Early Factor Dissociation Rate (${k}_{Eoff}$)**

The intrinsic binding affinity ${K}_{d,E}$ of an early CPA factor for a single phosphorylated Ser2-Pro3 CTD heptad repeat is measured in vitro at 1–10 μM (Mei et al., 2014; Meinhart & Cramer, 2004; Ni et al., 2014; Noble et al., 2005). These are affinities of the CID domains of CTD-binding factors (Pcf11, RPRD1A/1B, CREPT/p15RS), not of CPSF itself; we use them as a proxy for E–CTD binding. Using ${k}_{off}$ \= ${K}_{on}$ · ${K}_{d,E}$ for a single repeat:

$${k}_{Eoff,single} = (1–10) \times  {10}^{-6} M \times  5 \times  {10}^{5} M⁻¹ s⁻¹ = 0.5–5 s⁻¹$$

**Multivalent correction.** The E factor engages 2–4 tandem CTD repeats, as suggested by co-crystal structures of CID–CTD interactions (Meinhart & Cramer, 2004; Noble et al., 2005). Multivalent or bivalent engagement of a disordered CTD polymer substantially increases the effective binding affinity relative to the single-repeat value, which can push the effective ${k}_{Eoff}$ to the lower end of the range.

$${k}_{Eoff}^{eff}\approx 0.5{ s}^{-1}$$

This correction is qualitative rather than derived. ${k}_{Eoff}$ = 0.5 s⁻¹ corresponds to ${K}_{d,E}$ = ${k}_{Eoff}$/${K}_{on}$ = 1 μM, the lower end of the single-repeat range, and we adopt it as an effective multivalent value.

**Range evaluation:** 

${k}_{Eoff}$ in \[0.1, 5\] s⁻¹. 

### **Engaged-E Detachment Rate (${k}_{Eoff,engaged}$)**

In our model, we assumed that in a recognized complex, one CPA factor complex (one early factor (E)) bridges the CTD and the hexamer on the nascent RNA. Once they form the binding complex, the binding of E to the ser2P site on the CTD is assumed to be tighter. It therefore has a lower kEoff rate than an ordinary E (${k}_{Eoff}$ \= 0.5 s⁻¹). How firmly H holds the released E has not been measured. If E stayed bound to the hexamer, the rate would be far below ${k}_{Eoff}$. If E lost H right after release from the CTD, it would equal ${k}_{Eoff}$. We take a conservative value between these limits:

${k}_{Eoff,engaged}$ **\= 0.1 s⁻¹** (one-fifth of ${k}_{Eoff}$)

**Range evaluation:** ${k}_{Eoff,engaged}$ ∈ \[0.01, 0.5\] s⁻¹. The upper bound is close to ${k}_{Eoff}$. The lower bound allows a firm hold by H.

## **D. CPA Commitment Rate ${k}_{c}$**

The CPA commitment rate ${k}_{c}$ represents the productive, effectively irreversible assembly of the complete CPA machinery after early PAS recognition. This step includes diffusional recruitment of late CPA factors, such as CstF; productive engagement with the downstream sequence element (DSE); correct steric/orientational alignment; and assembly of the remaining components needed to form a complete cleavage-competent complex. Therefore, ${k}_{c}$ should be interpreted as an effective commitment rate rather than a raw diffusion-limited collision rate.

**Late CPA factor copy number.** We estimate the global nuclear copy number of late-processing factors to be approximately 300,000 molecules per nucleus (see Section G). Using our baseline nuclear volume of 3.74 × 10⁻¹³ L:

$$\left[{{L}_{free}}\right]=\frac{300,000}{6.022\times {10}^{23}\times 3.74\times {10}^{-13}}\approx 1.33 \mu M$$

(We take ${[L]}_{free}$ = ${L}_{total}$/(${N}_{A}{V}_{nucleus}$), so the late-factor pool is treated as non-limiting: we assume that assembly on other pre-mRNAs and binding at other cellular sites do not deplete it. (We over estimate it here))

Late CPA factors are multi-subunit assemblies built on dimerization: the CstF hexamer is a dimer of CstF-77/64/50 trimers, and the CFIm tetramer is a CFIm25 dimer clasped by two CFIm68/59 (Yang et al., 2011; Yang et al., 2018). Biochemically, assembly on the pre-mRNA is ordered and cooperative: CPSF recognizes the PAS first, CstF then joins, and the resulting complex is stable and committed to processing (Gilmartin & Nevins, 1989, 1991). Given this ordered, multi-component assembly, we assume that complete machinery assembly following early PAS recognition requires at least three sequential steps. The number of steps is a modelling assumption. The scheme is a coarse-grained representation, not a literal sequence of trimer binding events: L is a generic late-factor binding unit whose concentration is taken from the CstF pool (Section G), CFIm is not counted separately, and the third step (REHL₂ → commit) is given the same rate r for simplicity. Let R \= \[REH\], X \= \[REHL\], and Y \= \[REHL₂\]. Define r \= ${K}_{on}$${[L]}_{free}$ and q \= ${k}_{off}$ \= ${K}_{d,L}.{K}_{on}$ . We assume all three forward steps have rate r, and the first two steps reverse sequentially with rate q:

$$REH  {{⇌}^{r}}_{q}  REHL  {{⇌}^{r}}_{q}  REHL₂  {\rightarrow }^{r}  commit$$

The fraction of initial capture events that reaches commitment is:

$${f}_{lock} = \frac{{J}_{commit}}{{J}_{capture}} = \frac{{r}^{2}}{{q}^{2} + qr + {r}^{2}} = \frac{{[L]}^{2}}{{{K}_{d,L}}^{2} + {K}_{d,L}[L] + {[L]}^{2}}$$

The corresponding effective commitment rate is:

$${k}_{c} = r{f}_{lock} = {K}_{on}[L] \frac{{[L]}^{2}}{{{K}_{d,L}}^{2} + {K}_{d,L}[L] + {[L]}^{2}}$$

Using ${[L]}_{free}$ \= 1.332 μM and ${K}_{d,L}$ \= 2 μM gives ${k}_{c}$ \= 0.14 s⁻¹.

**Approximation.** ${k}_{c}$ = r${f}_{lock}$ is the commitment flux per REH when REH is held constant. It ignores the depletion of REH into REHL and REHL₂, and we keep it for simplicity. This underestimates the mean commitment time and so overestimates ${k}_{c}$.

**Range evaluation:** Varying the nuclear volume over \[300, 500\] μm³ and the late-factor copy number over \[200,000, 400,000\] gives ${[L]}_{free}$ ≈ 0.664–2.214 μM. The ${K}_{d,L}$ interval is supported by: \- **Tight end:** the CstF hexamer on GU-rich RNA, 0.12–0.9 μM, and CstF-77/64 on GU14, 0.16 μM (Yang et al., 2018); CFIm25/68 on two UGUA motifs, 0.22 μM (Yang et al., 2011). \- **Weak end:** the isolated CstF-64 RRM, 0.7–2.5 μM, and CstF-77/64 on a U-rich DSE, 3.7 μM (Yang et al., 2018; Grozdanov et al., 2018). Using ${K}_{d,L}$ ∈ \[0.1, 4\] μM at the baseline ${[L]}_{free}$ = 1.332 μM gives

${k}_{c}$ ≈ 0.05–0.62 s⁻¹.


## **E. Pol II Elongation Rates (${k}_{e}$, ${k}_{e2}$)**

### **Pre-PAS Elongation Rate (${k}_{e}$)**

Pol II elongation velocity in human cells has been measured by live-cell imaging of a reporter gene (Darzacq et al., 2007, 4.3 kb/min) and by metabolic labeling of nascent RNA with sequencing (Fuchs et al., 2014, most genes about 3.5 kb/min, range 2–6 kb/min). Together these support values of 3–5 kb/min. We use a central estimate of 4 kb/min, corresponding to 4,000 bp / 60 s ≈ 66.67 bp/s.

In the model, the gene is discretized into nodes of ${L}_{a}$ \= 100 bp, so the elongation rate in model units is:

$${k}_{e} = \frac{66.67 bp/s}{100 bp/node} = 0.667 node/s$$

The simulations use the rounded working value ${k}_{e}$ = 0.65 node s⁻¹ (65 bp/s, about 3.9 kb/min).

**Range: ${k}_{e}$** ∈ \[0.50, 0.83\] node s⁻¹, corresponding to elongation velocities of 3–5 kb/min.

### **Post-PAS Elongation Rate (${k}_{e2}$)**

After PAS recognition, Pol II slows down and pauses near the gene end: mNET-seq shows Pol II pausing at gene ends (Nojima et al., 2015), and the chromatin template slows transcription there (Proudfoot, 2016). We assume a post-PAS rate of approximately 1.8 kb/min (\~30 bp/s); this value is a modelling assumption, not a measured rate, giving:

$${k}_{e2}=\frac{30 bp/s}{100 bp/node}=0.3{ s}^{-1}$$

The standard parameter ${k}_{e2}$ \= 30/100 \= 0.3 node s⁻¹ reflects an approximately 2-fold slowdown relative to the pre-PAS rate. Range: ${k}_{e2}$ ∈ \[0.1, 0.5\] node s⁻¹.

# **F. CTD Ser2 Phosphorylation Parameters (${k}_{Pon,min}$, ${k}_{Pon,slope}$, ${k}_{Poff}$)**

In human cells, Ser2P enrichment develops during early elongation and persists through the gene body (Vihervaara et al., 2023). PAS function and Pol II pausing further enrich Ser2P near the gene end (Davidson et al., 2014). These observations support progressive phosphorylation during transcription.

**Formation.** We use a simple position-dependent rate:

${k}_{Pon}(x)={k}_{Pon,min}+{k}_{Pon,slope}x,$

where $x$ is the node index measured from the TSS, with 100 bp per node. We set ${k}_{Pon,min}=0$ and ${k}_{Pon,slope}=5\times 1{0}^{-4}$ **s−1 per node**. The linear form and numerical slope are modeling assumptions; the experimental enrichment profiles provide qualitative constraints rather than a quantitative calibration.

**Removal.** We use ${k}_{Poff}=0.01$ **s−1** per unoccupied site as a working assumption, corresponding to an unprotected-site half-life of approximately 69 s. CTD-binding proteins can inhibit dephosphorylation (Palancade et al., 2004); we approximate this effect by protecting E-bound sites.

**Sensitivity ranges:** ${k}_{Pon,slope}\in [1{0}^{-4},\ 2\times 1{0}^{-3}]$ s−1 per node and ${k}_{Poff}\in [0.001,\ 0.02]$ s−1. The removal range tests unprotected-site half-lives from approximately 35 s to 11.6 min. 

## **G. Total Factor Pool Sizes (${E}_{total}$, ${Pol}_{total}$)**

### **Total E-factor Pool (${E}_{total}$)**

The E-factor pool represents all early CPA assembly factors that associate with the Pol II CTD, principally the CPSF complex. Because published proteomics datasets do not uniformly provide direct absolute copy-number measurements for CPA factors, we estimated copy numbers by calibrating relative protein abundance data from the PaxDB database against the absolute measurements reported by Zeiler et al. (2012).

**Calibration procedure.** We retrieved protein abundance values (in parts per million, ppm) for HeLa cells from PaxDB (https\://pax-db.org/) for each CPA factor of interest, with up to 5 measurements per protein, ranked by dataset quality score. We then paired these relative abundances with absolute copy numbers per HeLa cell from Zeiler et al. (2012), who used a SILAC-based Protein Epitope Signature Tag (PrEST) approach to quantify 40 proteins spanning \~6,000 to \~20,000,000 copies per cell (typical SEM \~20%). We log₁₀-transformed both quantities and fitted a linear regression to the paired data. This log-log calibration yielded a conversion equation of the form:

$${log}_{10}\left({copies/cell}\right)=m\times {log}_{10}\left({PaxDB ppm}\right)+b$$

exploiting the well-established near-proportionality between iBAQ-derived ppm values and absolute protein abundances (Pearson R ≈ 0.95 across large calibration sets). We log-transformed the mean of the top five PaxDB ppm values for each CPA factor and substituted it into this equation to estimate absolute copy number.

For reproducibility, report the fitted coefficients m and b, calibration sample size, fit R², and prediction uncertainty with the analysis.

**Results.** Table G1 summarizes the resulting copy number estimates for the core CPA complex subunits, along with the active RNA Pol II reference value.

**Table G1.** Estimated copy numbers per HeLa cell for early CPA complex subunits and active RNA Pol II.

| Protein | Est. Copies/Cell | Mean PaxDB (ppm) |
| :---- | :---- | :---- |
| CPSF-160 | \~102,000 | 30.2 |
| WDR33 | \~214,000 | 57.9 |
| CPSF-30 | \~292,000 | 78.9 |
| CPSF-100 | \~241,000 | 65.1 |
| CPSF-73 | \~187,000 | 50.7 |
| CstF-50 | \~332,000 | 89.8 |
| CstF-64 | \~371,000 | 100.2 |
| CstF-77 | \~268,000 | 72.5 |
| Active RNA Pol II | \~75,000 | From BioNumbers |

The CPSF complex functions as a preassembled unit, so the stoichiometry-limiting subunit CPSF-160 sets an upper bound on the number of intact CPSF complexes, at \~102,000 copies per cell. Because CPSF is strongly nuclear-enriched and functions co-transcriptionally (Boisvert et al., 2012; Itzhak et al., 2016), this per-cell estimate is treated as approximately equal to the per-nucleus pool. We therefore adopt:

$${E}_{total}\approx 100,000 molecules per nucleus$$

This is internally consistent with the active Pol II estimate of \~75,000, yielding an ${E}_{total}$/${Pol}_{total}$ ratio of ≈ 1.3, which implies that E-factor availability is not grossly in excess of Pol II at baseline — a prerequisite for the resource competition effects central to the model. Range: ${E}_{total}$ ∈ \[50,000, 200,000\] molecules, representing an approximately twofold sensitivity range around the limiting-subunit baseline and uncertainty in nuclear enrichment.

### **Total Pol II Pool (${Pol}_{total}$)**

According to BioNumbers (BNID 111166), HeLa cells contain about 75,000 active Pol II/III complexes per cell, measured as nascent RNA chains, of which about 10% are Pol III. We use 75,000 as the Pol II pool:

$$Po{l}_{total}\approx 75,000 molecules per nucleus$$

**Range: ${Pol}_{total}$** ∈ \[25,000, 140,000\] molecules, obtained from 50,000–200,000 molecules per cell and a 50–70% nuclear fraction.

### **Late CPA Factor Pool (${L}_{total}$)**

The late CPA factor pool, principally the CstF heterotrimer (CstF-50/CstF-64/CstF-77), is estimated using the same PaxDB–Zeiler calibration procedure described above. The three CstF subunits span a narrow range of \~268,000–371,000 copies per cell (Table G1). Because CstF functions as a preformed heterotrimer, the stoichiometry-limiting subunit CstF-77 gives \~268,000 intact complexes per cell. We round this to a working baseline of \~300,000 molecules per nucleus, within the uncertainty of the copy-number and nuclear-enrichment estimates:

$${L}_{total}\approx 300,000 molecules per nucleus$$

Range: ${L}_{total}$ ∈ \[200,000, 400,000\] molecules.

This pool sets ${[L]}_{free}$ and therefore ${k}_{c}$ (Section D). With ${K}_{on}$ \= 5 × 10⁵ M⁻¹s⁻¹ and ${[L]}_{free}$ ≈ 1.332 μM, the baseline ${K}_{d,L}$ \= 2 μM gives ${k}_{c}$ ≈ 0.14 s⁻¹ (mean commitment time about 7 s), and a tighter ${K}_{d,L}$ \= 0.5 μM gives ${k}_{c}$ ≈ 0.44 s⁻¹ (about 2.3 s). The tighter case agrees with the ≲3 s commitment that Chao et al. (1999) measured at a strong PAS. This checks the resulting ${k}_{c}$ against in vivo timing. It is not an independent confirmation of ${L}_{total}$.

## **References**

1. Boisvert, F.-M., Ahmad, Y., Gierliński, M., Charrière, F., Lamont, D., Scott, M., Barton, G., & Lamond, A. I. (2012). A quantitative spatial proteomics analysis of proteome turnover in human cells. *Molecular & Cellular Proteomics*, 11(3), M111.011429. https\://doi.org/10.1074/mcp.M111.011429

2. Caliskan, G., Hyeon, C., Perez-Salas, U., Briber, R. M., Woodson, S. A., & Thirumalai, D. (2005). Persistence length changes dramatically as RNA folds. *Physical Review Letters*, 95(26), 268303\. https\://doi.org/10.1103/PhysRevLett.95.268303

3. Chao, L. C., Jamil, A., Kim, S. J., Huang, L., & Martinson, H. G. (1999). Assembly of the cleavage and polyadenylation apparatus requires about 10 seconds in vivo and is faster for strong than for weak poly(A) sites. *Molecular and Cellular Biology*, 19(8), 5588–5600. https\://doi.org/10.1128/MCB.19.8.5588

4. Corden, J. L., Cadena, D. L., Ahearn, J. M., & Dahmus, M. E. (1985). A unique structure at the carboxyl terminus of the largest subunit of eukaryotic RNA polymerase II. *Proceedings of the National Academy of Sciences*, 82(23), 7934–7938. https\://doi.org/10.1073/pnas.82.23.7934

5. Darzacq, X., Shav-Tal, Y., de Turris, V., Brody, Y., Shenoy, S. M., Phair, R. D., & Singer, R. H. (2007). *In vivo* dynamics of RNA polymerase II transcription. *Nature Structural & Molecular Biology*, 14(9), 796–806. https\://doi.org/10.1038/nsmb1280

6. Davidson, L., Muniz, L., & West, S. (2014). 3′ end formation of pre-mRNA and phosphorylation of Ser2 on the RNA polymerase II CTD are reciprocally coupled in human cells. *Genes & Development*, 28(4), 342–356. https\://doi.org/10.1101/gad.231274.113

7. Dhusia, K., & Wu, Y. (2021). Classification of protein–protein association rates based on biophysical informatics. *BMC Bioinformatics*, 22(1), 408\. https\://doi.org/10.1186/s12859-021-04323-0

8. Fuchs, G., Voichek, Y., Benjamin, S., Gilad, S., Amit, I., & Oren, M. (2014). 4sUDRB-seq: measuring genomewide transcriptional elongation rates and initiation frequencies within cells. *Genome Biology*, 15(5), R69. https\://doi.org/10.1186/gb-2014-15-5-r69

9. Gilmartin, G. M., & Nevins, J. R. (1989). An ordered pathway of assembly of components required for polyadenylation site recognition and processing. *Genes & Development*, 3(12B), 2180–2190. https\://doi.org/10.1101/gad.3.12b.2180

10. Gilmartin, G. M., & Nevins, J. R. (1991). Molecular analyses of two poly(A) site-processing factors that determine the recognition and efficiency of cleavage of the pre-mRNA. *Molecular and Cellular Biology*, 11(5), 2432–2438. https\://doi.org/10.1128/mcb.11.5.2432-2438.1991

11. Grozdanov, P. N., Masoumzadeh, E., Latham, M. P., & MacDonald, C. C. (2018). The structural basis of CstF-77 modulation of cleavage and polyadenylation through stimulation of CstF-64 activity. *Nucleic Acids Research*, 46(22), 12022–12039. https\://doi.org/10.1093/nar/gky862

12. Hamilton, K., Sun, Y., & Tong, L. (2019). Biophysical characterizations of the recognition of the AAUAAA polyadenylation signal. *RNA*, 25(12), 1673–1680. https\://doi.org/10.1261/rna.070870.119

13. Hsin, J. P., & Manley, J. L. (2012). The RNA polymerase II CTD coordinates transcription and RNA processing. *Genes & Development*, 26(19), 2119–2137. https\://doi.org/10.1101/gad.200303.112

14. Huang, L., Chu, H.-F., & Tong, L. (2025). Molecular basis for the recognition of low-frequency polyadenylation signals by mPSF. *Nucleic Acids Research*, 53(17), gkaf890. https\://doi.org/10.1093/nar/gkaf890

15. Hyeon, C., Dima, R. I., & Thirumalai, D. (2006). Size, shape, and flexibility of RNA structures. *The Journal of Chemical Physics*, 125(19), 194905\. https\://doi.org/10.1063/1.2364190

16. Itzhak, D. N., Tyanova, S., Cox, J., & Borner, G. H. H. (2016). Global, quantitative and dynamic mapping of protein subcellular localization. *eLife*, 5, e16950. https\://doi.org/10.7554/eLife.16950

17. Mei, K., Jin, Z., Ren, F., Wang, Y., Chang, Z., & Wang, X. (2014). Structural basis for the recognition of RNA polymerase II C-terminal domain by CREPT and p15RS. *Science China Life Sciences*, 57(1), 97–106. https\://doi.org/10.1007/s11427-013-4589-7

18. Meinhart, A., & Cramer, P. (2004). Recognition of RNA polymerase II carboxy-terminal domain by 3′\-RNA-processing factors. *Nature*, 430(6996), 223–226. https\://doi.org/10.1038/nature02679

19. Milo, R., Jorgensen, P., Moran, U., Weber, G., & Springer, M. (2010). BioNumbers—the database of key numbers in molecular and cell biology. *Nucleic Acids Research*, 38(Database issue), D750–D753. https\://doi.org/10.1093/nar/gkp889

20. Ni, Z., Xu, C., Guo, X., Hunter, G. O., Kuznetsova, O. V., Tempel, W., … & Greenblatt, J. F. (2014). RPRD1A and RPRD1B are human RNA polymerase II C-terminal domain scaffolds for Ser5 dephosphorylation. *Nature Structural & Molecular Biology*, 21(8), 686–695. https\://doi.org/10.1038/nsmb.2853

21. Noble, C. G., Hollingworth, D., Martin, S. R., Ennis-Adeniran, V., Smerdon, S. J., Kelly, G., Taylor, I. A., & Ramos, A. (2005). Key features of the interaction between Pcf11 CID and RNA polymerase II CTD. *Nature Structural & Molecular Biology*, 12(2), 144–151. https\://doi.org/10.1038/nsmb887

22. Nojima, T., Gomes, T., Grosso, A. R. F., Kimura, H., Dye, M. J., Dhir, S., Carmo-Fonseca, M., & Proudfoot, N. J. (2015). Mammalian NET-seq reveals genome-wide nascent transcription coupled to RNA processing. *Cell*, 161(3), 526–540. https\://doi.org/10.1016/j.cell.2015.03.027

23. Northrup, S. H., & Erickson, H. P. (1992). Kinetics of protein-protein association explained by Brownian dynamics computer simulation. *Proceedings of the National Academy of Sciences*, 89(8), 3338–3342. https\://doi.org/10.1073/pnas.89.8.3338

24. Palancade, B., Marshall, N. F., Tremeau-Bravard, A., Bensaude, O., Dahmus, M. E., & Dubois, M. F. (2004). Dephosphorylation of RNA polymerase II by CTD-phosphatase FCP1 is inhibited by phospho-CTD associating proteins. *Journal of Molecular Biology*, 335(2), 415–424. https\://doi.org/10.1016/j.jmb.2003.10.036

25. Phillips, R., Kondev, J., Theriot, J., & Garcia, H. G. (2013). *Physical biology of the cell* (2nd ed.). Garland Science.

26. Portz, B., Lu, F., Gibbs, E. B., Mayfield, J. E., Rachel Mehaffey, M., Zhang, Y. J., Brodbelt, J. S., Showalter, S. A., & Gilmour, D. S. (2017). Structural heterogeneity in the intrinsically disordered RNA polymerase II C-terminal domain. *Nature Communications*, 8, 15231\. https\://doi.org/10.1038/ncomms15231

27. Proudfoot, N. J. (2016). Transcriptional termination in mammals: stopping the RNA polymerase II juggernaut. *Science*, 352(6291), aad9926. https\://doi.org/10.1126/science.aad9926

28. Riback, J. A., Bowman, M. A., Zmyslowski, A. M., Plaxco, K. W., Clark, P. L., & Sosnick, T. R. (2019). Commonly used FRET fluorophores promote collapse of an otherwise disordered protein. *Proceedings of the National Academy of Sciences*, 116(18), 8889–8894. https\://doi.org/10.1073/pnas.1813038116

29. Rosonina, E., & Blencowe, B. J. (2004). Analysis of the requirement for RNA polymerase II CTD heptapeptide repeats in pre-mRNA splicing and 3′\-end cleavage. *RNA*, 10(4), 581–589. https\://doi.org/10.1261/rna.5207204

30. Schreiber, G., Haran, G., & Zhou, H.-X. (2009). Fundamental aspects of protein–protein association kinetics. *Chemical Reviews*, 109(3), 839–860. https\://doi.org/10.1021/cr800373w

31. Sdano, M. A., Fulcher, J. M., Palani, S., Chandrasekharan, M. B., Parnell, T. J., Whitby, F. G., Formosa, T., & Hill, C. P. (2017). A novel SH2 recognition mechanism recruits Spt6 to the doubly phosphorylated RNA polymerase II linker at sites of transcription. *eLife*, 6, e28723. https\://doi.org/10.7554/eLife.28723

32. Vihervaara, A., Versluis, P., Himanen, S. V., & Lis, J. T. (2023). PRO-IP-seq tracks molecular modifications of engaged Pol II complexes at nucleotide resolution. *Nature Communications*, 14, 7039\. https\://doi.org/10.1038/s41467-023-42715-3

33. Vos, S. M., Farnung, L., Boehning, M., Wigge, C., Linden, A., Urlaub, H., & Cramer, P. (2018). Structure of activated transcription complex Pol II–DSIF–PAF–SPT6. *Nature*, 560(7720), 607–612. https\://doi.org/10.1038/s41586-018-0440-4

34. Wang, M., Weiss, M., Simonovic, M., Haertinger, G., Schrimpf, S. P., Hengartner, M. O., & von Mering, C. (2012). PaxDb, a database of protein abundance averages across all three domains of life. *Molecular & Cellular Proteomics*, 11(8), 492–500. https\://doi.org/10.1074/mcp.O111.014704

35. Yang, Q., Coseno, M., Gilmartin, G. M., & Doublié, S. (2011). Crystal structure of a human cleavage factor CFIm25/CFIm68/RNA complex provides an insight into poly(A) site recognition and RNA looping. *Structure*, 19(3), 368–377. https\://doi.org/10.1016/j.str.2010.12.021

36. Yang, W., Hsu, P. L., Yang, F., Song, J.-E., & Varani, G. (2018). Reconstitution of the CstF complex unveils a regulatory role for CstF-50 in recognition of 3′\-end processing signals. *Nucleic Acids Research*, 46(2), 493–503. https\://doi.org/10.1093/nar/gkx1177

37. Zeiler, M., Straube, W. L., Lundberg, E., Uhlen, M., & Mann, M. (2012). A Protein Epitope Signature Tag (PrEST) library allows SILAC-based absolute quantification and multiplexed determination of protein copy numbers in cell lines. *Molecular & Cellular Proteomics*, 11(3), O111.009613. https\://doi.org/10.1074/mcp.O111.009613

38. Zhou, H.-X. (1993). Brownian dynamics study of the influences of electrostatic interaction and diffusion on protein-protein association kinetics. *Biophysical Journal*, 64(6), 1711–1726. https\://doi.org/10.1016/S0006-3495(93)81543-1
