# **Appendix: Parameter Estimation**

## **A. Molecular Association Rate ${K}_{on}$**

Throughout this estimation, we use a basal molecular association rate of ${K}_{on}$ \= $5\ ×\ {10}^{5}$ M⁻¹ s⁻¹ for bimolecular interactions. This value is derived from the Brownian diffusion-limited encounter rate in aqueous solution, which typically yields association rates on the order of \~$5\ ×\ {10}^{5}$ M⁻¹s⁻¹ for protein–protein interactions (Dhusia & Wu, 2021; Fu et al., 2026; Northrup & Erickson, 1992; Schreiber et al., 2009; Zhou, 1993).

Importantly, the exact value of ${K}_{on}$ does not affect the dimensionless equilibrium ratios that determine simulation behavior. For the hexamer-binding module, ${k}_{Hd}$ \= ${k}_{Hoff}$ / ${k}_{Hon}$, and for the E-factor CTD-binding module, ${k}_{E}$ \= ${k}_{Eoff}$ / ${k}_{Eon}$. Because the corresponding forward and reverse microscopic rates are scaled consistently from the same ${K}_{on}$ convention, ${K}_{on}$ cancels from these ratios. Thus, ${K}_{on}$ is retained for physical interpretability, unit consistency, and estimating absolute timescales, while the relevant simulation behavior is controlled primarily by ${k}_{Hd}$, ${k}_{Ed}$, and the downstream commitment rate ${k}_{c}$.

## **B. First-Principles Estimation of the Hexamer Binding Rate**

### **E-factor PAS Recognition Complex — Hexamer Binding (${k}_{Hon}$)**

The formation rate of the committed cleavage and polyadenylation (CPA) complex (${k}_{Hon}$) represents the tethered encounter between early CPA factors (E factors) pre-loaded on the Pol II CTD and the newly transcribed poly(A) signal hexamer. The relevant geometry includes two flexible tethers whose proximal ends attach to the Pol II body: the CTD arm and the nascent-RNA arm.

Because the CTD can be represented as an intrinsically disordered polymer, the worm-like chain (WLC) model approximates this problem well (Caliskan et al., 2005; Hyeon et al., 2006; Phillips et al., 2013).

**Physical parameterization.** We parameterize the model with a representative baseline and a physically motivated range:

* **Effective Contour Length (${L}_{eff}$):** Baseline 107 nm. The human Pol II CTD repeat array contains 52 heptad repeats, corresponding nominally to 364 residues and a fully extended contour length of approximately 131–138 nm (Corden et al., 1985; Rosonina & Blencowe, 2004; Hsin & Manley, 2012; Choi et al., 2011). An approximately 105-residue linker connects the repeat array to the Pol II core, adding roughly 38–40 nm when treated as a flexible tether (Vos et al., 2018). Assuming that E factors bind uniformly along the repeat array, their mean contour coordinate is therefore approximately

  ${L}_{eff}\ ≈\ 38\ +\ \frac{138}{2}\ ≈\ 107\ nm$

  We consequently use ${L}_{eff}\ =\ 107$ nm as the representative effective contour length

* **Adjusted Persistence Length (${l}_{p}$):** Baseline 2.0 nm. Range: \[1.0, 3.0\] nm. A naked, highly flexible IDP may have a short persistence length (\< 1.0 nm), whereas the phosphorylated CTD is expected to be more extended and stiffer due to its high negative charge (Portz et al., 2017).

* **Anchor Separation (d):** The distance between two anchors is estimated at d \= 10 nm, the distance from the CTD linker to the exit channel of the RNA .

**Effective concentration.** For independent, long, flexible chains, we approximate the two endpoint distributions as Gaussians. Their convolution gives a relative-coordinate Gaussian with mean-square span ${⟨{R}^{2}⟩}_{tot}$ \= ${⟨{R}^{2}⟩}_{CTD}$ \+ ${⟨{R}^{2}⟩}_{RNA}$. Using the 3D Gaussian probability density for a WLC end-to-end vector, the effective local concentration, J-density, at distance d from the anchor (Darabi et al., 2022; Phillips et al., 2013\) is:

${J}_{density}={\left({\frac{3}{2π⟨{R}^{2}{⟩}_{tot}}}\right)}^{3/2}exp​\left({\frac{-3{d}^{2}}{2⟨{R}^{2}{⟩}_{tot}}}\right)$

**CTD arm:** For the ideal WLC, ${⟨{R}^{2}⟩}_{CTD}$ ≈ 2 × ${L}_{eff}$ × ${l}_{p}$ \= 428 nm², corresponding to R ≈ 20.69 nm. This is of the same order as the SAXS measurements of the phosphorylated Drosophila CTD (Portz et al., 2017). In this report, the radius of gyration, ${R}_{g}$ \= 6.06 nm, is for the 42-repeat fly CTD. ${R}_{g}$ is the root-mean-square distance of all chain segments from the chain’s center of mass, and R is the end-to-end distance. Using the relationship R ≈ 2.5 ${R}_{g}$, which holds approximately for unfolded or disordered protein ensembles (Riback et al., 2019), the phosphorylated fly CTD gives R ≈ 15.2 nm. Scaling this value from the 42-repeat fly CTD to the 52-repeat human CTD as √(52/42) yields:

$R\ ≈\ 16.9\ nm\ (SAXS-scaled);\ baseline\ WLC\ estimate:\ 20.7\ nm$

**RNA arm:** For 200–600 nt, the contour-length range is ${L}_{RNA}$ \= 200 × 0.61 to 600 × 0.61 \= 122–366 nm. The baseline 400-nt segment has ${L}_{RNA}$ \= 244 nm. With ${l}_{p,RNA}$ \= 1.26 nm, ${⟨{R}^{2}⟩}_{RNA}$ ≈ 2 ${l}_{p}$ L \= 614.88 nm² at baseline.

For our baseline ${{R}^{2}}_{tot}$ \= 428 \+ 614.88 \= 1042.88 nm² and d \= 10 nm, this yields ${J}_{density}$ ≈ 8.484 × 10⁻⁶ nm⁻³. Using 1 nm⁻³ ≈ 1.66 M (Jeong et al., 2016):

${J}_{density}\ ≈\ 1.409\ ×\ 10⁻⁵\ M\ =\ 14.09\ μM$

Applying our baseline association rate (${K}_{on}$ \= 5 × 10⁵ M⁻¹s⁻¹):

${k}_{Hon}\ =\ {K}_{on}\ ·\ {J}_{eff}$

${k}_{Hon}\ ≈\ 5\ ×\ {10}^{5}\ M⁻¹\ s⁻¹\ ×\ 1.409\ ×\ {10}^{-5}\ M\ =\ 7\ s⁻¹$

**Range evaluation:**

${k}_{Hon}\ ∈\ [4,\ 17]\ s⁻¹$

### **E-factor PAS Recognition Complex — Hexamer Dissociation (${k}_{Hoff}$)**

The equilibrium dissociation constant ${K}_{d}$ for recognition of the AAUAAA polyadenylation signal (PAS) by mPSF (CPSF-160, WDR33, CPSF-30 and Fip1) depends on the surrounding RNA. In competition assays with unlabelled RNA it is about 3–40 nM for AAUAAA with flanking sequence (Hamilton et al., 2019; Huang, Chu & Tong, 2025) and about 2000 nM for the bare hexamer (Hamilton et al., 2019). Chao et al. (1999) found that strong (SV40 late) and weak (SV40 early) PASs with the same AAUAAA differ several-fold in assembly rate, so we treat ${k}_{Hoff}$ as an effective rate for the whole PAS region: a strong PAS lies in the 3–40 nM range and a weak PAS near 2000 nM. Using our nuclear association rate to derive the microscopic dissociation rate:

${k}_{Hoff}\ =\ {K}_{on}\ ·\ {K}_{d}\ =\ 5\ ×\ {10}^{5}\ M⁻¹\ s⁻¹\ ×\ 20\ ×\ {10}^{-9}\ M\ =\ 0.01\ s⁻¹$ (strong PAS, middle of 3–40 nM)

${k}_{Hoff}\ =\ {K}_{on}\ ·\ {K}_{d}\ =\ 5\ ×\ {10}^{5}\ M⁻¹\ s⁻¹\ ×\ 2000\ ×\ {10}^{-9}\ M\ =\ 1\ s⁻¹$ (weak PAS)

At the weak PAS the recognized complex lasts about 1 s, shorter than the mean commitment time 1/${k}_{c}$ ≈ 7 s, so loss of H competes directly with commitment. This fits the view of Chao et al. that PAS strength is kinetic. Variant hexamers bind more weakly than AAUAAA (79 nM–13 μM; Huang et al., 2025) and are enriched at proximal PASs (Tian et al., 2005; Beaudoing et al., 2000). They are not modelled separately.

**Range evaluation:**

${k}_{Hoff}\ ∈\ [0.0015,\ 1]\ s⁻¹$

## **C. E-factor Binding Parameters**

### **Early Factor Association Rate (${k}_{Eon}$)**

E-factor binding to the Pol II CTD is assumed to be diffusion-limited. We convert the molar rate constant into a per-molecule rate suitable for physiologically scaled simulations using:

${k}_{Eon}=\frac{{K}_{on}}{{N}_{A}\cdot {V}_{nucleus}}$

${k}_{Eon}\ =\ \frac{5\ ×\ {10}^{5}\ M⁻¹\ s⁻¹}{6.022\ ×\ {10}^{23}\ ×\ 3.74\ ×\ {10}^{-13}\ L}\ ≈\ 2.5\ ×\ {10}^{-6}\ molecule⁻¹\ s⁻¹$

where ${N}_{A}$ \= 6.022 × 10²³ molecules/mol and ${V}_{nucleus}$ \= 3.74 × 10⁻¹³ L corresponds to a representative 374 μm³ mammalian nucleus (BioNumbers). Range evaluation. Over the biological range of nuclear volumes \[200, 500\] μm³, ${k}_{Eon}$ ∈ \[1.66 × 10⁻⁶, 4.15 × 10⁻⁶\] molecule⁻¹s⁻¹.

### **Early Factor Dissociation Rate (${k}_{Eoff}$)**

The intrinsic binding affinity ${K}_{d,E}$ of an early CPA factor for a single phosphorylated Ser2-Pro3 CTD heptad repeat is measured in vitro at 1–10 μM (Lunde et al., 2010; Mei et al., 2014; Meinhart & Cramer, 2004; Moreno et al., 2025; Ni et al., 2014; Noble et al., 2005). Using ${k}_{off}$ \= ${K}_{on}$ · ${K}_{d,E}$ for a single repeat:

${k}_{Eoff,single}\ =\ (1–10)\ ×\ {10}^{-6}\ M\ ×\ 5\ ×\ {10}^{5}\ M⁻¹\ s⁻¹\ =\ 0.5–5\ s⁻¹$

**Multivalent correction.** The E factor engages 2–4 tandem CTD repeats, as suggested by co-crystal structures of CID–CTD interactions (Meinhart & Cramer, 2004; Noble et al., 2005). Multivalent or bivalent engagement of a disordered CTD polymer substantially increases the effective binding affinity relative to the single-repeat value, which can push the effective ${k}_{Eoff}$ to the lower end of the range.

${k}_{Eoff}^{eff}\approx 0.5{\ s}^{-1}$

**Supporting data.**
- **Single CIDs** on a doubly phosphorylated Ser2P diheptad bind with ${K}_{d}$ ≈ 6–160 μM. Examples: RPRD1A/1B/2 6.8–8.4 μM (Ni et al., 2014); Rtt103 6–15 μM; Pcf11 130–160 μM (Lunde et al., 2010; Jasnovidova et al., 2017).
- **Multivalent constructs** bind more tightly. Rtt103 on a four-repeat peptide has a first-site ${K}_{d}$ of 3 μM (Lunde et al., 2010), and dimeric full-length RPRD1B binds a two-site peptide at 1.7 μM (Ni et al., 2014).
- **In vivo residence.** The only direct kinetic measurement is single-molecule tracking in yeast: Pcf11, Cft1, Rtt103 and Ssu72 stay bound for 2.1–2.9 s (Ling et al., 2026). This gives ${k}_{off}$ ≈ 0.3–0.5 s⁻¹ for the whole complex.

Range evaluation: ${k}_{Eoff}$ in \[0.1, 5\] s⁻¹.
- The lower end allows tighter multivalent retention than the in vivo residence times. Rtt103 on the full phosphorylated CTD is reported at 0.3 μM (Jasnovidova et al., 2017), which gives 0.15 s⁻¹.
- The upper end is a single CID on one diheptad (${K}_{d}$ ≈ 10 μM).
- The working value 0.5 s⁻¹ lies between the multivalent in vitro values and the in vivo residence times.

### **Engaged-E Detachment Rate (${k}_{Eoff,engaged}$)**

In a recognized complex, one E bridges the CTD and the hexamer-bound H on the nascent RNA. When this engaged E releases the CTD, H keeps it near the RNA, where it can rebind a free Ser2P site. It therefore returns to the pool more slowly than an ordinary E (${k}_{Eoff}$ = 0.5 s⁻¹). How firmly H holds the released E has not been measured. If E stayed on H until H left the hexamer, the rate would be far below ${k}_{Eoff}$. If E lost H at once, it would equal ${k}_{Eoff}$. We take a conservative value between these limits:

${k}_{Eoff,engaged}$ **= 0.1 s⁻¹** (one-fifth of ${k}_{Eoff}$)

**Range:** ${k}_{Eoff,engaged}$ ∈ \[0.01, 0.5\] s⁻¹. The upper bound is ${k}_{Eoff}$ (no hold by H). The lower bound allows a firm hold by H.

## **D. CPA Commitment Rate ${k}_{c}$**

The CPA commitment rate ${k}_{c}$ represents the productive, effectively irreversible assembly of the complete CPA machinery after early PAS recognition. This step includes diffusional recruitment of late CPA factors, such as CstF; productive engagement with the downstream sequence element (DSE); correct steric/orientational alignment; and assembly of the remaining components needed to form a complete cleavage-competent complex. Therefore, ${k}_{c}$ should be interpreted as an effective commitment rate rather than a raw diffusion-limited collision rate.

**Late CPA factor copy number.** We estimate the global nuclear copy number of late-processing factors to be approximately 300,000 molecules per nucleus (see Section G). Using our baseline nuclear volume of 3.74 × 10⁻¹³ L:

$\left[{{L}_{free}}\right]=\frac{300,000}{6.022\times {10}^{23}\times 3.74\times {10}^{-13}}\approx 1.33\ μM$

Let R \= \[REH\], X \= \[REHL\], and Y \= \[REHL₂\]. Define r \= ${K}_{on}$${[L]}_{free}$ and q \= ${k}_{off}$ \= ${K}_{d,L}.{K}_{on}$ . We assume all three forward steps have rate r, and the first two steps reverse sequentially with rate q:

$REH\ \ {{⇌}^{r}}_{q}\ \ REHL\ \ {{⇌}^{r}}_{q}\ \ REHL₂\ \ {\rightarrow }^{r}\ \ commit$

The fraction of initial capture events that reaches commitment is:

${f}_{lock}\ =\ \frac{{J}_{commit}}{{J}_{capture}}\ =\ \frac{{r}^{2}}{{q}^{2}\ +\ qr\ +\ {r}^{2}}\ =\ \frac{{[L]}^{2}}{{{K}_{d,L}}^{2}\ +\ {K}_{d,L}[L]\ +\ {[L]}^{2}}$

The corresponding effective commitment rate is:

${k}_{c}\ =\ r{f}_{lock}\ =\ K\_on[L]\ \frac{{[L]}^{2}}{{{K}_{d,L}}^{2}\ +\ {K}_{d,L}[L]\ +\ {[L]}^{2}}$

Using ${[L]}_{free}$ \= 1.332 μM and ${K}_{d,L}$ \= 2 μM gives ${k}_{c}$ \= 0.14 s⁻¹.

Range evaluation. Varying the nuclear volume over \[300, 500\] μm³ and the late-factor copy number over \[200,000, 400,000\] gives ${[L]}_{free}$ ≈ 0.664–2.214 μM. Using ${K}_{d,L}$ ∈ \[0.1, 4\] μM gives ${k}_{c}$ ≈ 0.01–1 s⁻¹. We therefore use ${k}_{c}$ ≈ 0.14 s⁻¹ as the baseline effective commitment rate (Yang et al., 2018; Monarez et al., 2007).

**Measured late-factor affinities.** The ${K}_{d,L}$ interval is supported by:
- **Tight end:** the CstF hexamer on GU-rich RNA, 0.12–0.9 μM, and CstF-77/64 on GU14, 0.16 μM (Yang et al., 2018); CFIm25/68 on two UGUA motifs, 0.22 μM (Yang et al., 2011).
- **Weak end:** the isolated CstF-64 RRM, 0.7–2.5 μM, and CstF-77/64 on a U-rich DSE, 3.7 μM (Yang et al., 2018; Grozdanov et al., 2018).

All of these are in vitro values without competitor RNA, so the baseline ${K}_{d,L}$ = 2 μM is an effective in-nucleus value, not a measured one.

**In vivo timing.** Every measurement runs from PAS transcription to cleavage. That interval also includes recognition, so its rate is a lower bound on ${k}_{c}$.
- Chao et al. (1999) found commitment at the weak SV40 early PAS takes 10–20 s, about 0.05–0.1 s⁻¹ overall. The strong SV40 late PAS commits several times faster.
- Genome-wide cleavage half-lives in *Drosophila* S2 cells have medians of 36–42 s, about 0.016–0.019 s⁻¹ overall, and many sites are too fast to resolve (Torres-Ulloa et al., 2024).

**Range:** ${k}_{c}$ ∈ \[0.015, 0.6\] s⁻¹.
- The lower end is the slowest well-documented overall cleavage rate, the genome-wide median above, which bounds ${k}_{c}$ from below.
- The upper end is the scheme above with the tightest holo-CstF affinity (${K}_{d,L}$ = 0.12 μM) at baseline ${[L]}_{free}$. It gives ${k}_{c}$ ≈ 0.6 s⁻¹ (mean commitment time about 2 s), consistent with a strong PAS. The scheme cannot exceed ${K}_{on}{[L]}_{free}$ = 0.67 s⁻¹ at baseline.

## **E. Pol II Elongation Rates (${k}_{e}$, ${k}_{e2}$)**

### **Pre-PAS Elongation Rate (${k}_{e}$)**

Pol II elongation velocity in human cells has been measured by multiple approaches: metabolic labeling, live-cell imaging, and nascent RNA sequencing, and consistently yields values of 3–5 kb/min (Darzacq et al., 2007; Jonkers et al., 2014; Fuchs et al., 2014; Veloso et al., 2014; Krebs et al., 2017). We use a central estimate of 4 kb/min, corresponding to 4,000 bp / 60 s ≈ 66.67 bp/s.

In the model, the gene is discretized into nodes of ${L}_{a}$ \= 100 bp, so the elongation rate in model units is:

${k}_{e}\ =\ \frac{66.67\ bp/s}{100\ bp/node}\ =\ 0.65\ node/s$

**Range: ${k}_{e}$** ∈ \[0.50, 0.83\] node s⁻¹, corresponding to elongation velocities of 3–5 kb/min.

### **Post-PAS Elongation Rate (${k}_{e2}$)**

After PAS recognition, Pol II slows down, consistent with elongation-rate measurements showing reduced processivity downstream of cleavage sites (Fong et al., 2015; Jonkers et al., 2014). We estimate a post-PAS rate of approximately 1.8 kb/min (\~30 bp/s), giving:

${k}_{e2}=\frac{30\ bp/s}{100\ bp/node}=0.3{\ s}^{-1}$

The standard parameter ${k}_{e2}$ \= 30/100 \= 0.3 node s⁻¹ reflects approximately a 2-fold slowdown relative to the pre-PAS rate, consistent with the pausing and reduced elongation observed near poly(A) sites (Nojima et al., 2015; Proudfoot, 2016). Range: ${k}_{e2}$ ∈ \[0.1, 0.5\] node s⁻¹.

## **F. CTD Ser2 Phosphorylation Parameters (${k}_{Pon,min}$, ${k}_{Pon,slope}$, ${k}_{Poff}$)**

**What p counts.** One E-binding site is two adjacent phosphorylated Ser2P repeats. CID readers bind this unit most tightly (Lunde et al., 2010; Ni et al., 2014), and its 1–10 μM affinity is the one used for ${k}_{Eoff}$ in Section C. With $n$ marks scattered over the 52 repeats, the expected number of adjacent pairs is $n(n-1)/52$. Each phosphorylation event adds one phosphorylated Ser2P repeats, and E-bound sites are protected.

**Site formation (${k}_{Pon,min}$, ${k}_{Pon,slope}$).** Ser2 kinases act distributively (Baumli et al., 2012; Czudnochowski et al., 2012). If marks are deposited at random at a constant rate ${k}_{dep}$ per CTD, a polymerase that has been on the gene for time $t$ carries a Poisson number of marks with mean ${k}_{dep}t$, and on average $({k}_{dep}t)^{2}/52$ adjacent pairs. Pairs therefore form at rate $2{k}_{dep}^{2}t/52$, which rises linearly with distance from zero at the TSS:

${k}_{Pon}(x) = {k}_{Pon,min} + {k}_{Pon,slope}\cdot x, \quad {k}_{Pon,slope} = \frac{2{k}_{dep}^{2}}{52\,{k}_{e}}, \quad {k}_{Pon,min} = 0$

Here $x$ is the node index from the TSS. ${k}_{Pon}$ is the rate at which E sites (pairs) form, not the per-repeat phosphorylation rate ${k}_{dep}$. It is zero at the TSS because a site needs two marks.

Ser2P antibody signals plateau about 1.5 kb from the TSS (Vihervaara et al., 2023) and do not scale with gene length (Kim et al., 2010). We read the plateau as about 90% of Pol II carrying at least one mark ($n$ ≈ 2.3), reached after 23 s at 65 bp/s. This gives ${k}_{dep}$ ≈ 0.1 s⁻¹, well within kinase capacity: CDK9 turns over at 2.97 s⁻¹ in vitro (Baumli et al., 2012), and HEK293T cells contain roughly 0.6–0.8 CDK9 or CDK12 per Pol II (Cho et al., 2022). We therefore use ${k}_{Pon,slope}$ **= 5 × 10⁻⁴ s⁻¹ per node** and ${k}_{Pon,min}$ **= 0**. The Ser2P peak downstream of the PAS depends on PAS function (Davidson et al., 2014) and is not used in this calibration.

**Range:** ${k}_{Pon,slope}$ ∈ \[1 × 10⁻⁴, 2 × 10⁻³\] s⁻¹ per node, so that the mean of 2.3 marks is reached between about 3.6 and 0.8 kb from the TSS. The lower bound allows for slower elongation over the first kilobases (Jonkers et al., 2014).

**Site loss (${k}_{Poff}$).** ${k}_{Poff}$ is the loss rate of each unoccupied site. A pair is lost when either of its marks is removed, so ${k}_{Poff}$ ≈ 2${k}_{r}$, where ${k}_{r}$ is the removal rate per mark.
- **CDK12 inhibition.** Ser2P on chromatin-bound Pol II falls markedly within 15 min and keeps falling at 30 min (Tellier et al., 2020). This includes the 3′ ends of long genes, where Pol II is not replaced within that time.
- **Active removal.** The loss requires PP2A recruited by Integrator (Vervoort et al., 2021; Tellier et al., 2022), so it reflects dephosphorylation of engaged Pol II and not only run-off.
- **Rate.** These time courses give ${k}_{r}$ ≈ 8 × 10⁻⁴–3 × 10⁻³ s⁻¹ (half-life about 4–15 min). This is consistent with:
  - the slower Ser2P loss when CDK9 alone is inhibited (Laitem et al., 2015; Parua et al., 2020);
  - Ser2P depending on how long Pol II spends on the gene (Fong et al., 2017);
  - FCP1 capacity (Schwer et al., 2015; Cho et al., 2022).

The inhibitor time courses give ${k}_{Poff}$ ≈ 0.0016–0.006 s⁻¹. The two ends of this interval are not equally firm.
- **Lower end (firm).** At least half the signal is lost by 15 min. Any residual kinase activity only slows the observed decay, so this is a lower bound on ${k}_{r}$.
- **Upper end (soft).** The upper end rests on Ser2P still falling between 15 and 30 min. Other Ser2 kinases (CDK9, CDK13) stay active when CDK12 alone is inhibited, and kinase shut-off is not instantaneous. The observed decay is therefore slower than the true removal rate, so it does not cap ${k}_{r}$ from above.

We use ${k}_{Poff}$ **= 0.01 s⁻¹** (${k}_{r}$ = 0.005 s⁻¹, mark half-life ≈ 2.3 min). Removal is still slow compared with transit to a 20 kb PAS (about 5 min), as the linear ${k}_{Pon}(x)$ assumes.

**Range:** ${k}_{Poff}$ ∈ \[0.0016, 0.02\] s⁻¹, from the inhibitor-based lower bound to twice the working value.


## **G. Total Factor Pool Sizes (${E}_{total}$, ${Pol}_{total}$)**

### **Total E-factor Pool (${E}_{total}$)**

The E-factor pool represents all early CPA assembly factors that associate with the Pol II CTD, principally the CPSF complex. Because published proteomics datasets do not uniformly provide direct absolute copy-number measurements for CPA factors, we estimated copy numbers by calibrating relative protein abundance data from the PaxDB database against the absolute measurements reported by Zeiler et al. (2012).

**Calibration procedure.** We retrieved protein abundance values (in parts per million, ppm) for HeLa cells from PaxDB (https\://pax-db.org/) for each CPA factor of interest, with up to 5 measurements per protein, ranked by dataset quality score. We then paired these relative abundances with absolute copy numbers per HeLa cell from Zeiler et al. (2012), who used a SILAC-based Protein Epitope Signature Tag (PrEST) approach to quantify 40 proteins spanning \~6,000 to \~23,000,000 copies per cell (typical SEM \~20%). We log₁₀-transformed both quantities and fitted a linear regression to the paired data. This log-log calibration yielded a conversion equation of the form:

${log}_{10}\left({copies/cell}\right)=m\times {log}_{10}\left({PaxDB\ ppm}\right)+b$

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
| Active RNA Pol  II | \~70,000 | From BioNumbers |

The CPSF complex functions as a preassembled unit, so the stoichiometry-limiting subunit CPSF-160 sets an upper bound on the number of intact CPSF complexes, at \~102,000 copies per cell. Because CPSF is strongly nuclear-enriched and functions co-transcriptionally (Boisvert et al., 2012; Itzhak et al., 2016\), this per-cell estimate is treated as approximately equal to the per-nucleus pool. We therefore adopt:

${E}_{total}\approx 100,000\ molecules\ per\ nucleus$

This is internally consistent with the active Pol II estimate of \~70,000, yielding an ${E}_{total}$/${Pol}_{total}$ ratio of ≈ 1.4, which implies that E-factor availability is not grossly in excess of Pol II at baseline — a prerequisite for the resource competition effects central to the model. Range: ${E}_{total}$ ∈ \[50,000, 200,000\] molecules, representing an approximately twofold sensitivity range around the limiting-subunit baseline and uncertainty in nuclear enrichment.

### **Total Pol II Pool (${Pol}_{total}$)**

Total Pol II levels in human cells have been estimated by quantitative western blotting and mass spectrometry to be 50,000–200,000 molecules per HeLa cell, with the largest calibrated estimates placing RPB1 at \~150,000 molecules per cell (Bhatt et al., 2012; Liao & Bhatt, 2019; Mayer et al., 2010). A typical nuclear fraction of \~50–70% gives an effective nuclear Pol II pool of:

$Po{l}_{total}\approx 70,000\ molecules\ per\ nucleus$

**Range: ${Pol}_{total}$** ∈ \[25,000, 140,000\] molecules, obtained from 50,000–200,000 molecules per cell and a 50–70% nuclear fraction.

### **Late CPA Factor Pool (${L}_{total}$)**

The late CPA factor pool, principally the CstF heterotrimer (CstF-50/CstF-64/CstF-77), is estimated using the same PaxDB–Zeiler calibration procedure described above. The three CstF subunits span a narrow range of \~268,000–371,000 copies per cell (Table G1). Because CstF functions as a preformed heterotrimer, the stoichiometry-limiting subunit CstF-77 gives \~268,000 intact complexes per cell. We round this to a working baseline of \~300,000 molecules per nucleus, within the uncertainty of the copy-number and nuclear-enrichment estimates:

${L}_{total}\approx 300,000\ molecules\ per\ nucleus$

A complementary kinetic cross-check is consistent with this estimate: with ${K}_{on}$ \= 5 × 10⁵ M⁻¹s⁻¹, ${[L]}_{free}$ ≈ 1.332 μM, and ${K}_{d,L}$ \= 0.5 μM, the sequential mechanism gives ${k}_{c}$ ≈ 0.439 s⁻¹. Range: ${L}_{total}$ ∈ \[200,000, 400,000\] molecules.

## **References**

1. Baumli, S., Hole, A. J., Wang, L.-Z., Noble, M. E. M., & Endicott, J. A. (2012). The CDK9 tail determines the reaction pathway of positive transcription elongation factor b. *Structure*, 20(10), 1788–1795. https\://doi.org/10.1016/j.str.2012.08.011

2. Bhatt et al. (2012). \[UNVERIFIED — author to confirm: cited in Section G for total Pol II levels of 50,000–200,000 molecules per HeLa cell (quantitative western blotting/mass spectrometry; RPB1 ≈ 150,000 per cell)\].

3. Boisvert, F.-M., Ahmad, Y., Gierliński, M., Charrière, F., Lamont, D., Scott, M., Barton, G., & Lamond, A. I. (2012). A quantitative spatial proteomics analysis of proteome turnover in human cells. *Molecular & Cellular Proteomics*, 11(3), M111.011429. https\://doi.org/10.1074/mcp.M111.011429

4. Caliskan, G., Hyeon, C., Perez-Salas, U., Briber, R. M., Woodson, S. A., & Thirumalai, D. (2005). Persistence length changes dramatically as RNA folds. *Physical Review Letters*, 95(26), 268303\. https\://doi.org/10.1103/PhysRevLett.95.268303

5. Cho, N. H., Cheveralls, K. C., Brunner, A.-D., … Mann, M., & Leonetti, M. D. (2022). OpenCell: Endogenous tagging for the cartography of human cellular organization. *Science*, 375(6585), eabi6983. https\://doi.org/10.1126/science.abi6983

6. Choi et al. (2011). \[UNVERIFIED — author to confirm: cited in Section B for the fully extended contour length of the 52-repeat (364-residue) CTD, ≈131–138 nm\].

7. Corden, J. L., Cadena, D. L., Ahearn, J. M., & Dahmus, M. E. (1985). A unique structure at the carboxyl terminus of the largest subunit of eukaryotic RNA polymerase II. *Proceedings of the National Academy of Sciences*, 82(23), 7934–7938. https\://doi.org/10.1073/pnas.82.23.7934

8. Czudnochowski, N., Bösken, C. A., & Geyer, M. (2012). Serine-7 but not serine-5 phosphorylation primes RNA polymerase II CTD for P-TEFb recognition. *Nature Communications*, 3, 842. https\://doi.org/10.1038/ncomms1846

9. Darabi et al. (2022). \[UNVERIFIED — author to confirm: cited in Section B for the Gaussian end-to-end (J-density) expression for the effective local concentration of a tethered chain end\].

10. Darzacq, X., Shav-Tal, Y., de Turris, V., Brody, Y., Shenoy, S. M., Phair, R. D., & Singer, R. H. (2007). *In vivo* dynamics of RNA polymerase II transcription. *Nature Structural & Molecular Biology*, 14(9), 796–806. https\://doi.org/10.1038/nsmb1280

11. Davidson, L., Muniz, L., & West, S. (2014). 3′ end formation of pre-mRNA and phosphorylation of Ser2 on the RNA polymerase II CTD are reciprocally coupled in human cells. *Genes & Development*, 28(4), 342–356. https\://doi.org/10.1101/gad.231274.113

12. Dhusia, K., & Wu, Y. (2021). Classification of protein–protein association rates based on biophysical informatics. *BMC Bioinformatics*, 22(1), 408\. https\://doi.org/10.1186/s12859-021-04323-0

13. Fong, N., Brannan, K., Erickson, B., Kim, H., Cortazar, M. A., Sheridan, R. M., Nguyen, T., Karp, S., & Bentley, D. L. (2015). Effects of transcription elongation rate and Xrn2 exonuclease activity on RNA polymerase II termination suggest widespread kinetic competition. *Molecular Cell*, 60(2), 256–267. https\://doi.org/10.1016/j.molcel.2015.09.026

14. Fong, N., Saldi, T., Sheridan, R. M., Cortazar, M. A., & Bentley, D. L. (2017). RNA Pol II dynamics modulate co-transcriptional chromatin modification, CTD phosphorylation, and transcriptional direction. *Molecular Cell*, 66(4), 546–557.e3. https\://doi.org/10.1016/j.molcel.2017.04.016

15. Fu et al. (2026). \[UNVERIFIED — author to confirm: cited in Section A for diffusion-limited protein–protein association rates of order 5 × 10⁵ M⁻¹ s⁻¹\].

16. Fuchs, G., Voichek, Y., Benjamin, S., Gilad, S., Amit, I., & Oren, M. (2014). 4sUDRB-seq: measuring genomewide transcriptional elongation rates and initiation frequencies within cells. *Genome Biology*, 15(5), R69. https\://doi.org/10.1186/gb-2014-15-5-r69

17. Hamilton, K., Sun, Y., & Tong, L. (2019). Biophysical characterizations of the recognition of the AAUAAA polyadenylation signal. *RNA*, 25(12), 1673–1680. https\://doi.org/10.1261/rna.070870.119

18. Hsin, J. P., & Manley, J. L. (2012). The RNA polymerase II CTD coordinates transcription and RNA processing. *Genes & Development*, 26(19), 2119–2137. https\://doi.org/10.1101/gad.200303.112

19. Hyeon, C., Dima, R. I., & Thirumalai, D. (2006). Size, shape, and flexibility of RNA structures. *The Journal of Chemical Physics*, 125(19), 194905\. https\://doi.org/10.1063/1.2364190

20. Itzhak, D. N., Tyanova, S., Cox, J., & Borner, G. H. H. (2016). Global, quantitative and dynamic mapping of protein subcellular localization. *eLife*, 5, e16950. https\://doi.org/10.7554/eLife.16950

21. Jeong, J., Le, T. T., & Kim, H. D. (2016). Single-molecule fluorescence studies on DNA looping. *Methods*, 105, 34–43. https\://doi.org/10.1016/j.ymeth.2016.04.005

22. Jonkers, I., Kwak, H., & Lis, J. T. (2014). Genome-wide dynamics of Pol II elongation and its interplay with promoter proximal pausing, chromatin, and exons. *eLife*, 3, e02407. https\://doi.org/10.7554/eLife.02407

23. Kim, H., Erickson, B., Luo, W., Seward, D., Graber, J. H., Pollock, D. D., Megee, P. C., & Bentley, D. L. (2010). Gene-specific RNA polymerase II phosphorylation and the CTD code. *Nature Structural & Molecular Biology*, 17(10), 1279–1286. https\://doi.org/10.1038/nsmb.1913

24. Krebs, A. R., Imanci, D., Hoerner, L., Gaidatzis, D., Burger, L., & Schübeler, D. (2017). Genome-wide single-molecule footprinting reveals high RNA polymerase II turnover at paused promoters. *Molecular Cell*, 67(3), 411–422.e4. https\://doi.org/10.1016/j.molcel.2017.06.027

25. Laitem, C., Zaborowska, J., Isa, N. F., Kufs, J., Dienstbier, M., & Murphy, S. (2015). CDK9 inhibitors define elongation checkpoints at both ends of RNA polymerase II-transcribed genes. *Nature Structural & Molecular Biology*, 22(5), 396–403. https\://doi.org/10.1038/nsmb.3000

26. Liao & Bhatt (2019). \[UNVERIFIED — author to confirm: cited in Section G for total Pol II levels per HeLa cell (RPB1 ≈ 150,000 molecules per cell)\].

27. Lunde, B. M., Reichow, S. L., Kim, M., Suh, H., Leeper, T. C., Yang, F., Mutschler, H., Buratowski, S., Meinhart, A., & Varani, G. (2010). Cooperative interaction of transcription termination factors with the RNA polymerase II C-terminal domain. *Nature Structural & Molecular Biology*, 17(10), 1195–1201. https\://doi.org/10.1038/nsmb.1893

28. Mayer et al. (2010). \[UNVERIFIED — author to confirm: cited in Section G for total Pol II levels per HeLa cell; not to be confused with Mayer et al. (2012), a different paper on CTD Tyr1 phosphorylation\].

29. Mei, K., Jin, Z., Ren, F., Wang, Y., Chang, Z., & Wang, X. (2014). Structural basis for the recognition of RNA polymerase II C-terminal domain by CREPT and p15RS. *Science China Life Sciences*, 57(1), 97–106. https\://doi.org/10.1007/s11427-013-4589-7

30. Meinhart, A., & Cramer, P. (2004). Recognition of RNA polymerase II carboxy-terminal domain by 3′-RNA-processing factors. *Nature*, 430(6996), 223–226. https\://doi.org/10.1038/nature02679

31. Milo, R., Jorgensen, P., Moran, U., Weber, G., & Springer, M. (2010). BioNumbers—the database of key numbers in molecular and cell biology. *Nucleic Acids Research*, 38(Database issue), D750–D753. https\://doi.org/10.1093/nar/gkp889 (cited in the text as "BioNumbers").

32. Monarez, R. R., MacDonald, C. C., & Dass, B. (2007). Polyadenylation proteins CstF-64 and τCstF-64 exhibit differential binding affinities for RNA polymers. *Biochemical Journal*, 401(3), 651–658. https\://doi.org/10.1042/BJ20061097

33. Moreno, R. Y., Panina, S. B., & Zhang, Y. J. (2025). RPRD1B's direct interaction with phosphorylated RNA polymerase II regulates polyadenylation of cell cycle genes and drives cancer progression. *RSC Chemical Biology*, 6(3), 423–437. https\://doi.org/10.1039/d4cb00212a

34. Ni, Z., Xu, C., Guo, X., Hunter, G. O., Kuznetsova, O. V., Tempel, W., … & Greenblatt, J. F. (2014). RPRD1A and RPRD1B are human RNA polymerase II C-terminal domain scaffolds for Ser5 dephosphorylation. *Nature Structural & Molecular Biology*, 21(8), 686–695. https\://doi.org/10.1038/nsmb.2853

35. Noble, C. G., Hollingworth, D., Martin, S. R., Ennis-Adeniran, V., Smerdon, S. J., Kelly, G., Taylor, I. A., & Ramos, A. (2005). Key features of the interaction between Pcf11 CID and RNA polymerase II CTD. *Nature Structural & Molecular Biology*, 12(2), 144–151. https\://doi.org/10.1038/nsmb887

36. Nojima, T., Gomes, T., Grosso, A. R. F., Kimura, H., Dye, M. J., Dhir, S., Carmo-Fonseca, M., & Proudfoot, N. J. (2015). Mammalian NET-seq reveals genome-wide nascent transcription coupled to RNA processing. *Cell*, 161(3), 526–540. https\://doi.org/10.1016/j.cell.2015.03.027

37. Northrup, S. H., & Erickson, H. P. (1992). Kinetics of protein-protein association explained by Brownian dynamics computer simulation. *Proceedings of the National Academy of Sciences*, 89(8), 3338–3342. https\://doi.org/10.1073/pnas.89.8.3338

38. Parua, P. K., Kalan, S., Benjamin, B., Sansó, M., & Fisher, R. P. (2020). Distinct Cdk9-phosphatase switches act at the beginning and end of elongation by RNA polymerase II. *Nature Communications*, 11, 4338. https\://doi.org/10.1038/s41467-020-18173-6

39. Phillips, R., Kondev, J., Theriot, J., & Garcia, H. G. (2013). *Physical biology of the cell* (2nd ed.). Garland Science.

40. Portz, B., Lu, F., Gibbs, E. B., Mayfield, J. E., Rachel Mehaffey, M., Zhang, Y. J., Brodbelt, J. S., Showalter, S. A., & Gilmour, D. S. (2017). Structural heterogeneity in the intrinsically disordered RNA polymerase II C-terminal domain. *Nature Communications*, 8, 15231\. https\://doi.org/10.1038/ncomms15231

41. Proudfoot, N. J. (2016). Transcriptional termination in mammals: stopping the RNA polymerase II juggernaut. *Science*, 352(6291), aad9926. https\://doi.org/10.1126/science.aad9926

42. Riback, J. A., Bowman, M. A., Zmyslowski, A. M., Plaxco, K. W., Clark, P. L., & Sosnick, T. R. (2019). Commonly used FRET fluorophores promote collapse of an otherwise disordered protein. *Proceedings of the National Academy of Sciences*, 116(18), 8889–8894. https\://doi.org/10.1073/pnas.1813038116

43. Rosonina, E., & Blencowe, B. J. (2004). Analysis of the requirement for RNA polymerase II CTD heptapeptide repeats in pre-mRNA splicing and 3′-end cleavage. *RNA*, 10(4), 581–589. https\://doi.org/10.1261/rna.5207204

44. Schreiber, G., Haran, G., & Zhou, H.-X. (2009). Fundamental aspects of protein–protein association kinetics. *Chemical Reviews*, 109(3), 839–860. https\://doi.org/10.1021/cr800373w

45. Schüller, R., Forné, I., Straub, T., Schreieck, A., Texier, Y., Shah, N., Decker, T.-M., Cramer, P., Imhof, A., & Eick, D. (2016). Heptad-specific phosphorylation of RNA polymerase II CTD. *Molecular Cell*, 61(2), 305–314. https\://doi.org/10.1016/j.molcel.2015.12.003

46. Schwer, B., Ghosh, A., Sanchez, A. M., Lima, C. D., & Shuman, S. (2015). Genetic and structural analysis of the essential fission yeast RNA polymerase II CTD phosphatase Fcp1. *RNA*, 21(6), 1135–1146. https\://doi.org/10.1261/rna.050286.115

47. Suh, H., Ficarro, S. B., Kang, U.-B., Chun, Y., Marto, J. A., & Buratowski, S. (2016). Direct analysis of phosphorylation sites on the Rpb1 C-terminal domain of RNA polymerase II. *Molecular Cell*, 61(2), 297–304. https\://doi.org/10.1016/j.molcel.2015.12.021

48. Tellier, M., Zaborowska, J., Caizzi, L., Mohammad, E., Velychko, T., Schwalb, B., Ferrer-Vicens, I., Blears, D., Nojima, T., Cramer, P., & Murphy, S. (2020). CDK12 globally stimulates RNA polymerase II transcription elongation and carboxyl-terminal domain phosphorylation. *Nucleic Acids Research*, 48(14), 7712–7727. https\://doi.org/10.1093/nar/gkaa514

49. Tellier, M., Zaborowska, J., Neve, J., Nojima, T., Hester, S., Fournier, M., Furger, A., & Murphy, S. (2022). CDK9 and PP2A regulate RNA polymerase II transcription termination and coupled RNA maturation. *EMBO Reports*, 23(10), e54520. https\://doi.org/10.15252/embr.202154520

50. Veloso, A., Kirkconnell, K. S., Magnuson, B., Biewen, B., Paulsen, M. T., Wilson, T. E., & Ljungman, M. (2014). Rate of elongation by RNA polymerase II is associated with specific gene features and epigenetic modifications. *Genome Research*, 24(6), 896–905. https\://doi.org/10.1101/gr.171405.113

51. Vervoort, S. J., Welsh, S. A., Devlin, J. R., Barbieri, E., Knight, D. A., Offley, S., … Gardini, A., & Johnstone, R. W. (2021). The PP2A-Integrator-CDK9 axis fine-tunes transcription and can be targeted therapeutically in cancer. *Cell*, 184(12), 3143–3162.e32. https\://doi.org/10.1016/j.cell.2021.04.022

52. Vihervaara, A., Versluis, P., Himanen, S. V., & Lis, J. T. (2023). PRO-IP-seq tracks molecular modifications of engaged Pol II complexes at nucleotide resolution. *Nature Communications*, 14, 7039. https\://doi.org/10.1038/s41467-023-42715-3

53. Vos, S. M., Farnung, L., Boehning, M., Wigge, C., Linden, A., Urlaub, H., & Cramer, P. (2018). Structure of activated transcription complex Pol II–DSIF–PAF–SPT6. *Nature*, 560(7720), 607–612. https\://doi.org/10.1038/s41586-018-0440-4

54. Wang, M., Weiss, M., Simonovic, M., Haertinger, G., Schrimpf, S. P., Hengartner, M. O., & von Mering, C. (2012). PaxDb, a database of protein abundance averages across all three domains of life. *Molecular & Cellular Proteomics*, 11(8), 492–500. https\://doi.org/10.1074/mcp.O111.014704

55. Yang, W., Hsu, P. L., Yang, F., Song, J.-E., & Varani, G. (2018). Reconstitution of the CstF complex unveils a regulatory role for CstF-50 in recognition of 3′-end processing signals. *Nucleic Acids Research*, 46(2), 493–503. https\://doi.org/10.1093/nar/gkx1177

56. Zeiler, M., Straube, W. L., Lundberg, E., Uhlen, M., & Mann, M. (2012). A Protein Epitope Signature Tag (PrEST) library allows SILAC-based absolute quantification and multiplexed determination of protein copy numbers in cell lines. *Molecular & Cellular Proteomics*, 11(3), O111.009613. https\://doi.org/10.1074/mcp.O111.009613

57. Zhou, H.-X. (1993). Brownian dynamics study of the influences of electrostatic interaction and diffusion on protein-protein association kinetics. *Biophysical Journal*, 64(6), 1711–1726. https\://doi.org/10.1016/S0006-3495(93)81543-1

### **Additional references (range grounding, October 2026)**

These are cited in Sections B–D for the parameter ranges. Merge them into the numbered list above before submission.

- Beaudoing, E., Freier, S., Wyatt, J. R., Claverie, J.-M., & Gautheret, D. (2000). Patterns of variant polyadenylation signal usage in human genes. *Genome Research*, 10(7), 1001–1010. https://doi.org/10.1101/gr.10.7.1001
- Chao, L. C., Jamil, A., Kim, S. J., Huang, L., & Martinson, H. G. (1999). Assembly of the cleavage and polyadenylation apparatus requires about 10 seconds in vivo and is faster for strong than for weak poly(A) sites. *Molecular and Cellular Biology*, 19(8), 5588–5600. https://doi.org/10.1128/MCB.19.8.5588
- Grozdanov, P. N., et al. (2018). *Nucleic Acids Research*, 46(22), 12022–. https://doi.org/10.1093/nar/gky862
- Huang, L., Chu, H.-F., & Tong, L. (2025). Molecular basis for the recognition of low-frequency polyadenylation signals by mPSF. *Nucleic Acids Research*, 53(17), gkaf890. https://doi.org/10.1093/nar/gkaf890
- Jasnovidova, O., et al. (2017). *EMBO Reports*. https://doi.org/10.15252/embr.201643723
- Jasnovidova, O., et al. (2017). *PNAS*. https://doi.org/10.1073/pnas.1712450114
- Ling, Liang, Wang & Wu (2026). *Science*, 391, eads0960. https://doi.org/10.1126/science.ads0960
- Tian, B., Hu, J., Zhang, H., & Lutz, C. S. (2005). A large-scale analysis of mRNA polyadenylation of human and mouse genes. *Nucleic Acids Research*, 33(1), 201–212. https://doi.org/10.1093/nar/gki158
- Torres-Ulloa, L., Calvo-Roitberg, E., & Pai, A. A. (2024). *RNA*, 30, 256–270. https://doi.org/10.1261/rna.079783.123
- Yang, Q., Coseno, M., Gilmartin, G. M., & Doublié, S. (2011). Crystal structure of a human cleavage factor CFIm25/CFIm68/RNA complex provides an insight into poly(A) site recognition and RNA looping. *Structure*, 19(3), 368–377. https://doi.org/10.1016/j.str.2010.12.021

### **Note: sources to verify before submission**

References not cited in the text were removed from this list. The sources below still need the author's attention.

**Not yet identified (cited in the text, marked UNVERIFIED above)**

- Bhatt et al. (2012), Liao & Bhatt (2019) and Mayer et al. (2010), Section G: cited for 50,000–200,000 Pol II molecules per HeLa cell. No matching papers were found, and Mayer et al. (2010) in *Nat Struct Mol Biol* is a yeast occupancy study. A possible replacement is Kimura, H., Tao, Y., Roeder, R. G., & Cook, P. R. (1999). *Molecular and Cellular Biology*, 19(8), 5383–5392. https://doi.org/10.1128/MCB.19.8.5383. It reports about 140,000 soluble plus about 180,000 nuclear-bound Pol II molecules per HeLa cell (about 320,000 total), which is above the stated range.
- Choi et al. (2011), Section B: cited for the extended contour length of the 52-repeat CTD. Two 2011 Choi papers on disordered proteins fit; confirm which one is meant.
- Darabi et al. (2022), Section B: cited for the Gaussian end-to-end (J-density) expression. No match found.
- Fu et al. (2026), Section A: cited for diffusion-limited protein–protein association rates. No match found.

**Identified, but the source does not clearly support the claim**

- Jasnovidova et al. (2017), Section C. The 0.3 μM affinity of Rtt103 for the long phosphorylated CTD, which grounds the lower ${k}_{Eoff}$ bound, was seen only in an abstract or search snippet. Check it against the PNAS full text.
- Lunde et al. (2010), Section C. The PMC text gives the Rtt103 fluorescence-anisotropy ${K}_{d}$ for the doubly phosphorylated diheptad as "2.1"; it may be 21 μM. Check the PDF.
- Veloso et al. (2014), Section E: the paper reports a median elongation rate of about 1.5 kb/min. The text cites it for 3–5 kb/min.

**Matched by best judgement; confirm the intended paper**

- Moreno, Panina & Zhang (2025), *RSC Chemical Biology*: the only 2025 Moreno paper on CTD binding found.
- Vos et al. (2018), *Nature* 560:607: one of two 2018 Vos et al. papers; chosen because its abstract discusses the CTD linker.
