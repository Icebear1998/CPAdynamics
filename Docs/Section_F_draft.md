## **F. CTD Ser2 Phosphorylation Parameters (${k}_{Pon,min}$, ${k}_{Pon,slope}$, ${k}_{Poff}$)**

### **What p counts**

Each polymerase carries $p$ E-binding sites ($0 ≤ p ≤ M$). We define one site as two adjacent Ser2-phosphorylated CTD repeats (a Ser2P diheptad). CID-containing CTD readers bind Ser2P diheptads most tightly. For Rtt103, Ser2P on both repeats binds 5–30-fold more tightly than Ser2P on one repeat (15 vs 76–420 μM by NMR; 2.1 μM by fluorescence anisotropy), and neighbouring CIDs bind cooperatively on four-repeat peptides (Lunde et al., 2010). RPRD1A/1B bind Ser2P-containing segments as dimers that span two binding sites (Ni et al., 2014). The 1–10 μM affinities used for ${k}_{Eoff}$ in Section C are the affinities of such doubly phosphorylated segments. The extra mark matters less for Pcf11 (Lunde et al., 2010), so the diheptad definition is a coarse-grained choice rather than a property of every reader.

If $n$ Ser2P marks are scattered at random over the 52 human repeats, the expected number of adjacent pairs is $n(n-1)/52$. The default $M = 5$ sites therefore corresponds to $n ≈ 17$ Ser2P per CTD. This agrees with mass-spectrometry measurements: yeast Pol II carries 0.3–0.4 phosphates per repeat, mostly Ser2P and Ser5P (Suh et al., 2016), or about 16–21 on a 52-repeat CTD, and the human CTD is not heavily phosphorylated, with doubly phosphorylated repeats rare (Schüller et al., 2016).

Phosphorylation adds one site per event, at a rate ${k}_{Pon}$ that does not depend on $p$. Each of the $p-e$ unoccupied sites is lost independently at ${k}_{Poff}$, so a polymerase loses sites at $(p-e){k}_{Poff}$; E-bound sites are protected. The Ser2P profile in the model is therefore set by how sites accumulate as Pol II travels, not by a local two-state equilibrium.

### **Site formation (${k}_{Pon,min}$, ${k}_{Pon,slope}$)**

Ser2 kinases act distributively: P-TEFb releases the CTD after every phosphotransfer (Baumli et al., 2012; Czudnochowski et al., 2012), which favours phosphorylation of separated rather than consecutive repeats (Baumli et al., 2012). We therefore treat Ser2P deposition as random, at a constant rate ${k}_{dep}$ per CTD. A polymerase that has been on the gene for time $t$ then carries a Poisson number of marks $n$ with mean $\lambda = {k}_{dep}t$. Since $E[n(n-1)] = \lambda^{2}$ for a Poisson count, the expected number of pair sites is $({k}_{dep}t)^{2}/52$, and new pair sites form at rate

$\frac{d}{dt}\left[\frac{({k}_{dep}t)^{2}}{52}\right] = \frac{2{k}_{dep}^{2}t}{52}$.

This rises linearly with time from zero, and so linearly with distance at constant elongation rate. It gives the model's linear form

${k}_{Pon}(x) = {k}_{Pon,min} + {k}_{Pon,slope}\cdot x, \quad {k}_{Pon,slope} = \frac{2{k}_{dep}^{2}}{52\,{k}_{e}}, \quad {k}_{Pon,min} = 0$,

where $x$ is the node index from the TSS and ${k}_{e}$ = 0.65 node s⁻¹.

${k}_{Pon}$ is the rate at which E sites (pairs) form, not the per-repeat phosphorylation rate ${k}_{dep}$. A polymerase leaves the TSS with no marks and cannot carry a pair until it has acquired two, so ${k}_{Pon}$ is zero at the TSS even though ${k}_{dep}$ is not.

**Estimating ${k}_{dep}$.** Ser2P antibody signals rise from the promoter-proximal pause region and stop changing appreciably beyond about +1.5 kb in human K562 cells (Vihervaara et al., 2023), and do not scale with gene length in yeast (Kim et al., 2010). Antibody pull-down reports whether a polymerase carries recognizable Ser2P rather than how many marks it carries. We therefore read the plateau as about 90% of polymerases carrying at least one mark, a Poisson mean of $n ≈ 2.3$, by 1.5 kb. At 65 bp/s this takes 23 s, giving ${k}_{dep} ≈ 0.1$ s⁻¹. This is far below kinase capacity: CDK9/cyclin T turns over at 2.97 s⁻¹ in vitro (Baumli et al., 2012), and HEK293T cells contain roughly 0.6–0.8 CDK9 or CDK12 molecules per Pol II (Cho et al., 2022).

* ${k}_{Pon,slope}$ = 2(0.1)²/(52 × 0.65) ≈ 6 × 10⁻⁴; we use **${k}_{Pon,slope}$ = 5 × 10⁻⁴ s⁻¹ per node** (${k}_{dep}$ ≈ 0.09 s⁻¹).

* **${k}_{Pon,min}$ = 0**. A small positive value (0.001 s⁻¹) changes CAD₅₀ at the 20 kb PAS by less than 1 bp for M = 1, 5 and 10, and only slightly raises the number of E sites in the first few kilobases.

Without removal, about 17 marks would accumulate by 12 kb, so the Ser2P count keeps rising long after the antibody signal has plateaued. With ${k}_{Poff}$ = 0.003 s⁻¹ (below), E sites reach 50% of M at about 9 kb and 90% at about 14 kb for M = 5 (12.5 and 19 kb for M = 10). The number of E sites therefore keeps rising for roughly 10 kb after the antibody signal has plateaued.

The Ser2P peak downstream of the PAS is not used in this calibration. On human MYC in HeLa cells, Ser2P per Pol II rises mostly over the first ~1.6 kb and then roughly doubles just downstream of the pA site (Davidson et al., 2014, Fig. 1B, read by eye). That 3′ peak depends on a functional pA signal, on CDK12 and on Pol II pausing, so it is part of termination rather than of the pre-PAS phosphorylation rate. The model also cannot reproduce it, because it removes Pol II at commitment and so has no post-cleavage, paused Pol II.

**Range:** ${k}_{Pon,slope}$ ∈ \[1 × 10⁻⁴, 2 × 10⁻³\] s⁻¹ per node. The lower bound allows a slower first mark: the model elongates at a constant 65 bp/s, whereas Pol II accelerates over the first kilobases in vivo (Jonkers et al., 2014). At this bound the mean of 2.3 marks is reached at about 3.6 kb. At the upper bound it is reached at about 0.8 kb, close to the antibody plateau.

### **Site loss (${k}_{Poff}$)**

${k}_{Poff}$ is the loss rate of each unoccupied site. A pair site is lost when either of its two marks is removed. If each mark is removed independently at rate ${k}_{r}$, then ${k}_{Poff} ≈ 2{k}_{r}$, and a polymerase with $p-e$ free sites loses sites at $(p-e){k}_{Poff}$. E-bound sites are protected.

No study reports a fitted Ser2P half-life on engaged Pol II. The available evidence:

* **CDK12 inhibition (most direct).** In HEK293 cells with analog-sensitive CDK12, Ser2P on chromatin-bound Pol II is markedly reduced within 15 min of inhibition and falls further by 30 min (Tellier et al., 2020).
  * Run-off cannot explain all of the loss. At the 3′ ends of genes longer than about 34 kb, elongation factors are unchanged at 15 min, so pre-existing Pol II is still there.
  * Assuming the kinase is fully off and the signal tracks the number of marks: at least 50% loss by 15 min gives ${k}_{r}$ ≳ 8 × 10⁻⁴ s⁻¹. Continued loss between 15 and 30 min gives roughly ${k}_{r}$ ≲ 3 × 10⁻³ s⁻¹.
* **The loss is active dephosphorylation.**
  * Ser2P loss after CDK9 inhibition requires PP2A recruited by the Integrator subunit INTS6. Without INTS6, Ser2P stays on gene bodies after normalising to total Pol II (Vervoort et al., 2021).
  * The PP2A inhibitor calyculin A prevents Ser2P loss after DRB treatment (Tellier et al., 2022).
* **CDK9 inhibition alone.** Ser2P on Pol II and in the GAPDH gene body is greatly reduced within 30 min (Laitem et al., 2015), but can persist for at least 1 h while CDK12/13 remain active (Parua et al., 2020). These time courses bound ${k}_{r}$ from below only.
* **Dwell time.** Slower Pol II shifts Ser2P toward the 5′ ends of genes (Fong et al., 2017). Ser2P therefore depends on the time Pol II spends on the gene, which requires turnover slow compared with transit.
* **Phosphatase capacity.** Fcp1 turns over Ser2-PO₄ CTD substrate at 0.13 s⁻¹ (Schwer et al., 2015), and HEK293T cells contain about 0.17 FCP1 (CTDP1) per Pol II (Cho et al., 2022). FCP1 alone could remove roughly 10⁻³ s⁻¹ per mark. PP2A adds to this.

Together these give ${k}_{r}$ ≈ 8 × 10⁻⁴–3 × 10⁻³ s⁻¹ (half-life about 4–15 min), or ${k}_{Poff}$ ≈ 0.0016–0.006 s⁻¹. We take the centre of this range, ${k}_{r}$ = 1.5 × 10⁻³ s⁻¹ (half-life ≈ 7.7 min):

${k}_{Poff}$ **= 0.003 s⁻¹**

**Range:** ${k}_{Poff}$ ∈ \[0.0016, 0.006\] s⁻¹. This narrows an earlier, wider interval of \[0.002, 0.02\] in light of Tellier et al. (2020, 2022) and Vervoort et al. (2021).

Consistency check: deposition and removal balance at $n^{*} = {k}_{dep}/{k}_{r}$ ≈ 60 marks at the chosen value, and 30–110 across the range. Marks therefore keep accumulating along a 20 kb gene, as the derivation of ${k}_{Pon}(x)$ assumes. Because $p$ is capped at M, this ceiling does not change the number of E sites.

### **Sensitivity**

All results below use kHoff = 1, E_total = 100,000, Pol_total = 100,000, kc = 0.14 and kEoff_engaged = 0.05.

**At a 20 kb PAS:**

* **Baseline (${k}_{Pon,slope}$ = 5 × 10⁻⁴, ${k}_{Poff}$ = 0.003).**
  * M = 5 polymerases carry about 4.9 sites at the PAS, with about 0.9 bound E.
  * CAD₅₀ is 661, 444, 370, 332, 306 and 290 bp for M = 1–6.
  * Free E is 78% for M = 1 and 46% for M = 5.
* **${k}_{Poff}$ (slope fixed).**
  * Across the range (0.0016–0.006), CAD₅₀ is 654–675 bp for M = 1 and 306–308 bp for M = 5, because E sites are already near M at the PAS.
  * Faster loss starts to matter above the range: 742 and 335 bp at 0.02, and 888 and 446 bp at 0.05, where M = 5 carries only about 2.3 sites.
* **${k}_{Pon,slope}$ (${k}_{Poff}$ fixed).** Across 10⁻⁴–2 × 10⁻³, CAD₅₀ is 661–744 bp for M = 1 and 306–425 bp for M = 5. Both are highest at the lowest slope, where M = 5 carries about 2.4 sites at the PAS.

**PAS distance from the TSS** (downstream window 5 kb):

| TSS–PAS (kb) | 1 | 2 | 2.5 | 5 | 10 | 20 | 40 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CAD₅₀, M = 1 (bp) | 4311 | 3331 | 2856 | 1184 | 678 | 661 | 664 |
| CAD₅₀, M = 5 (bp) | 4246 | 3262 | 2781 | 995 | 362 | 306 | 335 |
| E sites at PAS, M = 5 | 0.03 | 0.13 | 0.21 | 0.86 | 3.1 | 4.9 | 4.9 |

At 1–2.5 kb only 62–80% of polymerases cleave within the 5 kb window, so those CAD₅₀ values are approximate. The slope mainly controls where along the gene E sites approach M, so it matters most for PASs within about 15 kb of the TSS.

### **Open issues**

* **Chao benchmark window.** With this calibration, M = 1 lands at 661 bp at a 20 kb PAS, inside the 400–800 bp window of Chao et al. (1999), so it does not separate from multi-site polymerases. Within the literature ${k}_{Poff}$ range, M = 1 exceeds 800 bp only with a slower ramp. For example, slope 7 × 10⁻⁵ with ${k}_{Poff}$ = 0.003 gives 849 bp for M = 1 and 507–605 bp for M = 2–6. The model code currently uses a provisional ${k}_{Poff}$ = 0.05 instead (888 bp for M = 1, 446 bp for M = 5). This file records the literature-based calibration.
* **Benchmark geometry.** Chao et al. used a pRSVcat reporter (RSV LTR → CAT → SV40 small-t intron and early poly(A) region) transiently transfected into COS cells. The paper does not state the TSS–PAS distance; from the vector layout it is about 1.5–2.5 kb. At that distance the pair-site model gives fewer than 0.2 E sites per polymerase at the PAS, and CAD₅₀ of about 3 kb for any M. A 20 kb benchmark assumes the sites have already accumulated.
* **Definition of an E site.** Antibody signals plateau early: by +1.5 kb in K562 cells and over the first ~1.6 kb on MYC in HeLa cells. If one Ser2P mark is enough to bind E, E sites would also saturate within a few kb. CAD would then be nearly independent of PAS distance beyond about 2–3 kb, and gene-length (PAS-distance) effects would be confined to very proximal PASs. Under the pair-site definition used here, E sites keep accumulating for about 10–15 kb.
