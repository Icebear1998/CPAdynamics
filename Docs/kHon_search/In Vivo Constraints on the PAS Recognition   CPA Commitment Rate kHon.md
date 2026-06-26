# In Vivo Constraints on the PAS Recognition / CPA Commitment Rate kHon

## Executive Overview

This report synthesizes in vivo and in vitro constraints on the **PAS recognition / CPA commitment rate constant (kHon)**, defined here as the effective first-order rate at which an elongating RNA Pol II complex, having just transcribed the PAS, commits to cleavage and polyadenylation via tethered encounter of pre-loaded CPSF factors with the AAUAAA hexamer on the nascent RNA.

The most informative direct measurement comes from the cis‑antisense rescue assay of Chao et al. on the SV40 early poly(A) site, which estimates a commitment/assembly time of roughly **10–20 s** with a 50% commitment distance of ~200 bp at 40 nt/s, implying an effective first-order rate on the order of **0.05–0.1 s⁻¹** for a moderately strong mammalian PAS. Genome-wide nascent-RNA kinetic profiling in Drosophila S2 cells yields median pre-mRNA 3′ cleavage half-lives of **~36–42 s**, corresponding to overall cleavage rate constants ~0.016–0.019 s⁻¹, with a tail of faster sites cleaved in ~30 s or less. Live imaging of transcription with MS2/PP7 in fly embryos reports cleavage times of **1.5–3 min** for a developmental reporter, implying slower, context-dependent CPA for that construct. Together, these data suggest that **kHon in vivo is unlikely to be slower than ~0.01 s⁻¹ and is plausibly in the range 0.05–0.1 s⁻¹ for strong canonical PASs in mammalian-like settings**, with substantial gene-to-gene variation.[^1][^2][^3][^4]

High-affinity binding of the human CPSF160–WDR33–CPSF30 ternary complex to AAUAAA (Kd ~3 nM) together with diffusion-limited association (kon ~10⁶ M⁻¹ s⁻¹) implies an intrinsic microscopic on-rate on the order of **1–10 s⁻¹** for complex formation in solution, substantially faster than the in vivo commitment time inferred from transcription-based assays. This gap is consistent with kHon reflecting not just initial CPSF–hexamer contact but a multistep commitment process (assembly, remodeling, engagement with additional CPA components, and transition to a cleavage-competent state), as emphasized by Chao et al.[^2][^5][^6][^7]

Finally, structural and biophysical work on the Pol II CTD indicates a phosphorylated CTD tail length of the order of **60–70 nm** when extended, compared to ~10 nm when compact, consistent with a flexible, disordered tether capable of sweeping a substantial volume around the RNA exit channel. Although direct J‑factor measurements for CTD‑tethered encounters with nascent RNA are not available, this geometry is compatible with effective local molarities in the sub‑millimolar range, which when multiplied by diffusion-limited kon reproduce the order of magnitude for the inferred in vivo commitment rate.[^8]

***

## 1. Direct in Vivo Kinetic Constraints on CPA Commitment

### 1.1 cis-Antisense Rescue Assay for SV40 Early PAS (Chao et al. 1999)

Chao et al. developed a **cis‑antisense rescue assay** in COS cells to measure the time required for commitment of the SV40 early poly(A) signal to cleavage and polyadenylation in vivo. An inverted copy of the early PAS placed immediately downstream of the authentic site produces an antisense RNA that forms a duplex with the PAS region and blocks cleavage/polyadenylation, drastically reducing reporter expression and leading to accumulation of uncleaved nuclear pre‑mRNA.[^2]

By moving the inverted antisense PAS progressively further downstream, they observed **distance-dependent rescue** of polyadenylation: once the polymerase had transcribed sufficiently far beyond the PAS, commitment to cleavage/polyadenylation occurred before the antisense sequence was synthesized, so antisense inhibition was relieved. Quantitative analysis of CAT expression and RNase protection showed that for the SV40 early site, **50% rescue occurred at ~200 bp separation**, and the rescue curve was well fit by a first-order process as a function of transcription distance.[^2]

Using a transcription elongation rate of **~40 nt/s** (2.4 kb/min) for mammalian Pol II, Chao et al. estimated that commitment and full assembly of the cleavage/polyadenylation apparatus for the SV40 early site required approximately **10–20 s**, with **50% commitment by ~5 s** (200 bp/40 nt/s) and a maximum domain of assembled RNA–protein interactions spanning ~200 nt around the PAS. Fitting the distance dependence of rescue with a first-order rate law yielded an exponential parameter **k ≈ 4.9 × 10⁻³ bp⁻¹**, which converts to an effective first-order rate in time of **k_eff ≈ 0.196 s⁻¹** when multiplied by 40 nt/s.[^2]

Importantly, this **0.2 s⁻¹** parameter should be interpreted as the apparent rate of **commitment/assembly of the full CPA complex** in vivo for a moderately strong viral PAS, not simply microscopic CPSF–AAUAAA binding. Chao et al. further showed that the strong SV40 late poly(A) site rescued itself from antisense inhibition much faster than the early or synthetic sites; in a replicating plasmid background, the late site was largely committed with an antisense sequence only 91 bp downstream, consistent with even larger effective k values for strong PASs.[^2]

For the kHon parameter in the kHon–kHoff PAS recognition module, an appropriate mapping from the Chao data is:

- **SV40 early PAS**: kHon on the order of **0.05–0.2 s⁻¹** (half-commitment in ~5–15 s, completion ~20 s), acknowledging uncertainties in exact elongation rate and the multistep nature of assembly.
- **Strong SV40 late PAS**: qualitatively several-fold faster commitment than early/synthetic sites at comparable factor abundance, consistent with effective kHon values **≥ 0.2 s⁻¹**.[^2]

### 1.2 Genome-wide 3′ Cleavage Half-Lives from 4sU Nascent RNA (Torres-Ulloa et al. 2024)

Torres-Ulloa et al. developed a genome-wide kinetic assay for pre-mRNA 3′ end cleavage using progressive 4sU metabolic labeling, nascent RNA-seq, and mathematical modeling in Drosophila S2 cells. Informative reads overlapping cleavage sites (CSs) or upstream windows were used to estimate the ratio of uncleaved to cleaved molecules as a function of labeling time, yielding site-specific **3′ cleavage half-lives** under steady-state conditions.[^4]

For ~2857 constitutive CSs (PAU > 95%), they reported median cleavage rate constants **k ≈ 0.019 min⁻¹**, corresponding to a median half-life **t₁/₂ ≈ 35.8 s**. For ~1601 alternative CSs (5% < PAU < 95%), median k ≈ 0.016 min⁻¹, giving **t₁/₂ ≈ 42.1 s**. Some CSs, such as that of PHGPx, had half-lives as short as **29.0 s**, while slower sites like the barricade (barc) CS had **t₁/₂ ≈ 2.1 min**. An additional subset of constitutive and many alternative sites had no detectable uncleaved reads at any time point, implying cleavage faster than the temporal resolution of the experiment (5 min labeling), and thus **t₁/₂ < several minutes**.[^4]

These 3′ cleavage half-lives report on the **overall kinetics from PAS recognition through endonucleolytic cleavage**, convolved over possible alternative sites within a gene, rather than the specific microscopic commitment step immediately after PAS transcription. Nevertheless, they place a robust **upper bound** on the in vivo CPA cycle time and suggest that for many sites cleavage is completed **within 30–60 s** of pre‑mRNA synthesis, consistent with earlier biochemical estimates that polyadenylation occurs within ~1 min of mRNA biogenesis in mammalian cells.[^4]

If kHon is the rate-limiting step for commitment, one expects kHon to be at least comparable to or somewhat faster than these half-life-derived rates (0.016–0.019 s⁻¹), implying **kHon ≳ 0.02 s⁻¹** on average, and higher for strong sites or high factor abundance.

### 1.3 Live Imaging of Cleavage Times in MS2/PP7 Reporters

Liu et al. applied a two-color MS2/PP7 live-imaging system and Bayesian inference to dissect the full transcription cycle of a hunchback reporter gene in Drosophila embryos, including initiation, elongation, and a **cleavage time parameter τ_cleave** representing the delay between Pol II reaching the end of the gene and release of the transcript.[^3]

Using dual-color traces and a model in which each polymerase traverses the gene at a constant elongation rate then experiences a deterministic cleavage delay, they inferred τ_cleave values on the order of **1.5–3 min**, with longer cleavage times in the anterior region (~3 min) and shorter (~1.5 min) posteriorly. The inferred mean elongation rate was ~1.72 kb/min, consistent with earlier Drosophila measurements.[^3]

These cleavage times are substantially slower than the median half-lives from the 4sU-based kinetic profiling, suggesting that in this particular reporter and developmental context, **post-elongation events (cleavage and release) can be rate limiting on the timescale of minutes**, possibly due to suboptimal PAS configuration, additional regulatory steps, or coupling to splicing and nuclear export. For the purpose of constraining kHon, these results reinforce that the **fast end of the CPA timescale is tens of seconds, but the slow end extends to minutes**, depending on gene architecture and regulatory context.[^3][^4]

### 1.4 Early Biochemical and Imaging Estimates

Earlier biochemical pulse-labeling work measuring the addition of radiolabeled adenines to nascent mRNA 3′ ends in mammalian cells suggested that polyadenylation (and thus cleavage) generally completed within **~1 min** of transcript synthesis, consistent with rapid, co-transcriptional CPA. More recent chromatin-associated RNA and imaging studies in mammalian systems support rapid 3′ end processing and demonstrate kinetic coupling between splicing and cleavage, although they often lack the temporal resolution to resolve sub-10-s events at single loci.[^9][^4]

Collectively, these in vivo measurements support an overall **CPA timescale from PAS transcription to completed cleavage** lying primarily between **~10 s and ~2–3 min**, with **strong constitutive PASs** at the fast end and alternative or regulated PASs often slower.[^9][^3][^4][^2]

***

## 2. CPSF–AAUAAA Recognition Kinetics and Binding Affinity

### 2.1 Equilibrium Affinity of CPSF Core for AAUAAA

Structural and biophysical studies have established that **a core CPSF complex comprising CPSF160, WDR33, and CPSF30, together with Fip1, directly recognizes the canonical AAUAAA PAS**. Cryo-EM structures show that CPSF30 and WDR33 make direct base-specific contacts with the AAUAAA hexamer within the core CPSF assembly.[^5][^10][^7][^11]

Hamilton, Sun, and Tong used fluorescence polarization to measure binding of the human CPSF160–WDR33–CPSF30 ternary complex to PAS RNA, reporting that AAUAAA is recognized with an affinity of **~3 nM** under their conditions. Variants of the RNA sequence or mutations in CPSF30 residues that engage the bases reduced affinity by more than an order of magnitude, consistent with PAS sequence playing a central role in PAS “strength”. Similarly, Clerici et al. reported equilibrium binding affinities in the **3–50 nM** range depending on CPSF complex composition and RNA sequence using FP assays.[^7][^5]

These high affinities imply that under nuclear protein concentrations (~0.1–1 μM for many factors), binding of CPSF core to a canonical AAUAAA sequence is **essentially saturated at thermodynamic equilibrium**, and that kinetic constraints on CPA are more likely to arise from **association rates, competitive binding, or downstream assembly steps** than from equilibrium occupancy per se.[^5][^7]

### 2.2 Association and Dissociation Rates

Direct measurements of the **kon and koff** for CPSF core binding to AAUAAA are sparse. Hamilton et al. focused on equilibrium FP and did not report explicit rate constants, but given Kd ≈ 3 nM and typical diffusion-limited association rates of **kon ~10⁶ M⁻¹ s⁻¹** for macromolecular interactions in the nucleus, one can infer that **koff ≈ Kd × kon ≈ 3 × 10⁻³ s⁻¹**, corresponding to a bound-state lifetime of **~300 s** (5 min) in vitro.[^12][^6][^7]

This simple estimate indicates that **once formed, the CPSF–PAS complex is relatively stable on the timescale of tens to hundreds of seconds**, consistent with the notion that PAS “strength” can depend on both high affinity and cooperative recruitment of additional factors. For weaker noncanonical PAS motifs with Kd in the tens of nM range, koff may be faster (~0.01–0.1 s⁻¹) but still slow compared to the sub-minute cleavage times observed at many sites.[^7][^9][^4]

In the context of a **tethered intramolecular encounter** between CTD-bound CPSF and a nascent AAUAAA hexamer, the **effective association rate kHon** should be approximated as

\(k_{\text{Hon}} \approx k_{\text{on}}^{\text{(diff)}} \times C_{\text{eff}}\)[^1]

where \(k_{\text{on}}^{\text{(diff)}}\) is the bimolecular diffusion-limited association rate (~10⁶ M⁻¹ s⁻¹) and \(C_{\text{eff}}\) is the effective local molar concentration (J‑factor) of CPSF relative to the PAS. For example, an effective molarity of **0.1 mM** would give kHon ≈ 100 s⁻¹, whereas **1 mM** would give kHon ≈ 1000 s⁻¹, both much faster than the 0.05–0.2 s⁻¹ commitment rates inferred in vivo.[^6][^2]

This discrepancy implies that **CPSF–AAUAAA recognition is not the dominant rate-limiting step in CPA in vivo**, and that kHon in the mechanistic model should represent a **coarse-grained rate that includes additional conformational and assembly transitions** (e.g., recruitment of CFIm/CstF/CF I/CF II, engagement of CPSF73, productive handoff among CTD, nascent RNA, and other cofactors), not just the initial hexamer contact.[^12][^9][^2]

### 2.3 Structural Basis for Recognition and Its Implications for kHon

Cryo-EM and crystallographic studies of the core CPSF complex show that CPSF30 zinc finger domains and an N‑terminal Lys/Arg-rich motif in WDR33 jointly contact AAUAAA bases, forming a compact recognition module that clamps the hexamer within a positively charged pocket. The CPSF160 scaffold organizes the positioning of CPSF30 and WDR33 and bridges to Fip1 and additional subunits.[^13][^10][^14]

The structural compactness and pre-organization of the PAS-binding site are consistent with **rapid association once the correct RNA sequence is presented within reach**, supporting the notion that the **intramolecular encounter per se could be very fast**, while the overall commitment step in vivo is slowed by recruitment and coordination of the full CPA machinery, as the Chao cis‑antisense data suggest.[^10][^2]

***

## 3. CTD Tethering, Polymer Properties, and Effective Local Concentration

### 3.1 CTD Length and Structural Range

The largest Pol II subunit Rpb1 contains a C‑terminal domain (CTD) composed of 26–52 tandem heptad repeats (YSPTSPS), with vertebrate CTDs containing 52 repeats, 21 of which are perfect consensus and the rest degenerate. Structural and biochemical data indicate that the **unphosphorylated yeast CTD can adopt a compact spiral of length ~100 Å (~10 nm)**, whereas the **phosphorylated CTD extends into a flexible tail ~650 Å (~65 nm) long**, roughly four times the diameter of the Pol II core.[^8]

This extended, intrinsically disordered tail provides multiple binding sites for RNA-processing factors, including CPSF, CstF, and other CPA components, and can reach over a large volume surrounding the RNA exit channel. Genome-wide ChIP studies show that Ser2 phosphorylation accumulates toward gene ends while Ser5 and Ser7 phosphorylation dominate promoters and early elongation regions, consistent with increasing recruitment of 3′-end processing factors as Pol II progresses into the 3′ flanking region.[^15][^8]

### 3.2 CTD-Mediated Tethering and J-Factor Considerations

Quantitative estimates of the **effective molarity or J‑factor** for CTD-tethered factors near the RNA exit channel are not directly available, but conceptual arguments can be made. The CTD behaves as a flexible, disordered polymer whose end-to-end distance distribution depends on repeat number, phosphorylation state, and interactions with other proteins.[^16][^17][^8]

If CPSF subunits are bound along multiple Ser2‑P heptads in a segment of length on the order of **tens of nanometers**, the instantaneous position of any given CPSF molecule relative to the RNA exit channel samples a **confined three-dimensional volume** whose radius is limited by the CTD contour length (~65 nm for a fully extended vertebrate CTD) and steric constraints near Pol II. Under such confinement, standard polymer-physics estimates for effective molarity of a tethered ligand in a reaction volume of tens to hundreds of nm³ give **sub‑millimolar to millimolar local concentrations**, orders of magnitude higher than bulk nuclear protein concentrations (~1–10 μM).[^15][^8]

While these numbers are approximate, they support the idea that **CTD tethering strongly enhances the local effective concentration of CPSF near the emerging PAS**, such that the **intramolecular encounter and binding step itself can be very fast** relative to the overall CPA cycle. In this regime, kHon is governed less by diffusion and more by **conformational search of the CTD and RNA, competitive binding, and multivalency of interacting complexes**.[^15][^8][^2]

### 3.3 CTD as an Organizer of CPA Factors

Reviews on CTD function emphasize that the CTD serves as a central platform that coordinates recruitment of capping enzymes, splicing factors, and **3′-end processing factors including CPSF and CstF**, with Ser2 phosphorylation being especially important for 3′ processing and termination. CTD phosphatase Ssu72, which is part of the cleavage and polyadenylation factor (CPF/CF) complex in yeast and has a human homolog Symplekin, dephosphorylates Ser5/Ser7-P near gene ends, linking CTD modification status to 3′-end processing.[^9][^8][^15]

ChIP and proteomics data show that CPSF and CstF can be recruited early and travel with Pol II, accumulating in a **3′ pause region 0.5–1.5 kb downstream of poly(A) sites**, consistent with early loading on the CTD and later engagement with the PAS. These observations support mechanistic models in which **Ser2P-enriched CTD segments carrying pre-loaded CPSF subunits provide a dynamic, high-local-concentration cloud of E-factors** that can rapidly engage a nascent PAS via tethered intramolecular collisions.[^18][^8]

***

## 4. Integrating Constraints into an In Vivo Range for kHon

### 4.1 Lower Bounds from Cleavage Half-Lives

Genome-wide kinetic profiling in S2 cells yields median 3′ cleavage half-lives of **~36–42 s** for constitutive and alternative sites, corresponding to overall cleavage rate constants **k_cleave ~0.016–0.019 s⁻¹**. If the PAS recognition/commitment step (kHon) is **rate-limiting or co‑limiting** in this overall process, then kHon must be **at least on the order of 10⁻² s⁻¹**. Faster sites (e.g., PHGPx, t₁/₂ ~29 s) imply k_cleave ~0.024 s⁻¹ and thus kHon ≳ 0.02–0.03 s⁻¹ in those contexts.[^4]

Earlier biochemical estimates that polyadenylation completes within ~1 min suggest that **k_cleave is rarely much slower than ~0.01 s⁻¹** for most active mRNA genes in mammalian-like systems, again implying **kHon ≳ 0.01 s⁻¹**.[^4]

### 4.2 Upper Bounds and Mechanistic Interpretation from Chao et al.

Chao et al.’s cis‑antisense rescue assay for the SV40 early poly(A) signal provides a more direct window on the commitment step itself, showing 50% commitment at ~200 bp (≈5 s) and completion in ~10–20 s, yielding an apparent first-order commitment rate of **k_eff ≈ 0.05–0.2 s⁻¹** depending on where in the curve one defines kHon. This rate incorporates **all steps between PAS emergence from Pol II and formation of a cleavage‑competent CPA apparatus** capable of resisting antisense occlusion.[^2]

Because the Chao system uses a strong viral PAS under high expression, in mammalian cells, and directly probes the competition between commitment and transcriptional progression, this **0.05–0.2 s⁻¹ window is highly relevant to mammalian kHon**. The observation that the SV40 late site assembles “several times faster” than early or synthetic sites suggests that **strong PAS motifs and favorable flanking elements can push kHon toward the upper end of this range or beyond**.[^2]

### 4.3 Reconciling CPSF Binding Kinetics with In Vivo Commitment

The high affinity of CPSF core for AAUAAA (Kd ~3 nM) and likely near-diffusion-limited kon (~10⁶ M⁻¹ s⁻¹) imply a microscopic association rate **kon·[CPSF] that is much larger than 0.1 s⁻¹** if the effective CPSF concentration near the PAS is in the μM range or higher through CTD tethering. For example, at 1 μM effective local concentration, the intrinsic CPSF association rate would be ~1 s⁻¹, and at 0.1 mM effective molarity, ~100 s⁻¹. These values are **one to three orders of magnitude faster than the in vivo commitment rates inferred from Chao’s work**.[^6][^7][^8][^2]

Thus, **kHon in the mechanistic model must be interpreted as a coarse-grained, effective rate constant for a multistep commitment process**, which includes not only the initial CPSF–AAUAAA encounter but also recruitment and positioning of other CPA factors, coordination with CTD phosphorylation state, and perhaps conformational changes in Pol II and the transcription complex. In this view, the CPSF binding step is fast and near-equilibrated, while one or more downstream assembly or remodeling steps are rate-limiting on the 5–20 s timescale.[^8][^9][^2]

### 4.4 Proposed Working In Vivo Range for kHon

Based on the assembled evidence, the following **working range for kHon in vivo** is consistent with current data:

- **Lower bound**: kHon **≳ 0.01 s⁻¹**
  - Justification: median cleavage half-lives of 36–42 s (k_cleave ~0.016–0.019 s⁻¹) in S2 cells and ~1 min polyadenylation in mammalian cells; commitment cannot be significantly slower than overall cleavage.[^4]

- **Plausible central range for canonical, reasonably strong PASs**: kHon **~0.05–0.1 s⁻¹**
  - Justification: Chao et al. SV40 early site shows 50% commitment in ~5 s and completion by 10–20 s (k_eff ~0.05–0.2 s⁻¹), with the early PAS representing a moderate-strength viral site.[^2]

- **Upper bound for very strong PASs in favorable contexts**: kHon **up to a few × 0.1 s⁻¹**
  - Justification: SV40 late PAS assembles “several times faster” than early and synthetic PASs in the same assay, implying significantly larger kHon; however, direct quantitative estimates are not provided, so **kHon ≲ 0.5 s⁻¹** is a conservative upper bound consistent with multi-second commitment.[^2]

Within this range, kHon should be modulated by:

- **PAS sequence and flanking cis-elements**, which influence CPSF affinity and recruitment of auxiliary factors.[^7][^9]
- **Local density and composition of PAS motifs**, which Torres-Ulloa et al. show increase cleavage rates cooperatively rather than additively.[^4]
- **E-factor loading level on the CTD** (⟨nE⟩), which will scale the effective collision cross-section and thereby kHon.
- **Elongation rate and 3′ Pol II pausing**, which shape the time window during which PAS recognition and commitment can occur before Pol II moves too far downstream.[^18][^9][^4]

***

## 5. Implications for Mechanistic Modeling and APA

### 5.1 Coupling kHon to CPA Assembly Distance (CAD)

In the mechanistic model where Pol II transits from R to REH at rate kHon after PAS transcription, the **mean CPA assembly distance (CAD)** downstream of the PAS is approximately \(v_{\text{elong}} / k_{\text{Hon}}\) when other steps are lumped into kHon. With v_elong ≈ 40 nt/s (mammalian) and kHon ≈ 0.05–0.1 s⁻¹, CAD would be **400–800 nt**, on the same order as empirically inferred pause and processing regions located 0.5–1.5 kb downstream of poly(A) sites in mammalian genes.[^18][^4][^2]

For stronger PASs (kHon approaching 0.2–0.5 s⁻¹), CAD would shrink to **80–200 nt**, consistent with Chao’s finding that the strong SV40 late site can assert commitment within <100 bp downstream of the PAS. This behavior matches the modeling intuition that **higher kHon produces shorter CAD and more efficient termination**.[^2]

### 5.2 kHon, kHoff, and PAS Strength in APA

In the APA module, the **ratio kHon/kHoff** defines the effective PAS “strength” at steady state in a rapid-equilibrium picture. The **in vitro CPSF affinity measurements (Kd ~3–50 nM)** and inferred koff (~10⁻³–10⁻² s⁻¹) imply that **kOn_micro/kOff_micro** is very large at the microscopic level for canonical AAUAAA sites.[^5][^7]

However, the effective in vivo **kHon** in the model acts on a different coarse-grained coordinate (commitment to cleavage and termination downstream of the PAS) and should be **calibrated against in vivo data such as the Chao antisense assay and genome-wide cleavage kinetics**, rather than the microscopic CPSF kon. kHoff, representing reversal of commitment, is likely small compared to kHon in strong PAS contexts but could be comparable for weak or regulated PASs where readthrough and delayed CPA are observed.[^19][^9][^4][^2]

Torres-Ulloa et al. find that alternative and upstream CSs are cleaved more slowly than constitutive and downstream CSs, and that genes with multiple PASs exhibit delayed cleavage at individual sites, consistent with **kinetic competition among sites and potential waiting until downstream PASs are transcribed**. In a multi-PAS model, **relative kHon values at proximal versus distal PASs, together with elongation rate and kHoff, will determine APA patterns**, and the in vivo ranges summarized here provide a quantitative starting point for such parameterization.[^4]

***

## 6. Gaps and Future Experimental Needs

Despite substantial progress, there remain important gaps in constraining kHon:

1. **Direct single-molecule measurements of CPSF–PAS association and dissociation rates** (kon, koff) in conditions mimicking nuclear crowding and CTD tethering are limited; current estimates rely on equilibrium FP and generic diffusion-limited kon values.[^6][^7]
2. **Live-cell measurements of the time between PAS transcription and irreversible commitment** (distinct from the time to physical cleavage) at single mammalian loci are scarce. Existing MS2/PP7 studies report overall cleavage delays but convolve multiple steps.[^20][^3]
3. **Quantitative J‑factor estimates for CTD-tethered factors** have not been experimentally determined. Polymer-physics modeling of the CTD based on recent structural heterogeneity and phase-separation studies could provide more precise effective molarity estimates.[^17][^16][^8]
4. **Systematic comparison of PAS variants and flanking elements** in matched live-cell reporters, with direct readouts of PAS-to-cleavage times and Pol II pausing, would help map how sequence-level PAS “strength” modulates kHon and CAD.

Addressing these gaps would allow tighter bounding of kHon and a more detailed mapping between microscopic binding steps and coarse-grained commitment in mechanistic models of transcription termination and APA.

---

## References

1. [Assembly of the Cleavage and Polyadenylation Apparatus Requires ...](https://pmc.ncbi.nlm.nih.gov/articles/PMC84411/) - We have devised a cis-antisense rescue assay of cleavage and polyadenylation to determine how long i...

2. [CPSF30 and Wdr33 directly bind to AAUAAA in mammalian mRNA 3](https://escholarship.org/uc/item/3w5260v8) - The CPSF subunit CPSF160 has been implicated in AAUAAA recognition, but direct evidence has been lac...

3. [[PDF] Real-time single-cell characterization of the eukaryotic transcription ...](https://garcialab.berkeley.edu/publications/Liu2021.pdf) - Precisely interpreting an MS2 or. PP7 signal therefore demands an integrated approach that accounts ...

4. [Genome-wide kinetic profiling of pre-mRNA 3′ end cleavage - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC10870368/) - We find that 3′ end cleavage is fast on average, with half-lives under a minute, but highly variable...

5. [Structural insights into the assembly and polyA signal recognition ...](https://pmc.ncbi.nlm.nih.gov/articles/PMC5760199/) - CPSF30 ZF3 domain and WDR33 N-terminal motif are principal determinants of AAUAAA motif. Equilibrium...

6. [end cleavage and polyadenylation of messenger RNA precursors.](https://www.semanticscholar.org/paper/The-biochemistry-of-3'-end-cleavage-and-of-RNA-Wahle-Keller/7a46a69bfba2fb4c676889d5e4e0dc5bea31526d) - The results confirmed that the polyadenylation-to-termination process is a straightforward two-step ...

7. [Biophysical characterizations of the recognition of the AAUAAA polyadenylation signal - PubMed](https://pubmed.ncbi.nlm.nih.gov/31462423/) - Most eukaryotic messenger RNA precursors must undergo 3'-end cleavage and polyadenylation for matura...

8. [The RNA polymerase II CTD coordinates transcription and ... - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3465734/) - The CTD functions to help couple transcription and processing of the nascent RNA and also plays role...

9. [Emerging Roles of RNA 3′-end Cleavage and Polyadenylation in ...](https://pmc.ncbi.nlm.nih.gov/articles/PMC7356254/) - A crucial feature of gene expression involves RNA processing to produce 3′ ends through a process te...

10. [Structural basis of AAUAAA polyadenylation signal recognition by ...](https://pubmed.ncbi.nlm.nih.gov/29358758/) - Mammalian mRNA biogenesis requires specific recognition of a hexanucleotide AAUAAA motif in the poly...

11. [CPSF30 and Wdr33 directly bind to AAUAAA in mammalian mRNA 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC4215182/) - Chan et al. found that CPSF subunits CPSF30 and Wdr33 directly contact AAUAAA. The CPSF30–RNA intera...

12. [Cleavage and polyadenylation: Ending the message expands gene ...](https://pmc.ncbi.nlm.nih.gov/articles/PMC5546720/) - Cleavage and polyadenylation (pA) is a fundamental step that is required for the maturation of prima...

13. [RCSB PDB - 6F9N: CRYSTAL STRUCTURE OF THE HUMAN CPSF160-WDR33 COMPLEX](https://www.rcsb.org/structure/6F9N) - CRYSTAL STRUCTURE OF THE HUMAN CPSF160-WDR33 COMPLEX

14. [Structural insights into the assembly and polyA signal ...](https://elifesciences.org/articles/33111.pdf)

15. [RNA Polymerase II C-Terminal Domain: Tethering Transcription to ...](https://pubs.acs.org/doi/10.1021/cr400158h) - The RNA polymerase II (Pol II) C-terminal domain (CTD) is a repetitive disordered domain that extend...

16. [RNA Pol II Length and Disorder Enable Cooperative Scaling of ...](https://www.biorxiv.org/content/10.1101/825299v1.full-text) - RNA Polymerase II contains a disordered C-terminal domain (CTD) whose length enigmatically correlate...

17. [Structural heterogeneity in the intrinsically disordered RNA ...](https://pubmed.ncbi.nlm.nih.gov/28497792/) - RNA polymerase II contains a repetitive, intrinsically disordered, C-terminal domain (CTD) composed ...

18. [RNA polymerase II pauses and associates with pre-mRNA ... - PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2836588/) - Transcription activation had little effect on CPSF at the start site, but enhanced its recruitment d...

19. [Elevated pre-mRNA 3′ end processing activity in cancer cells ...](https://www.nature.com/articles/s41467-023-39793-8) - Cleavage and polyadenylation (CPA) is responsible for 3′ end processing of eukaryotic poly(A)+ RNAs ...

20. [Kinetic competition during the transcription cycle results in stochastic ...](https://elifesciences.org/articles/03939) - We find that kinetic competition results in multiple competing pathways for pre-mRNA splicing. Splic...

