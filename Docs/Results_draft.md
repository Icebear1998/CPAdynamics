# CPAdynamics — Results (draft v1)

Oct 15, 2026 · @Kien

## Editor's note

This is a bullet-point draft of the Results; prose comes later. Values and figures come from `RunResultsAnalyses.m` (MATLAB R2018b, 2026-10-15) and `python/assemble_results_paper.py`; rerun both after any parameter change. One value still needs data: ⟦R5.iqr\_bp⟧.

- **Figures.** Figure 1, Figure 4A and Figure 5A are the original images from `Docs/Results.pdf`. All other figures are regenerated from the new runs in the style of the original MATLAB plots, with added detail (benchmarks, default markers, values, extra conditions). Kept from the originals: CAD vs M (M = 1–10, weak and strong PAS), assembly CDFs, M = 1/M = 5 maps (now part of Figure 2), Ser2P/AverageE profiles (now TSS to PAS), CAD vs gene length, CAD vs E\_total and usage vs inter-PAS distance. CAD vs Pol\_total went with the removed competition section.
- **Baselines.** `default_parameters.m` as of 2026-10-15: E\_total = Pol\_total = 100,000, k\_Poff = 0.01 s⁻¹, k\_Eoff\_engaged = 0.065 s⁻¹ at the weak PAS. R2–R3 use the weak PAS (k\_Hoff = 1 s⁻¹), R4–R5 a strong PAS (k\_Hoff = 0.0025 s⁻¹, with k\_Eoff\_engaged scaled by k\_Hoff/(k\_Hoff + k\_Hon) to 1.9 × 10⁻⁴ s⁻¹).
- **Benchmark criterion.** CAD₅₀ ≤ 600 bp at the weak PAS. The 200–600 bp band is drawn for reference only; commitment faster than 200 bp is not a failure.
- **Robustness.** Parameter-range maps for R2 (Figure 2D–E), range-end sweeps for R4–R5 (Figures S1–S2) and a k\_Eoff sweep for R5, all over the current `P.ranges`. Sweeps that vary k\_Hoff or k\_Eoff are shown two ways: k\_Eoff\_engaged held at its default (main figures), and k\_Eoff\_engaged scaled as k\_Eoff · k\_Hoff/(k\_Hoff + k\_Hon) (Figure S3; dashed lines in Figure 5D–E). Claims are paired comparisons or trend shapes; absolute values are at baseline.

Open items:

- [ ] Choose which k\_Eoff\_engaged treatment goes in the main figures (held fixed, or scaled with k\_Eoff and k\_Hoff). Both are shown for now.
- [ ] One site meets the 600 bp criterion at 1.3× the measured E or at k\_c ≥ 0.185 s⁻¹, both inside the appendix ranges, so the absolute R2 claim depends on E\_total and k\_c.
- [ ] k\_Hoff = 0.0025 s⁻¹ lies outside the appendix range (0.1–5 s⁻¹).
- [ ] Pol\_total is now 100,000; Appendix G still says about 70,000.
- [ ] Genome-wide free Pol II is only 82 molecules (R4); check the per-gene k\_in scaling.
- [ ] Add a citation for distal shifts after CPA factor knockdown (R5).
- [ ] Add the inter-PAS spread (R5), the cross-tissue abundance data (R1) and R6.

## Early CPA factors are about as abundant as Pol II

- Across 7,314 HeLa proteomics runs, CPSF/CstF subunits and Pol II subunits sit in the same abundance band (Figure 1a). Median abundances span 10^0.34–10^1.63 molar ppm for CPSF/CstF and 10^0.96–10^1.92 for Pol II.
- The CPSF subunits that limit intact complexes are the least abundant: CPSF1/CPSF160 at 10^0.60 and WDR33 at 10^0.34 ppm. That is within about 2–4-fold of RPB1 (10^0.96), the Pol II subunit that carries the CTD.
- Both sets lie mid-proteome: CPSF/CstF rank 1,582–4,122 and Pol II rank 1,141–2,966 of 4,393 proteins (Figure 1c). Detection is high for most subunits (Figure 1b), so the comparison is not driven by missing values.
- In the model this becomes E\_total = Pol\_total = 100,000 (E/Pol II = 1.0; 0.36–8 across the appendix bounds). Free E then sits at about 0.44 µM, half the effective per-site K\_d of 0.89 µM.
- Open question carried forward: more E would speed CPA (R2), so why is E not more abundant? R5 offers one answer.
- A cross-tissue dataset with the same ratio will be added.

![Figure 1: abundance of 3' end processing factors and Pol II across 7,314 HeLa runs](../Results/AssembleResultsPaper/Fig1_original.png)

**Figure 1.** Original figure. (a) Abundance (log10 molar ppm) of CPSF, CstF and Pol II subunits across 7,314 HeLa runs; numbers, medians. (b) Share of runs in which each subunit is quantified; dashed, 90%. (c) Whole-proteome abundance rank, with CPSF/CstF and Pol II subunits marked.

## One E site is not enough at the measured E level; five sites are

- **Constraint.** E is about as abundant as Pol II (R1), so E\_total is held at the measured 100,000 (E/Pol II = 1.0). The question is whether one E site per polymerase can reach the measured assembly distance under that limit, or whether several sites are needed.
- **Benchmark.** In the cis-antisense rescue assay, a downstream antisense copy of the PAS blocks processing unless commitment happens first, so the commitment CDF is the rescue curve (Chao et al., 1999). For the weak SV40 early PAS with a long antisense, rescue reaches 25% at about 435 bp and 50% at about 550 bp (about 480–500 bp after normalizing to its \~75% plateau), and is complete near 800 bp. The 400–800 bp span is the time from first protection to completion, not a CAD₅₀ range, so the criterion is CAD₅₀ ≤ 600 bp.
- **One site misses the criterion; a few sites meet it (Figure 2A–B).** At the weak PAS (k\_Hoff = 1 s⁻¹), CAD₅₀ is 720 bp with one E site. It falls to 472, 387 and 318 bp with two, three and five sites, and keeps falling slowly to 275 bp at M = 10. By 550 bp, 41% of single-site polymerases have committed, against 71% with five sites and 50% in Chao's data.
- **Strong PAS.** Same shape, about 25–35% shorter: 448, 246 and 224 bp for M = 1, 5 and 10. Chao et al. measured commitment only at the weak site, so the strong-PAS values are not benchmarked.
- **The shortfall is set by E abundance (Figure 2C).** One site reaches 600 bp only at 1.3× the measured E (595 bp at 130,000 molecules); five sites need 0.41×. At the measured E, one site misses the criterion and five sites meet it with room to spare.
- **Across the parameter ranges (Figure 2D–E).** Over the k\_Hd–k\_Ed, k\_c–k\_Hd and k\_c–k\_Ed planes (10 × 10 log-spaced points over the appendix ranges), CAD₅₀ is at most 600 bp in 30%, 30% and 26% of each plane with one site, against 82%, 57% and 54% with five; over all points, 29% against 64%. At baseline (star) one site lies just outside the 600 bp contour and five sites well inside it. The remaining failures with five sites sit at slow k\_c, where commitment, not E binding, is limiting, and at the weakest PAS end of the k\_Hd axis (k\_Hoff = 5 s⁻¹).
- **With k\_Eoff\_engaged scaled (Figure S3).** When k\_Eoff\_engaged moves with k\_Eoff and k\_Hoff along the map axes, the shares become 39%, 34% and 25% for one site and 77%, 58% and 50% for five; over all points, 33% against 62%. The two-fold advantage of five sites does not depend on this choice.
- **Where one site can pass.** One site fails across the whole k\_Poff range (679–771 bp for 0.0016–0.02 s⁻¹; five sites 311–341 bp). It passes with ≥ 1.3× the measured E or with k\_c ≥ 0.185 s⁻¹, both inside the appendix ranges. The paired gain holds everywhere; the absolute claim needs E\_total and k\_c near their measured values.

![Figure 2: one vs five E sites at the measured E: CAD vs M, assembly CDFs, CAD vs E abundance and CAD maps](../Results/AssembleResultsPaper/Fig2.png)

**Figure 2.** Weak PAS unless noted. (A) CAD₅₀ vs the maximum number of E sites per polymerase (M = 1–10) at the weak and strong PAS; band, 200–600 bp reference window. (B) CPA assembly completion (commitment CDF) vs Poly(A)–antiPoly(A) distance for M = 1, 2, 3 and 5; squares, Chao et al. 25% and 50% rescue; dashed drop lines, each curve's CAD₅₀; dash-dot line, rescue complete. (C) CAD₅₀ vs E\_total at Pol\_total = 100,000 for M = 1 and M = 5; red circles, measured E; band, 200–600 bp. (D1–D3, E1–E3) CAD₅₀ over k\_Hd–k\_Ed, k\_c–k\_Hd and k\_c–k\_Ed (10 × 10 log-spaced points over the appendix ranges) for M = 1 (D) and M = 5 (E), with k\_Eoff\_engaged held at its default. Red dashed, 200 and 600 bp contours; star, baseline; panel titles, share of each plane with CAD₅₀ ≤ 600 bp.

## Multiple sites work by concentrating E right before the PAS

- **Pre-loading (Figure 3A).** Just upstream of the PAS, 26% of polymerases carry at least one E with one site, against 62% with five; the share keeps rising slowly to 73% at M = 10.
- **E is concentrated right before the PAS (Figure 3B–C).** The five-site model holds less E near the TSS but more near the PAS. In the first 2 kb of loading, five-site polymerases carry 0.69–0.78× the E of single-site ones; the ratio crosses 1 about 15 kb before the PAS and reaches 3.4× just before it. With five sites 45% of the bound E on the gene body sits in the last 5 kb, against 34% with one. The same limited pool is spent where commitment happens.
- **Why less E near the TSS.** Near the TSS few polymerases carry more than one Ser2P site, so both models bind E the same way and bound E tracks free E. Free E is lower with five sites (51% against 79% of the total) because more E is held on polymerases near the PAS.
- **Profiles (Figure 3B).** Bound E and Ser2P sites rise from the TSS to the PAS. At the PAS, five-site polymerases carry 4.4 sites and 0.89 E, against 0.9 sites and 0.26 E with one site. Over the first \~8 kb, M = 3 and M = 5 carry the same E because few polymerases yet have more than three sites; the curves separate only nearer the PAS.

![Figure 3: E pre-loading at the PAS, profiles from TSS to PAS, bound E relative to one site](../Results/AssembleResultsPaper/Fig3.png)

**Figure 3.** Weak PAS. (A) Polymerases carrying at least one E at the node before the PAS, vs M. (B) Mean bound E (solid, AverageE) and Ser2P sites (dashed) from the TSS to the PAS for M = 1, 3 and 5. (C) Bound E per polymerase for M = 3 and M = 5 relative to M = 1 along the gene; grey, less E than with one site; dotted, where the ratio crosses 1.

## Gene length sets how much E polymerases bring to the PAS

- **Loading runway (Figure 4A).** On a short gene, Pol II reaches the PAS before Ser2P and E have accumulated, so commitment must finish farther downstream. A long gene gives the CTD a runway to load E before the PAS.
- **Genome-wide setup.** 10,000 genes share one pool of free Pol II and free E, with TSS-to-PAS lengths from the human log-normal distribution (median 22 kb); strong PAS, M = 5. At the shared steady state, free Pol II falls to 82 molecules while 39% of E stays free.
- **CAD vs length (Figure 4B).** CAD₅₀ is 3,022 bp at 2.5 kb, 1,256 bp at 5 kb and 271 bp at 20 kb, and levels off near 254 bp from about 20 kb. The median gene (22 kb) sits at the start of the plateau, so length matters for the shorter half: 47% of genes are under 20 kb and 27% under 10 kb, and CAD₅₀ drops below 600 bp only from about 8 kb.
- **What polymerases carry at the PAS (Figure 4D).** On a 5 kb gene they arrive with 0.74 Ser2P sites and 0.12 E; on a 20 kb gene with 4.40 sites and 0.72 E, about 6-fold more E.
- **Multivalency steepens the length dependence (Figure 4C).** At the same shared pools, five sites and one site start alike on short genes (3,022 against 3,117 bp at 2.5 kb), but five sites gain more from a longer runway. From 5 to 40 kb, CAD₅₀ falls 4.8-fold with five sites (1,256 to 259 bp) and 2.2-fold with one (1,492 to 684 bp); one site levels off near 661 bp, five sites near 254 bp. A long runway builds several Ser2P sites per polymerase, and only a multivalent CTD can fill them with E.
- **Robustness (Figure S1).** CAD₅₀ falls with length at every range end tested; the Ser2P slope moves the runway length most, and k\_c moves the plateau (Supplemental materials).

![Figure 4A: short gene vs long gene loading runway (original)](../Results/AssembleResultsPaper/Fig4A_original.png)

![Figure 4B-D: genome-wide CAD vs gene length, M = 5 vs M = 1, profiles from TSS to PAS for a short and a long gene](../Results/AssembleResultsPaper/Fig4.png)

**Figure 4.** (A) Original scheme: short genes lack a runway for CPA factor loading. (B–D) Strong PAS, 10,000 genes sharing free Pol II and E; M = 5 unless noted. (B) Genome-wide CAD₅₀ vs TSS-to-PAS distance under the human gene-length distribution (top strip); dotted, median 22 kb; dashed, plateau. (C) CAD₅₀ vs length for M = 5 and M = 1 at the same shared pools; circles, 5 and 40 kb. (D) Ser2P sites (green) and mean bound E (blue) from the TSS to the PAS for a 5 kb gene (dashed) and a 20 kb gene (solid); circles and box, values just before the PAS.

## E sits just above its most responsive level, so APA reacts strongly to losing E

- **Human PAS spacing (Figure 5A).** Distances between adjacent poly(A) sites in the 3′-most exon are roughly log-normal and peak at a few hundred nucleotides; the most common spacing is about 300 nt (spread ⟦R5.iqr\_bp⟧, to be quantified from the dataset).
- **Usage from the CDF.** A polymerase uses the proximal PAS only if it commits before reaching the distal one, so proximal usage equals the commitment CDF at the inter-PAS distance. At 300 bp and the strong PAS, usage is 58%: intermediate, not all-or-none.
- **PAS strength (Figure 5B).** Usage at 300 bp is 58% for k\_Hoff = 0.0025 and 0.05 s⁻¹, 56% for 0.2, 48% for 1, 29% for 5 and 11% for 25. Across the appendix range (0.1–5 s⁻¹) strength moves usage from about 56% to 29%; only k\_Hoff = 25 s⁻¹ lies outside it (dashed). Sites inside the range approach full usage once the spacing exceeds a few kilobases.
- **E abundance (Figure 5C).** At the measured ratio of 1.0 E per Pol II, usage at 300 bp is 58%, 10 points below the 68% ceiling set by k\_c. Halving E lowers usage by 14 points; doubling E raises it by only 7. Beyond 4-fold more E, the change falls to 0.4 points per doubling.
- **Why E may sit where it does.** At the measured level, more E buys little, but losing E shifts usage strongly toward the distal site. This fits the distal shifts reported after CPA factor knockdown (citation to add) and offers one answer to the question raised in R1: E is held where APA still responds to falling CPA factor levels.
- **Why E does not bind more tightly (Figure 5C–E).** Tighter binding also speeds CPA, but it erodes APA control. Tightening k\_Eoff 500-fold (0.5 to 0.001 s⁻¹) raises usage at 300 bp from 58% to 65%, near the ceiling, and sequesters E on polymerases: free E falls from 51% to 3%. The change per doubling of E falls from 11 to 6.4 points. The usage gap between a strong and a weak PAS (k\_Hoff = 0.0025 against 1 s⁻¹) narrows from 10 to 6 points, or to 4 with k\_Eoff\_engaged scaled (dashed). With one site the E response falls further, to 2.7 points per doubling.
- **Within the measured range.** Across the measured k\_Eoff range (0.1–5 s⁻¹) the five-site E response barely changes (10.5–11.6 points per doubling). Weaker binding widens the strong–weak PAS gap (20 points at 5 s⁻¹; 29 with k\_Eoff\_engaged scaled) but lengthens CAD₅₀ from 246 to 307 bp.
- **One tight site instead of five.** A single site tightened 25-fold to match the five-site CAD₅₀ (k\_Eoff = 0.020 s⁻¹) responds to E as strongly as five sites (13 against 11 points per doubling). Its PAS discrimination depends on how k\_Eoff\_engaged is treated: 15 points with it held fixed (better than five sites' 10), 8 points with it scaled (worse). APA therefore argues against tighter binding within M = 5, but does not by itself favour five weak sites over one tight site.
- **Robustness (Figure S2).** At every parameter range end but one, the E response at the measured ratio stays between 8 and 15 points per doubling. The most responsive E level lies between 0.35 and 1.0 times Pol II, at or below the measured ratio (Supplemental materials).

![Figure 5A: distance between adjacent poly(A) sites in the 3'-most exon (original)](../Results/AssembleResultsPaper/Fig5A_original.png)

![Figure 5B-E: proximal PAS usage vs spacing, E abundance and E binding affinity, with k\_Eoff\_engaged held or scaled](../Results/AssembleResultsPaper/Fig5.png)

**Figure 5.** (A) Original histogram of distances between adjacent poly(A) sites in the 3′-most exon (log scale), with a log-normal fit. (B–E) Strong PAS, M = 5 unless noted. (B) Proximal site usage vs inter-PAS distance for six k\_Hoff values; dashed, above the appendix range (0.1–5 s⁻¹); dotted, 300 bp. (C) Usage at 300 bp vs E\_total/Pol\_total at the baseline k\_Eoff and two tighter values; red circle, measured E; dashed, ceiling. (D) Change in usage at 300 bp per doubling of E, at the measured ratio, vs k\_Eoff for M = 5 and M = 1; solid, k\_Eoff\_engaged held at its default; dashed, k\_Eoff\_engaged scaled with k\_Eoff; white and grey diamonds, one site tightened to match the five-site CAD₅₀ under each treatment; red circle, baseline; blue band, measured k\_Eoff range. (E) Usage at 300 bp at a strong PAS (k\_Hoff = 0.0025 s⁻¹) minus a weak PAS (k\_Hoff = 1 s⁻¹) vs k\_Eoff; symbols as in D.

## Genome-wide findings

To be written. This section will test the model's predictions against genomic data: gene-length dependence of CPA (R4) and APA response to CPA factor levels (R5).

## Supplemental materials

### Gene-length dependence at the parameter range ends

- **What was varied.** Five-site CAD₅₀ vs TSS-to-PAS length at the shared genome-wide pools (strong PAS), with one runway parameter at a time set to its lower or upper appendix bound: the Ser2P slope (k\_Pon slope), k\_Poff, the elongation rate k\_e and the commitment rate k\_c.
- **Runway length.** The Ser2P slope moves where the curve levels off most. The steepest slope (2 × 10⁻³) brings the plateau in from about 20 kb to 8 kb; the shallowest (10⁻⁴) pushes it out to about 63 kb, and genes shorter than about 8 kb then do not reach 50% within 5 kb. Across its range (0.0016–0.02 s⁻¹), k\_Poff moves the plateau only between 16 and 32 kb.
- **Plateau level.** k\_c sets the level the curve settles to: about 116 bp at 0.6 s⁻¹ and 1.5 kb at 0.015 s⁻¹.
- **Elongation.** k\_e shifts short genes (799–1,591 bp at 5 kb against 1,256 at baseline) and leaves the plateau near 237–263 bp.
- **Conclusion.** CAD₅₀ falls with gene length at every range end; the parameters change where and how far it falls, not the direction.

![Figure S1: five-site CAD vs gene length at the range ends of the runway parameters](../Results/AssembleResultsPaper/FigS1.png)

**Figure S1.** Strong PAS, M = 5, shared genome-wide pools. CAD₅₀ vs TSS-to-PAS distance at baseline (teal) and with (A) the Ser2P slope, (B) k\_Poff, (C) k\_e and (D) k\_c at its lower (dotted) or upper (dashed) appendix bound; values in the legends. Missing points did not reach 50% within 5 kb.

### Sensitivity of APA to E at the parameter range ends

- **What was varied.** Proximal usage at 300 bp vs E\_total/Pol\_total (strong PAS, M = 5), with each of nine parameters (k\_Eon, k\_Eoff, k\_c, Ser2P slope, k\_Poff, Pol\_total, k\_e, k\_e2, k\_Hon) at its lower or upper appendix bound.
- **Result.** At every range end but one, the change per doubling of E at the measured ratio stays between 8 and 15 points. The exception is the slowest commitment rate (k\_c = 0.015 s⁻¹), where usage at 300 bp is only 11% and the change is 3 points.
- **Where the response peaks.** The E level of maximal response lies between 0.35 and 1.0 times Pol II at every range end, so the measured ratio is at or above the most responsive level throughout.

![Figure S2: change in proximal usage per doubling of E at the parameter range ends](../Results/AssembleResultsPaper/FigS2.png)

**Figure S2.** Strong PAS, M = 5. Change in usage at 300 bp per doubling of E vs E\_total/Pol\_total at baseline (teal) and with each of nine parameters at its lower or upper range end (grey); red dotted, measured ratio.

### CAD maps with k\_Eoff\_engaged scaled

- **What changes.** The maps of Figure 2D–E, but with k\_Eoff\_engaged recomputed at each point from its tethered-rebinding derivation, k\_Eoff · k\_Hoff/(k\_Hoff + k\_Hon), scaled from 0.065 s⁻¹ at baseline. In Figure 2 it is held at 0.065 s⁻¹ everywhere.
- **Result.** CAD₅₀ is at most 600 bp in 39%, 34% and 25% of the k\_Hd–k\_Ed, k\_c–k\_Hd and k\_c–k\_Ed planes with one site, against 77%, 58% and 50% with five; over all points, 33% against 62% (29% against 64% with it held).
- **Where it matters.** Scaling slows engaged-E loss at tight binding and strong PAS, which helps one site where it is close to the criterion. It speeds that loss at weak binding and weak PAS, which costs five sites a little at that corner.
- **Conclusion.** Five sites meet the criterion about twice as often as one under either treatment.

![Figure S3: CAD maps for one and five E sites with k\_Eoff\_engaged scaled with k\_Eoff and k\_Hoff](../Results/AssembleResultsPaper/FigS3.png)

**Figure S3.** Weak PAS baseline. As Figure 2D–E, with k\_Eoff\_engaged scaled with k\_Eoff and k\_Hoff at each map point, for M = 1 (A1–A3) and M = 5 (B1–B3). Red dashed, 200 and 600 bp contours; star, baseline; panel titles, share of each plane with CAD₅₀ ≤ 600 bp.
