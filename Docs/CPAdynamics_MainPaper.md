# **Introduction**

   
Biomolecular machineries must assemble efficiently to function in a highly heterogeneous physicochemical environment that is shaped by widely varying protein-protein and protein-RNA binding kinetics, molecular abundances, crowding-limited diffusion, and molecule-specific electrical charges (). 

These conditions are especially demanding for machineries that process nascent RNAs during transcription because these machineries assemble from many protein factors that bind onto transcribing RNA polymerase II (Pol II). 

The cleavage and polyadenylation (CPA) machinery, which cleaves nascent RNAs, adds the poly(A) tail to nascent RNAs, and is tightly coupled with transcription termination is a key RNA-processing machinery that must assemble efficiently. 

The CPA machinery assembly must complete within seconds to tens-of-seconds *in vivo* (Chao et al., 1999). This signifies that CPA assembly must be fast because transcription of a human gene of typical length (\~25 kilobases) takes minutes to tens-of-minutes given typical transcription elongation rates (2 to 6 kilobases/min) (Darzacq et al., 2007; Fuchs et al., 2014). This speed is striking because assembly is governed by a multitude of protein-protein interactions among CPA factors and between CPA factors and the Pol II carboxyl-terminal domain (CTD). Many protein-RNA interactions between RNA-recognizing CPA factors and 3’-end sequence elements of the nascent transcript, which includes the polyadenylation signal (PAS) and other auxiliary upstream and downstream sequence elements, also prevail (). 

This speed is particularly difficult to explain because the CPA factors involved in PAS recognition (CPSF160, CPSF30, WDR33, and Fip1) are comparable in abundance to or lower in abundance than actively transcribing Pol II (See Result 1, Bekker-Jensen et al., 2017; Jackson et al., 1998; Kimura et al., 1999). 

Consequently, it is likely that intense inter-Pol II competition for the limited pool of CPA factors exists within the nuclear environment making rapid CPA assembly and function unlikely.

Experiments also suggest that rapid CPA assembly is difficult to achieve: perturbing CPA-factor abundance has been shown to alter processing outcomes and shift poly(A) site usage genome-wide (Fusby et al., 2016; Grzechnik & Mischo, 2025; Jin et al., 2025; Kowalski et al., 2024; Mitschka & Mayr, 2022; Tang & Zhou, 2022).  Additionally, limiting CPA-factor availability likely reduces the probability of productive CPA assembly on nascent RNAs, thereby increasing the likelihood of delayed or unsuccessful cleavage and polyadenylation (Casañal et al., 2017; Davidson et al., 2014; Kowalski et al., 2024; Neve et al., 2017). Together, these observations indicate that CPA machinery assembly is sensitive to the local availability of factors within the nucleoplasmic pool and that rapid assembly will be harder still given intense inter-Pol II competition for CPA factors.

These findings raise the question of whether cells can assemble the CPA machinery reliably within seconds to tens-of-seconds despite the physically constraining conditions of the nucleoplasm, limited factor availability, and competition among elongating Pol II.

One way to speed assembly lies in the Pol II CTD, which can serve as a platform that carries CPA factors toward the PAS. Phosphorylation of CTD serine 2 (Ser2P) accumulates as Pol II elongates, and PAS function further enriches Ser2P near the gene end (Davidson et al., 2014). If Ser2P-marked repeats can each hold a CPA factor, a polymerase could arrive at the PAS already carrying several, so assembly would not have to begin from the free pool at the PAS. This preloading has a cost, because a factor held by one polymerase is unavailable to the others, which ties assembly at a single gene to the genome-wide competition described above. The same reasoning raises two further questions. Because phosphorylation accumulates during elongation, gene length should limit how much loading is possible before the PAS. And because commitment takes a finite distance, the spacing between alternative PASs should influence which site is used.

<!-- TODO: add citations for (i) CPA factors binding Ser2P-marked CTD repeats and (ii) existing quantitative models of CPA assembly/termination (gap statement). -->

In this study, we develop a model that quantifies how quickly CPA assembly completes and provides a predictive framework for how that speed depends on factor supply, gene length, and PAS strength. Since cleavage, polyadenylation, and termination each depend on the completion of the CPA machinery assembly, understanding the assembly dynamics is essential to understanding both 3′\-end processing and transcription termination. Specifically, we built a finite-rate kinetic model in which elongating Pol II acquires phosphorylated CTD binding sites, loads PAS-recognizing CPA factors from a shared nuclear pool, and commits to CPA downstream of the poly(A) signal (Results 2; Methods).

We quantify assembly speed with the CPA assembly distance (CAD₅₀), the distance downstream of the PAS by which half of the polymerases have committed to CPA. At the reference parameters, allowing each Pol II to carry up to five E factors rather than one shortened CAD₅₀ from 783 to 331 bp, by loading E onto Pol II before it reaches the PAS (Results 2–3). This loading draws on a shared pool, so CAD₅₀ falls as total E increases and rises as total Pol II increases (Results 3). Because loading builds during elongation, CAD₅₀ also depends strongly on gene length: for a strong PAS it fell from approximately 3,460 bp at a transcription start site (TSS)-to-PAS distance of 2.5 kb to approximately 288 bp at 200 kb, and about 31% of human protein-coding genes are shorter than the roughly 10.4 kb at which CAD₅₀ reaches 500 bp (Results 4). The same kinetics predict that proximal poly(A) site usage rises with the distance between sites and falls with weaker PAS recognition (Results 5).

# **Results**

# **1\. CPA factor abundance motivates a limited-resource model of assembly**

To establish the context for CPA-factor recruitment, we compared the abundances of selected CPA factors and RNA polymerase II subunits (Figure 1). The displayed abundance profiles varied substantially across CPA components, with several components less abundant than the Pol II subunits shown. These comparisons motivated a model in which the supply of factors available for CPA assembly is finite.

We represented this supply as a shared pool of effective assembly factors, denoted E, that can bind to elongating Pol II. E represents the PAS-recognizing CPA factors (CPSF160, CPSF30, WDR33, and Fip1).

 *\- Need to either include a conversion for Pol II core subunit (RPB1) based on the % of actively transcribing Pol II within the cell (based on the total Pol II pool) AND/OR utilize BIONUMBER’s actively transcribing Pol II quantity. See Discussion Section.*

![CPA factor and Pol II abundance][image2]

**Figure 1\.** CPA-factor and Pol II abundance comparison.

# **2\. Multiple E-binding sites accelerate CPA commitment**

Our kinetic model tests whether multiple binding sites on the Pol II carboxyl-terminal domain (CTD) can support efficient CPA assembly under finite factor availability (Figure 2A). As an elongating Pol II (state R) transcribes the gene, Ser2 phosphorylation (Ser2P) forms E-binding sites on its CTD, and each site can hold one E factor taken from the shared free pool. The CTD capacity M limits both the number of sites and the number of E factors a polymerase can carry. Downstream of the PAS, a polymerase carrying at least one E can form a reversible PAS-recognition complex (state REH), in which one E is engaged with the PAS hexamer (H); REH either dissociates back to R or commits irreversibly to CPA at rate kc. We quantified assembly efficiency using CAD₅₀, the distance traveled downstream of the PAS by which 50% of the modeled population has committed. All symbols are defined in Table M1 and reference values are listed in Table M2 (Methods).

We judged assembly against a working benchmark of CAD₅₀ ≤ 500 bp, a commitment distance of a few hundred base pairs motivated by the cis-antisense rescue experiments of Chao et al. (1999). The benchmark is a rough guide, not a measured threshold.

At the reference parameter set (Table M2), which corresponds to a weak PAS (kHoff \= 1 s⁻¹), increasing capacity from one to five sites reduced CAD₅₀ from 783 to 331 bp, a reduction of approximately 58% that moved CAD₅₀ from outside to inside the benchmark (Figure 2B, left). Further increases in capacity produced progressively smaller reductions in CAD₅₀, so we used M \= 5 as the working capacity (Results 3). The M \= 5 commitment curve rose more rapidly over the first 1,000 bp downstream of the PAS than the M \= 1 curve (Figure 2B, right).

We focused the paired maps on E binding (kEd \= kEoff/kEon), PAS recognition (kHd \= kHoff/kHon), and irreversible commitment (kc), three controls that strongly affect CAD₅₀ (Figure 2C). The maps compare kHd–kEd, kc–kHd, and kc–kEd, with association rates held fixed when varying the binding ratios. M \= 5 yielded shorter CAD₅₀ than M \= 1 wherever both conditions reached 50% commitment. Increasing PAS dissociation prolonged commitment distance, whereas increasing kc shortened it. Greater E dissociation also prolonged CAD₅₀ across the E/PAS-binding map. 

To assess the joint three-dimensional parameter space, we scored 10,000 random combinations of kHd, kEd, and kc against the benchmark, evaluating the same combinations for both capacities (Methods). One site met the benchmark in 23.3% of sampled conditions, compared with 57.0% for five sites. Thus, increasing capacity from one to five sites expanded the fraction of the sampled parameter space supporting commitment within 500 bp by approximately 34 percentage points.

![][image3]

**Figure 2A.** Model schematic linking E loading, reversible PAS recognition, and irreversible CPA commitment. Cartoon counts and distances are illustrative.

![Binding capacity and commitment curves][image4]

**Figure 2B.** Left, CAD₅₀ versus binding capacity M. Right, commitment curves for M \= 1 and M \= 5. Both use the reference parameter set. (Note: we can overlay Chao's data later)

![Pairwise parameter maps for one and five binding sites][image5]

**Figure 2C.** Paired parameter maps for M \= 1 and M \= 5\. Stars mark the baseline; dashed lines bound the 200–500 bp reference band. Benchmark attainment is scored as CAD₅₀ ≤ 500 bp over all 100 log-spaced conditions in each map, including conditions below 200 bp; an unreached 50% commitment threshold counts as not meeting the benchmark. Displayed colors are capped at 2,000 bp; grey regions and open circles mark an unreached threshold.

# **3\. Multiple binding sites concentrate E loading near the PAS while drawing from a shared factor pool**

To examine how multiple sites accelerate commitment, we followed site formation and E occupancy as Pol II approached the PAS (Figure 3). Under the assumed linear increase of the site-formation rate along the gene (Methods), mean site number and E occupancy increased along the gene. Compared with M \= 1, M \= 5 polymerases carried fewer E factors near the TSS but more immediately before the PAS, concentrating loading toward the recognition region. Immediately before the PAS, M \= 5 polymerases carried an average of 0.89 E factors across 4.44 formed sites, compared with 0.26 E factors across 0.93 sites at M \= 1\.

Increasing capacity from one to five sites raised the fraction of polymerases carrying at least one E before PAS recognition from approximately 26% to 62%. This fraction measures how many polymerases were preloaded; mean occupancy also counts additional factors on the same polymerase. Increasing M from 5 to 10 raised the preloaded fraction to 73% and reduced CAD₅₀ from 331 to 282 bp, showing progressively smaller gains from additional sites. We therefore used M \= 5 as the working capacity.

This local loading draws E factors from a finite shared factor pool (Figure 4A). E retained on one transcribing polymerase reduces the free E available for loading onto others across the gene population. Additional sites therefore increase recruitment capacity while increasing demand on the common pool; their occupancy depends on the E remaining available. 

The pool-size sweeps (Methods) show opposing effects of factor supply and polymerase abundance (Figure 4B). At M=5, increasing total E from 50,000 to 200,000 shortened CAD₅₀ from approximately 550 to 260 bp, with smaller reductions at higher supply. Increasing total Pol II from 40,000 to 150,000 increased CAD₅₀ approximately linearly, from 335 to 411 bp, at M \= 5\. These trends are consistent with commitment depending on the balance between factor supply and polymerase demand.

![E preloading and Ser2P profiles][image6]

**Figure 3\.** Fraction of polymerases carrying E immediately before the PAS and mean E and Ser2P-site profiles for selected binding capacities. Profiles exclude the PAS-recognized region.

![][image7]

**Figure 4A.** Conceptual link between local E loading and shared factor availability.

![Opposing effects of total E and total Pol II on commitment distance][image8]

**Figure 4B.** Pool-size sweeps: increasing total E decreases CAD₅₀ (left), whereas increasing total Pol II increases CAD₅₀ at M \= 5 (right). Red circles indicate the default condition. 

# **4\. Gene length limits the opportunity for E preloading**

Because E loading develops during elongation, we next examined whether the distance from the transcription start site (TSS) to the PAS affects commitment (Figure 5A). The simulated genes share the free E and Pol II pools left by a population of 10,000 genes at M \= 5, held fixed for every gene length and both capacities so that each comparison has the same factor supply. We used the strong-PAS value kHoff \= 0.01 s⁻¹, which represents the canonical AAUAAA hexamer, the most frequent signal among human poly(A) sites (about 53–58%; Beaudoing et al., 2000; Tian et al., 2005), rather than the weak-PAS reference value (Methods). Short genes reached the PAS with fewer formed sites and less bound E. This length dependence follows from the assumed linear increase of the site-formation rate along the gene: because sites accumulate with distance from the TSS, a shorter gene necessarily allows less loading. It should be read as a consequence of that assumption rather than an independent test of it, and its size depends on the ramp slope, which Ser2P enrichment profiles constrain only qualitatively (Appendix F).

<!-- TODO: consider a k_Pon,slope sensitivity sweep for the length dependence. -->


For a strong PAS with kHoff \= 0.01 s⁻¹ and M \= 5, CAD₅₀ decreased from approximately 3,460 bp at a TSS-to-PAS distance of 2.5 kb to approximately 288 bp at 200 kb (Figure 5B). CAD₅₀ reached the 500 bp benchmark at approximately 10.4 kb and changed little beyond ∼20 kb (313 bp at 20 kb versus 288 bp at 200 kb). In the same free pools, M \= 5 produced shorter commitment distances than M \= 1, with greater separation between the curves for longer genes. The accompanying profiles showed mean E occupancy immediately before the PAS rising from 0.12 on the 5 kb gene to 0.73 on the 20 kb gene.

To place this length dependence in a human genomic context, we calculated the genomic distance from the annotated transcript start to the transcript 3′ end for one MANE Select transcript per protein-coding gene on chromosomes 1–22, X, and Y (Figure 5C). The 19,229 retained genes had a median span of 22.5 kb, and approximately 31% were shorter than the ∼10.4 kb at which the modeled CAD₅₀ reaches 500 bp (approximately 16% were shorter than 5 kb).

The observed distribution therefore includes a substantial minority of genes in the length range where the model predicts that commitment depends strongly on preloading opportunity, while most genes are long enough to meet the benchmark.

We use these genomic spans, including introns, as a proxy for TSS-to-PAS distance. They end at the annotated transcript 3′ end rather than the upstream PAS hexamer and do not resolve alternative PAS usage (Methods).

![][image9]

![Gene length and assembly commitment][image10]

**Figure 5A–B.** Short- and long-gene schematic, CAD₅₀ versus TSS-to-PAS distance, and loading profiles. Green cartoon dots represent formed sites; grey positions have not yet formed available sites. Both capacities use the shared pools solved for M \= 5\. The 22 kb line marks the assumed population median.

![Human transcript spans aligned with modeled commitment distance][image11]

**Figure 5C.** Human transcript spans compared with the modeled length dependence of commitment. Upper panel: the existing M \= 5, kHoff \= 0.01 s⁻¹ curve, shown only over its simulated range of 2.5–200 kb. Lower panel: the full distribution of genomic TSS-to-transcript-3′-end distances for 19,229 MANE Select protein-coding genes on assembled chromosomes, with one transcript per gene. Both panels share a logarithmic distance axis; shading marks the simulated range, and the dashed line marks the observed median of 22.5 kb. Approximately 85.8% of spans fall within the simulated range; all shorter and longer spans remain in the histogram. Distances are coordinate differences, including introns, verified against the MANE mRNA annotation.

# **5\. Commitment kinetics provide a mechanism for distance-dependent APA regulation**

We next considered how the same commitment kinetics could influence the choice between proximal and distal polyadenylation sites. The interval between sites defines the distance over which proximal commitment can occur before the distal PAS is transcribed (Figure 6A). Assuming that reaching the distal site closes this opportunity, the cumulative commitment probability at the inter-PAS distance estimates proximal usage.

To represent differences in PAS recognition, we varied the effective dissociation rate kHoff from 0.001 to 10 s⁻¹ while holding kHon fixed (Methods); larger kHoff therefore corresponds to weaker effective binding.

Proximal commitment increased with inter-PAS distance, while larger kHoff values shifted commitment toward greater distances (Figure 6B). The curves for kHoff \= 0.001 and 0.1 s⁻¹ nearly overlapped; larger dissociation rates produced lower commitment at the same spacing. The aligned histogram shows how annotated inter-PAS spacings fall across these distances (Figure 6B \- lower).

![][image12]

**Figure 6A.** The interval between PASs provides an opportunity for proximal commitment.

![PAS strength and spacing with aligned inter-PAS distance histogram][image13]

**Figure 6B.** Proximal commitment for kHoff \= 0.001, 0.1, 1, 5, and 10 s⁻¹, with M \= 5 and a 50 kb simulated gene. Both panels share a logarithmic axis spanning 30–5,000 bp.

# **Discussion**

***Discussions derived from the RESULTS sections:***

- Chao’s results are used roughly for comparison; they should not be the absolute benchmark.


- Why we chose M \= 5: increasing M beyond 5 gains little efficiency


- The chosen parameter space is represented in a physiologically realistic manner. We utilized direct measurements and also derived all estimates based on published experimental data and literature.   
  - For parameters that we could not find published data, we carefully implemented model approximations to quantitatively estimate the unknown parameters (i.e., kc) or utilized parameter values that qualitatively simulate experimental data from which we do not have a well-defined quantitative metric (i.e., phosphorylation dynamics across gene length). 

- We need to state how the CPA factor/Pol II subunit abundance measurements translate into the numbers of E-factor versus the **number of actively transcribing Pol II** \- we can convert the abundance from Joshua into the predicted CPA numbers using the billion protein molecules per cell estimate. We also need to be clear that the abundance of Pol II core does not necessarily mean that this correlates with actively transcribing Pol II because only a fraction of Pol II will be transcribing within the cell.

  - [https\://doi.org/10.1083/jcb.200206019](https://doi.org/10.1083/jcb.200206019) \- Kimura et al. (2002) \- ∼20 (HeLa) to 25% (mammalian) of the total Pol II pool are actively transcribing or transiently immobilized (but engaged with DNA template \- aka paused) in human genes

  - [10.1091/mbc.9.6.1523](https://doi.org/10.1091/mbc.9.6.1523) \- Jackson et al. (1998) \- “each (HeLa) cell contains ∼75,000 active polymerases II/III”  

  -  [10.3389/fgene.2012.00221](https://doi.org/10.3389/fgene.2012.00221) \- Rieder et al. (2012) \~65,000 active polymerase molecules 

  - [https\://doi.org/10.1023/A:1024926710797](https://doi.org/10.1023/A:1024926710797) \- Martin & Pombo (2003) \- 65,000 to 85,000 in Hela cell. Also states that active pol II is ∼25% of total Pol II pool within HeLa cells based on studies Jackson et al. (1998) and Kimura et al. (1999).

    

- Simply raising the abundance of E-factors (PAS-recognizing CPA factors) does not necessarily allow the CPA machinery assembly to be more efficient/lower CAD within the cell.   
  - We show that even if the cell is able to raise E-factor abundance… , the CAD is decreased by … . raising this abundance leads to values that are likely to be physiologically unreasonable (that is, unattainable E-factor production by the cell.)

***Addressing Open Questions Section:***

- Why do certain CPA factors spatiotemporally associate with the Pol II CTD during certain portions of transcription?  
  -  Even though there is significant experimental evidence showing that E-factors have the tendency to associate with Pol II CTD early within transcription when compared to L-factors, there is no clear reason why L-factors cannot bind to the CTD before E-factors since both utilize Ser2P sites.   
  - How does this impact “intermediate” CPA factors?   
    - CPSF73 has the propensity to bind to the Pol II CTD intermediately within the transcription process and is vital for the cleavage function in the final assembly.

- How many Serine 2 phosphorylated repeats are required for CPA factors to effectively bind the Pol II CTD?   
    
- Is the Pol II CTD phosphorylated in a sequential or distributive manner? Is phosphorylation of the serine 2 positions processive or can these phosphorylation dynamics be described through unspecific position binding and unbinding events?

- Why do the majority/typical protein coding genes have a limited length distribution and not a larger range?  
  - Serine 2 level could potentially/partially explain the both minimum and maximum gene length limit within cells/organisms. (If a protein coding gene is too short, not enough Ser2P sites on Pol II CTD for proper CPA and post translational operations. If a protein coding gene is too long, perhaps there is a loss in APA and post translational modification fidelity \- that is, the proximal PAS will be selected as an uncontested majority by the CPA machinery.)

- Our model framework can be utilized to help answer even bigger questions like why do many organisms have a distribution of the number of chromosomes? (Could potentially be used to explain chromosome number variability of polyploids and different organisms.)

- Why does InterPAS distance vary greatly among different genes/transcripts?  
  - Median \= 1,521 nucleotides (when excluding below 30 nt)  
  - Range \[30, 10^6\] nt  
  - Regulators like RNA-binding proteins, antisense transcripts, etc. that have the potential to interact with or occlude PAS sites (especially the proximal \- in the canonical/our model’s predictions sense) that themselves can be regulated can act as another regulatory layer/mechanism that influences CPA assembly and CPA machinery “choice” for a PAS when multiple PAS are available via APA-associated transcripts.

- The commitment rate (kc) is a single lumped step in our model. We estimated it from a coarse-grained scheme in which a late-factor unit (L; CstF) binds REH in two reversible steps followed by a third, irreversible step (Appendix D). With the reference late-factor pool this gives kc = 0.14 s⁻¹ (a mean commitment time of about 7 s). A tighter late-factor affinity gives 0.44 s⁻¹ (about 2.3 s), in line with the ≲3 s commitment that Chao et al. (1999) measured at a strong PAS. This checks kc against in vivo timing; it is not an independent confirmation. In reality, more binding events occur than those involving CstF, and cleavage and polyadenylation take additional time. We fold both into kc rather than model them separately, which is an approximation: cleavage is assumed fast relative to commitment, but genome-wide cleavage half-lives are under a minute (Torres-Ulloa et al., 2024). The scheme also treats the late-factor pool as non-limiting, so kc does not respond to competition for L.
  <!-- TODO: add Torres-Ulloa, Calvo-Roitberg & Pai (2024), RNA 30:256-270, doi:10.1261/rna.079783.123 to References. -->


- Our Model makes a simplifying assumption that cleavage causes the transcribing Pol II to release from the DNA template. In reality this is not the case, but we use it here to quantify the assembly distance/time to quantify CPA assembly efficiency.  
  - Our model can be expanded to include the downstream processes that influence PAS-dependent transcription termination.

***Discussion Points on Model Limitations \+ Improving Model Results:***

- Binding kinetics of CPA factors to Pol II CTD and between CPA factors remain undetermined.   
  - These are critical values for creating more detailed/realistic future models that describe an accurate assembly process. Such values have the potential to significantly change the described CPA assembly efficiency provided by our model.

- Can Chip-chip be modified to quantify the number of repeats that are phosphorylated on transcribing Pol II in an accurate manner? This would be immensely helpful in quantifying Pol II CTD dynamics and for parameterizing this and future models.  
    
- **Does CHip-chip represent the % of serine 2 activated on the Pol II CTD at many repeats or a single repeat? How are multiple repeats counted?**  
- *From our understanding Chip-chip takes Pol II at a certain position of DNA, says if it is phosphorylated at a given CTD amino acid via a measured signal \- but does not provide the level nor number of phosphorylated CTD sites. This Ser2P signal “across the gene” is normalized via the “total Pol II at the gene positions” signal to account for Pol abundance changes in the gene due to pausing, elongation rate differences, etc. Thus, this Ser2P signal serves as a metric for a relative increase in Ser2P from the TSS to TTS.*

# **Methods**

## **A finite rate model of cotranscriptional CPA assembly**

We modeled how RNA polymerase II (Pol II) elongation, phosphorylation of the Pol II carboxyl-terminal domain (CTD), loading of cleavage and polyadenylation (CPA) factors, recognition of the poly(A) signal (PAS), and CPA commitment jointly determine the distance from the PAS to completed CPA assembly (Figure 2A). We represented these processes as a coupled system of finite-rate reactions and transport along a discretized gene.

The early, PAS-recognizing CPA factors (CPSF160, CPSF30, WDR33, and Fip1) were treated as a single effective factor, E (for early), drawn from a shared pool of free E (${E}_{f}$; the total pool is ${E}_{total}$). Phosphorylation of CTD serine 2 (Ser2P) forms E-binding sites, each a pair of adjacent Ser2P repeats, up to a capacity of $M$ sites per polymerase. An elongating polymerase carrying $p$ formed sites and $e$ bound E is in state $R(p,e)$. Downstream of the PAS, it can enter the recognition complex $REH(p,e)$, in which one bound E is engaged with H, the PAS hexamer recognition element on the nascent RNA.

CTD phosphorylation, CPA-factor loading, and PAS recognition are represented by reversible binding kinetics. 

Commitment, at rate $k_c$, is the effectively irreversible assembly of the complete, cleavage-competent CPA machinery, including recruitment of the late (L) factors such as CstF. Cleavage itself is not modeled as a separate step; it is assumed fast enough after commitment to be lumped into $k_c$. After commitment, E and Pol II are released after cleavage and return to the free pools (free E, ${E}_{f}$, and free Pol II, ${P}_{f}$; the total Pol II pool is ${P}_{total}$).

<!-- TODO: temporary wording on release of E and Pol II to the free pools after commitment/cleavage; revise. -->

We quantify CPA assembly efficiency with the CPA assembly distance (CAD), reported as ${CAD}_{50}$: the distance travelled downstream of the PAS by which 50% of the polymerases that reach the PAS have committed. At the typical transcription elongation rate, this distance converts to the time needed for CPA assembly to complete and commit.

![][image3]

**Figure 2A.** Model schematic linking E loading, reversible PAS recognition, and irreversible CPA commitment. Cartoon counts and distances are illustrative. (Repeated from Results 2.)

**Table M1.** Abbreviations and notation. Rate constants are written $k_{Hoff}$ in equations and kHoff in text and figures. Reference values are in Table M2.

| Symbol | Definition |
| :---- | :---- |
| Pol II | RNA polymerase II |
| CTD | carboxyl-terminal domain of the largest Pol II subunit |
| Ser2P | CTD repeat phosphorylated at serine 2 |
| CPA | cleavage and polyadenylation |
| PAS | poly(A) signal |
| TSS | transcription start site |
| APA | alternative polyadenylation |
| E | early (PAS-recognizing) CPA factors (CPSF160, CPSF30, WDR33, Fip1), treated as one effective factor drawn from a shared pool |
| H | PAS hexamer recognition element; in the recognition complex, one bound E is engaged with H |
| L | late CPA factors (such as CstF) whose assembly after PAS recognition is lumped into ${k}_{c}$ |
| $R(p,e)$ | elongating Pol II with $p$ formed E-binding sites and $e$ bound E, not in a recognition complex |
| $REH(p,e)$ | Pol II in the PAS recognition complex, with $p$ sites and $e$ bound E, exactly one of which is engaged with H |
| $p$ | number of formed E-binding sites; one site is a pair of adjacent Ser2P repeats |
| $e$ | number of bound E, with $0\leq e\leq p\leq M$ in R and $1\leq e\leq p\leq M$ in REH |
| $M$ | CTD capacity: maximum number of E-binding sites (and so of bound E) per polymerase |
| ${E}_{f}$, ${P}_{f}$ | free E and free Pol II; totals are ${E}_{total}$ and ${P}_{total}$ |
| $G$, ${x}_{PAS}$ | template length and PAS position, measured from the TSS |
| ${L}_{a}$, $N$, $q$, $i$ | node spacing (100 bp), number of nodes, PAS node, and node index |
| ${k}_{in}$ | Pol II initiation rate |
| ${k}_{e}$, ${k}_{e2}$ | elongation rate per node of R and of REH |
| ${k}_{Pon}(i)$, ${k}_{Pon,min}$, ${k}_{Pon,slope}$ | rate of E-binding site formation at node $i$, its value at the TSS, and its linear increase per node |
| ${k}_{Poff}$ | loss rate of a formed site not occupied by E |
| ${k}_{Eon}$, ${k}_{Eoff}$ | E association (per free E and open site) and dissociation (per bound E) rates |
| ${k}_{Eoff,engaged}$ | detachment rate of the E engaged with H in REH |
| ${k}_{Hon}$, ${k}_{Hoff}$ | PAS recognition rate (per bound E) and H dissociation rate |
| ${k}_{c}$ | CPA commitment rate from REH (see above) |
| kEd, kHd | dissociation ratios ${k}_{Eoff}/{k}_{Eon}$ and ${k}_{Hoff}/{k}_{Hon}$ |
| CAD, ${CAD}_{50}$ | CPA assembly distance; the downstream distance from the PAS at which cumulative commitment probability reaches 50% |

## **Gene discretization and polymerase microstates**

We represented each gene as a one-dimensional lattice with node spacing ${L}_{a}=100 bp$. For template length $G$ and a PAS located a distance ${x}_{PAS}$ downstream of the transcription start site (TSS), we used $N=⌊G/{L}_{a}⌋$ nodes (indexed $i=1,\dots,N$) and PAS node $q=⌊{x}_{PAS}/{L}_{a}⌋$. Recognition and commitment were permitted from node $q$ through $N$. The default single-gene template was 25 kb long, with its PAS 20 kb downstream of the TSS.

At node $i$, ${R}_{i}(p,e)$ denotes the population of polymerases that carry $p$ formed E-binding sites (pairs of phosphorylated Ser2 repeats) and $e$ bound E factors and have not formed a PAS recognition complex, with $0\leq e\leq p\leq M$. Polymerases in the recognition complex were represented by ${REH}_{i}(p,e)$, with $1\leq e\leq p\leq M$. 

In an REH state, exactly one of the $e$ factors is coupled with H. 

## **CTD phosphorylation and E factor exchange**

The rate of E-binding site formation (Ser2 phosphorylation), ${k}_{Pon}$, increased linearly along the gene from its value ${k}_{Pon,min}$ at the TSS with slope ${k}_{Pon,slope}$ per node, ${k}_{Pon}(i)={k}_{Pon,min}+{k}_{Pon,slope}(i-1)$. 

Site loss (dephosphorylation) decreased $p$ by one at rate ${k}_{Poff}$ only when $p>e$, protecting E-occupied sites.

E associated with either R or REH at $(p-e){k}_{Eon}{E}_{f}$, where ${k}_{Eon}$ is the association rate constant and $(p-e)$ is the number of unoccupied sites. Ordinary E is E that is not engaged with H in the recognition complex REH; its dissociation, with rate constant ${k}_{Eoff}$ per bound E, occurred at $e{k}_{Eoff}$ in R and $(e-1){k}_{Eoff}$ in REH. The recognition complex REH dissociates in two ways described below.

## **PAS recognition and alternative dissociation routes**

PAS recognition was modeled as the transition from an elongating Pol II state (R) to a PAS-bound REH state at the PAS. The effective PAS binding rate was obtained by scaling the basal PAS-recognition rate constant ${k}_{Hon}$ by the mean number of E factors loaded at the PAS. This coupling represents the assumption that Pol II molecules arriving with higher early-factor occupancy are more likely to enter a cleavage-competent, post-PAS complex. At and downstream of the PAS node, each bound E could initiate recognition, giving a transition $R(p,e)\rightarrow REH(p,e)$ at $e{k}_{Hon}$. 

Recognition complex REH could dissociate via two distinct routes. H dissociation, with rate constant ${k}_{Hoff}$, produced $REH(p,e)\rightarrow R(p,e)$, retaining every CTD-bound E factor. Alternatively, the engaged E detached at ${k}_{Eoff,engaged}$, returning the polymerase to $R(p,e-1)$ and releasing one E. This second route assumed rapid disassembly of the detached EH complex. Because exactly one E was engaged, ${k}_{Eoff,engaged}$ had no $e$ multiplier. Remaining bound E factors could subsequently support recognition again.

| Event | Transition and rate per polymerase |
| :---- | :---- |
| Phosphorylation | $(p,e)\rightarrow (p+1,e)$ at ${k}_{Pon}(i)$, if $p<M$; in R and REH |
| Dephosphorylation | $(p,e)\rightarrow (p-1,e)$ at ${k}_{Poff}$, if $p>e$; in R and REH |
| E association | $(p,e)\rightarrow (p,e+1)$ at $(p-e){k}_{Eon}{E}_{f}$; in R and REH |
| Ordinary E loss from R | $R(p,e)\rightarrow R(p,e-1)+{E}_{f}$ at $e{k}_{Eoff}$ |
| Ordinary E loss from REH | $REH(p,e)\rightarrow REH(p,e-1)+{E}_{f}$ at $(e-1){k}_{Eoff}$ |
| PAS recognition | $R(p,e)\rightarrow REH(p,e)$ at $e{k}_{Hon}$, for $i\geq q$ |
| H dissociation | $REH(p,e)\rightarrow R(p,e)$ at ${k}_{Hoff}$ |
| Engaged E detachment | $REH(p,e)\rightarrow R(p,e-1)+{E}_{f}$ at ${k}_{Eoff,engaged}$ |
| Commitment (then cleavage and release) | $REH(p,e)\rightarrow {P}_{f}+e{E}_{f}$ at ${k}_{c}$ |

## **Initiation, elongation, and removal**

New polymerases entered node 1 in state $R(0,0)$ at flux ${k}_{in}{P}_{f}$, where ${k}_{in}$ is the initiation rate. Elongation advanced a polymerase by one node without changing $p$ or $e$. R states moved at the elongation rate ${k}_{e}=(65 {bp s}^{-1})/{L}_{a}$, corresponding to 65 bp/s, and REH states at the slower rate ${k}_{e2}=(30 {bp s}^{-1})/{L}_{a}$, corresponding to 30 bp/s. Reduced velocity was associated with recognition, and polymerases in state R retained ${k}_{e}$ downstream of the PAS.

Commitment removed REH polymerases at ${k}_{c}$. At the final node, R and REH polymerases could also leave by terminal run-off at their respective elongation rates. Both commitment (after cleavage) and run-off returned one Pol II and all $e$ attached E factors to the free pools. Run-off represented escape from the finite simulated window and was not counted as commitment.

## **Coupled kinetic equations and resource conservation**

Let $X$ collect all R and REH microstate populations, ${E}_{f}$ denote free E, and ${P}_{f}$ denote free Pol II. The bound-state dynamics were written as

$$\frac{dX}{dt}=({A}_{0}+{E}_{f}{A}_{E})X+{P}_{f}b.$$

${A}_{0}$ contains phosphorylation, dephosphorylation, dissociation, recognition, commitment, and transport. ${A}_{E}$ contains E-association transitions per unit free E. The source $b={k}_{in}u$ inflows polymerases into ${R}_{1}(0,0)$, where $u$ is a unit vector.

At steady state, the total E pool equals free E plus E bound across all R and REH states, and the total Pol II pool equals free Pol II plus all bound Pol II. We solved the linear steady-state equations for $X$ at each trial ${E}_{f}$, scaled the occupancy to conserve Pol II, and updated ${E}_{f}$ until E conservation also held. No rapid-equilibrium approximation was used for any binding step. If $w$ contains the bound-E count $e$ of each microstate and $1$ denotes a vector of ones, conservation requires

$${E}_{total}={E}_{f}+{w}^{T}X,  {P}_{total}={P}_{f}+{1}^{T}X.$$

Free-pool derivatives were derived from the same reaction operator: $d{E}_{f}/dt=-{w}^{T}A({E}_{f})X$ and $d{P}_{f}/dt=-{1}^{T}A({E}_{f})X-{k}_{in}{P}_{f}$. Consequently, binding consumes free E, E detachment restores it, initiation consumes free Pol II, and commitment (after cleavage) and run-off recycle both resources.

## **CPA Assembly Distance**

We quantified CPA assembly efficiency using the CPA assembly distance (CAD). We report ${CAD}_{50}$, the downstream distance from the PAS at which the cumulative commitment probability reaches 50%. Conceptually, ${CAD}_{50}$ is the downstream distance required for half of polymerases that encounter the PAS to commit to CPA. We chose this definition because it maps naturally onto cis-antisense rescue assays, in which protection increases when CPA commitment occurs before downstream interference.

![][image1]

*Figure M1. CAD calculation from the commitment flux. The commitment flux at downstream node l is k\_c · REH(l). CAD₅₀ is obtained by normalizing the cumulative exit-flux distribution and interpolating the distance at which the cumulative probability reaches 50%.*

For each steady-state solution, the local commitment flux was computed as k\_c · REH(l). Run-off from the terminal node was treated as a terminal flux to ensure that polymerases reaching the end of the simulated template without commitment were included in the output distribution. The normalized cumulative sum of downstream commitment and terminal fluxes defined the commitment CDF F(l), and ${CAD}_{50}$ was obtained by linear interpolation at F(l) \= 0.5.

## **Gene-length and alternative polyadenylation analyses**

To test whether gene architecture affects CPA efficiency, we varied the TSS-to-PAS distance while holding the downstream model structure fixed and computed ${CAD}_{50}$ for each simulated gene. To evaluate genome-wide resource competition, we precomputed a three-dimensional lookup grid for R\_occupied and E\_occupied as functions of free Pol II, free E, and gene length. These functions were interpolated and combined with a log-normal gene-length distribution to solve conservation equations for total Pol II and total E. This procedure links the single-gene assembly model to genome-wide competition for limiting transcription and CPA factors.

In the gene-length analysis, 10,000 genes shared the free pools, and each gene's initiation rate was ${k}_{in}$/10,000. TSS-to-PAS distances spanned 2.5–200 kb with a fixed 5 kb downstream extension. Gene lengths followed a log-normal distribution with a median of 22 kb and a standard deviation of 0.68 in log10 space, normalized over the simulated range. The shared free E and free Pol II were solved once for $M=5$ and then held fixed for every gene length and both capacities, so each comparison had the same factor supply. This analysis used the strong-PAS value ${k}_{Hoff}$ = 0.01 s⁻¹, which represents the canonical AAUAAA hexamer, the most frequent signal among human poly(A) sites (Beaudoing et al., 2000; Tian et al., 2005).

To model alternative polyadenylation (APA), we treated the commitment CDF downstream of a proximal PAS as the probability that a polymerase commits to that PAS before reaching a distal PAS. Proximal PAS usage was estimated by evaluating the CDF at the inter-PAS distance. PAS strength was varied by sweeping the PAS dissociation parameter, allowing us to ask how intrinsic PAS strength and spacing jointly control proximal-site usage. ${k}_{Hoff}$ was swept from 0.001 to 10 s⁻¹ with ${k}_{Hon}$ fixed, using $M=5$ and a 50 kb simulated gene. This broad exploratory range was motivated by reported nanomolar-to-micromolar affinities of PAS-containing RNAs for CPSF/mPSF, which depend on both hexamer sequence and flanking RNA (Hamilton et al., 2019; Huang et al., 2025). These affinities motivate the sweep but do not directly determine its individual dissociation rates.

The histogram in Figure 6B (lower) shows GENCODE inter-PAS distances (14,212 consecutive intervals between 30 bp and 1 Mb; median 1,521 nt), digitized from a supplied plot.

<!-- TODO: the extraction procedure for this histogram is undocumented (feature type: polyA_site, polyA_signal or transcript ends; last-exon restriction; deduplication or clustering; interval filters). Analysis_writeup.md covers only the MANE spans (Figure 5C) and has no inter-PAS data; see "PAS-spacing interpretation" in Docs/Results_framework.md. -->

## **Parameter maps, benchmark, and pool-size sweeps**

We scored commitment against a benchmark of ${CAD}_{50}$ ≤ 500 bp, with 200–500 bp used as a reference band for commitment within a few hundred base pairs, motivated by the cis-antisense rescue experiments of Chao et al. (1999). A condition in which the 50% commitment threshold was not reached counted as not meeting the benchmark. The paired parameter maps (Figure 2C) used 10 × 10 log-spaced grids with the third parameter at its reference value and the association rates ${k}_{Eon}$ and ${k}_{Hon}$ held fixed while the dissociation ratios varied. The ranges were kHd from 8.8 × 10⁻⁵ to 0.25, kEd from 2.4 × 10⁴ to 3.0 × 10⁶ molecules, and ${k}_{c}$ from 0.015 to 0.6 s⁻¹. For the joint analysis, we drew 10,000 independent combinations of kHd, kEd, and ${k}_{c}$, sampling each uniformly in log space over the same ranges and evaluating identical combinations for $M=1$ and $M=5$, with all other parameters at their reference values. The pool-size sweeps (Figure 4B) varied ${E}_{total}$ from 50,000 to 200,000 with ${P}_{total}$ at its reference value, and ${P}_{total}$ from 40,000 to 150,000 with ${E}_{total}$ at its reference value, both at $M=5$.

## **Human gene-length data**

To place the model's length dependence in a genomic context, we used MANE Select v1.5 (GRCh38), which reconciles RefSeq and Ensembl/GENCODE transcript annotations and uses CAGE and polyadenylation evidence to help define transcript ends (Morales et al., 2022). We kept one MANE Select NM_ transcript per protein-coding gene on chromosomes 1–22, X, and Y (19,229 genes, with no length filtering). We computed the distance from the annotated transcript start to the transcript 3′ end from the summary-file coordinates and verified it against the MANE RefSeq GFF annotation. These spans include introns and serve as a proxy for TSS-to-PAS distance. Sources: NCBI MANE documentation (https://www.ncbi.nlm.nih.gov/refseq/MANE/) and MANE v1.5 files (https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/).

## **Implementation and parameters**

Simulations were implemented in MATLAB. The pipeline builds one sparse operator containing all finite-rate reactions and transport, solves the resulting linear microstate equations at each trial free E, scales occupancy to conserve Pol II, and uses fzero to conserve total E. Unless otherwise stated, figures used the reference parameter set in Table M2. The parameter-estimation appendix (Appendix) derives each value and the ranges explored.

**Table M2.** Reference parameter values. Rates are per second, and elongation rates are per 100 bp node.

| Parameter | Symbol | Reference value | Basis |
| :---- | :---- | :---- | :---- |
| Node spacing | ${L}_{a}$ | 100 bp | lattice discretization |
| Template length; PAS position | $G$; ${x}_{PAS}$ | 25 kb; 20 kb | single-gene default; the gene-length analysis varies ${x}_{PAS}$ with a 5 kb downstream extension |
| CTD capacity | $M$ | 5 | varied from 1 to 10 |
| Pol II initiation rate | ${k}_{in}$ | 2 s⁻¹ (divided by 10,000 per gene in the 10,000-gene population) | working value |
| Pre-PAS elongation (R) | ${k}_{e}$ | 0.65 node s⁻¹ (65 bp/s) | live-cell and nascent-RNA elongation rates, Appendix E |
| Post-PAS elongation (REH) | ${k}_{e2}$ | 0.30 node s⁻¹ (30 bp/s) | modeling assumption, Appendix E |
| Site formation | ${k}_{Pon,min}$; ${k}_{Pon,slope}$ | 0; 5 × 10⁻⁴ s⁻¹ per node | assumed linear ramp, Appendix F |
| Site loss | ${k}_{Poff}$ | 0.01 s⁻¹ | working assumption, Appendix F |
| E association | ${k}_{Eon}$ | 2.5 × 10⁻⁶ molecule⁻¹ s⁻¹ | basal association rate scaled by nuclear volume, Appendix C |
| E dissociation | ${k}_{Eoff}$ | 0.5 s⁻¹ | CTD-binding affinities of about 1 μM, Appendix C |
| Engaged-E detachment | ${k}_{Eoff,engaged}$ | 0.1 s⁻¹ | between the limiting cases, Appendix C |
| PAS recognition | ${k}_{Hon}$ | 7 s⁻¹ per bound E | tethered-encounter (worm-like chain) estimate, Appendix B |
| H dissociation | ${k}_{Hoff}$ | 1 s⁻¹ (weak PAS; reference); 0.01 s⁻¹ (strong PAS; gene-length analysis) | measured mPSF–PAS affinities, Appendix B |
| Commitment | ${k}_{c}$ | 0.14 s⁻¹ | three-step late-factor scheme, Appendix D |
| Total E | ${E}_{total}$ | 100,000 molecules | limiting CPSF subunit, Appendix G |
| Total Pol II | ${P}_{total}$ | 100,000 molecules | see Appendix G |
| Derived ratios | kEd; kHd | 2.0 × 10⁵ molecules; 0.14 | ${k}_{Eoff}/{k}_{Eon}$; ${k}_{Hoff}/{k}_{Hon}$ |

<!-- TODO: Appendix G states P_total = 75,000 molecules; the reference runs and default_parameters.m use 100,000. Reconcile. -->


# **References**

<!-- TODO: "Jin et al., 2025" (Introduction, paragraph on CPA-factor perturbation) could not be matched to a reference; add it here once confirmed. -->

1. Beaudoing, E., Freier, S., Wyatt, J. R., Claverie, J. M., & Gautheret, D. (2000). Patterns of variant polyadenylation signal usage in human genes. *Genome Research*, 10(7), 1001–1010. https\://doi.org/10.1101/gr.10.7.1001

2. Bekker-Jensen, D. B., Kelstrup, C. D., Batth, T. S., Larsen, S. C., Haldrup, C., Bramsen, J. B., Sørensen, K. D., Høyer, S., Ørntoft, T. F., Andersen, C. L., Nielsen, M. L., & Olsen, J. V. (2017). An Optimized Shotgun Strategy for the Rapid Generation of Comprehensive Human Proteomes. *Cell Systems*, 4(6), 587–599.e4. https\://doi.org/10.1016/j.cels.2017.05.009

3. Casañal, A., Kumar, A., Hill, C. H., Easter, A. D., Emsley, P., Degliesposti, G., Gordiyenko, Y., Santhanam, B., Wolf, J., Wiederhold, K., Dornan, G. L., Skehel, M., Robinson, C. V., & Passmore, L. A. (2017). Architecture of eukaryotic mRNA 3'-end processing machinery. *Science*, 358(6366), 1056–1059. https\://doi.org/10.1126/science.aao6535

4. Chao, L. C., Jamil, A., Kim, S. J., Huang, L., & Martinson, H. G. (1999). Assembly of the cleavage and polyadenylation apparatus requires about 10 seconds in vivo and is faster for strong than for weak poly(A) sites. *Molecular and Cellular Biology*, 19(8), 5588–5600. https\://doi.org/10.1128/MCB.19.8.5588

5. Darzacq, X., Shav-Tal, Y., de Turris, V., Brody, Y., Shenoy, S. M., Phair, R. D., & Singer, R. H. (2007). In vivo dynamics of RNA polymerase II transcription. *Nature Structural & Molecular Biology*, 14(9), 796–806. https\://doi.org/10.1038/nsmb1280

6. Davidson, L., Muniz, L., & West, S. (2014). 3' end formation of pre-mRNA and phosphorylation of Ser2 on the RNA polymerase II CTD are reciprocally coupled in human cells. *Genes & Development*, 28(4), 342–356. https\://doi.org/10.1101/gad.231274.113

7. Fuchs, G., Voichek, Y., Benjamin, S., Gilad, S., Amit, I., & Oren, M. (2014). 4sUDRB-seq: measuring genomewide transcriptional elongation rates and initiation frequencies within cells. *Genome Biology*, 15(5), R69. https\://doi.org/10.1186/gb-2014-15-5-r69

8. Fusby, B., Kim, S., Erickson, B., Kim, H., Peterson, M. L., & Bentley, D. L. (2016). Coordination of RNA Polymerase II Pausing and 3' End Processing Factor Recruitment with Alternative Polyadenylation. *Molecular and Cellular Biology*, 36(2), 295–303. https\://doi.org/10.1128/mcb.00898-15

9. Grzechnik, P. & Mischo, H. E. (2025). Fateful Decisions of Where to Cut the Line: Pathology Associated with Aberrant 3' End Processing and Transcription Termination. *Journal of Molecular Biology*, 437(1), 168802. https\://doi.org/10.1016/j.jmb.2024.168802

10. Hamilton, K., Sun, Y., & Tong, L. (2019). Biophysical characterizations of the recognition of the AAUAAA polyadenylation signal. *RNA*, 25(12), 1673–1680. https\://doi.org/10.1261/rna.070870.119

11. Huang, L., Chu, H. F., & Tong, L. (2025). Molecular basis for the recognition of low-frequency polyadenylation signals by mPSF. *Nucleic Acids Research*, 53(17), gkaf890. https\://doi.org/10.1093/nar/gkaf890

12. Jackson, D. A., Iborra, F. J., Manders, E. M., & Cook, P. R. (1998). Numbers and organization of RNA polymerases, nascent transcripts, and transcription units in HeLa nuclei. *Molecular Biology of the Cell*, 9(6), 1523–1536. https\://doi.org/10.1091/mbc.9.6.1523

13. Kimura, H., Tao, Y., Roeder, R. G., & Cook, P. R. (1999). Quantitation of RNA polymerase II and its transcription factors in an HeLa cell: little soluble holoenzyme but significant amounts of polymerases attached to the nuclear substructure. *Molecular and Cellular Biology*, 19(8), 5383–5392. https\://doi.org/10.1128/mcb.19.8.5383

14. Kimura, H., Sugaya, K., & Cook, P. R. (2002). The transcription cycle of RNA polymerase II in living cells. *Journal of Cell Biology*, 159(5), 777–782. https\://doi.org/10.1083/jcb.200206019

15. Kowalski, M. H., Wessels, H. H., Linder, J., Dalgarno, C., Mascio, I., Choudhary, S., Hartman, A., Hao, Y., Kundaje, A., & Satija, R. (2024). Multiplexed single-cell characterization of alternative polyadenylation regulators. *Cell*, 187(16), 4408–4425.e23. https\://doi.org/10.1016/j.cell.2024.06.005

16. Martin, S. & Pombo, A. (2003). Transcription factories: quantitative studies of nanostructures in the mammalian nucleus. *Chromosome Research*, 11(5), 461–470. https\://doi.org/10.1023/a:1024926710797

17. Mitschka, S. & Mayr, C. (2022). Context-specific regulation and function of mRNA alternative polyadenylation. *Nature Reviews Molecular Cell Biology*, 23(12), 779–796. https\://doi.org/10.1038/s41580-022-00507-5

18. Morales, J., Pujar, S., Loveland, J. E., Astashyn, A., Bennett, R., Berry, A., Cox, E., Davidson, C., Ermolaeva, O., Farrell, C. M., Fatima, R., Gil, L., Goldfarb, T., Gonzalez, J. M., Haddad, D., Hardy, M., Hunt, T., Jackson, J., Joardar, V. S., … Murphy, T. D. (2022). A joint NCBI and EMBL-EBI transcript set for clinical genomics and research. *Nature*, 604(7905), 310–315. https\://doi.org/10.1038/s41586-022-04558-8

19. Neve, J., Patel, R., Wang, Z., Louey, A., & Furger, A. M. (2017). Cleavage and polyadenylation: Ending the message expands gene regulation. *RNA Biology*, 14(7), 865–890. https\://doi.org/10.1080/15476286.2017.1306171

20. Rieder, D., Trajanoski, Z., & McNally, J. G. (2012). Transcription factories. *Frontiers in Genetics*, 3, 221. https\://doi.org/10.3389/fgene.2012.00221

21. Tang, P. & Zhou, Y. (2022). Alternative polyadenylation regulation: insights from sequential polyadenylation. *Transcription*, 13(4-5), 89–95. https\://doi.org/10.1080/21541264.2022.2114776

22. Tian, B., Hu, J., Zhang, H., & Lutz, C. S. (2005). A large-scale analysis of mRNA polyadenylation of human and mouse genes. *Nucleic Acids Research*, 33(1), 201–212. https\://doi.org/10.1093/nar/gki158
