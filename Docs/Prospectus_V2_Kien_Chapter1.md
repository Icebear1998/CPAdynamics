# Chapter 1: CPA Dynamics Regulation {#chapter-1:-cpa-dynamics-regulation}

## Section 1.1: Background {#section-1.1:-background}

### **1.1.1 Cleavage and Polyadenylation Process** {#1.1.1-cleavage-and-polyadenylation-process}

Cleavage and polyadenylation (CPA) is a fundamental and indispensable process in eukaryotic gene expression that facilitates transcription termination and adds a poly(A) tail to newly transcribed RNAs. This intricate mechanism is not merely an endpoint but a critical hub for post-transcriptional control, determining the fate of the messenger RNA (mRNA). The proper formation of the poly(A) tail is essential for several vital cellular functions, including mRNA stability, nuclear export, translation into protein, and even the molecule's specific localization within the cell. Without proper CPA, the entire cascade of gene expression would be severely compromised.

The CPA process is initiated when the elongating RNA Polymerase II (Pol II) transcribes a specific genetic sequence element known as the polyadenylation signal (PAS) on the nascent mRNA. The PAS signal, most commonly the hexamer AAUAAA, acts as a recognition site for assembling the multi-protein CPA machinery. Once assembled at the PAS, the machinery performs a precise cleavage of the pre-mRNA, typically 10-30 nucleotides downstream of the recognition site, and adds a poly(A)-tail to the 3’-end of the cleaved RNA. This CPA machinery formation event is the crucial first step that both defines the 3' end of the new mRNA and serves as a signal to trigger the termination of transcription by Pol II.

Furthermore, the CPA process provides a powerful layer of regulatory control. Most eukaryotic genes contain more than one PAS sequence, which allows the CPA machinery to choose between different sites. This process, known as Alternative Polyadenylation (APA), generates distinct mRNA isoforms from a single gene. These isoforms are identical in their protein-coding sequence but differ in the length and content of their 3' untranslated regions (3' UTRs). Since longer 3' UTRs generally contain more binding sites for regulatory molecules (such as microRNAs and RNA-binding proteins) than shorter ones, APA allows the cell to switch between stable, highly expressed isoforms and strictly regulated ones, thereby fine-tuning protein output.

### **1.1.2 The CPA Machinery and Physical Constraints** {#1.1.2-the-cpa-machinery-and-physical-constraints}

The CPA machinery is not a single protein but a large, dynamic conglomerate of multiple protein complexes that must assemble in a coordinated manner on the elongating Pol II. While many factors have been identified, the core machinery includes essential components like the Cleavage and Polyadenylation Specificity Factor (CPSF), which recognizes the PAS, and the Cleavage Stimulation Factor (CstF), recognizes the downstream signal in the PAS. The successful and timely recruitment of these and other factors is a prerequisite for the cleavage reaction and, ultimately, for proper transcription termination.

**Figure 1.1. Schematic of the Dynamic Assembly of CPA Machinery for the CPA and APA mechanisms.** **(A \- Upper) Assembly of CPA machinery:** CPA factors are recruited and pre-loaded onto the Pol II CTD scaffold while the polymerase is elongating. **(A \- Lower)** Upon transcription of a proximal polyadenylation signal (PAS), CPA factors continue to be recruited and docked CPA factors will interact with the PAS sequence, leading to the assembly of the full CPA machinery and cleavage of the nascent RNA. **(B)** **APA mechanism:** CPA can happen at any PAS site in the gene, resulting in mRNA isoforms with diferent 3' UTRs.

Experimental studies provide evidence that some CPA factors, like CPSF30, can bind to Pol II’s CTD during transcription initiation \[3\]. Many other CPA factors, like CPSF160, WDR33, and Fip1, will bind to Pol II’s CTD during transcription elongation before the PAS sequence is expressed. Yet another class of CPA factors can bind to Pol II’s CTD after the PAS sequence is expressed \[3\]. Based on this, we categorize CPA factors into two functional groups: '**Early' factors**, which can bind to the CTD prior to PAS transcription, and '**Late' factors**, which are recruited only after the PAS sequence is expressed. This categorization can be seen within Fig. 1.1 (B).

A key challenge to understanding this process comes from the physical constraints of the cellular environment. A central constraint is competition. Quantitative proteomic data reveals that core CPA factors are present at near-stoichiometric levels with actively transcribing Pol II, rather than in vast excess \[5\]. This means that CPA factors are a limited resource within the nucleoplasm. This creates a scenario of intense genome-wide competition, where thousands of Pol II complexes simultaneously vie for the same small pool of available CPA factors. This raises a critical question: how does any single polymerase manage to gather all the necessary components when so many others are competing for them at the same time?

Seminal experimental work has shown that the CPA machinery must fully assemble within a short "runway" of only **400 to 800 base pairs** past the PAS \[1\]. Given typical Pol II elongation rates, this timeframe translates to a window time of **10 to 20 seconds** after the PAS sequence emerges. If the machinery fails to assemble within this distance, the polymerase will continue transcribing, leading to transcriptional readthrough into downstream genes, which can have disruptive and often deleterious consequences. Therefore, the central problem that motivates this research is to understand how the cell resolves these competing pressures: how can the CPA machinery reliably and rapidly assemble within this tight spatiotemporal window when its components are limited and subject to intense genome-wide competition? 

### **1.1.3 The RNA Polymerase II CTD: A Dynamic Assembly Platform** {#1.1.3-the-rna-polymerase-ii-ctd:-a-dynamic-assembly-platform}

The key architectural feature that allows the cell to overcome the challenges of competition and time is the C-terminal domain (CTD) of RNA Polymerase II. The CTD is a long, intrinsically disordered domain extending from the core of the polymerase, composed of up to 52 tandem repeats of a heptad consensus sequence (YSPTSPS). This repetitive structure is not redundant; rather, it functions as a highly dynamic and programmable scaffold, or "landing pad," for a vast array of protein machineries involved in transcription and co-transcriptional RNA processing, including the CPA machinery.

Crucially, the Pol II CTD scaffold is dynamically regulated throughout transcription. The CTD undergoes extensive post-translational modifications, most notably phosphorylation, which creates a complex "CTD code" that dictates which factors can bind at which stage of the transcription cycle. For 3'-end processing, phosphorylation at the Serine-2 position of the heptad repeat is paramount. As shown by experimental ChIP-chip data, the level of Serine-2 phosphorylation ramps up progressively as Pol II travels from the 5' promoter region toward the 3' end of a gene.

### **1.1.4 Key Research Questions** {#1.1.4-key-research-questions}

Based on the existing knowledge about the CPA assembly process, we ask the following questions:

A. How does the CPA machinery efficiently assemble despite the competition for CPA factors among Pol-IIs? This is the core question addressing the kinetic challenge. We aim to determine if the proposed mechanism—the pre-loading of factors onto a progressively phosphorylated CTD scaffold—is quantitatively sufficient to explain the experimentally observed assembly time of 10-20 seconds.

B. What is the functional significance of the multiple heptameric repeats on the Pol II CTD? This question probes the architectural logic of the system. Why did a CTD with numerous potential binding sites evolve? We will investigate the quantitative relationship between the number of available binding sites, the efficiency of assembly, and the ability to overcome competition.

C. How does the physical assembly mechanism influence the selection between different PAS sites? This question explores the regulatory dimension of the process. If the assembly is a dynamic, distance-dependent process, how does this affect the cell's choice when faced with multiple, closely spaced PAS signals? We aim to determine if the physical constraints of assembly provide a mechanistic basis for regulating Alternative Polyadenylation (APA).

## Section 1.2: Methods {#section-1.2:-methods}

### **1.2.1 A Multi-Scale Model of CPA Machinery Assembly** {#1.2.1-a-multi-scale-model-of-cpa-machinery-assembly}

To capture the assembly and cleavage process of the CPA machinery, we developed a mechanistic, multi-scale model that maps a sequence of biochemical events into …. The model includes four key processes: *(i)* *Pol II Elongation*, RNA polymerase II (R) transcribes the gene and advances along a one-dimensional coordinate. *(ii) Site Activation,* progressive phosphorylation of Ser2 residues on *the C-terminal domain (CTD) of* the Pol II CTD turns into binding sites for CPA factors. *(iii) Early CPA factor assembly*, early CPA factors (E) bind reversibly to phosphorylated CTD sites during elongation. *(iv) PAS binding and commitment*, after transcription of the polyadenylation signal (PAS), the Pol II–CPA complex undergoes a commitment step leading to cleavage and termination.

**Model assumptions and formulation**

**Figure 1.2. The Pol II CTD as a dynamic, phosphorylation-dependent scaffold.** **(Top)** In the early stages of elongation, the CTD repeats are largely unphosphorylated at the Serine-2 position (light purple), providing very few high-affinity binding sites for early CPA factors (E), which remain free in the nucleoplasm. **(Bottom)** As Pol II transcribes the gene, progressive phosphorylation of the Serine-2 positions (dark purple) creates numerous high-affinity binding sites. This allows multiple early CPA factors (E) to be co-transcriptionally recruited and "pre-load" onto the CTD. This prepares the complex for rapid assembly with both late CPA factors (L) and the PAS once it is manufactured from the DNA template.

For each process, the key assumptions represented in the model are as follows. 

**(i) Transcriptional elongation**

Transcription is modeled as a one-dimensional reaction–advection system. The gene is discretized into nodes of length L\_a, and a coordinate l tell us where we are on the gene (position from transcription start site \= l\*L\_a). The number of polymerases R(l,t) is tracked at each node. Polymerases advance at an elongation rate k\_e​.

**(ii) Progressive CTD activation (Ser2 phosphorylation)**

Experimental data show that Ser2 phosphorylation gradually increases toward the 3′ end of genes. We represent this using a position-dependent activation rate k\_p(l), which increases as Pol II moves from the transcription start site toward the PAS.

As a result, the number of phosphorylated CTD repeats — and thus available binding sites for early CPA factors — increases with position. This establishes a spatial gradient of assembly competence along the gene (Figure 1.2).

**(iii) Early CPA factor recruitment**

Early-factor recruitment is modeled as reversible binding to phosphorylated Ser2 sites:

\\text{Ser2P} \+ E \\xrightleftharpoons\[k\_{E\\\_off}\]{k\_{E\\\_on}} RE

Each phosphorylated CTD repeat is treated as a potential independent binding site, allowing multiple E factors to bind simultaneously to a single polymerase.

Because phosphorylation and binding/unbinding events occur on a timescale (seconds) much faster than elongation (minutes), we assume these reactions reach quasi-steady-state equilibrium at each gene position. For each node, we construct a symbolic transition matrix describing all binding/unbinding states and compute its null space to obtain the equilibrium distribution of CTD occupancy.

From this distribution, we calculate the average number of bound early factors at position lll, denoted:

⟨nE(l)⟩\\langle n\_E(l) \\rangle⟨nE​(l)⟩

This quantity depends on:

* the local phosphorylation level (via kp(l)k\_p(l)kp​(l)), and

* the globally available pool of free factors EfreeE\_{free}Efree​.

**(iv) PAS recognition and commitment**

When a polymerase reaches the PAS node, it can transition to a PAS-bound committed state (denoted REH), initiates a commitment step toward cleavage and termination. This transition occurs at a rate proportional to a commitment constant k\_{Hon}​, which depends on the number of early CPA factors loaded onto the CTD:

R→REH  
with rate  
k\_{Hon}(l) \= f\\big(\\langle n\_E(l) \\rangle\\big)

Importantly, k\_{Hon}​ is not constant. It is parameterized by the average early-factor occupancy computed from the fast equilibrium module. Higher CTD occupancy increases the probability of commitment, creating a mechanistic coupling between co-transcriptional assembly and termination efficiency.

This model enables us to quantitatively measure assembly efficiency using experimentally observed data on CTD activation, factor recruitment, and PAS-dependent commitment (Figure 1.2).

#### **Time-scale separation and multi-scale coupling** {#time-scale-separation-and-multi-scale-coupling}

Directly simulating every binding microstate across the entire gene and genome is computationally prohibitive. We therefore exploit the separation of timescales in the model. More specifically, phosphorylation of Ser2P and the binding and unbinding of early CPA factors (“E factors”) to phosphorylated Ser2 residues of the Pol II CTD occur on a timescale much faster than transcriptional elongation, enabling the construction of a multi-scale framework that explicitly separates fast and slow reactions. Rapid site activation and factor binding (processes ii and iii) are treated as quasi-steady-state fast processes. Polymerase elongation and cleavage commitment (process iv) are modeled as slower, time-resolved processes. This separation keeps the model computationally tractable while preserving the mechanistic couplings needed. 

Building on this scheme, the multi-scale scheme is organized around a triangular relationship (Figure 1.3) that links three central components: (1) *Fast Processes—r*eversible binding of CPA factors to phosphorylated Ser2 sites on the CTD, reaching quasi-equilibrium. At each node, compute the equilibrium distribution of CTD binding states and extract ⟨n E (l)⟩ as a function of local Ser2P and the global free pool E\_free. (2) *Slow Processes*—progression of Pol II along the gene and commitment at a rate (kHon) that depends on the number of bound factors calculated in the Fast Process. Propagate polymerase densities along the gene using advection; when polymerases reach the PAS node, apply a commitment propensity k\_Hon(⟨nE⟩) to convert a fraction into the committed REH state and trigger cleavage/termination; and (3) Global Constraints—the limited pool of free CPA factors (Efree), which connects local molecular events to genome-wide competitive dynamics.

**Figure 1.3:** **Multi-scale coupling of fast molecular binding and slow transcriptional elongation.** The model is organized as a triangular scheme linking three key components: (1) *Fast Process*: Reversible binding of CPA factors to phosphorylated Ser2 sites on the CTD reaches a quasi-steady-state equilibrium. (2) *Slow Process*: Pol II elongates along the gene and commits to termination at a rate (kHon) dependent on the number of bound factors. (3) *Global Constraints*: The pool of free CPA factors (Efree) limits binding availability, coupling local assembly to genome-wide competition.

By coupling timescales, the framework provides a fully explicit description of all binding states across the gene, these states are crucial for determining how quickly the CPA machinery forms. The multi-scale strategy, therefore, provides a tractable means to probe termination efficiency, resource competition, and alternative polyadenylation outcomes.

## 1.2.2 Model Parameterization

Model parameters were derived using a combination of experimentally measured kinetic constants, quantitative proteomics databases, and first-principles biophysical calculations. Detailed derivations and sensitivity analyses for all physically estimated parameters are provided in the Appendix.

*Kinetic Rates from Literature and Biophysical Derivations*

The RNA polymerase II elongation rate (k\_e) and cleavage rate (k\_c) were selected within the experimentally observed in vivo range for mammalian transcription (BioNumbers database).

The PAS recognition and dissociation rates (k\_Hon and k\_Hoff) were estimated using a loop-mediated encounter framework based on polymer-physics–derived J-factor calculations and experimentally measured binding affinities. Full derivations are provided in the Appendix.

Association and dissociation rates for Early CPA factors (k\_Eon and k\_Eoff) were estimated from diffusion-limited binding principles under nuclear conditions and experimentally reported weak CTD-binding affinities. Detailed derivations and sensitivity bounds are provided in the Appendix.

CTD phosphorylation dynamics were parameterized to reproduce the experimentally observed 5′–3′ Ser2 phosphorylation gradient. The phosphorylation on-rate was modeled as a position-dependent linear ramp, k\_Pon(l) \= k\_Pon\_min \+ k\_Pon\_slope · l, while the dephosphorylation rate (k\_Poff) was chosen as a constant to maintain a dynamic steady-state gradient.

*Total Protein Abundances from Proteomic Databases*

Total Pol II abundance (Pol\_total) and total Early CPA factor abundance (E\_total) were estimated from consensus measurements in PaxDB for human cell lines. Conversion to nuclear copy numbers and parameter justification are described in the Appendix.

*Estimation of Maximum CPA Factor Binding Capacity*

A structural constraint in the model is the maximum number of CPA factors that can simultaneously occupy the CTD scaffold (N\_binding). Although the human CTD contains 52 heptad repeats, steric constraints and functional competition with other co-transcriptional machineries limit the number of repeats available for CPA assembly at any given time. The justification and supporting simulations for the chosen value are provided in the Appendix.

**Table 1.1. Key Model Parameters**			

| Parameter | Description | Value | Unit | Source |
| ----- | ----- | ----- | ----- | ----- |
| `Pol_total` | Total Pol II abundance | 70,000 | molecules | PaxDB / Model |
| `E_total` | Total CPA factor abundance | 70,000 | molecules | PaxDB / Model |
| `k_e` | Elongation rate | 2 | kb/min | Literature |
| `k_c` | Cleavage rate | 0.2 | s⁻¹ | Model fit |
| `N_binding` | Max CPA factor binding sites | 5 | sites | Structural constraint |
| `k_Hon` | PAS recognition rate | 0.05 | s⁻¹ | Appendix derivation |
| `k_Hoff` | PAS dissociation rate | 3 × 10⁻³ | s⁻¹ | Appendix derivation |
| `k_Eon` | Early factor association rate | 1.7 × 10⁻⁵ | molecule⁻¹ s⁻¹ | Appendix derivation |
| `k_Eoff` | Early factor dissociation rate | 50–500 | s⁻¹ | Appendix derivation |
| `k_Pon_slope` | Phosphorylation rate slope | 0.02 | s⁻¹ | Model fit |
| `k_Poff` | Dephosphorylation rate | 1 | s⁻¹ | Model assumption |

### **1.2.3 Quantification of CPA Assembly Efficiency**

### *CPA Assembly Distance (CAD) and the cis-antisense rescue assay.*

To quantify CPA assembly from model outputs we define the **CPA Assembly Distance (CAD)** as the downstream distance (bp) from the PAS at which 50% of polymerases have become irreversibly committed to cleavage/termination. The CAD is designed to directly comparable to the published cis-antisense rescue assay. In that assay, the authors systhesysed gene with antisense of the hexammer at different distance downstream of PAS. Transcription of the antisense element generates a complementary RNA that can base-pair with the newly transcribed PAS region and prevent CPA factors binding. When the antisense sequence is placed progressively farther downstream, inhibition is lost once the authentic PAS has already advanced through the necessary assembly steps and become committed before the antisense transcript can interfere. The reported rescue distance therefore reflects the transcriptional distance (and corresponding time) required for irreversible CPA commitment, which measures the downstream distance at which an inserted antisense copy of the PAS no longer blocks processing of the authentic site. For the canonical SV40 early site this 50% rescue distance is ≈200 bp (≈10–20 s when converted using typical Pol II speeds), providing an experimental benchmark for model calibration (/cite).

### *Mapping the model to the experimental CAD metric*

We compute a model-derived CAD that is directly comparable to the experimental distance above using the following pipeline.

1. **Steady-state solution (homeostatic assumption).** We solve the ODE system to steady state and obtain spatial densities of elongating polymerases 𝑅\_𝑠𝑜𝑙(𝑙) and termination-committed complexes 𝑅𝐸𝐻\_𝑠𝑜𝑙(𝑙) at each discrete node downstream of the PAS.  
2. **Exit flux vector.** Using the steady-state committed density, the cleavage exit flux at each downstream node is  
   			fluxcleavage(𝑙)  =  𝑘𝑐⋅𝑅𝐸𝐻𝑠𝑜𝑙(𝑙)  
   (we also include run-off flux at the terminal node to capture polymerases that reach the end of the discretized template without being cleaved).  
3. **Termination CDF.** Form the cumulative sum of the exit flux starting at the PAS and normalize by the total outflux to obtain the termination cumulative distribution function  
   		𝐹(𝑙)  =  ∑𝑙PAS𝑙fluxcleavage(𝑙′)/∑𝑙PAS𝐿endfluxcleavage(𝑙′),​  
   which gives the probability that a polymerase passing the PAS is terminated by the time it reaches position 𝑙.  
4. **CAD definition (model → experiment).** The model CAD is the downstream distance (in bp) where the termination CDF reaches 0.5:  
   		CAD\_model  =  𝑙, such that 𝐹(𝑙)=0.5,  
   computed by linear interpolation between discretized nodes that bracket the 50% point.

The cis-antisense assay probes the time required for a poly(A) site to become functionally resistant to an antisense antagonist; our CAD measures the spatial location where 50% of polymerases have undergone cleavage driven by assembled CPA machinery — these two quantities target the same underlying assembly kinetics and are therefore directly comparable when distance is converted to time by the elongation rate.

## Section 3: Results {#section-3:-results}

### **3.1 The CTD Scaffold with Multiple Binding Sites Explains Efficient CPA Assembly** {#3.1-the-ctd-scaffold-with-multiple-binding-sites-explains-efficient-cpa-assembly}

A central question is how the CPA machinery and the CTD's architectural features can quantitatively account for the rapid assembly of the CPA machinery despite the intense competition for its components. To address this question, we developed a mathematical model to depict the experimentally observed subprocesses of CPA machinery assembly and test the impact of progressive activation of multiple CPA binding sites on the assembly efficiency . We chose parameters based on biochemical measurements and typical physiological ranges (Table 1.1) and characterize the assembly efficiency via the **CPA Assembly Distance (CAD)**.

**Figure 1.4**: **Multi-site CTD scaffolding significantly improves CPA assembly efficiency.** The plot shows the predicted CPA Assembly Distance (CAD) as a function of the number of available binding sites on the Pol II CTD. A single binding site results in a prolonged assembly distance (\~1329 bp), exceeding the biologically observed window. Increasing the number of binding sites to five reduces the CAD to \~657 bp, consistent with experimental estimates (400-800 bp), demonstrating that the multi-valent scaffold mechanism is sufficient to explain rapid assembly kinetics.

The result in Figure 1.4 shows a clear and strong relationship between the CTD's binding capacity and the efficiency of CPA assembly. When the model was restricted to a single binding site—representing a simple, non-scaffold scenario—the predicted CAD was **1329 bp**. This distance is significantly longer than the 400-800 bp range estimated from experimental data, suggesting a single-site binding model is insufficient to explain the observed kinetics.

However, as the number of available binding sites increased, the CAD systematically decreased. With five available binding sites, the CAD was reduced to **657 bp**. This value falls within the experimentally observed range. This demonstrates that a multi-site "pre-loading" or "scaffolding" mechanism is sufficient to explain the known efficiency of CPA assembly. And the gradual increase in phosphorylated binding sites along the gene addresses the competition for CPA components and resources, which are concentrated toward the PAS site, where it needs to complete the machinery.

In our model, we chose five binding sites as a biophysically realistic upper limit, considering that the large physical size of some early CPA factors (e.g., CPSF160, CPSF30, CstF77, CstF64) and recruitment of other essential machineries (e.g., for splicing and capping) would significantly reduce the possible occupancy number of CPA factors . This result provides strong quantitative support for the hypothesis that the CTD's function as a multi-valent dynamic scaffold is the key mechanism that allows the cell to overcome kinetic barriers and ensure efficient transcription termination.

### **1.3.2 Model Prediction: CPA Efficiency is Strongly Dependent on Gene Length** {#1.3.2-model-prediction:-cpa-efficiency-is-strongly-dependent-on-gene-length}

Having established that our model can reproduce the known kinetics of CPA assembly, we next used it to make novel, testable predictions. We investigated how a fundamental architectural feature of a gene—its length—influences the efficiency of the CPA process. To do this, we simulated the CPA assembly process for genes with varying TSS-to-PAS distances and calculated the resulting CPA Assembly Distance (CAD).

The Figure 1.5 reveals a strong, non-linear relationship between gene length and CPA efficiency. The curve illustrates two distinct regimes:

**Figure 1.5: CPA assembly efficiency depends non-linearly on gene length.** The graph plots the CPA Assembly Distance (CAD) against the distance from the Transcription Start Site (TSS) to the PAS. For short genes (\<20 kb), the CAD is significantly elevated, indicating inefficient assembly due to insufficient time for factor recruitment (the "runway" effect). As gene length increases beyond 40 kb, the CAD decreases and plateaus at \~1350 bp, representing the maximal efficiency achievable when the CTD scaffold is fully saturated with CPA factors.

For **short genes** (with a TSS-to-PAS distance of less than 20 kb), the CAD is exceptionally long, as seen by the steep initial descent of the curve. The model suggests a clear physical mechanism for this: these genes provide an insufficient "runway" for the Pol II CTD to travel along and accumulate the necessary number of CPA factors before reaching the PAS.

For **long genes** (with a TSS-to-PAS distance greater than 40 kb), CAD shortens dramatically, reaching a stable plateau around 650 bp. This is visible as the flat part of the curve extending to the right. In this regime, the gene provides ample time and space for the CTD scaffold to fully load with CPA factors, so adding further length does not confer any additional benefit to the assembly process.

This prediction provides a powerful insight: gene architecture itself can be a key determinant of RNA processing efficiency. Furthermore, this result aligns with biological observations that very short Pol II transcripts, such as those for small nuclear RNAs (snRNAs), often employ distinct, specialized termination mechanisms, possibly because the standard CPA-dependent pathway is inherently inefficient for them. 

### **1.3.3 Mechanistic Insights into the Regulation of Alternative Polyadenylation (APA)** {#1.3.3-mechanistic-insights-into-the-regulation-of-alternative-polyadenylation-(apa)}

Finally, we applied our validated model to explore the regulatory logic of Alternative Polyadenylation (APA), a widespread mechanism of gene regulation where the choice between multiple PAS sites on a single transcript generates distinct mRNA isoforms. Genomic data reveal that, for human genes containing multiple PAS sites within the 3'-most exon, the distance between adjacent sites is not random; it shows a lognormal distribution peaking around 300 nucleotides and spanning about two orders of magnitude (Figure 1.6).  Is the observed range of inter-PAS distances important for the functioning of APA?

**Figure 1.6: Distribution of distances between adjacent poly(A) sites in human 3'-most exons.** The histogram displays the frequency of distances between adjacent PAS sites. The distribution is approximately log-normal, with a peak around 300 nucleotides and spanning roughly two orders of magnitude. This observed spacing provides a dynamic range that may facilitate regulation. (Tian et al. The RNA Polymerase II Carboxy-Terminal Domain (CTD) Code. \*Oxford Academic. Nucleic Acids Research.\* (2013))

![][image1]![][image2]  
**Figure 1.7:** To answer the above question, we used the model to estimate the probability of assembling CPA machinery at a proximal (first) PAS before a second, distal PAS is expressed as a function of their distance. The figure above presents the results of this analysis. Let's first focus on the zoomed-in plot (Figure 1.7), which examines a small range of the most frequent distances. The dashed line marks the typical 300 bp inter-PAS distance. Our model predicts that at this distance, the proximal PAS is used only approximately 25% of the time for a standard PAS strength. This result provides a crucial mechanistic insight: the system has an intermediate likelihood to assemble CPA machinery at the first PAS site it encounters. Such an intermediate value provides a proper dynamic range for regulation. Because the usage is not near 0% or 100%, the cell can use other molecular inputs, such as changes in CPA factor concentrations, to dial the usage of the proximal site up or down, thereby exerting precise control over which mRNA isoform is produced.

Figure 1.7 shows the "bigger picture," illustrating this same principle across a much wider range of inter-PAS distances and for different intrinsic PAS strengths (spanning a 1000-fold difference in `kHoff`). This range was chosen based on experimental studies showing that the binding affinity (`K_d`) of different PAS sequences can vary by up to 1000-fold \[4\]. This demonstrates that APA's tunable dynamic range is not a fluke of a single parameter set but a general and robust feature of the dynamic assembly process. Therefore, the physical constraints of CPA assembly appear to be exquisitely tuned to provide a flexible and highly regulatable system for controlling gene expression.

### **1.3.4 Discussion** {#1.3.4-discussion}

This multi-scale modeling approach enables the system-level exploration of CPA dynamics under realistic biophysical constraints. By analytically resolving fast binding equilibria and numerically simulating elongation–termination kinetics, the framework captures cooperative CPA assembly, accounts for resource competition, and provides a quantitative basis for predicting termination efficiency and alternative polyadenylation site usage.

The complete model integrates two local and global scales. Locally, the probability of termination commitment increases with higher E factor occupancy on the CTD. Globally, as more polymerases sequester E factors, the pool of unbound Efree shrinks, which reduces binding probabilities at other sites and creates competition among polymerases. This setup allows the model to capture both local cooperative assembly and genome-wide resource competition simultaneously. The triangular scheme (Figure 1.3) illustrates this multi-scale integration.

The resulting triangular scheme is illustrated in Figure 1.3. Coupling across scales and emergent behaviors produce the model’s biologically relevant outcomes: Local coupling: The probability of commitment at the PAS increases with CTD occupancy ⟨n E ​⟩; thus local phosphorylation and pre-loading accelerate termination. Global coupling (competition): As polymerases across the genome pre-load E factors, they reduce 𝐸 𝑓 𝑟 𝑒 𝑒 E free ​ , lowering binding probabilities elsewhere and creating genome-wide competition for limited CPA factors. This competition can produce non-linear effects (cooperativity, thresholding) and affects both termination efficiency and PAS choice (APA outcomes).

Our results resolve the central paradox of CPA assembly: how the machinery assembles rapidly (10-20s) despite intense genome-wide competition for limited factors. We demonstrate that the Pol II CTD functions as a dynamic "pre-loading" scaffold. This scaffolding mechanism acts as a molecular conveyor belt, triggered by the progressive phosphorylation gradient. It allows the polymerase to pre-load or accumulate multiple factors onto its tail before it even transcribes the PAS, effectively concentrating the necessary components locally and ensuring they are poised for immediate, cooperative assembly the moment the PAS signal emerges. This mechanism represents the cell's elegant solution to the kinetic problem, mitigating the effects of genome-wide competition by ensuring that factors are on-site and ready for action.

Furthermore, our model reveals that gene architecture is not merely a passive structural feature but an active determinant of regulation. We predict a non-linear dependence of termination efficiency on gene length, suggesting that short genes may require distinct termination mechanisms. Additionally, our analysis of Alternative Polyadenylation (APA) shows that the typical 300 bp spacing between PAS sites results in an intermediate usage probability for the proximal site. This "Goldilocks" regime—neither 0% nor 100% usage—provides the optimal dynamic range for regulatory control, allowing cells to tune isoform selection via subtle shifts in factor concentration or PAS strength. 

*Limitations and Future Directions*  
While our model provides significant mechanistic insights, it relies on a one-dimensional representation of the gene, neglecting the three-dimensional organization of the nucleus, such as transcriptional bursting or the clustering of polymerases in "transcription factories," which could further influence local factor availability. Additionally, we simplified the complex CPA machinery into "Early" and "Late" modules. Future work will aim to incorporate these spatial complexities and a more granular view of the molecular assembly. Moreover, while this deterministic framework captures population-average behavior, gene expression is inherently stochastic. In Chapter 3, we will propose extending this work into the stochastic regime to explore how cell-to-cell variability in CPA dynamics contributes to phenotypic heterogeneity.

**Appendix: Parameter Estimation**

**A. First-Principles Estimation of the Hexamer Binding Rate**

Throughout this estimation, we use a basal biomolecular association rate of *k*on \= 0.5 × 106 M−1s−1. Brownian dynamics computer simulations of protein-protein association estimate diffusion-limited collision rates in a standard aqueous environment to be ∼106 M−1s−1 (Northrup & Erickson, 1992). To account for the increased viscosity and macromolecular crowding within the nucleoplasm, we apply a 1/2 reduction ratio to this theoretical limit, consistent with measurements of nucleoplasmic diffusion constraints (Bancaud et al., 2009; Dix & Verkman, 2008). To capture biological variability, we define the full possible range for this nuclear association rate as *k*on ∈ \[0.25 × 106, 2.5 × 106\] M−1s−1.

To accurately model the formation rate of the committed cleavage and polyadenylation (CPA) complex (*k*Hon), we decompose the assembly process into two sequential physical events:

1. **Reaction 1 (k1):** The tethered encounter between the early CPA factors (E-factors) pre-loaded on the Pol II CTD and the newly transcribed AAUAAA hexamer.

2. **Reaction 2 (k2):** The 3D diffusional recruitment of late CPA factors (L-factors, e.g., CstF) from the nucleoplasm to the downstream sequence element (DSE) and the E-factor cloud.

**1\. Reaction 1: Tethered Encounter Rate (k₁)**

Because the nascent RNA tether holding the hexamer is short (∼20 nt) and relatively stiff compared to the massive Pol II CTD (Chen et al., 2012; Seol et al., 2007), we model this as a one-body problem: the E-factor is tethered to a highly flexible CTD sweeping the space surrounding the RNA exit channel.

We parameterize the physical dimensions with a representative baseline and a physical range:

* **Effective Contour Length (Leff):** Baseline 69.15 nm (middle of the CTD). Range: \[35, 138\] nm (representing an E-factor binding between the first quarter and the very tip of the 52-repeat mammalian CTD, calculating length via 7 amino acids per repeat at 0.38 nm per residue) (Corden, 1990; Meinhart et al., 2005; Rubinstein & Colby, 2003).

* **Adjusted Persistence Length (lp):** Baseline 1.2 nm. Range: \[0.8, 1.5\] nm. Single-stranded nucleic acids and intrinsically disordered proteins (IDPs) typically exhibit short persistence lengths on the order of 0.6–1.5 nm, depending on ionic strength and sequence. Chen et al. (2012) report ionic-strength-dependent persistence lengths for ssRNA and ssDNA that fall within this range, and Hofmann et al. (2012) and Marsh & Forman-Kay (2010) demonstrate that IDPs exhibit analogous polymeric scaling behavior with short effective segment lengths. The chosen baseline *l*p \= 1.2 nm and range \[0.8, 1.5\] nm are therefore consistent with published polymer-physics characterizations of disordered chains in aqueous solution. The upper portion of the range additionally captures localized chain stiffening caused by dense Ser2 phosphorylation and electrostatic repulsion (Borg et al., 2007).

* **Anchor Distance (d):** Baseline 8 nm. Range: \[6, 12\] nm (representing the physical distance across the Pol II core from the RNA exit channel to the RPB1 linker connecting the CTD) (Bernecky et al., 2016; Cramer et al., 2001).

Using the 3D Gaussian probability density function for the effective local concentration (*J*density) at the RNA exit channel:

*J*density \= ( 3 / (2π⟨R²⟩) )3/2 × exp( −3*d*2 / (2⟨R²⟩) )

For our baseline (⟨R²⟩ \= 2 × 69.15 × 1.2 \= 165.96 nm², *d* \= 8 nm), this yields *J*density \= 8.64 × 10−5 molecules/nm³, converting to a baseline *J*single \= 0.143 mM.

**Range evaluation:** A compact cloud (*L*eff \= 35, *l*p \= 0.8, *d* \= 6\) yields a maximum *J*single ≈ 0.50 mM. A diffuse cloud (*L*eff \= 138, *l*p \= 1.5, *d* \= 12\) yields a minimum *J*single ≈ 0.04 mM. Range: \[0.04, 0.50\] mM.

Applying our baseline in vivo association rate (*k*bind \= 0.5 × 106 M−1s−1) and a steric penalty ϕ \= 0.1 to account for the bulky Pol II core:

*k*1, single \= ϕ · *k*bind · *J*single ≈ 7.15 s−1

**Range evaluation:** Utilizing the full boundaries of *k*bind and *J*single, the tethered encounter rate ranges from *k*1, single ∈ \[1.0, 125\] s−1.

**2\. Reaction 2: Diffusional Search Rate (k₂)**

L-factors must diffuse from the bulk nucleoplasm. We estimate the global nuclear copy number of late processing factors at approximately 100,000 molecules (Beck et al., 2011; Hein et al., 2015). Nuclear volume is established with a baseline of 374 μm³ (3.74 × 10−13 L), consistent with quantitative measurements of human cell nuclear dimensions (Maul & Deaven, 1977), and a biological range of \[300, 500\] μm³ for typical mammalian cell lines (Milo et al., 2010).

\[*L*free\] \= 100,000 / ( 6.022 × 1023 mol−1 × 3.74 × 10−13 L ) ≈ 0.444 μM

**Range evaluation:** Varying the nuclear volume yields \[*L*free\] ∈ \[0.33, 0.55\] μM.

The basal recruitment rate for a single target is:

*k*2, basal \= *k*on · \[*L*free\] \= (0.5 × 106 M−1s−1) · (0.444 × 10−6 M) ≈ 0.22 s−1

**Range evaluation:** *k*2, basal ∈ \[0.08, 1.38\] s−1.

**3\. Resolution of the CPA Assembly Window (kᴴᵒⁿ)**

The experimental window for CPA assembly is roughly 10 seconds (requiring *k*Hon ≥ 0.1 s−1). While a single E-factor achieves a fast *k*1, relying on a single factor is biologically insufficient for two reasons:

* **In Vivo RNA Competition:** The nascent RNA is fiercely contested by abundant RNA-binding proteins (e.g., hnRNPs) which coat the transcript immediately upon exit (Dreyfuss et al., 2002; Müller-McNicoll & Neugebauer, 2013). This competition drastically increases the effective hexamer off-rate (*k*−1). Multiple E-factors are required to overwhelm this competition via avidity and maintain a high pre-equilibrium bound state (*P*bound → 1).

* **The “Antenna Effect” (Scaling k₂):** The true kinetic bottleneck is *k*2\. However, L-factors maintain a high, cooperative affinity for the E-factor complex itself (Kumar et al., 2019; Mandel et al., 2008). Therefore, the 52-repeat CTD acts as a spatial antenna. Multiple E-factors (⟨*n*E⟩) linearly increase the multivalent collision cross-section for the L-factor:

*k*2, effective \= ⟨*n*E⟩ · *k*2, basal

If the CTD maintains an average of ⟨*n*E⟩ \= 5 E-factors, the effective bottleneck rate jumps to a baseline of ≈1.11 s−1 (Range: \[0.40, 6.90\] s−1), comfortably accommodating the empirical 10-second CPA assembly window even at the lowest boundary conditions.

**B. Estimation of the Microscopic Dissociation Rates**

**1\. E-factor – Hexamer Dissociation (kᴸ-ᴴᵉˣᵃᵐᵉʳ\_ᵒᶠᶠ)**

The equilibrium dissociation constant *K*d for recognition of the canonical AAUAAA signal by the CPSF-160–WDR33–CPSF-30 ternary complex has been measured by fluorescence polarization at approximately 3 to 50 nM (Clerici et al., 2017). We use this experimentally determined affinity alongside our nuclear association rate to derive the microscopic dissociation rate:

*k*E-Hexamer\_off \= *K*d · *k*bind \= (3 to 50 × 10−9 M) · (0.5 × 106 M−1s−1)

**Baseline calculation:** ≈ 0.0015 to 0.025 s−1

**Full theoretical range:** By applying the complete *k*bind bounds (\[0.25 × 106, 2.5 × 106\]), the absolute range is *k*E-Hexamer\_off ∈ \[0.00075, 0.125\] s−1.

**2\. L-factor Dissociation (kᴸᵒᶠᶠ)**

Pérez Cañadillas and Varani (2003) and Takagaki and Manley (1997) demonstrated that CstF-64 binds to the downstream GU-rich element quite weakly on its own, with a *K*d generally in the 10 to 100 μM range.

*k*Loff \= *K*d · *k*bind \= (10 to 100 × 10−6 M) · (0.5 × 106 M−1s−1)

**Baseline calculation:** ≈ 5 to 50 s−1

**Full theoretical range:** *k*Loff ∈ \[2.5, 250\] s−1.

**C. Estimation of kᴸᵒⁿ and kᴸᵒᶠᶠ**

***Estimation of the Early Factor Association Rate (kᴸᵒⁿ)***

The binding of Early factors to the CTD is assumed to be diffusion-limited. Using our established nuclear association rate (*k*on \= 0.5 × 106 M−1s−1), we convert this molar rate constant into a molecule-based rate constant suitable for stochastic simulations:

*k*Eon \= *k*on / ( NA · Vnucleus )

*k*Eon \= 0.5 × 106 / ( 6.022 × 1023 × 3.74 × 10−13 ) ≈ 2.22 × 10−6 molecule−1s−1

**Full theoretical range:** Utilizing the full volume range (300–500 μm³) and association rate range (0.25–2.5 × 106 M−1s−1), the absolute boundary limits are *k*Eon ∈ \[0.83 × 10−6, 13.8 × 10−6\] molecule−1s−1.

***Estimation of the Early Factor Dissociation Rate (kᴸᵒᶠᶠ)***

The intrinsic binding affinity (*K*d) of Early factor adaptors (e.g., the Pcf11 CID) to a single phosphorylated Ser2-Pro3 CTD repeat is extremely weak, consistently measured in vitro via NMR and ITC at 10 to 50 μM (Lunde et al., 2010; Meinhart & Cramer, 2004; Noble et al., 2005). Using the relationship *k*off \= *K*d · *k*on:

*k*Eoff, intrinsic \= (10 to 50 × 10−6 M) · (0.5 × 106 M−1s−1) \= 5 to 25 s−1

**Full theoretical range:** *k*Eoff, intrinsic ∈ \[2.5, 125\] s−1.

This substantial dissociation rate mathematically highlights that a single E-factor connection breaks almost instantaneously. Stable tethering of the CPA cloud on the elongating polymerase strictly requires multivalent avidity across the 52 CTD repeats.

## D. Estimation of the Maximum Number of Early Factor (E) Binding Sites on the CTD

### Overview

The maximum number of Early factor (E) molecules that can simultaneously occupy the RNA Pol II C-terminal domain (CTD) is governed by two competing constraints: **steric exclusion** imposed by the physical footprint of CPSF160 on the CTD chain, and **kinetic competition** between binding and unbinding at individual phosphorylated heptad repeat sites. To obtain a principled estimate, we implemented a Monte Carlo simulation in which the CTD is represented as a three-dimensional random-walk polymer and CPSF160 binding events are subject to explicit steric clash detection. The simulation was run for 1,000 independent trajectories and the maximum simultaneous occupancy was recorded for each run.

### Physical Model and Parameters

The following parameters were used, grounded in structural data for the human CTD and the CPSF160 subunit of the cleavage and polyadenylation specificity factor:

| Parameter | Value | Basis |
| :---- | :---- | :---- |
| Total CTD heptad repeats (N) | 52 | Human Pol II CTD consensus |
| Repeat rise per residue (ℓ) | 30 Å | 7 aa × \~3.8 Å/aa ≈ 26.6 Å; rounded to 30 Å |
| CPSF160 exclusion diameter (d) | 80 Å (8 nm) | Structural estimate of CPSF160 footprint on CTD |
| Number of consensus repeats | 26 (repeats 1–26) | High-affinity, pSer2/pSer5-bearing heptads |
| Number of non-consensus repeats | 26 (repeats 27–52) | Lower-affinity divergent heptads |
| Per-step binding probability (pᵇ) | 0.04 (consensus & non-consensus) | Uniform in current model |
| Per-step unbinding probability (pᵘ) | 0.04 | Symmetric on/off kinetics |
| Number of simulation runs | 1,000 | Statistical convergence criterion |

***Table D1.** Simulation parameters for CTD steric occupancy model.*

### CTD Chain Geometry and Steric Constraint

The CTD is modeled as a freely jointed chain (random walk in 3D) with N \= 52 nodes separated by a step length ℓ \= 30 Å. The total contour length of the chain is therefore:

Lc \= N · l \= 52 × 30 Å \= 1560 Å \= 156 nm  
The root-mean-square (RMS) end-to-end distance of the random walk is:

Rrms \= √N · l \= √52 × 30 Å ≈ 216 Å ≈ 21.6 nm  
Because the RMS end-to-end distance (\~21.6 nm) is substantially smaller than the contour length (156 nm), the CTD is highly compact and coiled in 3D space. A naive one-dimensional packing argument — dividing contour length by the CPSF160 exclusion diameter — provides an upper bound on the number of sterically accommodated binders:

Nmax \= Lcd \= 1560 Å80 Å ≈ 19  
This value of \~19 represents an absolute upper bound assuming every exclusion zone along the chain could be simultaneously filled without any chain folding bringing distant sites into proximity. In reality, the three-dimensional conformation of the CTD means that sequentially distant sites can overlap in space, making the effective packing capacity considerably lower. The simulation captures this three-dimensional geometry explicitly via Euclidean distance-based steric clash detection.

### Kinetic Occupancy Model

At each simulation step, a previously unphosphorylated site is phosphorylated at random, simulating progressive CTD phosphorylation during transcription elongation. Binding and unbinding events occur stochastically with per-step probabilities pb \= 0.04 and pu \= 0.04, respectively, giving a symmetric on/off ratio. In the absence of steric constraints, the expected steady-state fractional occupancy per phosphorylated site would be:

focc \= pbpb \+ pu \= 0.040.04 \+ 0.04 \= 0.5  
Without steric constraints, the expected number of simultaneously bound CPSF160 molecules would approach N × focc \= 52 × 0.5 \= 26\. However, the explicit hard-sphere steric exclusion enforced in the simulation (no two CPSF160 molecules may bind within 80 Å of each other in 3D space) dramatically reduces this number. Steric clash rejection prevents binding at sites that, while phosphorylated, fall within the exclusion radius of an already-occupied neighbor in the current chain conformation.

### Simulation Results

Across 1,000 independent simulation runs, the distribution of maximum simultaneous CPSF160 occupancy was as follows:

| Max. Simultaneous Bound | Count (n \= 1,000) | Frequency (%) |
| ----- | ----- | ----- |
| 2 | 17 | 1.7 |
| 3 | 129 | 12.9 |
| 4 | 396 | 39.6 |
| 5 | 319 | 31.9 |
| 6 | 118 | 11.8 |
| 7 | 18 | 1.8 |
| 8 | 3 | 0.3 |

***Table D2.** Distribution of maximum simultaneous CPSF160 binding across 1,000 simulation runs.*

The mean maximum simultaneous binding was **4.46 ± 0.99 (mean ± SD)**, with a median of 4 and a 95% confidence interval of \[4.40, 4.52\]. The distribution was unimodal and approximately symmetric around its mode of 4, with values ranging from 2 to 8 across all runs. The modal outcome of 4 simultaneously bound CPSF160 molecules was observed in 39.6% of simulations, and values of 4 or 5 together accounted for 71.5% of all outcomes. Values ≥6 occurred in only 13.9% of runs, and the maximum observed across all simulations was 8\.

### Interpretation and Comparison with Analytical Bounds

The simulation result of a mean maximum occupancy of \~4–5 CPSF160 molecules is substantially lower than both the kinetic ceiling (26, from unconstrained occupancy) and the one-dimensional linear packing ceiling (\~19). This large reduction arises from the interplay of two effects. First, the compact three-dimensional random-walk conformation of the CTD brings many sequentially distant nodes into close spatial proximity, so steric clash rejection is far more frequent than a linear-chain argument would suggest. Second, the stochastic nature of phosphorylation means that at any given moment only a fraction of sites are available for binding; the transient, dynamic nature of binding events (pb \= pu \= 0.04) further limits how many molecules can accumulate simultaneously before some unbind.

The summary of bounds is as follows:

* **Kinetic upper bound (no steric constraint):** \~26 (N × focc \= 52 × 0.5)

* **Geometric upper bound (linear packing):** \~19 (Lc / d \= 1560 Å / 80 Å)

* **Simulation result (steric \+ kinetic, 3D chain):** 4.46 ± 0.99 (mean ± SD; n \= 1,000 runs)

We therefore adopt a working estimate of **\~4–5 Early factor molecules** as the maximum simultaneous CTD occupancy consistent with the physical constraints of the model. This value provides a biologically plausible upper bound on cooperative E-complex engagement with the CTD and is used to constrain the stoichiometric parameters of the full kinetic model.

### Sensitivity and Limitations

The estimated maximum occupancy is sensitive to the assumed CPSF160 exclusion diameter and the CTD repeat length. Increasing the exclusion diameter from 80 Å to 120 Å would reduce the linear packing ceiling from \~19 to \~13 and is expected to further decrease the simulated mean. Conversely, reducing the exclusion diameter or increasing the binding probability would shift the distribution upward. The current model treats binding probability as identical between consensus (repeats 1–26) and non-consensus (repeats 27–52) heptads; introducing differential affinity in future iterations would enable exploration of preferential occupancy at the consensus region. Additionally, the random-walk chain model does not account for RNA polymerase-associated structural constraints on CTD conformation, which may further restrict or orient the available binding surface.