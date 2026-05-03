# Comprehensive Synthesis: In Vivo Constraints on kHon

## 1. Definition of kHon in the Model

**kHon** is the effective first-order rate constant (s⁻¹) for PAS recognition / CPA commitment — the rate at which an elongating Pol II complex transitions from the elongating state (R) to the committed state (REH) upon transcribing the PAS. It represents the **tethered intramolecular encounter** between CTD-bound E-factors and the nascent AAUAAA hexamer.

The effective rate is:

$$k_{Hon} = \langle n_E \rangle \cdot k_{on} \cdot C_{eff} \cdot p_{gate}$$

where:
- **⟨nE⟩** = average number of E-factors pre-loaded on the CTD
- **k_on** = bimolecular association rate (~0.5–1 × 10⁶ M⁻¹s⁻¹ in nucleoplasm)
- **C_eff** = effective local concentration (J-factor) from CTD polymer physics
- **p_gate** (= ϕ in the model) = gating probability capturing steric occlusion, orientational constraints, conformational search, and RNA accessibility

---

## 2. Direct In Vivo Kinetic Constraints

### 2.1 Cis-Antisense Rescue Assay — Chao et al. (1999)
*Mol Cell Biol* 19(8): 5588–5600. PMID 10409748.

**The primary and most direct in vivo constraint on kHon.**

- Inverted PAS copy placed downstream blocks CPA via sense-antisense duplex formation
- Distance-dependent rescue: as antisense element moves further downstream, CPA commitment occurs before duplex forms
- **50% rescue distance for SV40 early PAS: ~200 bp**
- Fit to first-order kinetics: exponential parameter **k ≈ 4.9 × 10⁻³ bp⁻¹**
- Two distinct phases identified:
  1. **Protection phase** — CPA factors assemble around hexamer, blocking antisense access (= commitment, maps to kHon)
  2. **Cleavage phase** — actual endonucleolytic cut by CPSF73 (occurs after commitment)
- SV40 late PAS (strong): commits within **<100 bp** — several-fold faster than early PAS

**Converting distance to time:**
- At 40 nt/s (pre-PAS rate): 200 bp / 40 = **5 s** → kHon ≈ 0.14 s⁻¹
- At ~20 nt/s (post-PAS slowdown, Cortázar 2019): 200 bp / 20 = **10 s** → kHon ≈ 0.07 s⁻¹
- Total assembly completion: **10–20 s** → kHon ≈ 0.05–0.1 s⁻¹

**Quantitative constraint:**
| PAS | kHon estimate |
|---|---|
| SV40 early (moderate-strong) | **0.05–0.2 s⁻¹** |
| SV40 late (strong) | **≥ 0.2 s⁻¹** |

### 2.2 Genome-Wide 3′ Cleavage Half-Lives — Torres-Ulloa et al. (2024)
*PMC10870368*

- 4sU metabolic labeling + nascent RNA-seq + mathematical modeling in Drosophila S2 cells
- ~2857 constitutive cleavage sites: median **t₁/₂ ≈ 35.8 s** (k ≈ 0.019 s⁻¹)
- ~1601 alternative sites: median **t₁/₂ ≈ 42.1 s** (k ≈ 0.016 s⁻¹)
- Fast sites (e.g., PHGPx): **t₁/₂ ≈ 29 s** (k ≈ 0.024 s⁻¹)
- Slow sites (e.g., barc): **t₁/₂ ≈ 2.1 min**
- Alternative/upstream sites cleaved more slowly than constitutive/downstream sites
- Genes with multiple PASs show kinetic competition and delayed cleavage at individual sites

**Interpretation for kHon:** These half-lives measure the *entire* process (commitment + cleavage), not just the commitment step. Since commitment must precede cleavage, kHon must be **≥ these rates**. This provides a loose lower bound: **kHon ≳ 0.02 s⁻¹**.

### 2.3 Live-Cell Imaging Studies

**Liu et al. — MS2/PP7 dual-tagging, Drosophila embryos:**
- τ_cleave = **1.5–3 min** (hunchback reporter)
- Context-dependent; likely suboptimal PAS configuration
- Represents the "slow end" of CPA timing

**Darzacq et al. (2007) — HIV-1 live-cell reporter, mammalian cells:**
- **~55 s** pre-release step
- **~9 s** from polyadenylation onset to mRNA release
- Separates commitment-like timing from later release steps
- Consistent with commitment being faster than total 3'-end processing

**General MS2/PP7 studies:** Elongation at 2–4 kb/min with termination times of 15–60 s for reporter genes, consistent with kHon ~ 0.03–0.1 s⁻¹.

### 2.4 Pol II Post-PAS Slowdown — Cortázar et al. (2019)
- Pol II elongation rate drops from **>2 kb/min to <1 kb/min** downstream of PAS
- This ~2-fold slowdown directly affects the 200 bp → time conversion
- Supports kHon ≈ **0.07 s⁻¹** (using post-PAS speed) rather than 0.14 s⁻¹ (using pre-PAS speed)
- Already captured in the model by k_e2 = 30/100 (post-PAS elongation rate)

### 2.5 TT-seq Termination Zone — Schwalb et al. (2016)
*Science* 352(6290): 1225–1228.

- Median termination window: **~3,300 bp** downstream of PAS (K562 cells)
- Average of 4 termination sites per gene
- At 2–3 kb/min: total termination zone spans ~60–100 s
- This includes commitment + cleavage + XRN2 torpedo chase + PP1-mediated slowdown
- Commitment (kHon) occurs within the first few hundred bp; the rest is downstream termination machinery

### 2.6 CPSF73 as Gatekeeper — Eaton et al. (2020)
*Genes Dev* 34(1-2): 132–145. PMID 31805520.

- CPSF73 depletion → very extensive genome-wide readthrough (far beyond XRN2 depletion alone)
- Confirms two separable kinetic steps: commitment (kHon) and termination (XRN2)
- PP1 phosphatase decelerates Pol II downstream of PAS (dephosphorylates SPT5)
- Supports the model's kinetic separation of kHon from downstream termination rates

### 2.7 CPSF73 Cleavage Rate In Vitro
*RNA, 2021 (Martinson group)*

- Real-time fluorescence assay for CPSF73 endonuclease activity
- Cleavage rate: **~0.01–0.1 s⁻¹** in reconstituted systems
- Consistent with cleavage chemistry (kc) being rate-limiting *after* commitment
- Supports the model architecture where kHon and kc are separate parameters

---

## 3. CPSF–AAUAAA Binding Affinities

### 3.1 Equilibrium Dissociation Constants (Kd)

| Complex | RNA | Kd | Reference |
|---|---|---|---|
| CPSF160–WDR33–CPSF30–Fip1 (full core) | AAUAAA (16-mer) | **0.65 ± 0.09 nM** | Clerici et al. 2017 |
| CPSF160–WDR33–CPSF30 (ternary, no Fip1) | AAUAAA (17-mer) | **0.28 ± 0.07 nM** | Rios-Studer et al. 2019 |
| CPSF160–WDR33–CPSF30 (ternary) | AAUAAA (11-mer) | **0.32 ± 0.04 nM** | Rios-Studer et al. 2019 |
| CPSF160–WDR33–CPSF30–Fip1 | AAUAAA | **~3 nM** | Hamilton et al. 2019 |
| Reconstituted CPSF tetrameric | AAUAAA | **~2 nM** | Schönemann et al. 2014 |
| CPSF160–WDR33–CPSF30 (ternary) | AUUAAA | **~10 nM** | Rios-Studer et al. 2019 |
| CPSF160–WDR33–CPSF30 (ternary) | AUUAAA | **~17 nM** | Hamilton et al. 2019 |
| Core complex | AAGAAA | **120 ± 23 nM** | Clerici et al. 2017 |
| CPSF160–WDR33–CPSF30 (ternary) | AAGAAA | **≥ 50 nM** | Rios-Studer et al. 2019 |
| CPSF160–WDR33 (no CPSF30) | AAUAAA | **> 200 nM** | Clerici et al. 2017 |

**Key conclusions:**
- AAUAAA recognition is among the **highest-affinity** RNA–protein interactions known (0.3–3 nM)
- Fip1 inclusion has essentially no effect on RNA binding Kd
- Both WDR33 and CPSF30 are required; binary sub-complexes have Kd > 50 nM
- Structural basis: CPSF30 ZF2/ZF3 + WDR33 N-terminal basic motif jointly contact AAUAAA bases, including an unusual U3–A6 Hoogsteen base pair

### 3.2 Inferred Kinetic Rate Constants (kon, koff)

**Critical literature gap:** No published study has directly measured kon or koff for CPSF–AAUAAA by SPR, stopped-flow, or single-molecule kinetics. All estimates are inferred.

**Inferred koff values** (using kon = 0.5 × 10⁶ M⁻¹s⁻¹):

| PAS variant | Kd | Inferred koff | Residence time (1/koff) |
|---|---|---|---|
| AAUAAA (strong) | 0.3–1 nM | **1.5 × 10⁻⁴ – 5 × 10⁻⁴ s⁻¹** | ~30–100 min |
| AUUAAA (moderate) | ~10 nM | **~5 × 10⁻³ s⁻¹** | ~3 min |
| AAGAAA (weak) | ≥50 nM | **~2.5 × 10⁻² s⁻¹** | ~40 s |

**Key implication:** Once CPSF recognizes AAUAAA, the complex is **essentially irreversible** on the transcription timescale (seconds to minutes). For strong PAS, kHoff is negligible. For weak PAS variants (AUUAAA, AAGAAA), koff begins to compete with elongation, enabling APA regulation.

**Proxy measurement:** CstF2 (CstF64) kon has been directly measured at **4.35 × 10⁸ M⁻¹s⁻¹** — one of the few measured rates for any CPA factor. This is 2–3 orders of magnitude faster than the conservative nucleoplasmic estimate, suggesting electrostatic facilitation or fly-casting.

---

## 4. CTD Polymer Physics and J-Factor (Ceff)

### 4.1 CTD Physical Properties

| Property | Value | Source |
|---|---|---|
| Heptad repeats (human) | 52 (21 consensus + 31 degenerate) | Corden 1990 |
| Contour length | ~138 nm (364 aa × 0.38 nm/aa) | Calculated |
| Unphosphorylated CTD (compact) | ~10 nm (compact spiral, ~100 Å) | Structural data |
| Phosphorylated CTD (extended) | ~65 nm (~650 Å) | Structural data |
| Rg (Drosophila CTD, unphosphorylated, SAXS) | **~5–7 nm** | Portz et al. 2017 |
| Persistence length (smFRET/FCS, IDPs) | **0.4–0.5 nm** | smFRET literature |
| Persistence length (model current value) | **1.2 nm** (range 0.8–1.5) | Prospectus Appendix |
| Anchor distance (RNA exit channel to CTD linker) | **~8–10 nm** | Bernecky et al. 2016 cryo-EM |

**Discrepancy note:** The smFRET/FCS literature gives lp ≈ 0.4–0.5 nm, while the model uses 1.2 nm. If the shorter lp is correct, the chain is more compact, the J-factor is higher, but p_gate must be correspondingly larger to match observed kHon.

### 4.2 Effective Concentration (Ceff / J-Factor) Estimates

**Fuertes et al. (2019)** *PNAS* 116(46): 23124–23131:
- Experimental measurement of Ceff for IDR linkers of varying length
- **~100 residue GS-linker: Ceff ~ 0.1–10 mM**
- Scaling: Ceff ~ N^(−3ν), where ν = Flory exponent (0.4–0.6)
- Phosphorylated CTD (more extended, higher ν) → **lower Ceff** than unphosphorylated

**Sørensen & Kjaergaard (2019):**
- Disordered linker Ceff follows polymer-physics scaling
- Sequence-dependent compaction can shift Ceff by **~63-fold** for a 100-residue linker

**Freiburger et al. (2020)** *PNAS* 117(34): 20551–20562:
- Tethered kinase system validates that polymer-physics Ceff predicts actual intramolecular reaction rates
- For Ceff ~ 100–1,000 μM: 2–10-fold acceleration of tethered reactions

**Geometric estimates (single molecule in confined sphere):**

| Sphere radius | Ceff |
|---|---|
| 5 nm | ~3.2 mM |
| 10 nm | ~0.4 mM |
| 15 nm | ~0.12 mM |

**CTD length-dependent estimates** (using kon = 5 × 10⁵ M⁻¹s⁻¹):

| CTD segment | Ceff estimate | k_encounter (per E-factor) |
|---|---|---|
| ~4 heptads (30 residues) | ~10–50 mM | 5,000–25,000 s⁻¹ |
| ~14 heptads (100 residues) | ~0.5–5 mM | 250–2,500 s⁻¹ |
| ~28 heptads (200 residues) | ~0.05–0.5 mM | 25–250 s⁻¹ |
| Full CTD (365 residues) | ~10–100 μM | 5–50 s⁻¹ |

**Key conclusion:** The tethered encounter rate per E-factor (50–5,000 s⁻¹) is **orders of magnitude faster** than the observed commitment rate (~0.07 s⁻¹). This confirms that kHon is **not limited by the physical encounter** but by downstream gating factors (p_gate / ϕ).

### 4.3 LLPS Condensate Enhancement
- Pol II CTD phosphorylation promotes partitioning into liquid-liquid phase-separated transcriptional condensates
- Local concentration of RNA-binding proteins within condensates can reach **~9.6 mM** (SRSF1 estimate)
- Could further enhance Ceff but is difficult to incorporate quantitatively

---

## 5. The Gating Probability (p_gate / ϕ)

### 5.1 Definition and Physical Meaning

In the model: **kHon = ϕ · k_bind · J_single** (per E-factor), where ϕ is currently set to **0.1**.

The gating probability captures multiple effects that reduce the productive encounter rate below the diffusion limit:

1. **Steric occlusion by Pol II core** — RNA exit channel and CPSF binding pocket partially blocked
2. **Orientational constraints** — CPSF30 ZF2/ZF3 must adopt specific orientation relative to AAUAAA bases
3. **Conformational search cost** — disordered CTD must adopt the right conformation; nascent RNA must be accessible
4. **Competition from hnRNP proteins** — abundant RNA-binding proteins (hnRNP A1, hnRNP H) coat nascent transcript, occluding the hexamer
5. **Hexamer accessibility window** — the hexamer transits the "recognition window" near the exit channel in ~0.2 s out of a longer elongation cycle

### 5.2 Estimated Values

**From the gap between theory and observation:**
- k_encounter (theoretical, per E-factor) ≈ 50–5,000 s⁻¹
- kHon (observed, effective) ≈ 0.05–0.1 s⁻¹
- With ⟨nE⟩ = 5: kHon_theoretical = 250–25,000 s⁻¹
- **Implied p_gate ≈ 10⁻⁵ to 10⁻³**

However, much of this "gap" is already captured by other model parameters (⟨nE⟩ computed self-consistently, E-factor loading statistics). The residual gap attributable to ϕ specifically depends on which effects are modeled elsewhere.

**From the conservative Ceff = 1 μM estimate:**
- kHon_theoretical = 10⁶ × 10⁻⁶ × 2.5 = 2.5 s⁻¹
- kHon_observed ≈ 0.1 s⁻¹
- **Implied ϕ ≈ 0.04** — close to the current model value of 0.1

**From classical Brownian dynamics literature:**
- Northrup & Erickson (1992): κ ≈ 10⁻²–10⁻¹ for protein–protein recognition with geometric constraints
- Rogers et al. (2013): IDP binding at 10⁵–10⁶ M⁻¹s⁻¹ is not purely diffusion-limited; conformational search adds a penalty
- **Literature range for steric/orientational factors: ϕ ~ 0.01–0.1**

### 5.3 Key Literature Gap
No direct measurement of the steric/orientational penalty for a CTD-tethered protein encountering nascent RNA exists. Current estimates are all inferred.

---

## 6. Recommended Parameter Ranges

### 6.1 kHon

| PAS type | Example | Kd (CPSF) | Recommended kHon | kHoff |
|---|---|---|---|---|
| Strong (canonical) | AAUAAA, SV40 early | 0.3–1 nM | **0.1–0.2 s⁻¹** | ~10⁻³–10⁻⁴ s⁻¹ |
| Nominal single value | — | — | **0.07 s⁻¹** | — |
| Moderate | AUUAAA | ~10 nM | **0.01–0.05 s⁻¹** | ~10⁻² s⁻¹ |
| Weak | AAGAAA, AAUAAG | ≥50 nM | **<0.01 s⁻¹** | ~10⁻¹ s⁻¹ |
| Broad modeling prior | All sites | — | **0.02–0.2 s⁻¹** | — |

### 6.2 Gating Probability (ϕ / p_gate)

| Estimate source | ϕ value |
|---|---|
| Current model | 0.1 |
| Implied from Ceff = 1 μM | ~0.04 |
| Brownian dynamics literature | 0.01–0.1 |
| **Recommended range** | **0.01–0.1** |

---

## 7. Key Literature Gaps

1. **Direct kon/koff for CPSF–AAUAAA** — No SPR, stopped-flow, or smFRET measurements exist. All are inferred from Kd + diffusion limits.
2. **CTD-specific J-factor / Ceff** — No direct measurement for a CTD-tethered CPA factor encountering nascent RNA.
3. **Direct p_gate / ϕ measurement** — No experimental or computational estimate of the steric/orientational penalty for this specific system.
4. **Single-molecule CPSF dwell times at transcription sites** — FRAP/FCS of CPSF subunits at mammalian transcription sites has not produced published kon/koff values.
5. **Cis-antisense rescue for PAS variants** — Only SV40 early/late and a few synthetic sites characterized. Systematic mapping of kHon across PAS strength spectrum needed.

---

## 8. Complete Reference List

### Primary kHon constraint
1. **Chao LC, Jamil A, Kim SJ, Huang L, Martinson HG** (1999). Assembly of the cleavage and polyadenylation apparatus requires about 10 seconds in vivo and is faster for strong than for weak poly(A) sites. *Mol Cell Biol* 19(8): 5588–5600. PMID 10409748.

### Genome-wide cleavage kinetics
2. **Torres-Ulloa et al.** (2024). Genome-wide kinetic profiling of pre-mRNA 3′ end cleavage. *PMC10870368*.

### Live-cell imaging
3. **Liu et al.** Real-time single-cell characterization of the eukaryotic transcription cycle. *Garcia Lab, Berkeley*.
4. **Darzacq et al.** (2007). Live-cell mammalian HIV-1 reporter — 3'-end processing/release timing.

### CPSF–AAUAAA binding affinities
5. **Clerici M, Faini M, Aebersold R, Jinek M** (2017). Structural insights into the assembly and polyA signal recognition mechanism of the human CPSF complex. *eLife* 6: e33111. PMID 29274231. [Kd = 0.65 nM]
6. **Clerici M et al.** (2018). Structural basis of AAUAAA polyadenylation signal recognition by the human CPSF complex. *Nat Struct Mol Biol* 25: 135–138.
7. **Rios-Studer K et al.** (2019). Biophysical characterizations of the recognition of the AAUAAA polyadenylation signal. *RNA* 25(12): 1673–1685. PMID 31462423. [Kd = 0.28 nM; AUUAAA ~10 nM]
8. **Hamilton K, Sun Y, Tong L** (2019). Biophysical characterizations of the recognition of the AAUAAA polyadenylation signal. PMID 31462423.
9. **Schönemann L et al.** (2014). Reconstitution of CPSF active in polyadenylation. *Genes Dev* 28: 2381–2393.
10. **Chan SL et al.** (2014). CPSF30 and Wdr33 directly bind to AAUAAA in mammalian mRNA 3' processing. *Genes Dev* 28: 2370–2380.

### Termination mechanism
11. **Eaton JD, Francis L, Davidson L, West S** (2020). A unified allosteric/torpedo mechanism for transcriptional termination on human protein-coding genes. *Genes Dev* 34(1-2): 132–145. PMID 31805520.
12. **Schwalb B et al.** (2016). TT-seq maps the human transient transcriptome. *Science* 352(6290): 1225–1228.
13. **Cortázar et al.** (2019). Pol II speed near PAS — elongation slowdown downstream of PAS.

### CTD polymer physics
14. **Portz B et al.** (2017). Structural heterogeneity in the intrinsically disordered RNA polymerase II C-terminal domain. *Nat Commun* 8: 15231.
15. **Corden JL** (2013). RNA Polymerase II C-Terminal Domain: Tethering Transcription to Transcript and Template. *Chem Rev*.

### Effective concentration / J-factor
16. **Fuertes G, Iglesias-Bexiga M et al.** (2019). Effective concentrations enforced by intrinsically disordered linkers are governed by polymer physics. *PNAS* 116(46): 23124–23131.
17. **Freiburger L, Fischbach C, Bhatt B et al.** (2020). Intrinsically disordered linkers control tethered kinases via effective concentration. *PNAS* 117(34): 20551–20562.
18. **Sørensen CS, Kjaergaard M** (2019). Effective concentrations enforced by intrinsically disordered linkers. [Polymer-physics scaling of Ceff]

### Association kinetics
19. **Northrup SH, Erickson HP** (1992). Kinetics of protein-protein association explained by Brownian dynamics computer simulation. *PNAS*.
20. **Rogers JM et al.** (2013). IDP binding kinetics — fast binding need not be purely diffusion-limited.
21. **Bancaud A et al.** (2009); **Dix JA, Verkman AS** (2008). Nucleoplasmic diffusion constraints.

### Proteomics
22. **Beck M et al.** (2011); **Hein MY et al.** (2015). Quantitative proteomics of nuclear proteins.
