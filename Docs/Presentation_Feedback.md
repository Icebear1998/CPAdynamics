# Presentation Feedback Analysis — CPA Project

> Extracted from live-recorded comments ([Comments2.md](file:///Users/kienphan/WorkingSpace/ResearchProjs/Prospectus/Comments2.md)) on the [Prospectus Presentation](file:///Users/kienphan/WorkingSpace/ResearchProjs/Prospectus/Prospectus_Presentation.pptx).

---

## 1. Presentation Structure & Narrative Flow

These are the most impactful comments — they call for **reordering slides** so that the biology motivates the model, not the other way around.

| # | Feedback | Action |
|---|----------|--------|
| 1.1 | **State a "Central Hypothesis" explicitly.** Advisor noted you're already doing it implicitly ("You're just not explicitly saying *here is our central hypothesis*—you're just missing the word"). | Add a slide (or prominent callout) that says: *"Central Hypothesis: Multiple E-factor binding sites on the CTD are functionally required to enable efficient CPA assembly within the 400–800 bp window."* |
| 1.2 | **Key research questions should come *before* the model, not after.** ("Shouldn't this come before you introduce the model?" / "You don't raise questions after you already have the methods — doesn't logically make sense.") | Move the "Key Research Questions" slide (currently slide 15) to **before** the model introduction slides (slides 13–14). The flow should be: Biology → Questions → Model → Results. |
| 1.3 | **Introduce APA *before* raising APA research question (Q3).** You mentioned APA only briefly on one early slide and then the question appeared much later, surprising the audience. ("I was like, hey, where does that come from?") | Move the APA introduction to **right after** the first two CPA questions are posed, then raise Q3 about APA regulation immediately. This keeps the introduction close to the question so people don't forget. |
| 1.4 | **Frame it as "Central Hypothesis → Key Questions → Model."** Advisor suggested: list the central hypothesis, then the questions as aims, then say *"To validate/contradict this hypothesis, we build this model."* | Restructure the intro section to follow this arc explicitly. |
| 1.5 | **Include a cartoon/schematic alongside the research questions slide** so the audience can connect each question to a visual element of the biological process. ("If you have a graphical depiction side by side with these questions, then you can actually connect your questions to what it means in your model.") | Add the simplified CPA cartoon (showing Pol II, CTD, E-factors, PAS) next to or on the same slide as the research questions. |

### Suggested New Slide Order (Intro Section)

```
1. Why CPA matters (slides 1–2) — keep
2. How CPA works + CPA factors (slides 3–5) — keep
3. Physical constraints / competition (slide 6) — keep
4. Simplified CPA assembly summary (slides 7–8) — keep
5. CTD as dynamic platform + multiple E-binding (slides 10–11) — keep
6. Enhancement vs. Competition — the tension (slide 12) — keep
7. APA introduction (slide 4, moved/expanded) — MOVE here
8. ★ NEW: Central Hypothesis + Key Research Questions (with cartoon)
9. Model description (slides 13–14) — "To address these, we build..."
10. Results (slides 16–22)
```

---

## 2. Slide Titles — Make Them Statements

> *"Make the slide titles a statement… think of your slide titles as an outline of your whole presentation — if you read just the titles, it should read almost like an abstract."*

| Current Title | Suggested Revision |
|---|---|
| *"The Pol II CTD: A Dynamic Assembly Platform"* | **"Pol II CTD potentially provides multiple binding sites for CPA factors"** |
| *"CPA Assembly: Enhancement vs. Competition"* | **"Multiple E-binding may enhance local assembly but worsen genome-wide competition"** |
| *"Multi-Scale Model of CPA Assembly"* | **"A kinetic model captures CPA dynamics across fast binding and slow elongation"** |
| *"Key Research Questions"* | **"Central hypothesis and key questions"** |
| *"Quantifying CPA Assembly: The CAD Metric"* | **"CPA Assembly Distance (CAD₅₀) quantifies efficiency against experimental benchmarks"** |
| *"Result: CTD Scaffolding Enables Efficient CPA"* | **"Multi-site E pre-loading is sufficient to explain efficient CPA assembly"** |
| *"Result: CPA Efficiency Depends on Gene Length"* | **"Short genes lack sufficient runway for factor accumulation"** |
| *"Result: Mechanistic Basis for APA Regulation"* | **"The ~300 bp inter-PAS spacing maximizes APA regulatory flexibility"** |

---

## 3. Figure & Visual Improvements

| # | Feedback | Action |
|---|----------|--------|
| 3.1 | **Replace the music-note "pause" symbol** with a traffic-sign "slow" symbol. ("Not everyone recognizes them… pause literally means not moving, but you mean *slower*.") | Find/create a speed-limit or "slow" traffic sign icon for the elongation-rate change after PAS. |
| 3.2 | **Refine the main kinetic schematic (Figure 1.3).** Use a **triangular regime** layout: pre-PAS triangle (R states) and post-PAS section, with color-coded shading. Specifically, use **different shades of blue** for increasing E-binding number (white → light blue → dark blue) to reflect increasing exit flux. ("Almost like a gradient… the darker state would exit faster.") | Redesign Figure 1.3 with gradient-shaded triangles. |
| 3.3 | **Add a reminder cartoon on the E-binding profile result slide** (the average E vs. gene position plot). ("You may want to include a cartoon just reminding people what you're talking about.") | Insert a small inset cartoon showing CTD + E-factors next to the E-binding profile figure. |
| 3.4 | **On the inter-PAS distribution slide**, report the spread (e.g., "300 ± σ bp") alongside the peak value, since the distribution is right-skewed. | Add standard deviation or IQR to the reported 300 bp peak. |
| 3.5 | **Change "Multi-Scale" to "Kinetic"** in model title/description. ("Multi-scale… from *what*? It's not exactly clear. *Kinetic model* would be more accurate.") | Update slide 13 title and any references from "Multi-Scale Model" → "Kinetic Model of CPA Assembly." |

---

## 4. Technical / Scientific Clarifications Needed

| # | Issue Raised | What to Address |
|---|-------------|-----------------|
| 4.1 | **Time-scale separation between R and REH triangles may not hold.** Advisor pointed out that if you multiply KE_on × [E], it may be comparable to KH_on — so there's no clear separation. ("At least compare KE_on × [E] with KH_on… the first should be significantly faster.") | **Action:** Compute KE_on × [E_concentration] vs. KH_on for your parameter values. If they're comparable, either (a) justify the approximation or (b) implement the two-triangle coupling. At minimum, add a slide/footnote acknowledging this. |
| 4.2 | **The CAD metric should be labeled CAD₅₀** to remind the audience it's the 50% assembly distance, not 100%. ("You maybe want to name it as CAD subscript 50 to just remind people.") | Rename throughout slides and labels. |
| 4.3 | **Clarify what "n" (max E-binding sites) means on the result plot.** Labmate was confused — advisor suggested labeling it as *"Maximum E-factor binding sites on CTD"* rather than just "n". | Update axis label and legend on the CAD vs. n plot (slide 17). |
| 4.4 | **Clarify the deterministic model's relationship to stochasticity.** A labmate asked about fluctuations; the model is deterministic but estimates probabilities (the CAD₅₀ is a probability-based metric). Advisor: "The deterministic model essentially already has the stochasticity kind of built into it." | Add a brief note on the result slide explaining: *"Deterministic model; CAD₅₀ represents the distance at which 50% of polymerases have committed — fluctuations are implicit in the probability metric."* |
| 4.5 | **Explain the E-factor crossover** (n = 5 dips below n = 3 near the TSS, then surpasses it near PAS). Labmates asked for a physical explanation. Related to phosphorylation-dependent binding protection from dephosphorylation. | Prepare a clear explanation or supporting figure for why higher n means lower early accumulation but higher accumulation near PAS. |
| 4.6 | **Clarify what "conservative estimate" means.** Multiple people asked how you know the estimates are conservative, and whether conclusions change with different parameter values. | See Section 5 below — this feeds directly into the sensitivity analysis request. |

---

## 5. Parameter Sensitivity Analysis (New Work Required)

> **This was the single most emphasized action item. Both the advisor and labmates stressed it, and the advisor noted it was also raised at the committee meeting.**

| # | Task | Details |
|---|------|---------|
| 5.1 | **Perform parameter sensitivity analysis for the n = 1 case.** Show that even when varying parameters over their literature ranges, n = 1 still cannot achieve CAD within 400–800 bp. | This is the **foundation** of the argument — without it, "the whole foundation of this work is broken." (advisor's words) |
| 5.2 | **Rank parameters by sensitivity.** Identify which parameters the CAD is most sensitive to, and which are insensitive (where the conservative estimate doesn't matter). | Use local sensitivity analysis (one-at-a-time) or a multi-dimensional approach. |
| 5.3 | **Present it as supporting evidence** for the n = 1 result. The logic: (1) show conservative estimate → CAD too high for n = 1; (2) show sensitivity analysis → conclusion is robust; (3) *then* show n = 3, 5 bring CAD into range. | Add 1–2 slides after the first CAD result. |

---

## 6. Minor Slide-Level Fixes

| Slide | Fix |
|-------|-----|
| Enhancement vs. Competition (slide 12) | **Put competing effects on the same slide** — "locally enhances" vs. "globally worsens" — and make it visually clear (e.g., a balance/scale graphic). |
| Model overview (slide 13) | State that parameters are **conservative, literature-based estimates** (not fitted). Currently this is said verbally but not on the slide. |
| CAD result (slide 17) | Add a shaded band showing the 400–800 bp experimental range on the plot itself (already present, but make more prominent + add CAD₅₀ label). |
| E-binding profile (slide 18) | Specify that negative positions = before PAS, positive = after PAS. Some audience members were confused about the x-axis convention. |
| APA result — proximal PAS usage (slide 21) | The text says ~20%, but verbally you said ~40%. **Clarify and reconcile** the actual number. |

---

## 7. Follow-Up Research Tasks (From Discussion)

These emerged from the Q&A and are more research-level than presentation-level, but worth tracking:

| # | Task | Priority |
|---|------|----------|
| 7.1 | Verify time-scale separation: compute KE_on × [E] vs. KH_on for actual parameter values. | **High** — affects model validity |
| 7.2 | Parameter sensitivity analysis for n = 1 (CAD vs. parameter ranges). | **High** — critical for defense of results |
| 7.3 | Investigate whether phosphorylation-dependent E-binding protection from dephosphorylation explains the crossover. | Medium |
| 7.4 | Explore genomic data for correlation between inter-PAS distance and PAS strength. | Medium — future direction |
| 7.5 | Consider modeling without time-scale separation between R and REH triangles (two-triangle coupling). | Low — only if 7.1 shows separation doesn't hold |

---

## Summary of Top Priorities

1. **Restructure intro slides**: Biology → Central Hypothesis → Questions (with cartoon) → Model
2. **Make slide titles into statements** (abstract-like outline)
3. **Parameter sensitivity analysis for n = 1** — the most critical new work
4. **Verify time-scale separation** (KE_on × [E] vs. KH_on)
5. **Redesign kinetic schematic** with gradient-shaded triangles
6. **Rename CAD → CAD₅₀** and clarify axis labels
7. **Reconcile the proximal PAS usage number** (20% vs. 40%)
