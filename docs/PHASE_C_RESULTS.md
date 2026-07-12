# PHASE C — Integration & Target Scorecard

**Project:** Novel drug targets in CD4⁺ T-cell genome-scale CRISPRi Perturb-seq (Marson–Pritchard atlas)
**Scope:** Integrate the six Phase B directions into one ranked, evidence-backed target list, split by therapeutic direction and novelty, and hand off the primary deliverable.
**Compute:** `ssh:clust1-rocm-4`, env `perturb-seq`. **Outputs:** `phaseBC_outputs/TARGET_SCORECARD.csv` (301 candidates × 35 columns) + `top_targets_overview.png`.

## Why this phase

Phase B produced six complementary but separate views of the same DE matrix. Phase C answers the question the whole pipeline was built for: **taking all six lines of evidence together, which genes are the strongest, most tractable, most disease-relevant drug-target hypotheses?** The integration is a transparent composite score that keeps every supporting sub-score visible, rewards multi-direction convergence, and separates the candidates by therapeutic logic (suppress inflammation vs boost immunity) and by novelty (drug-naive vs already-drugged).

## Headline

- **301 convergent candidate genes** nominated by ≥2 independent Phase B directions; **61 by ≥3 directions**; **3 by all four** discovery directions (*LAT, PLCG1, SENP5*).
- **Two mechanistic axes** structure the top of the list: the **TCR-proximal signalosome** (validated pathway, several drug-naive entry points) and a **chromatin / transcriptional co-activator axis** (SAGA/Mediator + remodelers, largely undrugged).
- Therapeutic split: **185 suppress-inflammation, 64 boost-immunity, 49 mixed** (+3 hub/context, direction-dependent).
- Novelty split: **178 novel-druggable, 27 known-drug-target, 96 difficult**.
- **27 candidates are already drug targets** (CD3E, ITK, MALT1, CD28, IL4R, IL2RB, STAT3, PTPRC…) — recovered *de novo* by the CRISPRi-loss-of-function logic, validating the pipeline end-to-end.

![**Fig C1. Integrated target scorecard overview.** Six-panel summary of the 301-candidate scorecard — composite score, direction convergence, therapeutic direction, novelty, and disease-genetics support.](../results/phaseBC/top_targets_overview.png)

---

## Composite scoring method

Each candidate is scored at its **strongest condition** on six normalized (0–1) sub-scores, weighted and summed, then multiplied by a convergence bonus:

| Sub-score | Weight | Source |
|---|---|---|
| Immune-program biology | 0.20 | D1 cytokine + D2 polarization strength |
| Effect strength | 0.15 | on-target + downstream effect magnitude |
| Reproducibility | 0.15 | Phase A guide/donor agreement |
| Disease genetics | 0.15 | D5 immune-GWAS support |
| Druggability | 0.15 | D6 tractability + novelty |
| Context-specificity | 0.10 | D3 activation-dependence |
| Network hub | 0.10 | D4 out-degree |

**Composite = (weighted sum) × (1 + 0.10 · (n_directions − 1))** — a +10% bonus per additional nominating direction, so a 4-direction gene gets ×1.3. All seven sub-score columns are retained in `TARGET_SCORECARD.csv` alongside the composite, so any ranking can be re-derived or re-weighted. Composite scores span **0.925 (max) → 0.549 (median) → 0.340 (min)** across the 301 candidates.

**Novelty & direction.** `novel-druggable` = no approved/clinical/known drug AND small-molecule-or-antibody tractable. Therapeutic direction is inherited from the Phase B sign convention: *suppress-inflammation* if knockdown lowers effector programs (inhibit the target for autoimmune/inflammatory indications), *boost-immunity* if knockdown raises them (inhibit the target for immuno-oncology).

## Convergence statistics

Of the **1,021** genes nominated by any single discovery direction, convergence thins sharply: **720 at 1 direction, 240 at 2, 58 at 3, 3 at all 4**. The **301 candidate priority set** (≥2 directions) is what carries external annotation and the composite score. Within the 301, direction membership is **288 cytokine · 159 context · 114 network · 105 polarization**, and by best condition **140 Stim48hr · 132 Stim8hr · 29 Rest** — i.e. **90% of top candidates act in stimulated cells**, consistent with the therapeutic thesis of targeting activated T cells. Context-class composition: 129 constitutive, 96 activation-induced, 65 late-activation, 9 mixed, 2 rest-specific. Module composition among the annotated hubs: **14 TCR-proximal, 7 SAGA/Mediator, 4 mitochondrial, 1 DNA-repair**.

---

## Top-25 integrated targets

Full table with all sub-scores and annotations: `TARGET_SCORECARD.csv` (301 genes × 35 columns).

| Rank | Gene | Best cond. | Score | #Dir | Direction | Novelty | Module | Immune GWAS | Downstream |
|---|---|---|---|---|---|---|---|---|---|
| 1 | **LAT** | Stim8hr | 0.925 | 4 | suppress | novel-druggable | TCR-proximal | – | 5,535 |
| 2 | **SMARCE1** | Stim8hr | 0.909 | 3 | suppress | novel-druggable | – | 37 | 3,266 |
| 3 | **PLCG1** | Stim8hr | 0.908 | 4 | suppress | novel-druggable | TCR-proximal | – | 5,032 |
| 4 | **CD247** | Stim8hr | 0.905 | 3 | suppress | novel-druggable | TCR-proximal | 46 | 4,329 |
| 5 | **STAT6** | Stim48hr | 0.885 | 3 | boost | novel-druggable | – | 55 | 1,127 |
| 6 | CD3E | Stim8hr | 0.881 | 3 | suppress | known-drug-target | TCR-proximal | 3 | 5,710 |
| 7 | **ZAP70** | Stim8hr | 0.836 | 3 | suppress | novel-druggable | TCR-proximal | 1 | 5,021 |
| 8 | IL4R | Stim48hr | 0.832 | 3 | boost | known-drug-target | – | 48 | 728 |
| 9 | **VAV1** | Stim8hr | 0.830 | 2 | suppress | novel-druggable | TCR-proximal | – | 4,897 |
| 10 | **TMX1** | Stim8hr | 0.800 | 3 | suppress | novel-druggable | TCR-proximal | 1 | 2,346 |
| 11 | CD2 | Stim48hr | 0.790 | 3 | suppress | known-drug-target | – | 2 | 988 |
| 12 | IL2RB | Stim48hr | 0.781 | 3 | suppress | known-drug-target | – | 11 | 1,069 |
| 13 | CD28 | Stim8hr | 0.773 | 3 | suppress | known-drug-target | – | 29 | 1,505 |
| 14 | **IL7R** | Stim8hr | 0.769 | 2 | suppress | novel-druggable | – | 44 | 811 |
| 15 | ITK | Stim48hr | 0.751 | 3 | suppress | known-drug-target | TCR-proximal | 5 | 2,565 |
| 16 | **EMC2** | Stim48hr | 0.750 | 3 | boost | novel-druggable | – | 3 | 465 |
| 17 | **PKM** | Stim8hr | 0.749 | 2 | suppress | novel-druggable | – | – | 1,840 |
| 18 | TBX21 | Stim8hr | 0.749 | 2 | suppress | difficult | – | 16 | 307 |
| 19 | **SIK3** | Stim48hr | 0.742 | 3 | suppress | novel-druggable | – | – | 1,809 |
| 20 | **EIF4G1** | Stim8hr | 0.732 | 3 | boost | novel-druggable | – | – | 1,678 |
| 21 | **CBLB** | Stim8hr | 0.732 | 2 | boost | novel-druggable | – | 11 | 1,027 |
| 22 | **TRIP12** | Stim8hr | 0.729 | 3 | suppress | novel-druggable | – | – | 1,829 |
| 23 | **HEXIM1** | Stim48hr | 0.727 | 3 | suppress | novel-druggable | – | 2 | 1,357 |
| 24 | **CARMIL2** | Stim8hr | 0.722 | 3 | suppress | novel-druggable | TCR-proximal | – | 2,544 |
| 25 | **NSD1** | Stim48hr | 0.718 | 3 | suppress | novel-druggable | – | 1 | 2,174 |

*(Bold = novel-druggable. "#Dir" = independent nominating directions; "Immune GWAS" = count of autoimmune/allergic disease-level associations; "Downstream" = significant downstream DE genes at best condition.)* Just outside the top 25: MYC (#26), MED1 (#27), SKIC8 (#28), CPSF6 (#29), and **SENP5 (#30)** — a 4-direction hit but flagged *difficult* on druggability.

---

## The two mechanistic axes

**Axis 1 — TCR-proximal signalosome (suppress-inflammation).** The three 4-direction hits (LAT, PLCG1, SENP5) and much of the top-20 (CD3E, CD247, ZAP70, VAV1, ITK, TMX1, CARMIL2) form the TCR-proximal signaling module (D4 Module 81). Knockdown broadly reduces effector cytokines (IL2RA, IFNG, IL2) and is strongly activation-induced. Several are **novel-druggable** — **LAT and PLCG1 are drug-naive but structurally/antibody-tractable**, offering fresh entry points into a clinically validated pathway that already yields drugs at CD3E (anti-CD3) and ITK.

**Axis 2 — chromatin & transcriptional co-activators (drug-naive).** The SAGA/Mediator complex (D4 Module 86: TADA1, TADA2B, TAF6L, SGF29, SUPT20H, MED12/MED24) and chromatin remodelers (SMARCB1, SMARCE1, NSD1, CHD4) converge on ≥3 directions, several with strong autoimmune GWAS signal — **SMARCE1 alone carries 37 disease associations at −log₁₀P ≈ 68**. Structurally tractable but almost entirely undrugged, this axis is the highest-novelty output of the screen and is what Phase E later validates as mechanistically distinct from TCR signaling.

## Highest-confidence novel candidates

Genes combining **≥3-direction convergence + novel-druggability + immune-GWAS support** are the strongest genuinely novel hypotheses — functional perturbation evidence, population genetics, and a tractable protein all at once:

| Gene | #Dir | Immune GWAS | Score | Axis |
|---|---|---|---|---|
| SMARCE1 | 3 | 37 | 0.909 | chromatin remodeler |
| CD247 | 3 | 46 | 0.905 | TCR-proximal |
| STAT6 | 3 | 55 | 0.885 | Th2 master TF |
| ZAP70 | 3 | 1 | 0.836 | TCR-proximal |
| TMX1 | 3 | 1 | 0.800 | TCR-proximal |
| SGF29 | 3 | 5 | 0.702 | SAGA/Mediator |
| MED24 | 3 | 12 | 0.698 | SAGA/Mediator |
| NSD1 | 3 | 1 | 0.718 | chromatin (methyltransferase) |
| ARNT | 3 | 1 | 0.677 | bHLH-PAS TF |

Plus EMC2, HEXIM1, MTG1. These are the candidates carried into Phase E single-cell validation and Phase F translational dossiers.

## De-novo recovery of known drug targets

That the composite ranking surfaces **27 established immune drug targets** without any drug-database seeding is the pipeline's end-to-end validation: CD3E (OKT3/teplizumab class), ITK, MALT1, CD28 (abatacept), IL4R (dupilumab), IL2RB, CD2, PTPRC, STAT3 all appear among the top candidates with known-drug annotations, recovered purely from CRISPRi-loss-of-function biology. These are drug-target *recoveries*, distinct from — though overlapping (ZAP70, CD28) — the formal Phase A positive-control set that benchmarked the DE calls upstream.

## Directionality caveat

The evidence layers are **complementary, not redundant**. The three 4-direction hits (LAT, PLCG1, SENP5) have **no** common-variant disease GWAS signal — their value is functional convergence, not population genetics. Conversely, the GWAS-strongest genes (IL7R −log₁₀P 98.5, IL2RA, STAT6, IL4R) are canonical immune loci that validate the genetics layer but rank lower on functional convergence. A candidate strong on both — SMARCE1, CD247, STAT6 — is therefore more compelling than one strong on either alone, which is precisely why the composite rewards convergence.

## Caveats

- **Composite weights are a design choice.** The 0.15/0.15/0.20/0.10/0.10/0.15/0.15 weighting and the +10%/direction bonus are transparent but not learned; all sub-scores are retained so the ranking can be re-weighted for a different therapeutic priority.
- **Scored at best condition.** Each gene is scored at its single strongest context, which favors strongly activation-induced genes; a gene weak at its best condition is not rescued by breadth.
- **Annotation set = 301.** D5/D6 external annotation ran on the ≥2-direction priority set, so the ~720 single-direction genes carry no genetics/druggability score and cannot enter the top ranks even if biologically real.
- **Pseudobulk nominations.** The scorecard is built entirely on the pseudobulk DE layer; single-cell confirmation of the flagship targets is Phase D (donor D4) and Phase E (cross-donor + novel-axis).

## Deliverables

| File | Content |
|---|---|
| **`TARGET_SCORECARD.csv`** | **Primary output** — 301 candidates × 35 columns: composite + 7 sub-scores + all annotations |
| `top_targets_overview.png` | 6-panel integrated summary figure (Fig C1) |
| `candidate_union.csv` / `_full.csv` | 301 (≥2 directions) / 1,021 (≥1 direction) convergent genes |
| `TARGET_REPORT.md` | Top-25 table + top-20 per-target dossiers + therapeutic hypotheses (companion narrative) |

**Downstream:** Phase D single-cell validation (7 flagship targets, donor D4) → Phase E cross-donor + novel-axis validation (19 targets, 4 donors) → Phase F translational dossiers (structure, chemical matter, tractability) for the 19 validated targets.

*Data source: CZI Virtual Cells "Primary Human CD4⁺ T Cell Perturb-seq" (Zhu, Dann, Yan, Reyes Retana, Goto, Guitche, Petersen, Ota, Shy, Pritchard, Marson; bioRxiv 2025.12.23.696273). Unpublished; MIT license.*
