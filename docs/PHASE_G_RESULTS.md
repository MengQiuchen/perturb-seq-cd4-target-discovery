# Phase G — First-tier single-cell analyses on nominated targets

**Question this phase answers.** Phases A–F nominated targets from *pseudobulk*
(population-averaged) readouts. Phase G asks the three questions that **only single-cell
resolution can answer** — the readouts that disappear under averaging — for the 19
already-nominated targets (7 flagship TCR-proximal + SMARCE1 constitutive control + 11
novel chromatin/TF/kinase axis).

**Design.**
- **Pooled 4 donors** (D1–D4) per condition to reach adequate per-target power. Pooled KD-cell
  counts: Rest floor 115 (SEL1L), Stim8hr floor 142 (SMARCE1)/144 (ZAP70),
  Stim48hr floor 52 (ZAP70). Targets with n<150 are flagged `underpowered` in the tables
  (ZAP70, SMARCE1, SEL1L at Stim8hr; more at Stim48hr) and never silently averaged in.
- **Reused, not reinvented.** Per-target response uses the *exact* Phase A/E DE_stats
  150-gene signatures (with log2FC sign weights); the activation program is the canonical
  E3 marker set (IL2RA, CD69, TNFRSF9, MKI67, TNFRSF18, ICOS). Expression is the same
  normalize_total(1e4)+log1p space; cells filtered `low_quality==False` as in Phase E.
- **Per-cell metrics.** *dose* = residual target expression ÷ NTC mean (0 = full KD, 1 = none);
  *response R* = NTC-standardized, sign-weighted mean over the target's signature genes (NTC ≈ 0).
- **Circularity guard (important).** Every target appears in its own DE signature, and a few
  activation markers appear in some signatures. Both are **excluded** from R before any
  dose→R or baseline→R correlation, so "stronger KD ↔ stronger response" and "activation state ↔
  response" cannot be mechanically forced. Removing 1–3 genes from ~147 left every result
  essentially unchanged, which is itself evidence the signals are real biology, not artefact.

---

## Result 1 — Single-cell knockdown dose–response (`G1_*`, `G_dose_response_*.csv`)

Within KD cells, cells with **stronger residual knockdown show a stronger downstream
signature** — a true dose–response that pseudobulk (one number per population) cannot see.

- **Stim8hr: 18/19 targets** show a significant *negative* Spearman ρ (median −0.226).
  Strongest: CD247 (ρ=−0.43, p≈1e-19), ZAP70 (−0.39), SMARCB1 (−0.37), CD3E (−0.34),
  SMARCE1 (−0.33).
- Only **PLCG1** is null (ρ=−0.04, p=0.43) — consistent with its lower target detection
  (frac_expr 0.71) making the per-cell dose axis noisy.
- Dose–response is **strongest under stimulation** (Stim8hr 18/19 sig vs Rest 12/19), i.e. the
  quantitative KD→program coupling is itself activation-gated.

## Result 2 — Response heterogeneity / responder fraction (`G2_*`, `G_responder_*.csv`)

How many KD cells actually manifest the signature, and is the response a uniform shift or an
all-or-none subset?

- Response is a **coherent population shift**: each KD population's R distribution moves right
  of the NTC null as a whole (Fig G2b), not a bimodal responder/non-responder split.
- **Responder fraction** (KD cells above NTC 95th pct) is high for the most penetrant targets —
  Stim8hr ZAP70 0.97, SMARCE1 0.96, CD3E 0.95, MED24 0.92 — and lowest for TRIP12 (0.35),
  STAT6 (0.58), CHD4 (0.63), marking targets whose KD effect is carried by a smaller cell fraction.

## Result 3 — Baseline-state sensitivity (`G3_*`, `G_baseline_sensitivity_*.csv`)  ← headline

Stratifying KD cells by their **baseline activation score** reveals a clean mechanistic split
that the flagship-vs-novel-axis claim predicted:

- **TCR-proximal KD effects are larger in *less*-activated cells** — Stim8hr mean ρ=−0.54
  (PLCG1 −0.68, CD247 −0.68, VAV1 −0.66, LAT −0.53); e.g. CD247 R = 1.48 in low-activation
  vs 0.23 in high-activation cells.
- **Chromatin/TF KD effects are baseline-independent** — Stim8hr mean ρ=−0.04.
- The separation is significant (Mann–Whitney TCR < chromatin/TF, **p=4e-4**) and **only
  emerges under stimulation** (Rest TCR ρ=−0.10). Signaling nodes gate the *rate* of a cell
  entering activation, so their loss bites hardest where activation has not yet saturated;
  chromatin/TF regulators set programs independent of instantaneous activation level.

---

## Cross-condition summary (`G0_cross_condition_summary.png`)

Rest → Stim8hr → Stim48hr for all three readouts, split by axis. All three peak or separate at
**Stim8hr**, reinforcing the report's finding that the acute (8 hr) window is where TCR-proximal
biology is most legible. The 48 hr point shows partial decay of the signaling-specific signals.

## Files

Per condition (Rest / Stim8hr / Stim48hr):
- `G1_dose_response_<cond>.png` + `G_dose_response_<cond>.csv`
- `G2_responder_<cond>.png` + `G_responder_<cond>.csv`
- `G3_baseline_<cond>.png` + `G_baseline_sensitivity_<cond>.csv`
- `G0_cross_condition_summary.png` — three-readout × three-condition synthesis

Script: `phaseE_scripts/phaseG_firsttier.py` (parameterized by condition).

## Caveats

- Underpowered target×condition cells (n<150, esp. ZAP70/SMARCE1/SEL1L at Stim8hr and several
  at Stim48hr) are flagged in the CSVs; treat their point estimates as indicative.
- Dose axis is only reliable where the target is detectably expressed in NTC (frac_expr ≥ 0.30);
  low-expression targets are hatched in G1 and marked `dose_reliable=False`.
- Baseline "activation state" is a 6-gene score, not a full trajectory position; it separates
  terciles robustly but is not a substitute for a pseudotime/Milo analysis (second-tier).
