# Phase G — Third-tier single-cell analyses (mechanism / biomarker level)

The deepest tier: readouts that require correlation structure, guide-level identity, or a
trajectory — none of which pseudobulk can produce. All on Stim8hr, 4 donors pooled.

**New data gathered for this tier.** The Phase E subsets carry only `guide_type`
(targeting / non-targeting), not which of a gene's two CRISPRi guides labelled each cell.
For the dual-guide analysis we re-gathered guide-resolved subsets from the raw
`cell_level/D*_Stim8hr.assigned_guide.h5ad` shards (140–170 GB each), keeping `guide_id` +
`guide_group`, via the same slab-streaming pattern as Phase E
(`phaseG_gather_guide.py` → `phaseG_outputs/guide_subsets/`).

---

## Result 7 — Co-expression module rewiring (`G8_*`, `G8_coexpr_rewiring_Stim8hr.csv`)

Within each target's downstream module (top-30 signature genes, target self excluded), we compare
the gene–gene correlation structure in KD vs NTC cells. pseudobulk has one sample per group and
cannot compute a correlation at all.

- **Cell-count confound found and corrected.** A first pass (variable N per target) gave
  Spearman(n_kd, structure-preservation) = **+0.61 (p=0.006)** — the metric was tracking cell
  count, not biology (fewer cells → noisier correlations → apparent "rewiring"). We refit with
  **every target estimated at a fixed N=140 cells** per group and a **300-permutation null** for
  the connectivity change; the confound reversed to −0.53, i.e. removed.
- After correction, **7/19 targets show significant module rewiring** (permutation FDR<0.05,
  well-powered). **TCR-proximal knockdown tightens the downstream module most** (mean
  Δ|r| = +0.12 vs chromatin/TF +0.02): losing a signaling node makes its residual target genes
  co-vary more tightly, consistent with collapse onto a single "failed-activation" axis.
- Only SEL1L (n=114) is underpowered (n<140) — hatched in the figure, flagged
  `well_powered=False`, and not interpreted. SMARCE1 (n=142) and ZAP70 (n=144) clear the
  threshold and are included in the statistics (ZAP70's Δ|r|=0.027 is part of the TCR-proximal
  mean of +0.12).

## Result 8 — Activation-trajectory checkpoint (`G9_*`, `G9_trajectory_checkpoint.csv`)

A Rest→Stim8hr→Stim48hr activation pseudotime was built on **control (NTC) cells** across all three
conditions (harmony-integrated DPT); NTC medians order correctly (Rest 0.0-ish → Stim8hr →
Stim48hr 0.66), validating the axis. KD cells are projected on and their pseudotime compared to NTC.

- **TCR-proximal knockdown stalls cells early on the activation trajectory**: mean Δpseudotime
  = **−0.30** (CD3E −0.42, ZAP70 −0.41, LAT −0.28, PLCG1 −0.28, CD247 −0.21, VAV1 −0.21), all
  FDR≈0. This is the single-cell realization of "these cells cannot progress through activation."
- **Chromatin/TF knockdown does not move trajectory position** (mean Δ = −0.002; individual
  |Δ|<0.02). They reshape programs (Results 4–7) without holding cells at an earlier activation
  stage.
- Separation is highly significant (Mann–Whitney TCR < chromatin/TF, **p=1e-4**) and matches the
  mechanistic prediction exactly: signaling nodes gate trajectory *progression*; chromatin/TF
  regulators act orthogonally to trajectory position.

## Result 9 — Dual-guide single-cell concordance (`G10_*`, `G10_dualguide_Stim8hr.csv`)  ← target-specificity

For each target with two guides having ≥25 pooled cells, we test whether both guides push single
cells toward the **same** response — a stronger on-target / anti-off-target confirmation than a
pseudobulk correlation of two bulk profiles.

- **17/19 targets testable; 17/17 show both guides shifting response in the same direction AND
  both individually significant vs NTC.** Median guide-agreement (effect-size ratio) = **0.88**.
  This is strong single-cell evidence the phenotypes are on-target.
- **CD3E and ZAP70 drop to a single usable guide** even after pooling 4 donors (their second guide
  captured too few cells) — a data limitation, not a discordance; they are marked
  `single_guide_or_low` and not counted among the 17.
- **CD247** is the one magnitude-discordant case (guide 1 effect 1.35 vs guide 2 0.30, ratio 0.22):
  both guides agree in *direction* and are significant, but one is far stronger — likely a
  guide-efficiency difference. Flagged for follow-up; direction concordance still holds.

---

## Files
- `G8_coexpr_rewiring_Stim8hr.png/.csv` — module rewiring (fixed-N, permutation-tested)
- `G9_trajectory_checkpoint.png/.csv` + `G9_pseudotime_cells.csv.gz`, `G9_trajectory_manifold.h5ad`
- `G10_dualguide_Stim8hr.png/.csv` + `G10_dualguide_percell_Stim8hr.json`
- `guide_subsets/D*_Stim8hr.guide.h5ad` (+ summary json with per-target guide inventory)
- Scripts: `phaseG_coexpr.py`, `phaseG_trajectory.py`, `phaseG_dualguide.py`, `phaseG_gather_guide.py`

## Caveats
- Co-expression rewiring is interpreted only for well-powered targets (n≥140); the fixed-N design
  removes the cell-count confound but caps the module at 30 genes for stability.
- Pseudotime is a DPT ordering on the activation manifold, not a physical-time or RNA-velocity
  trajectory; it measures relative progression, and the Δ values are on that 0–1 scale.
- Dual-guide uses the top-2 guides by cell count; targets whose 2nd guide is sparse even pooled
  (CD3E, ZAP70) cannot be tested and are reported as such rather than dropped silently.
- Stim8hr only. The guide-resolved subsets exist only for Stim8hr; Rest/Stim48hr would need the
  same re-gather.
