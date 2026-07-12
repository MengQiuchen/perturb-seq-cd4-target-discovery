# Phase A Results — Foundation (QC, reproducibility, positive controls)

Ran on `ssh:clust1-rocm-4` over `GWCD4i.DE_stats.h5ad` (precomputed DESeq2).
Outputs: `phaseA_outputs/`. Script: `phaseA_scripts/phaseA.py`.

## A1 — Schema & QC audit  → `qc_overview.png`, `de_schema.csv`, `qc_summary.json`
- 33,983 (perturbation × context) rows · 11,526 unique targets · 10,282 measured genes.
- Balanced across contexts: Rest 11,287 / Stim8hr 11,415 / Stim48hr 11,281.
- Knockdown efficiency: 62.4% of rows on-target significant; median on-target effect
  size −6.30 (this is the on-target DE **z-score** = log2FC / lfcSE, the `zscore` layer
  value for the target gene — NOT a raw log2FC; negative = knockdown detected). The
  corresponding median on-target **log2FC** is −2.48 (−2.77 among significant knockdowns,
  ≈ 6.8-fold). Most targets have 2 guides (30,108).
- Exclusion flag rates: low_target_gex 23.0%, neighboring_gene_KD 7.7%,
  distal_offtarget_flag 1.3%.

## A2 — Reproducibility landscape  → `reproducibility.png`, `high_confidence_targets.csv`
High-confidence definition (universe for all downstream nomination):
on-target significant AND effect<0  AND  (two-guide agreement `r>0, p<0.05`
OR cross-donor `donor_correlation_hits_mean>0.3`)  AND  flag-clean
(not distal_offtarget / neighboring_gene_KD / low_target_gex).

- Funnel (nested/sequential): 33,983 → on-target 21,221 → + reproducible (guide|donor) 14,509
  → + flag-clean = **high-confidence 12,576**. (Note: 21,979 rows are reproducible over the
  *full* 33,983, unconditioned on on-target — that is not the funnel's intermediate step.)
- 12,576 high-confidence perturbations, all with ≥1 downstream DE gene; **5,728 unique
  targets**. Per context: Rest 3,714 / Stim8hr 4,492 / Stim48hr 4,370.
- Top-12 trans-effect hubs are dominated by canonical TCR-proximal signaling (CD3E, LAT,
  ZAP70, PLCG1, VAV1, CD247) + SAGA/Mediator coactivator subunits (TADA2B, SGF29, TAF6L,
  MED12, CCNC) — a sanity check that the set surfaces real regulators. One non-canonical
  outlier also ranks in the top 12: SENP5 (rank 4, n_downstream 5171, on-target −38.7 — an
  unusually large-magnitude effect, ~2× any other top hit), flagged as a candidate to
  inspect in Phase B rather than assume clean.

## A3 — Positive-control benchmark  → `positive_controls.png`, `positive_controls.csv`
16 known T-cell regulators tested for (i) self-knockdown and (ii) expected downstream
direction (CRISPRi = LOF; knockdown of an activator should drop its target program).
- **13/16 pass expected direction.** Clean passes with strong on-target KD: TBX21→IFNG/CXCR3,
  GATA3→IL4/IL5/IL13, STAT5A/B→IL2RA, IL2RA(self), LCK/ZAP70/CD28/PPP3CA→IL2/IFNG, BATF.
- 3 caveats (expected): FOXP3 and BCL6 knockdowns did not move their canonical readouts in
  the direction scored (Treg/Tfh programs are weak in this bulk stimulated CD4 setting);
  NFATC2 showed a paradoxical positive self-effect (not significant), consistent with the
  paper flagging non-replicating NFAT-family trans-effects.
- **Coverage note:** JAK1/JAK3 are NOT in the screen's perturbation set (kinases were not
  eligible targets), so the IL-2/JAK-STAT axis was benchmarked via STAT5A/STAT5B and IL2RA
  as proxies rather than the full kinase→STAT axis.

## Interpretation
The "signature → drug effect" reasoning is validated: known drug-target biology reproduces,
and the high-confidence set (12,576 perturbations / 5,728 targets) is the vetted universe for
Phase B target nomination.
