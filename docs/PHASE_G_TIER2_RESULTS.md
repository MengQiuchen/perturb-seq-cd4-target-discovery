# Phase G — Second-tier single-cell analyses (state / manifold level)

**What this tier adds.** First-tier (dose–response, responder fraction, baseline sensitivity)
worked on each target's own signature. Second-tier works on the **shared cell-state manifold**:
where KD cells sit relative to NTC cells in transcriptional-state space, whether that
displacement reproduces across donors, and whether it touches rare, safety-relevant states.
This upgrades the report's E3 analysis (single donor D4, 13 fixed clusters) to a
**4-donor harmony-integrated reference** and adds cross-donor statistics E3 could not do.

**Substrate — cross-donor reference manifold** (`G4_manifold_umap_Stim8hr.png`,
`G4_manifold_Stim8hr.h5ad` on cluster).
- 58,442 cells: all KD cells for the 19 targets + 12,000 NTC/donor, 4 donors pooled.
- E3 recipe (2000 HVG → PCA30 → UMAP → Leiden) **+ harmony(donor)** integration. Donors
  co-embed cleanly (Fig G4a), so residual structure is biological, not batch.
- 15 Leiden states, each labeled by its dominant program (naive_memory, activation, effector,
  cytokine, Treg-like, Tfh-like, exhaustion). Cluster 6/12 = naive/low-activation; 0/4/8 =
  Treg-like; 13 = Tfh-like (BCL6/CXCR5); 2/7 = effector/cytokine.

---

## Result 4 — Cross-donor neighborhood differential abundance (`G5_*`, `G5_crossdonor_DA_*.csv`)

For every target × state, the KD−NTC occupancy shift is tested across the 4 donors
(one-sample t on per-donor shifts) and required to be **direction-consistent in every donor**.

- **35 significant, donor-consistent DA events across 14 targets** (FDR<0.05).
- **TCR-proximal KD → cells pile into the naive / low-activation state and drain from
  activation states**: ZAP70 → cluster 6 (naive) shift **+0.78** and CD3E → cluster 6 **+0.73**
  are both FDR-significant and donor-consistent; CD3E additionally shows significant depletion of
  activation clusters 5/10/11 (ZAP70's shifts there are same-signed but do not clear FDR<0.05,
  padj≈0.06–0.15). This is the "failure to activate" phenotype, now shown reproducible across
  donors rather than in one.
- Chromatin/TF targets show weaker, program-specific shifts (below).

## Result 5 — Rare-state (safety) effects (`G6_*` panel a, `G6_rarestate_*.csv`)  ← key safety finding

The report flagged Treg/Tfh programs as **underpowered in bulk**; the pooled manifold isolates
them. Testing occupancy of the Treg-like (clusters 0/4/8) and Tfh-like (cluster 13) compartments:

- **Mediator-kinase module KD *expands* the Treg-like compartment**: MED12 **+0.41**
  (padj=0.013), MED24 **+0.20** (padj=0.036), TADA2B trend. If a drug against these targets
  likewise skewed cells toward a regulatory phenotype, that is an efficacy/safety consideration
  (immunosuppression) that only shows up at the rare-state level.
- **TCR-proximal + SMARCE1/SMARCB1/ARNT KD *deplete* Treg-like cells** (ZAP70 −0.22, CD3E −0.21,
  CD247 −0.15, SMARCE1 −0.13, SMARCB1 −0.13, ARNT −0.13) — consistent with impaired activation
  reducing all activation-dependent fates including induced-Treg.
- Tfh-like compartment is small (NTC occupancy ~0.9%); only CD3E reaches significance (depletion).

## Result 6 — Cross-donor concordance of state redistribution (`G6_*` panel b, `G7_*.csv`)

For each target, the full per-cluster shift **vector** is correlated between every donor pair.

- **TCR-proximal targets redistribute cells far more reproducibly across donors**
  (mean pairwise r = **0.84**) than chromatin/TF targets (**0.54**); SMARCE1 lowest (0.38).
- Most reproducible: CD3E r=0.99, MED12 0.98, ZAP70 0.98, TADA2B 0.95.
- Least reproducible (donor-variable redistribution): NSD1 0.04, SIK3 0.18, STAT6 0.26,
  SMARCB1 0.29 — these move cells, but *where* varies by donor.
- This independently confirms, at single-cell state resolution, the report's flagship claim that
  TCR-proximal biology is the more robust/reproducible axis.

---

## Files

- `G4_manifold_umap_Stim8hr.png` — reference manifold (donor / program / activation gradient)
- `G5_crossdonor_DA_Stim8hr.png` + `.csv` — target × program DA heatmap; CSV has per-cluster
  per-donor t, p, padj, direction-consistency
- `G6_rarestate_concordance_Stim8hr.png` — safety (Treg/Tfh) + concordance
- `G6_rarestate_Stim8hr.csv`, `G7_crossdonor_concordance_Stim8hr.csv`
- `G4_cluster_profiles_Stim8hr.csv`, `G4_occupancy_by_donor_Stim8hr.csv` — underlying tables
- `G4_manifold_Stim8hr.h5ad` (cluster) — the integrated manifold checkpoint
- Script: `phaseE_scripts/phaseG_manifold.py`

## Caveats

- Harmony integration removes donor as a linear batch covariate; strong donor-specific biology
  could still be partially absorbed. The concordance analysis (r per donor pair) is the honest
  check and shows the axis-level differences survive integration.
- Cluster labels are argmax over 7 program scores, not curated cell-type calls; "Treg-like"
  means highest Treg-program score, not a FOXP3⁺ gate. The direction of the MED12/MED24 Treg
  effect is robust; the absolute compartment definition is approximate.
- Underpowered target×donor cells (<10 KD cells in a donor) are dropped from that donor's DA row,
  so some targets are tested on 3 donors (noted in `n_donors`).
- Stim8hr only; Rest/Stim48hr manifolds can be built with the same script if a timecourse view
  of state redistribution is wanted.
