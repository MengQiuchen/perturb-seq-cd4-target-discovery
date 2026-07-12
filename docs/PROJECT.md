# PROJECT.md — Novel Drug Targets in CD4⁺ T Cell Perturb-seq

## Objective
Discover and prioritize **novel drug targets** — regulators of CD4⁺ T cell programs whose
inhibition (or activation) would be therapeutically valuable in autoimmune, allergic, and
immuno-oncology settings — by mining a genome-scale CRISPRi Perturb-seq atlas from the
labs of **Alex Marson** (Gladstone/UCSF) and **Jonathan Pritchard** (Stanford).

## The dataset
Zhu, Dann, Yan, Reyes Retana, Goto, Guitche, Brixi, Ota, Pritchard & Marson,
*Genome-scale perturb-seq in primary human CD4+ T cells maps context-specific regulators
of T cell programs and human immune traits*, bioRxiv **2025.12.23.696273** (2025).
Distributed on the CZI Virtual Cells Platform.

- **~22 million** primary human CD4⁺ T cells.
- **Genome-scale CRISPRi** — knockdown of essentially all expressed genes, 2 guides/gene.
- **4 donors** (D1–D4) × **3 contexts**: Rest, Stim8hr, Stim48hr (TCR/CD28 timecourse).
- Precomputed analysis layers (pseudobulk, genome-wide DESeq2 DE stats with per-guide and
  per-donor breakdowns, off-target/reproducibility flags) + full cell-level matrices.

## Why this enables target discovery
A CRISPRi knockdown is a genetic model of drug-induced loss of function. Mapping each
perturbation to its genome-wide transcriptional effect lets us ask, in reverse: *which
upstream regulators control a therapeutically desirable program?* Those regulators are
drug-target hypotheses. The screen's cross-guide and cross-donor reproducibility metrics,
plus its context (rest vs stimulation) axis, make it possible to nominate targets that are
robust and act specifically on activated, disease-driving T cells.

## Approach (see `ANALYSIS_PLAN.md` for detail)
- **Phase A** — QC, schema audit, reproducibility landscape, positive-control benchmark.
- **Phase B** — six parallel discovery directions: cytokine master regulators; helper
  lineage polarization; context-specific (stim-dependent) regulators; network hubs;
  human-genetics (GWAS) integration; druggability annotation.
- **Phase C** — composite ranked `TARGET_SCORECARD.csv` + `TARGET_REPORT.md` dossiers.
- **Phase D** (optional) — single-cell validation on shortlisted targets.
- **Phase E** — single-cell cross-donor reproduction + validation of the novel chromatin/TF axis
  (extends Phase D from 1 donor / 7 flagship targets to 4 donors / 19 targets).
- **Phase F** — translational target dossiers (structure, chemical matter, tractability, clinical
  precedent, human genetics) for the validated targets.

## Key deliverables
- `high_confidence_targets.csv` — reproducible, on-target, off-target-clean perturbation set.
- Per-direction candidate tables (cytokine / polarization / context / hub / GWAS / drug).
- `TARGET_SCORECARD.csv` — integrated ranked target list (**primary output**).
- `TARGET_REPORT.md` — top-target dossiers with figures and therapeutic hypotheses.
- `PHASE_E_RESULTS.md` — single-cell cross-donor + novel-axis validation report (19 targets, 4 donors).
- `TARGET_DOSSIERS.md` — Phase F translational dossiers for the 19 validated targets (structure,
  chemical matter, tractability, clinical precedent, human genetics; readiness ranking).

## Status
- Cluster env + storage set up (`ssh:clust1-rocm-4`). See `CLAUDE.md`.
- **Data download complete** — all analysis layers (pseudobulk + 3 DE objects) + suppl
  tables + metadata; one cell-level shard (D4_Rest, 118.6 GB) fetched for validation.
- **All data now lives inside this project dir** under `perturb-seq_data/` (the project
  dir is itself a symlink onto Lustre scratch, so this is not on the home quota).
- **Environment validated** on rocm-4 (scanpy/anndata/mudata import + h5ad round-trip pass).
- **Phase A Discover done** — schema audit + exploratory stats over the 3 analysis layers;
  see `discover_outputs/` (DATASET_README.md, 2 figures, 2 summary tables).
- **Phase A Foundation DONE** — QC audit, reproducibility landscape, positive-control benchmark over DE_stats.h5ad. High-confidence set = 12,576 perturbations / 5,728 targets (phaseA_outputs/high_confidence_targets.csv); 13/16 positive controls pass. See phaseA_outputs/PHASE_A_RESULTS.md.
- **Phase B DONE** — 6 discovery directions over the HC set: D1 cytokine regulators (3,900 gene×condition×cytokine effects), D2 polarization (1,109 shifts), D3 context-specific (5,728 classified; 985 activation-family), D4 network hubs (12,576 hubs, 4 co-regulation modules recovering TCR-proximal + SAGA/Mediator complexes), D5 GWAS integration (127/301 candidates with immune GWAS support), D6 druggability (178 novel-druggable). Ran as 4+2 parallel sub-agents + shared DE extraction.
- **Phase C DONE** — integrated composite scorecard. **301 convergent candidates** (≥2 directions; 61 at ≥3; LAT/PLCG1/SENP5 at all 4). Top targets: LAT, SMARCE1, PLCG1, CD247, STAT6, CD3E, ZAP70. Positive controls (CD3E/ITK/MALT1/CD28/IL4R/STAT3) recovered de novo. See `phaseBC_outputs/TARGET_SCORECARD.csv` + `TARGET_REPORT.md`.
- **Phase D DONE** — single-cell validation of 7 flagship targets on donor D4 (Rest vs Stim8hr; fetched the 134.6 GB `D4_Stim8hr` shard). **All 7 validate**: on-target knockdown 76-94% (KD vs NTC, MWU p<1e-8 for 6/7; SMARCE1 p=6.2e-4); single-cell signatures concordant with pseudobulk DE_stats (Pearson r 0.43-0.85); TCR-proximal targets (LAT/PLCG1/CD247/CD3E/ZAP70/VAV1) show 4.2-8.9x stronger effect in Stim8hr than Rest and shift KD cells from activated/effector toward naive/memory state; SMARCE1 confirmed constitutive (0.87x, equal in both) as internal specificity control. See `phaseD_outputs/PHASE_D_RESULTS.md`.
- **Phase E DONE** — single-cell cross-donor + novel-axis validation, **19 targets across up to 4 donors** (Stim8hr; 5 compact subsets of 82,398–82,715 cells each gathered by slab-streaming the 111–172 GB shards). **E1**: all **7 flagship** TCR-proximal targets reproduce across donors (KD>15% AND downstream concordance r>0.30 in every powered donor; D4 reproduces Phase D to ≤0.001). **E2**: all **12 novel chromatin/TF targets validate** single-cell — SAGA/Mediator co-activators (MED24 r≈0.80, MED12 0.76, TADA2B 0.74) concordance exceeds several flagship targets; remodelers CHD4/SMARCB1 weaker KD (48–52%) but concordant. **E3**: the novel axis is **mechanistically distinct** from TCR-proximal signaling — Ward clustering recovers the axis split at 89% (17/19), state-redistribution 61% (TCR) vs 26% (chromatin/TF) p=1e-4, effect magnitude 2.4× p=3e-4. Underpowered cells flagged (ZAP70/D3 n=2, VAV1/D3 n=13, SMARCE1/D4 n=12). See `PHASE_E_RESULTS.md` + 3 figures + 11 tables.
- **Phase F DONE** — translational dossiers for all **19 validated targets** (`TARGET_DOSSIERS.md`). Structure: AlphaFold models + per-residue pLDDT domain maps for all 19; **17/19 have experimental PDB structures** (LAT, TADA2B AlphaFold-only). Pharmacology: ChEMBL bioactivity + Open Targets tractability/known-drugs. **18/19 targets are clinically unprecedented** (only CD3E carries clinical-stage drugs — the anti-CD3 antibody class); **all 12 novel-axis targets have zero known drugs**. Cross-track readiness ranking: **STAT6** #1 (55 asthma/allergy GWAS, 552 ChEMBL bioactivities, max pChEMBL 9.15, druggable-family SM bucket), then ZAP70, SIK3 (novel kinase, 809 bioactivities), and the drug-naive SAGA/Mediator readers MED24/SGF29/TADA2B. CHD4/SMARCB1 flagged for essentiality/toxicity risk. See `TARGET_DOSSIERS.md` + `phaseF_master_druggability.csv` + 4 figures.

## Team framing
Marson lab: human immune-cell engineering / functional genomics. Pritchard lab: statistical
& population genetics, gene regulation. The GWAS-integration direction (Phase B5) is where
these two strengths meet and is central to the paper's own thesis ("regulators of T cell
programs *and human immune traits*").

## Files in this project
- `PROJECT.md` — this overview.
- `ANALYSIS_PLAN.md` — detailed scientific plan (phases, directions, deliverables, methods).
- `CLAUDE.md` — operational context (cluster, env, paths, storage rules, commands).
- `data/` — README, 12 metadata `.jsonld`, 3 suppl tables (small reference files).
- `discover_outputs/` — Phase A results: `DATASET_README.md`, `discover_overview.png`,
  `discover_de_landscape.png`, `discover_summary_layers.csv`, `discover_summary_metrics.csv`.
- `phaseBC_outputs/` — Phase B/C results: `TARGET_SCORECARD.csv` (301 candidates × 35 cols),
  `TARGET_REPORT.md`, `top_targets_overview.png`, the 6 per-direction tables
  (cytokine/polarization/context/network/gwas/druggability), candidate-union tables,
  and 6 direction figures.
- `phaseD_outputs/` — Phase D results: `PHASE_D_RESULTS.md`, single-cell subset checkpoints,
  flagship pseudobulk signatures, validation figures/tables (7 flagship targets, donor D4).
- `phaseE_outputs/` — Phase E results: 5 single-cell subset checkpoints (D1–D4 Stim8hr + D4 Stim48hr),
  57 pseudobulk signatures, cross-donor + novel-axis + state-manifold analysis tables and scripts.
- `PHASE_E_RESULTS.md` + 3 figures (`E1_cross_donor.png`, `E2_novel_axis.png`, `E3_state_manifold.png`)
  + 11 tables — cross-donor reproducibility, novel-axis KD/concordance/stim-dependence, state-manifold contrast.
- `TARGET_DOSSIERS.md` + Phase F artifacts: `phaseF_master_druggability.csv`, structure files
  (AlphaFold + representative PDB, render in Mol*), and 4 figures (druggability landscape, readiness
  ranking, per-track pLDDT domain maps).
- `perturb-seq_data/` — the actual data (pseudobulk, DE `.h5ad`/`.h5mu`, `cell_level/`,
  `suppl_tables/`); ~225 GB, physically on Lustre via the project-dir symlink.
- `s3_manifest.csv` (artifact) — full S3 file manifest with sizes.
