# Final Target List — CD4⁺ T-Cell Perturb-seq Drug-Target Discovery

The **19 regulators of CD4⁺ T-cell programs** nominated by this project and validated at single-cell resolution across up to 4 donors (Phases A–G), ranked by a transparent translational-readiness heuristic. This is the project's headline deliverable: a prioritized, evidence-backed set of therapeutic entry points.

Machine-readable version: [`FINAL_TARGET_LIST.csv`](FINAL_TARGET_LIST.csv). Full per-target narrative dossiers: [`docs/TARGET_DOSSIERS.md`](docs/TARGET_DOSSIERS.md) and the Phase F report [`docs/PHASE_F_RESULTS.md`](docs/PHASE_F_RESULTS.md).

## How this list was built

Each CRISPRi knockdown is read as a genetic model of drug-induced loss of function. Genome-scale nomination (Phase A–C) converged on 301 candidate regulators; the 19 below are those carried through **single-cell validation** (Phase D/E: on-target knockdown + pseudobulk concordance in every powered donor) and **translational annotation** (Phase F: structure, chemical matter, tractability, clinical precedent, human genetics). The **translational-readiness** score is a transparent heuristic (structure + chemical matter + genetics − essentiality risk); it intentionally rewards *existing* chemical matter and genetics, so it under-ranks the most novel drug-naive targets — read it together with the axis and novelty columns, not as a hard order.

**Headline.** 18/19 targets are clinically unprecedented (only **CD3E** carries approved/clinical-stage drugs against the target). All 12 chromatin/transcription-factor "novel-axis" targets have zero known drugs. Top-ranked opportunity: **STAT6**; standouts also include the novel kinase **SIK3** and the drug-naive SAGA/Mediator readers **MED24 / SGF29 / TADA2B**.

## Ranked targets

| Rank | Gene | Axis | Direction | Novelty | KD | Conc. r | Donors | pLDDT | #PDB | SM tractability | ChEMBL (max pChEMBL) | Known drugs | Immune GWAS | Readiness |
|---:|---|---|---|---|---:|---:|:--:|---:|---:|---|---|---:|---:|---:|
| 1 | **STAT6** | chromatin/TF | boost-immunity | novel-druggable | 82% | 0.566 | 4/4 | 76.51 | 7 | High-Quality Ligand | 552 (9.15) | 0 | 55 | **6.5** |
| 2 | **ZAP70** | signaling | suppress-inflammation | novel-druggable | 91% | 0.782 | 3/3 | 84.92 | 15 | High-Quality Pocket | 2390 (8.1) | 0 | 1 | **6.1** |
| 3 | **SIK3** | signaling | suppress-inflammation | novel-druggable | 87% | 0.523 | 3/3 | 50.65 | 5 | High-Quality Ligand | 809 (9.34) | 0 | 0 | **4.5** |
| 4 | **MED24** | chromatin/TF | mixed | novel-druggable | 96% | 0.796 | 4/4 | 84.14 | 10 | Structure with Ligand | 0 | 0 | 12 | **4.2** |
| 5 | **PLCG1** | signaling | suppress-inflammation | novel-druggable | 88% | 0.459 | 4/4 | 82.75 | 6 | High-Quality Ligand | 438 (6.75) | 0 | 0 | **4.0** |
| 6 | **VAV1** | signaling | suppress-inflammation | novel-druggable | 92% | 0.659 | 3/3 | 86.4 | 10 | Structure with Ligand | 1 (8.98) | 0 | 0 | **4.0** |
| 7 | **NSD1** | chromatin/TF | suppress-inflammation | novel-druggable | 92% | 0.653 | 3/3 | 44.68 | 4 | High-Quality Ligand | 147 (6.96) | 0 | 1 | **3.1** |
| 8 | **SEL1L** | chromatin/TF | suppress-inflammation | novel-druggable | 94% | 0.42 | 3/3 | 81.01 | 5 | — | 0 | 0 | 0 | **3.0** |
| 9 | **CD247** | signaling | suppress-inflammation | novel-druggable | 72% | 0.506 | 4/4 | 62.4 | 38 | Structure with Ligand | 0 | 0 | 46 | **2.5** |
| 10 | **SMARCE1** | chromatin/TF | suppress-inflammation | novel-druggable | 88% | 0.743 | 3/3 | 69.55 | 9 | Structure with Ligand | 8 (6.95) | 0 | 37 | **2.5** |
| 11 | **SGF29** | chromatin/TF | suppress-inflammation | novel-druggable | 94% | 0.692 | 4/4 | 91.75 | 8 | High-Quality Pocket | 1 | 0 | 5 | **2.5** |
| 12 | **CD3E** | signaling | suppress-inflammation | known-drug-target | 96% | 0.847 | 4/4 | 73.07 | 44 | Med-Quality Pocket | 0 | 22 | 3 | **1.8** |
| 13 | **ARNT** | chromatin/TF | mixed | novel-druggable | 89% | 0.691 | 4/4 | 55.5 | 46 | High-Quality Ligand | 25 | 0 | 1 | **1.6** |
| 14 | **TADA2B** | chromatin/TF | suppress-inflammation | difficult | 92% | 0.739 | 4/4 | 86.67 | 0 | — | 0 | 0 | 0 | **1.5** |
| 15 | **SMARCB1** | chromatin/TF | suppress-inflammation | novel-druggable | 48% | 0.587 | 4/4 | 80.77 | 18 | Structure with Ligand | 8 (8.03) | 0 | 0 | **1.5** |
| 16 | **TRIP12** | chromatin/TF | suppress-inflammation | novel-druggable | 84% | 0.696 | 4/4 | 66.77 | 5 | Structure with Ligand | 0 | 0 | 0 | **1.5** |
| 17 | **CHD4** | chromatin/TF | boost-immunity | novel-druggable | 52% | 0.381 | 4/4 | 64.62 | 12 | Structure with Ligand | 270 (8.52) | 0 | 0 | **1.0** |
| 18 | **LAT** | signaling | suppress-inflammation | novel-druggable | 85% | 0.671 | 4/4 | 59.35 | 0 | — | 2 | 0 | 0 | **0.0** |
| 19 | **MED12** | chromatin/TF | suppress-inflammation | difficult | 85% | 0.76 | 4/4 | 65.12 | 3 | — | 6 (6.89) | 0 | 0 | **0.0** |

## Column definitions

- **Axis** — mechanistic class: `signaling` (TCR-proximal signaling) or `chromatin/TF` (the novel chromatin/transcriptional-regulatory axis).
- **Direction** — therapeutic hypothesis inferred from the CRISPRi loss-of-function direction: `suppress-inflammation`, `boost-immunity`, or `mixed`. CRISPRi models loss-of-function only; directionality is an inference, not a gain-of-function measurement.
- **Novelty** — `novel-druggable` (tractable, no approved drug against the target), `known-drug-target`, or `difficult` (tractability-limited).
- **KD** — Phase E mean on-target knockdown (fraction), averaged across powered donors.
- **Conc. r** — Phase E mean single-cell-vs-pseudobulk signature concordance (Pearson r) across powered donors.
- **Donors** — donors with sufficient cell power for a verdict (e.g. 4/4, 3/3); under-powered target×donor cells were excluded, not pooled.
- **pLDDT** — AlphaFold mean per-residue confidence (model confidence, not verified pocket geometry).
- **#PDB** — number of experimental structures in the RCSB PDB.
- **SM tractability** — best small-molecule tractability bucket (Open Targets structural annotation).
- **ChEMBL (max pChEMBL)** — number of ChEMBL bioactivities and the representative max potency (ceiling over a 200-activity sample).
- **Known drugs** — count of clinical/approved drugs against the target (Open Targets).
- **Immune GWAS** — number of autoimmune/allergic/immune GWAS associations linked to the target (Phase B/C disease-genetics layer).
- **Readiness** — composite translational-readiness score (see above).

## Caveats

The readiness score is a judgment-weighted heuristic, not a fit; structure is monomeric AlphaFold (fold confidence, not a validated druggable pocket — no docking/MD was run); ChEMBL potency is a representative ceiling; and the therapeutic direction is inferred from the knockdown, not tested with a gain-of-function arm. See `docs/FULL_REPORT.md` §10 for the full boundaries/confidence discussion and `docs/PHASE_F_RESULTS.md` for per-axis detail.

*Source table: `results/phaseF/phaseF_master_druggability.csv` (19 × 25 columns), rank from `results/phaseF/phaseF_readiness_ranked.json`.*
