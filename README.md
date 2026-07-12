# CD4⁺ T-Cell Perturb-seq Drug-Target Discovery

Discovery and translational prioritization of **novel druggable regulators of CD4⁺ T-cell programs** from a genome-scale CRISPRi Perturb-seq atlas. A CRISPRi knockdown is read as a genetic model of drug-induced loss of function; this pipeline mines the perturbation→transcriptome map **backwards** — from therapeutically desirable T-cell programs to the upstream regulators that control them — and carries the resulting candidates from genome-scale nomination through single-cell validation to translational dossiers.

> **Scope note.** This repository contains the **analysis code, curated result tables, figures, and reports** produced on top of a public dataset. It does **not** redistribute the raw/processed expression matrices (hundreds of GB to ~2.3 TB); those are fetched from the public source (see [Data access](#data-access)).

---

## Key results

- **Foundation (Phase A).** From the precomputed DESeq2 layer (33,983 perturbation×context contrasts × 10,282 genes) a reproducible, on-target, off-target-clean **high-confidence set of 12,576 perturbations / 5,728 genes** was defined. Known drug-target biology reproduces: **13/16 positive controls** move their canonical readouts in the expected direction.
- **Discovery + integration (Phases B–C).** Four orthogonal functional lenses + two external evidence layers (human genetics, druggability) converged on **301 candidate genes** (≥2 directions; 61 at ≥3; LAT/PLCG1/SENP5 at all 4). **127/301** carry immune-disease GWAS support; **178** are novel-druggable (tractable but drug-naive). 27 known drug targets are recovered *de novo*.
- **Single-cell validation (Phases D–E).** All **7 flagship TCR-proximal targets** reproduce across **4 donors**, and — the central novelty result — all **12 drug-naive chromatin/transcriptional-axis targets validate** at single-cell resolution. The novel axis is **mechanistically distinct** from acute TCR signaling.
- **Translational dossiers (Phase F).** Structure + chemical matter + tractability + clinical precedent + genetics for all **19 validated targets**. **18/19 are clinically unprecedented** (only CD3E is drugged). Top-ranked opportunity: **STAT6**; other standouts SIK3 (novel kinase) and the drug-naive SAGA/Mediator readers MED24 / SGF29 / TADA2B.
- **Deep single-cell analysis (Phase G).** Nine analyses (three tiers, 19 targets, 4 donors pooled) answering questions only single-cell resolution can answer — activation-trajectory gating, 35 cross-donor-consistent differential-abundance events, a Treg-expansion rare-state safety signal for MED12/MED24, and 17/17 dual-guide-concordant targets.

**Bottom line.** The pipeline recovers established immune drug targets *de novo* (validating the logic) and nominates a genetically-supported, single-cell-validated, structurally-tractable set of **novel** CD4⁺ T-cell regulators — headed by STAT6 and a drug-naive chromatin/transcriptional axis — as fresh therapeutic entry points for autoimmune/allergic and immuno-oncology indications.

See **[`docs/FULL_REPORT.md`](docs/FULL_REPORT.md)** for the complete consolidated report (Phases A–G, methods, figures, and an explicit boundaries/confidence section).

---

## Repository layout

```
perturb-seq-cd4-target-discovery/
├── README.md                 # this file
├── LICENSE                   # MIT
├── CITATION.cff              # how to cite this work + the source dataset
├── docs/                     # English reports (start here)
│   ├── FULL_REPORT.md        # consolidated Phases A–G report — read this first
│   ├── PROJECT.md            # objective + status
│   ├── ANALYSIS_PLAN.md      # detailed scientific plan
│   ├── IMMUNOLOGY_PRIMER.md  # background primer for readers new to T-cell immunology
│   ├── PHASE_A_RESULTS.md    # QC, reproducibility, positive controls
│   ├── PHASE_B_RESULTS.md    # six discovery directions
│   ├── PHASE_C_RESULTS.md / TARGET_REPORT.md  # integrated scorecard + dossiers
│   ├── PHASE_D_RESULTS.md    # single-cell validation (1 donor)
│   ├── PHASE_E_RESULTS.md / PHASE_E_METHODS.md  # cross-donor + novel-axis validation
│   ├── TARGET_DOSSIERS.md    # Phase F translational dossiers (19 targets)
│   └── PHASE_G_*_RESULTS.md  # deep single-cell analysis (three tiers)
├── code/                     # analysis + visualization scripts, grouped by phase
│   ├── phaseA/               # QC / reproducibility / positive-control pipeline
│   ├── phaseE/               # single-cell cross-donor + novel-axis (+ scratch/ diagnostics)
│   ├── phaseG/               # deep single-cell analyses
│   └── viz/                  # UMAP embedding + plotting
├── results/                  # curated small outputs (tables, figures, structures)
│   ├── discover/ phaseA/ phaseBC/ phaseD/ phaseE/ phaseF/ phaseG/
│   └── phaseF/structures/    # AlphaFold + representative experimental PDB/CIF
└── data/                     # reference tables + how to obtain the full data
    ├── suppl_tables/         # sgRNA library, sample metadata, DE-stats supplementary table
    ├── metadata/             # per-shard dataset metadata (.jsonld)
    └── s3_manifest.csv       # full file manifest (paths + sizes) for the public S3 bucket
```

Each subdirectory has its own `README.md` describing its contents.

---

## The dataset

Zhu, Dann, Yan, Reyes Retana, Goto, Guitche, Petersen, Ota, Shy, Pritchard & Marson.
*Genome-scale perturb-seq in primary human CD4⁺ T cells maps context-specific regulators of T cell programs and human immune traits.* bioRxiv **2025.12.23.696273** (2025). Distributed on the CZI Virtual Cells Platform.

- **~22 million** primary human CD4⁺ T cells.
- **Genome-scale CRISPRi** — knockdown of essentially all expressed genes, 2 guides/gene.
- **4 donors** (D1–D4) × **3 contexts**: Rest, Stim8hr, Stim48hr (TCR/CD28 stimulation timecourse).
- Precomputed analysis layers (pseudobulk; genome-wide DESeq2 DE stats with per-guide and per-donor breakdowns; off-target/reproducibility flags) plus full cell-level matrices.

The **workhorse signal** is the `zscore` layer (= log2FC / lfcSE) of the precomputed genome-wide DESeq2 object (`GWCD4i.DE_stats.h5ad`), which avoids re-running differential expression on 22 M cells. Cell-level shards were used only for the single-cell validation (Phases D, E, G).

---

## Data access

The expression matrices are **not** included in this repository (they range from ~17 GB per analysis layer to 111–172 GB per cell-level shard, ~2.3 TB total). Obtain them from the public source:

- **Dataset home:** https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq
- **Public S3 bucket (no-sign-request):** `s3://genome-scale-tcell-perturb-seq/marson2025_data/`

```bash
# example: fetch one analysis layer
aws s3 cp --no-sign-request \
  s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad ./
```

The full file list with sizes is in [`data/s3_manifest.csv`](data/s3_manifest.csv). See [`data/README.md`](data/README.md) for a description of each layer and the built-in confidence flags used throughout.

---

## Reproducing the analysis

The analysis was run on a Slurm cluster with a dedicated conda environment (Python 3.11: scanpy 1.11, anndata 0.12, mudata 0.3, h5py, pandas, numpy, scipy, matplotlib, seaborn, leidenalg, igraph). External evidence (structures, tractability, genetics) was pulled from AlphaFold, RCSB PDB, ChEMBL, Open Targets, and the GWAS Catalog.

- **Phase A (foundation)** and **B/C (discovery + integration)** run on the precomputed DE-stats / pseudobulk objects — CPU, no cell-level data required.
- **Phases D / E / G (single-cell validation)** stream selected cell-level shards; these are memory- and IO-heavy (shards are 111–172 GB) and were run as detached jobs on the cluster.

See [`code/README.md`](code/README.md) for the environment details and how each script maps to a phase. Scripts contain absolute cluster paths from the original run; adapt the data-root and environment paths to your system.

---

## Analysis boundaries

The consolidated report is explicit about what the data can and cannot support (see `docs/FULL_REPORT.md` §10). In brief: **CRISPRi models loss-of-function only** (activation hypotheses are inferences from the opposite direction, not measurements); pseudobulk DESeq2 is a population-level contrast; single-cell validation is Stim8hr-centric and read at the KD-population vs NTC-population level; and the structural "druggability" layer is AlphaFold-model-based (pLDDT = fold confidence, not verified pocket geometry — no docking/MD was performed).

---

## License

Code and documentation in this repository are released under the **MIT License** (see [`LICENSE`](LICENSE)). The underlying dataset is distributed by its authors under the MIT license via the CZI Virtual Cells Platform; please cite the dataset (see [`CITATION.cff`](CITATION.cff)) and respect its terms.
