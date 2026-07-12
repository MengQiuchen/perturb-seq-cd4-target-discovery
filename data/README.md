# Data reference

This directory holds **small reference tables and metadata only**. The expression matrices themselves (analysis layers 17–45 GB each; cell-level shards 111–172 GB each, ~2.3 TB total) are **not** included — obtain them from the public source described below.

## Contents

| Path | Description |
|---|---|
| `suppl_tables/sgrna_library_metadata.suppl_table.csv` | sgRNA library: guide sequences and target-gene mapping (~10 MB). |
| `suppl_tables/DE_stats.suppl_table.csv` | Tabular form of the DE-stats `.obs` (per target×condition summary; ~5 MB). |
| `suppl_tables/sample_metadata.suppl_table.csv` | Per-sample (donor × condition × lane) metadata. |
| `metadata/marson_D*_*.assigned_guide.jsonld` | Croissant JSON-LD metadata for each of the 12 cell-level shards. |
| `s3_manifest.csv` | Full file manifest for the public S3 bucket — every object's key and size. |

## The dataset

Zhu, Dann, Yan, Reyes Retana, Goto, Guitche, Petersen, Ota, Shy, Pritchard & Marson.
*Genome-scale perturb-seq in primary human CD4⁺ T cells maps context-specific regulators of T cell programs and human immune traits.* bioRxiv **2025.12.23.696273** (2025). CZI Virtual Cells Platform, version 1.0.0 (released 2025-12-22). License: **MIT**. Citation status: unpublished at time of analysis.

- **~22 million** primary human CD4⁺ T cells, **4 donors (D1–D4) × 3 contexts** (Rest, Stim8hr, Stim48hr).
- Genome-scale CRISPRi (2 guides/gene); **12,731 perturbed genes**, 25,954 guides + non-targeting controls.
- Platform: Ultima sequencing; GEMX_flex_v2 (10x Genomics Flex).

## Getting the data

```bash
# public, no-sign-request S3 bucket
BUCKET=s3://genome-scale-tcell-perturb-seq/marson2025_data

# analysis layers (used by Phases A–C)
aws s3 cp --no-sign-request $BUCKET/GWCD4i.DE_stats.h5ad ./          # 16.8 GB — primary object
aws s3 cp --no-sign-request $BUCKET/GWCD4i.pseudobulk_merged.h5ad ./ # 44.6 GB
aws s3 cp --no-sign-request $BUCKET/GWCD4i.DE_stats.by_guide.h5mu ./ # 29.4 GB
aws s3 cp --no-sign-request $BUCKET/GWCD4i.DE_stats.by_donors.h5mu ./# 16.9 GB

# a cell-level shard (used by Phases D/E/G) — one donor × condition, 111–172 GB
aws s3 cp --no-sign-request $BUCKET/D4_Rest.assigned_guide.h5ad ./
```

See `s3_manifest.csv` for the complete file list with sizes.

## Data layers

| Layer / object | Size | Dimensions | Role |
|---|---|---|---|
| `GWCD4i.DE_stats.h5ad` | 16.8 GB | 33,983 (target×condition) × 10,282 genes | **Primary object** — genome-wide DESeq2. |
| `GWCD4i.pseudobulk_merged.h5ad` | 44.6 GB | 278,684 × 18,129 | Per (guide×donor×condition) pseudobulk. |
| `GWCD4i.DE_stats.by_guide.h5mu` | 29.4 GB | MuData | Per-guide DE (two-guide agreement). |
| `GWCD4i.DE_stats.by_donors.h5mu` | 16.9 GB | MuData | Per-donor-pair DE (cross-donor agreement). |
| `D*_*.assigned_guide.h5ad` | 111–172 GB × 12 | ~2.7 M cells × 18,130 genes each | Cell-level UMI counts (sparse CSR). |

### Key fields

**DE-stats `.layers`** (per target × gene): `log_fc`, `zscore` (= log2FC / lfcSE, the workhorse signal), `p_value`, `adj_p_value`, `baseMean`, `lfcSE`.

**DE-stats `.obs`**: `n_up_genes` / `n_down_genes` / `n_total_de_genes` (10% FDR), `ontarget_effect_size`, `ontarget_significant`, and the confidence/specificity flags used throughout the analysis —

- `guide_correlation_signif` — the two guides agree (reproducibility).
- `donor_correlation_*` — donors agree (reproducibility).
- `neighboring_gene_KD` — cis neighboring-gene knockdown (specificity flag).
- `distal_offtarget_flag` — distal off-target (exclusion flag).
- `low_target_gex` — target gene lowly expressed (exclusion flag).

The **high-confidence nomination set** (`../results/phaseA/high_confidence_targets.csv`) is defined as: on-target significant AND (two-guide OR cross-donor agreement) AND not flagged for off-target/low-expression. Use these built-in flags rather than re-deriving them.

**Significance threshold used throughout: 10% FDR.**
