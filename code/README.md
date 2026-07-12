# Code

Analysis and visualization scripts, grouped by project phase. Phases B, C, D, and F ran as interactive/agentic analyses over the precomputed DE-stats layer and external connectors (Open Targets, ChEMBL, AlphaFold, RCSB PDB, GWAS Catalog); their outputs are in `../results/` and their methods are documented in `../docs/`. The scripts here are the reusable single-cell and foundation pipelines.

## Layout

| Directory | Contents |
|---|---|
| `phaseA/` | Foundation pipeline: QC/schema audit, reproducibility landscape, positive-control benchmark over `GWCD4i.DE_stats.h5ad`. `phaseA.py` is the main run; `phaseA_regen.py` regenerates figures; `verify_oes.py` checks on-target effect-size fields. |
| `phaseE/` | Single-cell cross-donor + novel-axis validation. Gather scripts slab-stream the 111–172 GB cell-level shards into compact per-target subsets; analysis scripts compute knockdown efficiency, pseudobulk concordance, stimulation-dependence, and the state manifold. `.sbatch` files are the Slurm launchers. `scratch/` holds one-off diagnostic probes kept for provenance (not part of the pipeline). |
| `phaseG/` | Deep single-cell analyses on the 19 nominated targets: manifold construction, co-expression rewiring, activation-trajectory checkpoint, dual-guide concordance, guide-level subset gathering. |
| `viz/` | UMAP embedding + plotting utilities (per-target and per-module, by condition), manifest builder, and orchestration shell scripts. `*_config.json` are per-wave/per-donor plotting configs. |

## Environment

Scripts were run on a Slurm cluster in a dedicated conda environment:

```
python 3.11 · scanpy 1.11 · anndata 0.12 · mudata 0.3 · h5py · pandas · numpy · scipy · matplotlib · seaborn · leidenalg · igraph · aws-cli
```

## Adapting to your system

The scripts contain **absolute paths from the original cluster run** — a data root under `.../perturb-seq_data/` and a conda-env prefix. To reproduce elsewhere, edit those paths to point at your own dataset location (see `../data/README.md` for how to obtain the data) and environment. The cell-level scripts (`phaseE/`, `phaseG/`) are memory- and IO-heavy — each shard is 111–172 GB and is streamed in sequential slabs rather than loaded whole; run them as batch jobs, not interactively.
