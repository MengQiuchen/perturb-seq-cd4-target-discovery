"""Wave 1: joint D4 Rest + Stim8hr embedding."""
import os, sys, time, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from embed_common import log, load_shards, embed

PROJ = "/Users/meng01/qiuchen/project/hackathon/perturb-seq"
D = f"{PROJ}/perturb-seq_data/cell_level"
OUT = f"{PROJ}/viz_outputs"
SCORECARD = f"{PROJ}/phaseBC_outputs/TARGET_SCORECARD.csv"

PARAMS = dict(nt_cap=25_000, per_target_cap=400, bg_cap=90_000, chunk=200_000, seed=0)

def main():
    t0 = time.time()
    sc_df = pd.read_csv(SCORECARD)
    nominated = set(sc_df["gene"].dropna().astype(str))
    log(f"nominated targets from scorecard: {len(nominated)}")

    shards = [
        dict(donor="D4", condition="Rest",    path=f"{D}/D4_Rest.assigned_guide.h5ad"),
        dict(donor="D4", condition="Stim8hr",  path=f"{D}/D4_Stim8hr.assigned_guide.h5ad"),
    ]
    for sh in shards:
        assert os.path.exists(sh["path"]), f"missing {sh['path']}"

    adata = load_shards(shards, nominated, PARAMS)
    log(f"combined AnnData: {adata.shape}")
    log("obs condition counts:\n" + str(adata.obs["condition"].value_counts()))
    log("obs guide_type counts:\n" + str(adata.obs["guide_type"].value_counts()))
    log(f"nominated cells: {int(adata.obs['is_nominated'].sum()):,}")

    adata = embed(adata, integrate_donor=False, seed=PARAMS["seed"])

    ck = f"{OUT}/checkpoints/D4_Rest_Stim8hr.embedded.h5ad"
    log(f"writing checkpoint {ck}")
    adata.write_h5ad(ck)
    # also dump a small obs+umap table for quick plotting/inspection
    df = adata.obs.copy()
    df["UMAP1"] = adata.obsm["X_umap"][:, 0]
    df["UMAP2"] = adata.obsm["X_umap"][:, 1]
    df.to_parquet(f"{OUT}/checkpoints/D4_Rest_Stim8hr.obs_umap.parquet")
    log(f"DONE wave1 in {time.time()-t0:.0f}s; checkpoint cells={adata.n_obs:,}")

if __name__ == "__main__":
    main()
