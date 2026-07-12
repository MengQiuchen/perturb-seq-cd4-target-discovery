"""Generic config-driven embedding runner (waves 2/3).
Usage: python run_embed.py <config.json>
config = {tag, shards:[{donor,condition,path}], integrate_donor,
          nt_cap, per_target_cap, bg_cap, chunk, seed, n_hvg, n_pcs}
"""
import os, sys, time, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from embed_common import log, load_shards, embed

PROJ = "/Users/meng01/qiuchen/project/hackathon/perturb-seq"
OUT = f"{PROJ}/viz_outputs"
SCORECARD = f"{PROJ}/phaseBC_outputs/TARGET_SCORECARD.csv"

def main():
    cfg = json.load(open(sys.argv[1]))
    tag = cfg["tag"]
    t0 = time.time()
    nominated = set(pd.read_csv(SCORECARD)["gene"].dropna().astype(str))
    log(f"[{tag}] nominated targets: {len(nominated)}")
    for sh in cfg["shards"]:
        assert os.path.exists(sh["path"]), f"missing {sh['path']}"
    params = dict(nt_cap=cfg.get("nt_cap",25_000), per_target_cap=cfg.get("per_target_cap",400),
                  bg_cap=cfg.get("bg_cap",90_000), chunk=cfg.get("chunk",200_000), seed=cfg.get("seed",0))
    log(f"[{tag}] params={params} integrate_donor={cfg.get('integrate_donor',False)}")
    adata = load_shards(cfg["shards"], nominated, params)
    log(f"[{tag}] combined {adata.shape}; conditions={dict(adata.obs['condition'].value_counts())}; "
        f"donors={dict(adata.obs['donor'].value_counts())}; nominated={int(adata.obs['is_nominated'].sum()):,}")
    adata = embed(adata, n_hvg=cfg.get("n_hvg",2000), n_pcs=cfg.get("n_pcs",50),
                  integrate_donor=cfg.get("integrate_donor",False), seed=params["seed"])
    ck = f"{OUT}/checkpoints/{tag}.embedded.h5ad"
    log(f"[{tag}] writing {ck}")
    adata.write_h5ad(ck)
    df = adata.obs.copy()
    df["UMAP1"], df["UMAP2"] = adata.obsm["X_umap"][:,0], adata.obsm["X_umap"][:,1]
    try:
        df.to_parquet(f"{OUT}/checkpoints/{tag}.obs_umap.parquet")
    except Exception as e:
        log(f"[{tag}] parquet unavailable ({e}); writing CSV instead")
        df.to_csv(f"{OUT}/checkpoints/{tag}.obs_umap.csv.gz", compression="gzip")
    log(f"[{tag}] DONE embed in {time.time()-t0:.0f}s; cells={adata.n_obs:,}")

if __name__ == "__main__":
    main()
