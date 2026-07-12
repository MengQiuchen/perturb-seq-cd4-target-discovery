"""Generic config-driven plot driver (waves 2/3).
Usage: python run_plot.py <tag>
Loads viz_outputs/checkpoints/<tag>.embedded.h5ad and emits figures A/B/C/D/E.
"""
import os, sys, time, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import anndata as ad
import scipy.sparse as sp
import plot_umaps as P

PROJ = "/Users/meng01/qiuchen/project/hackathon/perturb-seq"
OUT = f"{PROJ}/viz_outputs"; FIG = f"{OUT}/figures"
SCORECARD = f"{PROJ}/phaseBC_outputs/TARGET_SCORECARD.csv"
ACT_MARKERS = ["IL2RA","CD69","MKI67","TNFRSF9","IFNG","SELL"]

def resolve_gene_rows(adata, names):
    out=[]; gn = adata.var["gene_name"].values if "gene_name" in adata.var else None
    vn = adata.var_names.values
    for nm in names:
        gi=None
        if gn is not None:
            h=np.where(gn==nm)[0]
            if len(h): gi=int(h[0])
        if gi is None:
            h=np.where(vn==nm)[0]
            if len(h): gi=int(h[0])
        out.append((nm,gi))
    return out

def main():
    tag = sys.argv[1]
    t0=time.time()
    ck = f"{OUT}/checkpoints/{tag}.embedded.h5ad"
    P.log(f"[{tag}] loading {ck}")
    adata = ad.read_h5ad(ck)
    xy = adata.obsm["X_umap"]; obs = adata.obs
    os.makedirs(FIG, exist_ok=True)
    sc_df = pd.read_csv(SCORECARD)
    present = set(obs["perturbed_gene_name"].unique())
    top_targets = [g for g in sc_df.sort_values("rank")["gene"].astype(str) if g in present][:25]
    module_map={}
    if "module_name" in sc_df.columns:
        for mod,sub in sc_df.dropna(subset=["module_name"]).groupby("module_name"):
            genes=[g for g in sub["gene"].astype(str) if g in present]
            if genes: module_map[str(mod)]=genes
    dir_map={}
    if "therapeutic_direction" in sc_df.columns:
        for d,sub in sc_df.dropna(subset=["therapeutic_direction"]).groupby("therapeutic_direction"):
            genes=[g for g in sub["gene"].astype(str) if g in present]
            if genes: dir_map[f"dir: {d}"]=genes
    lognorm = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
    def getcol(gi):
        col=lognorm[:,gi]; return np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
    markers = resolve_gene_rows(adata, ACT_MARKERS)

    P.fig_condition_comparison(xy, obs, FIG, tag, marker_genes=markers, lognorm_getter=getcol)
    P.fig_within_condition(xy, obs, FIG, tag)
    P.fig_cross_donor(xy, obs, FIG, tag)               # no-op if single donor
    P.fig_per_target(xy, obs, FIG, tag, top_targets, ncol=5)
    P.fig_per_module(xy, obs, FIG, tag, module_map, ncol=4)
    if dir_map:
        P.fig_per_module(xy, obs, FIG, f"{tag}_bydirection", dir_map, ncol=2)

    figs=sorted([f for f in os.listdir(FIG) if f.endswith(".png") and tag in f])
    json.dump({"tag":tag,"n_cells":int(adata.n_obs),
               "donors":sorted(set(obs['donor'])),"conditions":sorted(set(obs['condition'])),
               "figures":figs,"top_targets":top_targets,
               "modules":{k:len(v) for k,v in module_map.items()}},
              open(f"{OUT}/checkpoints/{tag}_manifest.json","w"), indent=2)
    P.log(f"[{tag}] DONE plotting in {time.time()-t0:.0f}s; {len(figs)} figures")

if __name__ == "__main__":
    main()
