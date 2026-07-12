"""Wave 1 plotting: load D4 Rest+Stim8hr embedded checkpoint, emit figures A/B/D/E."""
import os, sys, time, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import anndata as ad
import plot_umaps as P

PROJ = "/Users/meng01/qiuchen/project/hackathon/perturb-seq"
OUT = f"{PROJ}/viz_outputs"
FIG = f"{OUT}/figures"
CK = f"{OUT}/checkpoints/D4_Rest_Stim8hr.embedded.h5ad"
SCORECARD = f"{PROJ}/phaseBC_outputs/TARGET_SCORECARD.csv"
TAG = "D4_Rest_Stim8hr"

# activation markers to color (gene_name -> resolve to var row)
ACT_MARKERS = ["IL2RA","CD69","MKI67","TNFRSF9","IFNG","SELL"]

def resolve_gene_rows(adata, names):
    """Map gene_name -> integer var row index (using var['gene_name'] if present, else var_names)."""
    out = []
    gn = adata.var["gene_name"].values if "gene_name" in adata.var else None
    vn = adata.var_names.values
    for nm in names:
        gi = None
        if gn is not None:
            hit = np.where(gn == nm)[0]
            if len(hit): gi = int(hit[0])
        if gi is None:
            hit = np.where(vn == nm)[0]
            if len(hit): gi = int(hit[0])
        out.append((nm, gi))
    return out

def main():
    t0=time.time()
    P.log(f"loading checkpoint {CK}")
    adata = ad.read_h5ad(CK)
    P.log(f"loaded {adata.shape}")
    xy = adata.obsm["X_umap"]
    obs = adata.obs
    os.makedirs(FIG, exist_ok=True)

    sc_df = pd.read_csv(SCORECARD)
    # top targets by rank (composite) present in the data
    present = set(obs["perturbed_gene_name"].unique())
    top_targets = [g for g in sc_df.sort_values("rank")["gene"].astype(str) if g in present][:25]
    P.log(f"top targets present: {top_targets}")

    # module groups: group nominated targets by module_name (co-regulation modules)
    module_map = {}
    if "module_name" in sc_df.columns:
        for mod, sub in sc_df.dropna(subset=["module_name"]).groupby("module_name"):
            genes = [g for g in sub["gene"].astype(str) if g in present]
            if genes: module_map[str(mod)] = genes
    # therapeutic-direction groups too
    dir_map = {}
    if "therapeutic_direction" in sc_df.columns:
        for d, sub in sc_df.dropna(subset=["therapeutic_direction"]).groupby("therapeutic_direction"):
            genes=[g for g in sub["gene"].astype(str) if g in present]
            if genes: dir_map[f"dir: {d}"]=genes
    P.log(f"modules: {list(module_map)}; directions: {list(dir_map)}")

    # lognorm getter for marker coloring
    lognorm = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
    import scipy.sparse as sp
    def getcol(gi):
        col = lognorm[:, gi]
        return np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()

    markers = resolve_gene_rows(adata, ACT_MARKERS)
    P.log(f"marker rows: {markers}")

    P.fig_condition_comparison(xy, obs, FIG, TAG, marker_genes=markers, lognorm_getter=getcol)
    P.fig_within_condition(xy, obs, FIG, TAG)
    P.fig_per_target(xy, obs, FIG, TAG, top_targets, ncol=5)
    P.fig_per_module(xy, obs, FIG, TAG, module_map, ncol=4)
    if dir_map:
        P.fig_per_module(xy, obs, FIG, f"{TAG}_bydirection", dir_map, ncol=2)

    # manifest fragment
    figs = sorted([f for f in os.listdir(FIG) if f.endswith(".png") and TAG in f])
    json.dump({"tag":TAG,"n_cells":int(adata.n_obs),"figures":figs,
               "top_targets":top_targets,"modules":{k:len(v) for k,v in module_map.items()}},
              open(f"{OUT}/checkpoints/wave1_manifest.json","w"), indent=2)
    P.log(f"DONE wave1 plotting in {time.time()-t0:.0f}s; {len(figs)} figures")

if __name__ == "__main__":
    main()
