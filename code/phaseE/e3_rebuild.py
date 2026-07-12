#!/usr/bin/env python
# Rebuild the E3 state-manifold INDEPENDENTLY per donor and save per-cell cluster assignment,
# so cross-donor axis-recovery + a coverage-downsampling control can be computed identically.
# Manifold recipe identical to phaseE_manifold.py: 2000 HVG(seurat) -> scale(10) -> PCA(50,arpack)
#   -> neighbors(15,40pc) -> leiden(res1.0,igraph).  (UMAP skipped: not needed for occupancy.)
# Usage: e3_rebuild.py <subset_h5ad> <donor_tag> <outdir>
import sys, json, os, time, numpy as np, pandas as pd, scanpy as sc
SUB, TAG, OUTDIR = sys.argv[1:4]
os.makedirs(OUTDIR, exist_ok=True)
t0=time.time()
def log(*a): print("[%7.1fs]"%(time.time()-t0),TAG,*a,flush=True)
sc.settings.n_jobs=4

A=sc.read_h5ad(SUB); log("loaded",A.shape)
cls=A.obs["cell_class"].astype(str).values
ctrl_mask=np.isin(cls,["NTC","background"])
Ac=A[ctrl_mask].copy(); log("control cells",Ac.shape)

sc.pp.highly_variable_genes(Ac, n_top_genes=2000, flavor="seurat")
Ac.raw=Ac; Ac=Ac[:,Ac.var.highly_variable].copy()
_Xc=np.asarray(Ac.X.todense()) if hasattr(Ac.X,"todense") else np.asarray(Ac.X)
ctrl_mean=_Xc.mean(0).astype(np.float64); ctrl_std=_Xc.std(0).astype(np.float64)
sc.pp.scale(Ac, max_value=10)
sc.tl.pca(Ac, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(Ac, n_neighbors=15, n_pcs=40)
sc.tl.leiden(Ac, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
n_clusters=Ac.obs["leiden"].nunique()
log("manifold built, leiden clusters",n_clusters)

hvg=Ac.var_names
Xall=A[:,hvg].copy()
Xd=np.asarray(Xall.X.todense()) if hasattr(Xall.X,"todense") else np.asarray(Xall.X)
Xd=(Xd-ctrl_mean)/np.where(ctrl_std>0,ctrl_std,1.0)
Xd=np.clip(Xd,-10,10)
loadings=Ac.varm["PCs"][:,:40]
Xpca_all=Xd @ loadings
ctrl_leiden=Ac.obs["leiden"].values.astype(int)
Xpca_ctrl=Ac.obsm["X_pca"][:,:40]
cents=np.vstack([Xpca_ctrl[ctrl_leiden==k].mean(0) for k in range(n_clusters)])
from scipy.spatial.distance import cdist
assign=cdist(Xpca_all,cents).argmin(1)
log("assigned all cells")

# save compact per-cell assignment for downstream (occupancy, redistribution, downsampling)
out=pd.DataFrame({"cell_class":cls,
                  "perturbed_gene_name":A.obs["perturbed_gene_name"].astype(str).values,
                  "leiden_assigned":assign.astype(int)})
out.to_csv(f"{OUTDIR}/{TAG}.cell_assign.csv",index=False)
json.dump({"donor":TAG,"n_clusters":int(n_clusters),"n_control":int(ctrl_mask.sum()),
           "n_total":int(A.n_obs)}, open(f"{OUTDIR}/{TAG}.manifold_meta.json","w"))
log("DONE n_clusters",n_clusters); print("___E3_REBUILD_DONE___",TAG,n_clusters)
