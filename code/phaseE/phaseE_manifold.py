#!/usr/bin/env python
# Phase E E3 state-manifold. Build the leiden manifold on CONTROL cells (NTC+background) of the
# D4_Stim8hr subset, then compute per-target KD-vs-NTC cluster-occupancy shift + program scores.
# Manifold recipe (methods): 2000 HVG (seurat) -> scale -> PCA(50,arpack) -> neighbors(15,40pc)
#   -> leiden(res1.0, igraph) -> umap(min_dist0.3).
# Usage: phaseE_manifold.py <subset_h5ad> <outdir>
import sys, json, os, time, numpy as np, pandas as pd, scanpy as sc
SUB, OUTDIR = sys.argv[1:3]
os.makedirs(OUTDIR, exist_ok=True)
t0=time.time()
def log(*a): print("[%7.1fs]"%(time.time()-t0),*a,flush=True)
sc.settings.n_jobs=4

TARGETS=["LAT","PLCG1","CD247","CD3E","ZAP70","VAV1","SMARCE1",
 "SGF29","MED24","MED12","TADA2B","NSD1","CHD4","SMARCB1","STAT6","TRIP12","ARNT","SEL1L","SIK3"]
PROGRAMS={"naive_memory":["CCR7","SELL","TCF7","LEF1","IL7R"],
 "effector":["GZMB","GZMA","PRF1","IFNG","NKG7","GNLY","CCL5"],
 "activation":["IL2RA","CD69","TNFRSF9","MKI67","TNFRSF18","ICOS"],
 "treg":["FOXP3","IKZF2","CTLA4","IL2RA","TNFRSF18"],
 "cytokine":["IL2","IFNG","TNF","IL4","IL13","IL21","CSF2"],
 "exhaustion":["PDCD1","LAG3","HAVCR2","TIGIT","TOX","ENTPD1"]}

A=sc.read_h5ad(SUB); log("loaded",A.shape)
cls=A.obs["cell_class"].astype(str).values
# control cells = NTC + background
ctrl_mask=np.isin(cls,["NTC","background"])
Ac=A[ctrl_mask].copy(); log("control cells",Ac.shape)

# manifold on controls
sc.pp.highly_variable_genes(Ac, n_top_genes=2000, flavor="seurat")
Ac.raw=Ac; Ac=Ac[:,Ac.var.highly_variable].copy()
# capture control HVG log-norm mean/std BEFORE scaling (for projecting KD cells later)
_Xc=np.asarray(Ac.X.todense()) if hasattr(Ac.X,"todense") else np.asarray(Ac.X)
ctrl_mean=_Xc.mean(0).astype(np.float64); ctrl_std=_Xc.std(0).astype(np.float64)
sc.pp.scale(Ac, max_value=10)
sc.tl.pca(Ac, n_comps=50, svd_solver="arpack")
sc.pp.neighbors(Ac, n_neighbors=15, n_pcs=40)
sc.tl.leiden(Ac, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(Ac, min_dist=0.3)
n_clusters=Ac.obs["leiden"].nunique()
log("manifold built, leiden clusters",n_clusters)

# program scores on ALL cells (use full gene space, log-norm X)
for prog,glist in PROGRAMS.items():
    present=[g for g in glist if g in A.var["gene_name"].values]
    sc.tl.score_genes(A, present, score_name=f"prog_{prog}", use_raw=False)

# assign ALL cells (incl KD) to nearest control-cluster centroid in PCA space.
# Project A onto the CONTROL PCA: scale KD cells with the CONTROL mean/std (stored by sc.pp.scale
# on Ac), then multiply by the control PCA loadings. This keeps KD cells in the control basis.
hvg=Ac.var_names
Xall=A[:,hvg].copy()
Xd=np.asarray(Xall.X.todense()) if hasattr(Xall.X,"todense") else np.asarray(Xall.X)
Xd=(Xd-ctrl_mean)/np.where(ctrl_std>0,ctrl_std,1.0)
Xd=np.clip(Xd, -10, 10)  # match max_value=10
loadings=Ac.varm["PCs"][:,:40]  # genes x 40
Xpca_all=Xd @ loadings  # cells x 40
# control centroids
ctrl_leiden=Ac.obs["leiden"].values.astype(int)
Xpca_ctrl=Ac.obsm["X_pca"][:,:40]
cents=np.vstack([Xpca_ctrl[ctrl_leiden==k].mean(0) for k in range(n_clusters)])
# nearest centroid for all cells
from scipy.spatial.distance import cdist
d=cdist(Xpca_all, cents); assign=d.argmin(1)
A.obs["leiden_assigned"]=assign.astype(str)
log("assigned all cells to clusters")

# per-target KD-vs-NTC cluster occupancy shift
pgn=A.obs["perturbed_gene_name"].astype(str).values
ntc_mask=(cls=="NTC")
ntc_occ=np.bincount(assign[ntc_mask], minlength=n_clusters)/max(ntc_mask.sum(),1)
occ_rows=[]
for tgt in TARGETS:
    kd=(cls=="target")&(pgn==tgt)
    if kd.sum()==0: continue
    kd_occ=np.bincount(assign[kd], minlength=n_clusters)/kd.sum()
    row={"gene":tgt,"n_kd":int(kd.sum())}
    for k in range(n_clusters):
        row[f"c{k}_shift"]=round(float(kd_occ[k]-ntc_occ[k]),4)
    occ_rows.append(row)
pd.DataFrame(occ_rows).to_csv(f"{OUTDIR}/E3_cluster_occupancy_shift.csv",index=False)

# save control cluster profile (program means per cluster) + umap coords for controls
cl_prof=[]
for k in range(n_clusters):
    m=ctrl_leiden==k
    row={"cluster":k,"n_ctrl":int(m.sum())}
    # program means on control cells of this cluster
    ctrl_idx=np.where(ctrl_mask)[0][m]
    for p in PROGRAMS: row[f"prog_{p}"]=round(float(A.obs.iloc[ctrl_idx][f"prog_{p}"].mean()),4)
    cl_prof.append(row)
pd.DataFrame(cl_prof).to_csv(f"{OUTDIR}/E3_cluster_profiles.csv",index=False)

# ntc occupancy baseline
pd.DataFrame({"cluster":range(n_clusters),"ntc_occupancy":ntc_occ,
    "n_ctrl":[int((ctrl_leiden==k).sum()) for k in range(n_clusters)]}).to_csv(f"{OUTDIR}/E3_ntc_occupancy.csv",index=False)

# save umap coords (controls) with leiden for plotting
umap_df=pd.DataFrame(Ac.obsm["X_umap"],columns=["umap1","umap2"])
umap_df["leiden"]=ctrl_leiden
umap_df["cell_class"]=Ac.obs["cell_class"].values
umap_df.to_csv(f"{OUTDIR}/E3_ctrl_umap.csv",index=False)
summ={"n_clusters":int(n_clusters),"n_control":int(ctrl_mask.sum()),"n_total":int(A.n_obs),
      "ntc_occupancy":[round(float(x),4) for x in ntc_occ]}
json.dump(summ,open(f"{OUTDIR}/E3_manifold_summary.json","w"),indent=1)
log("DONE"); print("___MANIFOLD_DONE___ clusters",n_clusters)
