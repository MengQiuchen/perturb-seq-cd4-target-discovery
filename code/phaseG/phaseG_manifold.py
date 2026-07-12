#!/usr/bin/env python
# Phase G second-tier — pooled 4-donor reference manifold (harmony-integrated).
# Substrate for: (1) neighborhood DA, (2) rare-state (Treg/Tfh) effects, (3) cross-donor
# reproducibility of KD state redistribution.  Reuses E3 recipe (2000 HVG -> PCA30 -> UMAP ->
# Leiden) + harmony(donor) for cross-donor integration.  Saves a checkpoint h5ad.
import sys, os, json, time, glob
import numpy as np, pandas as pd
os.environ.setdefault("MPLCONFIGDIR","/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
import scanpy as sc, anndata as ad
from scipy import sparse
sc.settings.verbosity=1

COND = sys.argv[1] if len(sys.argv)>1 else "Stim8hr"
NTC_PER_DONOR = int(sys.argv[2]) if len(sys.argv)>2 else 12000
SEED=0
CKPT="/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/checkpoints"
OUT="/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseG_outputs"
os.makedirs(OUT, exist_ok=True)
t0=time.time()
def log(*a): print("[%7.1fs]"%(time.time()-t0),*a,flush=True)

TARGETS=["LAT","PLCG1","CD247","CD3E","ZAP70","VAV1","SMARCE1",
         "SGF29","MED24","MED12","TADA2B","NSD1","CHD4","SMARCB1",
         "STAT6","TRIP12","ARNT","SEL1L","SIK3"]
# E3 program marker sets + Tfh (rare-state) added for the safety/rare-state question
PROGRAMS={"naive_memory":["CCR7","SELL","TCF7","LEF1","IL7R"],
 "effector":["GZMB","GZMA","PRF1","IFNG","NKG7","GNLY","CCL5"],
 "activation":["IL2RA","CD69","TNFRSF9","MKI67","TNFRSF18","ICOS"],
 "treg":["FOXP3","IKZF2","CTLA4","IL2RA","TNFRSF18"],
 "cytokine":["IL2","IFNG","TNF","IL4","IL13","IL21","CSF2"],
 "exhaustion":["PDCD1","LAG3","HAVCR2","TIGIT","TOX","ENTPD1"],
 "tfh":["BCL6","CXCR5","PDCD1","ICOS","IL21","CD200"]}

rng=np.random.default_rng(SEED)
files=sorted(glob.glob(f"{CKPT}/D*_{COND}.subset.h5ad"))
assert files, f"no subsets for {COND}"
ads=[]
for f in files:
    donor=os.path.basename(f).split("_")[0]
    A=sc.read_h5ad(f); A.obs["donor"]=donor
    cls=A.obs["cell_class"].astype(str).values
    keep_kd = cls=="target"
    ntc_idx=np.where(cls=="NTC")[0]
    take=rng.choice(ntc_idx, size=min(NTC_PER_DONOR,ntc_idx.size), replace=False)
    sel=np.zeros(A.n_obs,bool); sel[keep_kd]=True; sel[take]=True
    sub=A[sel].copy(); del A
    ads.append(sub); log("loaded",donor,"kd",int(keep_kd.sum()),"ntc",take.size,"->",sub.shape)
var0=ads[0].var.copy(); vn0=ads[0].var_names.copy()
A=ad.concat(ads, join="inner", index_unique="-"); del ads
A=A[:, vn0].copy(); A.var=var0
A.layers["counts"]=A.layers["counts"] if "counts" in A.layers else A.X.copy()
log("pooled manifold set", A.shape, "| classes", dict(A.obs["cell_class"].value_counts()))

# --- E3 recipe + harmony --- (HVG on log-normalized X, cell_ranger flavor = E3's own recipe;
# avoids seurat_v3's skmisc dependency, which is absent in the fixed env)
sc.pp.highly_variable_genes(A, n_top_genes=2000, flavor="cell_ranger")
log("HVG done")
Ah=A[:, A.var.highly_variable].copy()
sc.pp.scale(Ah, max_value=10)
sc.tl.pca(Ah, n_comps=30, random_state=SEED)
log("PCA done")
# harmonypy 2.0.0: run directly; Z_corr is already (n_cells, n_pcs) — scanpy's wrapper
# wrongly transposes it, hence the shape error. Assign without transpose.
import harmonypy
ho=harmonypy.run_harmony(Ah.obsm["X_pca"], Ah.obs, ["donor"], max_iter_harmony=20, random_state=SEED)
Zc=np.asarray(ho.Z_corr)
if Zc.shape[0]!=Ah.n_obs: Zc=Zc.T          # be robust to either orientation
assert Zc.shape==(Ah.n_obs,30), f"harmony shape {Zc.shape}"
Ah.obsm["X_pca_harmony"]=Zc
log("harmony done", Zc.shape)
A.obsm["X_pca"]=Ah.obsm["X_pca"]; A.obsm["X_pca_harmony"]=Ah.obsm["X_pca_harmony"]; del Ah
sc.pp.neighbors(A, n_neighbors=30, use_rep="X_pca_harmony", random_state=SEED)
log("neighbors done")
sc.tl.leiden(A, resolution=1.0, random_state=SEED, flavor="igraph", n_iterations=2, directed=False)
sc.tl.umap(A, random_state=SEED)
log("leiden+umap done | n_clusters", A.obs["leiden"].nunique())

# --- program scores (symbol space, present genes only) ---
sym=A.var["gene_name"].astype(str).values
present={s:i for i,s in enumerate(sym)}
for prog,glist in PROGRAMS.items():
    pg=[g for g in glist if g in present]
    sc.tl.score_genes(A, pg, score_name=f"prog_{prog}", use_raw=False)
log("program scores done")

# --- cluster profiles: mean program score per cluster + size ---
prof=[]
for cl,idx in A.obs.groupby("leiden").groups.items():
    row={"leiden":cl,"n":len(idx)}
    sub=A.obs.loc[idx]
    for p in PROGRAMS: row[f"prog_{p}"]=round(float(sub[f"prog_{p}"].mean()),4)
    row["ntc_frac"]=round(float((sub["cell_class"]=="NTC").mean()),4)
    prof.append(row)
prof=pd.DataFrame(prof).sort_values("leiden")
prof.to_csv(f"{OUT}/G4_cluster_profiles_{COND}.csv",index=False)
log("cluster profiles:\n"+prof.to_string(index=False))

# --- per-target x per-donor x cluster occupancy (for cross-donor DA) ---
occ=[]
clusters=sorted(A.obs["leiden"].unique(), key=lambda x:int(x))
for donor in sorted(A.obs["donor"].unique()):
    dsel=A.obs["donor"]==donor
    ntc=A.obs[dsel & (A.obs["cell_class"]=="NTC")]
    ntc_occ=ntc["leiden"].value_counts(normalize=True)
    for t in TARGETS:
        kd=A.obs[dsel & (A.obs["cell_class"]=="target") & (A.obs["perturbed_gene_name"].astype(str)==t)]
        if len(kd)<10: continue
        kd_occ=kd["leiden"].value_counts(normalize=True)
        for cl in clusters:
            occ.append({"donor":donor,"gene":t,"leiden":cl,"n_kd":len(kd),
                        "kd_frac":round(float(kd_occ.get(cl,0.0)),5),
                        "ntc_frac":round(float(ntc_occ.get(cl,0.0)),5),
                        "shift":round(float(kd_occ.get(cl,0.0)-ntc_occ.get(cl,0.0)),5)})
occ=pd.DataFrame(occ)
occ.to_csv(f"{OUT}/G4_occupancy_by_donor_{COND}.csv",index=False)
log("occupancy rows",len(occ))

# --- save UMAP coords (subsample for plotting) + checkpoint ---
umap=pd.DataFrame(A.obsm["X_umap"],columns=["umap1","umap2"],index=A.obs_names)
umap["leiden"]=A.obs["leiden"].values; umap["cell_class"]=A.obs["cell_class"].values
umap["donor"]=A.obs["donor"].values; umap["perturbed_gene_name"]=A.obs["perturbed_gene_name"].astype(str).values
for p in PROGRAMS: umap[f"prog_{p}"]=A.obs[f"prog_{p}"].values
umap.to_csv(f"{OUT}/G4_umap_{COND}.csv.gz",index=False,compression="gzip")
log("umap saved",umap.shape)

# lean checkpoint: keep X (lognorm) + obs + obsm + var; drop counts to shrink
A.write(f"{OUT}/G4_manifold_{COND}.h5ad")
log("CHECKPOINT written", f"{OUT}/G4_manifold_{COND}.h5ad")
open(f"{OUT}/_DONE_MANIFOLD_{COND}","w").write("done\n")
log("ALL DONE manifold",COND)
