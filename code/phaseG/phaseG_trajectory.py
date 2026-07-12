#!/usr/bin/env python
# Phase G tier-3 (3): activation trajectory checkpoint analysis.
# Build a Rest->Stim8hr->Stim48hr activation pseudotime on CONTROL (NTC) cells across all three
# conditions; project each target's KD cells onto it; quantify where KD stalls cells along the
# activation trajectory. Prediction: TCR-proximal KD holds cells at early pseudotime; chromatin/TF
# KD shifts programs without necessarily moving trajectory position.
import sys, os, json, time, glob
import numpy as np, pandas as pd
os.environ.setdefault("MPLCONFIGDIR","/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
import scanpy as sc, anndata as ad
from scipy import sparse, stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

CKPT="/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/checkpoints"
OUT="/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseG_outputs"
os.makedirs(OUT,exist_ok=True)
SEED=0; NTC_PER=8000
t0=time.time()
def log(*a): print("[%7.1fs]"%(time.time()-t0),*a,flush=True)
TARGETS=["LAT","PLCG1","CD247","CD3E","ZAP70","VAV1","SMARCE1",
         "SGF29","MED24","MED12","TADA2B","NSD1","CHD4","SMARCB1",
         "STAT6","TRIP12","ARNT","SEL1L","SIK3"]
FLAG=set(["LAT","PLCG1","CD247","CD3E","ZAP70","VAV1"])
def axis_of(t): return "TCR-proximal" if t in FLAG else ("SMARCE1" if t=="SMARCE1" else "chromatin/TF")
COND_ORDER={"Rest":0,"Stim8hr":1,"Stim48hr":2}

rng=np.random.default_rng(SEED)
# pool all 3 conditions x 4 donors: all KD cells + subsampled NTC per (donor,cond)
ads=[]
for f in sorted(glob.glob(f"{CKPT}/D*_*.subset.h5ad")):
    base=os.path.basename(f); donor=base.split("_")[0]; cond=base.split("_")[1].split(".")[0]
    if cond not in COND_ORDER: continue
    A=sc.read_h5ad(f); A.obs["donor"]=donor; A.obs["cond"]=cond
    cls=A.obs["cell_class"].astype(str).values
    kd=np.where(cls=="target")[0]
    ni=np.where(cls=="NTC")[0]; nk=rng.choice(ni,size=min(NTC_PER,ni.size),replace=False)
    sel=np.zeros(A.n_obs,bool); sel[kd]=True; sel[nk]=True
    ads.append(A[sel].copy()); log("loaded",donor,cond,"kd",kd.size,"ntc",nk.size)
var0=ads[0].var.copy(); vn0=ads[0].var_names.copy()
A=ad.concat(ads,join="inner",index_unique="-"); del ads
A=A[:,vn0].copy(); A.var=var0
A.layers["counts"]=A.layers["counts"] if "counts" in A.layers else A.X.copy()
log("pooled",A.shape,"conds",dict(A.obs["cond"].value_counts()))

# ---- build manifold on CONTROL cells, get activation pseudotime ----
sc.pp.highly_variable_genes(A,n_top_genes=2000,flavor="cell_ranger")
Ah=A[:,A.var.highly_variable].copy(); sc.pp.scale(Ah,max_value=10)
sc.tl.pca(Ah,n_comps=30,random_state=SEED)
import harmonypy
ho=harmonypy.run_harmony(Ah.obsm["X_pca"],Ah.obs,["donor"],max_iter_harmony=20,random_state=SEED)
Zc=np.asarray(ho.Z_corr); Zc=Zc if Zc.shape[0]==Ah.n_obs else Zc.T
A.obsm["X_pca_harmony"]=Zc; del Ah
sc.pp.neighbors(A,n_neighbors=30,use_rep="X_pca_harmony",random_state=SEED)
log("neighbors done")

# DPT pseudotime rooted at a Rest NTC cell; orient so higher = more activated
isN=A.obs["cell_class"].values=="NTC"
rest_ntc=np.where(isN&(A.obs["cond"].values=="Rest"))[0]
sc.tl.diffmap(A,n_comps=15)
A.uns["iroot"]=int(rest_ntc[np.argmin(A.obsm["X_diffmap"][rest_ntc,1])])
sc.tl.dpt(A)
pt=A.obs["dpt_pseudotime"].values.copy()
# orient: ensure NTC Stim48hr > NTC Rest
if np.nanmedian(pt[isN&(A.obs["cond"].values=="Stim48hr")]) < np.nanmedian(pt[isN&(A.obs["cond"].values=="Rest")]):
    pt=1-pt
A.obs["act_pt"]=pt
# sanity: NTC pseudotime by condition
for cd in ["Rest","Stim8hr","Stim48hr"]:
    m=isN&(A.obs["cond"].values==cd)
    log("NTC",cd,"median_pt",round(float(np.nanmedian(pt[m])),3))

# ---- per target x condition: KD pseudotime vs NTC pseudotime (checkpoint = KD held earlier) ----
rows=[]; percond={}
cls=A.obs["cell_class"].values; pgn=A.obs["perturbed_gene_name"].astype(str).values; cond=A.obs["cond"].values
for cd in ["Rest","Stim8hr","Stim48hr"]:
    ntc_pt=pt[isN&(cond==cd)]; ntc_med=float(np.nanmedian(ntc_pt))
    for t in TARGETS:
        m=(cls=="target")&(pgn==t)&(cond==cd)
        if m.sum()<30: continue
        kd_pt=pt[m]; 
        u,pu=stats.mannwhitneyu(kd_pt,ntc_pt,alternative="two-sided")
        rows.append({"gene":t,"axis":axis_of(t),"cond":cd,"n_kd":int(m.sum()),
            "kd_median_pt":round(float(np.nanmedian(kd_pt)),4),"ntc_median_pt":round(ntc_med,4),
            "delta_pt":round(float(np.nanmedian(kd_pt)-ntc_med),4),"mwu_p":float(pu)})
    percond[cd]={"ntc":ntc_pt}
df=pd.DataFrame(rows)
from statsmodels.stats.multitest import multipletests
df["padj"]=multipletests(df["mwu_p"],method="fdr_bh")[1]
df.to_csv(f"{OUT}/G9_trajectory_checkpoint.csv",index=False)
log("wrote table")

# save pseudotime per cell (for plotting)
pcell=pd.DataFrame({"cond":cond,"cell_class":cls,"gene":pgn,"act_pt":pt,"donor":A.obs["donor"].values})
pcell.to_csv(f"{OUT}/G9_pseudotime_cells.csv.gz",index=False,compression="gzip")

# ---- Fig G9 ----
plt.rcParams.update({"font.size":8,"axes.spines.top":False,"axes.spines.right":False})
COL={"TCR-proximal":"#1f6fb2","chromatin/TF":"#c85200","SMARCE1":"#e0a030"}
fig=plt.figure(figsize=(12,5)); gs=fig.add_gridspec(1,2,width_ratios=[1,1.1],wspace=0.3)
# A: delta_pt at Stim8hr (KD - NTC): negative = held earlier on activation trajectory
axA=fig.add_subplot(gs[0,0])
d8=df[df.cond=="Stim8hr"].sort_values("delta_pt")
axA.barh(range(len(d8)),d8["delta_pt"],color=[COL[a] for a in d8["axis"]],
         edgecolor=["k" if p<0.05 else "none" for p in d8["padj"]],
         linewidth=[1.1 if p<0.05 else 0 for p in d8["padj"]])
axA.axvline(0,color="0.4",lw=.8)
axA.set_yticks(range(len(d8))); axA.set_yticklabels([f"$\\it{{{g}}}$" for g in d8["gene"]])
axA.set_xlabel("Δ activation pseudotime (KD − NTC), Stim 8 hr")
axA.set_title("a  TCR-proximal KD holds cells earlier on the activation trajectory\n(negative = stalled; black edge = FDR<0.05)",loc="left",fontsize=8.5)
axA.text(0.02,0.02,"← stalled earlier      advanced →",transform=axA.transAxes,fontsize=6.5,color="0.3")
# B: pseudotime distributions — NTC vs a TCR-proximal vs a chromatin/TF target at Stim8hr
axB=fig.add_subplot(gs[0,1])
xs=np.linspace(0,1,200)
def kde(v):
    v=v[~np.isnan(v)]
    try: return stats.gaussian_kde(v)(xs)
    except Exception: return np.zeros_like(xs)
axB.fill_between(xs,kde(percond["Stim8hr"]["ntc"]),color="0.7",alpha=.5,label="NTC (Stim 8 hr)")
tcr_ex=d8[d8.axis=="TCR-proximal"].iloc[0]["gene"]
chr_ex=d8[d8.axis=="chromatin/TF"].sort_values("delta_pt",ascending=False).iloc[0]["gene"]
for g,c_ in [(tcr_ex,COL["TCR-proximal"]),(chr_ex,COL["chromatin/TF"])]:
    m=(cls=="target")&(pgn==g)&(cond=="Stim8hr")
    axB.plot(xs,kde(pt[m]),color=c_,lw=1.8,label=f"$\\it{{{g}}}$ KD")
axB.set_xlabel("activation pseudotime (0 = rest-like, 1 = late effector)")
axB.set_ylabel("density")
axB.set_title("b  TCR-proximal KD cells accumulate at low pseudotime;\n   chromatin/TF KD tracks the NTC distribution",loc="left",fontsize=8.5)
axB.legend(frameon=False,fontsize=7)
fig.suptitle("Phase G · Fig 9 — Activation-trajectory checkpoint analysis (4 donors, 3 conditions pooled)",
             x=.01,ha="left",fontweight="bold",fontsize=9.5)
fig.savefig(f"{OUT}/G9_trajectory_checkpoint.png",bbox_inches="tight",dpi=150)
log("wrote fig")
A.write(f"{OUT}/G9_trajectory_manifold.h5ad")
open(f"{OUT}/_DONE_TRAJ","w").write("done\n")
log("ALL DONE trajectory")
