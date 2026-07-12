#!/usr/bin/env python
# Phase G tier-3 (1): co-expression network rewiring, KD vs NTC.
# For each target: within its signature-gene module, compute the gene-gene correlation matrix
# in KD cells vs NTC cells; quantify how much the co-regulation structure is disrupted.
# pseudobulk (one sample per group) cannot compute a correlation structure at all.
import sys, os, json, time, glob
import numpy as np, pandas as pd
os.environ.setdefault("MPLCONFIGDIR","/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
import scanpy as sc, anndata as ad
from scipy import sparse, stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

COND=sys.argv[1] if len(sys.argv)>1 else "Stim8hr"
CKPT="/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/checkpoints"
SIGF="/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/phaseE_signatures.json"
OUT="/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseG_outputs"
os.makedirs(OUT,exist_ok=True)
t0=time.time()
def log(*a): print("[%7.1fs]"%(time.time()-t0),*a,flush=True)

TARGETS=["LAT","PLCG1","CD247","CD3E","ZAP70","VAV1","SMARCE1",
         "SGF29","MED24","MED12","TADA2B","NSD1","CHD4","SMARCB1",
         "STAT6","TRIP12","ARNT","SEL1L","SIK3"]
FLAG=set(["LAT","PLCG1","CD247","CD3E","ZAP70","VAV1"])
def axis_of(t): return "TCR-proximal" if t in FLAG else ("SMARCE1" if t=="SMARCE1" else "chromatin/TF")
sig=json.load(open(SIGF)); SIG=sig["signatures"]; G2E=sig.get("gene_to_ensg",{})
NTOP=30        # module size: top-30 signature genes (excl. target self)
N_COMMON=140   # FIXED cell count for every target's corr matrices -> removes cell-count confound
NPERM=300      # permutation test for connectivity change

files=sorted(glob.glob(f"{CKPT}/D*_{COND}.subset.h5ad"))
ads=[]
for f in files:
    A=sc.read_h5ad(f); A.obs["donor"]=os.path.basename(f).split("_")[0]; ads.append(A)
var0=ads[0].var.copy(); vn0=ads[0].var_names.copy()
A=ad.concat(ads,join="inner",index_unique="-"); del ads
A=A[:,vn0].copy(); A.var=var0
log("pooled",A.shape)
ensg=A.var["ensg"].astype(str).values; ensg2col={}; [ensg2col.setdefault(e,i) for i,e in enumerate(ensg)]
Xc=A.X.tocsc() if sparse.issparse(A.X) else A.X
cls=A.obs["cell_class"].values; pgn=A.obs["perturbed_gene_name"].astype(str).values
ntc_idx=np.where(cls=="NTC")[0]

def dense(cols, rows):
    sub=Xc[:,cols][rows,:]
    return np.asarray(sub.todense()) if sparse.issparse(sub) else np.asarray(sub)

def corr_offdiag(M):
    # M: cells x genes ; return vector of off-diagonal Pearson r (upper triangle)
    C=np.corrcoef(M,rowvar=False)
    iu=np.triu_indices(C.shape[0],k=1)
    return C[iu], C

def mean_abs_offdiag(M):
    C=np.corrcoef(M,rowvar=False); iu=np.triu_indices(C.shape[0],k=1)
    v=C[iu]; return float(np.nanmean(np.abs(v))), v, C

rng=np.random.default_rng(0)
rows=[]; heat={}
for t in TARGETS:
    key=f"{t}_{COND}"
    if key not in SIG: continue
    genes=SIG[key]["top_genes"][:NTOP]; excl={G2E.get(t)}
    cols=[ensg2col[g] for g in genes if g in ensg2col and g not in excl]
    if len(cols)<10: log("skip",t,"few module genes"); continue
    kd_idx=np.where((cls=="target")&(pgn==t))[0]
    powered = kd_idx.size>=N_COMMON
    # FIXED N for both groups -> identical sampling noise across all targets
    n=N_COMMON if powered else kd_idx.size
    kd_s=rng.choice(kd_idx,size=n,replace=False) if kd_idx.size>=n else kd_idx
    ntc_s=rng.choice(ntc_idx,size=n,replace=False)
    Mkd=dense(cols,kd_s); Mntc=dense(cols,ntc_s)
    gg=(Mkd.std(0)>1e-6)&(Mntc.std(0)>1e-6)
    Mkd=Mkd[:,gg]; Mntc=Mntc[:,gg]
    makd,rkd,Ckd=mean_abs_offdiag(Mkd); mantc,rntc,Cntc=mean_abs_offdiag(Mntc)
    ok=~(np.isnan(rkd)|np.isnan(rntc))
    struct_preserved=float(np.corrcoef(rntc[ok],rkd[ok])[0,1]) if ok.sum()>3 else np.nan
    dC=makd-mantc
    # permutation null for dC: pool the 2n cells, split into two n-groups, recompute
    pool=np.vstack([Mkd,Mntc]); nulls=np.empty(NPERM)
    for b in range(NPERM):
        perm=rng.permutation(pool.shape[0])
        a1,a2=pool[perm[:n]],pool[perm[n:2*n]]
        m1,_,_=mean_abs_offdiag(a1); m2,_,_=mean_abs_offdiag(a2); nulls[b]=m1-m2
    p_perm=float((np.abs(nulls)>=abs(dC)).mean())
    rows.append({"gene":t,"axis":axis_of(t),"n_kd":int(kd_idx.size),"n_used":int(n),
        "well_powered":bool(powered),"n_module_genes":int(gg.sum()),
        "mean_abs_r_ntc":round(mantc,4),"mean_abs_r_kd":round(makd,4),
        "delta_connectivity":round(dC,4),"perm_p":p_perm,
        "struct_preserved_r":round(struct_preserved,4)})
    heat[t]={"Cntc":Cntc,"Ckd":Ckd}
    log("done",t,"n",n,"pow",powered,"dC",round(dC,3),"p",round(p_perm,3),"struct",round(struct_preserved,3))

df=pd.DataFrame(rows)
from statsmodels.stats.multitest import multipletests
df["perm_padj"]=multipletests(df["perm_p"],method="fdr_bh")[1]
df["sig_rewire"]=(df["perm_padj"]<0.05)&(df["well_powered"])
df.to_csv(f"{OUT}/G8_coexpr_rewiring_{COND}.csv",index=False)
log("wrote table")

# ---- Fig G8: rewiring summary + example matrices ----
plt.rcParams.update({"font.size":8,"axes.spines.top":False,"axes.spines.right":False})
COL={"TCR-proximal":"#1f6fb2","chromatin/TF":"#c85200","SMARCE1":"#e0a030"}
fig=plt.figure(figsize=(12,6)); gs=fig.add_gridspec(2,3,height_ratios=[1.2,1],hspace=0.5,wspace=0.35)
# A: connectivity CHANGE at fixed N (permutation-tested); positive = KD tightens co-expression
axA=fig.add_subplot(gs[0,:])
d2=df.sort_values("delta_connectivity")
axA.barh(range(len(d2)),d2["delta_connectivity"],color=[COL[a] for a in d2["axis"]],
         edgecolor=["k" if s else "none" for s in d2["sig_rewire"]],
         linewidth=[1.2 if s else 0 for s in d2["sig_rewire"]],
         hatch=["" if p else "///" for p in d2["well_powered"]])
axA.axvline(0,color="0.4",lw=.8)
axA.set_yticks(range(len(d2))); axA.set_yticklabels([f"$\\it{{{g}}}$" for g in d2["gene"]])
axA.set_xlabel(f"Δ mean |gene-gene r| in module (KD − NTC), fixed N={N_COMMON}   → KD increases co-expression coupling")
axA.set_title("a  Knockdown reshapes the downstream co-expression module at matched cell count\n(black edge = permutation FDR<0.05; hatched = <N, underpowered)",loc="left",fontsize=8.5)
from matplotlib.patches import Patch
axA.legend(handles=[Patch(color=COL[a],label=a) for a in ["TCR-proximal","chromatin/TF","SMARCE1"]],frameon=False,fontsize=6.5,loc="lower right")
# B-D: example correlation matrices (NTC vs KD) for the target with largest significant change
sigd=d2[d2.sig_rewire]
ex=(sigd.reindex(sigd["delta_connectivity"].abs().sort_values(ascending=False).index).iloc[0]["gene"]
    if len(sigd) else d2.iloc[-1]["gene"])
for j,(lab,key) in enumerate([("NTC (control)","Cntc"),(f"{ex} KD","Ckd")]):
    ax=fig.add_subplot(gs[1,j]); C=heat[ex][key]
    im=ax.imshow(C,cmap="RdBu_r",vmin=-1,vmax=1); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{lab}",fontsize=7.5)
# difference
axd=fig.add_subplot(gs[1,2]); D=heat[ex]["Ckd"]-heat[ex]["Cntc"]
imd=axd.imshow(D,cmap="PuOr",vmin=-1,vmax=1); axd.set_xticks([]); axd.set_yticks([])
axd.set_title(f"Δ (KD − NTC)",fontsize=7.5)
plt.colorbar(imd,ax=axd,shrink=0.7)
fig.suptitle(f"Phase G · Fig 8 — Co-expression module rewiring by knockdown ({COND}, 4 donors, matched N={N_COMMON}; example: $\\it{{{ex}}}$)",
             x=.01,ha="left",fontweight="bold",fontsize=9.5)
fig.savefig(f"{OUT}/G8_coexpr_rewiring_{COND}.png",bbox_inches="tight",dpi=150)
log("wrote fig")
open(f"{OUT}/_DONE_COEXPR_{COND}","w").write("done\n")
log("ALL DONE coexpr",COND)
