#!/usr/bin/env python
# Phase G tier-3 (2): dual-guide single-cell concordance.
# Each target has ~2 CRISPRi guides. At single-cell resolution, test whether the two guides push
# cells toward the SAME response state — a stronger on-target / anti-off-target confirmation than
# a pseudobulk correlation of two bulk profiles. Uses guide-aware re-gathered subsets (guide_id).
import sys, os, json, time, glob
import numpy as np, pandas as pd
os.environ.setdefault("MPLCONFIGDIR","/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
import scanpy as sc, anndata as ad
from scipy import sparse, stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

COND="Stim8hr"
GDIR="/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseG_outputs/guide_subsets"
SIGF="/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/phaseE_signatures.json"
OUT="/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseG_outputs"
MIN_CELLS_PER_GUIDE=25   # a guide needs >= this many cells (pooled) to enter the test
t0=time.time()
def log(*a): print("[%7.1fs]"%(time.time()-t0),*a,flush=True)
TARGETS=["LAT","PLCG1","CD247","CD3E","ZAP70","VAV1","SMARCE1",
         "SGF29","MED24","MED12","TADA2B","NSD1","CHD4","SMARCB1",
         "STAT6","TRIP12","ARNT","SEL1L","SIK3"]
FLAG=set(["LAT","PLCG1","CD247","CD3E","ZAP70","VAV1"])
def axis_of(t): return "TCR-proximal" if t in FLAG else ("SMARCE1" if t=="SMARCE1" else "chromatin/TF")
sig=json.load(open(SIGF)); SIG=sig["signatures"]; G2E=sig.get("gene_to_ensg",{})
NTOP=150

# pool 4 donors (guide-aware)
ads=[]
for f in sorted(glob.glob(f"{GDIR}/D*_{COND}.guide.h5ad")):
    A=sc.read_h5ad(f); A.obs["donor"]=os.path.basename(f).split("_")[0]; ads.append(A)
    log("loaded",os.path.basename(f),A.shape)
var0=ads[0].var.copy(); vn0=ads[0].var_names.copy()
A=ad.concat(ads,join="inner",index_unique="-"); del ads
A=A[:,vn0].copy(); A.var=var0
log("pooled",A.shape)
ensg=A.var["ensg"].astype(str).values; ensg2col={}; [ensg2col.setdefault(e,i) for i,e in enumerate(ensg)]
cls=A.obs["cell_class"].values; pgn=A.obs["perturbed_gene_name"].astype(str).values
gid=A.obs["guide_id"].astype(str).values
Xc=A.X.tocsc() if sparse.issparse(A.X) else A.X
ntc_idx=np.where(cls=="NTC")[0]

def dense(cols,rows):
    sub=Xc[:,cols][rows,:]; return np.asarray(sub.todense()) if sparse.issparse(sub) else np.asarray(sub)

# NTC baseline mean/sd per signature (for response R)
rows=[]; percell={}
for t in TARGETS:
    key=f"{t}_{COND}"
    if key not in SIG: continue
    s=SIG[key]; genes=s["top_genes"][:NTOP]; lfc=np.array(s["top_log_fc"][:NTOP],float)
    excl={G2E.get(t)}
    keep=[(g,l) for g,l in zip(genes,lfc) if g in ensg2col and g not in excl]
    cols=[ensg2col[g] for g,_ in keep]; w=np.sign([l for _,l in keep])
    # NTC stats on these genes
    Mntc=dense(cols,ntc_idx); mu=Mntc.mean(0); sd=Mntc.std(0)
    good=(sd>1e-6)&((Mntc>0).mean(0)>=0.05)
    cols=[c for c,gk in zip(cols,good) if gk]; w=w[good]; mu=mu[good]; sd=sd[good]
    def resp(rowidx):
        M=dense(cols,rowidx); Z=(M-mu)/sd; return (Z*w).mean(1)   # per-cell response
    # per-guide response for this target
    tmask=(cls=="target")&(pgn==t)
    guides=pd.Series(gid[tmask]).value_counts()
    good_guides=[g for g,cnt in guides.items() if cnt>=MIN_CELLS_PER_GUIDE]
    if len(good_guides)<2:
        rows.append({"gene":t,"axis":axis_of(t),"n_guides_tested":len(good_guides),
                     "status":"single_guide_or_low","guide_concordance":np.nan}); 
        log("skip",t,"guides>=%d:"%MIN_CELLS_PER_GUIDE,len(good_guides)); continue
    # take top-2 guides by cell count
    g1,g2=good_guides[0],good_guides[1]
    i1=np.where(tmask&(gid==g1))[0]; i2=np.where(tmask&(gid==g2))[0]
    r1=resp(i1); r2=resp(i2)
    # NTC response null
    rN=resp(ntc_idx[np.random.default_rng(0).choice(len(ntc_idx),size=min(2000,len(ntc_idx)),replace=False)])
    # 1) both guides shift response same direction vs NTC?
    d1=float(np.median(r1)-np.median(rN)); d2=float(np.median(r2)-np.median(rN))
    same_dir=bool(np.sign(d1)==np.sign(d2))
    # 2) are the two guides' response distributions statistically similar? (KS: high p = concordant)
    ks,ksp=stats.ks_2samp(r1,r2)
    # 3) guide-vs-guide effect-size ratio (closer to 1 = concordant)
    ratio=float(min(abs(d1),abs(d2))/max(abs(d1),abs(d2))) if max(abs(d1),abs(d2))>0 else np.nan
    # 4) each guide vs NTC significance
    _,p1=stats.mannwhitneyu(r1,rN,alternative="two-sided"); _,p2=stats.mannwhitneyu(r2,rN,alternative="two-sided")
    rows.append({"gene":t,"axis":axis_of(t),"n_guides_tested":2,"status":"tested",
        "n_g1":int(i1.size),"n_g2":int(i2.size),
        "d1_vs_ntc":round(d1,4),"d2_vs_ntc":round(d2,4),"same_direction":same_dir,
        "effect_ratio":round(ratio,4) if ratio==ratio else np.nan,
        "ks_stat":round(float(ks),4),"ks_p":float(ksp),
        "g1_sig":bool(p1<0.05),"g2_sig":bool(p2<0.05)})
    percell[t]={"r1":r1.tolist()[:2000],"r2":r2.tolist()[:2000],"rN":rN.tolist()[:2000],"g1":g1,"g2":g2}
    log("done",t,"d1",round(d1,3),"d2",round(d2,3),"same_dir",same_dir,"ratio",round(ratio,3) if ratio==ratio else None)

df=pd.DataFrame(rows)
df.to_csv(f"{OUT}/G10_dualguide_{COND}.csv",index=False)
json.dump(percell,open(f"{OUT}/G10_dualguide_percell_{COND}.json","w"))
log("wrote table; tested",int((df.status=='tested').sum()),"targets")

# ---- Fig G10 ----
plt.rcParams.update({"font.size":8,"axes.spines.top":False,"axes.spines.right":False})
COL={"TCR-proximal":"#1f6fb2","chromatin/TF":"#c85200","SMARCE1":"#e0a030"}
tt=df[df.status=="tested"].copy()
fig=plt.figure(figsize=(12,5)); gs=fig.add_gridspec(1,2,width_ratios=[1,1],wspace=0.32)
# A: guide1 vs guide2 effect (scatter) — on-diagonal = concordant
axA=fig.add_subplot(gs[0,0])
for _,r in tt.iterrows():
    axA.scatter(r["d1_vs_ntc"],r["d2_vs_ntc"],color=COL[r["axis"]],s=45,zorder=3,
                edgecolor="k" if (r["g1_sig"] and r["g2_sig"]) else "none",linewidth=0.8)
    axA.annotate(r["gene"],(r["d1_vs_ntc"],r["d2_vs_ntc"]),fontsize=6,xytext=(3,3),textcoords="offset points")
lim=np.nanmax(np.abs(np.r_[tt.d1_vs_ntc,tt.d2_vs_ntc]))*1.15
axA.plot([-lim,lim],[-lim,lim],"--",color="0.6",lw=.8)
axA.axhline(0,color="0.8",lw=.6); axA.axvline(0,color="0.8",lw=.6)
axA.set_xlim(-lim,lim); axA.set_ylim(-lim,lim)
axA.set_xlabel("guide 1 response shift vs NTC"); axA.set_ylabel("guide 2 response shift vs NTC")
axA.set_title("a  The two guides agree in direction & magnitude\n(on-diagonal = concordant; black edge = both guides sig)",loc="left",fontsize=8.5)
from matplotlib.patches import Patch
axA.legend(handles=[Patch(color=COL[a],label=a) for a in ["TCR-proximal","chromatin/TF","SMARCE1"]],frameon=False,fontsize=6.5,loc="lower right")
# B: example per-cell distributions for a flagship target
axB=fig.add_subplot(gs[0,1])
ex=[g for g in ["CD247","CD3E","LAT","PLCG1"] if g in percell]
ex=ex[0] if ex else list(percell.keys())[0]
xs=np.linspace(min(percell[ex]["rN"]+percell[ex]["r1"]+percell[ex]["r2"]),
               max(percell[ex]["rN"]+percell[ex]["r1"]+percell[ex]["r2"]),200)
def kde(v):
    v=np.array(v); 
    try: return stats.gaussian_kde(v)(xs)
    except Exception: return np.zeros_like(xs)
axB.fill_between(xs,kde(percell[ex]["rN"]),color="0.7",alpha=.5,label="NTC")
axB.plot(xs,kde(percell[ex]["r1"]),color=COL[axis_of(ex)],lw=1.8,label=f"guide 1 (n={df[df.gene==ex].n_g1.iloc[0]})")
axB.plot(xs,kde(percell[ex]["r2"]),color=COL[axis_of(ex)],lw=1.8,ls="--",label=f"guide 2 (n={df[df.gene==ex].n_g2.iloc[0]})")
axB.set_xlabel("per-cell response score"); axB.set_ylabel("density")
axB.set_title(f"b  Both guides of $\\it{{{ex}}}$ drive the same single-cell response",loc="left",fontsize=8.5)
axB.legend(frameon=False,fontsize=7)
fig.suptitle(f"Phase G · Fig 10 — Dual-guide single-cell concordance ({COND}, 4 donors, guide-resolved)",
             x=.01,ha="left",fontweight="bold",fontsize=9.5)
fig.savefig(f"{OUT}/G10_dualguide_{COND}.png",bbox_inches="tight",dpi=150)
log("wrote fig")
open(f"{OUT}/_DONE_DUALGUIDE","w").write("done\n")
log("ALL DONE dualguide")
