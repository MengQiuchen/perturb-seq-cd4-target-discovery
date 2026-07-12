#!/usr/bin/env python
# Phase G - First-tier single-cell analyses on already-nominated targets.
# (1) Knockdown dose-response  (2) Response heterogeneity / responder fraction
# (3) Baseline-state sensitivity.  Pools the 4 donors for one condition.
# Reuses the EXACT DE_stats 150-gene signatures + NTC standardization + 6-program marker sets.
import sys, os, json, time, glob
import numpy as np, pandas as pd
os.environ.setdefault("MPLCONFIGDIR","/mnt/scratche/slow/ghlab/qiuchen/tmp/mpl")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
import scanpy as sc, anndata as ad
from scipy import sparse, stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

COND = sys.argv[1] if len(sys.argv)>1 else "Stim8hr"
CKPT="/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/checkpoints"
SIGF="/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/phaseE_signatures.json"
OUT="/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseG_outputs"
os.makedirs(OUT, exist_ok=True)
t0=time.time()
def log(*a): print("[%7.1fs]"%(time.time()-t0),*a,flush=True)

TARGETS=["LAT","PLCG1","CD247","CD3E","ZAP70","VAV1","SMARCE1",
         "SGF29","MED24","MED12","TADA2B","NSD1","CHD4","SMARCB1",
         "STAT6","TRIP12","ARNT","SEL1L","SIK3"]
FLAGSHIP=set(["LAT","PLCG1","CD247","CD3E","ZAP70","VAV1"])
def axis_of(t):
    if t in FLAGSHIP: return "TCR-proximal"
    if t=="SMARCE1":  return "chromatin (SMARCE1)"
    return "chromatin/TF"
AXIS={t:axis_of(t) for t in TARGETS}
COL={"TCR-proximal":"#1f6fb2","chromatin/TF":"#c85200","chromatin (SMARCE1)":"#e0a030"}

# activation program (baseline-state axis) — canonical E3 marker set
PROG_ACT=["IL2RA","CD69","TNFRSF9","MKI67","TNFRSF18","ICOS"]

sig=json.load(open(SIGF)); SIG=sig["signatures"]; G2E=sig.get("gene_to_ensg",{})

# ---- load & pool the 4 donor subsets for this condition ----
files=sorted(glob.glob(f"{CKPT}/D*_{COND}.subset.h5ad"))
assert files, f"no subsets for {COND}"
log("files", [os.path.basename(f) for f in files])
ads=[]
for f in files:
    donor=os.path.basename(f).split("_")[0]
    A=sc.read_h5ad(f); A.obs["donor"]=donor
    ads.append(A); log("loaded",donor,A.shape)
var0=ads[0].var.copy()          # var identical across donors; concat drops it -> restore
vn0=ads[0].var_names.copy()
nv0=ads[0].n_vars
A=ad.concat(ads, join="inner", index_unique="-"); del ads
# concat preserves var order for identical var_names; reattach the annotation columns
A=A[:, vn0].copy()
A.var=var0
log("pooled", A.shape, "(var restored %d/%d)"%(A.n_vars,nv0),
    "| has gene_name:", "gene_name" in A.var, "| has ensg:", "ensg" in A.var)

sym=A.var["gene_name"].astype(str).values
ensg=A.var["ensg"].astype(str).values if "ensg" in A.var else A.var_names.astype(str).values
ensg2col={}; [ensg2col.setdefault(e,i) for i,e in enumerate(ensg)]
sym2col={};  [sym2col.setdefault(s,i) for i,s in enumerate(sym)]
sym2ensg={s:e for s,e in zip(sym,ensg)}
ACT_ENSG={sym2ensg[g] for g in PROG_ACT if g in sym2ensg}   # activation markers -> exclude from R
log("activation-marker ENSG for exclusion:", sorted(ACT_ENSG))

Xln=A.X.tocsc() if sparse.issparse(A.X) else A.X   # csc: fast column slicing
cls=A.obs["cell_class"].values
pgn=A.obs["perturbed_gene_name"].astype(str).values
donor=A.obs["donor"].values
is_ntc=cls=="NTC"; ntc_idx=np.where(is_ntc)[0]
log("cells: ntc",int(is_ntc.sum()),"| target",int((cls=='target').sum()),"| bg",int((cls=='background').sum()))

def colvec(col):
    v=Xln[:,col]
    return np.asarray(v.todense()).ravel() if sparse.issparse(v) else np.asarray(v).ravel()

# baseline activation score per cell (log-norm space), markers present only
act_cols=[sym2col[g] for g in PROG_ACT if g in sym2col]
ACT=np.zeros(A.n_obs)
if act_cols:
    M=np.column_stack([colvec(c) for c in act_cols])
    # z per gene across all cells, mean = score
    ACT=((M-M.mean(0))/(M.std(0)+1e-9)).mean(1)
log("activation score built from",len(act_cols),"markers")

rows_dr=[]; rows_resp=[]; rows_base=[]; percell={}
for t in TARGETS:
    key=f"{t}_{COND}"
    if key not in SIG or t not in sym2col:
        log("skip",t,"(no sig or not in var)"); continue
    s=SIG[key]; genes=s["top_genes"]; lfc=np.array(s["top_log_fc"],float)
    # circularity guard: drop the target's OWN ensg + activation markers from the response signature
    excl={G2E.get(t)} | ACT_ENSG
    n_excl=sum(1 for g in genes if g in excl)
    ok=np.array([(g not in excl) and (ensg2col.get(g,-1)>=0) for g in genes])
    cols=np.array([ensg2col.get(g,-1) for g in genes])[ok]; w=np.sign(lfc[ok])
    sub=Xln[:,cols]; sub=np.asarray(sub.todense()) if sparse.issparse(sub) else np.asarray(sub)
    mu=sub[ntc_idx].mean(0); sd=sub[ntc_idx].std(0); det=(sub[ntc_idx]>0).mean(0)
    good=(sd>1e-6)&(det>=0.05)                    # stable signature genes only
    subg=sub[:,good]; wg=w[good]; mug=mu[good]; sdg=sd[good]
    Z=(subg-mug)/sdg
    R=(Z*wg).mean(1)                              # per-cell response, NTC-centered ~0

    tcol=sym2col[t]; tgt=np.expm1(colvec(tcol))   # linear CP10k
    ntc_mean=tgt[ntc_idx].mean(); frac_ntc=float((tgt[ntc_idx]>0).mean())
    dose=tgt/(ntc_mean+1e-9)

    kd=(cls=="target")&(pgn==t); kd_idx=np.where(kd)[0]; nkd=kd_idx.size
    Rkd=R[kd_idx]; Rntc=R[ntc_idx]; dose_kd=dose[kd_idx]; act_kd=ACT[kd_idx]

    # (1) dose-response
    rho,prho=stats.spearmanr(dose_kd,Rkd)
    rows_dr.append({"gene":t,"axis":AXIS[t],"n_kd":int(nkd),"n_sig_used":int(good.sum()),
        "n_excl_circular":int(n_excl),
        "ntc_mean_cp10k":round(float(ntc_mean),3),"frac_expr_ntc":round(frac_ntc,3),
        "spearman_dose_R":round(float(rho),3),"p":float(prho),
        "dose_reliable":bool(frac_ntc>=0.30),
        "well_powered":bool(nkd>=150),"power_note":("ok" if nkd>=150 else "underpowered(<150)")})
    # (2) responder fraction / heterogeneity
    thr=np.percentile(Rntc,95); resp=float((Rkd>thr).mean())
    sk=stats.skew(Rkd); ku=stats.kurtosis(Rkd,fisher=True); bc=(sk**2+1)/(ku+3)
    mwu=stats.mannwhitneyu(Rkd,Rntc,alternative="greater").pvalue
    rows_resp.append({"gene":t,"axis":AXIS[t],"n_kd":int(nkd),
        "R_kd_median":round(float(np.median(Rkd)),3),"R_ntc_p95":round(float(thr),3),
        "responder_frac":round(resp,3),"R_kd_IQR":round(float(np.subtract(*np.percentile(Rkd,[75,25]))),3),
        "bimodality_coef":round(float(bc),3),"mwu_p_greater":float(mwu)})
    # (3) baseline-state sensitivity
    q1,q2=np.percentile(act_kd,[33.3,66.7]); lo=act_kd<=q1; hi=act_kd>q2
    R_lo=float(Rkd[lo].mean()) if lo.sum() else np.nan
    R_hi=float(Rkd[hi].mean()) if hi.sum() else np.nan
    rb,pb=stats.spearmanr(act_kd,Rkd)
    rows_base.append({"gene":t,"axis":AXIS[t],"n_kd":int(nkd),
        "spearman_baseAct_R":round(float(rb),3),"p":float(pb),
        "R_lowAct":round(R_lo,3),"R_highAct":round(R_hi,3),"delta_hi_lo":round(R_hi-R_lo,3)})

    percell[t]={"dose":dose_kd.astype(np.float32),"R":Rkd.astype(np.float32),
                "act":act_kd.astype(np.float32)}
    log("done",t,"nkd",nkd,"sig",int(good.sum()),"rho",round(rho,3),"resp",round(resp,3))

dr=pd.DataFrame(rows_dr); rp=pd.DataFrame(rows_resp); bs=pd.DataFrame(rows_base)
dr.to_csv(f"{OUT}/G_dose_response_{COND}.csv",index=False)
rp.to_csv(f"{OUT}/G_responder_{COND}.csv",index=False)
bs.to_csv(f"{OUT}/G_baseline_sensitivity_{COND}.csv",index=False)
# NTC-null R for figure reference
np.save(f"{OUT}/_Rntc_{COND}.npy", R[ntc_idx].astype(np.float32))
json.dump({t:{k:v.tolist() for k,v in d.items()} for t,d in percell.items()},
          open(f"{OUT}/_percell_{COND}.json","w"))
log("WROTE tables")

# ================= FIGURES =================
plt.rcParams.update({"font.size":8,"axes.spines.top":False,"axes.spines.right":False,
    "axes.titlesize":8,"axes.titleweight":"normal","figure.dpi":150})
def order_by_axis(df,val):
    return df.sort_values(["axis",val])

# ---- Fig G1: dose-response ----
fig=plt.figure(figsize=(11,6.2))
gs=fig.add_gridspec(2,3,height_ratios=[1.25,1],hspace=0.55,wspace=0.32)
# A: spearman bar
axA=fig.add_subplot(gs[0,:])
d2=dr.sort_values("spearman_dose_R").reset_index(drop=True)
cols_bar=[COL[a] for a in d2["axis"]]
axA.barh(range(len(d2)),d2["spearman_dose_R"],color=cols_bar,
         edgecolor=["none" if r else "k" for r in d2["dose_reliable"]],
         hatch=["" if r else "///" for r in d2["dose_reliable"]])
axA.axvline(0,color="0.4",lw=.8)
axA.set_yticks(range(len(d2))); axA.set_yticklabels([f"$\\it{{{g}}}$" for g in d2["gene"]])
axA.set_xlabel("Spearman ρ (per-cell residual dose → response)   ← stronger KD, more response")
axA.set_title("a  Single-cell dose–response: within KD cells, cells with stronger knockdown show a stronger downstream signature",loc="left")
for i,(v,p) in enumerate(zip(d2["spearman_dose_R"],d2["p"])):
    axA.text(v-0.01 if v<0 else v+0.01,i,("*" if p<0.05 else ""),va="center",ha="right" if v<0 else "left",fontsize=9)
from matplotlib.patches import Patch
axA.legend(handles=[Patch(color=COL["TCR-proximal"],label="TCR-proximal"),
                    Patch(color=COL["chromatin/TF"],label="chromatin/TF"),
                    Patch(color=COL["chromatin (SMARCE1)"],label="SMARCE1"),
                    Patch(facecolor="w",edgecolor="k",hatch="///",label="target lowly expressed (dose axis noisy)")],
           frameon=False,fontsize=6,loc="upper left",bbox_to_anchor=(0.01,0.99))
# B-D: example curves — pick strong flagship, strong chromatin/TF (reliable), and a lowly-expressed one
rel=dr[dr.dose_reliable]
ex=[]
f_ok=rel[rel.axis=="TCR-proximal"].sort_values("spearman_dose_R")
c_ok=rel[rel.axis=="chromatin/TF"].sort_values("spearman_dose_R")
if len(f_ok): ex.append(f_ok.iloc[0]["gene"])
if len(c_ok): ex.append(c_ok.iloc[0]["gene"])
low=dr[~dr.dose_reliable].sort_values("frac_expr_ntc")
ex.append(low.iloc[0]["gene"] if len(low) else dr.iloc[-1]["gene"])
for j,g in enumerate(ex[:3]):
    ax=fig.add_subplot(gs[1,j]); pc=percell[g]
    ax.scatter(pc["dose"],pc["R"],s=4,alpha=.18,color=COL[AXIS[g]],rasterized=True,linewidths=0)
    # decile means
    dq=np.clip(pc["dose"],0,np.percentile(pc["dose"],99))
    bins=np.linspace(0,max(dq.max(),1.0),9); idx=np.digitize(dq,bins)
    bx=[];by=[]
    for b in range(1,len(bins)):
        m=idx==b
        if m.sum()>=5: bx.append(dq[m].mean()); by.append(pc["R"][m].mean())
    ax.plot(bx,by,"-o",color="k",ms=3,lw=1.2,zorder=5)
    ax.axhline(0,color="0.6",lw=.7,ls=":")
    ax.set_title(f"$\\it{{{g}}}$  ({AXIS[g]}, n={pc['dose'].size})",loc="left",fontsize=7.5)
    ax.set_xlabel("residual target dose (1 = NTC)")
    if j==0: ax.set_ylabel("per-cell response R")
fig.suptitle(f"Phase G · Fig 1 — Single-cell knockdown dose–response ({COND}, 4 donors pooled)",
             x=.01,ha="left",fontweight="bold",fontsize=9.5)
fig.savefig(f"{OUT}/G1_dose_response_{COND}.png",bbox_inches="tight",dpi=150)
log("WROTE G1")

# ---- Fig G2: responder fraction + heterogeneity ----
fig2=plt.figure(figsize=(11,5.2))
gs2=fig2.add_gridspec(1,2,width_ratios=[1.15,1],wspace=0.28)
axA=fig2.add_subplot(gs2[0,0])
r2=rp.sort_values("responder_frac").reset_index(drop=True)
axA.barh(range(len(r2)),r2["responder_frac"],color=[COL[a] for a in r2["axis"]])
axA.axvline(0.05,color="0.4",lw=.9,ls="--"); axA.text(0.05,len(r2)-0.5,"NTC null 5%",fontsize=6,color="0.3")
axA.set_yticks(range(len(r2))); axA.set_yticklabels([f"$\\it{{{g}}}$" for g in r2["gene"]])
axA.set_xlabel("responder fraction (KD cells with R > NTC 95th pct)")
axA.set_title("a  What fraction of KD cells actually manifest the signature",loc="left")
# B: R distributions for a few targets vs NTC null
axB=fig2.add_subplot(gs2[0,1])
Rntc=np.load(f"{OUT}/_Rntc_{COND}.npy")
pick=[]
for cand in [rp.sort_values("responder_frac").iloc[-1]["gene"],
             rp[rp.axis=='chromatin/TF'].sort_values('responder_frac').iloc[-1]['gene'] if (rp.axis=='chromatin/TF').any() else None,
             rp.sort_values("responder_frac").iloc[0]["gene"]]:
    if cand and cand not in pick: pick.append(cand)
xs=np.linspace(min(Rntc.min(),-1),max([percell[g]['R'].max() for g in pick]+[Rntc.max()]),200)
def kde(v):
    try:
        k=stats.gaussian_kde(v); return k(xs)
    except Exception: return np.zeros_like(xs)
axB.fill_between(xs,kde(Rntc),color="0.7",alpha=.5,label="NTC (null)")
for g in pick:
    axB.plot(xs,kde(percell[g]["R"]),color=COL[AXIS[g]],lw=1.6,label=f"$\\it{{{g}}}$")
axB.axvline(np.percentile(Rntc,95),color="0.4",ls="--",lw=.8)
axB.set_xlabel("per-cell response R"); axB.set_ylabel("density")
axB.set_title("b  Each KD population's response distribution shifts right of the NTC null (larger for signaling than chromatin/TF)",loc="left",fontsize=7)
axB.legend(frameon=False,fontsize=6.5)
fig2.suptitle(f"Phase G · Fig 2 — Response heterogeneity across single KD cells ({COND}, 4 donors)",
              x=.01,ha="left",fontweight="bold",fontsize=9.5)
fig2.savefig(f"{OUT}/G2_responder_{COND}.png",bbox_inches="tight",dpi=150)
log("WROTE G2")

# ---- Fig G3: baseline-state sensitivity ----
fig3=plt.figure(figsize=(11,5.2))
gs3=fig3.add_gridspec(1,2,width_ratios=[1,1.05],wspace=0.3)
axA=fig3.add_subplot(gs3[0,0])
b2=bs.sort_values("spearman_baseAct_R").reset_index(drop=True)
axA.barh(range(len(b2)),b2["spearman_baseAct_R"],color=[COL[a] for a in b2["axis"]])
axA.axvline(0,color="0.4",lw=.8)
axA.set_yticks(range(len(b2))); axA.set_yticklabels([f"$\\it{{{g}}}$" for g in b2["gene"]])
axA.set_xlabel("Spearman ρ (baseline activation → response R)")
axA.set_title("a  Does a cell's activation state modulate the knockdown effect",loc="left")
# B: low/mid/high tercile response, mean over axis groups
axB=fig3.add_subplot(gs3[0,1])
for a,dfg in bs.groupby("axis"):
    axB.scatter(dfg["R_lowAct"],dfg["R_highAct"],color=COL[a],label=a,s=28)
lim=[min(bs[["R_lowAct","R_highAct"]].min())-.05,max(bs[["R_lowAct","R_highAct"]].max())+.05]
axB.plot(lim,lim,color="0.6",ls=":",lw=.8)
axB.set_xlim(lim); axB.set_ylim(lim)
axB.set_xlabel("response R in LOW-activation KD cells")
axB.set_ylabel("response R in HIGH-activation KD cells")
axB.set_title("b  Above the diagonal = effect larger in more-activated cells",loc="left")
axB.legend(frameon=False,fontsize=6.5)
# label only the notable points: largest |delta_hi_lo| plus the flagship extremes, to avoid pile-up
bs_lab=bs.reindex(bs["delta_hi_lo"].abs().sort_values(ascending=False).index)
label_set=set(bs_lab.head(6)["gene"]) | {"ZAP70","CD3E","PLCG1","CD247","CHD4"}
for _,rw in bs.iterrows():
    if rw["gene"] in label_set:
        axB.annotate(rw["gene"],(rw["R_lowAct"],rw["R_highAct"]),fontsize=6,style="italic",
                     xytext=(3,3),textcoords="offset points")
fig3.suptitle(f"Phase G · Fig 3 — Baseline-state sensitivity of the knockdown effect ({COND}, 4 donors)",
              x=.01,ha="left",fontweight="bold",fontsize=9.5)
fig3.savefig(f"{OUT}/G3_baseline_{COND}.png",bbox_inches="tight",dpi=150)
log("WROTE G3")
open(f"{OUT}/_DONE_{COND}","w").write("done\n")
log("ALL DONE",COND)
