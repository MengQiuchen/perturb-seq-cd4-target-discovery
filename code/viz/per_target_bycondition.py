"""
Per-target UMAP colored by CONDITION (not single highlight hue), + exact per-condition
cell-count crosstab for the nominated targets. Answers: for each knockdown target, how many
cells come from Rest vs Stim8hr vs Stim48hr, and where do each land on the manifold.

Usage: python per_target_bycondition.py <embedded_checkpoint.h5ad> <tag>
Writes:
  figures/D2_per_target_bycondition_<tag>.png   (top-25 targets, dots colored by condition)
  checkpoints/<tag>.per_target_condition_counts.csv  (target x condition count table)
"""
import os, sys, time, json
import numpy as np, pandas as pd
import anndata as ad
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

PROJ="/Users/meng01/qiuchen/project/hackathon/perturb-seq"
OUT=f"{PROJ}/viz_outputs"; FIG=f"{OUT}/figures"; CK=f"{OUT}/checkpoints"
SCORECARD=f"{PROJ}/phaseBC_outputs/TARGET_SCORECARD.csv"
COND_COLORS={"Rest":"#0072B2","Stim8hr":"#E69F00","Stim48hr":"#009E73"}
GREY="#D9D9D9"

STYLE={"figure.dpi":200,"savefig.dpi":300,"font.size":8.0,"axes.titlesize":8.0,
       "axes.labelsize":8.0,"xtick.labelsize":6.0,"ytick.labelsize":6.0,"legend.fontsize":7.0,
       "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":0.6,
       "axes.titlelocation":"left","font.family":"sans-serif","legend.frameon":False,
       "figure.facecolor":"white","savefig.bbox":"tight"}
mpl.rcParams.update(STYLE)

def _clean(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xlabel(""); ax.set_ylabel("")

def _arrows(ax):
    x0,y0,L=0.02,0.02,0.13
    ax.annotate("",xy=(x0+L,y0),xytext=(x0,y0),xycoords="axes fraction",arrowprops=dict(arrowstyle="-|>",color="#444",lw=0.8))
    ax.annotate("",xy=(x0,y0+L),xytext=(x0,y0),xycoords="axes fraction",arrowprops=dict(arrowstyle="-|>",color="#444",lw=0.8))
    ax.text(x0+L+0.01,y0,"UMAP1",transform=ax.transAxes,fontsize=5.5,ha="left",va="center",color="#444")
    ax.text(x0,y0+L+0.01,"UMAP2",transform=ax.transAxes,fontsize=5.5,ha="center",va="bottom",color="#444",rotation=90)

def main():
    ckpt, tag = sys.argv[1], sys.argv[2]
    log(f"loading {ckpt}")
    adata=ad.read_h5ad(ckpt)
    xy=adata.obsm["X_umap"]; obs=adata.obs
    conds=[c for c in ["Rest","Stim8hr","Stim48hr"] if c in set(obs["condition"])]
    n=xy.shape[0]; s=1.2 if n>100_000 else 2.5

    sc_df=pd.read_csv(SCORECARD)
    present=set(obs["perturbed_gene_name"].unique())
    top=[g for g in sc_df.sort_values("rank")["gene"].astype(str) if g in present][:25]

    # ---- exact per-condition counts for ALL present nominated targets ----
    nom=set(sc_df["gene"].astype(str))
    sub=obs[obs["perturbed_gene_name"].isin(nom)]
    ct=pd.crosstab(sub["perturbed_gene_name"], sub["condition"])
    ct["total"]=ct.sum(axis=1)
    # order by scorecard rank
    rank_map=dict(zip(sc_df["gene"].astype(str), sc_df["rank"]))
    ct["rank"]=[rank_map.get(g,9999) for g in ct.index]
    ct=ct.sort_values("rank")
    ct.to_csv(f"{CK}/{tag}.per_target_condition_counts.csv")
    log(f"wrote per-target x condition counts ({ct.shape[0]} targets)")

    # ---- state-colored per-target grid (top 25) ----
    ncol=5; nrow=int(np.ceil(len(top)/ncol))
    fig,axes=plt.subplots(nrow,ncol,figsize=(2.05*ncol,2.05*nrow),squeeze=False)
    for i,ax in enumerate(axes.flat):
        if i>=len(top): ax.axis("off"); continue
        t=top[i]; m=(obs["perturbed_gene_name"]==t).values
        ax.scatter(xy[:,0],xy[:,1],s=s,c=GREY,alpha=0.25,linewidths=0,rasterized=True,zorder=1)
        for cnd in conds:
            mm=m & (obs["condition"]==cnd).values
            if mm.any():
                ax.scatter(xy[mm,0],xy[mm,1],s=max(s*3,4),c=COND_COLORS[cnd],alpha=0.85,linewidths=0,rasterized=True,zorder=3)
        _clean(ax)
        parts=" ".join(f"{cnd[:1]}{int((m&(obs['condition']==cnd).values).sum())}" for cnd in conds)
        ax.set_title(f"{t}  (n={m.sum():,})", fontsize=7, style="italic")
    _arrows(axes[0,0])
    handles=[Line2D([0],[0],marker='o',ls='',mfc=COND_COLORS[c],mec='none',ms=5,label=c) for c in conds]
    fig.legend(handles=handles, loc="upper right", ncol=len(conds), fontsize=7, bbox_to_anchor=(0.99,1.0))
    fig.suptitle(f"Per-target cells colored by condition (Rest/Stim) — {tag}", fontsize=9, x=0.02, ha="left")
    fig.tight_layout(rect=[0,0,1,0.96])
    p=f"{FIG}/D2_per_target_bycondition_{tag}.png"; fig.savefig(p); plt.close(fig)
    log(f"wrote {p}")

if __name__=="__main__":
    main()
