"""
Request 2 (CPU SUBSET preview) — per-target UMAP from the FAST cached subsamples.

The raw-shard gather (all cells per target) is blocked by a degraded Lustre "slow"
OST (~1.4 MB/s cold). Meanwhile the 3 shardcache subsamples read at ~72 MB/s and
already contain every chosen target (nominated: 400/condition; non-nominated:
80-291 total). This script builds the per-target, condition-colored UMAPs from that
subset so we can see the distribution shape now; the full uncapped version runs on a
cuda batch node.

Reads:  viz_outputs/checkpoints/shardcache/D4_{cond}__...sub.h5ad  (raw counts, CSR)
Writes: phaseE_outputs/umap_within/
  checkpoints/pertarget_subset/{gene}.obs_umap.csv.gz
  figures/pertarget_subset/PTS_{gene}_bycondition.png
  figures/pertarget_subset/PTS_grid_bycondition.png
  per_target_subset_summary.json
"""
import os, sys, time, json, glob
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import scipy.sparse as sp
import anndata as ad

sys.path.insert(0, "/Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/code")
from embed_common import log

PROJ = "/Users/meng01/qiuchen/project/hackathon/perturb-seq"
CACHE = f"{PROJ}/viz_outputs/checkpoints/shardcache"
OUT = f"{PROJ}/phaseE_outputs/umap_within"
CKDIR = f"{OUT}/checkpoints/pertarget_subset"; FIG = f"{OUT}/figures/pertarget_subset"
os.makedirs(CKDIR, exist_ok=True); os.makedirs(FIG, exist_ok=True)
SCORECARD = f"{PROJ}/phaseBC_outputs/TARGET_SCORECARD.csv"
COND_COLORS = {"Rest":"#0072B2","Stim8hr":"#E69F00","Stim48hr":"#009E73"}

STYLE = {"figure.dpi":200,"savefig.dpi":300,"font.size":8.0,"axes.titlesize":8.0,
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

def _psize(n):
    if n>3000: return 5.0
    if n>1200: return 8.0
    if n>400:  return 12.0
    return 18.0

def load_targets_from_caches(targets):
    """Read only the chosen targets' cells from the 3 cache subsamples (raw counts)."""
    tset=set(targets); parts=[]
    for fn in sorted(glob.glob(f"{CACHE}/D4_*.sub.h5ad")):
        cond = fn.split("/")[-1].split("__")[0].replace("D4_","")
        log(f"[{cond}] loading cache obs")
        A = ad.read_h5ad(fn)   # fast pool ~72 MB/s
        pg = A.obs["perturbed_gene_name"].astype(str).values
        keep = np.isin(pg, list(tset))
        sub = A[keep].copy()
        sub.obs["condition"] = cond
        log(f"[{cond}] kept {sub.n_obs:,} target cells")
        parts.append(sub); del A
    U = ad.concat(parts, join="outer", index_unique="-")
    log(f"union subset {U.n_obs:,} cells x {U.n_vars:,} genes")
    return U

def embed_one(sub, seed=0):
    import scanpy as sc
    navail=len(os.sched_getaffinity(0)) if hasattr(os,"sched_getaffinity") else 8
    njobs=int(os.environ.get("EMBED_NJOBS", str(min(8,max(1,navail-1)))))
    sc.settings.n_jobs=njobs
    sc.pp.normalize_total(sub, target_sum=1e4); sc.pp.log1p(sub)
    sub.layers["lognorm"]=sub.X.copy()
    n=sub.n_obs
    nhvg=min(2000, max(200, sub.n_vars//5))
    sc.pp.highly_variable_genes(sub, n_top_genes=nhvg, flavor="seurat", subset=False)
    npcs=max(2, min(50, n-1, int(sub.var["highly_variable"].sum())))
    sc.pp.pca(sub, n_comps=npcs, mask_var="highly_variable", zero_center=True, random_state=seed)
    nnb=min(15, max(3, n//10))
    sc.pp.neighbors(sub, n_neighbors=nnb, n_pcs=npcs, use_rep="X_pca", random_state=seed)
    sc.tl.umap(sub, random_state=seed)
    return sub

def plot_one(sub, gene, path, rank, nominated):
    xy=sub.obsm["X_umap"]; cond=sub.obs["condition"].astype(str).values; n=xy.shape[0]; s=_psize(n)
    conds=[c for c in ["Rest","Stim8hr","Stim48hr"] if c in set(cond)]
    fig, ax = plt.subplots(figsize=(3.8,3.6))
    rng=np.random.default_rng(0); order=rng.permutation(n)
    colors=np.array([COND_COLORS[c] for c in cond],dtype=object)
    ax.scatter(xy[order,0],xy[order,1], s=s, c=list(colors[order]), alpha=0.8, linewidths=0, rasterized=True)
    _clean(ax); _arrows(ax)
    handles=[Line2D([0],[0],marker='o',ls='',mfc=COND_COLORS[c],mec='none',ms=5,
             label=f"{c} (n={(cond==c).sum():,})") for c in conds]
    ax.legend(handles=handles, loc="upper right", fontsize=6.5, handletextpad=0.2)
    tag = f"rank {int(rank)}" if (rank==rank and rank is not None) else "not nominated"
    ax.set_title(f"{gene} — per-target UMAP by condition (n={n:,} · {tag})", fontsize=8, style="italic")
    fig.tight_layout(); fig.savefig(path); plt.close(fig); log(f"  wrote {path}")

def main():
    t0=time.time()
    cfg=json.load(open(f"{PROJ}/phaseE_scripts/pt_targets.json"))
    targets=cfg["targets"]
    log(f"targets: {targets}")
    U = load_targets_from_caches(targets)
    sc_df=pd.read_csv(SCORECARD); sc_df["gene"]=sc_df["gene"].astype(str)
    rankmap=dict(zip(sc_df["gene"], sc_df["rank"]))
    nomset=set(sc_df.loc[sc_df["rank"].notna(),"gene"])
    summary={}; embedded={}
    for gene in targets:
        sub = U[U.obs["perturbed_gene_name"].astype(str)==gene].copy()
        if sub.n_obs < 30:
            log(f"[{gene}] only {sub.n_obs} cells — skipping"); continue
        log(f"[{gene}] embedding {sub.n_obs:,} cells")
        sub = embed_one(sub)
        rk = rankmap.get(gene, np.nan)
        odf=sub.obs.copy(); odf["UMAP1"]=sub.obsm["X_umap"][:,0]; odf["UMAP2"]=sub.obsm["X_umap"][:,1]
        odf.to_csv(f"{CKDIR}/{gene}.obs_umap.csv.gz", compression="gzip")
        plot_one(sub, gene, f"{FIG}/PTS_{gene}_bycondition.png", rk, gene in nomset)
        cc=sub.obs["condition"].value_counts().to_dict()
        summary[gene]={"n_cells":int(sub.n_obs),"rank":(None if pd.isna(rk) else int(rk)),
                       "nominated": gene in nomset,
                       "by_condition":{k:int(v) for k,v in cc.items()}}
        embedded[gene]=(sub.obsm["X_umap"], sub.obs["condition"].astype(str).values, rk)
        del sub
        log(f"[{gene}] done ({time.time()-t0:.0f}s)")

    genes=[g for g in targets if g in embedded]
    if genes:
        ncol=min(4,len(genes)); nrow=int(np.ceil(len(genes)/ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(2.6*ncol, 2.6*nrow), squeeze=False)
        for i,ax in enumerate(axes.flat):
            if i>=len(genes): ax.axis("off"); continue
            g=genes[i]; xy,cond,rk=embedded[g]; n=xy.shape[0]; s=_psize(n)
            rng=np.random.default_rng(0); order=rng.permutation(n)
            colors=np.array([COND_COLORS[c] for c in cond],dtype=object)
            ax.scatter(xy[order,0],xy[order,1], s=s, c=list(colors[order]), alpha=0.8, linewidths=0, rasterized=True)
            _clean(ax)
            rtxt="" if pd.isna(rk) else f" · r{int(rk)}"
            ax.set_title(f"{g} (n={n:,}{rtxt})", fontsize=7.5, style="italic")
        _arrows(axes[0,0])
        handles=[Line2D([0],[0],marker='o',ls='',mfc=COND_COLORS[c],mec='none',ms=5,label=c) for c in ["Rest","Stim8hr","Stim48hr"]]
        fig.legend(handles=handles, loc="upper right", ncol=3, fontsize=7, bbox_to_anchor=(0.99,1.0))
        fig.suptitle("Per-target UMAPs (CPU subset from caches) colored by condition — D4", fontsize=9, x=0.02, ha="left")
        fig.tight_layout(rect=[0,0,1,0.95]); p=f"{FIG}/PTS_grid_bycondition.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")
    json.dump(summary, open(f"{OUT}/per_target_subset_summary.json","w"), indent=2)
    log(f"ALL DONE in {time.time()-t0:.0f}s; {len(genes)} targets plotted")

if __name__=="__main__":
    main()
