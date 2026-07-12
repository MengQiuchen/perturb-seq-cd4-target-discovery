"""
Request 1 — WITHIN-CONDITION clustering & UMAP.

Unlike the Phase-D waves (which embed all conditions on ONE joint manifold, where the
dominant axis is Rest<->Stim activation and within-state substructure is squashed), this
embeds and Leiden-clusters EACH condition SEPARATELY. On a single-condition manifold the
activation axis is (mostly) removed, so finer substructure — naive/central-memory vs
effector, proliferating cells, cytotoxic/Treg-leaning subsets — can resolve.

Input: the per-shard stratified subsample caches already built for the joint waves
  viz_outputs/checkpoints/shardcache/D4_{cond}__nt25000_pt400_bg90000_s0_nom301.sub.h5ad
Each cache = one condition's subsample (~146k cells), raw counts in X, gene_name in var,
obs = guide_type / perturbed_gene_name / condition / is_nominated / QC.

Per condition it writes:
  checkpoints/within/{cond}.embedded.h5ad         (X_umap, leiden, lognorm layer)
  checkpoints/within/{cond}.obs_umap.csv.gz       (obs + UMAP coords + leiden)
  figures/within/W1_leiden_{cond}.png             UMAP colored by Leiden cluster (+labels)
  figures/within/W2_guidetype_{cond}.png          UMAP: targeting vs non-targeting control
  figures/within/W3_markers_{cond}.png            T-cell state marker genes (log-norm)
  figures/within/W4_composition_{cond}.png        cluster x guide_type composition bars
And a combined:
  figures/within/W0_leiden_all_conditions.png     3-panel side-by-side of the 3 manifolds
"""
import os, sys, time, json
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import scipy.sparse as sp
import anndata as ad

sys.path.insert(0, "/Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/code")
from embed_common import log, embed  # reuse the exact joint-wave embedding recipe

PROJ = "/Users/meng01/qiuchen/project/hackathon/perturb-seq"
CACHE = f"{PROJ}/viz_outputs/checkpoints/shardcache"
OUT = f"{PROJ}/phaseE_outputs/umap_within"
CKDIR = f"{OUT}/checkpoints/within"; FIG = f"{OUT}/figures/within"
os.makedirs(CKDIR, exist_ok=True); os.makedirs(FIG, exist_ok=True)

CONDS = ["Rest", "Stim8hr", "Stim48hr"]
COND_COLORS = {"Rest":"#0072B2","Stim8hr":"#E69F00","Stim48hr":"#009E73"}
GREY = "#D9D9D9"; HILITE = "#D55E00"
CACHE_TMPL = CACHE + "/D4_{cond}__nt25000_pt400_bg90000_s0_nom301.sub.h5ad"

# T-cell state marker panel (naive/memory, effector, activation, proliferation, Treg, lineage)
MARKERS = ["SELL","CCR7","TCF7","IL7R","GZMB","PRF1","IFNG","TNF",
           "IL2RA","CD69","TNFRSF9","MKI67","TOP2A","FOXP3","CTLA4","TBX21"]

# publication style (matches plot_umaps.py)
STYLE = {"figure.dpi":200,"savefig.dpi":300,"font.size":8.0,"axes.titlesize":8.0,
    "axes.labelsize":8.0,"xtick.labelsize":6.0,"ytick.labelsize":6.0,"legend.fontsize":7.0,
    "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":0.6,
    "axes.titlelocation":"left","font.family":"sans-serif","legend.frameon":False,
    "figure.facecolor":"white","savefig.bbox":"tight"}
mpl.rcParams.update(STYLE)
# 20-color CVD-aware-ish qualitative for leiden (tab20)
LEIDEN_CMAP = plt.get_cmap("tab20")

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
    if n>200_000: return 0.8
    if n>100_000: return 1.4
    if n>30_000:  return 3.0
    return 6.0

def resolve_gene_rows(adata, names):
    out=[]; gn = adata.var["gene_name"].values if "gene_name" in adata.var else None
    vn = adata.var_names.values
    for nm in names:
        gi=None
        if gn is not None:
            h=np.where(gn==nm)[0]
            if len(h): gi=int(h[0])
        if gi is None:
            h=np.where(vn==nm)[0]
            if len(h): gi=int(h[0])
        if gi is not None: out.append((nm,gi))
    return out

def plot_leiden(xy, leiden, cond, path, title=None):
    n=xy.shape[0]; s=_psize(n)
    cats=sorted(pd.unique(leiden), key=lambda x:int(x))
    fig, ax = plt.subplots(figsize=(4.0,3.8))
    for i,cl in enumerate(cats):
        m = leiden==cl
        ax.scatter(xy[m,0],xy[m,1], s=s, c=[LEIDEN_CMAP(i%20)], alpha=0.6, linewidths=0, rasterized=True)
    # cluster number labels at medians
    for cl in cats:
        m=leiden==cl; cx,cy=np.median(xy[m,0]),np.median(xy[m,1])
        ax.text(cx,cy,str(cl),fontsize=6.5,ha="center",va="center",fontweight="bold",
                color="#111",zorder=6,bbox=dict(boxstyle="round,pad=0.1",fc="white",ec="none",alpha=0.7))
    _clean(ax); _arrows(ax)
    ax.set_title(title or f"{cond} — within-condition Leiden clusters (n={n:,}, {len(cats)} clusters)", fontsize=8)
    fig.tight_layout(); fig.savefig(path); plt.close(fig); log(f"  wrote {path}")

def plot_guidetype(xy, gtype, cond, path):
    n=xy.shape[0]; s=_psize(n)
    fig, ax = plt.subplots(figsize=(3.8,3.6))
    tg = gtype=="targeting"; nt = gtype=="non-targeting"
    ax.scatter(xy[tg,0],xy[tg,1], s=s, c="#B9C7D6", alpha=0.35, linewidths=0, rasterized=True, zorder=1)
    ax.scatter(xy[nt,0],xy[nt,1], s=max(s,2.0), c=HILITE, alpha=0.7, linewidths=0, rasterized=True, zorder=3)
    _clean(ax); _arrows(ax)
    handles=[Line2D([0],[0],marker='o',ls='',mfc="#B9C7D6",mec='none',ms=5,label=f"targeting (n={tg.sum():,})"),
             Line2D([0],[0],marker='o',ls='',mfc=HILITE,mec='none',ms=5,label=f"non-targeting ctrl (n={nt.sum():,})")]
    ax.legend(handles=handles, loc="upper right", fontsize=6.5, handletextpad=0.2)
    ax.set_title(f"{cond} — perturbed vs non-targeting control", fontsize=8)
    fig.tight_layout(); fig.savefig(path); plt.close(fig); log(f"  wrote {path}")

def plot_markers(xy, lognorm, markers, cond, path):
    if not markers: log("  no markers present; skip"); return
    n=xy.shape[0]; s=_psize(n)
    ncol=4; nrow=int(np.ceil(len(markers)/ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.3*ncol, 2.3*nrow), squeeze=False)
    for i,ax in enumerate(axes.flat):
        if i>=len(markers): ax.axis("off"); continue
        g,gi = markers[i]
        col = lognorm[:,gi]
        vals = np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()
        order=np.argsort(vals)
        sc=ax.scatter(xy[order,0],xy[order,1], c=vals[order], s=s, cmap="viridis", alpha=0.8,
                      linewidths=0, rasterized=True, vmin=0, vmax=max(np.percentile(vals,99),1e-3))
        _clean(ax); ax.set_title(g, fontsize=8, style="italic")
        cb=fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02); cb.ax.tick_params(labelsize=5); cb.outline.set_visible(False)
    _arrows(axes[0,0])
    fig.suptitle(f"{cond} — T-cell state markers on within-condition manifold (log-norm)", fontsize=9, x=0.02, ha="left")
    fig.tight_layout(rect=[0,0,1,0.96]); fig.savefig(path); plt.close(fig); log(f"  wrote {path}")

def plot_composition(leiden, gtype, is_nom, cond, path):
    cats=sorted(pd.unique(leiden), key=lambda x:int(x))
    df = pd.DataFrame({"leiden":leiden,"gtype":gtype,"is_nom":is_nom})
    # fraction NTC per cluster and cluster sizes
    frac_nt = df.groupby("leiden")["gtype"].apply(lambda s:(s=="non-targeting").mean()).reindex(cats)
    sizes = df["leiden"].value_counts().reindex(cats)
    fig, (a1,a2) = plt.subplots(1,2, figsize=(7.0,3.0))
    a1.bar(range(len(cats)), sizes.values, color=COND_COLORS[cond], alpha=0.85)
    a1.set_xticks(range(len(cats))); a1.set_xticklabels(cats, fontsize=5.5)
    a1.set_xlabel("Leiden cluster"); a1.set_ylabel("cells"); a1.set_title(f"{cond} — cluster sizes", fontsize=8)
    a2.bar(range(len(cats)), 100*frac_nt.values, color=HILITE, alpha=0.85)
    a2.axhline(100*(gtype=="non-targeting").mean(), ls="--", lw=0.8, color="#444",
               label=f"overall {100*(gtype=='non-targeting').mean():.1f}%")
    a2.set_xticks(range(len(cats))); a2.set_xticklabels(cats, fontsize=5.5)
    a2.set_xlabel("Leiden cluster"); a2.set_ylabel("% non-targeting ctrl"); a2.set_title(f"{cond} — NTC fraction per cluster", fontsize=8)
    a2.legend(fontsize=6)
    fig.tight_layout(); fig.savefig(path); plt.close(fig); log(f"  wrote {path}")

def main():
    t0=time.time()
    umap_store={}   # cond -> (xy, leiden) for the combined panel
    summary={}
    for cond in CONDS:
        ck = f"{CKDIR}/{cond}.embedded.h5ad"
        if os.path.exists(ck):
            log(f"[{cond}] embedded checkpoint exists -> loading")
            A = ad.read_h5ad(ck)
        else:
            fn = CACHE_TMPL.format(cond=cond)
            assert os.path.exists(fn), f"missing cache {fn}"
            log(f"[{cond}] loading cache {fn}")
            A = ad.read_h5ad(fn)
            # keep only this condition's cells (cache is single-condition, but be safe)
            if "condition" in A.obs:
                A = A[A.obs["condition"]==cond].copy()
            log(f"[{cond}] {A.n_obs:,} cells x {A.n_vars:,} genes; embedding separately")
            A = embed(A, n_hvg=2000, n_pcs=50, n_neighbors=15, leiden_res=1.0,
                      integrate_donor=False, seed=0)
            A.write_h5ad(ck); log(f"[{cond}] wrote {ck}")
        xy = A.obsm["X_umap"]; leiden = A.obs["leiden"].astype(str).values
        gtype = A.obs["guide_type"].astype(str).values
        is_nom = A.obs["is_nominated"].values if "is_nominated" in A.obs else np.zeros(A.n_obs,bool)
        # obs+umap table
        odf = A.obs.copy(); odf["UMAP1"]=xy[:,0]; odf["UMAP2"]=xy[:,1]; odf["leiden"]=leiden
        odf.to_csv(f"{CKDIR}/{cond}.obs_umap.csv.gz", compression="gzip")
        # figures
        plot_leiden(xy, leiden, cond, f"{FIG}/W1_leiden_{cond}.png")
        plot_guidetype(xy, gtype, cond, f"{FIG}/W2_guidetype_{cond}.png")
        lognorm = A.layers["lognorm"] if "lognorm" in A.layers else A.X
        markers = resolve_gene_rows(A, MARKERS)
        plot_markers(xy, lognorm, markers, cond, f"{FIG}/W3_markers_{cond}.png")
        plot_composition(leiden, gtype, is_nom, cond, f"{FIG}/W4_composition_{cond}.png")
        umap_store[cond]=(xy, leiden)
        summary[cond]={"n_cells":int(A.n_obs),"n_clusters":int(len(set(leiden))),
                       "n_targeting":int((gtype=="targeting").sum()),
                       "n_nontargeting":int((gtype=="non-targeting").sum()),
                       "markers_present":[m for m,_ in markers]}
        del A
        log(f"[{cond}] done ({time.time()-t0:.0f}s elapsed)")

    # combined 3-panel leiden
    fig, axes = plt.subplots(1, len(CONDS), figsize=(3.4*len(CONDS), 3.4), squeeze=False); axes=axes[0]
    for ax,cond in zip(axes, CONDS):
        xy,leiden = umap_store[cond]; n=xy.shape[0]; s=_psize(n)
        cats=sorted(pd.unique(leiden), key=lambda x:int(x))
        for i,cl in enumerate(cats):
            m=leiden==cl
            ax.scatter(xy[m,0],xy[m,1], s=s, c=[LEIDEN_CMAP(i%20)], alpha=0.6, linewidths=0, rasterized=True)
        for cl in cats:
            m=leiden==cl; cx,cy=np.median(xy[m,0]),np.median(xy[m,1])
            ax.text(cx,cy,str(cl),fontsize=5.5,ha="center",va="center",fontweight="bold",color="#111",zorder=6,
                    bbox=dict(boxstyle="round,pad=0.1",fc="white",ec="none",alpha=0.7))
        _clean(ax); ax.set_title(f"{cond} (n={n:,}, {len(cats)} clusters)", fontsize=8)
    _arrows(axes[0])
    fig.suptitle("Within-condition Leiden clustering — each condition embedded separately (D4)", fontsize=9, x=0.02, ha="left")
    fig.tight_layout(rect=[0,0,1,0.95]); p=f"{FIG}/W0_leiden_all_conditions.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")

    json.dump(summary, open(f"{OUT}/within_condition_summary.json","w"), indent=2)
    log(f"ALL DONE in {time.time()-t0:.0f}s")

if __name__ == "__main__":
    main()
