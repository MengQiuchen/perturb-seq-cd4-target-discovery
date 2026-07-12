"""
Cross-donor replication (D1-D4, Stim8hr) of the WITHIN-CONDITION phenomenon
from the UMAP report: Leiden clusters correspond to T-cell states, and
non-targeting controls (NTC) spread across all clusters -> manifold structure
reflects cell STATE, not perturbation identity.

Uses the matched stratified Stim8hr subsets (40k NTC + 40k background +
~2.7k flagship-target cells) that exist for all 4 donors. Identical embed
pipeline per donor: normalize_total(1e4)->log1p->2000 HVG->PCA(50)->
neighbors(15)->UMAP->Leiden(res=1.0).

Reads:  phaseE_outputs/checkpoints/{D}_Stim8hr.subset.h5ad   (raw CSR counts)
Writes: phaseE_outputs/umap_within/crossdonor/
  checkpoints/{D}.obs_umap.csv.gz
  figures/CD1_leiden_{D}.png, CD2_guidetype_{D}.png, CD3_markers_{D}.png
  figures/CD0_leiden_grid.png, CD0_guidetype_grid.png, CD_ntc_fraction.png
  crossdonor_summary.json
"""
import os, sys, time, json
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import anndata as ad

sys.path.insert(0, "/Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/code")
from embed_common import log, embed

PROJ = "/Users/meng01/qiuchen/project/hackathon/perturb-seq"
CKIN = f"{PROJ}/phaseE_outputs/checkpoints"
OUT = f"{PROJ}/phaseE_outputs/umap_within/crossdonor"
CKD = f"{OUT}/checkpoints"; FIG = f"{OUT}/figures"
os.makedirs(CKD, exist_ok=True); os.makedirs(FIG, exist_ok=True)

DONORS = ["D1","D2","D3","D4"]
COND = "Stim8hr"
MARKERS = ["SELL","CCR7","TCF7","IL7R","GZMB","PRF1","IFNG","TNF",
           "IL2RA","CD69","TNFRSF9","MKI67","TOP2A","FOXP3","CTLA4","TBX21"]
GUIDE_COLORS = {"non-targeting":"#D55E00","targeting":"#999999"}  # NTC highlighted, targeting grey
DONOR_COLORS = {"D1":"#0072B2","D2":"#E69F00","D3":"#009E73","D4":"#CC79A7"}

STYLE = {"figure.dpi":200,"savefig.dpi":300,"font.size":8.0,"axes.titlesize":8.0,
    "axes.labelsize":8.0,"xtick.labelsize":6.0,"ytick.labelsize":6.0,"legend.fontsize":7.0,
    "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":0.6,
    "axes.titlelocation":"left","font.family":"sans-serif","legend.frameon":False,
    "figure.facecolor":"white","savefig.bbox":"tight"}
mpl.rcParams.update(STYLE)
try:
    import matplotlib.cm as cm
    TAB20 = [mpl.colors.to_hex(c) for c in cm.get_cmap("tab20").colors]
except Exception:
    TAB20 = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2",
             "#7f7f7f","#bcbd22","#17becf"]*2

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

def plot_leiden(sub, D, path):
    xy=sub.obsm["X_umap"]; lab=sub.obs["leiden"].astype(str).values
    cats=sorted(set(lab), key=lambda x:int(x)); n=xy.shape[0]
    fig, ax = plt.subplots(figsize=(4.2,4.0))
    rng=np.random.default_rng(0); order=rng.permutation(n)
    cmap={c:TAB20[i%len(TAB20)] for i,c in enumerate(cats)}
    colors=np.array([cmap[l] for l in lab],dtype=object)
    ax.scatter(xy[order,0],xy[order,1], s=3.0, c=list(colors[order]), alpha=0.75, linewidths=0, rasterized=True)
    # cluster labels at medians
    for c in cats:
        m=lab==c; mx,my=np.median(xy[m,0]),np.median(xy[m,1])
        ax.text(mx,my,c,fontsize=6,fontweight="bold",ha="center",va="center",
                bbox=dict(boxstyle="round,pad=0.1",fc="white",ec="none",alpha=0.6))
    _clean(ax); _arrows(ax)
    ax.set_title(f"{D} {COND} — Leiden clusters (n={n:,}, {len(cats)} clusters)", fontsize=8.5)
    fig.tight_layout(); fig.savefig(path); plt.close(fig); log(f"  wrote {path}")

def plot_guidetype(sub, D, path):
    xy=sub.obsm["X_umap"]; gt=sub.obs["guide_type"].astype(str).values; n=xy.shape[0]
    fig, ax = plt.subplots(figsize=(4.2,4.0))
    # draw targeting first (grey), NTC on top (orange) so spread is visible
    for g in ["targeting","non-targeting"]:
        m=gt==g
        ax.scatter(xy[m,0],xy[m,1], s=3.0, c=GUIDE_COLORS[g], alpha=0.5 if g=="targeting" else 0.55,
                   linewidths=0, rasterized=True, label=f"{g} (n={m.sum():,})")
    _clean(ax); _arrows(ax)
    ax.legend(loc="upper right", fontsize=6.8, markerscale=2.0, handletextpad=0.2)
    ax.set_title(f"{D} {COND} — NTC vs targeting", fontsize=8.5)
    fig.tight_layout(); fig.savefig(path); plt.close(fig); log(f"  wrote {path}")

def plot_markers(sub, D, path):
    import scanpy as sc
    present=[m for m in MARKERS if m in set(sub.var_names)]
    xy=sub.obsm["X_umap"]; n=xy.shape[0]
    ncol=4; nrow=int(np.ceil(len(present)/ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.5*ncol, 2.4*nrow), squeeze=False)
    L = sub.layers["lognorm"]
    for i,ax in enumerate(axes.flat):
        if i>=len(present): ax.axis("off"); continue
        g=present[i]; gi=list(sub.var_names).index(g)
        vals=np.asarray(L[:,gi].todense()).ravel() if hasattr(L,"todense") else np.asarray(L[:,gi]).ravel()
        vmax=np.quantile(vals,0.99) if vals.max()>0 else 1.0
        order=np.argsort(vals)  # high on top
        scx=ax.scatter(xy[order,0],xy[order,1], s=2.0, c=vals[order], cmap="viridis",
                       vmin=0, vmax=max(vmax,1e-6), alpha=0.85, linewidths=0, rasterized=True)
        _clean(ax); ax.set_title(g, fontsize=7.5, style="italic")
    _arrows(axes[0,0])
    fig.suptitle(f"{D} {COND} — canonical T-cell markers (lognorm)", fontsize=9, x=0.02, ha="left")
    fig.tight_layout(rect=[0,0,1,0.97]); fig.savefig(path); plt.close(fig); log(f"  wrote {path}")

def main():
    t0=time.time()
    results={}; embedded={}
    for D in DONORS:
        fn=f"{CKIN}/{D}_{COND}.subset.h5ad"
        if not os.path.exists(fn):
            log(f"[{D}] subset MISSING {fn}"); continue
        log(f"[{D}] loading {fn}")
        A=ad.read_h5ad(fn)
        # These subsets store LOG-NORMALIZED values in .X and RAW COUNTS in
        # layers['counts']. embed() runs normalize_total+log1p, so feed it the
        # raw counts (exactly matching the D4 within-condition run).
        if "counts" in A.layers:
            A.X = A.layers["counts"].copy()
            del A.layers["counts"]
        xs=A.X[:200].data if hasattr(A.X[:200],"data") else np.asarray(A.X[:200]).ravel()
        is_int=bool(np.allclose(xs, np.round(xs)))
        log(f"[{D}] n={A.n_obs:,} nvar={A.n_vars:,} raw_counts={is_int} (from counts layer)")
        assert is_int, f"{D}: X not raw counts after swap"
        A.obs["donor"]=D
        log(f"[{D}] embedding")
        A=embed(A, integrate_donor=False, seed=0)
        # metrics: NTC spread across clusters
        lab=A.obs["leiden"].astype(str).values
        gt=A.obs["guide_type"].astype(str).values
        cats=sorted(set(lab), key=lambda x:int(x))
        glob_ntc=(gt=="non-targeting").mean()
        rows=[]
        for c in cats:
            m=lab==c; ntc=(gt[m]=="non-targeting").mean(); ncell=int(m.sum())
            rows.append({"cluster":c,"n_cells":ncell,"ntc_frac":float(ntc),
                         "log2_obs_exp":float(np.log2((ntc+1e-6)/(glob_ntc+1e-6)))})
        cl=pd.DataFrame(rows)
        n_with_ntc=int((cl["ntc_frac"]>0).sum())
        results[D]={"n_cells":int(A.n_obs),"n_clusters":len(cats),
                    "global_ntc_frac":float(glob_ntc),
                    "n_clusters_with_ntc":n_with_ntc,
                    "ntc_frac_min":float(cl["ntc_frac"].min()),
                    "ntc_frac_max":float(cl["ntc_frac"].max()),
                    "ntc_frac_cv":float(cl["ntc_frac"].std()/max(cl["ntc_frac"].mean(),1e-9)),
                    "max_abs_log2_obs_exp":float(cl["log2_obs_exp"].abs().max()),
                    "per_cluster":cl.to_dict(orient="records")}
        # save obs+umap
        odf=A.obs.copy(); odf["UMAP1"]=A.obsm["X_umap"][:,0]; odf["UMAP2"]=A.obsm["X_umap"][:,1]
        odf.to_csv(f"{CKD}/{D}.obs_umap.csv.gz", compression="gzip")
        # per-donor figures
        plot_leiden(A, D, f"{FIG}/CD1_leiden_{D}.png")
        plot_guidetype(A, D, f"{FIG}/CD2_guidetype_{D}.png")
        plot_markers(A, D, f"{FIG}/CD3_markers_{D}.png")
        embedded[D]=(A.obsm["X_umap"].copy(), lab.copy(), gt.copy())
        del A
        log(f"[{D}] done ({time.time()-t0:.0f}s); {len(cats)} clusters, NTC in {n_with_ntc}/{len(cats)}, global NTC {glob_ntc:.3f}")

    donors=[d for d in DONORS if d in embedded]
    # combined Leiden grid
    if donors:
        fig, axes = plt.subplots(1, len(donors), figsize=(3.4*len(donors), 3.4), squeeze=False)
        for i,D in enumerate(donors):
            ax=axes[0,i]; xy,lab,gt=embedded[D]; n=xy.shape[0]
            cats=sorted(set(lab), key=lambda x:int(x))
            cmap={c:TAB20[j%len(TAB20)] for j,c in enumerate(cats)}
            rng=np.random.default_rng(0); order=rng.permutation(n)
            colors=np.array([cmap[l] for l in lab],dtype=object)
            ax.scatter(xy[order,0],xy[order,1], s=2.0, c=list(colors[order]), alpha=0.75, linewidths=0, rasterized=True)
            _clean(ax); ax.set_title(f"{D} ({len(cats)} clusters)", fontsize=8.5)
        _arrows(axes[0,0])
        fig.suptitle(f"Within-condition Leiden clustering — {COND}, all 4 donors", fontsize=9.5, x=0.02, ha="left")
        fig.tight_layout(rect=[0,0,1,0.95]); p=f"{FIG}/CD0_leiden_grid.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")
        # combined guidetype grid
        fig, axes = plt.subplots(1, len(donors), figsize=(3.4*len(donors), 3.4), squeeze=False)
        for i,D in enumerate(donors):
            ax=axes[0,i]; xy,lab,gt=embedded[D]
            for g in ["targeting","non-targeting"]:
                m=gt==g
                ax.scatter(xy[m,0],xy[m,1], s=2.0, c=GUIDE_COLORS[g], alpha=0.5 if g=="targeting" else 0.55,
                           linewidths=0, rasterized=True)
            _clean(ax); ax.set_title(D, fontsize=8.5)
        _arrows(axes[0,0])
        handles=[Line2D([0],[0],marker='o',ls='',mfc=GUIDE_COLORS[g],mec='none',ms=5,label=g) for g in ["non-targeting","targeting"]]
        fig.legend(handles=handles, loc="upper right", ncol=2, fontsize=8, bbox_to_anchor=(0.99,1.0))
        fig.suptitle(f"NTC vs targeting — {COND}, all 4 donors (NTC spreads across all clusters)", fontsize=9.5, x=0.02, ha="left")
        fig.tight_layout(rect=[0,0,1,0.95]); p=f"{FIG}/CD0_guidetype_grid.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")
        # NTC fraction per cluster (replication metric)
        fig, ax = plt.subplots(figsize=(6.2,3.8))
        for D in donors:
            cl=pd.DataFrame(results[D]["per_cluster"]); cl["ci"]=cl["cluster"].astype(int)
            cl=cl.sort_values("ci")
            ax.plot(range(len(cl)), cl["ntc_frac"], marker="o", ms=3.5, lw=1.0,
                    color=DONOR_COLORS[D], label=f"{D} (global {results[D]['global_ntc_frac']:.2f})")
            ax.axhline(results[D]["global_ntc_frac"], color=DONOR_COLORS[D], ls=":", lw=0.7, alpha=0.6)
        ax.set_xlabel("Leiden cluster (sorted)"); ax.set_ylabel("NTC fraction in cluster")
        ax.set_title("NTC fraction per cluster — roughly flat = structure driven by cell state, not perturbation", fontsize=8)
        ax.legend(fontsize=7, loc="best"); ax.set_ylim(0,1)
        fig.tight_layout(); p=f"{FIG}/CD_ntc_fraction.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")

    json.dump(results, open(f"{OUT}/crossdonor_summary.json","w"), indent=2)
    log(f"ALL DONE in {time.time()-t0:.0f}s; donors: {donors}")

if __name__=="__main__":
    main()
