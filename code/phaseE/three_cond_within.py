"""
Three-condition x 4-donor WITHIN-CONDITION replication of the UMAP-report
phenomenon: Leiden clusters correspond to T-cell states, and non-targeting
controls (NTC) spread across all clusters -> manifold structure reflects
cell STATE, not perturbation identity.

Generalizes cross_donor_within.py from Stim8hr-only to ALL THREE conditions
(Rest / Stim8hr / Stim48hr) x 4 donors (D1-D4). Uses the matched stratified
subsets (40k NTC + 40k background + flagship-target cells) built by
phaseE_gather.py (seed=1) for every donor x condition. Identical embed
pipeline per subset: normalize_total(1e4)->log1p->2000 HVG->PCA(50)->
neighbors(15)->UMAP->Leiden(res=1.0).

Reads:  phaseE_outputs/checkpoints/{D}_{COND}.subset.h5ad  (X=lognorm, layers['counts']=raw)
Writes: phaseE_outputs/umap_within/threecond/
  checkpoints/{D}_{COND}.obs_umap.csv.gz
  figures/TC1_leiden_{COND}_{D}.png, TC2_guidetype_{COND}_{D}.png, TC3_markers_{COND}_{D}.png
  figures/TC0_leiden_grid_{COND}.png, TC0_guidetype_grid_{COND}.png, TC_ntc_fraction_{COND}.png
  figures/TC_ntc_fraction_all.png
  threecond_summary.json
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
OUT = f"{PROJ}/phaseE_outputs/umap_within/threecond"
CKD = f"{OUT}/checkpoints"; FIG = f"{OUT}/figures"
os.makedirs(CKD, exist_ok=True); os.makedirs(FIG, exist_ok=True)

DONORS = ["D1","D2","D3","D4"]
CONDITIONS = ["Rest","Stim8hr","Stim48hr"]
MARKERS = ["SELL","CCR7","TCF7","IL7R","GZMB","PRF1","IFNG","TNF",
           "IL2RA","CD69","TNFRSF9","MKI67","TOP2A","FOXP3","CTLA4","TBX21"]
GUIDE_COLORS = {"non-targeting":"#D55E00","targeting":"#999999"}
DONOR_COLORS = {"D1":"#0072B2","D2":"#E69F00","D3":"#009E73","D4":"#CC79A7"}
COND_COLORS = {"Rest":"#0072B2","Stim8hr":"#E69F00","Stim48hr":"#009E73"}

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

def plot_leiden(sub, D, COND, path):
    xy=sub.obsm["X_umap"]; lab=sub.obs["leiden"].astype(str).values
    cats=sorted(set(lab), key=lambda x:int(x)); n=xy.shape[0]
    fig, ax = plt.subplots(figsize=(4.2,4.0))
    rng=np.random.default_rng(0); order=rng.permutation(n)
    cmap={c:TAB20[i%len(TAB20)] for i,c in enumerate(cats)}
    colors=np.array([cmap[l] for l in lab],dtype=object)
    ax.scatter(xy[order,0],xy[order,1], s=3.0, c=list(colors[order]), alpha=0.75, linewidths=0, rasterized=True)
    for c in cats:
        m=lab==c; mx,my=np.median(xy[m,0]),np.median(xy[m,1])
        ax.text(mx,my,c,fontsize=6,fontweight="bold",ha="center",va="center",
                bbox=dict(boxstyle="round,pad=0.1",fc="white",ec="none",alpha=0.6))
    _clean(ax); _arrows(ax)
    ax.set_title(f"{D} {COND} — Leiden clusters (n={n:,}, {len(cats)} clusters)", fontsize=8.5)
    fig.tight_layout(); fig.savefig(path); plt.close(fig); log(f"  wrote {path}")

def plot_guidetype(sub, D, COND, path):
    xy=sub.obsm["X_umap"]; gt=sub.obs["guide_type"].astype(str).values; n=xy.shape[0]
    fig, ax = plt.subplots(figsize=(4.2,4.0))
    for g in ["targeting","non-targeting"]:
        m=gt==g
        ax.scatter(xy[m,0],xy[m,1], s=3.0, c=GUIDE_COLORS[g], alpha=0.5 if g=="targeting" else 0.55,
                   linewidths=0, rasterized=True, label=f"{g} (n={m.sum():,})")
    _clean(ax); _arrows(ax)
    ax.legend(loc="upper right", fontsize=6.8, markerscale=2.0, handletextpad=0.2)
    ax.set_title(f"{D} {COND} — NTC vs targeting", fontsize=8.5)
    fig.tight_layout(); fig.savefig(path); plt.close(fig); log(f"  wrote {path}")

def plot_markers(sub, D, COND, path):
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
        order=np.argsort(vals)
        ax.scatter(xy[order,0],xy[order,1], s=2.0, c=vals[order], cmap="viridis",
                   vmin=0, vmax=max(vmax,1e-6), alpha=0.85, linewidths=0, rasterized=True)
        _clean(ax); ax.set_title(g, fontsize=7.5, style="italic")
    _arrows(axes[0,0])
    fig.suptitle(f"{D} {COND} — canonical T-cell markers (lognorm)", fontsize=9, x=0.02, ha="left")
    fig.tight_layout(rect=[0,0,1,0.97]); fig.savefig(path); plt.close(fig); log(f"  wrote {path}")

def main():
    t0=time.time()
    results={}     # results[COND][D] = metrics
    embedded={}    # embedded[COND][D] = (xy, lab, gt)
    for COND in CONDITIONS:
        results[COND]={}; embedded[COND]={}
        for D in DONORS:
            fn=f"{CKIN}/{D}_{COND}.subset.h5ad"
            if not os.path.exists(fn):
                log(f"[{D} {COND}] subset MISSING {fn}"); continue
            log(f"[{D} {COND}] loading {fn}")
            A=ad.read_h5ad(fn)
            if "counts" in A.layers:
                A.X = A.layers["counts"].copy()
                del A.layers["counts"]
            xs=A.X[:200].data if hasattr(A.X[:200],"data") else np.asarray(A.X[:200]).ravel()
            is_int=bool(np.allclose(xs, np.round(xs)))
            log(f"[{D} {COND}] n={A.n_obs:,} nvar={A.n_vars:,} raw_counts={is_int} (from counts layer)")
            assert is_int, f"{D} {COND}: X not raw counts after swap"
            A.obs["donor"]=D
            log(f"[{D} {COND}] embedding")
            A=embed(A, integrate_donor=False, seed=0)
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
            results[COND][D]={"n_cells":int(A.n_obs),"n_clusters":len(cats),
                        "global_ntc_frac":float(glob_ntc),
                        "n_clusters_with_ntc":n_with_ntc,
                        "ntc_frac_min":float(cl["ntc_frac"].min()),
                        "ntc_frac_max":float(cl["ntc_frac"].max()),
                        "ntc_frac_cv":float(cl["ntc_frac"].std()/max(cl["ntc_frac"].mean(),1e-9)),
                        "max_abs_log2_obs_exp":float(cl["log2_obs_exp"].abs().max()),
                        "per_cluster":cl.to_dict(orient="records")}
            odf=A.obs.copy(); odf["UMAP1"]=A.obsm["X_umap"][:,0]; odf["UMAP2"]=A.obsm["X_umap"][:,1]
            odf.to_csv(f"{CKD}/{D}_{COND}.obs_umap.csv.gz", compression="gzip")
            plot_leiden(A, D, COND, f"{FIG}/TC1_leiden_{COND}_{D}.png")
            plot_guidetype(A, D, COND, f"{FIG}/TC2_guidetype_{COND}_{D}.png")
            plot_markers(A, D, COND, f"{FIG}/TC3_markers_{COND}_{D}.png")
            embedded[COND][D]=(A.obsm["X_umap"].copy(), lab.copy(), gt.copy())
            del A
            log(f"[{D} {COND}] done ({time.time()-t0:.0f}s); {len(cats)} clusters, NTC in {n_with_ntc}/{len(cats)}, global NTC {glob_ntc:.3f}")

        # ---- per-condition combined grids ----
        donors=[d for d in DONORS if d in embedded[COND]]
        if donors:
            # Leiden grid
            fig, axes = plt.subplots(1, len(donors), figsize=(3.4*len(donors), 3.4), squeeze=False)
            for i,D in enumerate(donors):
                ax=axes[0,i]; xy,lab,gt=embedded[COND][D]; n=xy.shape[0]
                cats=sorted(set(lab), key=lambda x:int(x))
                cmap={c:TAB20[j%len(TAB20)] for j,c in enumerate(cats)}
                rng=np.random.default_rng(0); order=rng.permutation(n)
                colors=np.array([cmap[l] for l in lab],dtype=object)
                ax.scatter(xy[order,0],xy[order,1], s=2.0, c=list(colors[order]), alpha=0.75, linewidths=0, rasterized=True)
                _clean(ax); ax.set_title(f"{D} ({len(cats)} clusters)", fontsize=8.5)
            _arrows(axes[0,0])
            fig.suptitle(f"Within-condition Leiden clustering — {COND}, all 4 donors", fontsize=9.5, x=0.02, ha="left")
            fig.tight_layout(rect=[0,0,1,0.95]); p=f"{FIG}/TC0_leiden_grid_{COND}.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")
            # guidetype grid
            fig, axes = plt.subplots(1, len(donors), figsize=(3.4*len(donors), 3.4), squeeze=False)
            for i,D in enumerate(donors):
                ax=axes[0,i]; xy,lab,gt=embedded[COND][D]
                for g in ["targeting","non-targeting"]:
                    m=gt==g
                    ax.scatter(xy[m,0],xy[m,1], s=2.0, c=GUIDE_COLORS[g], alpha=0.5 if g=="targeting" else 0.55,
                               linewidths=0, rasterized=True)
                _clean(ax); ax.set_title(D, fontsize=8.5)
            _arrows(axes[0,0])
            handles=[Line2D([0],[0],marker='o',ls='',mfc=GUIDE_COLORS[g],mec='none',ms=5,label=g) for g in ["non-targeting","targeting"]]
            fig.legend(handles=handles, loc="upper right", ncol=2, fontsize=8, bbox_to_anchor=(0.99,1.0))
            fig.suptitle(f"NTC vs targeting — {COND}, all 4 donors (NTC spreads across all clusters)", fontsize=9.5, x=0.02, ha="left")
            fig.tight_layout(rect=[0,0,1,0.95]); p=f"{FIG}/TC0_guidetype_grid_{COND}.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")
            # NTC fraction per cluster
            fig, ax = plt.subplots(figsize=(6.2,3.8))
            for D in donors:
                cl=pd.DataFrame(results[COND][D]["per_cluster"]); cl["ci"]=cl["cluster"].astype(int)
                cl=cl.sort_values("ci")
                ax.plot(range(len(cl)), cl["ntc_frac"], marker="o", ms=3.5, lw=1.0,
                        color=DONOR_COLORS[D], label=f"{D} (global {results[COND][D]['global_ntc_frac']:.2f})")
                ax.axhline(results[COND][D]["global_ntc_frac"], color=DONOR_COLORS[D], ls=":", lw=0.7, alpha=0.6)
            ax.set_xlabel("Leiden cluster (sorted)"); ax.set_ylabel("NTC fraction in cluster")
            ax.set_title(f"NTC fraction per cluster — {COND} (flat = structure driven by state, not perturbation)", fontsize=8)
            ax.legend(fontsize=7, loc="best"); ax.set_ylim(0,1)
            fig.tight_layout(); p=f"{FIG}/TC_ntc_fraction_{COND}.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")

    # ---- overall NTC-fraction figure: 3 panels (one per condition) ----
    conds_done=[c for c in CONDITIONS if results.get(c)]
    if conds_done:
        fig, axes = plt.subplots(1, len(conds_done), figsize=(4.4*len(conds_done), 3.8), squeeze=False)
        for j,COND in enumerate(conds_done):
            ax=axes[0,j]
            for D in [d for d in DONORS if d in results[COND]]:
                cl=pd.DataFrame(results[COND][D]["per_cluster"]); cl["ci"]=cl["cluster"].astype(int)
                cl=cl.sort_values("ci")
                ax.plot(range(len(cl)), cl["ntc_frac"], marker="o", ms=3.0, lw=0.9,
                        color=DONOR_COLORS[D], label=D)
            gN=np.mean([results[COND][d]["global_ntc_frac"] for d in results[COND]])
            ax.axhline(gN, color="#333", ls=":", lw=0.8, alpha=0.7)
            ax.set_title(f"{COND}", fontsize=9, color=COND_COLORS.get(COND,"#333"))
            ax.set_xlabel("Leiden cluster (sorted)"); ax.set_ylim(0,1)
            if j==0: ax.set_ylabel("NTC fraction in cluster")
            ax.legend(fontsize=7, loc="upper right", ncol=2)
        fig.suptitle("NTC fraction per cluster across 3 conditions x 4 donors — flat lines = state-driven structure", fontsize=9.5, x=0.02, ha="left")
        fig.tight_layout(rect=[0,0,1,0.95]); p=f"{FIG}/TC_ntc_fraction_all.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")

    json.dump(results, open(f"{OUT}/threecond_summary.json","w"), indent=2)
    done=[(c,d) for c in CONDITIONS for d in DONORS if results.get(c,{}).get(d)]
    log(f"ALL DONE in {time.time()-t0:.0f}s; combos: {done}")

if __name__=="__main__":
    main()
