"""
Request 2 — PER-TARGET UMAP (one fresh embedding per target, cells pooled across
ALL conditions, colored by condition).

For each chosen knockdown target, this gathers *all* of that target's cells (not the
per-target-capped subsample used by the joint waves) from the 3 raw D4 shards, then
computes an INDEPENDENT UMAP on just that target's cells and colors by condition.
This asks, for a single perturbation: do its cells still separate by activation state
(Rest / Stim8hr / Stim48hr) on their own manifold, and where do the conditions sit?

Target list is read from a JSON file (phaseE_scripts/pt_targets.json):
  {"targets": ["GENE1","GENE2", ...]}

Gather step streams each raw shard's CSR ONCE (embed_common.read_masked_X) keeping only
rows whose perturbed_gene_name is in the chosen set AND low_quality==False (no caps),
then caches the union AnnData so re-plots skip the ~30-40 min gather.

Outputs (phaseE_outputs/umap_within):
  checkpoints/pertarget/union_gather.h5ad         all chosen targets' cells, raw counts
  checkpoints/pertarget/{gene}.obs_umap.csv.gz    per-target obs + UMAP coords
  figures/pertarget/PT_{gene}_bycondition.png     single per-target UMAP, colored by condition
  figures/pertarget/PT_grid_bycondition.png       grid: one panel per target, colored by condition
  per_target_umap_summary.json
"""
import os, sys, time, json
import numpy as np, pandas as pd
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import scipy.sparse as sp
import anndata as ad

sys.path.insert(0, "/Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/code")
from embed_common import log, read_obs_cols, read_masked_X

PROJ = "/Users/meng01/qiuchen/project/hackathon/perturb-seq"
CL = f"{PROJ}/perturb-seq_data/cell_level"
OUT = f"{PROJ}/phaseE_outputs/umap_within"
CKDIR = f"{OUT}/checkpoints/pertarget"; FIG = f"{OUT}/figures/pertarget"
os.makedirs(CKDIR, exist_ok=True); os.makedirs(FIG, exist_ok=True)
SCORECARD = f"{PROJ}/phaseBC_outputs/TARGET_SCORECARD.csv"

SHARDS = [("Rest",   f"{CL}/D4_Rest.assigned_guide.h5ad"),
          ("Stim8hr",f"{CL}/D4_Stim8hr.assigned_guide.h5ad"),
          ("Stim48hr",f"{CL}/D4_Stim48hr.assigned_guide.h5ad")]
COND_COLORS = {"Rest":"#0072B2","Stim8hr":"#E69F00","Stim48hr":"#009E73"}
GREY = "#D9D9D9"

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
    if n>20000: return 2.5
    if n>5000:  return 5.0
    if n>1500:  return 8.0
    return 14.0

def read_targeted_rows(fn, keep_mask, log_every=2000):
    """Read ONLY the CSR rows where keep_mask is True (targeted per-row slices).
    Avoids the full-array scan of read_masked_X; touches ~nnz(targets) not nnz(all)."""
    import h5py, numpy as np, scipy.sparse as sp, time
    kidx = np.nonzero(keep_mask)[0]
    with h5py.File(fn, "r", rdcc_nbytes=256*1024*1024, rdcc_nslots=100003) as f:
        X = f["X"]
        shp = X.attrs["shape"] if "shape" in X.attrs else None
        ncols = int(shp[1]) if shp is not None else int(f["var"][list(f["var"].keys())[0]].shape[0])
        t0=time.time(); indptr = X["indptr"][:]
        log(f"    indptr read {time.time()-t0:.0f}s ({indptr.shape[0]-1:,} rows)")
        data_ds = X["data"]; ind_ds = X["indices"]
        rows_data=[]; rows_ind=[]; new_indptr=[0]; total=0; t0=time.time()
        for j, r in enumerate(kidx):
            a=int(indptr[r]); b=int(indptr[r+1])
            if b>a:
                rows_data.append(data_ds[a:b]); rows_ind.append(ind_ds[a:b])
            total += (b-a); new_indptr.append(total)
            if (j+1) % log_every == 0:
                el=time.time()-t0
                log(f"    targeted {j+1:,}/{len(kidx):,} rows {total:,} nnz {el:.0f}s ({(j+1)/max(el,1e-9):.0f} rows/s)")
        data = np.concatenate(rows_data) if rows_data else np.empty(0, dtype=data_ds.dtype)
        indices = np.concatenate(rows_ind) if rows_ind else np.empty(0, dtype=ind_ds.dtype)
        Xc = sp.csr_matrix((data, indices, np.array(new_indptr, dtype=np.int64)), shape=(len(kidx), ncols))
        log(f"    built CSR {Xc.shape} nnz={Xc.nnz:,} in {time.time()-t0:.0f}s")
    return Xc, kidx

def gather_union(targets):
    """Stream each raw shard once; keep cells of chosen targets (no cap), low_quality-filtered."""
    ckf = f"{CKDIR}/union_gather.h5ad"
    if os.path.exists(ckf):
        log(f"union cache HIT {ckf}")
        A = ad.read_h5ad(ckf)
        have = set(A.obs["perturbed_gene_name"].unique())
        if set(targets).issubset(have):
            return A
        log(f"cache missing {set(targets)-have}; re-gathering")
    tset = set(targets); parts=[]; var_ref=None
    for cond, fn in SHARDS:
        log(f"[{cond}] reading obs for mask")
        f = __import__("h5py").File(fn,"r")
        obs = read_obs_cols(f, ["guide_type","perturbed_gene_name","low_quality","guide_id","lane_id"])
        # var (index = ensembl, gene_name = symbol)
        var = f["var"]; vidx = var.attrs.get("_index", b"_index")
        vidx = vidx.decode() if isinstance(vidx,(bytes,bytearray)) else vidx
        var_index = np.array([x.decode() if isinstance(x,(bytes,bytearray)) else x for x in var[vidx][:]], dtype=object)
        gene_name=None
        if "gene_name" in var:
            g=var["gene_name"]
            if hasattr(g,"keys") and "categories" in g:
                cats=np.array([x.decode() if isinstance(x,(bytes,bytearray)) else x for x in g["categories"][:]],dtype=object)
                gene_name=cats[g["codes"][:]]
            else:
                gene_name=np.array([x.decode() if isinstance(x,(bytes,bytearray)) else x for x in g[:]],dtype=object)
        f.close()
        pg = obs["perturbed_gene_name"]; lq = obs["low_quality"]
        ok = ~lq.astype(bool) if lq is not None else np.ones(len(pg),bool)
        keep = np.isin(pg, list(tset)) & ok
        log(f"[{cond}] keep {keep.sum():,}/{len(pg):,} cells (targets present: {sorted(set(pg[keep]))[:6]}...)")
        if keep.sum()==0: continue
        # Storage recovered (114 MB/s slab reads) -> use slab-streaming read_masked_X
        # (contiguous nnz scan, ~17 min/shard) instead of the per-row targeted reader,
        # which was a workaround for the degraded OST and is pathologically slow on
        # scattered small reads. Same (X, kept_idx) return signature.
        X, kidx = read_masked_X(fn, keep)
        odf = pd.DataFrame({"perturbed_gene_name": pg[kidx], "guide_type": obs["guide_type"][kidx],
                            "guide_id": obs["guide_id"][kidx] if obs["guide_id"] is not None else "NA",
                            "lane_id": obs["lane_id"][kidx] if obs["lane_id"] is not None else "NA"})
        odf["condition"]=cond; odf["donor"]="D4"
        odf.index=[f"D4_{cond}_{i}" for i in kidx]
        a=ad.AnnData(X=X, obs=odf); a.var_names=var_index
        if gene_name is not None: a.var["gene_name"]=gene_name
        if var_ref is None: var_ref=(var_index, gene_name)
        parts.append(a)
    A = ad.concat(parts, join="outer", index_unique=None)
    vi,gn = var_ref; A.var_names=vi
    if gn is not None: A.var["gene_name"]=gn
    tmp=ckf+".tmp"; A.write_h5ad(tmp); os.replace(tmp, ckf)
    log(f"union gathered {A.n_obs:,} cells; cached -> {ckf}")
    return A

def embed_one(sub, seed=0):
    """Independent embedding of one target's cells."""
    import scanpy as sc
    navail=len(os.sched_getaffinity(0)) if hasattr(os,"sched_getaffinity") else 8
    njobs=int(os.environ.get("EMBED_NJOBS", str(min(8,max(1,navail-1)))))
    try:
        import numba; njobs=min(njobs, numba.config.NUMBA_NUM_THREADS)
    except Exception: pass
    sc.settings.n_jobs=njobs
    sc.pp.normalize_total(sub, target_sum=1e4); sc.pp.log1p(sub)
    sub.layers["lognorm"]=sub.X.copy()
    n=sub.n_obs
    nhvg=min(2000, max(200, sub.n_vars//5))
    sc.pp.highly_variable_genes(sub, n_top_genes=nhvg, flavor="seurat", subset=False)
    npcs=min(50, n-1, int((sub.var["highly_variable"].sum())))
    npcs=max(2, min(npcs, 50))
    sc.pp.pca(sub, n_comps=npcs, use_highly_variable=True, zero_center=True, random_state=seed)
    nnb=min(15, max(3, n//10))
    sc.pp.neighbors(sub, n_neighbors=nnb, n_pcs=npcs, use_rep="X_pca", random_state=seed)
    sc.tl.umap(sub, random_state=seed)
    return sub

def plot_one(sub, gene, path):
    xy=sub.obsm["X_umap"]; cond=sub.obs["condition"].astype(str).values; n=xy.shape[0]; s=_psize(n)
    conds=[c for c in ["Rest","Stim8hr","Stim48hr"] if c in set(cond)]
    fig, ax = plt.subplots(figsize=(3.8,3.6))
    # fair overplot: shuffle draw order so no condition systematically hides another
    rng=np.random.default_rng(0); order=rng.permutation(n)
    colors=np.array([COND_COLORS[c] for c in cond],dtype=object)
    ax.scatter(xy[order,0],xy[order,1], s=s, c=list(colors[order]), alpha=0.75, linewidths=0, rasterized=True)
    _clean(ax); _arrows(ax)
    handles=[Line2D([0],[0],marker='o',ls='',mfc=COND_COLORS[c],mec='none',ms=5,
             label=f"{c} (n={(cond==c).sum():,})") for c in conds]
    ax.legend(handles=handles, loc="upper right", fontsize=6.5, handletextpad=0.2)
    ax.set_title(f"{gene} — per-target UMAP colored by condition (n={n:,})", fontsize=8, style="italic")
    fig.tight_layout(); fig.savefig(path); plt.close(fig); log(f"  wrote {path}")

def main():
    t0=time.time()
    cfg=json.load(open(f"{PROJ}/phaseE_scripts/pt_targets.json"))
    targets=cfg["targets"]
    log(f"targets: {targets}")
    A = gather_union(targets)
    sc_df=pd.read_csv(SCORECARD); sc_df["gene"]=sc_df["gene"].astype(str)
    rankmap=dict(zip(sc_df["gene"], sc_df["rank"]))
    summary={}; embedded={}
    for gene in targets:
        sub = A[A.obs["perturbed_gene_name"]==gene].copy()
        if sub.n_obs < 30:
            log(f"[{gene}] only {sub.n_obs} cells — skipping (too few for UMAP)"); continue
        log(f"[{gene}] embedding {sub.n_obs:,} cells")
        sub = embed_one(sub)
        odf=sub.obs.copy(); odf["UMAP1"]=sub.obsm["X_umap"][:,0]; odf["UMAP2"]=sub.obsm["X_umap"][:,1]
        odf.to_csv(f"{CKDIR}/{gene}.obs_umap.csv.gz", compression="gzip")
        plot_one(sub, gene, f"{FIG}/PT_{gene}_bycondition.png")
        cc=sub.obs["condition"].value_counts().to_dict()
        summary[gene]={"n_cells":int(sub.n_obs),"rank":(None if pd.isna(rankmap.get(gene,np.nan)) else int(rankmap.get(gene))),
                       "by_condition":{k:int(v) for k,v in cc.items()}}
        embedded[gene]=(sub.obsm["X_umap"], sub.obs["condition"].astype(str).values)
        del sub
        log(f"[{gene}] done ({time.time()-t0:.0f}s)")

    # combined grid
    genes=[g for g in targets if g in embedded]
    if genes:
        ncol=min(4,len(genes)); nrow=int(np.ceil(len(genes)/ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(2.6*ncol, 2.6*nrow), squeeze=False)
        for i,ax in enumerate(axes.flat):
            if i>=len(genes): ax.axis("off"); continue
            g=genes[i]; xy,cond=embedded[g]; n=xy.shape[0]; s=_psize(n)
            rng=np.random.default_rng(0); order=rng.permutation(n)
            colors=np.array([COND_COLORS[c] for c in cond],dtype=object)
            ax.scatter(xy[order,0],xy[order,1], s=s, c=list(colors[order]), alpha=0.75, linewidths=0, rasterized=True)
            _clean(ax)
            rk=rankmap.get(g,np.nan); rtxt="" if pd.isna(rk) else f" · rank {int(rk)}"
            ax.set_title(f"{g}  (n={n:,}{rtxt})", fontsize=7.5, style="italic")
        _arrows(axes[0,0])
        conds_all=[c for c in ["Rest","Stim8hr","Stim48hr"]]
        handles=[Line2D([0],[0],marker='o',ls='',mfc=COND_COLORS[c],mec='none',ms=5,label=c) for c in conds_all]
        fig.legend(handles=handles, loc="upper right", ncol=3, fontsize=7, bbox_to_anchor=(0.99,1.0))
        fig.suptitle("Per-target UMAPs (each target embedded on its own cells) colored by condition — D4", fontsize=9, x=0.02, ha="left")
        fig.tight_layout(rect=[0,0,1,0.95]); p=f"{FIG}/PT_grid_bycondition.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")

    json.dump(summary, open(f"{OUT}/per_target_umap_summary.json","w"), indent=2)
    log(f"ALL DONE in {time.time()-t0:.0f}s; {len(genes)} targets plotted")

if __name__=="__main__":
    main()
