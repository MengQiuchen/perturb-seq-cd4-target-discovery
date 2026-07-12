"""
Per-MODULE / per-DIRECTION UMAP grids colored by CONDITION (not single highlight hue).
The by-condition counterpart of Figure E (fig_per_module). Groups are built identically to
run_plot.py: by scorecard module_name (E1) and therapeutic_direction (E2). Within each panel,
the group's cells are colored by their source condition (Rest/Stim8hr[/Stim48hr]) over a grey
all-cell background; draw order is shuffled so no condition systematically occludes another.

Usage: python per_module_bycondition.py <embedded_checkpoint.h5ad> <tag>
Writes:
  figures/E_bycondition_per_module_<tag>.png            (co-regulation modules, by condition)
  figures/E_bycondition_per_module_<tag>_bydirection.png (therapeutic directions, by condition)
  checkpoints/<tag>.per_module_condition_counts.csv      (group x condition count table)
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
RNG=np.random.default_rng(0)

STYLE={"figure.dpi":200,"savefig.dpi":300,"font.size":8.0,"axes.titlesize":8.0,
       "axes.labelsize":8.0,"xtick.labelsize":6.0,"ytick.labelsize":6.0,"legend.fontsize":7.0,
       "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":0.6,
       "axes.titlelocation":"left","font.family":"sans-serif","legend.frameon":False,
       "figure.facecolor":"white","savefig.bbox":"tight"}
mpl.rcParams.update(STYLE)

def _point_size(n):
    if n>300_000: return 0.6
    if n>100_000: return 1.2
    if n>30_000:  return 2.5
    return 5.0
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

def build_groups(sc_df, present, col):
    gm={}
    if col in sc_df.columns:
        for k,sub in sc_df.dropna(subset=[col]).groupby(col):
            genes=[g for g in sub["gene"].astype(str) if g in present]
            if genes: gm[str(k)]=genes
    return gm

def plot_grid(xy, obs, group_map, conds, tag, out_png, ncol, suptitle):
    n=xy.shape[0]; s=_point_size(n)
    groups=[(g,ts) for g,ts in group_map.items()
            if any((obs["perturbed_gene_name"]==t).any() for t in ts)]
    if not groups: log(f"  no groups for {out_png}"); return None
    nrow=int(np.ceil(len(groups)/ncol))
    fig,axes=plt.subplots(nrow,ncol,figsize=(2.15*ncol,2.15*nrow),squeeze=False)
    rows=[]
    for i,ax in enumerate(axes.flat):
        if i>=len(groups): ax.axis("off"); continue
        g,ts=groups[i]; m=np.isin(obs["perturbed_gene_name"].values, ts)
        idx=np.where(m)[0]
        ax.scatter(xy[:,0],xy[:,1],s=s,c=GREY,alpha=0.25,linewidths=0,rasterized=True,zorder=1)
        cond_arr=obs["condition"].values[idx]
        cols=np.array([COND_COLORS.get(c,"#333") for c in cond_arr])
        perm=RNG.permutation(len(idx))  # fair overplot: shuffle draw order across conditions
        ax.scatter(xy[idx[perm],0],xy[idx[perm],1],s=max(s*2,3),c=cols[perm],
                   alpha=0.75,linewidths=0,rasterized=True,zorder=3)
        _clean(ax)
        ntg=len([t for t in ts if (obs['perturbed_gene_name']==t).any()])
        tw="target" if ntg==1 else "targets"
        ax.set_title(f"{g}\n({ntg} {tw}, n={m.sum():,})", fontsize=7)
        row={"group":g,"n_targets":ntg,"total":int(m.sum())}
        for cnd in conds: row[cnd]=int((m & (obs["condition"]==cnd).values).sum())
        rows.append(row)
    _arrows(axes[0,0])
    handles=[Line2D([0],[0],marker='o',ls='',mfc=COND_COLORS[c],mec='none',ms=5,label=c) for c in conds]
    # title on top line (left), legend stacked on a SECOND line just below it (left) — never collide
    if nrow>1:
        title_y, leg_y, top = 1.0, 0.965, 0.90
    else:
        title_y, leg_y, top = 1.06, 0.99, 0.90
    fig.suptitle(suptitle, fontsize=9, x=0.02, y=title_y, ha="left")
    fig.legend(handles=handles, loc="upper left", ncol=len(conds), fontsize=7,
               bbox_to_anchor=(0.02, leg_y), frameon=False)
    fig.tight_layout(rect=[0,0,1,top]); fig.savefig(out_png, bbox_inches="tight"); plt.close(fig)
    log(f"  wrote {out_png}")
    return pd.DataFrame(rows)

def main():
    ckpt, tag = sys.argv[1], sys.argv[2]
    log(f"loading {ckpt}")
    adata=ad.read_h5ad(ckpt)
    xy=adata.obsm["X_umap"]; obs=adata.obs
    conds=[c for c in ["Rest","Stim8hr","Stim48hr"] if c in set(obs["condition"])]
    sc_df=pd.read_csv(SCORECARD)
    present=set(obs["perturbed_gene_name"].unique())

    module_map=build_groups(sc_df, present, "module_name")
    dir_map={f"dir: {k[5:] if k.startswith('dir: ') else k}":v for k,v in build_groups(sc_df, present, "therapeutic_direction").items()}
    log(f"modules: {list(module_map)}; directions: {list(dir_map)}")

    dfm=plot_grid(xy,obs,module_map,conds,tag,
                  f"{FIG}/E_bycondition_per_module_{tag}.png", 4,
                  f"Module cells colored by condition (Rest/Stim) — {tag}")
    dfd=plot_grid(xy,obs,dir_map,conds,tag,
                  f"{FIG}/E_bycondition_per_module_{tag}_bydirection.png", 2,
                  f"Therapeutic-direction cells colored by condition — {tag}")
    # counts table
    parts=[]
    if dfm is not None: dfm["kind"]="module"; parts.append(dfm)
    if dfd is not None: dfd["kind"]="direction"; parts.append(dfd)
    if parts:
        allc=pd.concat(parts, ignore_index=True)
        allc.to_csv(f"{CK}/{tag}.per_module_condition_counts.csv", index=False)
        log(f"wrote {CK}/{tag}.per_module_condition_counts.csv ({len(allc)} groups)")

if __name__=="__main__":
    main()
