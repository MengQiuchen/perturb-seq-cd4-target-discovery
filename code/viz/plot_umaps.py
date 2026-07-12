"""
UMAP figure generator for the CD4+ Perturb-seq atlas embeddings.

Reads an embedded .h5ad checkpoint (X_umap in obsm, obs with condition/donor/
guide_type/perturbed_gene_name/is_nominated, lognorm layer for gene coloring)
and writes a standard family of figures:

  A. condition-comparison : UMAP faceted by condition (+ colored-by-condition overlay
                            + activation-marker gene panels)
  B. within-condition     : per condition, targeting vs non-targeting overlay
  C. cross-donor          : (multi-donor only) faceted by donor; targeting vs NT
  D. per-target highlight : top-N nominated targets, one panel each (target cells
                            highlighted over a grey background of all cells)
  E. per-module grid      : nominated targets grouped by co-regulation module /
                            therapeutic direction (the target-group figures)

Figure-style rules baked in (skill figure-style):
  §6.6 embeddings: no ticks/labels, corner axis-arrow pair, leader-line cluster labels
  §4 color: threaded CVD-safe palette, grey background context, one alarm hue reserved
  §2 label economy: title = takeaway, n stated, minimal annotations
"""
import os, sys, json, time
import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrow

def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

# ---- figure style (captured from apply_figure_style) ----
STYLE = {
    "figure.dpi": 200, "savefig.dpi": 300, "font.size": 8.0,
    "axes.titlesize": 8.0, "axes.labelsize": 8.0, "xtick.labelsize": 6.0,
    "ytick.labelsize": 6.0, "legend.fontsize": 7.0, "axes.spines.top": False,
    "axes.spines.right": False, "axes.linewidth": 0.6, "axes.titlelocation": "left",
    "axes.titleweight": "normal", "font.family": "sans-serif",
    "legend.frameon": False, "figure.facecolor": "white", "savefig.bbox": "tight",
}
mpl.rcParams.update(STYLE)

# CVD-safe qualitative palette (Okabe-Ito), alarm hue (vermillion) reserved for perturbation
OKABE = {"blue":"#0072B2","orange":"#E69F00","green":"#009E73","yellow":"#F0E442",
         "skyblue":"#56B4E9","vermillion":"#D55E00","purple":"#CC79A7","black":"#000000"}
COND_COLORS = {"Rest":"#0072B2", "Stim8hr":"#E69F00", "Stim48hr":"#009E73"}
DONOR_COLORS = {"D1":"#0072B2","D2":"#E69F00","D3":"#009E73","D4":"#CC79A7"}
GREY = "#D9D9D9"
HILITE = "#D55E00"  # alarm/focal hue for highlighted cells

def _clean_embed_ax(ax):
    """§6.6: strip ticks/labels/spines for an embedding scatter."""
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xlabel(""); ax.set_ylabel("")

def _axis_arrows(ax, label="UMAP"):
    """§6.6: small corner arrow-pair naming the axes."""
    x0 = 0.02; y0 = 0.02; L = 0.13
    ax.annotate("", xy=(x0+L, y0), xytext=(x0, y0), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="#444", lw=0.8))
    ax.annotate("", xy=(x0, y0+L), xytext=(x0, y0), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color="#444", lw=0.8))
    ax.text(x0+L+0.01, y0, f"{label}1", transform=ax.transAxes, fontsize=5.5,
            ha="left", va="center", color="#444")
    ax.text(x0, y0+L+0.01, f"{label}2", transform=ax.transAxes, fontsize=5.5,
            ha="center", va="bottom", color="#444", rotation=90)

def _point_size(n):
    if n > 300_000: return 0.6
    if n > 100_000: return 1.2
    if n > 30_000:  return 2.5
    return 5.0

def _rasterized_scatter(ax, xy, mask=None, color=None, s=None, alpha=0.5, zorder=1):
    pts = xy if mask is None else xy[mask]
    ax.scatter(pts[:,0], pts[:,1], s=s, c=color, alpha=alpha, linewidths=0,
               rasterized=True, zorder=zorder)

def leiden_labels(ax, xy, leiden):
    """Label leiden clusters with leader lines to text in whitespace (light-touch)."""
    for cl in pd.unique(leiden):
        m = leiden == cl
        cx, cy = xy[m,0].mean(), xy[m,1].mean()
        ax.text(cx, cy, str(cl), fontsize=5.5, ha="center", va="center",
                color="#222", zorder=6,
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.6))

# ---------- Figure A: condition comparison ----------
def fig_condition_comparison(ad_umap, obs, outdir, tag, marker_genes=None, lognorm_getter=None):
    conds = [c for c in ["Rest","Stim8hr","Stim48hr"] if c in set(obs["condition"])]
    xy = ad_umap
    n = xy.shape[0]
    s = _point_size(n)

    # A1: faceted by condition
    ncol = len(conds)
    fig, axes = plt.subplots(1, ncol, figsize=(2.6*ncol, 2.8), squeeze=False)
    axes = axes[0]
    for ax, cnd in zip(axes, conds):
        m = (obs["condition"] == cnd).values
        _rasterized_scatter(ax, xy, color=GREY, s=s, alpha=0.35, zorder=1)  # context
        _rasterized_scatter(ax, xy, mask=m, color=COND_COLORS[cnd], s=s, alpha=0.55, zorder=2)
        _clean_embed_ax(ax)
        ax.set_title(f"{cnd}  (n={m.sum():,})", fontsize=8)
    _axis_arrows(axes[0])
    fig.suptitle(f"CD4$^+$ T-cell states shift with TCR stimulation — {tag}", fontsize=9, x=0.02, ha="left")
    fig.tight_layout(rect=[0,0,1,0.94])
    p = f"{outdir}/A1_condition_facets_{tag}.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")

    # A2: single overlay colored by condition
    fig, ax = plt.subplots(figsize=(3.6,3.4))
    for cnd in conds:
        m = (obs["condition"]==cnd).values
        _rasterized_scatter(ax, xy, mask=m, color=COND_COLORS[cnd], s=s, alpha=0.5, zorder=2)
    _clean_embed_ax(ax); _axis_arrows(ax)
    handles = [Line2D([0],[0], marker='o', ls='', mfc=COND_COLORS[c], mec='none', ms=5,
               label=f"{c} (n={(obs['condition']==c).sum():,})") for c in conds]
    ax.legend(handles=handles, loc="upper right", fontsize=6.5, handletextpad=0.2)
    ax.set_title(f"Joint embedding colored by condition — {tag}", fontsize=8)
    fig.tight_layout(); p=f"{outdir}/A2_condition_overlay_{tag}.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")

    # A3: activation-marker gene panels (if genes available)
    if marker_genes and lognorm_getter is not None:
        avail = [(g,gi) for g,gi in marker_genes if gi is not None]
        if avail:
            ncol = len(avail)
            fig, axes = plt.subplots(1, ncol, figsize=(2.3*ncol, 2.5), squeeze=False); axes=axes[0]
            for ax,(g,gi) in zip(axes, avail):
                vals = lognorm_getter(gi)
                order = np.argsort(vals)  # plot high-expressing on top
                scat = ax.scatter(xy[order,0], xy[order,1], c=vals[order], s=s, cmap="viridis",
                                  alpha=0.7, linewidths=0, rasterized=True,
                                  vmin=0, vmax=max(np.percentile(vals,99),1e-3))
                _clean_embed_ax(ax); ax.set_title(g, fontsize=8, style="italic")
                cb = fig.colorbar(scat, ax=ax, fraction=0.046, pad=0.02)
                cb.ax.tick_params(labelsize=5); cb.outline.set_visible(False)
            _axis_arrows(axes[0])
            fig.suptitle(f"Activation markers on the joint embedding — {tag}  (log-norm expr)", fontsize=9, x=0.02, ha="left")
            fig.tight_layout(rect=[0,0,1,0.92]); p=f"{outdir}/A3_activation_markers_{tag}.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")

# ---------- Figure B: within-condition targeting vs non-targeting ----------
def fig_within_condition(ad_umap, obs, outdir, tag):
    conds = [c for c in ["Rest","Stim8hr","Stim48hr"] if c in set(obs["condition"])]
    xy = ad_umap; n=xy.shape[0]; s=_point_size(n)
    ncol=len(conds)
    fig, axes = plt.subplots(1, ncol, figsize=(2.6*ncol,2.9), squeeze=False); axes=axes[0]
    for ax,cnd in zip(axes,conds):
        cm = (obs["condition"]==cnd).values
        nt = cm & (obs["guide_type"]=="non-targeting").values
        tg = cm & (obs["guide_type"]=="targeting").values
        _rasterized_scatter(ax, xy, mask=tg, color="#B9C7D6", s=s, alpha=0.35, zorder=1)
        _rasterized_scatter(ax, xy, mask=nt, color=HILITE, s=max(s,2.0), alpha=0.7, zorder=3)
        _clean_embed_ax(ax); ax.set_title(f"{cnd}", fontsize=8)
    _axis_arrows(axes[0])
    handles=[Line2D([0],[0],marker='o',ls='',mfc="#B9C7D6",mec='none',ms=5,label="targeting (perturbed)"),
             Line2D([0],[0],marker='o',ls='',mfc=HILITE,mec='none',ms=5,label="non-targeting control")]
    axes[-1].legend(handles=handles, loc="upper right", fontsize=6.5, handletextpad=0.2)
    fig.suptitle(f"Perturbed vs non-targeting control cells within each condition — {tag}", fontsize=9, x=0.02, ha="left")
    fig.tight_layout(rect=[0,0,1,0.93]); p=f"{outdir}/B_within_condition_targeting_{tag}.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")

# ---------- Figure C: cross-donor ----------
def fig_cross_donor(ad_umap, obs, outdir, tag):
    donors = sorted(set(obs["donor"]))
    if len(donors) < 2: return
    xy=ad_umap; n=xy.shape[0]; s=_point_size(n)
    ncol=len(donors)
    fig, axes = plt.subplots(1, ncol, figsize=(2.5*ncol,2.8), squeeze=False); axes=axes[0]
    for ax,dn in zip(axes,donors):
        m=(obs["donor"]==dn).values
        _rasterized_scatter(ax, xy, color=GREY, s=s, alpha=0.3, zorder=1)
        _rasterized_scatter(ax, xy, mask=m, color=DONOR_COLORS.get(dn,"#333"), s=s, alpha=0.55, zorder=2)
        _clean_embed_ax(ax); ax.set_title(f"{dn}  (n={m.sum():,})", fontsize=8)
    _axis_arrows(axes[0])
    fig.suptitle(f"Cross-donor reproducibility of the embedding — {tag}", fontsize=9, x=0.02, ha="left")
    fig.tight_layout(rect=[0,0,1,0.93]); p=f"{outdir}/C_cross_donor_{tag}.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")

# ---------- Figure D: per-target highlight ----------
def fig_per_target(ad_umap, obs, outdir, tag, targets, ncol=5):
    xy=ad_umap; n=xy.shape[0]; s=_point_size(n)
    targets=[t for t in targets if (obs["perturbed_gene_name"]==t).any()]
    if not targets: log("  no targets present for per-target panel"); return
    nrow=int(np.ceil(len(targets)/ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.05*ncol, 2.05*nrow), squeeze=False)
    for i,ax in enumerate(axes.flat):
        if i>=len(targets): ax.axis("off"); continue
        t=targets[i]; m=(obs["perturbed_gene_name"]==t).values
        _rasterized_scatter(ax, xy, color=GREY, s=s, alpha=0.3, zorder=1)
        ax.scatter(xy[m,0], xy[m,1], s=max(s*3,4), c=HILITE, alpha=0.8, linewidths=0, rasterized=True, zorder=3)
        _clean_embed_ax(ax); ax.set_title(f"{t}  (n={m.sum():,})", fontsize=7.5, style="italic")
    _axis_arrows(axes[0,0])
    fig.suptitle(f"Cells by knockdown target (CRISPRi) — {tag}", fontsize=9, x=0.02, ha="left")
    fig.tight_layout(rect=[0,0,1,0.965]); p=f"{outdir}/D_per_target_{tag}.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")

# ---------- Figure E: per-module grid ----------
def fig_per_module(ad_umap, obs, outdir, tag, module_map, ncol=4):
    """module_map: dict group_name -> list of target genes."""
    xy=ad_umap; n=xy.shape[0]; s=_point_size(n)
    groups=[(g,ts) for g,ts in module_map.items() if any((obs["perturbed_gene_name"]==t).any() for t in ts)]
    if not groups: log("  no modules present"); return
    nrow=int(np.ceil(len(groups)/ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.15*ncol, 2.15*nrow), squeeze=False)
    for i,ax in enumerate(axes.flat):
        if i>=len(groups): ax.axis("off"); continue
        g,ts=groups[i]; m=np.isin(obs["perturbed_gene_name"].values, ts)
        _rasterized_scatter(ax, xy, color=GREY, s=s, alpha=0.3, zorder=1)
        ax.scatter(xy[m,0], xy[m,1], s=max(s*2,3), c=HILITE, alpha=0.7, linewidths=0, rasterized=True, zorder=3)
        _clean_embed_ax(ax)
        ntg=len([t for t in ts if (obs['perturbed_gene_name']==t).any()])
        tw = "target" if ntg==1 else "targets"
        ax.set_title(f"{g}\n({ntg} {tw}, n={m.sum():,})", fontsize=7)
    _axis_arrows(axes[0,0])
    fig.suptitle(f"Cells grouped by co-regulation module / therapeutic direction — {tag}", fontsize=9, x=0.02, ha="left")
    fig.tight_layout(rect=[0,0,1,0.95]); p=f"{outdir}/E_per_module_{tag}.png"; fig.savefig(p); plt.close(fig); log(f"  wrote {p}")
