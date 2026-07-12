#!/usr/bin/env python
"""Phase A — Foundation: schema/QC audit, reproducibility landscape, positive-control benchmark.
Runs on rocm-4 in the perturb-seq conda env. Operates on GWCD4i.DE_stats.h5ad (precomputed DESeq2).
All outputs -> $PROJ/phaseA_outputs/.
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import TwoSlopeNorm

PROJ = "/Users/meng01/qiuchen/project/hackathon/perturb-seq"
D    = PROJ + "/perturb-seq_data"
OUT  = "/mnt/scratche/slow/ghlab/qiuchen/tmp/regen_out"
os.makedirs(OUT, exist_ok=True)

# ---- publication figure style (correctness rules baked in; no house style) ----
mpl.rcParams.update({
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlelocation": "left", "axes.titleweight": "normal",
    "figure.facecolor": "white", "axes.facecolor": "white",
    "font.family": "DejaVu Sans", "svg.fonttype": "none",
})
FOCAL = "#2166ac"; ACCENT = "#b2182b"; GREY = "#9e9e9e"
DIVERGE = "RdBu_r"  # diverging, centered at 0 for signed zscore
SEQ = "viridis"

print("[load] reading DE_stats.h5ad into memory ...", flush=True)
A = ad.read_h5ad(f"{D}/GWCD4i.DE_stats.h5ad")
obs = A.obs.copy()
gene_name = A.var["gene_name"].astype(str).values
gene_to_col = {g: i for i, g in enumerate(gene_name)}
zscore = A.layers["zscore"]  # (33983 x 10282), dense float
logfc  = A.layers["log_fc"]
adjp   = A.layers["adj_p_value"]
print(f"[load] done. shape={A.shape}", flush=True)

CONDS = ["Rest", "Stim8hr", "Stim48hr"]
COND_ORDER = {c: i for i, c in enumerate(CONDS)}

# =====================================================================
# A1 — SCHEMA & QC AUDIT
# =====================================================================
print("[A1] schema & QC audit ...", flush=True)
schema_rows = []
for col in obs.columns:
    s = obs[col]
    dt = str(s.dtype)
    n_na = int(s.isna().sum())
    if s.dtype == bool:
        desc = f"True={int(s.sum())} ({100*s.mean():.1f}%)"
    elif str(s.dtype) == "category":
        nu = s.nunique()
        top = s.value_counts().head(3)
        desc = f"{nu} categories; top: " + ", ".join(f"{k}={v}" for k, v in top.items())
    else:
        sn = pd.to_numeric(s, errors="coerce")
        desc = f"min={sn.min():.3g} med={sn.median():.3g} max={sn.max():.3g}"
    schema_rows.append({"field": col, "dtype": dt, "n_missing": n_na, "summary": desc})
de_schema = pd.DataFrame(schema_rows)
de_schema.to_csv(f"{OUT}/de_schema.csv", index=False)
print(de_schema.to_string(index=False), flush=True)

# knockdown efficiency: on-target effect (negative = knockdown worked)
n_rows = A.n_obs
n_targets = obs["target_contrast_gene_name"].nunique()
ontarget_sig_rate = float(obs["ontarget_significant"].mean())
eff = obs["ontarget_effect_size"].astype(float)
# per-condition perturbation counts
per_cond = obs["culture_condition"].value_counts().reindex(CONDS)
# exclusion flags
flag_cols = ["distal_offtarget_flag", "neighboring_gene_KD", "low_target_gex"]
flag_rates = {c: float(obs[c].mean()) for c in flag_cols}

qc_summary = {
    "n_perturbation_condition_rows": int(n_rows),
    "n_unique_targets": int(n_targets),
    "n_genes_measured": int(A.n_vars),
    "rows_per_condition": {c: int(per_cond[c]) for c in CONDS},
    "ontarget_significant_rate": ontarget_sig_rate,
    "ontarget_effect_size_median": float(eff.median()),
    "flag_rates": flag_rates,
    "n_guides_dist": {str(k): int(v) for k, v in obs["n_guides"].value_counts().items()},
}
json.dump(qc_summary, open(f"{OUT}/qc_summary.json", "w"), indent=2)
print("[A1]", json.dumps(qc_summary, indent=2), flush=True)

# ---- FIGURE: qc_overview.png (2x3) ----
fig, axs = plt.subplots(2, 3, figsize=(11, 6.2))

ax = axs[0, 0]
bars = ax.bar(range(len(CONDS)), [per_cond[c] for c in CONDS], color=FOCAL, width=0.62)
ax.set_xticks(range(len(CONDS))); ax.set_xticklabels(CONDS)
ax.set_ylabel("perturbation×condition rows"); ax.set_title("a  Tested perturbations per context", loc="left", fontweight="bold")
for b, c in zip(bars, CONDS):
    ax.text(b.get_x()+b.get_width()/2, b.get_height()+80, f"{int(per_cond[c]):,}", ha="center", va="bottom", fontsize=6)
ax.margins(y=0.12)

ax = axs[0, 1]
ax.hist(eff.clip(-30, 8), bins=60, color=FOCAL)
ax.axvline(0, color=GREY, lw=0.8, ls="--")
ax.axvline(eff.median(), color=ACCENT, lw=1.2)
ax.text(eff.median()-1.5, ax.get_ylim()[1]*0.9, f"median\n{eff.median():.1f}", color=ACCENT, ha="right", fontsize=6)
ax.set_xlabel("on-target z-score (log2FC / lfcSE of target gene)")
ax.set_ylabel("perturbations")
ax.set_title("b  Knockdown efficiency (negative = KD worked)", loc="left", fontweight="bold")

ax = axs[0, 2]
sig = obs.groupby("culture_condition", observed=True)["ontarget_significant"].mean().reindex(CONDS)
bars = ax.bar(range(len(CONDS)), sig.values*100, color=FOCAL, width=0.62)
ax.set_xticks(range(len(CONDS))); ax.set_xticklabels(CONDS)
ax.set_ylabel("% on-target significant"); ax.set_ylim(0, 100)
ax.set_title("c  Knockdown detected (on-target significant)", loc="left", fontweight="bold")
for b, v in zip(bars, sig.values):
    ax.text(b.get_x()+b.get_width()/2, v*100+1.5, f"{v*100:.0f}%", ha="center", va="bottom", fontsize=6)

ax = axs[1, 0]
ncells = obs["n_cells_target"].astype(float)
ax.hist(np.log10(ncells.clip(lower=1)), bins=60, color=FOCAL)
med = np.log10(ncells.median())
ax.axvline(med, color=ACCENT, lw=1.2)
ax.text(med+0.05, ax.get_ylim()[1]*0.9, f"median\n{int(ncells.median())} cells", color=ACCENT, fontsize=6)
ax.set_xlabel("cells per target (log10)")
ax.set_ylabel("perturbations")
ax.set_title("d  Cells per target", loc="left", fontweight="bold")

ax = axs[1, 1]
ndl = obs["n_total_de_genes"].astype(float)
ax.hist(np.log10(ndl.clip(lower=0.5)), bins=60, color=FOCAL)
ax.set_xlabel("downstream DE genes per perturbation (log10, 10% FDR)")
ax.set_ylabel("perturbations")
ax.set_title("e  Trans-effect magnitude (hub-ness)", loc="left", fontweight="bold")

ax = axs[1, 2]
fr = pd.Series(flag_rates)*100
names = ["distal\noff-target", "neighbor\ngene KD", "low target\nexpression"]
bars = ax.barh(range(len(fr)), fr.values, color=ACCENT)
ax.set_yticks(range(len(fr))); ax.set_yticklabels(names)
ax.invert_yaxis()
ax.set_xlabel("% of perturbations flagged")
ax.set_title("f  Exclusion flags", loc="left", fontweight="bold")
for i, v in enumerate(fr.values):
    ax.text(v+0.3, i, f"{v:.1f}%", va="center", fontsize=6)
ax.margins(x=0.12)

fig.suptitle("Phase A · QC audit — GWCD4i genome-scale CRISPRi Perturb-seq (33,983 perturbation×context rows; 11,526 targets; 10,282 genes)",
             x=0.01, ha="left", fontsize=9, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(f"{OUT}/qc_overview.png")
plt.close(fig)
print("[A1] saved qc_overview.png + de_schema.csv", flush=True)

# =====================================================================
# A2 — REPRODUCIBILITY LANDSCAPE + HIGH-CONFIDENCE SET
# =====================================================================
print("[A2] reproducibility landscape ...", flush=True)
# two-guide agreement: guides significantly positively correlated
guide_ok = (obs["guide_correlation_all"].astype(float) > 0) & \
           (obs["guide_correlation_all_pval"].astype(float) < 0.05)
# cross-donor agreement: mean hit correlation across donor pairs positive & present
donor_ok = obs["donor_correlation_hits_mean"].astype(float) > 0.3
# on-target: knockdown detected AND meaningful effect
ontarget_ok = obs["ontarget_significant"].astype(bool) & (obs["ontarget_effect_size"].astype(float) < 0)
# exclusion flags clean
clean = (~obs["distal_offtarget_flag"].astype(bool)) & \
        (~obs["neighboring_gene_KD"].astype(bool)) & \
        (~obs["low_target_gex"].astype(bool))
reproducible = guide_ok.fillna(False) | donor_ok.fillna(False)
high_conf = ontarget_ok & reproducible & clean

obs["_guide_ok"] = guide_ok.fillna(False)
obs["_donor_ok"] = donor_ok.fillna(False)
obs["_ontarget_ok"] = ontarget_ok
obs["_clean"] = clean
obs["_high_confidence"] = high_conf

hc = obs[high_conf].copy()
# has downstream trans-effect
hc_active = hc[hc["n_total_de_genes"] >= 1]
hc_out = hc.reset_index().rename(columns={"index": "perturbation_id"})
keep_cols = ["perturbation_id", "target_contrast_gene_name", "target_contrast", "culture_condition",
             "ontarget_effect_size", "ontarget_significant",
             "guide_correlation_all", "guide_correlation_all_pval",
             "donor_correlation_hits_mean", "n_downstream", "n_total_de_genes",
             "n_up_genes", "n_down_genes", "n_cells_target", "target_baseMean",
             "_guide_ok", "_donor_ok"]
keep_cols = [c for c in keep_cols if c in hc_out.columns]
hc_out[keep_cols].to_csv(f"{OUT}/high_confidence_targets.csv", index=False)

repro_summary = {
    "total_rows": int(n_rows),
    "ontarget_ok": int(ontarget_ok.sum()),
    "guide_agreement": int(guide_ok.fillna(False).sum()),
    "donor_agreement": int(donor_ok.fillna(False).sum()),
    "reproducible_either": int(reproducible.sum()),
    "clean_flags": int(clean.sum()),
    "high_confidence": int(high_conf.sum()),
    "high_confidence_with_transeffect": int((high_conf & (obs["n_total_de_genes"] >= 1)).sum()),
    "high_confidence_unique_targets": int(hc["target_contrast_gene_name"].nunique()),
    "hc_per_condition": {c: int((hc["culture_condition"] == c).sum()) for c in CONDS},
}
json.dump(repro_summary, open(f"{OUT}/reproducibility_summary.json", "w"), indent=2)
print("[A2]", json.dumps(repro_summary, indent=2), flush=True)

# ---- FIGURE: reproducibility.png (2x2) ----
fig, axs = plt.subplots(2, 2, figsize=(9.5, 7.4))

ax = axs[0, 0]
gc = obs["guide_correlation_all"].astype(float).dropna()
ax.hist(gc, bins=60, color=FOCAL)
ax.axvline(0, color=GREY, lw=0.8, ls="--")
frac_pos = float((guide_ok.fillna(False)).sum())/n_rows
ax.set_xlabel("two-guide correlation (all genes)")
ax.set_ylabel("perturbations")
ax.set_title(f"a  Guide-1 vs guide-2 agreement\n{repro_summary['guide_agreement']:,} pass (r>0, p<0.05)", loc="left", fontweight="bold")

ax = axs[0, 1]
dc = obs["donor_correlation_hits_mean"].astype(float).dropna()
ax.hist(dc, bins=40, color=FOCAL)
ax.axvline(0.3, color=ACCENT, lw=1.2)
ax.text(0.32, ax.get_ylim()[1]*0.85, "threshold\nr>0.3", color=ACCENT, fontsize=6)
ax.set_xlabel("cross-donor correlation (mean over hit genes)")
ax.set_ylabel("perturbations")
ax.set_title(f"b  Cross-donor agreement\n{repro_summary['donor_agreement']:,} pass; {int(dc.shape[0]):,} evaluable", loc="left", fontweight="bold")

# funnel
ax = axs[1, 0]
steps = ["all rows", "on-target\nKD", "+ reproducible\n(guide|donor)", "+ flag-clean\n= high-conf"]
vals = [n_rows, int(ontarget_ok.sum()), int((ontarget_ok & reproducible).sum()), int(high_conf.sum())]
bars = ax.bar(range(len(steps)), vals, color=[GREY, FOCAL, FOCAL, ACCENT], width=0.66)
ax.set_xticks(range(len(steps))); ax.set_xticklabels(steps, fontsize=6)
ax.set_ylabel("perturbations")
ax.set_title("c  High-confidence filter funnel", loc="left", fontweight="bold")
for b, v in zip(bars, vals):
    ax.text(b.get_x()+b.get_width()/2, v+250, f"{v:,}", ha="center", va="bottom", fontsize=6)
ax.margins(y=0.14)

# high-conf per condition, split active/inactive
ax = axs[1, 1]
act = [int(((hc["culture_condition"] == c) & (hc["n_total_de_genes"] >= 1)).sum()) for c in CONDS]
ina = [int(((hc["culture_condition"] == c) & (hc["n_total_de_genes"] < 1)).sum()) for c in CONDS]
x = np.arange(len(CONDS))
ax.bar(x, act, color=FOCAL, width=0.62, label="with trans-effect (≥1 DE gene)")
ax.bar(x, ina, bottom=act, color=GREY, width=0.62, label="no trans-effect")
ax.set_xticks(x); ax.set_xticklabels(CONDS)
ax.set_ylabel("high-confidence perturbations")
ax.set_title("d  High-confidence set by context", loc="left", fontweight="bold")
ax.legend(frameon=False, fontsize=6, loc="upper left")
for i, (a, b) in enumerate(zip(act, ina)):
    ax.text(i, a+b+30, f"{a+b:,}", ha="center", va="bottom", fontsize=6)
ax.margins(y=0.12)

fig.suptitle(f"Phase A · Reproducibility landscape — high-confidence set = {high_conf.sum():,} perturbations ({hc['target_contrast_gene_name'].nunique():,} unique targets)",
             x=0.01, ha="left", fontsize=9, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(f"{OUT}/reproducibility.png")
plt.close(fig)
print("[A2] saved reproducibility.png + high_confidence_targets.csv", flush=True)

# =====================================================================
# A3 — POSITIVE-CONTROL BENCHMARK
# =====================================================================
print("[A3] positive-control benchmark ...", flush=True)
# curated regulator -> expected downstream readouts & direction.
# CRISPRi = loss-of-function. For an ACTIVATOR of a program, KD => downstream DOWN (negative z).
PC = {
    # Th subset master TFs (activators of their lineage cytokines)
    "TBX21":  {"program": "Th1",  "activates": ["IFNG", "CXCR3"]},
    "GATA3":  {"program": "Th2",  "activates": ["IL4", "IL5", "IL13"]},
    "RORC":   {"program": "Th17", "activates": ["IL23R", "CCR6", "IL17F"]},
    "FOXP3":  {"program": "Treg", "activates": ["IL2RA", "CTLA4", "IKZF2"]},
    "BCL6":   {"program": "Tfh",  "activates": ["CXCR5"]},
    # IL-2 / JAK-STAT axis
    "STAT5A": {"program": "IL2-STAT5", "activates": ["IL2RA"]},
    "STAT5B": {"program": "IL2-STAT5", "activates": ["IL2RA"]},
    "IL2RA":  {"program": "IL2R",  "activates": ["IL2RA"]},
    # TCR / costim / calcineurin-NFAT (activators of activation cytokines)
    "LCK":    {"program": "TCR",   "activates": ["IL2", "IFNG"]},
    "ZAP70":  {"program": "TCR",   "activates": ["IL2", "IFNG"]},
    "CD28":   {"program": "costim","activates": ["IL2"]},
    "PPP3CA": {"program": "calcineurin","activates": ["IL2", "IFNG"]},
    "NFATC1": {"program": "NFAT",  "activates": ["IL2", "IFNG"]},
    "NFATC2": {"program": "NFAT",  "activates": ["IL2", "IFNG"]},
    # broad effector TFs
    "IRF4":   {"program": "effector", "activates": ["IL2", "IFNG"]},
    "BATF":   {"program": "effector", "activates": ["IFNG", "IL2"]},
}
# readout gene panel (columns of the heatmap)
READOUTS = ["IFNG", "CXCR3", "IL4", "IL5", "IL13", "IL23R", "CCR6", "IL17F",
            "IL2", "IL2RA", "CTLA4", "IKZF2", "CXCR5", "TNF", "IL21", "GZMB"]
READOUTS = [g for g in READOUTS if g in gene_to_col]

# choose condition per target: prefer Stim48hr (cytokine biology active), else Stim8hr, else Rest
def row_for(target, cond):
    m = (obs["target_contrast_gene_name"].astype(str) == target) & (obs["culture_condition"].astype(str) == cond)
    idx = np.where(m.values)[0]
    return int(idx[0]) if len(idx) else None

records = []
heat = []  # rows = PC targets present, cols = readouts, value = zscore
heat_labels = []
selfz = []
for tgt, meta in PC.items():
    # pick condition
    chosen = None
    for cond in ["Stim48hr", "Stim8hr", "Rest"]:
        r = row_for(tgt, cond)
        if r is not None:
            chosen = (cond, r); break
    if chosen is None:
        continue
    cond, ri = chosen
    # on-target self knockdown
    self_col = gene_to_col.get(tgt)
    self_z = float(zscore[ri, self_col]) if self_col is not None else np.nan
    ote = float(obs["ontarget_effect_size"].iloc[ri])
    osig = bool(obs["ontarget_significant"].iloc[ri])
    # expected downstream readouts
    exp = meta["activates"]
    exp_present = [g for g in exp if g in gene_to_col]
    exp_z = {g: float(zscore[ri, gene_to_col[g]]) for g in exp_present}
    # PASS if majority of expected activation targets go DOWN (negative z) under KD
    downs = [v for v in exp_z.values() if not np.isnan(v)]
    n_expected_down = sum(1 for v in downs if v < 0)
    direction_pass = (len(downs) > 0) and (n_expected_down >= (len(downs)+1)//2)
    records.append({
        "target": tgt, "program": meta["program"], "condition_used": cond,
        "ontarget_effect_size": ote, "ontarget_significant": osig, "self_zscore": self_z,
        "n_expected_readouts": len(downs), "n_readouts_down": n_expected_down,
        "expected_direction_pass": bool(direction_pass),
        **{f"z[{g}]": exp_z.get(g, np.nan) for g in exp_present},
        "n_total_de_genes": int(obs["n_total_de_genes"].iloc[ri]),
    })
    # heatmap row across full readout panel
    heat.append([float(zscore[ri, gene_to_col[g]]) for g in READOUTS])
    heat_labels.append(f"{tgt} ({meta['program']}, {cond[:5]})")
    selfz.append(self_z)

pc_df = pd.DataFrame(records)
pc_df.to_csv(f"{OUT}/positive_controls.csv", index=False)
print(pc_df.to_string(index=False), flush=True)
n_pass = int(pc_df["expected_direction_pass"].sum())
print(f"[A3] direction PASS: {n_pass}/{len(pc_df)} controls", flush=True)

# ---- FIGURE: positive_controls.png ----
heat = np.array(heat)
fig = plt.figure(figsize=(11, 6.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 2.6], wspace=0.28)

# left: on-target self-knockdown (validates the on-target axis)
ax0 = fig.add_subplot(gs[0, 0])
order = np.argsort(pc_df["ontarget_effect_size"].values)
lab = pc_df["target"].values[order]
val = pc_df["ontarget_effect_size"].values[order]
sig = pc_df["ontarget_significant"].values[order]
colors = [FOCAL if s else GREY for s in sig]
ax0.hlines(range(len(val)), 0, val, color=colors, lw=1)
ax0.scatter(val, range(len(val)), c=colors, s=22, zorder=3)
ax0.set_yticks(range(len(lab))); ax0.set_yticklabels(lab, fontsize=6)
ax0.axvline(0, color=GREY, lw=0.8)
ax0.set_xlabel("on-target z-score\n(log2FC / lfcSE of KD'd gene)")
ax0.set_title("a  Knockdown works on positive controls", loc="left", fontweight="bold")
ax0.text(0.02, -0.14, "blue = on-target significant · negative = knockdown detected",
         transform=ax0.transAxes, fontsize=6, color=GREY)

# right: downstream signature heatmap
ax1 = fig.add_subplot(gs[0, 1])
vmax = np.nanpercentile(np.abs(heat), 98)
vmax = max(vmax, 3.0)
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
im = ax1.imshow(heat, aspect="auto", cmap=DIVERGE, norm=norm)
ax1.set_xticks(range(len(READOUTS))); ax1.set_xticklabels(READOUTS, rotation=90, fontsize=6)
ax1.set_yticks(range(len(heat_labels))); ax1.set_yticklabels(heat_labels, fontsize=6)
ax1.set_title("b  Downstream signature of each knockdown (DE zscore)", loc="left", fontweight="bold")
# annotate cells that are expected activation targets
for i, tgt in enumerate(pc_df["target"].values):
    for g in PC[tgt]["activates"]:
        if g in READOUTS:
            j = READOUTS.index(g)
            ax1.add_patch(plt.Rectangle((j-0.5, i-0.5), 1, 1, fill=False, edgecolor="black", lw=1.1))
cb = fig.colorbar(im, ax=ax1, fraction=0.025, pad=0.02)
cb.set_label("DE zscore (KD vs control)\nnegative = downstream gene down", fontsize=6)
ax1.text(0.0, 1.03, "boxed = known target activated by this regulator → expected to drop under KD",
         transform=ax1.transAxes, fontsize=6, color="black")

fig.suptitle(f"Phase A · Positive-control benchmark — known T-cell regulators reproduce expected biology  ({n_pass}/{len(pc_df)} pass expected direction)",
             x=0.01, ha="left", fontsize=9, fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(f"{OUT}/positive_controls.png")
plt.close(fig)
print("[A3] saved positive_controls.png + positive_controls.csv", flush=True)

print("[DONE] Phase A complete. Outputs in", OUT, flush=True)
for f in sorted(os.listdir(OUT)):
    sz = os.path.getsize(f"{OUT}/{f}")
    print(f"  {f}  ({sz/1024:.1f} KB)", flush=True)
