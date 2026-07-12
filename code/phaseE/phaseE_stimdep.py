#!/usr/bin/env python
# Phase E stim-dependence: for each target, mean |KD-vs-NTC log-FC| over the target's
# Stim8hr top-150 signature genes, computed within each condition subset. Ratio of the
# effect magnitude between conditions quantifies stimulation-dependence.
# Compares D4_Stim8hr vs D4_Stim48hr (both gathered). Stim-dependent = ratio>1.5.
# Usage: phaseE_stimdep.py <sig_json> <out_csv> <s8_h5ad> <s48_h5ad>
import sys, json, numpy as np, pandas as pd, scanpy as sc
from scipy import sparse
SIG, OUT, S8, S48 = sys.argv[1:5]
TARGETS=["LAT","PLCG1","CD247","CD3E","ZAP70","VAV1","SMARCE1",
 "SGF29","MED24","MED12","TADA2B","NSD1","CHD4","SMARCB1","STAT6","TRIP12","ARNT","SEL1L","SIK3"]
sig=json.load(open(SIG))["signatures"]

def mean_abs_delta(h5, cond_tag):
    A=sc.read_h5ad(h5)
    ensg=A.var["ensg"].astype(str).values if "ensg" in A.var else A.var_names.astype(str).values
    ensg2col={}
    for i,e in enumerate(ensg): ensg2col.setdefault(e,i)
    cls=A.obs["cell_class"].astype(str).values; pgn=A.obs["perturbed_gene_name"].astype(str).values
    ntc=np.where(cls=="NTC")[0]
    Xln=A.X
    ntc_mean=np.asarray(Xln[ntc].mean(axis=0)).ravel()
    out={}
    for tgt in TARGETS:
        kd=np.where((cls=="target")&(pgn==tgt))[0]
        if kd.size==0: out[tgt]=(np.nan,0); continue
        kd_mean=np.asarray(Xln[kd].mean(axis=0)).ravel()
        sc_lfc=kd_mean-ntc_mean
        # use the target's STIM8HR signature genes as the reference gene set (per methods)
        s=sig.get(f"{tgt}_Stim8hr") or sig.get(f"{tgt}_Stim48hr") or sig.get(f"{tgt}_Rest")
        genes=s["top_genes"]
        cols=np.array([ensg2col.get(g,-1) for g in genes]); ok=cols>=0
        vals=np.abs(sc_lfc[cols[ok]])
        out[tgt]=(float(np.nanmean(vals)), int(kd.size))
    return out

d8=mean_abs_delta(S8,"Stim8hr"); d48=mean_abs_delta(S48,"Stim48hr")
rows=[]
for t in TARGETS:
    m8,n8=d8[t]; m48,n48=d48[t]
    ratio8=m8/m48 if (m48 and m48>0) else np.nan
    rows.append({"gene":t,"mean_abs_delta_Stim8hr":round(m8,4) if np.isfinite(m8) else np.nan,
                 "n_kd_Stim8hr":n8,"mean_abs_delta_Stim48hr":round(m48,4) if np.isfinite(m48) else np.nan,
                 "n_kd_Stim48hr":n48,"ratio_8hr_over_48hr":round(ratio8,3) if np.isfinite(ratio8) else np.nan,
                 "stim8_dominant":bool(np.isfinite(ratio8) and ratio8>1.5),
                 "stim48_dominant":bool(np.isfinite(ratio8) and ratio8<(1/1.5))})
pd.DataFrame(rows).to_csv(OUT,index=False)
print("___STIMDEP_DONE___", OUT)
for r in rows: print(r["gene"], r["mean_abs_delta_Stim8hr"], r["mean_abs_delta_Stim48hr"], r["ratio_8hr_over_48hr"])
