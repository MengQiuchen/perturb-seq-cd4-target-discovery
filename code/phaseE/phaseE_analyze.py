#!/usr/bin/env python
# Phase E per-shard analysis. Reads a compact normalized subset h5ad (from phaseE_gather.py),
# computes: (1) on-target KD efficiency (linear-space % reduction + MWU p), (2) pseudobulk
# concordance (Pearson r of single-cell KD-vs-NTC log-FC vs DE_stats log_fc over top-150 sig),
# (3) 6-program score shifts (KD vs NTC).  Writes per-shard CSVs + JSON.
# Usage: phaseE_analyze.py <subset_h5ad> <tag> <cond> <signatures_json> <outdir>
import sys, json, os, time, numpy as np, pandas as pd, h5py
import scanpy as sc, anndata as ad
from scipy import sparse
from scipy.stats import mannwhitneyu, pearsonr, spearmanr

SUB, TAG, COND, SIGJSON, OUTDIR = sys.argv[1:6]
os.makedirs(OUTDIR, exist_ok=True)
t0=time.time()
def log(*a): print("[%7.1fs]"%(time.time()-t0),*a,flush=True)

TARGETS=["LAT","PLCG1","CD247","CD3E","ZAP70","VAV1","SMARCE1",
 "SGF29","MED24","MED12","TADA2B","NSD1","CHD4","SMARCB1",
 "STAT6","TRIP12","ARNT","SEL1L","SIK3"]
PROGRAMS={
 "naive_memory":["CCR7","SELL","TCF7","LEF1","IL7R"],
 "effector":["GZMB","GZMA","PRF1","IFNG","NKG7","GNLY","CCL5"],
 "activation":["IL2RA","CD69","TNFRSF9","MKI67","TNFRSF18","ICOS"],
 "treg":["FOXP3","IKZF2","CTLA4","IL2RA","TNFRSF18"],
 "cytokine":["IL2","IFNG","TNF","IL4","IL13","IL21","CSF2"],
 "exhaustion":["PDCD1","LAG3","HAVCR2","TIGIT","TOX","ENTPD1"],
}

A=sc.read_h5ad(SUB); log("loaded",A.shape, "conds", A.obs["condition"].unique().tolist())
sig=json.load(open(SIGJSON))["signatures"]

# --- var: symbol + ensg maps (X is log-normalized; layers['counts'] raw) ---
sym=A.var["gene_name"].astype(str).values
ensg=A.var["ensg"].astype(str).values if "ensg" in A.var else A.var_names.astype(str).values
sym2col={}
for i,s in enumerate(sym):
    sym2col.setdefault(s,i)  # first occurrence
ensg2col={}
for i,e in enumerate(ensg):
    ensg2col.setdefault(e,i)

cls=A.obs["cell_class"].astype(str).values
pgn=A.obs["perturbed_gene_name"].astype(str).values
ntc_mask=(cls=="NTC")
Xln=A.X  # csr log-normalized
Xlin=A.layers["counts"] if "counts" in A.layers else None  # raw counts

# precompute NTC mean log-normalized per gene (dense vector) once
ntc_idx=np.where(ntc_mask)[0]
def colmean_ln(rows):
    if len(rows)==0: return np.zeros(Xln.shape[1])
    return np.asarray(Xln[rows].mean(axis=0)).ravel()
ntc_mean_ln=colmean_ln(ntc_idx)
log("NTC cells",ntc_idx.size)

# expm1 for linear-space KD %
def target_expr_vec(colidx, rows):
    # returns per-cell log-norm expression of one gene for given rows
    v=np.asarray(Xln[rows, colidx].todense()).ravel() if sparse.issparse(Xln) else Xln[rows, colidx]
    return v

# ---------- KD efficiency + concordance + present targets ----------
kd_rows=[]; conc_rows=[]
present=[]
for tgt in TARGETS:
    kd_idx=np.where((cls=="target")&(pgn==tgt))[0]
    n_kd=kd_idx.size
    if n_kd==0:
        continue
    present.append(tgt)
    # on-target gene column (by symbol)
    col=sym2col.get(tgt)
    if col is None:
        # fall back: skip KD calc
        kd_rows.append({"gene":tgt,"cond":COND,"n_kd":n_kd,"note":"target gene not in var"})
    else:
        kd_ln=target_expr_vec(col, kd_idx)
        ntc_ln=target_expr_vec(col, ntc_idx)
        kd_lin=np.expm1(kd_ln); ntc_lin=np.expm1(ntc_ln)
        m_kd_lin=float(kd_lin.mean()); m_ntc_lin=float(ntc_lin.mean())
        rel_kd=1.0 - (m_kd_lin/m_ntc_lin) if m_ntc_lin>0 else np.nan
        try:
            p=float(mannwhitneyu(kd_ln, ntc_ln, alternative="less").pvalue)
        except Exception:
            p=np.nan
        kd_rows.append({"gene":tgt,"cond":COND,"n_kd":int(n_kd),
            "mean_kd_lognorm":round(float(kd_ln.mean()),4),"mean_ntc_lognorm":round(float(ntc_ln.mean()),4),
            "mean_kd_linear":round(m_kd_lin,4),"mean_ntc_linear":round(m_ntc_lin,4),
            "rel_knockdown":round(float(rel_kd),4),
            "frac_expr_kd":round(float((kd_ln>0).mean()),3),"frac_expr_ntc":round(float((ntc_ln>0).mean()),3),
            "mwu_p_less":p})
    # concordance vs signature for this cond
    s=sig.get(f"{tgt}_{COND}")
    if s is not None and n_kd>0:
        kd_mean_ln=colmean_ln(kd_idx)
        sc_lfc=kd_mean_ln - ntc_mean_ln  # per-gene single-cell log-FC (log-norm space)
        genes=s["top_genes"]; ref_lfc=np.array(s["top_log_fc"],dtype=float)
        cols=np.array([ensg2col.get(g,-1) for g in genes])
        ok=cols>=0
        sc_vals=np.where(ok, sc_lfc[np.where(ok,cols,0)], np.nan)
        m=ok & np.isfinite(ref_lfc) & np.isfinite(sc_vals)
        n_match=int(m.sum())
        if n_match>=10 and np.std(sc_vals[m])>0 and np.std(ref_lfc[m])>0:
            r=float(pearsonr(sc_vals[m], ref_lfc[m])[0]); rho=float(spearmanr(sc_vals[m], ref_lfc[m])[0])
        else:
            r=np.nan; rho=np.nan
        conc_rows.append({"gene":tgt,"cond":COND,"n_kd":int(n_kd),"n_sig_genes":len(genes),
            "n_matched":n_match,"pearson_r":round(r,4) if np.isfinite(r) else np.nan,
            "spearman_r":round(rho,4) if np.isfinite(rho) else np.nan})
log("KD + concordance done; present targets", len(present))

pd.DataFrame(kd_rows).to_csv(f"{OUTDIR}/{TAG}.kd.csv",index=False)
pd.DataFrame(conc_rows).to_csv(f"{OUTDIR}/{TAG}.concordance.csv",index=False)

# ---------- program score shifts ----------
# score on symbol space (X log-normalized). Use only genes present.
for prog,glist in PROGRAMS.items():
    present_g=[g for g in glist if g in sym2col]
    sc.tl.score_genes(A, present_g, score_name=f"prog_{prog}", use_raw=False)
prog_rows=[]
ntc_scores={p:A.obs.loc[ntc_mask,f"prog_{p}"].mean() for p in PROGRAMS}
for tgt in present:
    kd_idx=np.where((cls=="target")&(pgn==tgt))[0]
    row={"gene":tgt,"cond":COND,"n_kd":int(kd_idx.size)}
    for p in PROGRAMS:
        kd_s=A.obs.iloc[kd_idx][f"prog_{p}"].mean()
        ntc_s=ntc_scores[p]
        row[f"{p}_kd"]=round(float(kd_s),4)
        row[f"{p}_ntc"]=round(float(ntc_s),4)
        row[f"{p}_shift"]=round(float(kd_s-ntc_s),4)
    prog_rows.append(row)
pd.DataFrame(prog_rows).to_csv(f"{OUTDIR}/{TAG}.progshift.csv",index=False)
log("program shifts done")

summ={"tag":TAG,"cond":COND,"subset":os.path.basename(SUB),"n_cells":int(A.n_obs),
      "present_targets":present,"n_present":len(present)}
json.dump(summ,open(f"{OUTDIR}/{TAG}.analysis_summary.json","w"),indent=1)
print("___ANALYZE_DONE___",TAG,"present",len(present))
