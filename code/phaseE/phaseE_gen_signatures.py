#!/usr/bin/env python
# Generate top-150 pseudobulk signatures (by |zscore|) for all 19 Phase E targets x 3 conditions,
# from DE_stats.h5ad. Same structure as Phase D flagship_pseudobulk_signatures.json:
#   {gene_to_ensg, signatures: {"SYM_cond": {ensg_cond_row, ontarget_significant,
#       ontarget_effect_size, n_downstream, top_genes[150 ENSG], top_zscore[150], top_log_fc[150]}}}
import json, numpy as np, h5py
PROJ="/Users/meng01/qiuchen/project/hackathon/perturb-seq"
DE=PROJ+"/perturb-seq_data/GWCD4i.DE_stats.h5ad"
OUT=PROJ+"/phaseE_outputs/phaseE_signatures.json"

ENSG={"LAT":"ENSG00000213658","PLCG1":"ENSG00000124181","CD247":"ENSG00000198821","CD3E":"ENSG00000198851",
 "ZAP70":"ENSG00000115085","VAV1":"ENSG00000141968","SMARCE1":"ENSG00000073584","SGF29":"ENSG00000176476",
 "MED24":"ENSG00000008838","MED12":"ENSG00000184634","TADA2B":"ENSG00000173011","NSD1":"ENSG00000165671",
 "CHD4":"ENSG00000111642","SMARCB1":"ENSG00000099956","STAT6":"ENSG00000166888","TRIP12":"ENSG00000153827",
 "ARNT":"ENSG00000143437","SEL1L":"ENSG00000071537","SIK3":"ENSG00000160584"}
TARGETS=list(ENSG.keys())
CONDS=["Rest","Stim8hr","Stim48hr"]
NTOP=150

f=h5py.File(DE,"r")
def rc(name):
    g=f["obs"][name]
    if isinstance(g,h5py.Group):
        cats=g["categories"][:]; cats=cats.astype(str) if cats.dtype.kind in "SUO" else cats
        return np.asarray(cats)[g["codes"][:]]
    arr=g[:]; return arr.astype(str) if arr.dtype.kind in "SUO" else arr

tcg=rc("target_contrast_gene_name")
cc=rc("culture_condition")
ontarget_sig=rc("ontarget_significant")
ontarget_eff=f["obs/ontarget_effect_size"][:]
n_down=rc("n_downstream")
var_ensg=f["var/gene_ids"][:].astype(str) if "gene_ids" in f["var"] else f["var/_index"][:].astype(str)

zscore_ds=f["layers/zscore"]; logfc_ds=f["layers/log_fc"]

sig={}
for sym in TARGETS:
    for cond in CONDS:
        mask=(tcg==sym)&(cc==cond)
        rows=np.where(mask)[0]
        if rows.size==0:
            continue
        row=int(rows[0])  # one (perturbation x condition) row
        z=zscore_ds[row,:].astype(np.float64)
        lfc=logfc_ds[row,:].astype(np.float64)
        finite=np.isfinite(z)
        order=np.argsort(-np.abs(np.where(finite,z,0)))[:NTOP]
        sig[f"{sym}_{cond}"]={
            "ensg_cond_row":row,
            "ontarget_significant":bool(ontarget_sig[row]=="True") if ontarget_sig.dtype.kind in "SUO" else bool(ontarget_sig[row]),
            "ontarget_effect_size":float(ontarget_eff[row]),
            "n_downstream":float(n_down[row]) if n_down.dtype.kind in "iuf" else float(n_down[row]),
            "top_genes":[str(x) for x in var_ensg[order]],
            "top_zscore":[round(float(x),4) for x in z[order]],
            "top_log_fc":[round(float(x),4) for x in lfc[order]],
        }
f.close()
out={"gene_to_ensg":ENSG,"signatures":sig}
json.dump(out,open(OUT,"w"))
print("WROTE",OUT,"n_signatures",len(sig))
print("targets_with_all3:",sum(1 for t in TARGETS if all(f"{t}_{c}" in sig for c in CONDS)))
