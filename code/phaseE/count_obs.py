
import h5py, numpy as np, pandas as pd, time, sys, json
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
def dec(x): return x.decode() if isinstance(x,(bytes,bytearray)) else x
PROJ="/Users/meng01/qiuchen/project/hackathon/perturb-seq"
CL=f"{PROJ}/perturb-seq_data/cell_level"
shards=[("Rest",f"{CL}/D4_Rest.assigned_guide.h5ad"),
        ("Stim8hr",f"{CL}/D4_Stim8hr.assigned_guide.h5ad"),
        ("Stim48hr",f"{CL}/D4_Stim48hr.assigned_guide.h5ad")]

def read_cat(obs, col):
    o=obs[col]
    if isinstance(o,h5py.Group) and "categories" in o:
        cats=np.array([dec(x) for x in o["categories"][:]],dtype=object)
        codes=o["codes"][:]
        vals=np.empty(codes.shape[0],dtype=object); good=codes>=0
        vals[good]=cats[codes[good]]; vals[~good]="NA"
        return vals
    return np.array([dec(x) for x in o[:]],dtype=object)

per_cond={}   # cond -> Series (gene -> count) for targeting, low_quality-filtered
gt_counts={}  # cond -> dict guide_type->count (post low_quality filter)
n_obs_tot={}
for cond,fn in shards:
    log(f"reading obs: {cond}")
    f=h5py.File(fn,"r"); obs=f["obs"]
    n=obs["_index"].shape[0]; n_obs_tot[cond]=int(n)
    gt=read_cat(obs,"guide_type")
    pg=read_cat(obs,"perturbed_gene_name")
    lq=obs["low_quality"][:].astype(bool) if "low_quality" in obs else np.zeros(n,bool)
    f.close()
    ok=~lq
    gt_counts[cond]={"non-targeting":int(((gt=="non-targeting")&ok).sum()),
                     "targeting":int(((gt=="targeting")&ok).sum()),
                     "low_quality_dropped":int(lq.sum()), "n_total":int(n)}
    tmask=(gt=="targeting")&ok
    s=pd.Series(pg[tmask]).value_counts()
    per_cond[cond]=s
    log(f"  {cond}: n={n:,} tgt={tmask.sum():,} nt={gt_counts[cond]['non-targeting']:,} lq={lq.sum():,} genes={s.shape[0]:,}")

# combine into a target x condition table (ALL targeting genes, uncapped)
allg=sorted(set().union(*[set(s.index) for s in per_cond.values()]))
tab=pd.DataFrame(index=allg)
for cond,_ in shards:
    tab[cond]=per_cond[cond].reindex(allg).fillna(0).astype(int)
tab["total"]=tab[[c for c,_ in shards]].sum(axis=1)
tab=tab.sort_values("total",ascending=False)
tab.index.name="perturbed_gene_name"

# annotate with scorecard nominated + rank + module + direction
sc=pd.read_csv(f"{PROJ}/phaseBC_outputs/TARGET_SCORECARD.csv")
sc["gene"]=sc["gene"].astype(str)
rankmap=dict(zip(sc["gene"],sc["rank"]))
modmap=dict(zip(sc["gene"],sc["module_name"]))
dirmap=dict(zip(sc["gene"],sc["therapeutic_direction"]))
nomset=set(sc["gene"])
tab["is_nominated"]=[g in nomset for g in tab.index]
tab["rank"]=[rankmap.get(g,np.nan) for g in tab.index]
tab["module_name"]=[modmap.get(g,None) for g in tab.index]
tab["therapeutic_direction"]=[dirmap.get(g,None) for g in tab.index]

OUT=f"{PROJ}/phaseE_outputs/umap_within"
import os; os.makedirs(OUT,exist_ok=True)
tab.to_csv(f"{OUT}/D4_uncapped_per_target_counts.csv")
json.dump({"gt_counts":gt_counts,"n_obs_tot":n_obs_tot,
           "n_targeting_genes":int(tab['is_nominated'].shape[0]),
           "n_nominated_present":int(tab['is_nominated'].sum())},
          open(f"{OUT}/D4_obs_summary.json","w"),indent=2)

log("=== TOP 30 targeting genes by TOTAL cells (uncapped, low_quality-filtered) ===")
print(tab.head(30).to_string())
log("=== TOP 20 NOMINATED targets by TOTAL cells ===")
print(tab[tab["is_nominated"]].head(20).to_string())
log(f"DONE. wrote {OUT}/D4_uncapped_per_target_counts.csv")
