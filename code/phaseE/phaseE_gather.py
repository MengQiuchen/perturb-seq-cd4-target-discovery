#!/usr/bin/env python
# Phase E slab-streaming gather (adapted from proven Phase D gather_reference.py, frame 49de3a80).
# Keeps: all 19-target KD cells + 40k NTC + 40k background (other targeting genes), low_quality==False.
# low_quality is a BOOL array in these shards (spec's "True"/"False" string is wrong) -> use ~lq.
# Usage: phaseE_gather.py <shard_path> <cond_tag> <out_h5ad> <summary_json> <seed>
import sys, time, json, os, numpy as np, h5py
from scipy import sparse
import scanpy as sc, anndata as ad, pandas as pd

SHARD, COND, OUT_H5, OUT_SUM, SEED = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
t0=time.time()
def log(*a): print("[%7.1fs]"%(time.time()-t0),*a,flush=True)

TARGETS=["LAT","PLCG1","CD247","CD3E","ZAP70","VAV1","SMARCE1",
         "SGF29","MED24","MED12","TADA2B","NSD1","CHD4","SMARCB1",
         "STAT6","TRIP12","ARNT","SEL1L","SIK3"]
TARGET_SET=set(TARGETS)
N_NTC=40000; N_BG=40000
rng=np.random.default_rng(SEED)

f=h5py.File(SHARD,"r",rdcc_nbytes=0); og=f["obs"]
def rc(n):
    g=og[n]
    if isinstance(g,h5py.Group):
        cats=g["categories"][:]; cats=cats.astype(str) if cats.dtype.kind in "SUO" else cats
        return np.asarray(cats)[g["codes"][:]]
    arr=g[:]
    return arr.astype(str) if arr.dtype.kind in "SUO" else arr

gt=rc("guide_type")
lq=rc("low_quality")           # BOOL array
pgn=rc("perturbed_gene_name")  # target symbol per cell
lane=rc("lane_id")
ng=og["n_genes_by_counts"][:]; tc=og["total_counts"][:]; mt=og["pct_counts_mt"][:]

# low_quality may arrive as bool or as str; normalize to a "is good quality" mask
if lq.dtype==bool:
    goodq = ~lq
else:
    goodq = (lq.astype(str)!="True")
log("n_cells",gt.size,"goodq",int(goodq.sum()),"lq_dtype",str(lq.dtype))

is_tgt = (gt=="targeting") & goodq & np.isin(pgn, TARGETS)
is_ntc = (gt=="non-targeting") & goodq
is_bg_pool = (gt=="targeting") & goodq & (~np.isin(pgn, TARGETS))
log("target-KD cells", int(is_tgt.sum()), "| NTC pool", int(is_ntc.sum()), "| bg pool", int(is_bg_pool.sum()))

# sample NTC and background
ntc_idx=np.where(is_ntc)[0]
bg_idx =np.where(is_bg_pool)[0]
ntc_keep = rng.choice(ntc_idx, size=min(N_NTC,ntc_idx.size), replace=False)
bg_keep  = rng.choice(bg_idx,  size=min(N_BG, bg_idx.size),  replace=False)
tgt_keep = np.where(is_tgt)[0]

keep=np.zeros(gt.size, dtype=bool)
keep[tgt_keep]=True; keep[ntc_keep]=True; keep[bg_keep]=True
# class code: 1 target, 2 ntc, 3 background
cls=np.zeros(gt.size, dtype=np.int8)
cls[tgt_keep]=1; cls[ntc_keep]=2; cls[bg_keep]=3
idx=np.where(keep)[0]
log("total kept", idx.size, "(target %d + ntc %d + bg %d)"%(tgt_keep.size, ntc_keep.size, bg_keep.size))

# --- slab-streaming CSR read (proven pattern) ---
indptr=f["X/indptr"][:].astype(np.int64); nnz=int(indptr[-1])
ncols=int(f["X"].attrs["shape"][1]) if "shape" in f["X"].attrs else 18130
data_ds=f["X/data"]; ind_ds=f["X/indices"]
rl=(indptr[idx+1]-indptr[idx]).astype(np.int64); tot=int(rl.sum())
od=np.empty(tot,dtype=np.float32); oi=np.empty(tot,dtype=np.int32); wp=0
SLAB=500_000_000; nsl=(nnz+SLAB-1)//SLAB
log("nnz",nnz,"ncols",ncols,"kept_nnz",tot,"slabs",nsl)
for si in range(nsl):
    lo=si*SLAB; hi=min(lo+SLAB,nnz); d=data_ds[lo:hi]; ii=ind_ds[lo:hi]
    rs=max(np.searchsorted(indptr,lo,side="right")-1,0); re=np.searchsorted(indptr,hi,side="left")
    rr=np.arange(rs,re); sl=np.maximum(indptr[rr],lo); sh=np.minimum(indptr[rr+1],hi)
    counts=np.clip((sh-sl),0,None).astype(np.int64); km=np.repeat(keep[rr],counts)
    kd=d[km]; ki=ii[km]; od[wp:wp+kd.size]=kd; oi[wp:wp+ki.size]=ki; wp+=kd.size
    log(" slab",si+1,"/",nsl,"wp",wp)
assert wp==tot, (wp,tot)
nip=np.zeros(idx.size+1,dtype=np.int64); nip[1:]=np.cumsum(rl)
X=sparse.csr_matrix((od,oi,nip),shape=(idx.size,ncols))
gene_name=f["var/gene_name"][:].astype(str)
# ENSG var index (for concordance matching); prefer gene_ids else _index
if "gene_ids" in f["var"]:
    ensg=f["var/gene_ids"][:].astype(str)
else:
    ensg=f["var/_index"][:].astype(str)
f.close()

# class label + per-cell target symbol aligned to idx
clsi=cls[idx]
class_name=np.where(clsi==1,"target",np.where(clsi==2,"NTC","background"))
pgn_i=pgn[idx].copy()
pgn_i[clsi==2]="non-targeting"   # NTC label
obs=pd.DataFrame({
    "cell_class":class_name,
    "perturbed_gene_name":pgn_i,
    "guide_type":gt[idx],
    "lane_id":lane[idx],
    "n_genes_by_counts":ng[idx],
    "total_counts":tc[idx],
    "pct_counts_mt":mt[idx],
    "condition":COND,
})
A=ad.AnnData(X=X,obs=obs)
A.var["gene_name"]=gene_name
A.var["ensg"]=ensg
A.var_names=pd.Index(gene_name).astype(str); A.var_names_make_unique()
log("built AnnData",A.shape)

# normalize
A.layers["counts"]=A.X.copy()
sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
log("normalized")
A.write_h5ad(OUT_H5)
log("WROTE",OUT_H5)

# summary
per_tgt={t:int(np.sum((clsi==1)&(pgn_i==t))) for t in TARGETS}
summ={"cond":COND,"shard":os.path.basename(SHARD),"n_cells":int(A.n_obs),"n_genes":int(A.n_vars),
      "n_target_cells":int((clsi==1).sum()),"n_ntc":int((clsi==2).sum()),"n_bg":int((clsi==3).sum()),
      "per_target_KD":per_tgt,"seed":SEED,
      "ntc_pool":int(is_ntc.sum()),"bg_pool":int(is_bg_pool.sum())}
json.dump(summ,open(OUT_SUM,"w"),indent=1)
log("WROTE",OUT_SUM)
print("___GATHER_DONE___",json.dumps(per_tgt))
