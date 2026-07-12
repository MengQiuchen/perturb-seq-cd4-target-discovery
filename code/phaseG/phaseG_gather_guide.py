#!/usr/bin/env python
# Phase G tier-3 guide-aware gather: keeps 19-target KD cells (with guide_id/guide_group) +
# N_NTC non-targeting cells, low_quality==False.  Adapted verbatim from phaseE_gather.py
# slab-streaming; only adds guide_id + guide_group to obs and drops the background pool.
# Usage: phaseG_gather_guide.py <shard> <cond> <out_h5ad> <summary_json> <seed>
import sys, time, json, os, numpy as np, h5py
from scipy import sparse
import scanpy as sc, anndata as ad, pandas as pd

SHARD, COND, OUT_H5, OUT_SUM, SEED = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
t0=time.time()
def log(*a): print("[%7.1fs]"%(time.time()-t0),*a,flush=True)

TARGETS=["LAT","PLCG1","CD247","CD3E","ZAP70","VAV1","SMARCE1",
         "SGF29","MED24","MED12","TADA2B","NSD1","CHD4","SMARCB1",
         "STAT6","TRIP12","ARNT","SEL1L","SIK3"]
N_NTC=40000
rng=np.random.default_rng(SEED)

f=h5py.File(SHARD,"r",rdcc_nbytes=0); og=f["obs"]
def rc(n):
    g=og[n]
    if isinstance(g,h5py.Group):
        cats=g["categories"][:]; cats=cats.astype(str) if cats.dtype.kind in "SUO" else cats
        return np.asarray(cats)[g["codes"][:]]
    arr=g[:]
    return arr.astype(str) if arr.dtype.kind in "SUO" else arr

gt=rc("guide_type"); lq=rc("low_quality"); pgn=rc("perturbed_gene_name"); lane=rc("lane_id")
guide_id=rc("guide_id"); guide_group=rc("guide_group")
ng=og["n_genes_by_counts"][:]; tc=og["total_counts"][:]; mt=og["pct_counts_mt"][:]
goodq = ~lq if lq.dtype==bool else (lq.astype(str)!="True")
log("n_cells",gt.size,"goodq",int(goodq.sum()))

is_tgt = (gt=="targeting") & goodq & np.isin(pgn, TARGETS)
is_ntc = (gt=="non-targeting") & goodq
ntc_idx=np.where(is_ntc)[0]
ntc_keep=rng.choice(ntc_idx, size=min(N_NTC,ntc_idx.size), replace=False)
tgt_keep=np.where(is_tgt)[0]
log("target-KD",tgt_keep.size,"| NTC kept",ntc_keep.size)

keep=np.zeros(gt.size,bool); keep[tgt_keep]=True; keep[ntc_keep]=True
cls=np.zeros(gt.size,np.int8); cls[tgt_keep]=1; cls[ntc_keep]=2
idx=np.where(keep)[0]

indptr=f["X/indptr"][:].astype(np.int64); nnz=int(indptr[-1])
ncols=int(f["X"].attrs["shape"][1]) if "shape" in f["X"].attrs else 18130
data_ds=f["X/data"]; ind_ds=f["X/indices"]
rl=(indptr[idx+1]-indptr[idx]).astype(np.int64); tot=int(rl.sum())
od=np.empty(tot,np.float32); oi=np.empty(tot,np.int32); wp=0
SLAB=500_000_000; nsl=(nnz+SLAB-1)//SLAB
log("nnz",nnz,"kept_nnz",tot,"slabs",nsl)
for si in range(nsl):
    lo=si*SLAB; hi=min(lo+SLAB,nnz); d=data_ds[lo:hi]; ii=ind_ds[lo:hi]
    rs=max(np.searchsorted(indptr,lo,side="right")-1,0); re=np.searchsorted(indptr,hi,side="left")
    rr=np.arange(rs,re); sl=np.maximum(indptr[rr],lo); sh=np.minimum(indptr[rr+1],hi)
    counts=np.clip((sh-sl),0,None).astype(np.int64); km=np.repeat(keep[rr],counts)
    kd=d[km]; ki=ii[km]; od[wp:wp+kd.size]=kd; oi[wp:wp+ki.size]=ki; wp+=kd.size
    log(" slab",si+1,"/",nsl,"wp",wp)
assert wp==tot,(wp,tot)
nip=np.zeros(idx.size+1,np.int64); nip[1:]=np.cumsum(rl)
X=sparse.csr_matrix((od,oi,nip),shape=(idx.size,ncols))
gene_name=f["var/gene_name"][:].astype(str)
ensg=f["var/gene_ids"][:].astype(str) if "gene_ids" in f["var"] else f["var/_index"][:].astype(str)
f.close()

clsi=cls[idx]
class_name=np.where(clsi==1,"target","NTC")
pgn_i=pgn[idx].copy(); pgn_i[clsi==2]="non-targeting"
obs=pd.DataFrame({"cell_class":class_name,"perturbed_gene_name":pgn_i,
    "guide_id":guide_id[idx],"guide_group":guide_group[idx],"guide_type":gt[idx],
    "lane_id":lane[idx],"n_genes_by_counts":ng[idx],"total_counts":tc[idx],
    "pct_counts_mt":mt[idx],"condition":COND})
A=ad.AnnData(X=X,obs=obs)
A.var["gene_name"]=gene_name; A.var["ensg"]=ensg
A.var_names=pd.Index(gene_name).astype(str); A.var_names_make_unique()
A.layers["counts"]=A.X.copy()
sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
A.write_h5ad(OUT_H5)
log("WROTE",OUT_H5,A.shape)

# per-target guide inventory: how many distinct guides & cells/guide (for dual-guide analysis)
inv={}
for t in TARGETS:
    m=(clsi==1)&(pgn_i==t)
    gids=guide_id[idx][m]
    uq,cnt=np.unique(gids,return_counts=True)
    inv[t]={"n_cells":int(m.sum()),"n_guides":int(len(uq)),
            "guides":{str(g):int(c) for g,c in zip(uq,cnt)}}
summ={"cond":COND,"shard":os.path.basename(SHARD),"n_cells":int(A.n_obs),
      "n_target":int((clsi==1).sum()),"n_ntc":int((clsi==2).sum()),
      "seed":SEED,"guide_inventory":inv}
json.dump(summ,open(OUT_SUM,"w"),indent=1)
log("WROTE",OUT_SUM)
print("___GATHER_DONE___")
