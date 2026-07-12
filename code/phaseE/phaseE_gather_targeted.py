#!/usr/bin/env python
# Phase E TARGETED gather for shards on a degraded OST (slow raw bandwidth).
# Instead of scanning all X nnz sequentially, read ONLY the kept rows' CSR spans,
# coalescing adjacent spans within COALESCE_NNZ to bound the number of Lustre RPCs.
# Same output contract as phaseE_gather.py.
# Usage: phaseE_gather_targeted.py <shard> <cond> <out_h5ad> <out_sum> <seed> [coalesce_nnz]
import sys, time, json, os, numpy as np, h5py
from scipy import sparse
import scanpy as sc, anndata as ad, pandas as pd

SHARD, COND, OUT_H5, OUT_SUM, SEED = sys.argv[1:6]
SEED=int(SEED)
COALESCE=int(sys.argv[6]) if len(sys.argv)>6 else 200_000  # merge gaps < this many nnz
t0=time.time()
def log(*a): print("[%7.1fs]"%(time.time()-t0),*a,flush=True)

TARGETS=["LAT","PLCG1","CD247","CD3E","ZAP70","VAV1","SMARCE1",
         "SGF29","MED24","MED12","TADA2B","NSD1","CHD4","SMARCB1",
         "STAT6","TRIP12","ARNT","SEL1L","SIK3"]
N_NTC=40000; N_BG=40000
rng=np.random.default_rng(SEED)

f=h5py.File(SHARD,"r",rdcc_nbytes=0); og=f["obs"]
def rc(n):
    g=og[n]
    if isinstance(g,h5py.Group):
        cats=g["categories"][:]; cats=cats.astype(str) if cats.dtype.kind in "SUO" else cats
        return np.asarray(cats)[g["codes"][:]]
    arr=g[:]; return arr.astype(str) if arr.dtype.kind in "SUO" else arr

gt=rc("guide_type"); lq=rc("low_quality"); pgn=rc("perturbed_gene_name"); lane=rc("lane_id")
ng=og["n_genes_by_counts"][:]; tc=og["total_counts"][:]; mt=og["pct_counts_mt"][:]
goodq = (~lq) if lq.dtype==bool else (lq.astype(str)!="True")
log("n_cells",gt.size,"goodq",int(goodq.sum()),"lq_dtype",str(lq.dtype))

is_tgt=(gt=="targeting")&goodq&np.isin(pgn,TARGETS)
is_ntc=(gt=="non-targeting")&goodq
is_bg=(gt=="targeting")&goodq&(~np.isin(pgn,TARGETS))
log("target-KD",int(is_tgt.sum()),"| NTC pool",int(is_ntc.sum()),"| bg pool",int(is_bg.sum()))
ntc_idx=np.where(is_ntc)[0]; bg_idx=np.where(is_bg)[0]
ntc_keep=rng.choice(ntc_idx,size=min(N_NTC,ntc_idx.size),replace=False)
bg_keep =rng.choice(bg_idx, size=min(N_BG, bg_idx.size), replace=False)
tgt_keep=np.where(is_tgt)[0]
keep=np.zeros(gt.size,bool); keep[tgt_keep]=True; keep[ntc_keep]=True; keep[bg_keep]=True
cls=np.zeros(gt.size,np.int8); cls[tgt_keep]=1; cls[ntc_keep]=2; cls[bg_keep]=3
idx=np.sort(np.where(keep)[0])
log("total kept",idx.size,"(t %d ntc %d bg %d)"%(tgt_keep.size,ntc_keep.size,bg_keep.size))

indptr=f["X/indptr"][:].astype(np.int64); nnz=int(indptr[-1])
ncols=int(f["X"].attrs["shape"][1]) if "shape" in f["X"].attrs else 18130
data_ds=f["X/data"]; ind_ds=f["X/indices"]

# ---- build coalesced read segments over kept rows ----
# each kept row r spans nnz [indptr[r], indptr[r+1]). Merge consecutive kept rows into a
# segment when the byte-gap between them is < COALESCE nnz. Read the whole segment span once,
# then slice out the kept rows within it.
starts=indptr[idx]; ends=indptr[idx+1]
segs=[]  # (row_lo_pos, row_hi_pos, nnz_lo, nnz_hi)  positions into idx[]
seg_lo=0; cur_nnz_hi=ends[0]
for p in range(1,idx.size):
    if starts[p]-cur_nnz_hi <= COALESCE:
        cur_nnz_hi=ends[p]
    else:
        segs.append((seg_lo,p-1,indptr[idx[seg_lo]],cur_nnz_hi)); seg_lo=p; cur_nnz_hi=ends[p]
segs.append((seg_lo,idx.size-1,indptr[idx[seg_lo]],cur_nnz_hi))
seg_span=sum(hi-lo for _,_,lo,hi in segs)
kept_nnz=int((ends-starts).sum())
log("segments",len(segs),"seg_span_nnz",seg_span,"kept_nnz",kept_nnz,"over-read x%.2f"%(seg_span/max(kept_nnz,1)))

# ---- read each segment, extract kept rows ----
rl=(ends-starts).astype(np.int64); tot=int(rl.sum())
od=np.empty(tot,np.float32); oi=np.empty(tot,np.int32); wp=0
report_every=max(1,len(segs)//20)
for si,(plo,phi,nlo,nhi) in enumerate(segs):
    dseg=data_ds[nlo:nhi]; iseg=ind_ds[nlo:nhi]
    for p in range(plo,phi+1):
        a=starts[p]-nlo; b=ends[p]-nlo
        n=b-a
        od[wp:wp+n]=dseg[a:b]; oi[wp:wp+n]=iseg[a:b]; wp+=n
    if si%report_every==0 or si==len(segs)-1:
        log(" seg",si+1,"/",len(segs),"wp",wp,"/",tot)
assert wp==tot,(wp,tot)
nip=np.zeros(idx.size+1,np.int64); nip[1:]=np.cumsum(rl)
X=sparse.csr_matrix((od,oi,nip),shape=(idx.size,ncols))
gene_name=f["var/gene_name"][:].astype(str)
ensg=f["var/gene_ids"][:].astype(str) if "gene_ids" in f["var"] else f["var/_index"][:].astype(str)
f.close()

clsi=cls[idx]
class_name=np.where(clsi==1,"target",np.where(clsi==2,"NTC","background"))
pgn_i=pgn[idx].copy(); pgn_i[clsi==2]="non-targeting"
obs=pd.DataFrame({"cell_class":class_name,"perturbed_gene_name":pgn_i,"guide_type":gt[idx],
    "lane_id":lane[idx],"n_genes_by_counts":ng[idx],"total_counts":tc[idx],
    "pct_counts_mt":mt[idx],"condition":COND})
A=ad.AnnData(X=X,obs=obs); A.var["gene_name"]=gene_name; A.var["ensg"]=ensg
A.var_names=pd.Index(gene_name).astype(str); A.var_names_make_unique()
log("built AnnData",A.shape)
A.layers["counts"]=A.X.copy(); sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
A.write_h5ad(OUT_H5); log("WROTE",OUT_H5)
per_tgt={t:int(np.sum((clsi==1)&(pgn_i==t))) for t in TARGETS}
summ={"cond":COND,"shard":os.path.basename(SHARD),"n_cells":int(A.n_obs),"n_genes":int(A.n_vars),
    "n_target_cells":int((clsi==1).sum()),"n_ntc":int((clsi==2).sum()),"n_bg":int((clsi==3).sum()),
    "per_target_KD":per_tgt,"seed":SEED,"method":"targeted","n_segments":len(segs),
    "overread":round(seg_span/max(kept_nnz,1),3)}
json.dump(summ,open(OUT_SUM,"w"),indent=1); log("WROTE",OUT_SUM)
print("___GATHER_DONE___",json.dumps(per_tgt))
