
import h5py, numpy as np, time
P="/Users/meng01/qiuchen/project/hackathon/perturb-seq/perturb-seq_data/cell_level/D1_Stim8hr.assigned_guide.h5ad"
TARGETS=["LAT","PLCG1","CD247","CD3E","ZAP70","VAV1","SMARCE1","SGF29","MED24","MED12","TADA2B","NSD1","CHD4","SMARCB1","STAT6","TRIP12","ARNT","SEL1L","SIK3"]
rng=np.random.default_rng(1)
t=time.time()
f=h5py.File(P,"r",rdcc_nbytes=0); og=f["obs"]
def rc(n):
    g=og[n]
    if isinstance(g,h5py.Group):
        cats=g["categories"][:]; cats=cats.astype(str) if cats.dtype.kind in "SUO" else cats
        return np.asarray(cats)[g["codes"][:]]
    a=g[:]; return a.astype(str) if a.dtype.kind in "SUO" else a
gt=rc("guide_type"); lq=rc("low_quality"); pgn=rc("perturbed_gene_name")
goodq=(~lq) if lq.dtype==bool else (lq.astype(str)!="True")
is_tgt=(gt=="targeting")&goodq&np.isin(pgn,TARGETS)
is_ntc=(gt=="non-targeting")&goodq
is_bg=(gt=="targeting")&goodq&(~np.isin(pgn,TARGETS))
ntc=rng.choice(np.where(is_ntc)[0],40000,replace=False)
bg=rng.choice(np.where(is_bg)[0],40000,replace=False)
keep=np.zeros(gt.size,bool); keep[np.where(is_tgt)[0]]=True; keep[ntc]=True; keep[bg]=True
idx=np.sort(np.where(keep)[0])
print("obs read %.1fs, kept=%d of %d (%.2f%%)"%(time.time()-t,idx.size,gt.size,100*idx.size/gt.size))
t=time.time(); indptr=f["X/indptr"][:].astype(np.int64); print("indptr read %.1fs"%(time.time()-t))
nnz=int(indptr[-1]); starts=indptr[idx]; ends=indptr[idx+1]
kept_nnz=int((ends-starts).sum())
print("total nnz=%d, kept_nnz=%d (%.2f%%)"%(nnz,kept_nnz,100*kept_nnz/nnz))
# data bytes: float32 (4B) + indices int64 (8B) = 12 B/nnz
print("kept data+idx = %.2f GB"%(kept_nnz*12/1e9))
print("full data+idx = %.2f GB"%(nnz*12/1e9))
for COAL in [0,10_000,50_000,200_000,1_000_000,5_000_000]:
    segs=0; span=0; seg_lo_nnz=starts[0]; cur_hi=ends[0]
    for p in range(1,idx.size):
        if starts[p]-cur_hi<=COAL: cur_hi=ends[p]
        else: span+=cur_hi-seg_lo_nnz; segs+=1; seg_lo_nnz=starts[p]; cur_hi=ends[p]
    span+=cur_hi-seg_lo_nnz; segs+=1
    print("COAL=%9d -> segs=%7d span_nnz=%12d overread=%.2fx read=%.1fGB"%(COAL,segs,span,span/kept_nnz,span*12/1e9))
f.close()
