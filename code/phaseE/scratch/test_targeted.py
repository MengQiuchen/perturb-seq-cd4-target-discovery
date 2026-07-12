
import h5py, numpy as np, time
def log(m): print("[%s] %s"%(time.strftime("%H:%M:%S"),m), flush=True)
fn="/Users/meng01/qiuchen/project/hackathon/perturb-seq/perturb-seq_data/cell_level/D4_Rest.assigned_guide.h5ad"
log("import done, opening")
f=h5py.File(fn,"r")
X=f["X"]; data=X["data"]; indices=X["indices"]; indptr=X["indptr"]
log("data chunks %s indices chunks %s indptr %d"%(data.chunks, indices.chunks, indptr.shape[0]))
t0=time.time(); ip=indptr[:]; log("indptr read %.1fs"%(time.time()-t0))
rng=np.random.default_rng(0)
rows=np.sort(rng.choice(ip.shape[0]-1, size=300, replace=False))
t0=time.time(); nnz=0
for r in rows:
    a,b=int(ip[r]),int(ip[r+1]); d=data[a:b]; ix=indices[a:b]; nnz+=d.shape[0]
dt=time.time()-t0
log("300 scattered rows: %.1fs nnz=%d per-row=%.1fms => est 12835 rows/shard = %.1f min"%(dt,nnz,1000*dt/300, 12835*(dt/300)/60))
f.close(); log("DONE")
