
import h5py, numpy as np, time
def log(m): print("[%s] %s"%(time.strftime("%H:%M:%S"),m),flush=True)
fn="/Users/meng01/qiuchen/project/hackathon/perturb-seq/perturb-seq_data/cell_level/D4_Rest.assigned_guide.h5ad"
f=h5py.File(fn,"r", rdcc_nbytes=256*1024*1024, rdcc_nslots=100003)
X=f["X"]; data=X["data"]
log("opened; data shape %d chunks %s"%(data.shape[0], data.chunks))
# (1) large contiguous slice 40MB
t0=time.time(); a=data[0:10_000_000]; dt=time.time()-t0
log("40MB contiguous slice: %.1fs = %.1f MB/s"%(dt, 40/dt))
# (2) another fresh offset
t0=time.time(); b=data[500_000_000:510_000_000]; dt=time.time()-t0
log("40MB @2GB offset: %.1fs = %.1f MB/s"%(dt, 40/dt))
# (3) indptr full read
t0=time.time(); ip=X["indptr"][:]; dt=time.time()-t0
log("indptr 21MB: %.1fs = %.1f MB/s"%(dt, 21/dt))
f.close(); log("DONE")
