
import h5py, time, numpy as np
fn="/Users/meng01/qiuchen/project/hackathon/perturb-seq/perturb-seq_data/cell_level/D1_Rest.assigned_guide.h5ad"
f=h5py.File(fn,"r")
X=f["X"]; d=X["data"]
n=d.shape[0]; print("nnz=",n)
# read a 200M-element slab from a fresh offset (mimics 200k-row chunk contiguous read)
off=n//2
t=time.time(); a=d[off:off+50_000_000]; dt=time.time()-t
mb=a.nbytes/1e6
print(f"h5py slab 50M float32 @mid = {mb:.0f}MB in {dt:.1f}s = {mb/dt:.0f} MB/s")
f.close()
