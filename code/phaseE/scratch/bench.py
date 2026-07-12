
import h5py,time,numpy as np
P="/Users/meng01/qiuchen/project/hackathon/perturb-seq/perturb-seq_data/cell_level/D3_Stim8hr.assigned_guide.h5ad"
f=h5py.File(P,"r",rdcc_nbytes=0)
print("data dtype",f["X/data"].dtype,"chunks",f["X/data"].chunks,"shape",f["X/data"].shape)
print("indices dtype",f["X/indices"].dtype,"chunks",f["X/indices"].chunks)
t=time.time(); d=f["X/data"][0:100_000_000]; dt=time.time()-t
print("read 100M float32 (%.2f GB) in %.1fs = %.1f MB/s"%(d.nbytes/1e9,dt,d.nbytes/1e6/dt))
t=time.time(); ii=f["X/indices"][0:100_000_000]; dt=time.time()-t
print("read 100M indices (%.2f GB) in %.1fs = %.1f MB/s"%(ii.nbytes/1e9,dt,ii.nbytes/1e6/dt))
f.close()
