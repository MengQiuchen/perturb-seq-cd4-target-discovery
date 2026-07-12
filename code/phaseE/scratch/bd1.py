
import h5py,time,numpy as np
P="/Users/meng01/qiuchen/project/hackathon/perturb-seq/perturb-seq_data/cell_level/D1_Stim8hr.assigned_guide.h5ad"
f=h5py.File(P,"r",rdcc_nbytes=0)
n=f["X/data"].shape[0]; print("D1 X/data n=",n)
for lbl,off in [("start",0),("mid",n//2),("late",int(n*0.8))]:
    hi=min(off+50_000_000,n)
    t=time.time(); d=f["X/data"][off:hi]; dt=time.time()-t
    print("D1 data %-5s off=%d read %.2fGB in %.1fs = %.1f MB/s"%(lbl,off,d.nbytes/1e9,dt,d.nbytes/1e6/dt))
f.close()
