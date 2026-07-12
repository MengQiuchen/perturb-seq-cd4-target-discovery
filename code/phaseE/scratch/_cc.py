
import h5py, numpy as np, collections
def dec(x): return x.decode() if isinstance(x,(bytes,bytearray)) else x
for D in ["D1","D2","D3","D4"]:
    f=h5py.File(f"/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/checkpoints/{D}_Stim8hr.subset.h5ad","r"); obs=f["obs"]
    o=obs["cell_class"]
    if isinstance(o,h5py.Group) and "categories" in o:
        cats=np.array([dec(x) for x in o["categories"][:]],dtype=object); codes=o["codes"][:]
        v=np.empty(codes.shape[0],dtype=object); g=codes>=0; v[g]=cats[codes[g]]; v[~g]="NA"
    else:
        v=np.array([dec(x) for x in o[:]],dtype=object)
    print(D, "cell_class:", dict(collections.Counter(v)))
    f.close()
