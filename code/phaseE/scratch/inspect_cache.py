
import h5py, numpy as np, glob, collections
def dec(x): return x.decode() if isinstance(x,(bytes,bytearray)) else x
CACHE="/Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/checkpoints/shardcache"
for fn in sorted(glob.glob(CACHE+"/D4_*.sub.h5ad")):
    f=h5py.File(fn,"r"); obs=f["obs"]
    n=obs["_index"].shape[0]
    def cat(col):
        o=obs[col]
        if isinstance(o,h5py.Group) and "categories" in o:
            cats=np.array([dec(x) for x in o["categories"][:]],dtype=object); codes=o["codes"][:]
            v=np.empty(codes.shape[0],dtype=object); g=codes>=0; v[g]=cats[codes[g]]; v[~g]="NA"; return v
        return np.array([dec(x) for x in o[:]],dtype=object)
    gt=cat("guide_type"); cond=cat("condition")
    print(fn.split("/")[-1], "n="+str(n),
          "cond:",dict(collections.Counter(cond)),
          "gtype:",dict(collections.Counter(gt)), flush=True)
    var=f["var"]; vi=dec(var.attrs.get("_index",b"_index"))
    print("   var keys:", list(var.keys())[:6], "n_var=", var[vi].shape[0], flush=True)
    f.close()
