
import h5py, numpy as np, collections
def dec(x): return x.decode() if isinstance(x,(bytes,bytearray)) else x
for D in ["D1","D2","D3"]:
    fn=f"/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/checkpoints/{D}_Stim8hr.subset.h5ad"
    import os
    if not os.path.exists(fn):
        print(D, "MISSING"); continue
    f=h5py.File(fn,"r"); obs=f["obs"]
    n=obs["_index"].shape[0]
    def cat(col):
        if col not in obs: return None
        o=obs[col]
        if isinstance(o,h5py.Group) and "categories" in o:
            cats=np.array([dec(x) for x in o["categories"][:]],dtype=object); codes=o["codes"][:]
            v=np.empty(codes.shape[0],dtype=object); g=codes>=0; v[g]=cats[codes[g]]; v[~g]="NA"; return v
        return np.array([dec(x) for x in o[:]],dtype=object)
    gt=cat("guide_type"); cond=cat("condition"); pg=cat("perturbed_gene_name")
    X=f["X"]; xkind="CSR" if hasattr(X,"keys") else "dense"
    nvar=f["var"]["_index"].shape[0] if "_index" in f["var"] else "?"
    print(f"{D}_Stim8hr: n={n:,} nvar={nvar} Xkind={xkind}")
    print(f"   obs cols: {list(obs.keys())[:14]}")
    if gt is not None: print(f"   guide_type: {dict(collections.Counter(gt))}")
    if cond is not None: print(f"   condition: {dict(collections.Counter(cond))}")
    if pg is not None:
        c=collections.Counter(pg); print(f"   n_genes={len(c)}, top5={c.most_common(5)}")
    f.close()
