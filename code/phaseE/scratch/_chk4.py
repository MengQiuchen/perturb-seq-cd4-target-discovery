
import h5py, numpy as np, collections
def dec(x): return x.decode() if isinstance(x,(bytes,bytearray)) else x
fn=f"/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/checkpoints/D4_Stim8hr.subset.h5ad"
f=h5py.File(fn,"r"); obs=f["obs"]; n=obs["_index"].shape[0]
def cat(col):
    o=obs[col]
    if isinstance(o,h5py.Group) and "categories" in o:
        cats=np.array([dec(x) for x in o["categories"][:]],dtype=object); codes=o["codes"][:]
        v=np.empty(codes.shape[0],dtype=object); g=codes>=0; v[g]=cats[codes[g]]; v[~g]="NA"; return v
    return np.array([dec(x) for x in o[:]],dtype=object)
gt=cat("guide_type"); cond=cat("condition"); pg=cat("perturbed_gene_name")
print("D4_Stim8hr.subset: n=%d nvar=%d"%(n, f["var"]["_index"].shape[0]))
print("obs cols:", list(obs.keys()))
print("guide_type:", dict(collections.Counter(gt)))
print("condition:", dict(collections.Counter(cond)))
c=collections.Counter(pg); print("n_genes=%d top8=%s"%(len(c), c.most_common(8)))
# var index kind — symbols or ensembl?
vi=np.array([dec(x) for x in f["var"]["_index"][:8]],dtype=object)
print("var_index[:8]:", list(vi))
print("var keys:", list(f["var"].keys()))
f.close()
