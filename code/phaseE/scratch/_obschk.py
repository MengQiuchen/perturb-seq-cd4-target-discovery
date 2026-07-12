
import h5py, numpy as np
fn="/Users/meng01/qiuchen/project/hackathon/perturb-seq/perturb-seq_data/cell_level/D1_Rest.assigned_guide.h5ad"
f=h5py.File(fn,"r")
og=f["obs"]
print("obs keys:", list(og.keys()))
for k in ["guide_type","low_quality","perturbed_gene_name","lane_id","n_genes_by_counts","total_counts","pct_counts_mt"]:
    if k in og:
        g=og[k]
        if isinstance(g,h5py.Group):
            cats=g["categories"][:]
            print(f"  {k}: categorical, {len(cats)} cats, e.g. {cats[:4].astype(str).tolist() if cats.dtype.kind in 'SUO' else cats[:4].tolist()}")
        else:
            print(f"  {k}: array dtype={g.dtype}, shape={g.shape}")
    else:
        print(f"  {k}: *** MISSING ***")
print("var keys:", list(f["var"].keys())[:8])
print("X shape attr:", dict(f["X"].attrs) if hasattr(f["X"],"attrs") else "n/a")
f.close()
