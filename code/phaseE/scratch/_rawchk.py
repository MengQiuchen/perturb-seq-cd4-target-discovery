
import h5py, numpy as np
for D in ["D1","D4"]:
    f=h5py.File(f"/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/checkpoints/{D}_Stim8hr.subset.h5ad","r")
    X=f["X"]
    d=X["data"][:50]
    print(f"{D}: X[data][:20] =", np.round(d[:20],4).tolist())
    print(f"    all_integer(first 5000)?", bool(np.allclose(X["data"][:5000], np.round(X["data"][:5000]))))
    print(f"    min={float(X['data'][:5000].min()):.4f} max={float(X['data'][:5000].max()):.4f}")
    # check layers
    print(f"    layers:", list(f["layers"].keys()) if "layers" in f else "none")
    # row sums for a couple cells (raw counts -> integer totals ~ thousands; lognorm -> ~ hundreds non-integer)
    ip=X["indptr"][:6]
    for r in range(3):
        a,b=int(ip[r]),int(ip[r+1]); rs=float(X["data"][a:b].sum())
        print(f"    cell{r} nnz={b-a} sum={rs:.3f}")
    f.close()
