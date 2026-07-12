
import h5py, numpy as np
f = h5py.File("/Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/checkpoints/shardcache/D4_Rest__nt25000_pt400_bg90000_s0_nom301.sub.h5ad", "r")
print("keys:", list(f.keys()))
print("layers:", list(f["layers"].keys()) if "layers" in f else "none")
X = f["X"]
if hasattr(X, "keys"):
    print("X kind: group/CSR; subkeys", list(X.keys()))
    d = X["data"][:300]
else:
    print("X kind: dense")
    d = np.asarray(X[:5]).ravel()[:300]
print("sample data[:12]:", np.round(d[:12], 3).tolist())
print("integer_like:", bool(np.allclose(d, np.round(d))))
print("min/max(300):", float(d.min()), float(d.max()))
print("obs cols:", list(f["obs"].keys())[:20])
f.close()
