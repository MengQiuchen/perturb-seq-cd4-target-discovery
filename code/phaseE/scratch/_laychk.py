
import h5py, numpy as np
for D in ["D1","D4"]:
    f=h5py.File(f"/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/checkpoints/{D}_Stim8hr.subset.h5ad","r")
    L=f["layers"]["counts"]
    kind = "group/CSR" if hasattr(L,"keys") else "dense"
    print(f"{D}: counts layer kind={kind}", list(L.keys()) if hasattr(L,"keys") else "")
    d=L["data"][:20] if hasattr(L,"keys") else np.asarray(L[:3]).ravel()[:20]
    print(f"    counts[:20] =", np.round(d,3).tolist())
    dd = L["data"][:5000] if hasattr(L,"keys") else np.asarray(L[:50]).ravel()
    print(f"    all_integer? {bool(np.allclose(dd, np.round(dd)))} min={float(dd.min())} max={float(dd.max())}")
    f.close()
