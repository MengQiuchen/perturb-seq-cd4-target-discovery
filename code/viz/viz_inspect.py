
import h5py, sys, numpy as np
fn = sys.argv[1]
f = h5py.File(fn, "r")

def kind(o):
    if isinstance(o, h5py.Group):
        enc = o.attrs.get("encoding-type", b"")
        enc = enc.decode() if isinstance(enc, bytes) else enc
        return f"group[{enc}]" if enc else "group"
    return f"dataset {o.shape} {o.dtype}"

print("TOP-LEVEL KEYS:", list(f.keys()))
# X
if "X" in f:
    X = f["X"]
    if isinstance(X, h5py.Group):
        et = X.attrs.get("encoding-type", b"")
        et = et.decode() if isinstance(et, bytes) else et
        shp = X.attrs.get("shape", None)
        print(f"\nX: {et}, shape attr={shp}, members={list(X.keys())}")
        if "data" in X: print("   X/data:", X["data"].shape, X["data"].dtype)
    else:
        print(f"\nX: dense {X.shape} {X.dtype}")

# obs
if "obs" in f:
    obs = f["obs"]
    cols = list(obs.keys())
    print(f"\nOBS columns ({len(cols)}):", cols)
    idx = obs.attrs.get("_index", b"_index")
    idx = idx.decode() if isinstance(idx, bytes) else idx
    if idx in obs:
        n = obs[idx].shape[0]
        print("  n_obs:", n)
    # For each obs column, if categorical, report n categories + sample
    for col in cols:
        o = obs[col]
        if isinstance(o, h5py.Group) and "categories" in o:
            cats = o["categories"][:]
            cats = [x.decode() if isinstance(x, bytes) else x for x in cats]
            print(f"  [cat] {col}: {len(cats)} categories; e.g. {cats[:8]}")
        elif isinstance(o, h5py.Dataset):
            print(f"  [num] {col}: {o.shape} {o.dtype}")

# var
if "var" in f:
    var = f["var"]
    idx = var.attrs.get("_index", b"_index")
    idx = idx.decode() if isinstance(idx, bytes) else idx
    if idx in var:
        genes = var[idx][:20]
        genes = [x.decode() if isinstance(x, bytes) else x for x in genes]
        nvar = var[idx].shape[0]
        print(f"\nVAR: n_var={nvar}, index='{idx}', first genes={genes[:10]}")
    print("  var columns:", list(var.keys()))

# obsm / layers / uns
for grp in ["obsm", "layers", "varm", "obsp", "uns"]:
    if grp in f:
        g = f[grp]
        items = {}
        for k in g.keys():
            o = g[k]
            items[k] = (o.shape if isinstance(o, h5py.Dataset) else kind(o))
        print(f"\n{grp.upper()}:", items)
f.close()
