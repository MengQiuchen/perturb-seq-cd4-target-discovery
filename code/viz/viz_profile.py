
import h5py, sys, numpy as np
for fn in sys.argv[1:]:
    f = h5py.File(fn, "r")
    name = fn.split("/")[-1]
    obs = f["obs"]
    n = obs["_index"].shape[0]
    print(f"\n=== {name}: n_obs={n:,} ===")
    Xg = f["X"]; shp = Xg.attrs.get("shape", None); nnz = Xg["data"].shape[0]
    print(f"  X shape={list(shp)}, nnz={nnz:,}, density={nnz/(shp[0]*shp[1]):.4f}")
    # guide_type
    gt = obs["guide_type"]
    cats = [x.decode() if isinstance(x,bytes) else x for x in gt["categories"][:]]
    codes = gt["codes"][:]
    for i,ct in enumerate(cats):
        print(f"  guide_type '{ct}': {(codes==i).sum():,}")
    # low_quality
    if "low_quality" in obs:
        lq = obs["low_quality"][:]
        print(f"  low_quality True: {lq.sum():,} ({100*lq.mean():.1f}%)")
    # cells per target (perturbed_gene_name)
    pg = obs["perturbed_gene_name"]
    pcats = pg["categories"][:]
    pcodes = pg["codes"][:]
    counts = np.bincount(pcodes[pcodes>=0], minlength=len(pcats))
    print(f"  n targets (perturbed_gene_name): {len(pcats):,}")
    print(f"  cells/target: min={counts.min()}, median={int(np.median(counts))}, "
          f"mean={counts.mean():.0f}, max={counts.max():,}")
    for thr in [50,100,200,500,1000]:
        print(f"    targets with >={thr} cells: {(counts>=thr).sum():,}")
    # guides per... guide_id count
    gid = obs["guide_id"]; print(f"  n guide_id: {len(gid['categories']):,}")
    f.close()
