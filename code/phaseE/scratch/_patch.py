
import io, sys
p = "/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_scripts/per_target_umap.py"
src = open(p).read()
reader = r'''
def read_targeted_rows(fn, keep_mask, log_every=2000):
    """Read ONLY the CSR rows where keep_mask is True (targeted per-row slices).
    Avoids the full-array scan of read_masked_X; touches ~nnz(targets) not nnz(all)."""
    import h5py, numpy as np, scipy.sparse as sp, time
    kidx = np.nonzero(keep_mask)[0]
    with h5py.File(fn, "r", rdcc_nbytes=256*1024*1024, rdcc_nslots=100003) as f:
        X = f["X"]
        shp = X.attrs["shape"] if "shape" in X.attrs else None
        ncols = int(shp[1]) if shp is not None else int(f["var"][list(f["var"].keys())[0]].shape[0])
        t0=time.time(); indptr = X["indptr"][:]
        log(f"    indptr read {time.time()-t0:.0f}s ({indptr.shape[0]-1:,} rows)")
        data_ds = X["data"]; ind_ds = X["indices"]
        rows_data=[]; rows_ind=[]; new_indptr=[0]; total=0; t0=time.time()
        for j, r in enumerate(kidx):
            a=int(indptr[r]); b=int(indptr[r+1])
            if b>a:
                rows_data.append(data_ds[a:b]); rows_ind.append(ind_ds[a:b])
            total += (b-a); new_indptr.append(total)
            if (j+1) % log_every == 0:
                el=time.time()-t0
                log(f"    targeted {j+1:,}/{len(kidx):,} rows {total:,} nnz {el:.0f}s ({(j+1)/max(el,1e-9):.0f} rows/s)")
        data = np.concatenate(rows_data) if rows_data else np.empty(0, dtype=data_ds.dtype)
        indices = np.concatenate(rows_ind) if rows_ind else np.empty(0, dtype=ind_ds.dtype)
        Xc = sp.csr_matrix((data, indices, np.array(new_indptr, dtype=np.int64)), shape=(len(kidx), ncols))
        log(f"    built CSR {Xc.shape} nnz={Xc.nnz:,} in {time.time()-t0:.0f}s")
    return Xc, kidx
'''
# insert reader before 'def gather_union'
anchor = "def gather_union(targets):"
assert anchor in src, "anchor not found"
assert "def read_targeted_rows" not in src, "already patched"
src = src.replace(anchor, reader.strip("\n") + "\n\n" + anchor, 1)
# swap the read call
old = "X, kidx = read_masked_X(fn, keep, chunk=200_000)"
assert old in src, "read call not found"
src = src.replace(old, "X, kidx = read_targeted_rows(fn, keep)", 1)
open(p,"w").write(src)
print("patched OK; new length", len(src), "lines", src.count(chr(10)))
