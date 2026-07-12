"""
Shared cell-level embedding pipeline for the CD4+ Perturb-seq atlas.

Builds a JOINT embedding across a set of (donor x condition) shards by:
  1) reading .obs from each raw-CSR h5ad via h5py (no full-matrix load),
  2) building a stratified keep-mask (all/ capped non-targeting controls +
     full coverage of nominated targets + a background targeting sample),
  3) streaming the CSR X once in row-chunks and keeping only masked rows,
  4) concatenating shards, then normalize -> log1p -> HVG -> PCA -> neighbors
     -> UMAP -> Leiden (Harmony over donor when >1 donor present),
  5) writing an embedded .h5ad checkpoint (lognorm X retained for gene coloring).

Standard scanpy on CPU. Condition/donor are preserved as biology; only donor is
integrated (Harmony) when multiple donors are present.
"""
import os, sys, time, json
import numpy as np
import h5py
import scipy.sparse as sp

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def _dec(x):
    return x.decode() if isinstance(x, (bytes, bytearray)) else x

def read_obs_cols(f, cols):
    """Read a set of .obs columns from an open h5ad h5py file.
    Returns dict col-> (numpy array of decoded values). Categorical columns are
    materialized to string arrays; numeric/bool returned as-is."""
    obs = f["obs"]
    out = {}
    for c in cols:
        if c not in obs:
            out[c] = None
            continue
        o = obs[c]
        if isinstance(o, h5py.Group) and "categories" in o and "codes" in o:
            cats = np.array([_dec(x) for x in o["categories"][:]], dtype=object)
            codes = o["codes"][:]
            vals = np.empty(codes.shape[0], dtype=object)
            good = codes >= 0
            vals[good] = cats[codes[good]]
            vals[~good] = "NA"
            out[c] = vals
        else:
            out[c] = o[:]
    return out

def build_keep_mask(obs, nominated_set, nt_cap, per_target_cap, bg_cap, seed):
    """Stratified subsample mask over one shard's cells.
    - drop low_quality True
    - non-targeting: random up to nt_cap
    - targeting & perturbed_gene_name in nominated_set: up to per_target_cap each
    - other targeting (background): random up to bg_cap
    """
    rng = np.random.default_rng(seed)
    n = obs["guide_type"].shape[0]
    gtype = obs["guide_type"]
    pgene = obs["perturbed_gene_name"]
    lq = obs.get("low_quality")
    ok = np.ones(n, dtype=bool)
    if lq is not None:
        ok &= ~lq.astype(bool)
    keep = np.zeros(n, dtype=bool)

    is_nt = (gtype == "non-targeting") & ok
    nt_idx = np.where(is_nt)[0]
    if len(nt_idx) > nt_cap:
        nt_idx = rng.choice(nt_idx, nt_cap, replace=False)
    keep[nt_idx] = True

    is_tgt = (gtype == "targeting") & ok
    nom_mask = is_tgt & np.isin(pgene, list(nominated_set))
    # per-target cap for nominated
    for g in np.unique(pgene[nom_mask]):
        gi = np.where(nom_mask & (pgene == g))[0]
        if len(gi) > per_target_cap:
            gi = rng.choice(gi, per_target_cap, replace=False)
        keep[gi] = True

    # background: targeting but not nominated
    bg_idx = np.where(is_tgt & ~np.isin(pgene, list(nominated_set)))[0]
    if len(bg_idx) > bg_cap:
        bg_idx = rng.choice(bg_idx, bg_cap, replace=False)
    keep[bg_idx] = True
    return keep

def read_masked_X(fn, keep, chunk=200_000):
    """Stream the CSR X of a raw h5ad once, keeping only rows where keep=True.
    Returns (scipy.csr matrix of kept rows in original order, kept_row_indices)."""
    f = h5py.File(fn, "r")
    Xg = f["X"]
    shape = tuple(int(x) for x in Xg.attrs["shape"])
    n_obs, n_var = shape
    indptr = Xg["indptr"][:]
    data_ds = Xg["data"]; ind_ds = Xg["indices"]
    keep_idx = np.where(keep)[0]
    blocks = []
    for s in range(0, n_obs, chunk):
        e = min(s + chunk, n_obs)
        local = keep_idx[(keep_idx >= s) & (keep_idx < e)]
        if len(local) == 0:
            continue
        d0, d1 = int(indptr[s]), int(indptr[e])
        data = data_ds[d0:d1]
        indices = ind_ds[d0:d1]
        lip = indptr[s:e+1] - d0
        blk = sp.csr_matrix((data, indices, lip), shape=(e - s, n_var))
        blocks.append(blk[local - s])
        log(f"    chunk {s:,}-{e:,}: kept {len(local):,} rows (cum {sum(b.shape[0] for b in blocks):,})")
    f.close()
    X = sp.vstack(blocks).tocsr()
    return X, keep_idx

CACHE_DIR = "/Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/checkpoints/shardcache"

def _cache_key(sh, params, nominated_set):
    return (f"{sh['donor']}_{sh['condition']}__nt{params['nt_cap']}_pt{params['per_target_cap']}"
            f"_bg{params['bg_cap']}_s{params['seed']}_nom{len(nominated_set)}")

def load_one_shard(sh, nominated_set, params, use_cache=True):
    """Read + subsample ONE shard, cache the subsampled AnnData to disk, reuse if present.
    Returns (AnnData, (var_index, gene_name))."""
    import anndata as ad, pandas as pd
    os.makedirs(CACHE_DIR, exist_ok=True)
    ckf = f"{CACHE_DIR}/{_cache_key(sh, params, nominated_set)}.sub.h5ad"
    if use_cache and os.path.exists(ckf):
        log(f"  shard {sh['donor']}_{sh['condition']}: CACHE HIT {ckf}")
        a = ad.read_h5ad(ckf)
        gn = a.var["gene_name"].values if "gene_name" in a.var else None
        return a, (a.var_names.values, gn)
    log(f"  shard {sh['donor']}_{sh['condition']}: reading obs (cache miss)")
    f = h5py.File(sh["path"], "r")
    obs_cols = ["guide_type","perturbed_gene_name","guide_id","lane_id",
                "n_genes_by_counts","pct_counts_mt","total_counts","low_quality"]
    obs = read_obs_cols(f, obs_cols)
    var = f["var"]; vidx = _dec(var.attrs.get("_index", b"_index"))
    var_index = np.array([_dec(x) for x in var[vidx][:]], dtype=object)
    gene_name = None
    if "gene_name" in var:
        gng = var["gene_name"]
        if isinstance(gng, h5py.Group) and "categories" in gng:
            cats = np.array([_dec(x) for x in gng["categories"][:]], dtype=object)
            gene_name = cats[gng["codes"][:]]
        else:
            gene_name = np.array([_dec(x) for x in gng[:]], dtype=object)
    f.close()
    keep = build_keep_mask(obs, nominated_set, params["nt_cap"],
                           params["per_target_cap"], params["bg_cap"], params["seed"])
    log(f"  shard {sh['donor']}_{sh['condition']}: keep {keep.sum():,}/{keep.shape[0]:,}")
    X, kidx = read_masked_X(sh["path"], keep, params["chunk"])
    sub = {k: (v[kidx] if v is not None else None) for k, v in obs.items()}
    obs_df = pd.DataFrame({
        "guide_type": sub["guide_type"], "perturbed_gene_name": sub["perturbed_gene_name"],
        "guide_id": sub["guide_id"], "lane_id": sub["lane_id"],
        "n_genes_by_counts": sub["n_genes_by_counts"], "pct_counts_mt": sub["pct_counts_mt"],
        "total_counts": sub["total_counts"]})
    obs_df["donor"] = sh["donor"]; obs_df["condition"] = sh["condition"]
    obs_df["is_nominated"] = np.isin(sub["perturbed_gene_name"], list(nominated_set))
    obs_df.index = [f"{sh['donor']}_{sh['condition']}_{i}" for i in kidx]
    a = ad.AnnData(X=X, obs=obs_df)
    a.var_names = var_index
    if gene_name is not None: a.var["gene_name"] = gene_name
    tmp = ckf + ".tmp"
    a.write_h5ad(tmp); os.replace(tmp, ckf)
    log(f"  shard {sh['donor']}_{sh['condition']}: cached -> {ckf} ({a.n_obs:,} cells)")
    return a, (var_index, gene_name)

def load_shards(shards, nominated_set, params, use_cache=True):
    """shards: list of dict(donor, condition, path). Returns AnnData (per-shard cached)."""
    import anndata as ad
    parts = []; var_ref = None
    for sh in shards:
        a, vr = load_one_shard(sh, nominated_set, params, use_cache=use_cache)
        if var_ref is None: var_ref = vr
        parts.append(a)
    log(f"  concatenating {len(parts)} shard(s)")
    adata = ad.concat(parts, join="outer", index_unique=None)
    vi, gn = var_ref
    adata.var_names = vi
    if gn is not None: adata.var["gene_name"] = gn
    return adata

def embed(adata, n_hvg=2000, n_pcs=50, n_neighbors=15, leiden_res=1.0,
          integrate_donor=False, seed=0):
    import scanpy as sc
    # This node exposes only ~10 usable CPUs (cgroup); numba requires threads <= that.
    # Cap conservatively so UMAP's numba parallel layer does not raise.
    navail = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else 8
    njobs = int(os.environ.get("EMBED_NJOBS", str(min(8, max(1, navail - 1)))))
    njobs = min(njobs, max(1, navail))
    try:
        import numba
        njobs = min(njobs, numba.config.NUMBA_NUM_THREADS)
    except Exception:
        pass
    log(f"  using n_jobs={njobs} (affinity={navail})")
    sc.settings.n_jobs = njobs
    log(f"  normalize_total + log1p  (n_obs={adata.n_obs:,}, n_var={adata.n_vars:,})")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()
    log("  highly_variable_genes")
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat", subset=False)
    log("  PCA")
    sc.pp.pca(adata, n_comps=n_pcs, use_highly_variable=True, zero_center=True, random_state=seed)
    rep = "X_pca"
    if integrate_donor and adata.obs["donor"].nunique() > 1:
        log("  Harmony integration over donor")
        import scanpy.external as sce
        sce.pp.harmony_integrate(adata, key="donor", basis="X_pca", adjusted_basis="X_pca_harmony")
        rep = "X_pca_harmony"
    log(f"  neighbors (rep={rep})")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=rep, random_state=seed)
    log("  UMAP")
    sc.tl.umap(adata, random_state=seed)
    log("  Leiden")
    sc.tl.leiden(adata, resolution=leiden_res, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
    return adata
