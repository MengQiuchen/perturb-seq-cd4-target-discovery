
import h5py, numpy as np
f = h5py.File("/Users/meng01/qiuchen/project/hackathon/perturb-seq/perturb-seq_data/GWCD4i.DE_stats.h5ad","r")

def cat(col):
    g = f[f"obs/{col}"]
    codes = g["codes"][:]
    cats  = g["categories"][:].astype(str)
    return codes, cats

tg_codes, tg_cats = cat("target_contrast_gene_name")
cc_codes, cc_cats = cat("culture_condition")
oes = f["obs/ontarget_effect_size"][:].astype(float)

# var gene_name
vg = f["var/gene_name"]
if isinstance(vg, h5py.Group):
    gcodes = vg["codes"][:]; gcats = vg["categories"][:].astype(str)
    gene_name = gcats[gcodes]
else:
    gene_name = vg[:].astype(str)
gcol = {g:i for i,g in enumerate(gene_name)}

print("layers:", list(f["layers"].keys()))
print("oes range: min=%.2f med=%.2f max=%.2f" % (oes.min(), np.median(oes), oes.max()))

# pick rows spanning the range whose target gene is measured
order = np.argsort(oes)
picks = list(order[:2]) + list(order[len(order)//2-1:len(order)//2+1]) + list(order[-2:])
print("\nrow | target | cond | ontarget_effect_size | log_fc[self] | zscore[self] | adj_p[self]")
for r in picks:
    tgt = tg_cats[tg_codes[r]]
    j = gcol.get(tgt, None)
    lfc = f["layers/log_fc"][r, j] if j is not None else np.nan
    z   = f["layers/zscore"][r, j] if j is not None else np.nan
    ap  = f["layers/adj_p_value"][r, j] if (j is not None and "adj_p_value" in f["layers"]) else np.nan
    print("%6d | %-10s | %-8s | %10.3f | %10.3f | %10.3f | %8.2e" % (r, tgt, cc_cats[cc_codes[r]], oes[r], lfc, z, ap))
f.close()
