
import h5py
P="/Users/meng01/qiuchen/project/hackathon/perturb-seq/perturb-seq_data/GWCD4i.DE_stats.h5ad"
f=h5py.File(P,"r")
for L in ["zscore","log_fc"]:
    d=f["layers/"+L]; print(L,"shape",d.shape,"chunks",d.chunks,"dtype",d.dtype)
f.close()
