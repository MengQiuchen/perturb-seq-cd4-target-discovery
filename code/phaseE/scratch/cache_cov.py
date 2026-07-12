
import h5py, numpy as np, glob, json, time, collections
def dec(x): return x.decode() if isinstance(x,(bytes,bytearray)) else x
def log(m): print("[%s] %s"%(time.strftime("%H:%M:%S"),m),flush=True)
CACHE="/Users/meng01/qiuchen/project/hackathon/perturb-seq/viz_outputs/checkpoints/shardcache"
targets=["ZNF613","TRAF3","UBXN1","TADA2B","ATAD5","CBLB","CLCC1","ATF7IP2","LRP2BP","NFKBIL1","EDC3","PARP14","ELOF1","CYP20A1","CASP3","PIGB"]
tset=set(targets)
agg=collections.defaultdict(lambda: collections.defaultdict(int))
for fn in sorted(glob.glob(CACHE+"/D4_*.sub.h5ad")):
    cond=fn.split("/")[-1].split("__")[0].replace("D4_","")
    f=h5py.File(fn,"r"); obs=f["obs"]
    o=obs["perturbed_gene_name"]
    if isinstance(o,h5py.Group) and "categories" in o:
        cats=np.array([dec(x) for x in o["categories"][:]],dtype=object); codes=o["codes"][:]
        v=np.empty(codes.shape[0],dtype=object); g=codes>=0; v[g]=cats[codes[g]]; v[~g]="NA"
    else:
        v=np.array([dec(x) for x in o[:]],dtype=object)
    cnt=collections.Counter(v)
    for t in targets: agg[t][cond]=int(cnt.get(t,0))
    f.close()
    log(f"{cond}: scanned {len(v):,} cells")
print("GENE\tRest\tStim8hr\tStim48hr\tTOTAL_in_cache", flush=True)
for t in targets:
    r=agg[t]; tot=r["Rest"]+r["Stim8hr"]+r["Stim48hr"]
    print(f"{t}\t{r['Rest']}\t{r['Stim8hr']}\t{r['Stim48hr']}\t{tot}", flush=True)
json.dump({t:dict(agg[t]) for t in targets}, open("/Users/meng01/qiuchen/project/hackathon/perturb-seq/phaseE_outputs/umap_within/cache_target_coverage.json","w"), indent=2)
log("DONE")
