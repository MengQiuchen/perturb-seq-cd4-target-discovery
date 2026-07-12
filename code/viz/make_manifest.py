"""Build viz_outputs/FIGURES_README.md cataloguing every figure produced.
Reads per-wave manifest JSONs + the figures dir; writes a human-readable manifest.
Usage: python make_manifest.py
"""
import os, json, glob, time
PROJ = "/Users/meng01/qiuchen/project/hackathon/perturb-seq"
OUT = f"{PROJ}/viz_outputs"; FIG=f"{OUT}/figures"; CK=f"{OUT}/checkpoints"

FIG_DESC = {
    "A1_condition_facets": "Joint UMAP faceted by condition; each panel highlights one condition's cells over a grey all-cell background. Shows how the T-cell state manifold shifts Rest→Stim.",
    "A2_condition_overlay": "Single joint UMAP with all cells colored by condition. Direct overlay of the condition axis.",
    "A3_activation_markers": "Joint UMAP colored by log-norm expression of canonical activation markers (IL2RA, CD69, MKI67, TNFRSF9, IFNG, SELL). Annotates the resting→activated direction.",
    "B_within_condition_targeting": "Per condition, non-targeting control cells (alarm hue) overlaid on targeting/perturbed cells (grey). Within-condition comparison of perturbed vs control.",
    "C_cross_donor": "Joint UMAP faceted by donor (multi-donor embeddings only); each panel highlights one donor. Shows cross-donor reproducibility after Harmony integration.",
    "D_per_target": "Grid of per-target highlight panels: for each nominated target gene, cells carrying a knockdown of that gene are highlighted over the full-embedding grey background.",
    "E_per_module": "Grid of per-group highlight panels: nominated targets grouped by co-regulation module (Phase B4) or therapeutic direction; all cells in the group highlighted together (reduces plot count vs one-per-gene).",
}

def desc_for(fn):
    for k, v in FIG_DESC.items():
        if fn.startswith(k): return v
    return "—"

def main():
    manifests = {}
    for mf in sorted(glob.glob(f"{CK}/*_manifest.json")) + sorted(glob.glob(f"{CK}/wave*_manifest.json")):
        try: manifests[os.path.basename(mf)] = json.load(open(mf))
        except Exception: pass
    figs = sorted(os.path.basename(p) for p in glob.glob(f"{FIG}/*.png"))

    lines = []
    lines.append("# Cell-level UMAP visualizations — figure manifest\n")
    lines.append(f"_Generated {time.strftime('%Y-%m-%d %H:%M')} on ssh:clust1-rocm-4._\n")
    lines.append("Visualizations of the Marson–Pritchard genome-scale CRISPRi Perturb-seq atlas "
                 "in primary human CD4⁺ T cells (cell-level `*.assigned_guide.h5ad` shards).\n")
    lines.append("## Method\n")
    lines.append("- **Subsampling** (stratified, per shard): all non-targeting controls (capped), "
                 "full coverage of the 301 nominated scorecard targets (per-target cap), plus a "
                 "background sample of other targeting cells. Low-quality cells dropped.\n")
    lines.append("- **Embedding**: raw CSR counts streamed via h5py (no full-matrix load) → "
                 "`normalize_total(1e4)` → `log1p` → 2000 HVG → PCA(50) → neighbors(15) → UMAP → "
                 "Leiden. Multi-donor embeddings add **Harmony integration over donor** "
                 "(condition preserved as biology). scanpy 1.11.5, CPU.\n")
    lines.append("- **Figure style**: embedding scatters drop ticks (corner axis arrows instead), "
                 "CVD-safe threaded palette, grey background context, rasterized points.\n")
    lines.append("\n## Waves / embeddings\n")
    for mfname, m in manifests.items():
        tag = m.get("tag", mfname)
        lines.append(f"### `{tag}`")
        donors = ",".join(m.get("donors", [])) if m.get("donors") else "D4"
        conds = ",".join(m.get("conditions", [])) if m.get("conditions") else "?"
        lines.append(f"- cells embedded: **{m.get('n_cells','?'):,}**" if isinstance(m.get('n_cells'), int) else f"- cells embedded: {m.get('n_cells','?')}")
        lines.append(f"- donors: {donors} · conditions: {conds}")
        if m.get("modules"): lines.append(f"- module groups: {m['modules']}")
        if m.get("top_targets"): lines.append(f"- per-target panels: {', '.join(m['top_targets'][:25])}")
        lines.append(f"- checkpoint: `checkpoints/{tag}.embedded.h5ad`  (+ `{tag}.obs_umap.parquet`)")
        lines.append("")
    lines.append("## All figures\n")
    lines.append("| file | description |")
    lines.append("|---|---|")
    for fn in figs:
        lines.append(f"| `figures/{fn}` | {desc_for(fn)} |")
    lines.append("\n## Code (the plotting 'track')\n")
    for cf in sorted(os.path.basename(p) for p in glob.glob(f"{OUT}/code/*.py")):
        lines.append(f"- `code/{cf}`")
    lines.append("\n### Config files\n")
    for cf in sorted(os.path.basename(p) for p in glob.glob(f"{OUT}/code/*_config.json")):
        lines.append(f"- `code/{cf}`")
    open(f"{OUT}/FIGURES_README.md", "w").write("\n".join(lines))
    print(f"wrote {OUT}/FIGURES_README.md with {len(figs)} figures, {len(manifests)} manifests")

if __name__ == "__main__":
    main()
