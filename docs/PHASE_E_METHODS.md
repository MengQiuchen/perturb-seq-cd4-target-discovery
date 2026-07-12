# Phase E methods reference (adapt the proven Phase D pipeline)

## Data facts (verified on cluster)
- Cell-level shards: `$PROJ/perturb-seq_data/cell_level/D{1..4}_{Rest,Stim8hr,Stim48hr}.assigned_guide.h5ad` (all 12 on Lustre, 111-161 GB each).
  - X = CSR, UNCOMPRESSED: `X/data` (float32), `X/indices` (int32/64), `X/indptr`. `X.attrs["shape"]`=(n_cells, 18130).
  - obs categorical cols: `perturbed_gene_name` (target symbol; "non-targeting" NOT used here — see guide_type), `perturbed_gene_id` (ENSG),
    `guide_type` (values: "targeting" / "non-targeting"), `guide_id`, `guide_group`, `lane_id`, `low_quality` (str "True"/"False").
  - obs numeric: `n_genes_by_counts`, `total_counts`, `pct_counts_mt`, `top_guide_UMI_counts`, `PuroR`.
  - var: `gene_name` (symbol), `gene_ids`/`_index` (ENSG). Use `var/gene_name` for canonical markers.
- DE_stats (pseudobulk reference): `$PROJ/perturb-seq_data/GWCD4i.DE_stats.h5ad`.
  - 33,983 (perturbation×condition) rows × 10,282 genes. obs uses key `index` (NOT `_index`).
  - obs `target_contrast_gene_name` (11,526 gene cats) × `culture_condition` (Rest/Stim8hr/Stim48hr) identify each row.
  - layers: `log_fc`, `zscore`, `p_value`, `adj_p_value`, `baseMean`, `lfcSE`. var has `gene_name`, `gene_ids`, `_index`.
  - obs QC/repro flags: `ontarget_significant` (bool), `ontarget_effect_size`, `guide_correlation_signif`, `donor_correlation_*`, `n_downstream`, `distal_offtarget_flag`, `neighboring_gene_KD`, `low_target_gex`.

## Phase D validation method (reproduce, extend to novel targets + all donors)
1. GATHER (per shard): keep cells = {all target-KD cells for the target list} ∪ {40k NTC sampled from guide_type=="non-targeting"} ∪ {40k background = random targeting cells for other genes}, all with low_quality!="True".
   - Use the slab-streaming CSR read in `gather_reference.py` (SLAB=500M nnz, per-nnz keep-mask via np.repeat(keep[rows],counts)). ~30-35s/slab, ~12 min/shard. DO NOT use backed .to_memory() on a scattered mask.
   - Cache each subset to a compact normalized h5ad + summary json (per-target cell counts).
2. NORMALIZE: layers["counts"]=X.copy(); sc.pp.normalize_total(target_sum=1e4); sc.pp.log1p.
3. ON-TARGET KD: per target, mean linear expression of the target gene in KD cells vs the 40k NTC pool; report % reduction + Mann-Whitney p. KD confirmed = >15% reduction.
4. CONCORDANCE: single-cell KD-vs-NTC log-FC per gene (log-normalized) vs DE_stats `log_fc` over the target's top-150 signature genes (by |zscore|). Pearson r. Concordant = r>0.30.
   - Flagship signatures already in `phaseD_outputs/flagship_pseudobulk_signatures.json`. For the 12 novel targets, regenerate the same structure from DE_stats (row = target_contrast_gene_name==sym & culture_condition==cond; rank genes by |zscore|, take top 150).
5. STIM-DEPENDENCE: mean |Δ| over the target's Stim8hr signature (KD vs NTC) in each condition. For E1 use Stim8hr vs (Rest, if gathered). For E2 (novel) compare Stim8hr vs Stim48hr, and vs Rest where the target's best_condition warrants. Stim-dependent = ratio>1.5.
6. STATE MANIFOLD (E3): build on control cells (NTC+background) of the D4_Stim8hr subset — 2000 HVG (seurat) -> scale -> PCA(50, arpack) -> neighbors(15,40pc) -> leiden(res1.0, igraph flavor) -> umap(min_dist0.3). Score 6 programs (below) with sc.tl.score_genes. Compare each target's KD-vs-NTC program-score shift (Δ). Then CLUSTER the per-target Δ-vectors (all 19 targets) and test whether the chromatin/TF axis separates from the TCR-proximal naive/memory-reversion cluster (e.g. hierarchical clustering + silhouette, or a simple 2-group contrast).

## Program gene sets (verbatim from Phase D)
PROGRAMS = {
 "naive_memory":["CCR7","SELL","TCF7","LEF1","IL7R"],
 "effector":["GZMB","GZMA","PRF1","IFNG","NKG7","GNLY","CCL5"],
 "activation":["IL2RA","CD69","TNFRSF9","MKI67","TNFRSF18","ICOS"],
 "treg":["FOXP3","IKZF2","CTLA4","IL2RA","TNFRSF18"],
 "cytokine":["IL2","IFNG","TNF","IL4","IL13","IL21","CSF2"],
 "exhaustion":["PDCD1","LAG3","HAVCR2","TIGIT","TOX","ENTPD1"],
}

## Operational (host)
- Env: source /mnt/scratche/slow/ghlab/qiuchen/miniforge3/etc/profile.d/conda.sh && conda activate /mnt/scratche/slow/ghlab/qiuchen/miniforge3/envs/perturb-seq
- Direct call_command works ONLY while the user's Slurm allocation adopts the node (pam_slurm_adopt). login_shell=False always. python3 not on PATH without conda.
- Gathers are long: launch as detached `nohup setsid ... > log 2>&1 &` writing a `.done` marker; poll the marker. Keep <=2 concurrent gathers (I/O saturation risk). NO c.download (SFTP wedges) — pull small results via base64, tar large ones.
- perturb-seq env has NO pyarrow -> write CSV not parquet.
- Write everything under $PROJ (Lustre symlink); large subsets under perturb-seq_data or phaseE_outputs/checkpoints.
