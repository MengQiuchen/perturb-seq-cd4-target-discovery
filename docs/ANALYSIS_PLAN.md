# Analysis Plan — Novel Drug Targets from CD4⁺ T Cell Genome-Scale Perturb-seq

**Project goal:** Nominate novel, druggable, disease-relevant regulators of CD4⁺ T cell
programs using the Marson–Pritchard genome-scale CRISPRi Perturb-seq atlas, and hand off
a ranked, evidence-backed target list.

**Guiding logic.** Each perturbation in this screen knocks down one gene and measures the
genome-wide transcriptional consequence. A perturbation therefore *is* a candidate target:
if knocking down gene *X* moves a therapeutically desirable program (e.g. suppresses a
pathogenic cytokine, expands a regulatory program), then a drug that inhibits *X* is
predicted to reproduce that effect. Target discovery here is **reading the perturbation →
program map backwards**: from the program we want to change, to the upstream regulators that
control it, filtered to those that are reproducible, disease-linked, and druggable.

---

## 0. Dataset in one screen

- **Assay:** genome-scale CRISPR interference (CRISPRi), probe-based Perturb-seq (10x Flex).
- **Scale:** ~22 M primary human CD4⁺ T cells, 4 donors (D1–D4), 3 contexts
  (**Rest, Stim8hr, Stim48hr** = TCR/CD28 re-stimulation timecourse).
- **Reference:** Zhu, Dann, Yan, Reyes Retana, Goto, Guitche, Brixi, Ota, Pritchard,
  Marson. *Genome-scale perturb-seq in primary human CD4+ T cells maps context-specific
  regulators of T cell programs and human immune traits.* bioRxiv 2025.12.23.696273
  (v1, Dec 2025); also SSRN 6137047.
- **Analysis layers available (S3 `s3://genome-scale-tcell-perturb-seq/marson2025_data/`):**
  - `GWCD4i.pseudobulk_merged.h5ad` (44.6 GB) — per (guide × donor × condition) pseudobulk,
    18,129 genes, with QC/eligibility flags (`keep_for_DE`, `keep_test_genes`, `n_cells`…).
  - `GWCD4i.DE_stats.h5ad` (16.8 GB) — **primary analysis object.** 33,983 (perturbation ×
    condition) rows × 10,282 gene columns. `.obs` carries on-target effect, trans-effect
    counts, reproducibility and off-target flags; `.layers` carry per-gene `log_fc`,
    `zscore`, `p_value`, `adj_p_value`, `baseMean`, `lfcSE`.
  - `GWCD4i.DE_stats.by_guide.h5mu` (29.4 GB) — per-guide DE (guide_1 / guide_2 modalities).
  - `GWCD4i.DE_stats.by_donors.h5mu` (16.9 GB) — per-donor-pair DE.
  - `suppl_tables/` — `DE_stats.suppl_table.csv` (tabular `.obs`), `sgrna_library_metadata`,
    `sample_metadata`, per-lane QC.
  - Cell-level `D*_*.assigned_guide.h5ad` — 12 files, ~1.8 TB total (fetch selectively).

**Built-in confidence filters (use these throughout, do not re-derive):**
`ontarget_significant` + `ontarget_effect_size` (did the knockdown work),
`guide_correlation_signif` (the two guides agree), `donor_correlation_*` (donors agree),
`n_downstream` / `n_total_de_genes` (trans-effect magnitude = "hub-ness"),
and exclusion flags `distal_offtarget_flag`, `neighboring_gene_KD`, `low_target_gex`.

---

## Phase A — Foundation (data QC, schema, positive controls)

Deliverables that everything downstream depends on. Runs on the DE-stats + pseudobulk
objects (no cell-level data needed).

1. **Schema & QC audit.** Load `DE_stats.h5ad`; confirm dimensions, enumerate `.obs`/`.var`
   fields, tabulate perturbations per condition, cells-per-target distribution, fraction
   passing `keep_for_DE`/`keep_test_genes`, and knockdown efficiency
   (`ontarget_significant` rate, `ontarget_effect_size` distribution). → `qc_overview.png`,
   `de_schema.csv`.
2. **Reproducibility landscape.** Distributions of `guide_correlation_signif` and
   `donor_correlation_hits_mean`; define a **"high-confidence" perturbation set** =
   on-target significant AND (two-guide agreement OR cross-donor agreement) AND not flagged
   for off-target/low-expression. This set is the universe for all target nomination.
   → `reproducibility.png`, `high_confidence_targets.csv`.
3. **Positive-control benchmark.** Confirm known T-cell drug-target biology reproduces:
   e.g. knockdown of the IL-2/JAK-STAT axis, NFAT/calcineurin pathway, costimulation
   regulators, and TF master regulators (TBX21, GATA3, FOXP3, RORC) produce the expected
   downstream signatures. Establishes that "signature → drug effect" reasoning is valid
   here before we trust novel hits. → `positive_controls.png`, `positive_controls.csv`.

## Phase B — Target-discovery directions (parallel workstreams)

Each direction produces a candidate list keyed by (gene, condition) with an effect summary.
They are complementary lenses on the same DE matrix; a target surfacing in several is
stronger.

### Direction 1 — Master regulators of effector cytokine programs
Rank perturbations by their signed effect on therapeutically central cytokines and their
receptors: **IL2, IFNG, IL4, IL13, IL17A, IL21, TNF, IL10, IL2RA, CTLA4**. Knockdowns that
*suppress* pathogenic cytokines → candidate **anti-inflammatory / autoimmunity** targets;
knockdowns that *boost* effector output or IL2RA → candidate **immuno-oncology /
immunodeficiency** targets. → `cytokine_regulators.csv`, per-cytokine ranked bar plots.

### Direction 2 — Regulators of helper-lineage polarization
Score each perturbation on lineage-defining programs — Th1 (TBX21, IFNG, CXCR3),
Th2 (GATA3, IL4, IL5, IL13), Th17 (RORC, IL17A/F, IL23R), Treg (FOXP3, IL2RA, IKZF2),
Tfh (BCL6, CXCR5) — using curated program gene sets scored on the DE `zscore` layer.
Nominate regulators that shift the Th1/Th2 or effector/Treg balance. This mirrors the
paper's Th1/Th2 validation axis. → `polarization_regulators.csv`, program-score heatmap.

### Direction 3 — Context-specific (stimulation-dependent) regulators
Contrast each perturbation's effect across Rest → Stim8hr → Stim48hr. Prioritize targets
active **only in stimulated cells** — attractive because a drug would act on activated,
disease-driving T cells while sparing the resting repertoire. Classify each target as
constitutive / activation-induced / rest-specific. → `context_specific_targets.csv`,
condition-interaction heatmap.

### Direction 4 — Network hubs / master regulators (topology)
Build the perturbation → downstream-gene bipartite network from significant DE edges;
compute out-degree (`n_downstream`), identify hub regulators and co-regulated gene
**modules** (clustered DE signatures). Hubs upstream of many disease-relevant genes are
high-value, mechanism-defining targets. → `regulator_network_hubs.csv`, module heatmap,
network figure.

### Direction 5 — Human-genetics integration (disease relevance)
Intersect the high-confidence regulator set (and the genes each controls) with autoimmune/
allergic/immune GWAS loci — RA, IBD/Crohn's/UC, MS, T1D, SLE, asthma, allergy, psoriasis.
A regulator is prioritized when it is a network hub controlling disease-program genes AND
sits at / colocalizes with a disease locus. Pull associations via **Open Targets** and
**GWAS Catalog** (through the literature/omics connectors + Open Targets API).
→ `gwas_regulator_overlap.csv`, disease × regulator dot plot.

### Direction 6 — Druggability & tractability annotation
Annotate every candidate with: protein class (kinase / GPCR / enzyme / TF / surface),
small-molecule & antibody **tractability**, existing drugs / tool compounds, and known
mechanism — via **Open Targets tractability** and **ChEMBL** (bioactivities, mechanisms).
Flag *novel* (no approved drug for T-cell / autoimmune indications) but *tractable*
regulators as the highest-interest set; TF hubs are noted as high-impact-but-hard.
→ `druggability_annotation.csv`.

## Phase C — Integration & hand-off

7. **Composite target scorecard.** Integrate the six directions into one ranked table:
   per candidate, an effect-strength score, reproducibility score, context specificity,
   network hub score, disease-genetics score, and druggability tier — combined into a
   transparent composite rank with the supporting evidence columns retained. Separate the
   list by therapeutic direction (suppress-inflammation vs enhance-immunity).
   → `TARGET_SCORECARD.csv` (primary deliverable), `top_targets_overview.png`.
8. **Target dossiers + report.** For the top ~15–25 candidates, a one-paragraph dossier
   (what it does in T cells, which program its knockdown moves and in which context,
   reproducibility, disease link, druggability, and the proposed therapeutic hypothesis
   incl. inhibit-vs-activate direction). Compile into `TARGET_REPORT.md` with embedded
   figures. Explicitly flag which are known vs genuinely novel.

## Phase D — (Optional) Cell-level validation
For a shortlist, pull the relevant cell-level `D*_*.assigned_guide.h5ad` shards to confirm
effects at single-cell resolution and map perturbation effects onto T-cell state manifolds
(effector / memory / Treg / exhaustion). Gated on shortlist size given ~150 GB per shard.

---

## Method & rigor notes
- **Primary object is `DE_stats.h5ad`** — precomputed genome-wide DE (DESeq2) with the
  `zscore` layer as the workhorse signal; avoids re-running DE on 22 M cells.
- **Nominate only reproducible, on-target, off-target-clean perturbations** (Phase A2 set).
  Reproducibility (two-guide + cross-donor) is the main defense against artifacts — the
  paper itself flags non-replicating trans-effects (e.g. IL2RA, NFAT5, TNFAIP3).
- **Direction matters:** always record whether the therapeutic hypothesis is to *inhibit*
  or *activate* the target, since CRISPRi knockdown models loss-of-function only.
- **External data via connectors:** Open Targets (disease association + tractability),
  ChEMBL (bioactivity/mechanism), literature graph, GWAS Catalog — all read-only lookups.
- **Compute:** DE-stats/pseudobulk analysis fits the connected rocm-4 node (CPU, <200 GB);
  cell-level shards, if needed, warrant heavier jobs. See `CLAUDE.md` for env & paths.
