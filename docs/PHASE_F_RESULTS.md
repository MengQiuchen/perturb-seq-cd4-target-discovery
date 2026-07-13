# PHASE F — Translational Target Dossiers (Structure, Chemistry, Tractability, Genetics)

**Project:** Novel drug targets in CD4⁺ T-cell genome-scale CRISPRi Perturb-seq (Marson–Pritchard atlas)
**Scope:** For the **19 single-cell-validated targets** (Phase D/E), assemble the translational evidence a target-selection committee needs — protein structure (AlphaFold + experimental PDB), chemical matter (ChEMBL), modality tractability + clinical precedent (Open Targets), and human genetics — and rank them by a transparent translational-readiness heuristic.
**Compute:** `ssh:clust1-rocm-4`, env `perturb-seq`, external evidence via connectors (AlphaFold DB, RCSB PDB, ChEMBL, Open Targets, GWAS Catalog). **Outputs:** `phaseF_outputs/` — master table `phaseF_master_druggability.csv`, per-track dossiers `phaseF_trackA_signaling.json` / `phaseF_trackB_chromatin_tf.json`, readiness ranking `phaseF_readiness_ranked.json`, structure profiles `phaseF_plddt_profiles.json`, and four figures. Full prose dossiers in `TARGET_DOSSIERS.md`.

## Why this phase

Phases A–E answered *whether* these 19 genes are real, reproducible, stimulation-dependent regulators of CD4⁺ T-cell programs. Phase F asks the orthogonal, translational question: **if you wanted to drug one of them, where would you start, and which are the least-crowded opportunities?** It reads the CRISPRi phenotype as a genetic model of drug-induced loss of function and layers on the four axes that decide whether a validated regulator is a viable drug target — is there a folded pocket to bind, is there chemical matter already, what modality is tractable, is the target genetically tied to the disease, and is anyone already there. No new perturbation data is generated; this is an evidence-integration and prioritization phase.

## Headline

**18 of 19 targets are clinically unprecedented** — only **CD3E** carries clinical-stage drugs against the target (the anti-CD3 antibody class). **All 12 chromatin/transcription-factor "novel-axis" targets have zero known drugs.** 17/19 have experimental PDB structures (LAT and TADA2B are AlphaFold-only). The top-ranked opportunity is **STAT6** — a genetically anchored, chemically mature, drug-naive master Th2 transcription factor — followed by the novel kinase **SIK3** and a cluster of entirely undrugged SAGA/Mediator readers (**MED24, SGF29, TADA2B**).

---

## Cross-track translational-readiness ranking (all 19)

| Gene | Axis | KD% | AF pLDDT | #PDB | SM buckets | ChEMBL acts | max pChEMBL | Known drugs | Immune GWAS | Readiness |
|---|---|---|---|---|---|---|---|---|---|---|
| **STAT6** | chr/TF | 82 | 77 | 7 | 3 | 552 | 9.15 | 0 | 55 | **6.5** |
| **ZAP70** | sig | 91 | 85 | 15 | 4 | 2390 | 8.10 | 0 | 1 | **6.1** |
| **SIK3** | sig | 87 | 51 | 5 | 3 | 809 | 9.34 | 0 | 0 | **4.5** |
| MED24 | chr/TF | 96 | 84 | 10 | 1 | 0 | — | 0 | 12 | 4.2 |
| PLCG1 | sig | 88 | 83 | 6 | 1 | 438 | 6.75 | 0 | 0 | 4.0 |
| VAV1 | sig | 92 | 86 | 10 | 1 | 1 | 8.98 | 0 | 0 | 4.0 |
| NSD1 | chr/TF | 92 | 45 | 4 | 2 | 147 | 6.96 | 0 | 1 | 3.1 |
| SEL1L | chr/TF | 94 | 81 | 5 | 0 | 0 | — | 0 | 0 | 3.0 |
| CD247 | sig | 72 | 62 | 38 | 1 | 0 | — | 0 | 46 | 2.5 |
| SMARCE1 | chr/TF | 88 | 70 | 9 | 1 | 8 | 6.95 | 0 | 37 | 2.5 |
| SGF29 | chr/TF | 94 | 92 | 8 | 1 | 1 | — | 0 | 5 | 2.5 |
| CD3E | sig | 96 | 73 | 44 | 2 | 0 | — | 22 | 3 | 1.8 |
| ARNT | chr/TF | 89 | 56 | 46 | 2 | 25 | — | 0 | 1 | 1.6 |
| TADA2B | chr/TF | 92 | 87 | 0 | 0 | 0 | — | 0 | 0 | 1.5 |
| SMARCB1 | chr/TF | 48 | 81 | 18 | 1 | 8 | 8.03 | 0 | 0 | 1.5 |
| TRIP12 | chr/TF | 84 | 67 | 5 | 1 | 0 | — | 0 | 0 | 1.5 |
| CHD4 | chr/TF | 52 | 65 | 12 | 1 | 270 | 8.52 | 0 | 0 | 1.0 |
| LAT | sig | 85 | 59 | 0 | 0 | 2 | — | 0 | 0 | 0.0 |
| MED12 | chr/TF | 85 | 65 | 3 | 0 | 6 | 6.89 | 0 | 0 | 0.0 |

*Axis: sig = TCR-proximal signaling; chr/TF = chromatin/transcription-factor novel axis. KD% = Phase E mean on-target knockdown. **The readiness score is a prioritization aid, not a cutoff** — it deliberately rewards existing chemical matter and genetics, so it under-scores genuinely novel but chemically virgin targets (e.g. the SAGA/Mediator readers), which are read from the structure and axis columns instead.*

### Readiness score (transparent heuristic)

Readiness = **structure** (pLDDT ≥ 80 → +1.5; longest confident domain ≥ 150 aa → +1.0; ≥ 5 experimental PDB → +0.5) + **chemical matter** (≥ 2 SM tractability buckets → +1.0; ≥ 100 ChEMBL bioactivities → +1.0; max pChEMBL ≥ 8 → +1.0) + **genetics** (min(immune-GWAS / 20, 1) × 2.0) − **essentiality risk** (Phase E KD% < 60% → −1.5). The weights are a judgment, not a fit; every component column is retained in `phaseF_master_druggability.csv` so the ranking can be re-derived or re-weighted for a different translational priority.

---

## 1. Structure — AlphaFold confidence + experimental coverage

Per UniProt accession: the AlphaFold DB canonical model (v6), per-residue pLDDT read from the B-factor column, confident folded domain = a contiguous run ≥ 30 residues at pLDDT ≥ 70; experimental coverage = RCSB PDB match count. The 19 split cleanly into three structural classes:

- **Crisp, high-confidence single-domain targets** — clean starting points for structure-based design. **SGF29** (mean pLDDT 92; Tudor reader domain res 201–288 at pLDDT 96.6), **ZAP70** (85; tandem SH2 res 3–256 + kinase domains, all pLDDT 92–93), **VAV1** (86), **TADA2B** (87), **MED24** (84), **SMARCB1** (81), **SEL1L** (81; a single 370-aa Sel1-repeat solenoid at pLDDT 92.9), and **STAT6** (77 overall but with a 236-aa DNA-binding/SH2 module res 398–633 at pLDDT 94.9).
- **Large scaffolds with discrete ordered modules** — the drug target is the module, not the chain. **NSD1** (2,696 aa, overall pLDDT 45 but a confident SET methyltransferase-region domain res 1631–1820 at pLDDT 89.5), **CHD4** (1,912 aa), **TRIP12** (1,992 aa), **MED12** (2,177 aa), and **SIK3** (1,321 aa overall pLDDT 51, but a crisp catalytic kinase domain res 60–336 at pLDDT 94.0). For these, the confident-domain coordinates in `phaseF_plddt_profiles.json` define the actual druggable surface.
- **Disordered adaptors / receptor chains** — no confident folded domain (**LAT** mean pLDDT 59, **CD247**/CD3ζ 62), pointing to PPI-disruption, degrader, or antibody logic rather than a classical pocket. CD247 nonetheless has **38 experimental PDB** entries (captured in TCR/CD3-complex structures), and CD3E **44** — abundant complex structures despite modest monomer confidence.

![**Fig F1. Cross-track translational-readiness ranking.** Score = structure + chemical matter + genetics − essentiality risk. ○ = clinically unprecedented (0 known drugs); ⚠ = essentiality/toxicity risk (KD < 60%). STAT6 leads; nearly every target is clinically unprecedented.](../results/phaseF/phaseF_readiness_ranking.png)

![**Fig F3. Track A (signaling) AlphaFold confidence maps.** Per-residue pLDDT with confident folded domains (blue bars, ≥ 30 aa at pLDDT ≥ 70). ZAP70 and SIK3 show crisp catalytic domains; LAT and CD247 are largely disordered adaptors (PPI / degrader logic).](../results/phaseF/phaseF_trackA_structure_maps.png)

![**Fig F4. Track B (chromatin/TF) AlphaFold confidence maps.** SGF29 and TADA2B are compact high-confidence readers; NSD1, CHD4, MED12 and TRIP12 are large scaffolds where only discrete catalytic/reader modules are ordered — those modules are the actual drug targets.](../results/phaseF/phaseF_trackB_structure_maps.png)

## 2. Chemical matter — ChEMBL bioactivities

ChEMBL bioactivity count, max pChEMBL (over a 200-activity representative sample), and sub-µM compound count per single-protein human target. Chemical maturity is concentrated in a handful of targets:

- **ZAP70** — the richest signaling target: 2,390 bioactivities, 20 sub-µM compounds, max pChEMBL 8.1 (potent ATP-competitive kinase inhibitors documented since the late 1990s).
- **SIK3** — 809 bioactivities, 18 sub-µM, max pChEMBL 9.34; the ligand set is dominated by multi-kinase inhibitors that bind its catalytic domain (staurosporine Kd 0.46 nM, crenolanib 2 nM, dasatinib 28 nM) — chemical proof the pocket is drug-accessible even though no SIK3-selective drug exists.
- **STAT6** — 552 bioactivities and the deepest potent-compound bench of any target (**70 sub-µM**, max pChEMBL 9.15), unusual for a transcription factor and reflecting an SH2-domain-directed medicinal-chemistry history.
- **CHD4** (270 acts, max 8.52; dinaciclib binds), **PLCG1** (438, max 6.75), **NSD1** (147, max 6.96) — moderate benches. The remaining targets have thin or no ChEMBL matter (MED24, TADA2B, SEL1L, TRIP12, CD247 have zero) — expected for the novel chromatin/scaffold axis and the reason the readiness score under-weights them.

## 3. Tractability & clinical precedent — Open Targets

Small-molecule (SM), antibody (AB), and protein-degradation (PR) tractability buckets, plus known-drug / clinical-candidate counts:

- **CD3E is the only precedented target** — 22 known drugs against the target: the approved anti-CD3 antibodies **teplizumab** (type-1 diabetes) and **otelixizumab** (phase 3), **foralumab** (phase 2), plus a large family of CD3-engaging bispecifics (blinatumomab, tarlatamab, epcoritamab, mosunetuzumab, …). It is `known-drug-target`, not a novel opportunity, and appears here as the pipeline's *de novo* positive control — the logic recovered a real, drugged immune target.
- **All 18 other targets carry zero known drugs.** Small-molecule tractability is strongest for the kinases and TFs (ZAP70 4 SM buckets; SIK3 and STAT6 3 each; NSD1, ARNT 2), while the disordered adaptors lean antibody/degrader (LAT 4 AB buckets, 0 SM; CD247 3 AB; SEL1L AB-only). Protein-degradation tractability is broad (most targets 2–4 PR buckets), consistent with the CRISPRi loss-of-function readout translating naturally to a degrader modality.

![**Fig F2. Druggability landscape.** Structural readiness (pLDDT + longest folded domain), modality tractability (SM/AB buckets), and chemical matter (ChEMBL bioactivities + max potency) for all 19 targets, colored by axis.](../results/phaseF/phaseF_druggability_landscape.png)

## 4. Two mechanistic axes and the top opportunities

The 19 targets fall into the two axes that ran through the whole project, and they present very differently for drug discovery:

**Axis 1 — TCR-proximal signaling (7 targets: LAT, PLCG1, CD247, CD3E, ZAP70, VAV1, SIK3).** A clinically validated pathway (CD3E, ITK already drugged) with several drug-naive entry points. ZAP70 and SIK3 are the standouts — folded catalytic domains, real chemical matter, clean clinical whitespace. LAT and CD247 are disordered adaptors better suited to antibody or degrader approaches. All are `suppress-inflammation` (knockdown blocks the activation program), so the therapeutic logic is inhibition for autoimmune/allergic disease.

**Axis 2 — chromatin / transcriptional co-activators (12 targets: SMARCE1, SGF29, MED24, MED12, TADA2B, NSD1, CHD4, SMARCB1, STAT6, TRIP12, ARNT, SEL1L).** Structurally tractable, largely undrugged, several with strong autoimmune GWAS signal — the program's most novel output. The SAGA/Mediator readers (MED24, SGF29, TADA2B) are compact, high-confidence, entirely drug-naive folds — ideal structure-based-design starts (reader-domain antagonism, degrader, or PPI disruption).

**Top opportunities:**

1. **STAT6 (rank #1)** — master Th2 transcription factor; **55 asthma/allergy GWAS associations** (max −log₁₀P ≈ 38.5), the richest chemical matter of the set (552 bioactivities, 70 sub-µM, max pChEMBL 9.15), a druggable-family SM bucket and a confident SH2/DNA-binding module — yet **no approved STAT6 drug**. Therapeutic logic: STAT6 **inhibition** for type-2 inflammation (the same axis dupilumab targets upstream at IL4R).
2. **SIK3 (rank #3)** — the one novel kinase; a crisp catalytic domain (res 60–336, pLDDT ≈ 94), 809 bioactivities, ligand-bound experimental structures (8R4O series), and clean clinical whitespace. An ATP-competitive or allosteric inhibitor is the obvious modality.
3. **MED24 / SGF29 / TADA2B** — structurally crisp, entirely drug-naive SAGA/Mediator readers and scaffolds; the highest-novelty starts for reader-domain antagonism, degrader, or PPI-disruption chemistry. MED24 additionally carries 12 asthma/allergy GWAS associations.
4. **Highest-confidence novel candidates** (≥ 3-direction convergence + novel-druggable + immune-GWAS support, from Phase B/C): **SMARCE1, CD247, STAT6, ZAP70, TRIP12, NSD1, ARNT, SEL1L** — the set that is strong on functional perturbation evidence, population genetics, and a tractable protein simultaneously.

### Essentiality-flagged (proceed with a window assessment)

**CHD4 (KD% 52%) and SMARCB1 (48%)** have the weakest Phase E knockdowns and are core, broadly essential remodeler-complex subunits. Their readiness carries the −1.5 essentiality penalty; a narrow therapeutic window is likely, and any program against them needs a careful selectivity/window assessment before it is worth chemical investment.

---

## Caveats

- **Readiness is a heuristic, not a fit.** The weights (structure/chemistry/genetics/essentiality) are a transparent judgment call. It intentionally rewards *existing* chemical matter and genetics, so it systematically under-ranks the most novel targets (drug-naive readers with no ChEMBL history). Read the ranking together with the axis and structure columns, not as a hard order.
- **Max pChEMBL is a representative ceiling.** It is taken over a 200-activity sample per target, not the full activity set, so it indicates whether potent chemical matter *exists*, not the single best compound ever reported.
- **Structure is monomeric AlphaFold.** pLDDT confidence is for the isolated chain; disordered adaptors (LAT, CD247) fold in complex, and their many experimental PDB entries are complex structures. "No confident domain" means *no monomeric pocket*, not "undruggable" — it redirects modality (PPI/degrader/antibody), it does not disqualify.
- **Tractability buckets are Open Targets annotations**, i.e. genome-informatics predictions of modality feasibility, not experimental confirmation of a bindable site.
- **Genetics coverage is uneven.** Immune-GWAS counts come from the Phase B/C disease-association layer; a zero (e.g. LAT, PLCG1, VAV1, SIK3) means no common-variant disease signal was found, which lowers readiness but does not argue against the functionally convergent evidence that nominated the target.
- **Phase F does not test the therapeutic direction.** Whether to inhibit or activate each target is inferred from the CRISPRi knockdown direction (Phase B/C); confirming it needs a CRISPRa arm or an arrayed functional assay (a Phase-H/follow-up item).

## Conclusion

Across the 19 single-cell-validated regulators, the pipeline recovers exactly one already-drugged immune target (**CD3E**) — end-to-end validation of the logic — and nominates **18 clinically unprecedented targets**, all 12 chromatin/transcriptional novel-axis genes among them with zero known drugs. The most translationally ready opportunity is **STAT6** (genetics + chemistry + structure + drug-naive), with **SIK3** the standout novel kinase and **MED24/SGF29/TADA2B** the highest-novelty structure-based-design starts. **CHD4 and SMARCB1** are flagged for essentiality/window risk. Together with the per-target prose dossiers in `TARGET_DOSSIERS.md`, this phase converts the validated regulator list into a prioritized, evidence-backed set of therapeutic entry points for autoimmune/allergic and immuno-oncology indications.

## Files (`phaseF_outputs/`)

- `PHASE_F_RESULTS.md` — this report
- `phaseF_master_druggability.csv` — **master table**: all 19 targets × structure / chemistry / tractability / clinical / genetics columns
- `phaseF_readiness_ranked.json` — readiness score + all component sub-scores, ranked
- `phaseF_trackA_signaling.json`, `phaseF_trackB_chromatin_tf.json` — full per-target dossiers (AlphaFold, PDB, ChEMBL compounds, Open Targets drugs, tractability, Phase E validation)
- `phaseF_plddt_profiles.json` — per-target pLDDT profile + confident-domain coordinates
- `phaseF_readiness_ranking.png` (Fig F1), `phaseF_druggability_landscape.png` (Fig F2), `phaseF_trackA_structure_maps.png` (Fig F3), `phaseF_trackB_structure_maps.png` (Fig F4)
- `structures/` — downloaded AlphaFold `.cif`/`.pdb` models for the 19 targets
- `TARGET_DOSSIERS.md` (project root) — full narrative per-target dossiers for the top targets

**Downstream:** Phase G — single-cell analyses (dose–response, state redistribution, trajectory, dual-guide) that only single-cell resolution can produce, deepening the mechanistic case for the same 19 targets.
