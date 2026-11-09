# ECM gene program in clear-cell renal cell carcinoma

Analysis for the Opeongo computational biology home assignment: where a fixed 10-gene
extracellular-matrix (ECM) panel is expressed across cell populations in three ccRCC single-cell
cohorts, whether it differs between tumor and adjacent non-tumor tissue in one population, whether
combining cohorts changes that, and whether linked bulk RNA and protein data agree.

Panel: COL1A1, COL1A2, COL3A1, COL5A1, COL6A1, COL6A2, FN1, SPARC, DCN, LUM.

**Author:** Leviel Fluhr · **Time spent:** about 6 hours · **AI tools:** used, with every output checked
(see [AI assistance](#ai-assistance)).

**For reviewers, in reading order:**
1. the 4-page report, [`report/report.pdf`](report/report.pdf);
2. the analysis notebook, [`notebooks/ecm_analysis.ipynb`](notebooks/ecm_analysis.ipynb);
3. [How decisions were made](#how-decisions-were-made), for the reasoning behind each stage.

Everything can be regenerated with `uv sync && uv run opengo all` (see [Setup](#setup)).

**Deliverables**

- **Report (4 pages):** [`report/report.pdf`](report/report.pdf) (source `report/report.html`,
  figures `report/figures/`).
- **Analysis notebook:** [`notebooks/ecm_analysis.ipynb`](notebooks/ecm_analysis.ipynb), executed,
  with paired source `notebooks/ecm_analysis.py`. It walks through every stage's key numbers and figures.
- **Decision record:** "How decisions were made" below, plus
  [`outputs/tables/decisions_log.tsv`](outputs/tables/decisions_log.tsv) (inclusion/exclusion
  record) and [`PLAN.md`](PLAN.md) (plan, review gates, change log).

**Main findings**

1. Cohort C was excluded: its counts are capped at ≤ 409 non-mitochondrial counts per cell.
2. A (single-nucleus) and B (single-cell) were combined with Harmony (batch = cohort).
3. The ECM panel is a stromal program. In the fibroblast + mural compartment, tumor exceeds adjacent
   in 9/9 cohort B donors (+1.04 SD [0.63, 1.46]), robust to the fibroblast:pericyte mix.
4. Cohort A (3 pairs) is too imprecise to confirm or refute this (its 95 % CI, −1.77 to +1.81,
   contains B's estimate). The pooled estimate is therefore dominated by B, and how it is pooled
   matters (fixed effect +0.88 [0.56, 1.21]; random effects +0.60 [−0.39, 1.59]).
5. Full CPTAC bulk RNA (72 pairs) and protein (80 pairs) show the same gene pattern: collagens, FN1
   and SPARC up; DCN/LUM not. This may reflect per-cell change, expansion of mural and vascular
   tissue, or both. The 3 matched multimodal donors are too few to decide.

## Inputs (not committed)

The supplied assignment package is not in this repository. Place its contents at the repository
root so that the layout is:

```
data/rna/cohort_A_GSE227898/      34 GEO sample archives (10x Matrix Market, raw droplets)
data/rna/cohort_B_GSE178481/      GSE178481_RAW.tar (26 count CSVs)
data/rna/cohort_C_GSE159115/      13 legacy 10x HDF5 files
data/bulk/                        CPTAC bulk RNA and protein matrices (tumor / normal)
metadata/                         manifests, label mapping, ECM panel, multimodal crosswalk
metadata/original/                author annotations and specimen metadata
FILE_MANIFEST.tsv                 sizes + SHA-256 of every supplied file (committed)
```

The supplied package `README.md` is kept, byte-identical, at `docs/supplied/README.md`, so that
this README can sit at the root. `opengo verify` accounts for that move.

## Setup

Requires [uv](https://docs.astral.sh/uv/). The environment is fully specified by
`pyproject.toml` + `uv.lock` (Python 3.11).

```bash
uv sync
```

```bash
uv run opengo verify
```

`verify` checks every supplied file against `FILE_MANIFEST.tsv` and writes
`outputs/tables/input_verification.tsv`.

## Execution order

The whole analysis runs from the supplied inputs with one command (≈ 25 min on a laptop, peak
≈ 2.6 GB RAM, as measured in the verified rebuild; headless Chrome or Chromium is needed only for
the final PDF step):

```bash
uv run opengo all
```

The steps, which can also be run individually:

| Step | Command | Output | Status |
|---|---|---|---|
| 0 | `uv run opengo verify` | `outputs/tables/input_verification.tsv` | done |
| 1 | `uv run opengo load all` (or `A` / `B` / `C`), then `uv run opengo genes` | raw-count checkpoints in `outputs/checkpoints/`; `outputs/tables/stage1_*.tsv`, `gene_correspondence*.tsv`, `ecm_panel_presence.tsv` | done (~10 min, peak ~1.9 GB RAM) |
| 1 fig | `uv run opengo figures 1` | `outputs/step_1_output/*.png` | done |
| 2 | `uv run opengo qc`, then `uv run opengo explore`, then `uv run opengo figures 2` | per-cell QC flags, pseudobulk checkpoint, `outputs/tables/stage2_*.tsv`, `outputs/step_2_output/*.png` | done (~20 min) |
| 3 | `uv run opengo integrate`, then `uv run opengo figures 3` | embeddings + metrics, `outputs/tables/stage3_*.tsv`, `outputs/step_3_output/*.png` | done (~6 min, peak ~2.7 GB) |
| 4 | `uv run opengo overview`, then `uv run opengo figures 4` | ECM overview, panel coherence, focus coverage, `outputs/tables/stage4_*.tsv`, `outputs/step_4_output/*.png` | done (~1 min) |
| 5 | `uv run opengo pseudobulk`, then `uv run opengo figures 5` | paired tumor vs adjacent effects, meta-analysis, robustness, `outputs/tables/stage5_*.tsv`, `outputs/step_5_output/*.png` | done (< 1 min) |
| 6 | `uv run opengo multimodal`, then `uv run opengo figures 6` | crosswalk validation, cross-assay agreement, composition proxy, `outputs/tables/stage6_*.tsv`, `outputs/step_6_output/*.png` | done (< 1 min) |
| 7 | `uv run opengo figures report`, then `uv run opengo report` | `report/figures/*.png`, `report/report.pdf` (fails if > 4 pages) | done |
| nb | `uv run jupytext --to ipynb --execute notebooks/ecm_analysis.py` | `notebooks/ecm_analysis.ipynb` (reads the tables; seconds) | done |

Each command reads the supplied inputs or earlier checkpoints and can be re-run from a clean session.

**Verified clean rebuild (2026-11-07, 17:11–17:36).** `opengo all` was run end to end from the raw
inputs into a separate folder:

```bash
OPENGO_OUTPUTS_DIR=rebuild_check/outputs OPENGO_REPORT_DIR=rebuild_check/report uv run opengo all
```

- It finished in 25 min (exit 0, peak 2.6 GB RAM) and produced a 4-page report.
- Every result table was compared with the committed outputs
  ([`outputs/tables/rebuild_comparison.tsv`](outputs/tables/rebuild_comparison.tsv)):
  - 63 of 65 tables are identical;
  - one differs at the 1e-9 level (floating point);
  - in one, only the Leiden-based ARI/NMI of the Harmony embeddings differ (see Stage 3).
- No number quoted in the report changed.

## How decisions were made

Every stage follows the same loop: **spot an issue → decide how to test it → look at the evidence
(a table or figure) → act, and log the decision with its evidence**. Nothing downstream runs until
the stage has been reviewed. The full record, with alternatives considered, is in
[`outputs/tables/decisions_log.tsv`](outputs/tables/decisions_log.tsv); the rationale for each
stage's choices is in [`PLAN.md`](PLAN.md).

### Stage 1: loading and reconciliation

Figures: [`outputs/step_1_output/`](outputs/step_1_output/) (regenerate with `uv run opengo figures 1`).

**1. Did we load exactly the cells the authors analysed?**
- *Issue:* cohort A is deposited as raw droplet matrices (6.8 M barcodes per specimen); B and C as
  processed matrices. A silent join error would change every downstream number.
- *Check:* join A barcodes to the author annotation by specimen + barcode; compare cells per
  specimen with the package manifest; pytest assertions on the documented totals.
- *Evidence:* `02_load_reconciliation.png`: all 73 specimens on the identity line; A 141,950,
  B 157,881, C 29,499 cells, every author annotation matched, all counts non-negative integers.
- *Action:* A restricted to annotated nuclei (no re-calling, as instructed); B and C keep all
  deposited cells until QC.

**2. Can genes be compared across cohorts?**
- *Issue:* A and C use Ensembl IDs from different releases; B has symbols only, from an older
  annotation. Matching by symbol alone can silently pick the wrong gene.
- *Check:* A↔C by Ensembl ID; each B symbol assigned only if it matches exactly one A (else C)
  symbol; collisions and ambiguities counted rather than guessed.
- *Evidence:* `05_gene_correspondence.png`: 29,531 genes shared by all three cohorts; 532 B
  symbols (old and new names of the same 266 genes, mostly lncRNA clone names) left out; all 10
  ECM genes matched by Ensembl ID in every cohort.
- *Action:* integration uses the shared set; a gene absent from a cohort is treated as unmeasured,
  never as zero.

**3. Is a suggested cell label wrong? (TRM)**
- *Issue:* the supplied `label_mapping.tsv` maps cohort A "TRM" (510 nuclei) to Lymphoid, but the
  authors' own grouping column places all 510 under Macrophages. "TRM" is ambiguous: tissue-resident
  *macrophages* or tissue-resident *memory T cells*.
- *First instinct, and why it was not enough:* the initial proposal was to remap from the label
  columns alone. The Gemini gatekeeper review (independent of the working agent) pointed out that labels
  can carry clerical errors and that lineage marker genes should decide. We agreed and tested it.
- *Check:* fraction of nuclei expressing myeloid vs T-cell markers in TRM, compared with reference
  author groups (macrophages, monocytes, CD4/CD8 T cells, Tregs, NK, tumor).
- *Evidence:* `06_cohortA_label_check.png` (the label conflict) and `07_TRM_marker_check.png`
  (markers). TRM nuclei: CD3D 0 %, CD3E 1 %, CD2 1 %, CD8A 2 % (T cells 11–66 %); CSF1R 50 %
  (macrophages 51 %); PTPRC 94 %. Resident-T markers CD69/ITGAE are at macrophage, not T-cell, levels.
- *Action:* TRM remapped to **Myeloid** (tissue-resident macrophages), logged as a deliberate
  deviation from the suggested guide; the original suggestion is kept in `population_suggested`.
  Caveat noted: TRM show a weaker classic macrophage profile (CD163 13 % vs 55 %), so they are
  re-examined in QC. Impact on the ECM question is negligible (0.36 % of A, immune cells).

**4. Which populations can support a tumor-vs-adjacent comparison?**
- *Issue:* donors, not cells, are the replicates; a population needs enough cells in *both*
  conditions in the *same* donors.
- *Check:* count paired donors per population at minimum-cell thresholds 1–100.
- *Evidence:* `01_study_design.png`, `04_paired_donor_coverage.png`: at ≥ 10 cells the
  fibroblast/myofibroblast compartment has only 3 paired donors in A and 2 in B, whereas
  endothelial cells have 3 (A) + 9 (B) and B pericytes 9.
- *Action:* no exclusion yet; this shapes Stage 2 (re-examine B's fibroblast/pericyte/endothelial
  label conflicts and A myofibroblasts with non-panel markers) and the Stage 4 population choice.
  Panel genes (DCN, LUM, collagens) are never used to define cells, to avoid circularity.

**5. Are the cohorts comparable in composition and depth?**
- *Evidence:* `03_population_composition.png`: A tumors are mostly malignant nuclei, B tumors
  mostly immune cells, C unlabelled; mean counts per cell differ ~10-fold (C adjacent lowest).
- *Interpretation:* consistent with known single-nucleus vs single-cell capture biases
  (Wu et al. 2019; Denisenko et al. 2020).
- *Action:* population proportions are never compared across cohorts; tumor-vs-adjacent effects
  are estimated within cohort, with cohort as a model term; QC thresholds are set per specimen and
  per assay.

### Stage 2: quality control and exploratory structure

Figures: [`outputs/step_2_output/`](outputs/step_2_output/) (regenerate with `uv run opengo qc`,
`uv run opengo explore`, `uv run opengo figures 2`). Decisions approved at Check-in 2.

**1. Which cells are low quality, and is a flag technical or biological?**
- *Approach:* per-specimen, assay-aware flags (no fixed global cut-offs): counts/genes below
  median − 5 MAD, counts above median + 5 MAD, mitochondrial % above median + 3 MAD *and* above
  5 % (nuclei, where mitochondrial reads indicate ambient RNA) or 25 % (whole kidney cells, where
  proximal-tubule cells are naturally high), Scrublet doublets. Nothing is removed at this step.
- *Evidence:* `01_qc_distributions.png`, `02_flag_rates.png`. A and B were already filtered by
  their authors (A counts sit between 10³ and 10⁴; B's author QC fields match our recomputation
  exactly, `03_B_author_qc_vs_recomputed.png`). Flags are rare in the ECM-relevant populations
  (A stroma 1 %, A endothelium 4 %, B ≤ 1 %); they concentrate in A immune nuclei (high
  mitochondrial fraction, 10–22 %) and B non-malignant epithelium (12 %, proximal tubule biology).
- *Action (approved):* flagged cells (union of flags) removed in A (4.2 %) and B (0.9 %); no A or B
  specimen excluded. Removal is uneven across populations (A lymphoid 22 %, A myeloid 10 %, B
  non-malignant epithelium 13 %), which is why population proportions are never compared.

**2. Is cohort C usable?**
- *Issue noticed:* C's median mitochondrial fraction is 48 % and ~175 genes are detected per cell.
- *Check:* non-mitochondrial counts per cell, mitochondrial fraction per specimen, ECM panel
  detection by cohort.
- *Evidence:* `08_cohortC_integrity.png`: every C cell has **≤ 409 non-mitochondrial counts**
  (a hard ceiling; A and B reach 10⁴–10⁵); mitochondrial reads are 76 % of all C counts; the two
  worst specimens (median 92–96 % mitochondrial) are the tumors of two of C's four paired donors;
  DCN and LUM are detected in 0.2 % of C cells. Marker-based labelling still recovers broad
  populations in C (`06_cohortC_marker_labels.png`, after validating the procedure at 81 % / 95 %
  agreement with author labels in A / B, `05_marker_label_validation.png`), but per-cell ECM
  signal is too sparse to support a comparison.
- *Action (approved):* C excluded from integration and from all ECM analyses; the evidence is
  reported. A C-based sensitivity analysis (two clean pairs) was considered and not run, because
  per-cell panel detection is too sparse to be informative.

**3. What dominates variation between specimens?**
- *Check:* PCA of whole-specimen and endothelial pseudobulk, with R² of cohort / condition / donor
  for each component.
- *Evidence:* `04_pseudobulk_pca.png`: PC1 (83 %) is cohort (R² 0.99, i.e. single-nucleus vs
  single-cell); tumor vs adjacent appears only on later components, within cohort (endothelial PC4,
  R² 0.58).
- *Action:* confirms the plan: effects are estimated within cohort and combined with a cohort
  term or meta-analysis; never by pooling raw expression across cohorts.

**4. What is the stromal compartment made of?** (sets up the Stage 4 choice)
- *Check:* non-panel fibroblast (PDGFRA, C7, FBLN1, MFAP5), mural (RGS5, NOTCH3, MYH11, PDGFRB)
  and endothelial markers in A and B stromal labels, plus marker-based sub-clustering.
- *Evidence:* `07_stromal_identity.png`: A "Fibroblasts" are mostly mural-like (RGS5 50 %,
  PDGFRA 2 %); A "Myofibroblasts" are a fibroblast (C7 42 %) + mural mix; B "Fibroblast" is also
  RGS5-high; B fibroblast/pericyte label conflicts fall mostly with mural cells. In ccRCC the
  fibroblast–pericyte boundary is a continuum, not two clean types.
- *Implication:* a combined **fibroblast + mural (stromal) compartment** is defensible across A and
  B and has 3 (A) + 9 (B) paired donors at ≥ 10 cells, the same coverage as endothelium, but it is
  the compartment that actually produces the ECM panel. Decided at Stage 4 from the ECM dot plot.

### Stage 3: combining cohorts: unintegrated baseline vs Harmony

Figures: [`outputs/step_3_output/`](outputs/step_3_output/) (regenerate with
`uv run opengo integrate`, `uv run opengo figures 3`). Cohorts A + B (C excluded at Stage 2).
Embeddings use a stratified subsample (87,489 cells: ≤ 1,000 per specimen plus every endothelial,
fibroblast/myofibroblast and pericyte cell). **Embeddings are used for visualisation, evaluation and
label transfer only, never for testing** (tests use pseudobulk counts, Stages 4–5).

**1. Did the integration method actually run? (a tool problem caught by checking)**
- *Issue:* the first Harmony run reproduced the unintegrated metrics to three decimals.
- *Check:* Harmony on identical PCs with harmonypy 2.1.0 vs the long-established 0.0.10.
- *Evidence:* 2.1.0 stopped after one iteration (cohort silhouette 0.308 → 0.307); 0.0.10 ran
  2 × 19 iterations (0.31 → 0.25) and improved population separation (0.25 → 0.28).
- *Action:* harmonypy pinned to 0.0.10; a test now fails if Harmony stops correcting.

**2. Which embedding, and which batch key?**
- *Compared:* A alone; A + B unintegrated; A + B Harmony with batch = cohort; A + B Harmony with
  batch = specimen.
- *Evidence:* `02_integration_metrics.png`, `03_condition_structure_and_purity.png`,
  `04_label_transfer.png`.

| | A + B unintegrated | Harmony (cohort) | Harmony (specimen) |
|---|---|---|---|
| Within-population cohort mixing (batch ASW) | 0.49 | **0.81** | 0.43 |
| Cell-type separation (ASW) / ARI* | 0.62 / 0.29 | **0.65 / 0.34–0.39** | 0.63 / 0.56–0.57 |
| Label transfer A→B / B→A (balanced acc.) | 0.78 / 0.65 | **0.88 / 0.88** | 0.81 / 0.55 |
| Tumor/adjacent structure, stroma A / B (ASW) | 0.26 / 0.23 | 0.20 / 0.21 | **0.00 / 0.05** |

\* ARI/NMI of Leiden clusters on the Harmony embeddings vary slightly between runs (ARI 0.39 vs
0.34 in the clean rebuild). Harmony's output differs only at the ~1e-9 level, which is enough to change
near-tied neighbours, and Leiden amplifies that. All other metrics are identical across runs, and the
report does not use ARI.

- *Decision (approved):* **Harmony with batch = cohort** is the integrated embedding. Batch =
  specimen is rejected: specimen is confounded with condition, and it erased the tumor/adjacent
  structure in stroma and endothelium, exactly the signal Stage 5 tests.

**3. What was preserved, weakened or lost?**
- *Preserved:* every major population stays distinct (kNN purity 0.89–1.00); tumor vs adjacent
  structure within stroma and endothelium (reduced ~10–20 %, not erased).
- *Weakened:* the malignant vs non-malignant epithelium boundary across cohorts (both derive from
  proximal tubule); global cohort mixing stays low (median iLISI ≈ 1) because the two cohorts have
  very different compositions (A mostly tumor nuclei, B mostly immune cells) and assays; Harmony
  aligns matching identities rather than forcing everything to mix.
- *Lost:* nothing with batch = cohort; with batch = specimen, the condition structure.

**4. What does combining cohorts add?**
- B's stromal cells (fibroblasts + pericytes) are predicted as stroma from A's labels in 99 % of
  cases after Harmony (baseline: 69 %), supporting a **combined fibroblast + mural compartment**
  across cohorts.
- Paired donors for that compartment rise from 3 (A) to 12 (A + B) at ≥ 10 cells.
- B adds an explicit pericyte label and an independent single-cell assay to test whether a
  tumor-vs-adjacent effect seen in nuclei replicates.

### Stage 4: where is the ECM panel expressed, and which population to test?

Figures: [`outputs/step_4_output/`](outputs/step_4_output/) (regenerate with
`uv run opengo overview`, `uv run opengo figures 4`). Values: cell level = log1p(CP10k) and
fraction expressing (overview only); donor level = pseudobulk log2-CPM per specimen (≥ 10 cells),
library size from all genes of the cohort.

**1. Where is the program expressed?**
- *Evidence:* `01_ecm_dotplot.png`, `02_ecm_donor_level.png`.
- *Observation:* in both cohorts the stromal compartment has the highest expression of every panel
  gene (median log2-CPM 6–9 in A, 4–10 in B). Endothelial cells express SPARC and FN1 strongly
  (log2-CPM 7–10) but fibrillar collagens weakly and DCN/LUM not at all. Immune and epithelial
  populations are low throughout.
- *Within stroma:* DCN and LUM are **fibroblast-specific** (B fibroblasts log2-CPM ≈ 9; B pericytes
  3.5 and 0; A "Fibroblasts", which are mostly mural cells, ≈ 0). The collagens, FN1 and SPARC are
  shared by fibroblasts and mural cells.

**2. Which summary of the 10 genes?**
- *Check:* do the genes co-vary across specimens within stroma? Gene–gene Spearman, PC1, and
  agreement between candidate summaries (`03_panel_coherence.png`).
- *Evidence:* median pairwise ρ 0.48 (A) and 0.70 (B); PC1 explains 49 % / 64 %. The mean of
  per-gene z-scores agrees with PC1 (ρ 1.00 / 0.97) and with log of summed CPM (0.96 / 0.98).
- *Decision:* **mean of per-gene z-scored log2-CPM** (z within cohort × population): equal gene
  weights, simple, and equivalent to the alternatives here.
- *Genes that challenge it:* in B, **DCN and LUM do not co-vary with the other eight** (PC1 loadings
  ≈ 0; ρ with the rest ≈ −0.3 to 0.5). Because they mark fibroblasts rather than mural cells, they
  track the fibroblast-to-pericyte mix inside the compartment. They are always reported separately,
  and the summary is also reported without them.

**3. Which focus population? (approved at Check-in 4)**
- *Evidence:* `04_focus_population_choice.png`.

| Candidate | Panel expression | Paired donors ≥ 10 cells (A + B) | Verdict |
|---|---|---|---|
| **Stroma: fibroblast + mural** | highest | 3 + 9 = **12** | **chosen** |
| Fibroblast/myofibroblast only | highest | 3 + 2 = 5 | too few pairs in B |
| Pericyte only | high | 0 + 9 = 9 | B only; sensitivity |
| Endothelial | SPARC/FN1 only | 3 + 9 = 12 | weak for collagens/DCN/LUM; contrast only |

- *Justification:* the stromal compartment is where the panel is expressed, it is supported as one
  compartment across cohorts by markers (Stage 2) and label transfer (Stage 3, 99 %), and it has
  the best paired coverage.
- *Known risk, carried into Stage 5:* a tumor-vs-adjacent difference in this compartment can come
  from a change in its fibroblast:mural mix rather than a change in expression per cell. Stage 5
  therefore reports the pericyte-only (B) result and the fibroblast fraction within the
  compartment alongside the main estimate.

### Stage 5: tumor vs adjacent in the stromal compartment

Figures: [`outputs/step_5_output/`](outputs/step_5_output/) (regenerate with
`uv run opengo pseudobulk`, `uv run opengo figures 5`).

**What is tested.** Pseudobulk raw counts of QC-passing stromal (fibroblast + mural) cells per
specimen → log2-CPM (library size = all genes of the cohort). **Donors are the replicates:**
repeated tumor regions are averaged within donor; the primary analysis is paired (tumor − adjacent
within donor); specimens with < 10 compartment cells are dropped. Panel summary = mean of per-gene
standardized differences (the paired form of the Stage 4 mean-z score). Per cohort: mean, 95 %
t-interval, exact sign-flip p (all 2ⁿ sign patterns; with n = 3 the smallest possible p is 0.25),
BH correction over the 10 genes. Cohorts combined by inverse-variance meta-analysis
(DerSimonian–Laird random effects, I²). The integration embedding is not used.

**Observations** (`01_paired_panel_score.png`, `02_gene_effects.png`)

| | n paired donors | Panel score, tumor − adjacent [95 % CI] | Donors up | Sign-flip p |
|---|---|---|---|---|
| Cohort B (scRNA) | 9 | **+1.04 [+0.63, +1.46]** | 9 / 9 | 0.004 |
| Cohort A (snRNA) | 3 | +0.02 [−1.77, +1.81] | 1 / 3 | 1.0 |
| A + B, fixed effect | 12 | +0.88 [+0.56, +1.21] | 10 / 12 | 0.003 (pooled sign-flip) |
| A + B, random effects (k = 2; τ² poorly estimated) | 12 | +0.60 [−0.39, +1.59] | | |

- *B, per gene:* COL1A2, COL3A1, FN1 and SPARC increase in all 9 donors (log2FC +1.7 to +2.6,
  BH q = 0.01); COL1A1, COL6A1, COL6A2 also increase (q ≤ 0.03); COL5A1 is borderline (q = 0.09).
- *Genes that challenge the summary:* **DCN and LUM do not increase** (B log2FC −1.7 and −0.7,
  q ≥ 0.17) and vary wildly between donors (±5 log2 units), consistent with their being
  fibroblast-specific in a compartment of which ~97 % of cells are author-labelled pericytes in B.
- *A:* no consistent direction; one donor up, two down; intervals span ±2.

**Robustness check, targeted at the main uncertainty, compartment composition**
(`03_robustness.png`). If the effect came from more fibroblasts in tumor stroma, it should vanish in
pure mural cells and without the fibroblast genes, and the fibroblast share should rise in tumors.
- Pericytes only (B): **+1.27 [+0.85, +1.70], 9 / 9 up.**
- Panel without DCN/LUM (B): **+1.41 [+0.96, +1.86], 9 / 9 up.**
- Fibroblast share of the compartment, tumor − adjacent (B): +0.04 [−0.06, +0.14]; only 2 / 9 donors
  up (the mean is driven by RCC-PR1).
- → composition does not explain the B effect.

**Sensitivities** (same figure): B is unchanged with ≥ 20 or ≥ 50 cells (+1.04, +1.07), with one
tumor region per donor (+1.04), in a mixed model using all 26 specimens with a donor random
intercept (+1.11 [+0.85, +1.36]), and when any single donor is left out (+0.93 to +1.15). A stays
null in a mixed model that adds the unpaired tumor-only donors (24 donors, 32 specimens; −0.18
[−0.85, +0.48]). Those donors add only between-donor information; the contrast still rests on 4
adjacent specimens. One of A's 25 donors has < 10 compartment cells and is dropped.
*Specificity:* endothelial cells show no change (B +0.13 [−0.26, +0.52]), so this is not a global
tumor shift in ECM transcripts.

**Interpretation (separate from the observations above)**
- In single-cell data (B), tumor-associated mural/pericyte cells express more fibrillar collagen,
  FN1 and SPARC transcript than their adjacent counterparts in the same patients; the fibroblast
  proteoglycans DCN/LUM do not follow. Higher transcript levels are not, by themselves, evidence of
  active matrix remodelling.
- **Cohort A neither confirms nor contradicts this.** Its CI (−1.77 to +1.81) contains B's estimate:
  3 paired donors and 4 adjacent specimens are too few. It also uses a different assay (nuclei) and has
  a compartment of different make-up (fibroblast-like + mural, versus mostly pericytes in B). The
  apparent heterogeneity (I² = 80 %) is not reliable evidence of disagreement: it comes from two
  cohorts, treats A's n = 3 standard error as known, and each cohort is standardized by its own SD.
  The B result is robust within B; whether it generalises to snRNA data and to fibroblasts is not
  established here.

### Stage 6: linked modalities: snRNA vs bulk RNA vs protein

Figures: [`outputs/step_6_output/`](outputs/step_6_output/) (regenerate with
`uv run opengo multimodal`, `uv run opengo figures 6`). Bulk values are used as supplied (log2,
already normalised); missing protein values are not imputed (none occur for the 10 panel genes).

**1. Is the snRNA ↔ bulk correspondence what the crosswalk says?** (`01_crosswalk_validation.png`)
- *Checks per snRNA specimen:* link given → resolves through `bulk_samples.tsv` → column exists in
  the condition-specific RNA and protein file → protein aliquot equals the authors' sample table →
  snRNA parent specimen equals protein parent specimen.
- *Evidence:* 29 / 34 linked; all 29 resolve to real RNA and protein columns; all protein aliquots
  match the author table; the 2 parent mismatches are exactly the 2 flagged links. 27 primary links
  (23 tumor, 4 adjacent), as documented.
- *Matched subset:* donors with primary-linked tumor **and** adjacent = **C3L-00079, C3L-00088,
  C3N-01200 (n = 3)**. C3N-00242 (tumor protein from two parent specimens) is sensitivity only;
  C3L-00908 (different parent) has no adjacent anyway. The 5 unlinked snRNA specimens get no bulk
  partner: no pairing is invented from donor identity. Bulk RNA is a donor/condition column whose
  exact aliquot is not reconciled; each assay used a different portion of the parent specimen.

**2. Do the assays agree on direction?** (`02_cross_assay_agreement.png`, `03_matched_subset.png`;
each effect estimated within its own assay, tumor − adjacent, paired by donor)

| Evidence | Panel score | Collagens / FN1 / SPARC | DCN / LUM |
|---|---|---|---|
| B scRNA stroma (9 pairs, Stage 5) | **+1.04** | up (7/8 genes, q ≤ 0.03) | not up |
| A snRNA stroma (3 pairs) | +0.02 | mixed | down / mixed |
| Matched subset, 3 donors: snRNA stroma / whole specimen / bulk RNA / protein | −0.22 / −0.14 / +0.48 / +0.26 | FN1, SPARC up in all four; collagens mixed | down |
| **Bulk RNA, all CPTAC pairs (72)** | **+0.52 [0.31, 0.73], 54/72 up** | up (FN1 72/72, SPARC 72/72) | **down** (DCN 4/72 up) |
| **Protein, all CPTAC pairs (80)** | **+0.68 [0.49, 0.87], 67/80 up** | up (SPARC 78/80, COL3A1 72/80) | DCN down, LUM ≈ 0 |

- *Observation:* the well-powered bulk cohorts and the single-cell cohort B show the **same gene
  pattern** (fibrillar collagens, FN1 and SPARC higher in tumor, DCN/LUM lower or unchanged) at
  both RNA and protein level. COL6A1 is the one discrepancy (down in bulk RNA, up in protein and B).
- *Observation:* the 3-donor matched subset cannot settle direction (donors disagree, all p ≥ 0.5).
  Adding the flagged C3N-00242 link (sensitivity) moves bulk RNA from +0.48 to +0.07 and protein
  from +0.26 to −0.17, which reinforces this;
  on average FN1 and SPARC are higher in tumor in every matched assay (not in every donor). The
  matched snRNA-stroma value (−0.22) differs from Stage 5's A estimate (+0.02) because Stage 6 uses
  only the crosswalk-linked tumor specimen, whereas Stage 5 averages all tumor regions per donor.

**3. Could tissue composition explain the bulk signal?** (`04_composition_proxy.png`; non-panel
marker genes in bulk RNA, 72 pairs)
- Fibroblast markers (PDGFRA, C7, FBLN1) **fall** in 97 % of tumors (mean −1.54 SD); mural (+0.56)
  and endothelial (+1.63) markers rise, the hypervascular ccRCC pattern.
- Across donors the ECM change tracks fibroblast (ρ 0.77; 0.68 without DCN/LUM) and mural (ρ 0.62)
  marker change: composition strongly shapes the bulk signal.
- Yet the 8-gene panel (without DCN/LUM) rises on average (+0.96 SD) **despite** the fibroblast
  loss, so the bulk increase is not explained by more fibroblasts. It may still reflect the
  expansion of mural and endothelial compartments (both rise; ECM tracks mural change, ρ 0.62).
  Bulk data cannot separate per-cell change from composition. (A regression extrapolated to
  "no marker change" gives +1.46; shown but not relied on, because almost no donor has zero change.)
- In snRNA, the matched specimens' stromal fraction differs strongly between tumor and adjacent
  (e.g. C3L-00079 8 % vs 1.5 %); snRNA captures a different aliquot and nuclei capture is biased, so
  it is not used as a quantitative composition estimate for bulk.

**Interpretation and limitations**
- Bulk RNA, protein and single-cell B agree that collagen I/III, FN1 and SPARC are higher in
  ccRCC tissue than in adjacent kidney while DCN/LUM are not. Single-cell B places at least part of
  the RNA increase within mural/pericyte cells (Stage 5). **Cell-specific RNA does not establish the
  cellular source of bulk protein**: endothelial cells, which expand in tumors, also express FN1 and
  SPARC strongly (Stage 4).
- The matched snRNA–bulk comparison (n = 3) is too small to confirm or refute anything; the
  full-cohort bulk analysis is context, not a matched validation. It comes from the same CPTAC study
  as cohort A, its 72/80 pairs include the matched donors, and most donors have no snRNA data.
- Bulk values were used as supplied; the panel genes have no missing values, but upstream
  imputation cannot be ruled out (no imputation flag is provided).
- Composition is a major driver of bulk differences; the marker-based proxy is approximate and
  RNA-based.

## Parameters, seeds and records

- Parameters and the random seed (`RANDOM_SEED = 42`) are in
  [`src/opengo/config.py`](src/opengo/config.py).
- Every inclusion/exclusion and analysis decision is logged with its evidence in
  [`outputs/tables/decisions_log.tsv`](outputs/tables/decisions_log.tsv).

## Repository layout

```
src/opengo/          analysis package: config, io (loaders), genes, qc, explore, integrate,
                     stage4-6, plots_stage1-6, plots_report, report, cli
tests/               pytest checks (counts, joins, label override, Harmony not a no-op,
                     exact sign-flip p, crosswalk, 4-page report)
outputs/tables/      every number in the report + decisions_log.tsv
outputs/step_N_output/  per-stage diagnostic figures
outputs/logs/        run logs
outputs/checkpoints/ large intermediates (regenerated, not committed)
notebooks/           analysis notebook (jupytext pair, executed)
report/              4-page report (HTML source, figures, PDF)
PLAN.md              analysis plan, review gates, change log
```

Development checks:

```bash
uv run ruff check src tests
```

```bash
uv run pytest
```

## AI assistance

**Two separate AI systems, with different roles:**

- **Working agent: Claude Code** (Anthropic, Claude Opus 5.5). It planned the analysis with the
  author, wrote the code, ran the pipeline and drafted the text.
- **Gatekeeper: Gemini** (Google). The author used it as an independent reviewer, kept separate from
  Claude and unbiased by its reasoning. At every stage it evaluated the conclusions and the QC before
  the work moved on. Where it disagreed, the disagreement was settled by data, not by preference. For
  example, it challenged the first TRM relabelling, which was then decided by lineage markers. Its
  advice on population choice, statistics and the composition bridge to bulk data is recorded in
  `PLAN.md` (rows marked "Gemini review").
- **Fresh-context reviews:** in addition, two independent Claude sub-agents with no access to the
  working conversation reviewed the finished work. One checked compliance with the assignment and
  every number in the report against the tables; the other scored the final submission. Their
  findings were applied.

**Workflow.** Work proceeded in stages, each closed by a gate:
- Claude stopped and reported what was run, the key numbers and figures, any surprises, and its
  proposed decisions with alternatives.
- Gemini independently evaluated the stage's conclusions and QC.
- The author then approved the stage, redirected it, or asked for a redo before the next stage
  started.

These steps are recorded in the `PLAN.md` change log and in `decisions_log.tsv`.

**How consequential outputs were checked**

- *Inputs:* SHA-256 of all 70 supplied files against the package manifest.
- *Joins and counts:* pytest assertions against the package documentation: cells and specimens
  per cohort, every author annotation matched, integer counts, ECM genes present by Ensembl ID,
  crosswalk counts (29 linked, 27 primary), no multimodal pair invented from donor identity.
- *Figures:* every figure was inspected before a decision was based on it; several were redrawn
  after inspection.
- *Tool failures caught by checking, not assumed away:*
  - harmonypy 2.1.0 silently made no correction. It was detected because the metrics equalled the
    baseline, confirmed against 0.0.10 on identical input, then pinned and guarded by a test.
  - scanpy's neighbour search failed on a sparse type above ~8k cells; a wrapper was added after
    reproducing the failure.
- *Gatekeeper challenges settled by data:* e.g. the TRM relabelling was challenged by the Gemini
  review and then decided by
  lineage-marker expression, not by label text.
- *Statistics:* the exact sign-flip minimum p is unit-tested. With n = 3 pairs no p below 0.25 is
  possible, so effect sizes and confidence intervals lead.
- *Report:* every number in the report comes from `outputs/tables/`. The PDF build fails if the
  report exceeds four pages.

## Time spent

| Date | Time | Activity | Hours (approx.) |
|---|---|---|---|
| 2026-11-07 | 11:00–14:00 | Reading the assignment and data documentation, best-practice review, plan, uv environment, CLI, GitHub setup | 3.0 |
| 2026-11-07 | 14:00–16:00 | Stages 1–3: loading and reconciliation, QC, cohort C assessment, integration (incl. two tool issues) | 2.0 |
| 2026-11-07 | 16:00–16:45 | Stages 4–6: ECM overview, tumor vs adjacent, robustness, multimodal comparison | 0.75 |
| 2026-11-07 | 16:45–17:00 | Stage 7: report, notebook, README | 0.25 |
| | | **Total** | **≈ 6** |

Split by activity estimated from commit timestamps; the total is the author's own record (11:00–17:00).
