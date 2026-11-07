# ECM gene program in clear-cell renal cell carcinoma

Analysis for the Opeongo computational biology home assignment: where a fixed 10-gene
extracellular-matrix (ECM) panel is expressed across cell populations in three ccRCC single-cell
cohorts, whether it differs between tumor and adjacent non-tumor tissue in one population, whether
combining cohorts changes that, and whether linked bulk RNA and protein data agree.

Panel: COL1A1, COL1A2, COL3A1, COL5A1, COL6A1, COL6A2, FN1, SPARC, DCN, LUM.

> **Status:** work in progress. Stages 0–1 (environment, input verification, loading and
> reconciliation of all three cohorts) are done;
> the analysis stages are added one at a time after review. See [`PLAN.md`](PLAN.md) for the
> full plan, decision gates and change log.

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

| Step | Command | Output | Status |
|---|---|---|---|
| 0 | `uv run opengo verify` | `outputs/tables/input_verification.tsv` | done |
| 1 | `uv run opengo load all` (or `A` / `B` / `C`), then `uv run opengo genes` | raw-count checkpoints in `outputs/checkpoints/`; `outputs/tables/stage1_*.tsv`, `gene_correspondence*.tsv`, `ecm_panel_presence.tsv` | done (~10 min, peak ~1.9 GB RAM) |
| 1 fig | `uv run opengo figures 1` | `outputs/step_1_output/*.png` | done |
| 2 | `uv run opengo qc`, then `uv run opengo explore`, then `uv run opengo figures 2` | per-cell QC flags, pseudobulk checkpoint, `outputs/tables/stage2_*.tsv`, `outputs/step_2_output/*.png` | done (~20 min); decisions pending review |
| 3 | `uv run opengo integrate` | baseline vs Harmony comparison | planned |
| 4–5 | `uv run opengo pseudobulk` | ECM overview, tumor vs adjacent tests | planned |
| 6 | `uv run opengo multimodal` | crosswalk validation, bulk RNA / protein comparison | planned |
| 7 | notebook in `notebooks/` | report figures and numbers | planned |

Each command reads the supplied inputs or earlier checkpoints and can be re-run from a clean session.

## How decisions were made

Every stage follows the same loop: **spot an issue → decide how to test it → look at the evidence
(a table or figure) → act, and log the decision with its evidence**. Nothing downstream runs until
the stage has been reviewed. The full record, with alternatives considered, is in
[`outputs/tables/decisions_log.tsv`](outputs/tables/decisions_log.tsv); the rationale for each
stage's choices is in [`PLAN.md`](PLAN.md).

### Stage 1 — loading and reconciliation

Figures: [`outputs/step_1_output/`](outputs/step_1_output/) (regenerate with `uv run opengo figures 1`).

**1. Did we load exactly the cells the authors analysed?**
- *Issue:* cohort A is deposited as raw droplet matrices (6.8 M barcodes per specimen); B and C as
  processed matrices. A silent join error would change every downstream number.
- *Check:* join A barcodes to the author annotation by specimen + barcode; compare cells per
  specimen with the package manifest; pytest assertions on the documented totals.
- *Evidence:* `02_load_reconciliation.png` — all 73 specimens on the identity line; A 141,950,
  B 157,881, C 29,499 cells, every author annotation matched, all counts non-negative integers.
- *Action:* A restricted to annotated nuclei (no re-calling, as instructed); B and C keep all
  deposited cells until QC.

**2. Can genes be compared across cohorts?**
- *Issue:* A and C use Ensembl IDs from different releases; B has symbols only, from an older
  annotation. Matching by symbol alone can silently pick the wrong gene.
- *Check:* A↔C by Ensembl ID; each B symbol assigned only if it matches exactly one A (else C)
  symbol; collisions and ambiguities counted rather than guessed.
- *Evidence:* `05_gene_correspondence.png` — 29,531 genes shared by all three cohorts; 532 B
  symbols (old and new names of the same 266 genes, mostly lncRNA clone names) left out; all 10
  ECM genes matched by Ensembl ID in every cohort.
- *Action:* integration uses the shared set; a gene absent from a cohort is treated as unmeasured,
  never as zero.

**3. Is a suggested cell label wrong? (TRM)**
- *Issue:* the supplied `label_mapping.tsv` maps cohort A "TRM" (510 nuclei) to Lymphoid, but the
  authors' own grouping column places all 510 under Macrophages. "TRM" is ambiguous: tissue-resident
  *macrophages* or tissue-resident *memory T cells*.
- *First instinct, and why it was not enough:* the initial proposal was to remap from the label
  columns alone. A second review (another AI agent, consulted by the author) pointed out that labels
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
- *Evidence:* `01_study_design.png`, `04_paired_donor_coverage.png` — at ≥ 10 cells the
  fibroblast/myofibroblast compartment has only 3 paired donors in A and 2 in B, whereas
  endothelial cells have 3 (A) + 9 (B) and B pericytes 9.
- *Action:* no exclusion yet; this shapes Stage 2 (re-examine B's fibroblast/pericyte/endothelial
  label conflicts and A myofibroblasts with non-panel markers) and the Stage 4 population choice.
  Panel genes (DCN, LUM, collagens) are never used to define cells, to avoid circularity.

**5. Are the cohorts comparable in composition and depth?**
- *Evidence:* `03_population_composition.png` — A tumors are mostly malignant nuclei, B tumors
  mostly immune cells, C unlabelled; mean counts per cell differ ~10-fold (C adjacent lowest).
- *Interpretation:* consistent with known single-nucleus vs single-cell capture biases
  (Wu et al. 2019; Denisenko et al. 2020).
- *Action:* population proportions are never compared across cohorts; tumor-vs-adjacent effects
  are estimated within cohort, with cohort as a model term; QC thresholds are set per specimen and
  per assay.

### Stage 2 — quality control and exploratory structure

Figures: [`outputs/step_2_output/`](outputs/step_2_output/) (regenerate with `uv run opengo qc`,
`uv run opengo explore`, `uv run opengo figures 2`). Status: **proposals awaiting review**.

**1. Which cells are low quality — and is a flag technical or biological?**
- *Approach:* per-specimen, assay-aware flags (no fixed global cut-offs): counts/genes below
  median − 5 MAD, counts above median + 5 MAD, mitochondrial % above median + 3 MAD *and* above
  5 % (nuclei, where mitochondrial reads indicate ambient RNA) or 25 % (whole kidney cells, where
  proximal-tubule cells are naturally high), Scrublet doublets. Nothing is removed at this step.
- *Evidence:* `01_qc_distributions.png`, `02_flag_rates.png`. A and B were already filtered by
  their authors (A counts sit between 10³ and 10⁴; B's author QC fields match our recomputation
  exactly, `03_B_author_qc_vs_recomputed.png`). Flags are rare in the ECM-relevant populations
  (A stroma 1 %, A endothelium 4 %, B ≤ 1 %); they concentrate in A immune nuclei (high
  mitochondrial fraction, 10–22 %) and B non-malignant epithelium (12 %, proximal tubule biology).
- *Proposed action:* exclude flagged cells (union of flags) in A and B; no A or B specimen excluded.

**2. Is cohort C usable?**
- *Issue noticed:* C's median mitochondrial fraction is 48 % and ~175 genes are detected per cell.
- *Check:* non-mitochondrial counts per cell, mitochondrial fraction per specimen, ECM panel
  detection by cohort.
- *Evidence:* `08_cohortC_integrity.png` — every C cell has **≤ 409 non-mitochondrial counts**
  (a hard ceiling; A and B reach 10⁴–10⁵); mitochondrial reads are 76 % of all C counts; the two
  worst specimens (median 92–96 % mitochondrial) are the tumors of two of C's four paired donors;
  DCN and LUM are detected in 0.2 % of C cells. Marker-based labelling still recovers broad
  populations in C (`06_cohortC_marker_labels.png`, after validating the procedure at 81 % / 95 %
  agreement with author labels in A / B, `05_marker_label_validation.png`), but per-cell ECM
  signal is too sparse to support a comparison.
- *Proposed action:* exclude C from integration and from the primary ECM analysis; report the
  evidence; optionally use C (non-mitochondrial pseudobulk, two clean pairs) only as a labelled
  sensitivity check.

**3. What dominates variation between specimens?**
- *Check:* PCA of whole-specimen and endothelial pseudobulk, with R² of cohort / condition / donor
  for each component.
- *Evidence:* `04_pseudobulk_pca.png` — PC1 (83 %) is cohort (R² 0.99, i.e. single-nucleus vs
  single-cell); tumor vs adjacent appears only on later components, within cohort (endothelial PC4,
  R² 0.58).
- *Action:* confirms the plan — effects are estimated within cohort and combined with a cohort
  term or meta-analysis; never by pooling raw expression across cohorts.

**4. What is the stromal compartment made of?** (sets up the Stage 4 choice)
- *Check:* non-panel fibroblast (PDGFRA, C7, FBLN1, MFAP5), mural (RGS5, NOTCH3, MYH11, PDGFRB)
  and endothelial markers in A and B stromal labels, plus marker-based sub-clustering.
- *Evidence:* `07_stromal_identity.png` — A "Fibroblasts" are mostly mural-like (RGS5 50 %,
  PDGFRA 2 %); A "Myofibroblasts" are a fibroblast (C7 42 %) + mural mix; B "Fibroblast" is also
  RGS5-high; B fibroblast/pericyte label conflicts fall mostly with mural cells. In ccRCC the
  fibroblast–pericyte boundary is a continuum, not two clean types.
- *Implication:* a combined **fibroblast + mural (stromal) compartment** is defensible across A and
  B and has 3 (A) + 9 (B) paired donors at ≥ 10 cells — the same coverage as endothelium, but it is
  the compartment that actually produces the ECM panel. Decided at Stage 4 from the ECM dot plot.

## Parameters, seeds and records

- Parameters and the random seed (`RANDOM_SEED = 42`) are in
  [`src/opengo/config.py`](src/opengo/config.py).
- Every inclusion/exclusion and analysis decision is logged with its evidence in
  [`outputs/tables/decisions_log.tsv`](outputs/tables/decisions_log.tsv).

## Repository layout

```
src/opengo/      importable analysis code (config, verify, CLI; loaders and stages to follow)
tests/           pytest checks for joins, counts and metadata wiring
outputs/         tables and figures (checkpoints are regenerated, not committed)
notebooks/       report notebook (jupytext)
PLAN.md          analysis plan, decision gates, change log
```

Development checks:

```bash
uv run ruff check src tests
```

```bash
uv run pytest
```

## AI assistance

Claude Code (Anthropic, Claude Opus 5.5) was used to plan the analysis, write code and draft text.
How consequential outputs were checked:

- Input integrity: SHA-256 of all 70 supplied files against the package manifest.
- Joins and counts: pytest checks against counts stated in the package documentation
  (specimens and cells per cohort, annotation matches, integer counts, ECM panel presence).
- Each analysis stage is reviewed by the author before the next stage starts (see `PLAN.md`);
  every figure is inspected before a decision is taken.
- Contested decisions get a second opinion and an empirical test: e.g. the TRM relabelling was
  challenged by a second AI agent and settled by lineage-marker expression (see
  "How decisions were made").

This section is updated as the analysis proceeds.

## Time spent

| Date | Activity | Hours |
|---|---|---|
| 2026-11-07 | Reading materials, planning, environment and repository setup | _to fill_ |
