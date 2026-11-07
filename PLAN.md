# Working plan — ECM program in ccRCC

This is a living document. The analysis is stepwise and adaptive: each stage ends with a
**check-in** (below) and what we learn can change the next stage. Changes are logged at the bottom.

**Hard constraint:** ~8.6 GB RAM / 8 cores. ~330k cells will not fit as one dense-ish object.
Load specimens one at a time → one sparse raw-count `.h5ad` per cohort in `outputs/checkpoints/`;
integration on an HVG subset only; pseudobulk summed from raw counts per cohort.

---

## Working protocol — check-in after every stage

At the end of each stage Claude stops and reports:
1. **What was run** (command, parameters, runtime).
2. **Key numbers and 1–3 diagnostic figures/tables** (paths in `outputs/`).
3. **Surprises / problems** (anything that contradicts expectations or the metadata).
4. **Proposed decision(s)** with alternatives and what each would change downstream.
5. **Effect on later stages** — does the plan need to change?

Then we decide together: **proceed / redo with changes / revise plan**. No stage starts without
an explicit go. Every decision is appended to `outputs/tables/decisions_log.tsv`
(stage, decision, evidence, alternatives considered) — this doubles as the
inclusion/exclusion record required by the assignment.

Decisions that are deliberately **left open until we see data** (flagged ⟶ in the stages):
QC thresholds · whether C is integrated or used only for robustness · Harmony batch key ·
focus population · min cells per pseudobulk · repeated-region handling · robustness check.

---

## Stage 0 — Environment, CLI, repo
- `.venv` exists (uv, Python 3.11) but no `pyproject.toml`/`uv.lock` → convert to a uv project, pin
  installed packages. All installs via `uv add` (no pip).
- Code layout: logic in `src/opengo/` (importable), heavy steps through a CLI, exploration and
  report figures in a notebook that imports the same functions:
  `uv run opengo verify | load A|B|C | qc | integrate | pseudobulk | multimodal | all`
  (built with `click`, already installed; each command reads/writes checkpoints, logs params).
- Dev tools (dependency group `dev`): `ruff` (lint + format), `pytest` (tests for joins/counts),
  `ipykernel` (notebook kernel). Analysis additions: `igraph` + `leidenalg` (clustering for C).
  LISI comes from `harmonypy` (already installed). `pydeseq2` only if a gate calls for it.
- `git init`; `.gitignore` excludes `data/`, `.venv/`, `outputs/checkpoints/`.
- `opengo verify`: SHA-256 of every input vs `FILE_MANIFEST.tsv`.
- **Check-in 0:** env locked, CLI runs, hashes all match.

## Stage 1 — Load & reconcile (no filtering)
| Cohort | Load | Cell key | Cells kept at load |
|---|---|---|---|
| A snRNA | extract each tar to tmp, `read_10x_mtx` | `orig.ident` + `individual_barcode` | barcodes in author annotation (expect 141,950) |
| B scRNA | stream each CSV from the tar | `specimen_barcode` | all; join annotations (expect 157,881), report unmatched both ways |
| C scRNA | `read_10x_h5(genome='GRCh38')` | `sample_barcode` | all (expect 29,499) |
- Gene correspondence: A/C by Ensembl ID; B symbols mapped to Ensembl via A/C features; log
  duplicates/unmatched. Absent gene ≠ zero → integration on the shared gene set; the 10 ECM genes
  verified by Ensembl ID in every cohort.
- Checks (as `pytest` tests): integer counts, orientation, per-specimen cell counts vs
  `paper_annotated_cells`, donor/condition join against `sample_manifest.tsv`.
- **Check-in 1:** specimen table (cells, annotation match %, donor, condition, repeats).

## Stage 2 — Exploratory QC & structure (all populations, all cohorts)
- Per-specimen metrics: UMI, genes, % mito, % ribo. B source QC fields recomputed and compared.
  Missing acquisition diagnostics are not treated as pass/fail.
- Assay-aware criteria: nuclei (A) — mito should be ~0, high mito = ambient/stripping, not dying
  cells; kidney cells (B/C) — mito can be biologically high → per-specimen MAD-based outliers, not a
  fixed %. ⟶ thresholds decided from the plots.
- Coverage table: cells per specimen × major population incl. `Unresolved` (kept in coverage, not
  in contrasts).
- **Exploratory structure before cell-level work:**
  - per-specimen pseudobulk PCA (log-CPM, shared genes), coloured by cohort / condition / donor →
    how much does cohort dominate vs tumor/adjacent?
  - per-cohort cell-level PCA/UMAP coloured by author label, specimen, QC metrics.
- C has no labels → per-cohort Leiden + canonical markers (PTPRC, CD3E, CD68, VWF/PECAM1,
  RGS5/PDGFRB, DCN/LUM, CA9/NDUFA4L2, LRP2/UMOD …) for coarse labels; later cross-checked by kNN label
  transfer in the integrated space. Never use ECM panel genes alone to define fibroblasts (circular).
- **Check-in 2:** cell/specimen exclusions; ⟶ C in integration or only as robustness cohort.

## Stage 3 — Section 1: A first, then what B/C add; baseline vs Harmony
- Sequence the assignment asks for: **A alone** → **A+B** → **A+B+C** (if C passes Check-in 2).
- Combined HVG object; normalize_total → log1p; HVGs with `batch_key`; ECM genes kept in object
  (not forced into HVGs); PCA.
- **Baseline:** PCA → kNN → UMAP, no correction. **Integration:** Harmony on PCA.
  ⟶ batch key = cohort (primary) vs specimen — specimen is confounded with condition; decided by
  the evaluation.
- One concise comparison (table + figure):
  - mixing: iLISI by cohort; identity: cLISI + silhouette by major population; Leiden vs labels ARI/NMI;
  - biology: is tumor vs adjacent structure within one population still visible; do rare populations
    (fibroblasts, pericytes) stay distinct or get merged.
- Write up explicitly: **signals preserved / weakened / lost**, and **what combining adds**
  (donors, adjacent samples, population coverage, scRNA vs snRNA capture differences).
- The embedding is used only for visualisation and label transfer — **never for testing**.
- **Check-in 3.**

## Stage 4 — ECM overview across populations
- Per-gene dotplot (fraction expressing + mean log-norm) across major populations × cohort.
- **Summary (justified):** per pseudobulk sample, mean of per-gene z-scored log-CPM over the 10 genes
  (equal gene weights, interpretable as an average log-fold shift; z-scoring within cohort removes
  cohort scale). Cross-checked against `sc.tl.score_genes` at cell level. Genes that disagree with
  the summary are shown individually.
- ⟶ **Focus population** from a coverage table: donors with ≥ N cells in BOTH tumor and adjacent, per
  cohort. Candidates: fibroblast/stromal (main ECM source, but only 460 cells in B; A fibro + myofibro
  ≠ B pericytes), endothelial (good A+B coverage), pericyte (B only). Malignant has no adjacent
  counterpart → reported as non-estimable.
- **Check-in 4:** chosen population + justification.

## Stage 5 — Tumor vs adjacent in the focus population
- Values tested: **raw counts summed per specimen** (pseudobulk) → log-CPM (edgeR-style TMM or CPM).
  ⟶ min cells per pseudobulk (e.g. ≥ 20); repeated tumor regions averaged per donor (primary) — or
  donor random effect.
- Primary: paired within-donor differences (tumor − adjacent) for the panel score and each gene.
  - per cohort: mean difference, 95 % CI, per-donor spaghetti plot;
  - combined: `~ cohort/donor + condition` (donor nested in cohort) and, as a check, inverse-variance
    meta-analysis of per-cohort effects with heterogeneity → **does adding cohorts change it?**
  - small n (≈4–9 pairs per cohort): emphasise effect size + CI and sign consistency; exact
    sign-flip permutation p-value; BH across the 10 genes. Unpaired donors used only in a sensitivity
    model.
- ⟶ Robustness check aimed at the main uncertainty seen by then — likely candidates: leave-one-cohort /
  leave-one-donor-out; restrict to paired donors; vary min-cells; alternative score (rank-based);
  cell-number/composition sensitivity.
- **Check-in 5.**

## Stage 6 — Section 3: linked modalities (cohort A)
- Validate crosswalk: 34 A specimens → 27 primary (23 T, 4 Adj), 2 flagged (C3L-00908 T different parent;
  **C3N-00242 T multiple protein IDs**), 5 with no bulk link. Confirm every `bulk_sample_id` resolves
  via `bulk_samples.tsv` to a real column in the condition-specific file.
- ⚠ Consequence: of the 4 donors with snRNA tumor+adjacent, **only 3 have both conditions on primary
  links** (C3N-00242 tumor is flagged) → primary paired n = 3, n = 4 in sensitivity. Bulk RNA is
  donor/condition-level only (aliquot not reconciled).
- Within each assay, same direction question for the panel and each gene:
  (a) snRNA focus-population pseudobulk; (b) snRNA **whole-specimen** pseudobulk — a bridge that mixes
  composition the same way bulk does; (c) bulk RNA; (d) protein — matched subset, paired where possible.
- Context (clearly labelled *not* the matched subset): full bulk cohort has 72 paired RNA and 80 paired
  protein tumor/adjacent donors → much better powered direction estimate.
- Composition: stromal/endothelial fraction per specimen from snRNA, to explain bulk vs cell-type
  discrepancies. Protein missingness per gene reported; no zero-filling.
- **Check-in 6:** agreement table, discrepancies, limitations.

## Stage 7 — Report & repo
- Report ≤ 4 pages: observations separated from interpretations, uncertainty stated, negative results OK.
- README (inputs, `uv sync`, run order `opengo all` + notebook), params + seed (42) in `src/opengo/config.py`,
  decisions log, AI-assistance note (what was checked and how), time log.

---

## Stage 1 findings → changes to later stages (proposed at Check-in 1)
Figures: `outputs/step_1_output/` (`uv run opengo figures 1`).

| Finding | Best-practice basis | Change |
|---|---|---|
| TRM (510 A nuclei) mapped to Lymphoid, authors group it with Macrophages | Harmonise labels to a shared ontology/reference, log deviations | Remap TRM → Myeloid (pending OK); log as deviation from the guide |
| Stromal paired coverage thin: A ≤ 3, B ≤ 2 donors (≥ 10 cells) | Pseudobulk floor ≈ 10 cells / 1,000 counts (50 conservative) | Stage 4: focus population probably endothelial (A 3 + B 9 paired) or a mural/stromal compartment; decided from coverage + QC |
| B `Unresolved` = fibroblast/pericyte/endothelial label conflicts (2,087 cells) | Fibroblast vs pericyte not cleanly separable by single markers | Stage 2: re-check these with non-panel markers (fibroblast: PDGFRA, C7, FBLN1; pericyte: RGS5, MYH11, NOTCH3; endothelial: VWF, CDH5, RAMP2 — PECAM1 absent in B). **Never use DCN/LUM/COL genes to define cells (they are the panel).** |
| A has no pericyte label; A "Myofibroblasts" may contain mural cells | Same | Stage 2: check RGS5/MYH11 in A myofibroblasts before treating A stroma ≈ B fibroblasts |
| Composition differs strongly by cohort (B tumors mostly immune; A mostly tumor nuclei) | snRNA vs scRNA capture biases (Wu 2019; Denisenko 2020) | Never compare population proportions across cohorts; compare within-population expression only, cohort as a fixed term |
| Depth differs ~10× (C adjacent specimens ~300–1,000 counts/cell) | Per-specimen QC | Stage 2: per-specimen/assay-aware QC; C may be robustness-only |
| Gene names: B older annotation, 2,204 unmatched, 532 rename pairs left out | Ensembl ID first, symbol fallback, log ambiguities | Keep; check every marker gene used later is in the shared set |
| Paired n small (A 3–4, C 4, B 9) | Sign-flip minimum two-sided p = 2/2ⁿ (n = 3 → 0.25) | Lead with effect size + CI and direction consistency; pool across cohorts (cohort-stratified model + meta-analysis) |

## Stage 2 findings and second-review input → changes to later stages (Check-in 2)
| Source | Point | Change |
|---|---|---|
| Stage 2 QC | C non-mitochondrial counts capped at ≤ 409/cell; 76 % mito reads; DCN/LUM in 0.2 % of cells | Proposed: C excluded from integration + primary analysis; optional sensitivity only |
| Stage 2 QC | A/B already author-filtered; flags rare in stroma/endothelium | Proposed: drop flagged cells (union) in A and B |
| Stage 2 PCA | PC1 83 % = cohort (sn vs sc) | Within-cohort effects + cohort term / meta-analysis (unchanged plan, now evidenced) |
| Stage 2 stroma | Fibroblast/mural boundary blurred in A and B | Stage 4 candidate: combined fibroblast + mural compartment (A 3 + B 9 paired donors ≥ 10 cells); split by marker sub-label as sensitivity |
| Second review | ECM panel is fibroblast-driven; don't choose endothelium on counts alone | Agreed: ECM dot plot across populations first (Stage 4); if endothelium is tested, report sparse baseline collagen |
| Second review | Combined model needs cohort term | Agreed, refined: donor nested in cohort (keeps pairing) |
| Second review | Stromal fraction as bridge to bulk | Agreed with caveats (captured-nuclei proxy, different aliquot); add bulk-internal stromal marker genes in the full bulk cohort (72 pairs) |
| Second review | 3-figure report layout | Agreed in principle (Stage 7); add ambient-RNA limitation |
| Stage 2 fix | scanpy 1.11.5 kNN fails on pynndescent `csr_array` for > ~8k cells | `_CsrKnn` wrapper (same algorithm); scipy pinned < 1.17 (not the cause, kept for consistency of QC outputs) |

## Requirements traceability (assignment → plan)
| Assignment requirement | Stage |
|---|---|
| Assess all 3 cohorts, QC, criteria for cells/specimens/cohorts, all populations | 1, 2 |
| A cells identified via author annotations, no re-calling | 1 |
| Start with A, assess what B/C add | 3 |
| Unintegrated baseline vs one integration method; mixing + identity + biology | 3 |
| Signals preserved/weakened/lost; what combining adds; one concise comparison | 3 |
| Compact ECM overview across populations | 4 |
| Concise, justified summary; individual genes where they challenge it | 4, 5 |
| Choose one population with adequate donor + T/Adj coverage, justify | 4 |
| Direction, size, donor variability, effect of adding cohorts | 5 |
| Donors as replicates, repeated specimens, pairing; state tested values (not embedding) | 5 |
| One robustness check on the main conclusion | 5 |
| Validate crosswalk, matched subset, exclusions/ambiguity (aliquot ≠ specimen) | 6 |
| Direction vs bulk RNA and protein, within assay, paired | 6 |
| Agreement, discrepancies, limitations (composition, coverage, missingness) | 6 |
| ≤ 4-page report, figures/tables, uncertainty, observations vs interpretations | 7 |
| README, clean-session scripts, env spec, params/seeds, incl/excl records | 0, 7 |
| AI-assistance description, time spent | 7 |

## Known risks to watch
- B has few fibroblasts → focus population may be dominated by A if stromal.
- snRNA (A) vs scRNA (B/C) capture different cell types and genes (dissociation bias) → compare effects
  within cohort, don't pool raw counts across cohorts without a cohort term.
- C labels are ours, not the authors' → C results only as supporting evidence.
- Tumor/adjacent is confounded with specimen → specimen-level batch correction can erase the signal of
  interest (another reason not to test on the embedding).

## References (practice)
- Heumos et al. 2023 *Nat Rev Genet* — best practices; per-sample QC, ambient RNA, doublets.
- Squair et al. 2021 *Nat Commun* — pseudobulk with donor replicates.
- Luecken et al. 2022 *Nat Methods* (scIB) — batch removal vs bio conservation metrics.
- OSCA "single-nuclei processing"; Hippen et al. 2021 miQC — assay-aware mito QC.
- Clark et al. 2019 *Cell* — CPTAC ccRCC proteogenomics (bulk source).

## Plan changes log
| Date | Stage | Change | Reason |
|---|---|---|---|
| 2026-11-07 | all | Added check-in protocol, traceability table, CLI, multimodal n = 3 finding | Plan review |
| 2026-11-07 | 7 | Done: 4-page report (report/report.pdf), executed notebook, `opengo all`, README finalised (AI note, time log to confirm) | Stage 7 |
| 2026-11-07 | 6 | Stage 6 findings accepted | Check-in 6 |
| 2026-11-07 | 6 | Done: crosswalk validated (27 primary, 3 matched pairs + C3N-00242). Full bulk RNA +0.52 (54/72), protein +0.68 (67/80); same gene pattern as B scRNA (collagens/FN1/SPARC up, DCN/LUM not). Fibroblast markers fall in 97 % of tumors; ECM tracks composition but rises despite fibroblast loss. Matched n = 3 inconclusive | Check-in 6 |
| 2026-11-07 | 5 | Stage 5 findings accepted | Check-in 5 |
| 2026-11-07 | 5 | Done: B +1.04 [0.63, 1.46] 9/9 donors; A +0.02 (n=3); meta RE +0.60 [−0.39, 1.59], I² 80 %. Robust to composition (pericytes only, no DCN/LUM, fibroblast share unchanged) and to all sensitivities; endothelial null | Check-in 5 |
| 2026-11-07 | 4 | Stromal focus + mean-z approved | Check-in 4 |
| 2026-11-07 | 4 | Done: ECM highest in stroma (both cohorts); summary = mean z (agrees with PC1, ρ ≥ 0.97); DCN/LUM fibroblast-specific, discordant in B. Proposed focus: stroma (fibro + mural), 3 + 9 paired donors; sensitivity pericyte-only and fibroblast fraction | Check-in 4 |
| 2026-11-07 | 3 | Harmony by cohort approved | Check-in 3 |
| 2026-11-07 | 3 | Done: A, A+B baseline, Harmony(cohort), Harmony(specimen). harmonypy pinned 0.0.10 (2.1.0 no-op). Proposed: Harmony by cohort; specimen key rejected (erases condition). Combined stromal compartment supported by label transfer (99 %) | Check-in 3 |
| 2026-11-07 | 2 | Done: C excluded, flagged cells removed (approved) | Check-in 2 |
| 2026-11-07 | 1 | Done: 329,330 cells loaded (A 141,950 / B 157,881 / C 29,499), all counts match docs; 29,531 shared genes. Early signal: stromal paired coverage is thin (A ≤3, B ≤2 donors) — affects Stage 4 choice | Check-in 1 |
| 2026-11-07 | 0 | Done: uv project (96 locked pkgs), `src/opengo/` package, CLI with `verify`, pytest + ruff, git init (no commit). 70/70 inputs verified | Check-in 0 |
