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
| 2 | `uv run opengo qc` | QC tables and figures | planned |
| 3 | `uv run opengo integrate` | baseline vs Harmony comparison | planned |
| 4–5 | `uv run opengo pseudobulk` | ECM overview, tumor vs adjacent tests | planned |
| 6 | `uv run opengo multimodal` | crosswalk validation, bulk RNA / protein comparison | planned |
| fig | `uv run opengo figures 1` | `outputs/step_1_output/*.png` | done |
| 7 | notebook in `notebooks/` | report figures and numbers | planned |

Each command reads the supplied inputs or earlier checkpoints and can be re-run from a clean session.

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
  (e.g. specimens per cohort; cells per cohort to follow at Stage 1).
- Each analysis stage is reviewed by the author before the next stage starts (see `PLAN.md`).

This section is updated as the analysis proceeds.

## Time spent

| Date | Activity | Hours |
|---|---|---|
| 2026-11-07 | Reading materials, planning, environment and repository setup | _to fill_ |
