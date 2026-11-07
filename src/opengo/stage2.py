"""Stage 2 driver: QC metrics, pseudobulk, exploratory PCA, marker labels, stromal checks."""

from __future__ import annotations

import pandas as pd

from opengo import config, explore, qc

STAGE2_DIR = config.CHECKPOINTS_DIR / "stage2"


def run_qc(progress=None) -> dict[str, pd.DataFrame]:
    cells = pd.concat([qc.run_qc(c, progress=progress) for c in config.COHORTS])
    tables = qc.summarise(cells)
    if "author_doublet_score" in cells:
        b = cells[cells["cohort"] == "B"]
        tables["stage2_B_author_vs_recomputed"] = pd.DataFrame(
            {
                "metric": ["total_counts", "n_genes", "pct_mt", "doublet_score"],
                "author_field": ["TotalUMI", "gene.num", "MT.ratio (x100)", "doublest.score"],
                "spearman": [
                    b["total_counts"].corr(b["author_total_umi"], method="spearman"),
                    b["n_genes"].corr(b["author_n_genes"], method="spearman"),
                    b["pct_mt"].corr(100 * b["author_mt_ratio"], method="spearman"),
                    b["doublet_score"].corr(b["author_doublet_score"], method="spearman"),
                ],
                "exact_match_frac": [
                    (b["total_counts"] == b["author_total_umi"]).mean(),
                    (b["n_genes"] == b["author_n_genes"]).mean(),
                    ((b["pct_mt"] - 100 * b["author_mt_ratio"]).abs() < 0.01).mean(),
                    float("nan"),
                ],
            }
        )
    for name, df in tables.items():
        df.to_csv(config.TABLES_DIR / f"{name}.tsv", sep="\t", index=False)
    return tables


def run_explore(progress=None) -> dict[str, pd.DataFrame]:
    STAGE2_DIR.mkdir(parents=True, exist_ok=True)
    out: dict[str, pd.DataFrame] = {}

    # 1-2. pseudobulk + specimen-level PCA (whole specimen, and endothelial cells only)
    pb = explore.build_pseudobulk(progress=progress)
    for label, sel in {
        "ALL": pb.obs["population"] == "ALL",
        "Endothelial": (pb.obs["population"] == "Endothelial")
        & (pb.obs["n_cells"] >= config.PSEUDOBULK_MIN_CELLS),
    }.items():
        res = explore.pseudobulk_pca(pb[sel.to_numpy()].copy())
        res["scores"].to_csv(STAGE2_DIR / f"pca_scores_{label}.tsv", sep="\t")
        out[f"stage2_pseudobulk_pca_r2_{label}"] = res["r2"]

    # 3. marker-label validation on A and B subsamples, then C
    acc_rows = []
    for cohort in ("A", "B"):
        sub = explore._subsample(cohort, config.MARKER_VALIDATION_N)
        labelled, clusters = explore.marker_label(sub, explore.MARKERS)
        conf = explore.compare_to_author(labelled)
        conf.to_csv(config.TABLES_DIR / f"stage2_marker_vs_author_{cohort}.tsv", sep="\t")
        acc_rows.append(
            {"cohort": cohort, "cells": labelled.n_obs, "accuracy": explore.label_accuracy(conf)}
        )
    out["stage2_marker_label_accuracy"] = pd.DataFrame(acc_rows)

    c = explore._subsample("C", 0)
    labelled, clusters = explore.marker_label(c, explore.MARKERS)
    keep = ["sample_id", "condition", "leiden", "marker_label", "label_margin"]
    keep += [k for k in labelled.obs if k.startswith("score_")]
    obs_c = labelled.obs[keep].copy()
    obs_c[["umap_1", "umap_2"]] = labelled.obsm["X_umap"]
    obs_c.to_csv(STAGE2_DIR / "C_marker_labels.tsv", sep="\t")
    clusters.to_csv(config.TABLES_DIR / "stage2_C_marker_clusters.tsv", sep="\t")
    out["stage2_C_label_counts"] = (
        obs_c.groupby(["marker_label", "condition"], observed=True).size().unstack(fill_value=0)
    ).reset_index()
    dot = _marker_dot(labelled, "marker_label")
    dot.to_csv(STAGE2_DIR / "C_marker_dot.tsv", sep="\t", index=False)

    # 4. stromal sub-analysis (A: fibro/myofibro/endothelial; B: + pericytes and conflicts)
    stromal_pops = {"Stromal_fibro_myo", "Pericyte", "Endothelial", "Unresolved"}
    b_conflict = {"Fibroblast", "Pericytes", "Endothelial"}
    for cohort in ("A", "B"):

        def mask(obs, cohort=cohort):
            pop = obs["population"].astype(str)
            m = pop.isin(stromal_pops - {"Unresolved"})
            if cohort == "B":
                m |= (pop == "Unresolved") & obs["author_broad"].astype(str).isin(b_conflict)
            return m.to_numpy()

        sub = explore._subsample(cohort, 0, mask_fn=mask)
        labelled, _ = explore.marker_label(sub, explore.STROMAL_MARKERS, resolution=0.8)
        author = labelled.obs["author_broad"].astype(str)
        if cohort == "B":
            author = author.where(
                labelled.obs["population"].astype(str) != "Unresolved",
                "Unresolved: " + author + " / " + labelled.obs["author_detailed"].astype(str),
            )
        ct = pd.crosstab(author, labelled.obs["marker_label"])
        ct.to_csv(config.TABLES_DIR / f"stage2_stromal_{cohort}_author_vs_marker.tsv", sep="\t")
        dot = _marker_dot(labelled, "author_broad", genes=_flat(explore.STROMAL_MARKERS))
        dot.to_csv(STAGE2_DIR / f"stromal_{cohort}_marker_dot.tsv", sep="\t", index=False)

    for name, df in out.items():
        df.to_csv(config.TABLES_DIR / f"{name}.tsv", sep="\t", index=False)
    return out


def _flat(markers: dict[str, list[str]]) -> list[str]:
    return [g for gs in markers.values() for g in gs]


def _marker_dot(adata, group: str, genes: list[str] | None = None) -> pd.DataFrame:
    """Fraction expressing + mean log-norm expression per group (adata is log-normalised)."""
    import numpy as np
    import scipy.sparse as sp

    genes = [g for g in (genes or _flat(explore.MARKERS)) if g in adata.var_names]
    x = adata[:, genes].layers["counts"]
    x = x.toarray() if sp.issparse(x) else np.asarray(x)
    lx = adata[:, genes].X
    lx = lx.toarray() if sp.issparse(lx) else np.asarray(lx)
    rows = []
    for g in adata.obs[group].astype(str).unique():
        m = (adata.obs[group].astype(str) == g).to_numpy()
        for j, gene in enumerate(genes):
            rows.append(
                {
                    "group": g,
                    "gene": gene,
                    "n": int(m.sum()),
                    "frac": float((x[m, j] > 0).mean()),
                    "mean": float(lx[m, j].mean()),
                }
            )
    return pd.DataFrame(rows)
