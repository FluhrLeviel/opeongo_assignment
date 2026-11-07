"""Stage 1 figures: study design, load reconciliation, composition, coverage, gene matching.

Reads only Stage 1 tables (plus cohort A obs for the label cross-tab) and writes PNGs to
outputs/step_1_output/. Run with `uv run opengo figures 1`.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

from opengo import config  # noqa: E402

OUT_DIR = config.OUTPUTS_DIR / "step_1_output"
COND_COLORS = {"Tumor": "#c0392b", "Adjacent": "#2471a3"}
POP_ORDER = [
    "Malignant_epithelial",
    "Nonmalignant_epithelial",
    "Endothelial",
    "Stromal_fibro_myo",
    "Pericyte",
    "Myeloid",
    "Lymphoid",
    "Unresolved",
    "Unlabelled",
]
POP_COLORS = dict(
    zip(
        POP_ORDER,
        [
            "#7b241c",
            "#e59866",
            "#1f618d",
            "#196f3d",
            "#82e0aa",
            "#7d3c98",
            "#d4ac0d",
            "#808b96",
            "#d5d8dc",
        ],
        strict=True,
    )
)


def _tables() -> dict[str, pd.DataFrame]:
    t = config.TABLES_DIR
    man = pd.read_csv(config.SAMPLE_MANIFEST_PATH, sep="\t")
    spec = pd.read_csv(t / "stage1_specimens.tsv", sep="\t")
    spec = spec.merge(man[["sample_id", "donor_id"]], on="sample_id")
    cov = pd.read_csv(t / "stage1_population_by_specimen.tsv", sep="\t")
    cov = cov.merge(man[["sample_id", "donor_id"]], on="sample_id")
    return {
        "spec": spec,
        "cov": cov,
        "genes": pd.read_csv(t / "gene_correspondence.tsv", sep="\t", low_memory=False),
    }


def _save(fig: plt.Figure, name: str) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / name
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return path


# ------------------------------------------------------------------ 1. study design


def fig_design(spec: pd.DataFrame) -> Path:
    """Donor x specimen layout: which donors have tumor, adjacent, repeated regions."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 7.5), gridspec_kw={"width_ratios": [25, 11, 8]})
    for ax, (cohort, d) in zip(axes, spec.groupby("cohort"), strict=True):
        counts = d.groupby(["donor_id", "condition"]).size().unstack(fill_value=0)
        for c in ("Tumor", "Adjacent"):
            if c not in counts:
                counts[c] = 0
        counts["paired"] = (counts["Tumor"] > 0) & (counts["Adjacent"] > 0)
        counts = counts.sort_values(["paired", "Adjacent", "Tumor"], ascending=False)
        y = np.arange(len(counts))
        ax.barh(y, counts["Adjacent"], color=COND_COLORS["Adjacent"], label="Adjacent")
        ax.barh(
            y, counts["Tumor"], left=counts["Adjacent"], color=COND_COLORS["Tumor"], label="Tumor"
        )
        ax.set_yticks(y, [f"{i}{' ●' if p else ''}" for i, p in counts["paired"].items()])
        ax.invert_yaxis()
        ax.set_xlabel("specimens")
        ax.xaxis.get_major_locator().set_params(integer=True)
        n_p = int(counts["paired"].sum())
        assay = d["assay"].iloc[0]
        ax.set_title(f"Cohort {cohort} ({assay})\n{len(counts)} donors, {n_p} paired ●")
        ax.tick_params(axis="y", labelsize=8)
    axes[0].legend(loc="lower right")
    fig.suptitle(
        "Study design: specimens per donor (● = donor has both tumor and adjacent)",
        fontweight="bold",
    )
    fig.tight_layout()
    return _save(fig, "01_study_design.png")


# ------------------------------------------------------------------ 2. reconciliation


def fig_reconciliation(spec: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))

    ax = axes[0]
    for cohort, d in spec.groupby("cohort"):
        ax.scatter(d["manifest_cells"], d["cells_loaded"], label=cohort, s=30, alpha=0.8)
    lim = spec[["manifest_cells", "cells_loaded"]].to_numpy().max() * 1.05
    ax.plot([0, lim], [0, lim], "k--", lw=0.8)
    ax.set(
        xlabel="cells expected (sample_manifest)",
        ylabel="cells loaded",
        title="Cells per specimen: loaded vs expected\n(all 73 on the identity line)",
    )
    ax.legend(title="cohort")

    ax = axes[1]
    a = spec[spec["cohort"] == "A"].sort_values("raw_barcodes")
    ax.bar(np.arange(len(a)), a["raw_barcodes"], color="#bbbbbb", label="raw droplet barcodes")
    ax.bar(np.arange(len(a)), a["cells_loaded"], color="#1f618d", label="author-annotated nuclei")
    ax.set_yscale("log")
    ax.set(
        xticks=[],
        xlabel="34 cohort A specimens",
        ylabel="barcodes (log)",
        title="Cohort A: annotated nuclei kept from raw droplet matrices",
    )
    ax.legend(loc="center right", fontsize=8)

    ax = axes[2]
    spec = spec.assign(mean_umi=spec["total_counts"] / spec["cells_loaded"])
    order = spec.sort_values(["cohort", "condition", "mean_umi"])
    colors = order["condition"].map(COND_COLORS)
    ax.bar(np.arange(len(order)), order["mean_umi"], color=colors)
    for cohort, d in order.reset_index(drop=True).groupby("cohort"):
        ax.text(d.index.to_numpy().mean(), ax.get_ylim()[1] * 0.95 or 1, cohort, ha="center")
    ax.set_yscale("log")
    ax.set(
        xticks=[],
        xlabel="specimens (grouped by cohort, condition)",
        ylabel="mean counts per cell (log)",
        title="Sequencing depth differs by cohort\n(preview of Stage 2 QC)",
    )
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in COND_COLORS.values()]
    ax.legend(handles, COND_COLORS.keys(), fontsize=8)
    fig.tight_layout()
    return _save(fig, "02_load_reconciliation.png")


# ------------------------------------------------------------------ 3. composition


def fig_composition(cov: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(
        1, 3, figsize=(18, 5.5), gridspec_kw={"width_ratios": [34, 26, 13]}, sharey=True
    )
    for ax, (cohort, d) in zip(axes, cov.groupby("cohort"), strict=True):
        wide = d.pivot_table(
            index=["condition", "donor_id", "sample_id"],
            columns="population",
            values="cells",
            fill_value=0,
        )
        wide = wide[[p for p in POP_ORDER if p in wide.columns]]
        frac = wide.div(wide.sum(axis=1), axis=0)
        frac = frac.sort_index(level=["condition", "donor_id"])
        bottom = np.zeros(len(frac))
        for pop in frac.columns:
            ax.bar(
                np.arange(len(frac)),
                frac[pop],
                bottom=bottom,
                color=POP_COLORS[pop],
                label=pop,
                width=0.9,
            )
            bottom += frac[pop].to_numpy()
        conds = frac.index.get_level_values("condition")
        ax.set_xticks(
            np.arange(len(frac)),
            [f"{d_}" for d_ in frac.index.get_level_values("donor_id")],
            rotation=90,
            fontsize=6,
        )
        for i, c in enumerate(conds):
            ax.get_xticklabels()[i].set_color(COND_COLORS[c])
        ax.set_title(f"Cohort {cohort}: {len(frac)} specimens")
    axes[0].set_ylabel("fraction of cells (author labels → suggested populations)")
    handles, labels = [], []
    for ax in axes:
        for h, lab in zip(*ax.get_legend_handles_labels(), strict=True):
            if lab not in labels:
                handles.append(h)
                labels.append(lab)
    fig.legend(handles, labels, loc="upper center", ncol=9, bbox_to_anchor=(0.5, 1.06))
    fig.text(0.5, -0.04, "x labels: donor (red = tumor, blue = adjacent)", ha="center")
    fig.tight_layout()
    return _save(fig, "03_population_composition.png")


# ------------------------------------------------------------------ 4. paired coverage


def fig_paired_coverage(cov: pd.DataFrame) -> Path:
    d = cov.groupby(["cohort", "donor_id", "condition", "population"])["cells"].sum()
    d = d.unstack("condition", fill_value=0)
    pops = [
        p
        for p in POP_ORDER
        if p in d.index.get_level_values("population") and p != "Malignant_epithelial"
    ]
    thresholds = [1, 10, 20, 50, 100]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), sharey=True)
    for ax, cohort in zip(axes, config.COHORTS, strict=True):
        mat = pd.DataFrame(0, index=pops, columns=thresholds)
        sub = d.xs(cohort, level="cohort")
        for t in thresholds:
            ok = sub[(sub.get("Tumor", 0) >= t) & (sub.get("Adjacent", 0) >= t)]
            vc = ok.index.get_level_values("population").value_counts()
            mat[t] = vc.reindex(pops, fill_value=0)
        sns.heatmap(mat, annot=True, fmt="d", cmap="Greens", vmin=0, vmax=9, cbar=False, ax=ax)
        ax.set(title=f"Cohort {cohort}", xlabel="min cells in BOTH tumor & adjacent")
    axes[0].set_ylabel("population")
    fig.suptitle(
        "Paired donors usable for a tumor-vs-adjacent comparison, by min-cells threshold\n"
        "(malignant cells omitted: no adjacent counterpart)",
        fontweight="bold",
    )
    fig.tight_layout()
    return _save(fig, "04_paired_donor_coverage.png")


# ------------------------------------------------------------------ 5. genes


def fig_genes(genes: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    ax = axes[0]
    combo = (
        genes["in_A"].map({True: "A", False: ""})
        + genes["in_B"].map({True: "B", False: ""})
        + genes["in_C"].map({True: "C", False: ""})
    )
    vc = combo.value_counts().sort_values()
    ax.barh(vc.index, vc.values, color="#5d6d7e")
    for i, v in enumerate(vc.values):
        ax.text(v, i, f" {v:,}", va="center", fontsize=8)
    ax.set(
        xlabel="genes (Ensembl ID, or B symbol if unmapped)",
        ylabel="present in cohorts",
        title="Gene overlap across cohorts\n(shared A∩B∩C = integration gene set)",
    )
    ax.set_xlim(0, vc.max() * 1.2)

    ax = axes[1]
    b = genes.loc[genes["in_B"], "B_match"].value_counts()
    b = b.rename(
        {
            "A_symbol": "unique A symbol",
            "C_symbol": "unique C symbol",
            "collision": "old+new name of one ID\n(both in B; left out)",
            "unmatched": "no match",
            "ambiguous": "symbol → several IDs",
        }
    )
    ax.barh(b.index, b.values, color=["#1e8449", "#52be80", "#e67e22", "#c0392b", "#7f8c8d"])
    for i, v in enumerate(b.values):
        ax.text(v, i, f" {v:,}", va="center", fontsize=8)
    ax.invert_yaxis()
    ax.set(
        xlabel="cohort B gene symbols",
        title="How B symbols were mapped to Ensembl IDs\n(all 10 ECM genes: unique A symbol)",
    )
    ax.set_xlim(0, b.max() * 1.25)
    fig.tight_layout()
    return _save(fig, "05_gene_correspondence.png")


# ------------------------------------------------------------------ 6. label check


def fig_label_check() -> Path:
    a = ad.read_h5ad(config.cohort_raw_path("A"), backed="r")
    obs = a.obs[["author_broad", "author_group", "population"]].astype(str)
    a.file.close()
    ct = pd.crosstab(obs["author_broad"] + "  →  " + obs["population"], obs["author_group"])
    ct = ct.loc[ct.sum(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(10, 9))
    sns.heatmap(
        np.log10(ct + 1),
        annot=ct,
        fmt="d",
        cmap="Blues",
        cbar_kws={"label": "log10(nuclei + 1)"},
        annot_kws={"fontsize": 7},
        ax=ax,
    )
    trm = [i for i, lab in enumerate(ct.index) if lab.startswith("TRM")]
    for i in trm:
        ax.add_patch(plt.Rectangle((0, i), ct.shape[1], 1, fill=False, ec="red", lw=2))
    ax.set(
        xlabel="author group column (Cell_group_w_epithelialcelltypes)",
        ylabel="author broad label  →  suggested population (label_mapping.tsv)",
        title="Cohort A label check: TRM (red) is grouped with Macrophages by the authors\n"
        "but mapped to Lymphoid by the suggested guide",
    )
    fig.tight_layout()
    return _save(fig, "06_cohortA_label_check.png")


def make_all() -> list[Path]:
    sns.set_theme(style="white", context="notebook")
    t = _tables()
    return [
        fig_design(t["spec"]),
        fig_reconciliation(t["spec"]),
        fig_composition(t["cov"]),
        fig_paired_coverage(t["cov"]),
        fig_genes(t["genes"]),
        fig_label_check(),
        fig_trm_markers(),
    ]


# ------------------------------------------------------------------ 7. TRM marker check

TRM_GROUPS = {
    "TRM": None,
    "Macrophages": 1500,
    "Monocytes": None,
    "CD8 CTL": None,
    "CD4 T-cells": 1500,
    "Tregs": None,
    "NK cells strong": None,
    "Tumor cells": 1500,
}
TRM_MARKERS = {
    "pan-immune": ["PTPRC"],
    "myeloid": ["CD68", "CD163", "C1QA", "C1QC", "CSF1R", "MRC1", "CD4"],
    "T cell": ["CD3D", "CD3E", "CD2", "CD8A", "IL7R"],
    "tissue-resident T": ["ITGAE", "CD69"],
    "NK": ["NKG7"],
}


def trm_marker_table() -> pd.DataFrame:
    """Fraction of nuclei with >0 UMI per marker, TRM vs reference author labels (cohort A)."""
    rng = np.random.default_rng(config.RANDOM_SEED)
    a = ad.read_h5ad(config.cohort_raw_path("A"), backed="r")
    labels = a.obs["author_broad"].astype(str).to_numpy()
    idx = []
    for g, n in TRM_GROUPS.items():
        ii = np.where(labels == g)[0]
        if n and len(ii) > n:
            ii = rng.choice(ii, n, replace=False)
        idx.extend(ii.tolist())
    sub = a[np.sort(idx)].to_memory()
    a.file.close()
    sym = sub.var["symbol"].astype(str).to_numpy()
    grp = sub.obs["author_broad"].astype(str).to_numpy()
    rows = []
    for lineage, genes in TRM_MARKERS.items():
        for m in genes:
            j = np.where(sym == m)[0]
            if len(j) != 1:
                continue
            x = np.asarray(sub.X[:, j[0]].todense()).ravel()
            for g in TRM_GROUPS:
                rows.append((lineage, m, g, int((grp == g).sum()), float((x[grp == g] > 0).mean())))
    out = pd.DataFrame(rows, columns=["lineage", "marker", "author_label", "n_nuclei", "frac_pos"])
    out.to_csv(config.TABLES_DIR / "stage1_TRM_marker_check.tsv", sep="\t", index=False)
    return out


def fig_trm_markers() -> Path:
    t = trm_marker_table()
    mat = t.pivot_table(index=["lineage", "marker"], columns="author_label", values="frac_pos")
    mat = mat.reindex(
        [(lin, m) for lin, ms in TRM_MARKERS.items() for m in ms if (lin, m) in mat.index]
    )[list(TRM_GROUPS)]
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.heatmap(
        mat * 100,
        annot=True,
        fmt=".0f",
        cmap="Purples",
        vmin=0,
        vmax=100,
        ax=ax,
        cbar_kws={"label": "% nuclei with ≥1 UMI"},
    )
    ax.set_yticklabels([f"{m}  ({lin})" for lin, m in mat.index], rotation=0)
    ax.add_patch(plt.Rectangle((0, 0), 1, len(mat), fill=False, ec="red", lw=2))
    ax.set(
        xlabel="cohort A author label",
        ylabel="",
        title="TRM nuclei express myeloid (CSF1R, CD4) but not T-cell markers (CD3D/E, CD2, CD8A)\n"
        "→ supports remapping TRM to Myeloid (tissue-resident macrophages)",
    )
    fig.tight_layout()
    return _save(fig, "07_TRM_marker_check.png")
